/* input_ofd.c — OpenFDTD 互換 .ofd パーサ (音響ソルバーが使う部分集合)
 *
 * 書式は本家 OpenFDTD/sol/input_data.c と OpenFDTD-X/src/io/OfdIO.cpp が正。
 * 読むキー : xmesh/ymesh/zmesh, geometry, feed, point, frequency1/frequency2,
 * title。その他の未知キーは本家と同様に無視する (前方互換)。
 *
 * 音響での解釈 (v1) :
 *   - メッシュ : 単一の一様格子 dx = 最小刻み。非一様なら警告して一様化。
 *   - geometry : 材質によらず剛体 (法線速度 0) としてボクセル化。
 *   - feed     : 位置のみ使用 (方向・電圧は電磁用)。複数あれば #1 のみ。
 *   - point    : 受音点。行末の "# 名前" を WAV ファイル名に使う
 *                (本家トークナイザは余剰トークンとして無視するので互換)。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "acoustic.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <dirent.h>
#endif

#define MAXTOKEN 1024
#define MAXLINE  8192

static int ends_with(const char *s, const char *suf)
{
	size_t ls = strlen(s), lf = strlen(suf);
	return (ls >= lf) && (strcmp(s + ls - lf, suf) == 0);
}

/* working_dir 直下の唯一の .ofd を探す (.ofdx は除外)。
 * 0 個・複数個は曖昧なので非零終了 — 数値を捏造した合成 RIR は出さない。 */
int ac_find_input(ac_t *ac)
{
	char found[AC_NAME_MAX * 4] = "";
	int count = 0;
#ifdef _WIN32
	WIN32_FIND_DATAA fd;
	HANDLE h;
	char pat[AC_PATH_MAX + 8];
	snprintf(pat, sizeof(pat), "%s\\*", ac->workdir);
	h = FindFirstFileA(pat, &fd);
	if (h != INVALID_HANDLE_VALUE) {
		do {
			if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) &&
			    ends_with(fd.cFileName, ".ofd")) {
				count++;
				snprintf(found, sizeof(found), "%s", fd.cFileName);
			}
		} while (FindNextFileA(h, &fd));
		FindClose(h);
	}
#else
	DIR *d = opendir(ac->workdir);
	struct dirent *e;
	if (!d) {
		ac_err(ac, "cannot open working dir '%s'", ac->workdir);
		return 1;
	}
	while ((e = readdir(d)) != NULL) {
		if (ends_with(e->d_name, ".ofd")) {
			count++;
			snprintf(found, sizeof(found), "%s", e->d_name);
		}
	}
	closedir(d);
#endif
	if (count == 0) {
		ac_err(ac, "no input given and no .ofd file found in '%s'", ac->workdir);
		return 1;
	}
	if (count > 1) {
		ac_err(ac, "no input given and %d .ofd files found in '%s' "
		       "(ambiguous — pass the input file name explicitly)",
		       count, ac->workdir);
		return 1;
	}
	if (snprintf(ac->ofd_path, sizeof(ac->ofd_path), "%s/%s",
	             ac->workdir, found) >= (int)sizeof(ac->ofd_path)) {
		ac_err(ac, "input path too long");
		return 1;
	}
	return 0;
}

/* 空白区切りトークナイザ (本家 input_data.c と同じ流儀。str は破壊される) */
static int tokenize(char *str, char **tok, int maxtok)
{
	int n = 0;
	char *s = str;
	while (n < maxtok) {
		while (*s == ' ' || *s == '\t' || *s == '\r' || *s == '\n') s++;
		if (*s == '\0') break;
		tok[n++] = s;
		while (*s && *s != ' ' && *s != '\t' && *s != '\r' && *s != '\n') s++;
		if (*s) *s++ = '\0';
	}
	return n;
}

/* mesh 1 軸分 : "x0 d1 x1 d2 x2 ..." → 範囲と最小刻み。
 * 戻り値 0 = OK。min/max/dmin を更新する。 */
static int parse_mesh(ac_t *ac, char **tok, int ntoken,
                      double *lo, double *hi, const char *key)
{
	double prev, next, d;
	int i, ndiv, div;
	if (ntoken < 5 || ntoken % 2 == 0) {
		ac_err(ac, "invalid %s data", key);
		return 1;
	}
	ndiv = (ntoken - 3) / 2;
	prev = atof(tok[2]);
	*lo = prev;
	for (i = 0; i < ndiv; i++) {
		div  = atoi(tok[2 * i + 3]);
		next = atof(tok[2 * i + 4]);
		if (div <= 0 || next <= prev) {
			ac_err(ac, "invalid %s data (non-increasing nodes or bad division)", key);
			return 1;
		}
		d = (next - prev) / div;
		if (ac->dxmin <= 0.0 || d < ac->dxmin) {
			if (ac->dxmin > 0.0 && (ac->dxmin - d) > 1e-9 * ac->dxmin)
				ac->nonuniform = 1;
			ac->dxmin = d;
		}
		else if ((d - ac->dxmin) > 1e-9 * ac->dxmin) {
			ac->nonuniform = 1;
		}
		prev = next;
	}
	*hi = prev;
	return 0;
}

/* 行末 "# 名前" を取り出す (無ければ空文字)。line は '#' の手前で切られる */
static void take_trailing_name(char *line, char *name, size_t cap)
{
	char *hash = strchr(line, '#');
	size_t n;
	name[0] = '\0';
	if (!hash) return;
	*hash = '\0';
	hash++;
	while (*hash == ' ' || *hash == '\t') hash++;
	n = strlen(hash);
	while (n > 0 && (hash[n-1] == ' ' || hash[n-1] == '\t' ||
	                 hash[n-1] == '\r' || hash[n-1] == '\n')) n--;
	if (n >= cap) n = cap - 1;
	memcpy(name, hash, n);
	name[n] = '\0';
}

int ac_read_ofd(ac_t *ac)
{
	char line[MAXLINE], save[MAXLINE];
	char *tok[MAXTOKEN];
	int ntoken, nline = 0;
	int have_mesh[3] = { 0, 0, 0 };
	FILE *fp = fopen(ac->ofd_path, "r");

	if (!fp) {
		ac_err(ac, "cannot open input '%s'", ac->ofd_path);
		return 1;
	}

	while (fgets(line, sizeof(line), fp) != NULL) {
		if (!strncmp(line, "end", 3)) break;
		snprintf(save, sizeof(save), "%s", line);
		ntoken = tokenize(line, tok, MAXTOKEN);
		if (ntoken == 0) continue;

		if (nline == 0) {
			/* ヘッダ行 : "OpenFDTD <major> <minor>" (OpenTHFD / OpenRCWA も許容
			 * — OpenFDTD-X は光 RCWA 有効時にヘッダを切り替えるため) */
			if (strcmp(tok[0], "OpenFDTD") && strcmp(tok[0], "OpenTHFD") &&
			    strcmp(tok[0], "OpenRCWA")) {
				ac_err(ac, "not OpenFDTD/OpenTHFD data: %s", ac->ofd_path);
				fclose(fp);
				return 1;
			}
			nline++;
			continue;
		}
		/* "key = value..." 以外の行 (コメント等) は本家と同様に読み飛ばす */
		if (ntoken < 3 || strcmp(tok[1], "=")) continue;

		if (!strcmp(tok[0], "title")) {
			const char *eq = strchr(save, '=');
			size_t n;
			if (!eq) continue;
			eq++;
			while (*eq == ' ') eq++;
			n = strlen(eq);
			while (n > 0 && (eq[n-1] == '\r' || eq[n-1] == '\n')) n--;
			if (n >= sizeof(ac->title)) n = sizeof(ac->title) - 1;
			memcpy(ac->title, eq, n);
			ac->title[n] = '\0';
		}
		else if (!strcmp(tok[0], "xmesh")) {
			if (parse_mesh(ac, tok, ntoken, &ac->x0, &ac->x1, "xmesh")) { fclose(fp); return 1; }
			have_mesh[0] = 1;
		}
		else if (!strcmp(tok[0], "ymesh")) {
			if (parse_mesh(ac, tok, ntoken, &ac->y0, &ac->y1, "ymesh")) { fclose(fp); return 1; }
			have_mesh[1] = 1;
		}
		else if (!strcmp(tok[0], "zmesh")) {
			if (parse_mesh(ac, tok, ntoken, &ac->z0, &ac->z1, "zmesh")) { fclose(fp); return 1; }
			have_mesh[2] = 1;
		}
		else if (!strcmp(tok[0], "geometry")) {
			int shape, npar, n;
			ac_geom_t *g;
			if (ntoken < 4) {
				ac_err(ac, "invalid geometry data #%d", ac->ngeom + 1);
				fclose(fp);
				return 1;
			}
			shape = atoi(tok[3]);
			switch (shape) {
				case 1: case 2: case 11: case 12: case 13:
					npar = 6; break;
				case 31: case 32: case 33:
				case 41: case 42: case 43:
				case 51: case 52: case 53:
					npar = 8; break;
				default:
					npar = 0; break;   /* 未知 shape : setup で警告してスキップ */
			}
			if (ntoken < 4 + npar) {
				ac_err(ac, "invalid geometry data #%d (shape %d needs %d params)",
				       ac->ngeom + 1, shape, npar);
				fclose(fp);
				return 1;
			}
			g = (ac_geom_t *)realloc(ac->geom, (ac->ngeom + 1) * sizeof(ac_geom_t));
			if (!g) { ac_err(ac, "out of memory (geometry)"); fclose(fp); return 1; }
			ac->geom = g;
			memset(&ac->geom[ac->ngeom], 0, sizeof(ac_geom_t));
			ac->geom[ac->ngeom].shape = shape;
			for (n = 0; n < npar; n++)
				ac->geom[ac->ngeom].g[n] = atof(tok[4 + n]);
			ac->ngeom++;
		}
		else if (!strcmp(tok[0], "feed")) {
			if (ntoken < 6) {
				ac_err(ac, "invalid feed data #%d", ac->nfeed + 1);
				fclose(fp);
				return 1;
			}
			if (ac->nfeed == 0) {
				ac->srcx = atof(tok[3]);
				ac->srcy = atof(tok[4]);
				ac->srcz = atof(tok[5]);
			}
			ac->nfeed++;
		}
		else if (!strcmp(tok[0], "point")) {
			ac_recv_t *r;
			if (ntoken < 6) {
				ac_err(ac, "invalid point data #%d", ac->nrecv + 1);
				fclose(fp);
				return 1;
			}
			r = (ac_recv_t *)realloc(ac->recv, (ac->nrecv + 1) * sizeof(ac_recv_t));
			if (!r) { ac_err(ac, "out of memory (point)"); fclose(fp); return 1; }
			ac->recv = r;
			memset(&ac->recv[ac->nrecv], 0, sizeof(ac_recv_t));
			ac->recv[ac->nrecv].x = atof(tok[3]);
			ac->recv[ac->nrecv].y = atof(tok[4]);
			ac->recv[ac->nrecv].z = atof(tok[5]);
			take_trailing_name(save, ac->recv[ac->nrecv].name,
			                   sizeof(ac->recv[ac->nrecv].name));
			ac->nrecv++;
		}
		else if (!strcmp(tok[0], "frequency1") || !strcmp(tok[0], "frequency2")) {
			/* 電磁解析用の周波数リスト。物理には使わないが記録する */
			int n;
			if (!strcmp(tok[0], "frequency1")) {
				double *f = (double *)realloc(ac->freq1,
				            (ac->nfreq1 + ntoken - 2) * sizeof(double));
				if (!f) { ac_err(ac, "out of memory (frequency1)"); fclose(fp); return 1; }
				ac->freq1 = f;
				for (n = 2; n < ntoken; n++)
					ac->freq1[ac->nfreq1++] = atof(tok[n]);
			}
		}
		/* name / material / abc / pml / その他 : 音響では未使用 — 無視 (前方互換) */
	}
	fclose(fp);

	if (nline == 0) {
		ac_err(ac, "empty input '%s'", ac->ofd_path);
		return 1;
	}
	if (!have_mesh[0] || !have_mesh[1] || !have_mesh[2]) {
		ac_err(ac, "no %smesh data in '%s'",
		       !have_mesh[0] ? "x" : (!have_mesh[1] ? "y" : "z"), ac->ofd_path);
		return 1;
	}
	if (ac->nfeed == 0) {
		ac_err(ac, "no feed (acoustic source position) in '%s'", ac->ofd_path);
		return 1;
	}
	if (ac->nfeed > 1)
		ac_log(ac, "warning: %d feeds found — using feed #1 only "
		       "(multi-source is not supported in v1)", ac->nfeed);
	if (ac->nrecv == 0) {
		ac_err(ac, "no point (receiver position) in '%s' — nothing to record",
		       ac->ofd_path);
		return 1;
	}
	if (ac->nonuniform)
		ac_log(ac, "warning: non-uniform mesh — using uniform grid with "
		       "dx = %.6g m (the minimum spacing)", ac->dxmin);

	ac_log(ac, "title: %s", ac->title[0] ? ac->title : "(none)");
	ac_log(ac, "domain: [%.6g, %.6g] x [%.6g, %.6g] x [%.6g, %.6g] m",
	       ac->x0, ac->x1, ac->y0, ac->y1, ac->z0, ac->z1);
	ac_log(ac, "mesh: min spacing %.6g m, %d geometries, %d feed(s), %d point(s)",
	       ac->dxmin, ac->ngeom, ac->nfeed, ac->nrecv);
	if (ac->nfreq1 > 0)
		ac_log(ac, "frequency1: %d values (unused by the acoustic solver)",
		       ac->nfreq1);
	return 0;
}
