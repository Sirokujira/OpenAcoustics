/* ga_input.c — OpenFDTD 互換 .ofd パーサ (幾何音響ソルバーが使う部分集合)
 *
 * 書式は本家 OpenFDTD/sol/input_data.c と OpenFDTD-X/src/io/OfdIO.cpp が正。
 * 読むキー : xmesh/ymesh/zmesh, geometry, feed, point, title。
 * その他の未知キーは本家と同様に無視する (前方互換)。
 *
 * 幾何音響での解釈 :
 *   - メッシュ : 「室」の直方体 [x0,x1]x[y0,y1]x[z0,z1] を決めるためだけに読む。
 *                幾何音響に格子は無いので刻みは使わない (最小刻みは .ofdx の
 *                受音点名を座標一致で引き当てる際の許容値としてのみ使う)。
 *   - geometry : 剛体障害物。可視性判定 (遮蔽) と光線追跡の反射に使う。
 *   - feed     : 音源位置 (#1 のみ)。
 *   - point    : 受音点。行末の "# 名前" を WAV ファイル名に使う。
 *
 * ofdx_acoustic_fdtd の src/input_ofd.c と同じ書式を読むが、別バイナリとして
 * 独立させる方針 (FDTD 側のソースには手を入れない) のため実装は共有しない。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ga.h"

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

/* working_dir 直下の唯一の .ofd を探す。0 個・複数個は曖昧なので非零終了 —
 * 数値を捏造した合成 RIR は出さない。 */
int ga_find_input(ga_t *g)
{
	char found[GA_NAME_MAX * 4] = "";
	int count = 0;
#ifdef _WIN32
	WIN32_FIND_DATAA fd;
	HANDLE h;
	char pat[GA_PATH_MAX + 8];
	snprintf(pat, sizeof(pat), "%s\\*", g->workdir);
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
	DIR *d = opendir(g->workdir);
	struct dirent *e;
	if (!d) {
		ga_err(g, "cannot open working dir '%s'", g->workdir);
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
		ga_err(g, "no input given and no .ofd file found in '%s'", g->workdir);
		return 1;
	}
	if (count > 1) {
		ga_err(g, "no input given and %d .ofd files found in '%s' "
		       "(ambiguous — pass the input file name explicitly)",
		       count, g->workdir);
		return 1;
	}
	if (snprintf(g->ofd_path, sizeof(g->ofd_path), "%s/%s",
	             g->workdir, found) >= (int)sizeof(g->ofd_path)) {
		ga_err(g, "input path too long");
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

/* mesh 1 軸分 : "x0 d1 x1 d2 x2 ..." → 範囲と最小刻み */
static int parse_mesh(ga_t *g, char **tok, int ntoken,
                      double *lo, double *hi, const char *key)
{
	double prev, next, d;
	int i, ndiv, div;
	if (ntoken < 5 || ntoken % 2 == 0) {
		ga_err(g, "invalid %s data", key);
		return 1;
	}
	ndiv = (ntoken - 3) / 2;
	prev = atof(tok[2]);
	*lo = prev;
	for (i = 0; i < ndiv; i++) {
		div  = atoi(tok[2 * i + 3]);
		next = atof(tok[2 * i + 4]);
		if (div <= 0 || next <= prev) {
			ga_err(g, "invalid %s data (non-increasing nodes or bad division)", key);
			return 1;
		}
		d = (next - prev) / div;
		if (g->dxmin <= 0.0 || d < g->dxmin) g->dxmin = d;
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

int ga_read_ofd(ga_t *g)
{
	char line[MAXLINE], save[MAXLINE];
	char *tok[MAXTOKEN];
	int ntoken, nline = 0;
	int have_mesh[3] = { 0, 0, 0 };
	FILE *fp = ac_fopen(g->ofd_path, "r");

	if (!fp) {
		ga_err(g, "cannot open input '%s'", g->ofd_path);
		return 1;
	}

	while (fgets(line, sizeof(line), fp) != NULL) {
		if (!strncmp(line, "end", 3)) break;
		snprintf(save, sizeof(save), "%s", line);
		ntoken = tokenize(line, tok, MAXTOKEN);
		if (ntoken == 0) continue;

		if (nline == 0) {
			if (strcmp(tok[0], "OpenFDTD") && strcmp(tok[0], "OpenTHFD") &&
			    strcmp(tok[0], "OpenRCWA")) {
				ga_err(g, "not OpenFDTD/OpenTHFD data: %s", g->ofd_path);
				fclose(fp);
				return 1;
			}
			nline++;
			continue;
		}
		if (ntoken < 3 || strcmp(tok[1], "=")) continue;

		if (!strcmp(tok[0], "title")) {
			const char *eq = strchr(save, '=');
			size_t n;
			if (!eq) continue;
			eq++;
			while (*eq == ' ') eq++;
			n = strlen(eq);
			while (n > 0 && (eq[n-1] == '\r' || eq[n-1] == '\n')) n--;
			if (n >= sizeof(g->title)) n = sizeof(g->title) - 1;
			memcpy(g->title, eq, n);
			g->title[n] = '\0';
		}
		else if (!strcmp(tok[0], "xmesh")) {
			if (parse_mesh(g, tok, ntoken, &g->x0, &g->x1, "xmesh")) { fclose(fp); return 1; }
			have_mesh[0] = 1;
		}
		else if (!strcmp(tok[0], "ymesh")) {
			if (parse_mesh(g, tok, ntoken, &g->y0, &g->y1, "ymesh")) { fclose(fp); return 1; }
			have_mesh[1] = 1;
		}
		else if (!strcmp(tok[0], "zmesh")) {
			if (parse_mesh(g, tok, ntoken, &g->z0, &g->z1, "zmesh")) { fclose(fp); return 1; }
			have_mesh[2] = 1;
		}
		else if (!strcmp(tok[0], "geometry")) {
			int shape, npar, n;
			ga_geom_t *gg;
			if (ntoken < 4) {
				ga_err(g, "invalid geometry data #%d", g->ngeom + 1);
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
				ga_err(g, "invalid geometry data #%d (shape %d needs %d params)",
				       g->ngeom + 1, shape, npar);
				fclose(fp);
				return 1;
			}
			gg = (ga_geom_t *)realloc(g->geom, (g->ngeom + 1) * sizeof(ga_geom_t));
			if (!gg) { ga_err(g, "out of memory (geometry)"); fclose(fp); return 1; }
			g->geom = gg;
			memset(&g->geom[g->ngeom], 0, sizeof(ga_geom_t));
			g->geom[g->ngeom].shape = shape;
			for (n = 0; n < npar; n++)
				g->geom[g->ngeom].g[n] = atof(tok[4 + n]);
			g->ngeom++;
		}
		else if (!strcmp(tok[0], "feed")) {
			if (ntoken < 6) {
				ga_err(g, "invalid feed data #%d", g->nfeed + 1);
				fclose(fp);
				return 1;
			}
			if (g->nfeed == 0) {
				g->srcx = atof(tok[3]);
				g->srcy = atof(tok[4]);
				g->srcz = atof(tok[5]);
			}
			g->nfeed++;
		}
		else if (!strcmp(tok[0], "point")) {
			ga_recv_t *r;
			if (ntoken < 6) {
				ga_err(g, "invalid point data #%d", g->nrecv + 1);
				fclose(fp);
				return 1;
			}
			r = (ga_recv_t *)realloc(g->recv, (g->nrecv + 1) * sizeof(ga_recv_t));
			if (!r) { ga_err(g, "out of memory (point)"); fclose(fp); return 1; }
			g->recv = r;
			memset(&g->recv[g->nrecv], 0, sizeof(ga_recv_t));
			g->recv[g->nrecv].x = atof(tok[3]);
			g->recv[g->nrecv].y = atof(tok[4]);
			g->recv[g->nrecv].z = atof(tok[5]);
			take_trailing_name(save, g->recv[g->nrecv].name,
			                   sizeof(g->recv[g->nrecv].name));
			g->nrecv++;
		}
		/* frequency1/2 / name / material / abc / pml ほか : 幾何音響では未使用 */
	}
	fclose(fp);

	if (nline == 0) {
		ga_err(g, "empty input '%s'", g->ofd_path);
		return 1;
	}
	if (!have_mesh[0] || !have_mesh[1] || !have_mesh[2]) {
		ga_err(g, "no %smesh data in '%s'",
		       !have_mesh[0] ? "x" : (!have_mesh[1] ? "y" : "z"), g->ofd_path);
		return 1;
	}
	if (g->nfeed == 0) {
		ga_err(g, "no feed (acoustic source position) in '%s'", g->ofd_path);
		return 1;
	}
	if (g->nfeed > 1)
		ga_log(g, "warning: %d feeds found — using feed #1 only "
		       "(multi-source is not supported in v1)", g->nfeed);
	if (g->nrecv == 0) {
		ga_err(g, "no point (receiver position) in '%s' — nothing to record",
		       g->ofd_path);
		return 1;
	}

	ga_log(g, "title: %s", g->title[0] ? g->title : "(none)");
	ga_log(g, "room: [%.6g, %.6g] x [%.6g, %.6g] x [%.6g, %.6g] m",
	       g->x0, g->x1, g->y0, g->y1, g->z0, g->z1);
	ga_log(g, "input: %d geometries, %d feed(s), %d point(s)",
	       g->ngeom, g->nfeed, g->nrecv);
	return 0;
}
