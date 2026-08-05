/* input_ofdx.c — .ofdx (JSON サイドカー) の acoustic.absorption 読み
 *
 * 必要キーのみを読む自前の tolerant scanner (外部 JSON ライブラリ禁止)。
 * 対象は OpenFDTD-X/src/io/OfdIO.cpp が書く
 *   root.acoustic.absorption = [ { enabled, role, name, area, alpha[6], air_a }, ... ]
 * のみで、未知キー・未知構造はすべて構文的に読み飛ばす (前方互換)。
 *
 * role の対応 (指示書 / AcousticTab の吸音表) :
 *   4 = Floor    → z-        1 = Ceiling → z+
 *   2 = SideWall → y- と y+   3 = RearWall → x- と x+
 * 他の role と欠落した壁は既定 α = AC_ALPHA_DEFAULT (0.1)。
 * 6 バンド alpha 配列は v1 では帯域平均を使う (周波数依存境界は将来課題)。
 * 同一 role の行が複数あるときは最初の enabled 行を使う。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "acoustic.h"

/* ── 最小 JSON スキャナ (カーソル前進型) ───────────────────────── */

typedef struct {
	const char *s;    /* 現在位置 */
	const char *end;
} js_t;

static void jws(js_t *j)
{
	while (j->s < j->end && (*j->s == ' ' || *j->s == '\t' ||
	                         *j->s == '\r' || *j->s == '\n'))
		j->s++;
}

static int jpeek(js_t *j)
{
	jws(j);
	return (j->s < j->end) ? (unsigned char)*j->s : -1;
}

static int jexpect(js_t *j, char c)
{
	if (jpeek(j) != c) return 1;
	j->s++;
	return 0;
}

/* 文字列 : 開始の '"' から終端まで進め、out (cap>0 なら) へコピー */
static int jstring(js_t *j, char *out, size_t cap)
{
	size_t n = 0;
	if (jexpect(j, '"')) return 1;
	while (j->s < j->end && *j->s != '"') {
		char c = *j->s++;
		if (c == '\\' && j->s < j->end) {
			c = *j->s++;   /* エスケープは素通し (キー比較にのみ使う) */
			if (c == 'u') {  /* \uXXXX */
				int k;
				for (k = 0; k < 4 && j->s < j->end; k++) j->s++;
				c = '?';
			}
		}
		if (out && n + 1 < cap) out[n++] = c;
	}
	if (out && cap > 0) out[n] = '\0';
	return jexpect(j, '"');
}

static int jnumber(js_t *j, double *v)
{
	char *ep = NULL;
	double d;
	jws(j);
	d = strtod(j->s, &ep);
	if (ep == j->s || ep > j->end) return 1;
	if (v) *v = d;
	j->s = ep;
	return 0;
}

/* 任意の値 1 個を構文的に読み飛ばす (再帰) */
static int jskip(js_t *j)
{
	int c = jpeek(j);
	if (c < 0) return 1;
	if (c == '"') return jstring(j, NULL, 0);
	if (c == '{' || c == '[') {
		char close = (c == '{') ? '}' : ']';
		j->s++;
		if (jpeek(j) == close) { j->s++; return 0; }
		for (;;) {
			if (c == '{') {
				if (jstring(j, NULL, 0)) return 1;
				if (jexpect(j, ':')) return 1;
			}
			if (jskip(j)) return 1;
			if (jpeek(j) == ',') { j->s++; continue; }
			return jexpect(j, close);
		}
	}
	if (c == 't') { if (j->end - j->s >= 4 && !strncmp(j->s, "true", 4))  { j->s += 4; return 0; } return 1; }
	if (c == 'f') { if (j->end - j->s >= 5 && !strncmp(j->s, "false", 5)) { j->s += 5; return 0; } return 1; }
	if (c == 'n') { if (j->end - j->s >= 4 && !strncmp(j->s, "null", 4))  { j->s += 4; return 0; } return 1; }
	return jnumber(j, NULL);
}

static int jbool(js_t *j, int *v)
{
	int c = jpeek(j);
	if (c == 't') { *v = 1; return jskip(j); }
	if (c == 'f') { *v = 0; return jskip(j); }
	return 1;   /* true/false 以外 */
}

/* ── absorption 1 行の解釈 ─────────────────────────────────────── */

typedef struct {
	int    enabled;
	int    role;
	double asum;
	int    acount;
} row_t;

static int parse_alpha_array(js_t *j, row_t *row)
{
	if (jexpect(j, '[')) return 1;
	if (jpeek(j) == ']') { j->s++; return 0; }
	for (;;) {
		double v;
		if (jnumber(j, &v)) return 1;
		row->asum += v;
		row->acount++;
		if (jpeek(j) == ',') { j->s++; continue; }
		return jexpect(j, ']');
	}
}

static int parse_row(js_t *j, row_t *row)
{
	char key[64];
	row->enabled = 1;      /* 欠落キーは既定値 (旧ファイル互換) */
	row->role = 0;
	row->asum = 0.0;
	row->acount = 0;
	if (jexpect(j, '{')) return 1;
	if (jpeek(j) == '}') { j->s++; return 0; }
	for (;;) {
		double v;
		if (jstring(j, key, sizeof(key))) return 1;
		if (jexpect(j, ':')) return 1;
		if (!strcmp(key, "enabled")) {
			if (jbool(j, &row->enabled)) return 1;
		}
		else if (!strcmp(key, "role")) {
			if (jnumber(j, &v)) return 1;
			row->role = (int)v;
		}
		else if (!strcmp(key, "alpha")) {
			if (parse_alpha_array(j, row)) return 1;
		}
		else {
			if (jskip(j)) return 1;   /* name / area / air_a / 未知キー */
		}
		if (jpeek(j) == ',') { j->s++; continue; }
		return jexpect(j, '}');
	}
}

/* role → 壁インデックス列。戻り値 = 壁数 */
static int role_walls(int role, int walls[2])
{
	switch (role) {
		case 4: walls[0] = AC_ZM; return 1;                  /* Floor */
		case 1: walls[0] = AC_ZP; return 1;                  /* Ceiling */
		case 2: walls[0] = AC_YM; walls[1] = AC_YP; return 2; /* SideWall */
		case 3: walls[0] = AC_XM; walls[1] = AC_XP; return 2; /* RearWall */
		default: return 0;
	}
}

static int parse_absorption(ac_t *ac, js_t *j)
{
	int applied[5] = { 0, 0, 0, 0, 0 };   /* role 1..4 の採用済みフラグ */
	if (jexpect(j, '[')) return 1;
	if (jpeek(j) == ']') { j->s++; return 0; }
	for (;;) {
		row_t row;
		if (parse_row(j, &row)) return 1;
		if (row.enabled && row.role >= 1 && row.role <= 4 && !applied[row.role]) {
			int walls[2], nw, w;
			double a = (row.acount > 0) ? row.asum / row.acount : AC_ALPHA_DEFAULT;
			if (a < 0.0) a = 0.0;
			if (a > 1.0) a = 1.0;
			nw = role_walls(row.role, walls);
			for (w = 0; w < nw; w++) ac->alpha[walls[w]] = a;
			applied[row.role] = 1;
			ac_log(ac, ".ofdx: absorption role %d -> alpha = %.4g (band average of %d)",
			       row.role, a, row.acount);
		}
		if (jpeek(j) == ',') { j->s++; continue; }
		return jexpect(j, ']');
	}
}

/* オブジェクトを歩き、目的キーだけ target で処理し他は読み飛ばす */
static int walk_object(ac_t *ac, js_t *j, const char *target, int depth)
{
	char key[64];
	if (jexpect(j, '{')) return 1;
	if (jpeek(j) == '}') { j->s++; return 0; }
	for (;;) {
		if (jstring(j, key, sizeof(key))) return 1;
		if (jexpect(j, ':')) return 1;
		if (!strcmp(key, target)) {
			if (depth == 0) {
				/* root.acoustic → その中の absorption へ */
				if (walk_object(ac, j, "absorption", 1)) return 1;
			}
			else {
				if (parse_absorption(ac, j)) return 1;
			}
		}
		else {
			if (jskip(j)) return 1;
		}
		if (jpeek(j) == ',') { j->s++; continue; }
		return jexpect(j, '}');
	}
}

int ac_read_ofdx(ac_t *ac)
{
	FILE *fp = fopen(ac->ofdx_path, "rb");
	char *buf;
	long len;
	js_t j;
	int rc;

	if (!fp) {
		ac_log(ac, "no .ofdx sidecar (%s) — all walls use default alpha = %g",
		       ac->ofdx_path, AC_ALPHA_DEFAULT);
		return 0;
	}
	ac->have_ofdx = 1;
	if (fseek(fp, 0, SEEK_END) != 0 || (len = ftell(fp)) < 0 ||
	    fseek(fp, 0, SEEK_SET) != 0) {
		ac_err(ac, "cannot read .ofdx '%s'", ac->ofdx_path);
		fclose(fp);
		return 1;
	}
	buf = (char *)malloc((size_t)len + 1);
	if (!buf) {
		ac_err(ac, "out of memory (.ofdx, %ld bytes)", len);
		fclose(fp);
		return 1;
	}
	if (fread(buf, 1, (size_t)len, fp) != (size_t)len) {
		ac_err(ac, "cannot read .ofdx '%s'", ac->ofdx_path);
		free(buf);
		fclose(fp);
		return 1;
	}
	fclose(fp);
	buf[len] = '\0';

	j.s = buf;
	j.end = buf + len;
	/* UTF-8 BOM は読み飛ばす */
	if (len >= 3 && (unsigned char)buf[0] == 0xEF &&
	    (unsigned char)buf[1] == 0xBB && (unsigned char)buf[2] == 0xBF)
		j.s += 3;

	ac_log(ac, "reading .ofdx sidecar: %s", ac->ofdx_path);
	rc = walk_object(ac, &j, "acoustic", 0);
	free(buf);
	if (rc != 0) {
		/* 壊れた JSON で黙って既定値を使うのは「静かに狂う」ので明確に失敗させる */
		ac_err(ac, "malformed JSON in .ofdx '%s' (fix or remove the sidecar)",
		       ac->ofdx_path);
		return 1;
	}
	ac_log(ac, "wall absorption (band-averaged): "
	       "x- %.4g, x+ %.4g, y- %.4g, y+ %.4g, z- %.4g, z+ %.4g",
	       ac->alpha[AC_XM], ac->alpha[AC_XP], ac->alpha[AC_YM],
	       ac->alpha[AC_YP], ac->alpha[AC_ZM], ac->alpha[AC_ZP]);
	return 0;
}
