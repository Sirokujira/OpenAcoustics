/* ga_ofdx.c — .ofdx (JSON サイドカー) の読み取り (幾何音響ソルバー版)
 *
 * 必要キーのみを読む自前の tolerant scanner (外部 JSON ライブラリ禁止)。
 * 読むのは root.acoustic の下の 3 つだけで、未知キー・未知構造はすべて
 * 構文的に読み飛ばす (前方互換) :
 *
 *   absorption[] = [ { enabled, role, name, area, alpha[6], air_a }, ... ]
 *       role の対応 (OpenFDTD-X の AcousticTab 吸音表) :
 *         4 = Floor    → z-        1 = Ceiling → z+
 *         2 = SideWall → y- と y+   3 = RearWall → x- と x+
 *       FDTD 側は alpha[6] を帯域平均に潰すが、幾何音響は**バンド別に使う**
 *       (これが幾何音響側の利点)。alpha[] は 125/250/500/1k/2k/4k Hz の
 *       6 バンドなので、8 kHz バンドは 4 kHz の値を外挿する (log に明示)。
 *
 *   receivers[] = [ { enabled, name, x, y, z, ... }, ... ]
 *       GUI が持つ受音点名を **座標一致** (許容 = .ofd メッシュ最小刻み) で
 *       .ofd の point 行に引き当てる (FDTD 側と同じ規則)。
 *
 *   ga = { image_order, rays, temperature_c, humidity_percent, pressure_kpa,
 *          air_absorption, receiver_radius_m, duration_s }
 *       幾何音響ソルバー固有の設定。すべて省略可 (省略時は既定値)。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ga.h"

/* ── 最小 JSON スキャナ (カーソル前進型) ───────────────────────── */

typedef struct {
	const char *s;
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

static int jstring(js_t *j, char *out, size_t cap)
{
	size_t n = 0;
	if (jexpect(j, '"')) return 1;
	while (j->s < j->end && *j->s != '"') {
		char c = *j->s++;
		if (c == '\\' && j->s < j->end) {
			c = *j->s++;
			if (c == 'u') {
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
	return 1;
}

/* ── absorption 1 行 ───────────────────────────────────────────── */

typedef struct {
	int    enabled;
	int    role;
	double a[GA_NBAND];
	int    acount;      /* 読めた alpha の個数 */
} arow_t;

static int parse_alpha_array(js_t *j, arow_t *row)
{
	if (jexpect(j, '[')) return 1;
	if (jpeek(j) == ']') { j->s++; return 0; }
	for (;;) {
		double v;
		if (jnumber(j, &v)) return 1;
		if (v < 0.0) v = 0.0;
		if (v > 1.0) v = 1.0;
		if (row->acount < GA_NBAND) row->a[row->acount] = v;
		row->acount++;
		if (jpeek(j) == ',') { j->s++; continue; }
		return jexpect(j, ']');
	}
}

static int parse_arow(js_t *j, arow_t *row)
{
	char key[64];
	int b;
	row->enabled = 1;      /* 欠落キーは既定値 (旧ファイル互換) */
	row->role = 0;
	row->acount = 0;
	for (b = 0; b < GA_NBAND; b++) row->a[b] = GA_ALPHA_DEFAULT;
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

static int role_walls(int role, int walls[2])
{
	switch (role) {
		case 4: walls[0] = GA_ZM; return 1;                  /* Floor */
		case 1: walls[0] = GA_ZP; return 1;                  /* Ceiling */
		case 2: walls[0] = GA_YM; walls[1] = GA_YP; return 2; /* SideWall */
		case 3: walls[0] = GA_XM; walls[1] = GA_XP; return 2; /* RearWall */
		default: return 0;
	}
}

static int parse_absorption(ga_t *g, js_t *j)
{
	int applied[5] = { 0, 0, 0, 0, 0 };   /* role 1..4 の採用済みフラグ */
	if (jexpect(j, '[')) return 1;
	if (jpeek(j) == ']') { j->s++; return 0; }
	for (;;) {
		arow_t row;
		if (parse_arow(j, &row)) return 1;
		if (row.enabled && row.role >= 1 && row.role <= 4 && !applied[row.role]) {
			int walls[2], nw, w, b;
			double a[GA_NBAND];
			for (b = 0; b < GA_NBAND; b++) {
				if (row.acount <= 0)              a[b] = GA_ALPHA_DEFAULT;
				else if (b < row.acount)          a[b] = row.a[b];
				else                              a[b] = row.a[row.acount - 1];
			}
			nw = role_walls(row.role, walls);
			for (w = 0; w < nw; w++)
				for (b = 0; b < GA_NBAND; b++) g->alpha[walls[w]][b] = a[b];
			applied[row.role] = 1;
			g->have_band_alpha = 1;
			if (row.acount > 0 && row.acount < GA_NBAND)
				ga_log(g, ".ofdx: absorption role %d has %d bands — band %d..%d "
				       "(%.0f Hz and above) extrapolated from the last value %.4g",
				       row.role, row.acount, row.acount + 1, GA_NBAND,
				       GA_BAND_F0 * (1 << row.acount), row.a[row.acount - 1]);
			ga_log(g, ".ofdx: absorption role %d -> alpha = "
			       "[%.4g %.4g %.4g %.4g %.4g %.4g %.4g] (125 Hz .. 8 kHz)",
			       row.role, a[0], a[1], a[2], a[3], a[4], a[5], a[6]);
		}
		if (jpeek(j) == ',') { j->s++; continue; }
		return jexpect(j, ']');
	}
}

/* ── receivers ─────────────────────────────────────────────────── */

static void assign_recv_name(ga_t *g, double x, double y, double z,
                             const char *name)
{
	double tol = (g->dxmin > 0.0) ? g->dxmin : 1e-6;
	double best = -1.0;
	int bi = -1, i;

	if (name[0] == '\0') return;
	for (i = 0; i < g->nrecv; i++) {
		double dx = g->recv[i].x - x;
		double dy = g->recv[i].y - y;
		double dz = g->recv[i].z - z;
		double d2 = dx * dx + dy * dy + dz * dz;
		if (g->recv[i].name[0] != '\0') continue;   /* .ofd の名前を優先 */
		if (d2 > tol * tol) continue;
		if (bi < 0 || d2 < best) { best = d2; bi = i; }
	}
	if (bi < 0) return;
	snprintf(g->recv[bi].name, sizeof(g->recv[bi].name), "%s", name);
	ga_log(g, ".ofdx: receiver #%d (%g %g %g) -> name \"%s\"",
	       bi + 1, g->recv[bi].x, g->recv[bi].y, g->recv[bi].z, name);
}

static int parse_receiver_row(ga_t *g, js_t *j)
{
	char key[64];
	char name[GA_NAME_MAX];
	double x = 0.0, y = 0.0, z = 0.0, v;
	int enabled = 1;

	name[0] = '\0';
	if (jexpect(j, '{')) return 1;
	if (jpeek(j) == '}') { j->s++; return 0; }
	for (;;) {
		if (jstring(j, key, sizeof(key))) return 1;
		if (jexpect(j, ':')) return 1;
		if (!strcmp(key, "enabled")) {
			if (jbool(j, &enabled)) return 1;
		}
		else if (!strcmp(key, "x")) { if (jnumber(j, &v)) return 1; x = v; }
		else if (!strcmp(key, "y")) { if (jnumber(j, &v)) return 1; y = v; }
		else if (!strcmp(key, "z")) { if (jnumber(j, &v)) return 1; z = v; }
		else if (!strcmp(key, "name")) {
			if (jstring(j, name, sizeof(name))) return 1;
		}
		else {
			if (jskip(j)) return 1;
		}
		if (jpeek(j) == ',') { j->s++; continue; }
		if (jexpect(j, '}')) return 1;
		if (enabled) assign_recv_name(g, x, y, z, name);
		return 0;
	}
}

static int parse_receivers(ga_t *g, js_t *j)
{
	if (jexpect(j, '[')) return 1;
	if (jpeek(j) == ']') { j->s++; return 0; }
	for (;;) {
		if (parse_receiver_row(g, j)) return 1;
		if (jpeek(j) == ',') { j->s++; continue; }
		return jexpect(j, ']');
	}
}

/* ── acoustic.ga (幾何音響ソルバー固有の設定) ───────────────────── */

static int parse_ga(ga_t *g, js_t *j)
{
	char key[64];
	if (jexpect(j, '{')) return 1;
	if (jpeek(j) == '}') { j->s++; return 0; }
	for (;;) {
		double v;
		int bv;
		if (jstring(j, key, sizeof(key))) return 1;
		if (jexpect(j, ':')) return 1;
		if (!strcmp(key, "image_order")) {
			if (jnumber(j, &v)) return 1;
			g->order = (int)v;
			if (g->order < GA_ORDER_MIN || g->order > GA_ORDER_MAX) {
				ga_err(g, ".ofdx acoustic.ga.image_order = %d is out of range "
				       "(%d..%d)", g->order, GA_ORDER_MIN, GA_ORDER_MAX);
				return 1;
			}
		}
		else if (!strcmp(key, "rays")) {
			if (jnumber(j, &v)) return 1;
			g->nrays = (int)v;
			if (g->nrays < GA_RAYS_MIN || g->nrays > GA_RAYS_MAX) {
				ga_err(g, ".ofdx acoustic.ga.rays = %d is out of range (%d..%d)",
				       g->nrays, GA_RAYS_MIN, GA_RAYS_MAX);
				return 1;
			}
		}
		else if (!strcmp(key, "scattering")) {
			if (jnumber(j, &v)) return 1;
			if (v < 0.0 || v > 1.0) {
				ga_err(g, ".ofdx acoustic.ga.scattering = %g is out of range "
				       "(0 = specular only .. 1 = fully diffuse)", v);
				return 1;
			}
			g->scatter = v;
		}
		else if (!strcmp(key, "temperature_c")) {
			if (jnumber(j, &v)) return 1;
			if (v < -20.0 || v > 50.0) {
				ga_err(g, ".ofdx acoustic.ga.temperature_c = %g is out of range "
				       "(-20 .. 50 degC — ISO 9613-1 の適用範囲)", v);
				return 1;
			}
			g->temp_c = v;
		}
		else if (!strcmp(key, "humidity_percent")) {
			if (jnumber(j, &v)) return 1;
			if (v < 0.0 || v > 100.0) {
				ga_err(g, ".ofdx acoustic.ga.humidity_percent = %g is out of "
				       "range (0 .. 100)", v);
				return 1;
			}
			g->humid = v;
		}
		else if (!strcmp(key, "pressure_kpa")) {
			if (jnumber(j, &v)) return 1;
			if (v < 50.0 || v > 120.0) {
				ga_err(g, ".ofdx acoustic.ga.pressure_kpa = %g is out of range "
				       "(50 .. 120 kPa)", v);
				return 1;
			}
			g->press_kpa = v;
		}
		else if (!strcmp(key, "air_absorption")) {
			if (jbool(j, &bv)) return 1;
			g->air_on = bv;
		}
		else if (!strcmp(key, "receiver_radius_m")) {
			if (jnumber(j, &v)) return 1;
			if (v <= 0.0 || v > 10.0) {
				ga_err(g, ".ofdx acoustic.ga.receiver_radius_m = %g is out of "
				       "range (0 < r <= 10 m)", v);
				return 1;
			}
			g->rsphere_user = v;
		}
		else if (!strcmp(key, "duration_s")) {
			if (jnumber(j, &v)) return 1;
			if (v <= 0.0 || v > 10.0) {
				ga_err(g, ".ofdx acoustic.ga.duration_s = %g is out of range "
				       "(0 < T <= 10 s)", v);
				return 1;
			}
			g->duration_user = v;
		}
		else {
			if (jskip(j)) return 1;   /* 未知キー (前方互換) */
		}
		if (jpeek(j) == ',') { j->s++; continue; }
		return jexpect(j, '}');
	}
}

/* root.acoustic の中身 */
static int walk_acoustic(ga_t *g, js_t *j)
{
	char key[64];
	if (jexpect(j, '{')) return 1;
	if (jpeek(j) == '}') { j->s++; return 0; }
	for (;;) {
		if (jstring(j, key, sizeof(key))) return 1;
		if (jexpect(j, ':')) return 1;
		if (!strcmp(key, "absorption")) {
			if (parse_absorption(g, j)) return 1;
		}
		else if (!strcmp(key, "receivers")) {
			if (parse_receivers(g, j)) return 1;
		}
		else if (!strcmp(key, "ga")) {
			if (parse_ga(g, j)) return 1;
		}
		else {
			if (jskip(j)) return 1;
		}
		if (jpeek(j) == ',') { j->s++; continue; }
		return jexpect(j, '}');
	}
}

static int walk_object(ga_t *g, js_t *j, const char *target)
{
	char key[64];
	if (jexpect(j, '{')) return 1;
	if (jpeek(j) == '}') { j->s++; return 0; }
	for (;;) {
		if (jstring(j, key, sizeof(key))) return 1;
		if (jexpect(j, ':')) return 1;
		if (!strcmp(key, target)) {
			if (walk_acoustic(g, j)) return 1;
		}
		else {
			if (jskip(j)) return 1;
		}
		if (jpeek(j) == ',') { j->s++; continue; }
		return jexpect(j, '}');
	}
}

int ga_read_ofdx(ga_t *g)
{
	FILE *fp = ac_fopen(g->ofdx_path, "rb");
	char *buf;
	long len;
	js_t j;
	int rc, w;

	if (!fp) {
		ga_log(g, "no .ofdx sidecar (%s) — all walls use default alpha = %g "
		       "in every band", g->ofdx_path, GA_ALPHA_DEFAULT);
		return 0;
	}
	g->have_ofdx = 1;
	if (fseek(fp, 0, SEEK_END) != 0 || (len = ftell(fp)) < 0 ||
	    fseek(fp, 0, SEEK_SET) != 0) {
		ga_err(g, "cannot read .ofdx '%s'", g->ofdx_path);
		fclose(fp);
		return 1;
	}
	buf = (char *)malloc((size_t)len + 1);
	if (!buf) {
		ga_err(g, "out of memory (.ofdx, %ld bytes)", len);
		fclose(fp);
		return 1;
	}
	if (fread(buf, 1, (size_t)len, fp) != (size_t)len) {
		ga_err(g, "cannot read .ofdx '%s'", g->ofdx_path);
		free(buf);
		fclose(fp);
		return 1;
	}
	fclose(fp);
	buf[len] = '\0';

	j.s = buf;
	j.end = buf + len;
	if (len >= 3 && (unsigned char)buf[0] == 0xEF &&
	    (unsigned char)buf[1] == 0xBB && (unsigned char)buf[2] == 0xBF)
		j.s += 3;

	ga_log(g, "reading .ofdx sidecar: %s", g->ofdx_path);
	rc = walk_object(g, &j, "acoustic");
	free(buf);
	if (rc != 0) {
		/* 壊れた JSON や値域外の設定で黙って既定値を使うのは「静かに狂う」ので
		 * 明確に失敗させる。値域エラーはどのキーが悪いかを先に ga_err 済み。 */
		ga_err(g, "malformed JSON or bad value in .ofdx '%s' "
		       "(fix or remove the sidecar)", g->ofdx_path);
		return 1;
	}
	for (w = 0; w < GA_NWALL; w++)
		ga_log(g, "wall %s: alpha = [%.4g %.4g %.4g %.4g %.4g %.4g %.4g]",
		       (w == GA_XM) ? "x-" : (w == GA_XP) ? "x+" :
		       (w == GA_YM) ? "y-" : (w == GA_YP) ? "y+" :
		       (w == GA_ZM) ? "z-" : "z+",
		       g->alpha[w][0], g->alpha[w][1], g->alpha[w][2], g->alpha[w][3],
		       g->alpha[w][4], g->alpha[w][5], g->alpha[w][6]);
	return 0;
}
