/* ga_setup.c — 幾何音響ソルバーの前処理 (室・境界・空気吸収・時間軸・確保)
 *
 * ここで決まる物理量 :
 *   - 室 : .ofd のメッシュ範囲がそのまま直方体の室。格子は作らない
 *          (幾何音響は連続座標で計算する — 受音点も音源も丸めない)。
 *   - 境界 : バンド別の圧力反射係数 R_b = sqrt(1 - alpha_b)。
 *          幾何音響の標準的な「エネルギー反射率 = 1 - alpha」の定義そのもの
 *          (垂直入射・局所反応。斜入射の角度依存は v1 では扱わない)。
 *   - 空気吸収 : ISO 9613-1 の減衰係数 [dB/m] をバンド中心周波数で評価する。
 *   - 計算時間 : T = clamp(1.5 * max_b T_Eyring(b), 0.5, 3.0) [s]。
 *          T_Eyring = 0.161 V / (-S ln(1 - alpha_mean) + 4 m V)
 *          (Eyring–Norris + 空気吸収項。m はエネルギー減衰係数 [1/m])。
 *          直接音が必ず入るよう「最遠受音点までの伝搬時間 + 50 ms」も下限にする。
 *
 * 有効帯域 (metadata の valid_band_hz) :
 *   下限 f_lo = Schroeder 周波数 2000 sqrt(T60 / V) [Hz]。
 *     これより下ではモード密度が足りず、鏡像法・光線追跡の前提 (波長 <<
 *     部屋・壁の寸法、干渉を無視した強度加算) が成り立たない。
 *   上限 f_hi = min(fs/2, 8 kHz オクターブバンドの上端 8000*sqrt(2) Hz)。
 *     吸音率も空気吸収も 8 kHz バンドまでしか定義していないため、
 *     それより上は外挿になる (11.3 kHz 〜 24 kHz は 8 kHz バンドの値を延長)。
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ga.h"

#define GA_MAX_SAMPLES (10 * GA_FS)

/* ── ISO 9613-1 : 大気による純音減衰係数 [dB/m] ──────────────────
 * ISO 9613-1:1993 の式 (3)〜(5)。T は絶対温度 [K]、p_a は気圧、
 * h は水蒸気のモル濃度 [%]。20 ℃・50 %・101.325 kPa で
 * 4 kHz -> 約 0.0328 dB/m (= 32.8 dB/km、規格 Table 1 と一致)。 */
double ga_air_alpha_db_m(double f_hz, double temp_c, double humid_pct,
                         double press_kpa)
{
	const double T0  = 293.15;    /* 基準温度 (20 ℃) [K] */
	const double T01 = 273.16;    /* 水の三重点 [K] */
	const double pr  = 101.325;   /* 基準大気圧 [kPa] */
	double T   = temp_c + 273.15;
	double pa  = press_kpa / pr;            /* 相対気圧 p_a/p_r */
	double Tr  = T / T0;                    /* 相対温度 T/T0 */
	double psat, h, frO, frN, f2;

	if (T <= 0.0 || pa <= 0.0) return 0.0;
	/* 飽和水蒸気圧 p_sat/p_r (ISO 9613-1 式 (B.3)) */
	psat = pow(10.0, -6.8346 * pow(T01 / T, 1.261) + 4.6151);
	h    = humid_pct * psat / pa;           /* モル濃度 [%] */
	/* 酸素 / 窒素の緩和周波数 [Hz] */
	frO  = pa * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h));
	frN  = pa * pow(Tr, -0.5)
	     * (9.0 + 280.0 * h * exp(-4.170 * (pow(Tr, -1.0 / 3.0) - 1.0)));
	f2   = f_hz * f_hz;
	return 8.686 * f2
	     * (1.84e-11 / pa * sqrt(Tr)
	        + pow(Tr, -2.5)
	          * (0.01275 * exp(-2239.1 / T) / (frO + f2 / frO)
	             + 0.1068 * exp(-3352.0 / T) / (frN + f2 / frN)));
}

/* ── 幾何 : AABB (shape 1 は厳密に同じ箱)。戻り値 0 = 対応しない shape ──
 * 並びは本家 OpenFDTD sol/ingeometry.c と同じ (FDTD 側 fdtd.c と同一規則)。 */
static int geom_aabb(const ga_geom_t *g, double lo[3], double hi[3])
{
	const double *p = g->g;
	int i;
	double a, b;
	switch (g->shape) {
		case 1: case 2: case 11: case 12: case 13:
			for (i = 0; i < 3; i++) {
				a = p[2*i]; b = p[2*i + 1];
				lo[i] = (a < b) ? a : b;
				hi[i] = (a < b) ? b : a;
			}
			return 1;
		case 31: case 32: case 33: {
			int ax = g->shape - 31;
			int u = (ax + 1) % 3, w = (ax + 2) % 3;
			lo[ax] = (p[0] < p[1]) ? p[0] : p[1];
			hi[ax] = (p[0] < p[1]) ? p[1] : p[0];
			lo[u] = hi[u] = p[2];
			lo[w] = hi[w] = p[5];
			for (i = 1; i < 3; i++) {
				if (p[2+i] < lo[u]) lo[u] = p[2+i];
				if (p[2+i] > hi[u]) hi[u] = p[2+i];
				if (p[5+i] < lo[w]) lo[w] = p[5+i];
				if (p[5+i] > hi[w]) hi[w] = p[5+i];
			}
			return 1;
		}
		case 41: case 42: case 43:
		case 51: case 52: case 53: {
			int ax = (g->shape % 10) - 1;
			int u = (ax + 1) % 3, w = (ax + 2) % 3;
			double hu = ((p[4] > p[6]) ? p[4] : p[6]) / 2;
			double hw = ((p[5] > p[7]) ? p[5] : p[7]) / 2;
			lo[ax] = (p[0] < p[1]) ? p[0] : p[1];
			hi[ax] = (p[0] < p[1]) ? p[1] : p[0];
			lo[u] = p[2] - hu; hi[u] = p[2] + hu;
			lo[w] = p[3] - hw; hi[w] = p[3] + hw;
			return 1;
		}
		default:
			return 0;
	}
}

/* 点を含む剛体形状の番号 (0 起点。無ければ -1) */
static int geom_containing(const ga_t *g, double x, double y, double z)
{
	int n;
	for (n = 0; n < g->ngeom; n++) {
		const ga_geom_t *o = &g->geom[n];
		if (!o->ok) continue;
		if (x >= o->lo[0] && x <= o->hi[0] && y >= o->lo[1] && y <= o->hi[1] &&
		    z >= o->lo[2] && z <= o->hi[2])
			return n;
	}
	return -1;
}

/* 点から AABB 表面までの距離 (内部なら 0) */
static double geom_distance(const ga_geom_t *o, double x, double y, double z)
{
	double p[3], d = 0.0;
	int i;
	p[0] = x; p[1] = y; p[2] = z;
	for (i = 0; i < 3; i++) {
		double e = 0.0;
		if (p[i] < o->lo[i]) e = o->lo[i] - p[i];
		else if (p[i] > o->hi[i]) e = p[i] - o->hi[i];
		d += e * e;
	}
	return sqrt(d);
}

/* 名前 → ファイル名に使える形 (FDTD 側 fdtd.c と同じ規則) */
static void safe_name(const char *name, char *out, size_t cap)
{
	size_t i;
	for (i = 0; name[i] != '\0' && i + 1 < cap; i++) {
		unsigned char c = (unsigned char)name[i];
		int bad = (c < 0x20) || c == '/' || c == '\\' || c == ':' ||
		          c == '*' || c == '?' || c == '"' || c == '<' ||
		          c == '>' || c == '|' || c == ' ';
		out[i] = bad ? '_' : (char)c;
	}
	out[i] = '\0';
}

static void recv_filename(ga_recv_t *r, int index)
{
	char safe[GA_NAME_MAX];
	r->alias[0] = '\0';
	if (index == 0) {
		strcpy(r->file, "rir.wav");
		if (r->name[0] != '\0') {
			safe_name(r->name, safe, sizeof(safe));
			snprintf(r->alias, sizeof(r->alias), "rir_%s.wav", safe);
		}
		return;
	}
	if (r->name[0] == '\0') {
		snprintf(r->file, sizeof(r->file), "rir_%d.wav", index + 1);
	}
	else {
		safe_name(r->name, safe, sizeof(safe));
		snprintf(r->file, sizeof(r->file), "rir_%s.wav", safe);
	}
}

static int next_pow2(int n)
{
	int m = 1;
	while (m < n) m <<= 1;
	return m;
}

int ga_setup(ga_t *g)
{
	const double c = GA_C0;
	double lx, ly, lz, swall[GA_NWALL], stot, rmax, tsolve;
	double rbase;
	int w, b, n, r;

	/* ── 室 ── */
	lx = g->x1 - g->x0;
	ly = g->y1 - g->y0;
	lz = g->z1 - g->z0;
	if (lx <= 0.0 || ly <= 0.0 || lz <= 0.0) {
		ga_err(g, "degenerate room %.6g x %.6g x %.6g m — check xmesh/ymesh/zmesh",
		       lx, ly, lz);
		return 1;
	}
	g->vol = lx * ly * lz;
	swall[GA_XM] = swall[GA_XP] = ly * lz;
	swall[GA_YM] = swall[GA_YP] = lx * lz;
	swall[GA_ZM] = swall[GA_ZP] = lx * ly;
	stot = 2.0 * (ly * lz + lx * lz + lx * ly);
	g->surf = stot;

	/* ── 剛体障害物 (AABB) ── */
	for (n = 0; n < g->ngeom; n++) {
		ga_geom_t *o = &g->geom[n];
		o->ok = geom_aabb(o, o->lo, o->hi);
		o->exact = (o->shape == 1);
		if (!o->ok) {
			ga_log(g, "warning: geometry #%d has unknown shape %d — ignored",
			       n + 1, o->shape);
			continue;
		}
		if (!o->exact)
			ga_log(g, "warning: geometry #%d shape %d approximated by its AABB "
			       "(only shape 1 = box is exact in v1) — occlusion and ray "
			       "reflection use the box [%.4g %.4g] x [%.4g %.4g] x [%.4g %.4g]",
			       n + 1, o->shape, o->lo[0], o->hi[0], o->lo[1], o->hi[1],
			       o->lo[2], o->hi[2]);
	}
	if (g->ngeom > 0)
		ga_log(g, "obstacles: %d geometry entries treated as rigid (alpha = 0) — "
		       "they occlude image-source paths and reflect rays, but their own "
		       "1st/2nd order reflections are only in the ray (late) part",
		       g->ngeom);

	/* ── 音源・受音点の妥当性 (捏造しない : 室外・剛体内は失敗させる) ── */
	if (g->srcx < g->x0 || g->srcx > g->x1 || g->srcy < g->y0 || g->srcy > g->y1 ||
	    g->srcz < g->z0 || g->srcz > g->z1) {
		ga_err(g, "feed position (%.4g, %.4g, %.4g) is outside the room "
		       "[%.4g,%.4g] x [%.4g,%.4g] x [%.4g,%.4g]",
		       g->srcx, g->srcy, g->srcz, g->x0, g->x1, g->y0, g->y1,
		       g->z0, g->z1);
		return 1;
	}
	n = geom_containing(g, g->srcx, g->srcy, g->srcz);
	if (n >= 0) {
		ga_err(g, "feed position (%.4g, %.4g, %.4g) is inside rigid geometry "
		       "#%d (shape %d) — move the source out of the object "
		       "(e.g. 1.5 m above the stage floor)",
		       g->srcx, g->srcy, g->srcz, n + 1, g->geom[n].shape);
		return 1;
	}
	rmax = 0.0;
	for (r = 0; r < g->nrecv; r++) {
		ga_recv_t *rv = &g->recv[r];
		double dx, dy, dz, d;
		if (rv->x < g->x0 || rv->x > g->x1 || rv->y < g->y0 || rv->y > g->y1 ||
		    rv->z < g->z0 || rv->z > g->z1) {
			ga_err(g, "point #%d (%.4g, %.4g, %.4g) is outside the room "
			       "[%.4g,%.4g] x [%.4g,%.4g] x [%.4g,%.4g]",
			       r + 1, rv->x, rv->y, rv->z, g->x0, g->x1, g->y0, g->y1,
			       g->z0, g->z1);
			return 1;
		}
		n = geom_containing(g, rv->x, rv->y, rv->z);
		if (n >= 0) {
			ga_err(g, "point #%d (%.4g, %.4g, %.4g) is inside rigid geometry "
			       "#%d (shape %d) — move it out of the object, e.g. to ear "
			       "height above the floor/audience block "
			       "(ISO 3382-1: 1.2 m seated)",
			       r + 1, rv->x, rv->y, rv->z, n + 1, g->geom[n].shape);
			return 1;
		}
		dx = rv->x - g->srcx; dy = rv->y - g->srcy; dz = rv->z - g->srcz;
		d = sqrt(dx * dx + dy * dy + dz * dz);
		if (d < 1e-6) {
			ga_err(g, "point #%d coincides with the feed position "
			       "(%.4g, %.4g, %.4g) — the direct sound 1/(4 pi r) diverges",
			       r + 1, rv->x, rv->y, rv->z);
			return 1;
		}
		if (d > rmax) rmax = d;
		recv_filename(rv, r);
	}

	/* ── 境界 : バンド別の圧力反射係数 R = sqrt(1 - alpha) ── */
	for (w = 0; w < GA_NWALL; w++)
		for (b = 0; b < GA_NBAND; b++) {
			double a = g->alpha[w][b];
			if (a < 0.0) a = 0.0;
			if (a > 1.0) a = 1.0;
			g->alpha[w][b] = a;
			g->refl[w][b] = sqrt(1.0 - a);
		}

	/* ── 空気吸収 (ISO 9613-1) ── */
	for (b = 0; b < GA_NBAND; b++) {
		double fc = GA_BAND_F0 * (double)(1 << b);
		g->air_db_m[b] = g->air_on
		               ? ga_air_alpha_db_m(fc, g->temp_c, g->humid, g->press_kpa)
		               : 0.0;
		/* エネルギー減衰 exp(-m s) と音圧レベル減衰 10^(-a s/20) の対応 :
		 * 10^(-a s/10) = exp(-m s) より m = a ln(10)/10 [1/m] */
		g->air_e[b] = g->air_db_m[b] * log(10.0) / 10.0;
	}
	if (g->air_on)
		ga_log(g, "air absorption (ISO 9613-1, %.4g degC, %.4g %%RH, %.6g kPa) "
		       "[dB/m]: %.6g %.6g %.6g %.6g %.6g %.6g %.6g",
		       g->temp_c, g->humid, g->press_kpa,
		       g->air_db_m[0], g->air_db_m[1], g->air_db_m[2], g->air_db_m[3],
		       g->air_db_m[4], g->air_db_m[5], g->air_db_m[6]);
	else
		ga_log(g, "air absorption disabled (.ofdx acoustic.ga.air_absorption = false)");

	/* ── 残響時間の推定 (バンド別)。障害物は算入しない (室のみの推定値で、
	 * 計算時間の決定と metrics の参考値にしか使わない — 実際の減衰は
	 * 光線追跡の結果そのもの) ──
	 *   Eyring : T = 0.161 V / (-S ln(1 - alpha_mean) + 4 m V)
	 *   Sabine : T = 0.161 V / ( S alpha_mean        + 4 m V)  (参考値) */
	g->t60max = 0.0;
	g->tsabmax = 0.0;
	for (b = 0; b < GA_NBAND; b++) {
		double sa = 0.0, amean, den_e, den_s, air4;
		for (w = 0; w < GA_NWALL; w++) sa += swall[w] * g->alpha[w][b];
		amean = sa / stot;
		air4 = 4.0 * g->air_e[b] * g->vol;
		if (amean >= 1.0 - 1e-12) {
			g->teyring[b] = g->tsabine[b] = 0.0;   /* 全吸音 : 残響なし */
			continue;
		}
		den_e = -stot * log(1.0 - amean) + air4;
		den_s = stot * amean + air4;
		if (den_e <= 1e-12 || den_s <= 1e-12) {
			g->teyring[b] = g->tsabine[b] = -1.0;  /* 無損失 */
			continue;
		}
		g->teyring[b] = 0.161 * g->vol / den_e;
		g->tsabine[b] = 0.161 * g->vol / den_s;
		if (g->teyring[b] > g->t60max)  g->t60max = g->teyring[b];
		if (g->tsabine[b] > g->tsabmax) g->tsabmax = g->tsabine[b];
	}
	/* どれか 1 バンドでも無損失なら「残響時間は無限大」(-1) として扱う。
	 * 逆に全バンド alpha = 1 (完全吸音) は T = 0 で、両者を取り違えると
	 * 計算時間が 3 s に張り付く / 0 になるので明確に分ける。 */
	{
		int lossless = 0;
		for (b = 0; b < GA_NBAND; b++) if (g->teyring[b] < 0.0) lossless = 1;
		if (lossless) { g->t60max = -1.0; g->tsabmax = -1.0; }
	}
	ga_log(g, "T_Eyring per band [s]: %.4g %.4g %.4g %.4g %.4g %.4g %.4g",
	       g->teyring[0], g->teyring[1], g->teyring[2], g->teyring[3],
	       g->teyring[4], g->teyring[5], g->teyring[6]);

	/* ── 計算時間 ── */
	if (g->duration_user > 0.0) {
		tsolve = g->duration_user;
		ga_log(g, "duration: %.4g s (forced by .ofdx acoustic.ga.duration_s)",
		       tsolve);
	}
	else {
		tsolve = (g->t60max >= 0.0) ? 1.5 * g->t60max : GA_TMAX;
		if (tsolve < GA_TMIN) tsolve = GA_TMIN;
		if (tsolve > GA_TMAX) tsolve = GA_TMAX;
		/* 直接音が必ず入るようにする (長大な自由空間ケース対策) */
		if (tsolve < rmax / c + 0.05) tsolve = rmax / c + 0.05;
		if (g->t60max >= 0.0)
			ga_log(g, "duration: max_band T_Eyring = %.4g s "
			       "(V = %.4g m^3, S = %.4g m^2) -> T = %.4g s",
			       g->t60max, g->vol, g->surf, tsolve);
		else
			ga_log(g, "duration: lossless room (alpha = 0, no air absorption) "
			       "-> T = %.4g s (upper clamp)", tsolve);
	}
	g->nsamples = (int)ceil(tsolve * GA_FS);
	if (g->nsamples < 16) g->nsamples = 16;
	if (g->nsamples > GA_MAX_SAMPLES) {
		ga_err(g, "requested duration %.4g s exceeds the limit of %g s",
		       tsolve, (double)GA_MAX_SAMPLES / GA_FS);
		return 1;
	}
	g->duration = (double)g->nsamples / GA_FS;

	/* ── 受音検出球の半径 ──────────────────────────────────────────
	 * 光線追跡は「球を通過したレイの本数 x 断面積の逆数」で強度を推定する。
	 * 半径は統計 (本数) と時間分解能 (通過時刻のぼけ +-R/c) のトレードオフ。
	 * 既定は室の代表寸法 V^(1/3) の 1/10 を 0.3〜1.5 m に収めた値とし、
	 * 壁・障害物・音源に食い込まないよう受音点ごとに切り詰める。 */
	rbase = (g->rsphere_user > 0.0) ? g->rsphere_user
	                                : pow(g->vol, 1.0 / 3.0) / 10.0;
	if (g->rsphere_user <= 0.0) {
		if (rbase < GA_RSPHERE_MIN) rbase = GA_RSPHERE_MIN;
		if (rbase > GA_RSPHERE_MAX) rbase = GA_RSPHERE_MAX;
	}
	for (r = 0; r < g->nrecv; r++) {
		ga_recv_t *rv = &g->recv[r];
		double lim = rbase, d;
		double dw[6];
		int i;
		dw[0] = rv->x - g->x0; dw[1] = g->x1 - rv->x;
		dw[2] = rv->y - g->y0; dw[3] = g->y1 - rv->y;
		dw[4] = rv->z - g->z0; dw[5] = g->z1 - rv->z;
		for (i = 0; i < 6; i++) if (0.9 * dw[i] < lim) lim = 0.9 * dw[i];
		d = sqrt((rv->x - g->srcx) * (rv->x - g->srcx)
		       + (rv->y - g->srcy) * (rv->y - g->srcy)
		       + (rv->z - g->srcz) * (rv->z - g->srcz));
		if (0.5 * d < lim) lim = 0.5 * d;
		for (n = 0; n < g->ngeom; n++) {
			if (!g->geom[n].ok) continue;
			d = 0.9 * geom_distance(&g->geom[n], rv->x, rv->y, rv->z);
			if (d < lim) lim = d;
		}
		if (lim < 0.02) {
			ga_err(g, "point #%d (%.4g, %.4g, %.4g) is too close to a wall, an "
			       "obstacle or the source — the ray detection sphere would be "
			       "smaller than 2 cm (move the receiver at least a few cm away)",
			       r + 1, rv->x, rv->y, rv->z);
			return 1;
		}
		rv->radius = lim;
		ga_log(g, "point #%d%s%s: (%.6g, %.6g, %.6g) m, detection sphere "
		       "r = %.4g m -> %s", r + 1, rv->name[0] ? " " : "", rv->name,
		       rv->x, rv->y, rv->z, rv->radius, rv->file);
	}
	ga_log(g, "source: (%.6g, %.6g, %.6g) m (t = 0 は音源発火時刻 — 直接音は "
	       "t = r/c に立つ)", g->srcx, g->srcy, g->srcz);

	/* ── 有効帯域 ── */
	{
		double t60ref = (g->t60max > 0.0) ? g->t60max : g->duration;
		g->flo = 2000.0 * sqrt(t60ref / g->vol);      /* Schroeder 周波数 */
		g->fhi = GA_BAND_F0 * (double)(1 << (GA_NBAND - 1)) * sqrt(2.0);
		if (g->fhi > 0.5 * GA_FS) g->fhi = 0.5 * GA_FS;
	}

	/* ── 配列確保 ── */
	g->nbins = (int)ceil(g->duration / GA_BIN_S) + 1;
	g->fftn  = next_pow2(g->nsamples + 4096);
	g->echo  = (double *)calloc((size_t)g->nrecv * g->nbins * GA_NBAND,
	                            sizeof(double));
	g->rir   = (double *)calloc((size_t)g->nrecv * g->nsamples, sizeof(double));
	g->band  = (double *)calloc((size_t)GA_NBAND * g->nsamples, sizeof(double));
	g->fr    = (double *)calloc((size_t)g->fftn, sizeof(double));
	g->fi    = (double *)calloc((size_t)g->fftn, sizeof(double));
	g->yr    = (double *)calloc((size_t)g->fftn, sizeof(double));
	g->yi    = (double *)calloc((size_t)g->fftn, sizeof(double));
	g->twr   = (double *)calloc((size_t)g->fftn / 2, sizeof(double));
	g->twi   = (double *)calloc((size_t)g->fftn / 2, sizeof(double));
	g->brev  = (int *)calloc((size_t)g->fftn, sizeof(int));
	if (!g->echo || !g->rir || !g->band || !g->fr || !g->fi || !g->yr ||
	    !g->yi || !g->twr || !g->twi || !g->brev) {
		ga_err(g, "out of memory (%d samples x %d receivers)",
		       g->nsamples, g->nrecv);
		return 1;
	}

	ga_log(g, "scattering coefficient: %.4g (0 = specular only, 1 = fully "
	       "diffuse Lambert). 鏡像法の像は (1-s)^(n/2) で減じ、拡散した分は "
	       "レイ側が受け持つ (二重計上なし)", g->scatter);
	ga_log(g, "geometric acoustics: image order %d, %d rays, "
	       "%d samples at %d Hz (%.4g s), echogram %d bins of %.4g ms",
	       g->order, g->nrays, g->nsamples, GA_FS, g->duration,
	       g->nbins, GA_BIN_S * 1e3);
	ga_log(g, "valid band: %.4g .. %.4g Hz "
	       "(下限 = Schroeder 周波数、上限 = 8 kHz バンド上端)",
	       g->flo, g->fhi);
	return 0;
}
