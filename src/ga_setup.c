/* ga_setup.c — 幾何音響ソルバーの前処理 (室・境界・空気吸収・時間軸・確保)
 *
 * ここで決まる物理量 :
 *   - 室 : .ofd のメッシュ範囲がそのまま直方体の室。格子は作らない
 *          (幾何音響は連続座標で計算する — 受音点も音源も丸めない)。
 *   - 境界 : バンド別の圧力反射係数 R_b = sqrt(1 - alpha_b)。幾何音響の
 *          標準的な「エネルギー反射率 = 1 - alpha」の定義そのもの。
 *          .ofdx で角度依存吸音を on にすると、局所反応境界として
 *          R(theta) = (zeta cos - 1)/(zeta cos + 1) を使う。zeta は
 *          吸音表の**ランダム入射** alpha から Paris の式を逆に解いて決める。
 *   - 反射面 : 室の 6 面と障害物 (AABB) の各 6 面をまとめた surf[]。
 *          鏡像法も光線追跡もこのリストだけを見る。
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

/* ── 局所反応境界の統計入射吸音率と、その逆問題 ──────────────────
 * 実数の規格化インピーダンス zeta の面に平面波が角度 theta で入るとき
 *   R(theta) = (zeta cos(theta) - 1) / (zeta cos(theta) + 1)
 *   alpha(theta) = 1 - R^2 = 4 zeta cos / (zeta cos + 1)^2
 * Paris の式 (ランダム入射 = cos 重み平均) に入れて x = cos(theta) で積分すると
 *   alpha_stat = 8 zeta int_0^1 x^2/(zeta x + 1)^2 dx
 *              = (8/zeta^2) [ zeta + 1 - 2 ln(1+zeta) - 1/(1+zeta) ]
 * という閉形式になる。zeta -> 0 と zeta -> infinity で 0、zeta = 1.55 付近で
 * 最大 0.951 を取る山型なので、逆問題は「山より右 (zeta が大きい = 空気より
 * 硬い、通常の吸音材)」の枝で二分法を使う。最大値を超える吸音率は局所反応の
 * 実インピーダンスでは表現できないので、呼び出し側で警告してクランプする。 */
double ga_alpha_stat(double zeta)
{
	if (zeta <= 0.0) return 0.0;
	return (8.0 / (zeta * zeta))
	     * (zeta + 1.0 - 2.0 * log(1.0 + zeta) - 1.0 / (1.0 + zeta));
}

double ga_zeta_from_alpha(double alpha_stat)
{
	double lo = 1.0, hi = 1.0e9, peak;
	int i;
	/* 山の位置を三分探索で求める (決定的) */
	{
		double a = 0.1, b = 10.0;
		for (i = 0; i < 200; i++) {
			double m1 = a + (b - a) / 3.0, m2 = b - (b - a) / 3.0;
			if (ga_alpha_stat(m1) < ga_alpha_stat(m2)) a = m1; else b = m2;
		}
		peak = 0.5 * (a + b);
	}
	if (alpha_stat >= ga_alpha_stat(peak)) return peak;
	if (alpha_stat <= 0.0) return hi;
	lo = peak;
	for (i = 0; i < 200; i++) {          /* alpha_stat は [peak, inf) で単調減少 */
		double mid = 0.5 * (lo + hi);
		if (ga_alpha_stat(mid) > alpha_stat) lo = mid; else hi = mid;
	}
	return 0.5 * (lo + hi);
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
	int w, b, n, r, si2;

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
	g->area = stot;

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
	{
		int nmat = 0;
		for (n = 0; n < g->ngeom; n++)
			if (g->geom[n].ok && g->geom[n].has_mat) nmat++;
		ga_log(g, "obstacles: %d geometry entries (%d with .ofdx materials, "
		       "the rest rigid alpha = 0). Their faces are part of the "
		       "image-source surface set, so their own low-order specular "
		       "reflections are included, and they also occlude other paths "
		       "and reflect rays.", g->ngeom, nmat);
	}

	/* ── 使用する音源集合 (複数音源の契約 : ADR-0010) ──
	 * 既定は feed #1 のみ (従来動作と完全一致)。multi_source なら全 feed を
	 * 強度 1 で t = 0 に同時発火し、RIR は重ね合わせになる (1/N 正規化なし)。 */
	g->nsrc = (g->multi_source && g->nfeed > 0) ? g->nfeed : 1;
	if (g->nfeed > 1 && !g->multi_source)
		ga_log(g, "warning: %d feeds found — using feed #1 only "
		       "(set acoustic.multi_source in the .ofdx to sum all sources)",
		       g->nfeed);
	if (g->multi_source)
		ga_log(g, "multi_source: %d feed(s) fire simultaneously at t = 0 — "
		       "rir.wav is the superposition (unit strength each, no 1/N)",
		       g->nsrc);

	/* ── 音源ごとのゲイン・遅延 (ADR-0010 Decision 7) ──
	 * .ofdx acoustic.sources[] が無ければ既定値 (gain = 1, delay = 0) で
	 * 確保する — 以降は常に配列越しに参照でき、既定値なら従来動作と
	 * 完全一致する (1.0 * x と x - 0.0 は IEEE でビット等価)。 */
	if (!g->srcgain) {
		g->srcgain  = (double *)malloc((size_t)(g->nfeed > 0 ? g->nfeed : 1)
		                               * sizeof(double));
		g->srcdelay = (double *)calloc((size_t)(g->nfeed > 0 ? g->nfeed : 1),
		                               sizeof(double));
		if (!g->srcgain || !g->srcdelay) {
			ga_err(g, "out of memory (%d source gains)", g->nfeed);
			return 1;
		}
		for (n = 0; n < g->nfeed; n++) g->srcgain[n] = 1.0;
	}
	g->delay_max = 0.0;
	for (si2 = 0; si2 < g->nsrc; si2++)
		if (g->srcdelay[si2] > g->delay_max) g->delay_max = g->srcdelay[si2];

	/* ── 音源・受音点の妥当性 (捏造しない : 室外・剛体内は失敗させる) ── */
	for (si2 = 0; si2 < g->nsrc; si2++) {
		double sx = g->feedpos[3 * si2 + 0];
		double sy = g->feedpos[3 * si2 + 1];
		double sz = g->feedpos[3 * si2 + 2];
		if (sx < g->x0 || sx > g->x1 || sy < g->y0 || sy > g->y1 ||
		    sz < g->z0 || sz > g->z1) {
			ga_err(g, "feed #%d position (%.4g, %.4g, %.4g) is outside the "
			       "room [%.4g,%.4g] x [%.4g,%.4g] x [%.4g,%.4g]",
			       si2 + 1, sx, sy, sz, g->x0, g->x1, g->y0, g->y1,
			       g->z0, g->z1);
			return 1;
		}
		n = geom_containing(g, sx, sy, sz);
		if (n >= 0) {
			ga_err(g, "feed #%d position (%.4g, %.4g, %.4g) is inside rigid "
			       "geometry #%d (shape %d) — move the source out of the "
			       "object (e.g. 1.5 m above the stage floor)",
			       si2 + 1, sx, sy, sz, n + 1, g->geom[n].shape);
			return 1;
		}
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
		for (si2 = 0; si2 < g->nsrc; si2++) {
			dx = rv->x - g->feedpos[3 * si2 + 0];
			dy = rv->y - g->feedpos[3 * si2 + 1];
			dz = rv->z - g->feedpos[3 * si2 + 2];
			d = sqrt(dx * dx + dy * dy + dz * dz);
			if (d < 1e-6) {
				ga_err(g, "point #%d coincides with feed #%d "
				       "(%.4g, %.4g, %.4g) — the direct sound 1/(4 pi r) "
				       "diverges", r + 1, si2 + 1, rv->x, rv->y, rv->z);
				return 1;
			}
			if (d > rmax) rmax = d;
		}
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

	/* ── 反射面リスト : 室の 6 面 + 障害物 (AABB) 1 個あたり 6 面 ──
	 * 鏡像法も光線追跡もこのリストだけを見る。障害物の面が入っているので、
	 * 障害物による低次の鏡面反射も鏡像法が扱える (v1 では遮蔽のみだった)。
	 * 障害物は剛体 (alpha = 0) — FDTD 側の geometry の扱いと揃えてある。 */
	{
		int nmax = GA_NWALL + 6 * g->ngeom;
		double rlo[3], rhi[3];
		int si = 0;
		rlo[0] = g->x0; rlo[1] = g->y0; rlo[2] = g->z0;
		rhi[0] = g->x1; rhi[1] = g->y1; rhi[2] = g->z1;
		g->surf = (ga_surf_t *)calloc((size_t)nmax, sizeof(ga_surf_t));
		if (!g->surf) {
			ga_err(g, "out of memory (%d reflecting surfaces)", nmax);
			return 1;
		}
		for (w = 0; w < GA_NWALL; w++) {
			ga_surf_t *s = &g->surf[si++];
			int ax = w / 2, hiside = w % 2;      /* GA_XM=0 -> axis 0, lo 側 */
			int u = (ax + 1) % 3, v = (ax + 2) % 3;
			s->axis  = ax;
			s->coord = hiside ? rhi[ax] : rlo[ax];
			s->nrm   = hiside ? -1.0 : 1.0;      /* 室内 (音場側) を向く法線 */
			s->lo[0] = rlo[u]; s->hi[0] = rhi[u];
			s->lo[1] = rlo[v]; s->hi[1] = rhi[v];
			for (b = 0; b < GA_NBAND; b++) {
				s->alpha[b] = g->alpha[w][b];
				s->refl[b]  = g->refl[w][b];
			}
			s->scatter = (g->wall_scatter[w] >= 0.0)
			           ? g->wall_scatter[w] : g->scatter;
			s->wall = w;
			s->geom = -1;
		}
		for (n = 0; n < g->ngeom; n++) {
			ga_geom_t *o = &g->geom[n];
			int ax;
			o->surf0 = -1;
			if (!o->ok) {
				if (o->has_mat)
					ga_log(g, "warning: .ofdx acoustic.ga.obstacles has a "
					       "material for geometry #%d, but that geometry has "
					       "an unknown shape and was ignored", n + 1);
				continue;
			}
			o->surf0 = si;      /* 面の並びは 2*axis + side (光線追跡が使う) */
			for (ax = 0; ax < 3; ax++) {
				int side, u = (ax + 1) % 3, v = (ax + 2) % 3;
				for (side = 0; side < 2; side++) {
					ga_surf_t *s = &g->surf[si++];
					s->axis  = ax;
					s->coord = side ? o->hi[ax] : o->lo[ax];
					s->nrm   = side ? 1.0 : -1.0;   /* 障害物の外向き = 音場側 */
					s->lo[0] = o->lo[u]; s->hi[0] = o->hi[u];
					s->lo[1] = o->lo[v]; s->hi[1] = o->hi[v];
					/* 既定は剛体 (alpha = 0)。.ofdx acoustic.ga.obstacles で
					 * 材質があればバンド別に使う。FDTD 側は常に剛体なので、
					 * 材質を与えると障害物の扱いは高域のみ変わる (ReadMe)。 */
					for (b = 0; b < GA_NBAND; b++) {
						double a = o->has_mat ? o->mat_alpha[b] : 0.0;
						s->alpha[b] = a;
						s->refl[b]  = sqrt(1.0 - a);
					}
					s->scatter = (o->has_mat && o->mat_scatter >= 0.0)
					           ? o->mat_scatter : g->scatter;
					s->wall = -1;
					s->geom = n;
				}
			}
		}
		g->nsurf = si;
	}

	/* ── 角度依存吸音 (局所反応) 用の規格化インピーダンス ──
	 * 吸音表の alpha は**ランダム入射**の値なので、Paris の式を逆に解いて
	 * zeta を決める (垂直入射の値として使うと二重に効いてしまう)。
	 * 局所反応の実インピーダンスで表せる上限は alpha_stat = 0.951 なので、
	 * それを超える値は正直にクランプして警告する。 */
	{
		int clamped = 0;
		double amax = ga_alpha_stat(ga_zeta_from_alpha(2.0));   /* 山の高さ */
		for (n = 0; n < g->nsurf; n++) {
			ga_surf_t *s = &g->surf[n];
			for (b = 0; b < GA_NBAND; b++) {
				double a = s->alpha[b];
				if (a > amax) { a = amax; clamped = 1; }
				s->zeta[b] = ga_zeta_from_alpha(a);
			}
		}
		if (g->angle_dep) {
			ga_log(g, "angle-dependent absorption: locally reacting boundary "
			       "R(theta) = (zeta cos - 1)/(zeta cos + 1); zeta solved from "
			       "the random-incidence alpha via Paris' formula "
			       "(max representable alpha_stat = %.4f)", amax);
			if (clamped)
				ga_log(g, "warning: some alpha exceed %.4f, the largest "
				       "random-incidence absorption a locally reacting real "
				       "impedance can have — clamped (the surface is treated "
				       "as the most absorbing such boundary, not as alpha = 1)",
				       amax);
			for (w = 0; w < GA_NWALL; w++)
				ga_log(g, "wall %s: zeta = [%.4g %.4g %.4g %.4g %.4g %.4g %.4g]",
				       (w == GA_XM) ? "x-" : (w == GA_XP) ? "x+" :
				       (w == GA_YM) ? "y-" : (w == GA_YP) ? "y+" :
				       (w == GA_ZM) ? "z-" : "z+",
				       g->surf[w].zeta[0], g->surf[w].zeta[1], g->surf[w].zeta[2],
				       g->surf[w].zeta[3], g->surf[w].zeta[4], g->surf[w].zeta[5],
				       g->surf[w].zeta[6]);
		}
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
		if (tsolve < g->delay_max + rmax / c)
			ga_log(g, "warning: forced duration %.4g s is shorter than "
			       "max source delay + longest direct path (%.4g s) — "
			       "some arrivals will fall outside the RIR",
			       tsolve, g->delay_max + rmax / c);
	}
	else {
		tsolve = (g->t60max >= 0.0) ? 1.5 * g->t60max : GA_TMAX;
		if (tsolve < GA_TMIN) tsolve = GA_TMIN;
		if (tsolve > GA_TMAX) tsolve = GA_TMAX;
		/* 音源遅延 (ADR-0010 Decision 7) : 最後の発火から数えて T 残すよう
		 * max(delay) を足す (遅延で残響が窓の外にこぼれない)。delay = 0 なら
		 * 従来と完全一致。 */
		tsolve += g->delay_max;
		/* 直接音が必ず入るようにする (長大な自由空間ケース対策) */
		if (tsolve < g->delay_max + rmax / c + 0.05)
			tsolve = g->delay_max + rmax / c + 0.05;
		if (g->delay_max > 0.0)
			ga_log(g, "duration: extended by max source delay %.4g s",
			       g->delay_max);
		if (g->t60max >= 0.0)
			ga_log(g, "duration: max_band T_Eyring = %.4g s "
			       "(V = %.4g m^3, S = %.4g m^2) -> T = %.4g s",
			       g->t60max, g->vol, g->area, tsolve);
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
		for (si2 = 0; si2 < g->nsrc; si2++) {
			d = sqrt((rv->x - g->feedpos[3 * si2 + 0])
			           * (rv->x - g->feedpos[3 * si2 + 0])
			       + (rv->y - g->feedpos[3 * si2 + 1])
			           * (rv->y - g->feedpos[3 * si2 + 1])
			       + (rv->z - g->feedpos[3 * si2 + 2])
			           * (rv->z - g->feedpos[3 * si2 + 2]));
			if (0.5 * d < lim) lim = 0.5 * d;
		}
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
	for (si2 = 0; si2 < g->nsrc; si2++)
		ga_log(g, "source #%d: (%.6g, %.6g, %.6g) m, gain = %.4g, "
		       "delay = %.4g s%s", si2 + 1,
		       g->feedpos[3 * si2 + 0], g->feedpos[3 * si2 + 1],
		       g->feedpos[3 * si2 + 2], g->srcgain[si2], g->srcdelay[si2],
		       (si2 == 0) ? " (発火は t = delay — 直接音は t = delay + r/c)"
		                  : "");

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

	ga_log(g, "scattering coefficient: default %.4g, per wall "
	       "[x- %.4g, x+ %.4g, y- %.4g, y+ %.4g, z- %.4g, z+ %.4g] "
	       "(0 = specular only, 1 = fully diffuse Lambert). 鏡像法の像は "
	       "面ごとの (1-s)^(1/2) を掛けて減じ、抜けた拡散分はレイ側が受け持つ "
	       "(二重計上なし)", g->scatter,
	       g->surf[GA_XM].scatter, g->surf[GA_XP].scatter,
	       g->surf[GA_YM].scatter, g->surf[GA_YP].scatter,
	       g->surf[GA_ZM].scatter, g->surf[GA_ZP].scatter);
	ga_log(g, "reflecting surfaces: %d (room 6 + %d obstacle faces)",
	       g->nsurf, g->nsurf - GA_NWALL);
	ga_log(g, "geometric acoustics: image order %d, %d rays, "
	       "%d samples at %d Hz (%.4g s), echogram %d bins of %.4g ms",
	       g->order, g->nrays, g->nsamples, GA_FS, g->duration,
	       g->nbins, GA_BIN_S * 1e3);
	ga_log(g, "valid band: %.4g .. %.4g Hz "
	       "(下限 = Schroeder 周波数、上限 = 8 kHz バンド上端)",
	       g->flo, g->fhi);
	return 0;
}
