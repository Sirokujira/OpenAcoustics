/* ga_shape.c — 幾何音響の障害物の厳密形状 (本家 OpenFDTD の全 shape)
 *
 * パラメータの並びは本家 OpenFDTD sol/ingeometry.c と同じで、FDTD 側
 * (fdtd.c の geom_inside) と**同じ規則**で内外を決める (境界を含む)。
 * 低域 (FDTD) と高域 (幾何音響) が同じ障害物を見ることがハイブリッド合成の
 * クロスオーバー整合の前提なので、片方だけ規則を変えないこと。
 *
 *   shape 1        直方体 [x1,x2] x [y1,y2] x [z1,z2]
 *   shape 2        その箱に内接する楕円体
 *   shape 11..13   x/y/z 軸の楕円柱 (断面は箱に内接する楕円)
 *   shape 31..33   x/y/z 軸の三角柱 : g[0..1] 軸の範囲、頂点の (u, w) は
 *                  (g[2], g[5]), (g[3], g[6]), (g[4], g[7])
 *   shape 41..43   x/y/z 軸の角錐台 : g[0..1] 軸の両端、g[2..3] 断面中心
 *                  (u, w)、g[4..5] は g[0] 側の全幅、g[6..7] は g[1] 側の全幅
 *   shape 51..53   同じ並びの円錐台 (楕円断面)
 *   断面の 2 軸は u = (軸+1)%3, w = (軸+2)%3。
 *
 * 交差の求め方:
 *   - 楕円体・楕円柱 : 2 次方程式 (閉形式)。楕円柱は側面と両端の平面の
 *     区間の共通部分。
 *   - 三角柱・角錐台 : 平面の凸多面体なので、各面の半空間で直線を切る
 *     (Cyrus–Beck)。入った面が反射面 (鏡像法と同じ面データ)。
 *   - 円錐台 : 断面の半幅が軸方向に線形に変わる楕円なので、境界は
 *     G(t) = U^2 B^2 + W^2 A^2 - A^2 B^2 = 0 (U, W は断面中心からのずれ、
 *     A, B は半幅 — どれも t の 1 次式) の 4 次方程式になる。両端の平面で
 *     切った区間を 64 等分して符号の変わり目を探し、二分法で詰める
 *     (区間の 1/64 より短い「かすめる」交差は取りこぼし得る — 面積 0 に
 *     近い寄与なので許容)。法線は F = (U/A)^2 + (W/B)^2 の勾配。
 *     断面の縦横比が両端で違う円錐台は凸とは限らないので、最初の入口を
 *     この符号探索で求める (凸を仮定した 2 根の公式は使わない)。
 * どれも乱数を使わない決定的な計算。
 */
#include <math.h>
#include <string.h>
#include "ga.h"

#define SH_TINY   1e-15
#define SH_CHORD  1e-9     /* これより短い弦 (かすめる接触) は交差としない [m] */
#define SH_SCAN   64       /* 円錐台の符号探索の分割数 */
#define SH_BISECT 64       /* 円錐台の二分法の回数 */
#define SH_NDIR   1024     /* 距離の推定に使う方向の数 (球面フィボナッチ格子) */

static double dot3(const double a[3], const double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static int unit3(double a[3])
{
	double l = sqrt(dot3(a, a));
	if (!(l > 0.0)) return 0;
	a[0] /= l; a[1] /= l; a[2] /= l;
	return 1;
}

/* (軸 ax の座標 a, 断面の u 座標, w 座標) から 3 次元の点 */
static void axpt(int ax, double a, double cu, double cw, double P[3])
{
	P[ax] = a;
	P[(ax + 1) % 3] = cu;
	P[(ax + 2) % 3] = cw;
}

/* ── 外接箱 (従来の AABB と同じ — 交差判定の前段と距離の下限に使う) ── */
static int solid_aabb(ga_geom_t *o)
{
	const double *p = o->g;
	double *lo = o->lo, *hi = o->hi;
	int i;
	double a, b;
	switch (o->shape) {
		case 1: case 2: case 11: case 12: case 13:
			for (i = 0; i < 3; i++) {
				a = p[2*i]; b = p[2*i + 1];
				lo[i] = (a < b) ? a : b;
				hi[i] = (a < b) ? b : a;
			}
			return 1;
		case 31: case 32: case 33: {
			int ax = o->shape - 31;
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
			int ax = (o->shape % 10) - 1;
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

/* ── 多面体の面を 1 枚足す ──
 * 法線は Newell の方法 (縮退した 4 角形 = 3 角形でも安定)、向きは立体の
 * 内点 cen から遠ざかる側 (外向き) にそろえ、頂点はその法線のまわりに
 * 反時計回りに並べ直す。面積 0 の面は fallback (外向きの単位法線) で
 * 半空間だけを持ち、鏡像法には入れない (fnv = 0)。 */
static void poly_face(ga_geom_t *o, int nv, double v[4][3], const double cen[3],
                      const double fallback[3])
{
	int f = o->nface++, i;
	double N[3] = { 0.0, 0.0, 0.0 }, vc[3] = { 0.0, 0.0, 0.0 }, rel[3];
	for (i = 0; i < nv; i++) {
		const double *a = v[i], *b = v[(i + 1) % nv];
		N[0] += (a[1] - b[1]) * (a[2] + b[2]);
		N[1] += (a[2] - b[2]) * (a[0] + b[0]);
		N[2] += (a[0] - b[0]) * (a[1] + b[1]);
		vc[0] += a[0] / nv; vc[1] += a[1] / nv; vc[2] += a[2] / nv;
	}
	rel[0] = vc[0] - cen[0]; rel[1] = vc[1] - cen[1]; rel[2] = vc[2] - cen[2];
	if (sqrt(dot3(N, N)) > 1e-12) {
		if (dot3(N, rel) < 0.0) {             /* 内向き -> 反転 */
			double tmp[4][3];
			N[0] = -N[0]; N[1] = -N[1]; N[2] = -N[2];
			for (i = 0; i < nv; i++) memcpy(tmp[i], v[nv - 1 - i], sizeof(tmp[i]));
			for (i = 0; i < nv; i++) memcpy(v[i], tmp[i], sizeof(tmp[i]));
		}
		unit3(N);
		o->fnv[f] = nv;
	}
	else {
		N[0] = fallback[0]; N[1] = fallback[1]; N[2] = fallback[2];
		o->fnv[f] = 0;                        /* 面積 0 : 鏡像法には入れない */
	}
	memcpy(o->fn[f], N, sizeof(N));
	o->fd[f] = dot3(N, vc);
	for (i = 0; i < nv; i++) memcpy(o->fv[f][i], v[i], sizeof(v[i]));
}

/* 三角柱 (31..33) の 5 面 */
static int poly_prism(ga_geom_t *o)
{
	const double *p = o->g;
	int ax = o->shape - 31, k;
	double a0 = (p[0] < p[1]) ? p[0] : p[1];
	double a1 = (p[0] < p[1]) ? p[1] : p[0];
	double area2 = (p[3] - p[2]) * (p[7] - p[5]) - (p[6] - p[5]) * (p[4] - p[2]);
	double cen[3], fb[3], v[4][3];
	if (area2 == 0.0 || !(a1 > a0)) return 0;          /* 体積 0 */
	axpt(ax, 0.5 * (a0 + a1), (p[2] + p[3] + p[4]) / 3.0,
	     (p[5] + p[6] + p[7]) / 3.0, cen);
	o->nface = 0;
	fb[0] = fb[1] = fb[2] = 0.0;
	for (k = 0; k < 3; k++) axpt(ax, a0, p[2 + k], p[5 + k], v[k]);
	fb[ax] = -1.0; poly_face(o, 3, v, cen, fb);
	for (k = 0; k < 3; k++) axpt(ax, a1, p[2 + k], p[5 + k], v[k]);
	fb[ax] = 1.0;  poly_face(o, 3, v, cen, fb);
	fb[ax] = 0.0;
	for (k = 0; k < 3; k++) {
		int k2 = (k + 1) % 3;
		axpt(ax, a0, p[2 + k],  p[5 + k],  v[0]);
		axpt(ax, a0, p[2 + k2], p[5 + k2], v[1]);
		axpt(ax, a1, p[2 + k2], p[5 + k2], v[2]);
		axpt(ax, a1, p[2 + k],  p[5 + k],  v[3]);
		poly_face(o, 4, v, cen, fb);
	}
	return 1;
}

/* 角錐台 (41..43) の 6 面 */
static int poly_frustum(ga_geom_t *o)
{
	const double *p = o->g;
	int ax = (o->shape % 10) - 1, e, su, sw;
	double hu[2], hw[2], ea[2], cen[3], fb[3], v[4][3];
	ea[0] = p[0]; ea[1] = p[1];
	hu[0] = 0.5 * p[4]; hw[0] = 0.5 * p[5];
	hu[1] = 0.5 * p[6]; hw[1] = 0.5 * p[7];
	if (p[1] == p[0] || hu[0] < 0.0 || hw[0] < 0.0 || hu[1] < 0.0 || hw[1] < 0.0 ||
	    (hu[0] == 0.0 && hu[1] == 0.0) || (hw[0] == 0.0 && hw[1] == 0.0))
		return 0;                                       /* 体積 0 / 不正な幅 */
	axpt(ax, 0.5 * (p[0] + p[1]), p[2], p[3], cen);
	o->nface = 0;
	for (e = 0; e < 2; e++) {                           /* 両端の矩形 */
		fb[0] = fb[1] = fb[2] = 0.0;
		fb[ax] = (ea[e] > ea[1 - e]) ? 1.0 : -1.0;
		axpt(ax, ea[e], p[2] - hu[e], p[3] - hw[e], v[0]);
		axpt(ax, ea[e], p[2] + hu[e], p[3] - hw[e], v[1]);
		axpt(ax, ea[e], p[2] + hu[e], p[3] + hw[e], v[2]);
		axpt(ax, ea[e], p[2] - hu[e], p[3] + hw[e], v[3]);
		poly_face(o, 4, v, cen, fb);
	}
	fb[0] = fb[1] = fb[2] = 0.0;
	for (su = -1; su <= 1; su += 2) {                   /* +-u 側面 */
		axpt(ax, ea[0], p[2] + su * hu[0], p[3] - hw[0], v[0]);
		axpt(ax, ea[0], p[2] + su * hu[0], p[3] + hw[0], v[1]);
		axpt(ax, ea[1], p[2] + su * hu[1], p[3] + hw[1], v[2]);
		axpt(ax, ea[1], p[2] + su * hu[1], p[3] - hw[1], v[3]);
		poly_face(o, 4, v, cen, fb);
	}
	for (sw = -1; sw <= 1; sw += 2) {                   /* +-w 側面 */
		axpt(ax, ea[0], p[2] - hu[0], p[3] + sw * hw[0], v[0]);
		axpt(ax, ea[0], p[2] + hu[0], p[3] + sw * hw[0], v[1]);
		axpt(ax, ea[1], p[2] + hu[1], p[3] + sw * hw[1], v[2]);
		axpt(ax, ea[1], p[2] - hu[1], p[3] + sw * hw[1], v[3]);
		poly_face(o, 4, v, cen, fb);
	}
	return 1;
}

int ga_solid_prepare(ga_geom_t *o)
{
	const double *p = o->g;
	int a;
	o->nface = 0;
	if (!solid_aabb(o)) return 0;
	switch (o->shape) {
		case 1:
			o->kind = GA_SOLID_BOX;
			return 1;
		case 2: case 11: case 12: case 13:
			for (a = 0; a < 3; a++)
				if (!(fabs(p[2*a + 1] - p[2*a]) > 0.0)) return 0;   /* 体積 0 */
			o->kind = GA_SOLID_CURVED;
			return 1;
		case 31: case 32: case 33:
			o->kind = GA_SOLID_POLY;
			return poly_prism(o);
		case 41: case 42: case 43:
			o->kind = GA_SOLID_POLY;
			return poly_frustum(o);
		case 51: case 52: case 53:
			if (p[1] == p[0] || p[4] < 0.0 || p[5] < 0.0 || p[6] < 0.0 ||
			    p[7] < 0.0 || (p[4] == 0.0 && p[6] == 0.0) ||
			    (p[5] == 0.0 && p[7] == 0.0))
				return 0;
			o->kind = GA_SOLID_CURVED;
			return 1;
		default:
			return 0;
	}
}

/* ── 内外判定 (境界を含む) — FDTD 側 fdtd.c の geom_inside と同じ規則 ── */
int ga_solid_inside(const ga_geom_t *o, double x, double y, double z)
{
	const double *p = o->g;
	double q[3];
	if (x < o->lo[0] || x > o->hi[0] || y < o->lo[1] || y > o->hi[1] ||
	    z < o->lo[2] || z > o->hi[2])
		return 0;
	q[0] = x; q[1] = y; q[2] = z;
	switch (o->shape) {
		case 1:
			return 1;
		case 2: {
			double sum = 0.0;
			int a;
			for (a = 0; a < 3; a++) {
				double c = 0.5 * (p[2*a] + p[2*a + 1]);
				double r = 0.5 * fabs(p[2*a + 1] - p[2*a]);
				if (r <= 0.0) return 0;
				sum += ((q[a] - c) / r) * ((q[a] - c) / r);
			}
			return sum <= 1.0;
		}
		case 11: case 12: case 13: {
			int ax = o->shape - 11, a;
			double sum = 0.0;
			for (a = 0; a < 3; a++) {
				double c, r;
				if (a == ax) continue;
				c = 0.5 * (p[2*a] + p[2*a + 1]);
				r = 0.5 * fabs(p[2*a + 1] - p[2*a]);
				if (r <= 0.0) return 0;
				sum += ((q[a] - c) / r) * ((q[a] - c) / r);
			}
			return sum <= 1.0;
		}
		case 31: case 32: case 33: {
			int ax = o->shape - 31;
			int u = (ax + 1) % 3, w = (ax + 2) % 3, k;
			double area2 = (p[3] - p[2]) * (p[7] - p[5])
			             - (p[6] - p[5]) * (p[4] - p[2]);
			int neg = 0, pos = 0;
			if (area2 == 0.0) return 0;
			for (k = 0; k < 3; k++) {
				int k2 = (k + 1) % 3;
				double cr = (p[2+k2] - p[2+k]) * (q[w] - p[5+k])
				          - (p[5+k2] - p[5+k]) * (q[u] - p[2+k]);
				if (cr < 0.0) neg = 1;
				if (cr > 0.0) pos = 1;
			}
			return !(neg && pos);
		}
		case 41: case 42: case 43:
		case 51: case 52: case 53: {
			int ax = (o->shape % 10) - 1;
			int u = (ax + 1) % 3, w = (ax + 2) % 3;
			double len = p[1] - p[0], t, hu, hw, du, dw;
			if (len == 0.0) return 0;
			t = (q[ax] - p[0]) / len;
			if (t < 0.0 || t > 1.0) return 0;
			hu = 0.5 * (p[4] + (p[6] - p[4]) * t);
			hw = 0.5 * (p[5] + (p[7] - p[5]) * t);
			du = q[u] - p[2];
			dw = q[w] - p[3];
			if (o->shape < 50)
				return fabs(du) <= hu && fabs(dw) <= hw;
			if (hu <= 0.0 || hw <= 0.0) return 0;
			return (du / hu) * (du / hu) + (dw / hw) * (dw / hw) <= 1.0;
		}
		default:
			return 0;
	}
}

/* ── 直線との交差 ───────────────────────────────────────────────── */

/* 2 次方程式 a t^2 + b t + c = 0 の 2 根 (t0 <= t1)。実根が無い / 接するだけ
 * なら 0 (数値的に安定な形 : 桁落ちの少ない側から求める)。 */
static int quad_roots(double a, double b, double c, double *t0, double *t1)
{
	double disc = b * b - 4.0 * a * c, sq, qq, r0, r1;
	if (!(a > 0.0) || !(disc > 0.0)) return 0;
	sq = sqrt(disc);
	qq = -0.5 * (b + ((b >= 0.0) ? sq : -sq));
	r0 = qq / a;
	r1 = (qq != 0.0) ? c / qq : -r0;
	if (r0 > r1) { double t = r0; r0 = r1; r1 = t; }
	*t0 = r0; *t1 = r1;
	return 1;
}

/* 区間 [tE, tL] (形状の内側) から「t in (tmin, tmax) で入る」を判定 */
static int accept_entry(double tE, double tL, double tmin, double tmax, double dl)
{
	if (!(tE > tmin) || !(tE < tmax)) return 0;      /* 窓の外 / 始点が内側 */
	if (!((tL - tE) * dl > SH_CHORD)) return 0;       /* かすめるだけ */
	return 1;
}

/* 箱 (shape 1) : 面の半空間 (slab) で切る。面番号は 2*axis + side */
static int box_entry(const ga_geom_t *o, const double p[3], const double d[3],
                     double tmin, double tmax, double dl,
                     double *t, double nrm[3], int *face)
{
	double tE = -1e300, tL = 1e300;
	int i, fa = -1, fs = 0;
	for (i = 0; i < 3; i++) {
		if (fabs(d[i]) < SH_TINY) {
			if (p[i] < o->lo[i] || p[i] > o->hi[i]) return 0;
		}
		else {
			double t1 = (o->lo[i] - p[i]) / d[i];
			double t2 = (o->hi[i] - p[i]) / d[i];
			if (t1 > t2) { double tt = t1; t1 = t2; t2 = tt; }
			if (t1 > tE) { tE = t1; fa = i; fs = (d[i] > 0.0) ? 0 : 1; }
			if (t2 < tL) tL = t2;
		}
	}
	if (fa < 0 || tE > tL || !accept_entry(tE, tL, tmin, tmax, dl)) return 0;
	*t = tE;
	nrm[0] = nrm[1] = nrm[2] = 0.0;
	nrm[fa] = fs ? 1.0 : -1.0;
	*face = 2 * fa + fs;
	return 1;
}

/* 凸多面体 (Cyrus–Beck) */
static int poly_entry(const ga_geom_t *o, const double p[3], const double d[3],
                      double tmin, double tmax, double dl,
                      double *t, double nrm[3], int *face)
{
	double tE = -1e300, tL = 1e300;
	int f, fE = -1;
	for (f = 0; f < o->nface; f++) {
		double den = dot3(o->fn[f], d);
		double num = o->fd[f] - dot3(o->fn[f], p);    /* >= 0 : 半空間の内側 */
		if (fabs(den) < SH_TINY) {
			if (num < 0.0) return 0;                    /* 平行で外側 */
			continue;
		}
		if (den < 0.0) {
			double tt = num / den;
			if (tt > tE) { tE = tt; fE = f; }
		}
		else {
			double tt = num / den;
			if (tt < tL) tL = tt;
		}
		if (tE > tL) return 0;
	}
	if (fE < 0 || !accept_entry(tE, tL, tmin, tmax, dl)) return 0;
	*t = tE;
	memcpy(nrm, o->fn[fE], 3 * sizeof(double));
	*face = fE;
	return 1;
}

/* 楕円体 (2) */
static int ellipsoid_entry(const ga_geom_t *o, const double p[3], const double d[3],
                           double tmin, double tmax, double dl,
                           double *t, double nrm[3], int *face)
{
	const double *g = o->g;
	double a = 0.0, b = 0.0, c = -1.0, c3[3], r2[3], t0, t1, X;
	int i;
	for (i = 0; i < 3; i++) {
		double r = 0.5 * fabs(g[2*i + 1] - g[2*i]);
		double q;
		c3[i] = 0.5 * (g[2*i] + g[2*i + 1]);
		r2[i] = r * r;
		q = p[i] - c3[i];
		a += d[i] * d[i] / r2[i];
		b += 2.0 * q * d[i] / r2[i];
		c += q * q / r2[i];
	}
	if (!quad_roots(a, b, c, &t0, &t1)) return 0;
	if (!accept_entry(t0, t1, tmin, tmax, dl)) return 0;
	*t = t0;
	for (i = 0; i < 3; i++) {
		X = p[i] + t0 * d[i];
		nrm[i] = (X - c3[i]) / r2[i];
	}
	unit3(nrm);
	*face = -1;
	return 1;
}

/* 楕円柱 (11..13) : 側面 (2 次) と両端の平面の区間の共通部分 */
static int cylinder_entry(const ga_geom_t *o, const double p[3], const double d[3],
                          double tmin, double tmax, double dl,
                          double *t, double nrm[3], int *face)
{
	const double *g = o->g;
	int ax = o->shape - 11, k, cap = 0;
	int ab[2];
	double a = 0.0, b = 0.0, c = -1.0, c3[3], r2[3];
	double l0 = -1e300, l1 = 1e300, s0 = -1e300, s1 = 1e300, tE, tL;
	ab[0] = (ax + 1) % 3; ab[1] = (ax + 2) % 3;
	for (k = 0; k < 2; k++) {
		int i = ab[k];
		double r = 0.5 * fabs(g[2*i + 1] - g[2*i]), q;
		c3[i] = 0.5 * (g[2*i] + g[2*i + 1]);
		r2[i] = r * r;
		q = p[i] - c3[i];
		a += d[i] * d[i] / r2[i];
		b += 2.0 * q * d[i] / r2[i];
		c += q * q / r2[i];
	}
	if (a > 1e-30) {
		if (!quad_roots(a, b, c, &l0, &l1)) return 0;
	}
	else if (c > 0.0) {
		return 0;                                        /* 軸に平行で外側 */
	}
	if (fabs(d[ax]) < SH_TINY) {
		if (p[ax] < o->lo[ax] || p[ax] > o->hi[ax]) return 0;
	}
	else {
		s0 = (o->lo[ax] - p[ax]) / d[ax];
		s1 = (o->hi[ax] - p[ax]) / d[ax];
		if (s0 > s1) { double tt = s0; s0 = s1; s1 = tt; }
	}
	tE = l0; tL = l1;
	if (s0 > tE) { tE = s0; cap = 1; }
	if (s1 < tL) tL = s1;
	if (tE > tL || !accept_entry(tE, tL, tmin, tmax, dl)) return 0;
	*t = tE;
	nrm[0] = nrm[1] = nrm[2] = 0.0;
	if (cap) {
		nrm[ax] = (d[ax] > 0.0) ? -1.0 : 1.0;
	}
	else {
		for (k = 0; k < 2; k++) {
			int i = ab[k];
			nrm[i] = (p[i] + tE * d[i] - c3[i]) / r2[i];
		}
		unit3(nrm);
	}
	*face = -1;
	return 1;
}

/* 円錐台 (51..53) : 楕円断面の半幅 A, B が軸方向に線形 */
typedef struct {
	int    ax, u, w;
	double len, s0, s1;          /* s(t) = s0 + s1 t (0..1 が両端) */
	double U0, U1, W0, W1;       /* 断面中心からのずれ */
	double A0, A1, B0, B1;       /* 半幅 */
} cone_t;

static int cone_inside_t(const cone_t *k, double t)
{
	double s = k->s0 + k->s1 * t;
	double A = k->A0 + k->A1 * t, B = k->B0 + k->B1 * t;
	double U = k->U0 + k->U1 * t, W = k->W0 + k->W1 * t;
	if (s < 0.0 || s > 1.0 || !(A > 0.0) || !(B > 0.0)) return 0;
	return U * U * B * B + W * W * A * A - A * A * B * B <= 0.0;
}

static int cone_entry(const ga_geom_t *o, const double p[3], const double d[3],
                      double tmin, double tmax, double dl,
                      double *t, double nrm[3], int *face)
{
	const double *g = o->g;
	cone_t k;
	double sa = -1e300, sb = 1e300, lo, hi, tE = 0.0;
	int i, j, it, src = 0, how = 0;   /* how : 1 = 端の平面から、2 = 側面から */
	(void)dl;
	k.ax = (o->shape % 10) - 1;
	k.u = (k.ax + 1) % 3; k.w = (k.ax + 2) % 3;
	k.len = g[1] - g[0];
	k.s0 = (p[k.ax] - g[0]) / k.len;
	k.s1 = d[k.ax] / k.len;
	k.A0 = 0.5 * (g[4] + (g[6] - g[4]) * k.s0); k.A1 = 0.5 * (g[6] - g[4]) * k.s1;
	k.B0 = 0.5 * (g[5] + (g[7] - g[5]) * k.s0); k.B1 = 0.5 * (g[7] - g[5]) * k.s1;
	k.U0 = p[k.u] - g[2]; k.U1 = d[k.u];
	k.W0 = p[k.w] - g[3]; k.W1 = d[k.w];
	if (fabs(k.s1) < SH_TINY) {
		if (k.s0 < 0.0 || k.s0 > 1.0) return 0;
	}
	else {
		sa = -k.s0 / k.s1;
		sb = (1.0 - k.s0) / k.s1;
		if (sa > sb) { double tt = sa; sa = sb; sb = tt; }
	}
	/* 探索区間 [lo, hi] : 窓 (tmin, tmax)・両端の平面・外接箱の側面の共通部分。
	 * 外接箱でも絞るのは、大半の直線がここで外れて符号探索を省けるのと、
	 * 探索の刻みが細かくなるため (外接箱は形状を含むので入口を失わない)。
	 * lo を決めたもの (src : 0 = 窓の始点、1 = 端の平面、2 = 外接箱の側面) を
	 * 覚えておく — 入口の直後で既に内側なら、そこが入った面になる。 */
	lo = tmin; hi = tmax;
	if (sa > lo) { lo = sa; src = 1; }
	if (sb < hi) hi = sb;
	for (i = 0; i < 3; i++) {
		double t1, t2;
		if (i == k.ax) continue;                /* 軸方向は端の平面で済んでいる */
		if (fabs(d[i]) < SH_TINY) {
			if (p[i] < o->lo[i] || p[i] > o->hi[i]) return 0;
			continue;
		}
		t1 = (o->lo[i] - p[i]) / d[i];
		t2 = (o->hi[i] - p[i]) / d[i];
		if (t1 > t2) { double tt = t1; t1 = t2; t2 = tt; }
		if (t1 > lo) { lo = t1; src = 2; }
		if (t2 < hi) hi = t2;
	}
	if (!(lo < hi)) return 0;

	/* ちょうど lo は s の丸めで [0,1] の外に出得るので僅かに先を見る */
	if (cone_inside_t(&k, lo + (hi - lo) * 1e-9)) {
		if (src == 0) return 0;                 /* 始点が内側 : 交差なし */
		tE = lo;
		how = (src == 1) ? 1 : 2;
	}
	else {
		/* 区間を刻んで最初に内側へ入る点を探し、二分法で詰める */
		double tprev = lo;
		for (j = 1; j <= SH_SCAN && !how; j++) {
			double tj = lo + (hi - lo) * (double)j / (double)SH_SCAN;
			if (cone_inside_t(&k, tj)) {
				double ta = tprev, tb = tj;       /* 外 / 内 */
				for (it = 0; it < SH_BISECT; it++) {
					double tm = 0.5 * (ta + tb);
					if (cone_inside_t(&k, tm)) tb = tm; else ta = tm;
				}
				tE = tb;
				how = 2;
			}
			tprev = tj;
		}
		if (!how) return 0;
	}
	if (!(tE > tmin) || !(tE < tmax)) return 0;

	nrm[0] = nrm[1] = nrm[2] = 0.0;
	if (how == 1) {
		/* s が 0 の端 (g[0] 側) から入ったか 1 の端 (g[1] 側) から入ったか */
		if (k.s1 > 0.0) nrm[k.ax] = (g[0] > g[1]) ? 1.0 : -1.0;
		else            nrm[k.ax] = (g[1] > g[0]) ? 1.0 : -1.0;
	}
	else {
		/* 側面 : 法線 = F = (U/A)^2 + (W/B)^2 の勾配 */
		double s = k.s0 + k.s1 * tE;
		double A = 0.5 * (g[4] + (g[6] - g[4]) * s);
		double B = 0.5 * (g[5] + (g[7] - g[5]) * s);
		double U = k.U0 + k.U1 * tE, W = k.W0 + k.W1 * tE;
		double dA = 0.5 * (g[6] - g[4]) / k.len, dB = 0.5 * (g[7] - g[5]) / k.len;
		nrm[k.u]  = 2.0 * U / (A * A);
		nrm[k.w]  = 2.0 * W / (B * B);
		nrm[k.ax] = -2.0 * U * U * dA / (A * A * A) - 2.0 * W * W * dB / (B * B * B);
		if (!unit3(nrm)) return 0;
	}
	*t = tE;
	*face = -1;
	return 1;
}

int ga_solid_entry(const ga_geom_t *o, const double p[3], const double d[3],
                   double tmin, double tmax, double *t, double nrm[3], int *face)
{
	double dl = sqrt(dot3(d, d));
	if (!o->ok || !(dl > 0.0)) return 0;
	switch (o->shape) {
		case 1:
			return box_entry(o, p, d, tmin, tmax, dl, t, nrm, face);
		case 2:
			return ellipsoid_entry(o, p, d, tmin, tmax, dl, t, nrm, face);
		case 11: case 12: case 13:
			return cylinder_entry(o, p, d, tmin, tmax, dl, t, nrm, face);
		case 31: case 32: case 33:
		case 41: case 42: case 43:
			return poly_entry(o, p, d, tmin, tmax, dl, t, nrm, face);
		case 51: case 52: case 53:
			return cone_entry(o, p, d, tmin, tmax, dl, t, nrm, face);
		default:
			return 0;
	}
}

/* ── 点から形状までの距離の推定 (受音球の半径の上限用) ──
 * 箱は外接箱までの距離そのもの (厳密 — 従来と同じ値)。その他の形状は
 * 球面フィボナッチ格子の 1024 方向へ直線を伸ばして最初に当たる距離の最小値
 * (真の距離より僅かに大きく出得る — 呼び出し側が 0.9 倍の余裕を取る)。
 * 内側なら 0。 */
double ga_solid_distance(const ga_geom_t *o, double x, double y, double z)
{
	const double golden = GA_PI * (3.0 - sqrt(5.0));
	double p[3], best = 1e300;
	int i;
	p[0] = x; p[1] = y; p[2] = z;
	if (o->kind == GA_SOLID_BOX) {
		double dsum = 0.0;
		for (i = 0; i < 3; i++) {
			double e = 0.0;
			if (p[i] < o->lo[i]) e = o->lo[i] - p[i];
			else if (p[i] > o->hi[i]) e = p[i] - o->hi[i];
			dsum += e * e;
		}
		return sqrt(dsum);
	}
	if (ga_solid_inside(o, x, y, z)) return 0.0;
	for (i = 0; i < SH_NDIR; i++) {
		double zz = 1.0 - (2.0 * i + 1.0) / (double)SH_NDIR;
		double rxy = sqrt((zz * zz < 1.0) ? (1.0 - zz * zz) : 0.0);
		double d[3], t, nrm[3];
		int face;
		d[0] = rxy * cos(golden * i);
		d[1] = rxy * sin(golden * i);
		d[2] = zz;
		if (ga_solid_entry(o, p, d, 0.0, best, &t, nrm, &face) && t < best)
			best = t;
	}
	return best;
}
