/* ga_trace.c — 早期反射 (鏡像法) と後期残響 (光線追跡)
 *
 * ── 鏡像法 (ga_images) ────────────────────────────────────────────
 * 室は直方体なので鏡像音源は閉形式で書ける。x 軸について整数 k の像は
 *   k 偶数 : x_k = x0 + k Lx + (xs - x0)
 *   k 奇数 : x_k = x0 + k Lx + (Lx - (xs - x0))
 * (k = 0 が実音源、k = ±1 が 1 回反射)。反射回数は |k|、内訳は
 * k > 0 なら x+ 壁を ceil(k/2) 回・x- 壁を floor(k/2) 回 (k < 0 は逆)。
 * y, z も同様で、次数 = |a| + |b| + |c| <= order の像だけを使う。
 *
 * 振幅 (バンド別) :
 *   A_b = 1/(4 pi r) * prod_w R_w,b^{n_w} * 10^(-a_b r / 20)
 *   - 1/(4 pi r) : 自由音場の直接音を 1/(4 pi r) に正規化する規約
 *     (FDTD 側と合成するときの共通規約)。像でも経路長 r で同じ形になる。
 *   - R_w,b = sqrt(1 - alpha_w,b) : エネルギー反射率 1 - alpha に対応する
 *     圧力反射係数。
 *   - 10^(-a_b r/20) : ISO 9613-1 の空気吸収 (音圧レベル a_b [dB/m])。
 *   到達時刻は t = r/c。t = 0 は音源発火時刻なので遅延は入れない。
 *
 * 可視性判定 : 直方体の室では像は必ず「幾何学的に到達可能」だが、室内の
 * 剛体障害物に遮られる経路は棄却しなければならない。展開空間の直線を
 * 各セル境界で分割し、セルごとに室内へ折り返して実際の経路 (折れ線) を
 * 復元し、その全区間を障害物 AABB と交差判定する (近似なしの厳密判定)。
 *
 * ── 光線追跡 (ga_rays) ────────────────────────────────────────────
 * 方向は球面フィボナッチ格子 (決定的な準一様分布 — 乱数は使わない)。
 *   z_i = 1 - (2i+1)/N,  phi_i = i * pi (3 - sqrt(5))
 * 各レイはバンド別エネルギー w_b を運び、壁で (1 - alpha_b) を、
 * 空気中で exp(-m_b s) を掛ける。障害物は剛体 (alpha = 0) として鏡面反射。
 *
 * 検出は半径 R の受音球の通過本数で行う。自由音場で N 本のレイのうち
 * 球を通過するのは N R^2/(4 r^2) 本、1 本あたりのエネルギーは 1/N なので
 * 通過エネルギー和 = I * pi R^2 (I = 1/(4 pi r^2))。本ソルバーの振幅規約
 * A = 1/(4 pi r) では A^2 = I/(4 pi) なので、通過 1 本あたりの寄与は
 *   dE = (1/N) * w_b / (4 pi^2 R^2)
 * となる (自由音場で dE の総和が厳密に 1/(4 pi r)^2 になることで検算できる)。
 *
 * 二重計上の回避 : 鏡像法が厳密に受け持つ経路 — 反射回数が order 以下で、
 * 障害物に当たらず、一度も拡散反射していない鏡面経路 — はレイ側で検出しない。
 * 障害物に当たったレイ・拡散したレイは鏡像法に含まれないので回数によらず
 * 検出する。散乱係数 s のぶんは鏡像の振幅から (1-s)^(n/2) で抜いてあるので、
 * 「鏡面成分 = 鏡像法、拡散成分 = レイ」で過不足なく 1 回ずつ数えられる。
 *
 * 拡散反射 (Lambert) : 散乱係数 s は .ofdx acoustic.ga.scattering (既定 0.1)。
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ga.h"

#define GA_EPS        1e-9
#define GA_MAXBOUNCE  200000   /* 退化した幾何での無限ループ止め */
#define GA_PROG_RAYS  45       /* "progress k/50" のうちレイ追跡に割く分 */
#define GA_PROG_TOTAL 50

/* ── 線分 vs AABB (交差すれば 1) ────────────────────────────────
 * 接触 (面上をかすめる) を遮蔽と誤判定しないよう箱をごく僅か縮める。 */
static int seg_hits_box(const double p0[3], const double p1[3],
                        const double lo[3], const double hi[3])
{
	double tmin = 0.0, tmax = 1.0;
	int i;
	for (i = 0; i < 3; i++) {
		double a = lo[i] + GA_EPS, b = hi[i] - GA_EPS;
		double d = p1[i] - p0[i];
		if (b < a) { a = b = 0.5 * (lo[i] + hi[i]); }
		if (fabs(d) < 1e-15) {
			if (p0[i] < a || p0[i] > b) return 0;
		}
		else {
			double t1 = (a - p0[i]) / d;
			double t2 = (b - p0[i]) / d;
			if (t1 > t2) { double t = t1; t1 = t2; t2 = t; }
			if (t1 > tmin) tmin = t1;
			if (t2 < tmax) tmax = t2;
			if (tmin > tmax) return 0;
		}
	}
	return 1;
}

static int any_obstacle_hits(const ga_t *g, const double p0[3], const double p1[3])
{
	int n;
	for (n = 0; n < g->ngeom; n++) {
		if (!g->geom[n].ok) continue;
		if (seg_hits_box(p0, p1, g->geom[n].lo, g->geom[n].hi)) return 1;
	}
	return 0;
}

/* ── 鏡像法 ─────────────────────────────────────────────────────── */

/* 軸 i について k 回折り返した像の座標 */
static double image_coord(double s, double lo, double len, int k)
{
	double sl = s - lo;
	double v = ((k % 2) == 0) ? (k * len + sl) : (k * len + (len - sl));
	return lo + v;
}

/* 展開空間の座標 X を室内へ折り返す (セル番号 k を与える) */
static double fold_coord(double X, double lo, double len, int k)
{
	double frac = (X - lo) - (double)k * len;
	return ((k % 2) == 0) ? (lo + frac) : (lo + len - frac);
}

/* 挿入ソート (決定的。要素数は最大 11) */
static void sort_asc(double *t, int n)
{
	int i, j;
	for (i = 1; i < n; i++) {
		double v = t[i];
		for (j = i - 1; j >= 0 && t[j] > v; j--) t[j + 1] = t[j];
		t[j + 1] = v;
	}
}

/* 像音源 A (展開座標) から受音点 B (室内) への実経路が障害物に遮られるか。
 * 展開直線をセル境界で分割し、区間ごとに室内へ折り返して判定する。 */
static int image_blocked(const ga_t *g, const double A[3], const double B[3])
{
	double lo[3], len[3];
	double ts[24];   /* 交点は最大 3 軸 x order(<=3) + 端点 2 = 11 */
	int nt = 0, i, s;

	if (g->ngeom <= 0) return 0;
	lo[0] = g->x0; lo[1] = g->y0; lo[2] = g->z0;
	len[0] = g->x1 - g->x0; len[1] = g->y1 - g->y0; len[2] = g->z1 - g->z0;

	ts[nt++] = 0.0;
	ts[nt++] = 1.0;
	for (i = 0; i < 3; i++) {
		double d = B[i] - A[i];
		int m, m0, m1;
		if (fabs(d) < 1e-15) continue;
		m0 = (int)floor((A[i] - lo[i]) / len[i]);
		m1 = (int)floor((B[i] - lo[i]) / len[i]);
		if (m0 > m1) { int tmp = m0; m0 = m1; m1 = tmp; }
		for (m = m0; m <= m1 + 1; m++) {
			double plane = lo[i] + (double)m * len[i];
			double t = (plane - A[i]) / d;
			if (t > 1e-12 && t < 1.0 - 1e-12 && nt < 24) ts[nt++] = t;
		}
	}
	sort_asc(ts, nt);

	for (s = 0; s + 1 < nt; s++) {
		double t0 = ts[s], t1 = ts[s + 1], tm = 0.5 * (t0 + t1);
		double p0[3], p1[3];
		int k[3];
		if (t1 - t0 < 1e-12) continue;
		for (i = 0; i < 3; i++) {
			double xm = A[i] + tm * (B[i] - A[i]);
			k[i] = (int)floor((xm - lo[i]) / len[i]);
		}
		for (i = 0; i < 3; i++) {
			p0[i] = fold_coord(A[i] + t0 * (B[i] - A[i]), lo[i], len[i], k[i]);
			p1[i] = fold_coord(A[i] + t1 * (B[i] - A[i]), lo[i], len[i], k[i]);
		}
		if (any_obstacle_hits(g, p0, p1)) return 1;
	}
	return 0;
}

int ga_images(ga_t *g, int ri, double *band)
{
	ga_recv_t *rv = &g->recv[ri];
	double L[3], O[3], S[3], B[3];
	double amp[GA_NBAND];
	int a, b, c, w, nb;
	int order = g->order;

	O[0] = g->x0; O[1] = g->y0; O[2] = g->z0;
	L[0] = g->x1 - g->x0; L[1] = g->y1 - g->y0; L[2] = g->z1 - g->z0;
	S[0] = g->srcx; S[1] = g->srcy; S[2] = g->srcz;
	B[0] = rv->x; B[1] = rv->y; B[2] = rv->z;

	rv->nimage = 0;
	rv->nblocked = 0;

	for (a = -order; a <= order; a++) {
		int na = (a < 0) ? -a : a;
		for (b = -order; b <= order; b++) {
			int nbb = (b < 0) ? -b : b;
			if (na + nbb > order) continue;
			for (c = -order; c <= order; c++) {
				int nc = (c < 0) ? -c : c;
				int cnt[GA_NWALL];
				double A[3], dx, dy, dz, r, g0;
				if (na + nbb + nc > order) continue;

				A[0] = image_coord(S[0], O[0], L[0], a);
				A[1] = image_coord(S[1], O[1], L[1], b);
				A[2] = image_coord(S[2], O[2], L[2], c);
				dx = A[0] - B[0]; dy = A[1] - B[1]; dz = A[2] - B[2];
				r = sqrt(dx * dx + dy * dy + dz * dz);
				if (r < 1e-9) continue;
				if (r / GA_C0 >= g->duration) continue;

				/* 反射した壁の回数 */
				for (w = 0; w < GA_NWALL; w++) cnt[w] = 0;
				if (a > 0) { cnt[GA_XP] = (a + 1) / 2; cnt[GA_XM] = a / 2; }
				else if (a < 0) { cnt[GA_XM] = (na + 1) / 2; cnt[GA_XP] = na / 2; }
				if (b > 0) { cnt[GA_YP] = (b + 1) / 2; cnt[GA_YM] = b / 2; }
				else if (b < 0) { cnt[GA_YM] = (nbb + 1) / 2; cnt[GA_YP] = nbb / 2; }
				if (c > 0) { cnt[GA_ZP] = (c + 1) / 2; cnt[GA_ZM] = c / 2; }
				else if (c < 0) { cnt[GA_ZM] = (nc + 1) / 2; cnt[GA_ZP] = nc / 2; }

				/* 散乱係数 s のぶんは鏡面成分から抜ける (エネルギーで (1-s)^n、
				 * 振幅で (1-s)^(n/2))。抜けた拡散成分は光線追跡側が
				 * 「1 回でも拡散したレイ」として受け持つので二重計上しない。 */
				g0 = 1.0 / (4.0 * GA_PI * r)
				   * pow(1.0 - g->scatter, 0.5 * (na + nbb + nc));
				for (nb = 0; nb < GA_NBAND; nb++) {
					double v = g0 * pow(10.0, -g->air_db_m[nb] * r / 20.0);
					for (w = 0; w < GA_NWALL; w++) {
						int q;
						for (q = 0; q < cnt[w]; q++) v *= g->refl[w][nb];
					}
					amp[nb] = v;
				}
				/* 全バンドで消えている像は置かない (alpha = 1 の壁など) */
				{
					int alive = 0;
					for (nb = 0; nb < GA_NBAND; nb++)
						if (amp[nb] > 1e-300) { alive = 1; break; }
					if (!alive) continue;
				}
				if (image_blocked(g, A, B)) {
					rv->nblocked++;
					continue;
				}
				ga_deposit(g, band, r / GA_C0, amp);
				rv->nimage++;
			}
		}
	}
	ga_log(g, "point #%d: image sources up to order %d -> %d visible, "
	       "%d blocked by obstacles", ri + 1, order, rv->nimage, rv->nblocked);
	return 0;
}

/* ── 光線追跡 ───────────────────────────────────────────────────── */

/* ── 拡散反射 (Lambert) ─────────────────────────────────────────
 * 実在の壁は完全な鏡ではなく、凹凸・什器で入射エネルギーの一部を拡散する。
 * 拡散の割合 (散乱係数 s) は .ofdx acoustic.ga.scattering で与える。
 * これは見た目の細部ではなく**残響時間そのもの**に効く : 鏡面反射だけの
 * 直方体ではレイの方向余弦が保存されるため反射回数が方向ごとに固定され、
 * 場が混合しない。混合しない場の減衰は指数関数の重ね合わせになり、対数軸で
 * 凸に折れて Eyring 式 (拡散場の仮定) から late で数 % 〜 10 % ずれる。
 * 拡散反射を入れると自由行程が更新過程になり、拡散場に収束する。
 *
 * 方向の抽選は**乱数生成器を使わない**。決定的カウンタ g->qidx に固定の
 * 32 bit 整数ハッシュ (splitmix 系のアバランチ) を掛けて [0,1) を作るだけで、
 * 種も状態も外部に依存しない — 再実行・スレッド数・OS によらず同じ列になる
 * (レイのループは直列で、qidx の進み方も入力だけで決まる)。
 * cos 重み付き半球サンプリング :
 *   cos(theta) = sqrt(u1),  phi = 2 pi u2
 * これは Lambert の余弦則そのもの (放射輝度が方向によらない)。
 *
 * 準乱数列 (Halton) を使わない理由 : 基数 2,3 の基数逆列を連番で引くと、
 * 直方体の反射列と共鳴して平均自由行程が 4V/S から系統的にずれる
 * (10 m 立方体で 6.99 m vs 理論 6.6667 m = +4.8%、残響時間がそのまま
 * 数 % 狂う)。ハッシュ列では 6.666 m と理論値に一致する。 */
static double hash01(uint32_t x)
{
	x ^= x >> 16; x *= 0x7feb352du;
	x ^= x >> 15; x *= 0x846ca68bu;
	x ^= x >> 16;
	return (double)(x >> 8) * (1.0 / 16777216.0);   /* [0,1) — 24 bit */
}

/* 法線 n (軸番号 ax、向き sgn = +1/-1) の半球へ cos 重みで散乱させる */
static void lambert_dir(double d[3], int ax, double sgn, double u1, double u2)
{
	int t1 = (ax + 1) % 3, t2 = (ax + 2) % 3;
	double ct = sqrt(u1);                    /* cos(theta) = sqrt(u1) */
	double st = sqrt((u1 < 1.0) ? (1.0 - u1) : 0.0);
	double ph = 2.0 * GA_PI * u2;
	d[ax] = sgn * ct;
	d[t1] = st * cos(ph);
	d[t2] = st * sin(ph);
}

/* 室内から見た壁までの距離 (最短) と壁番号 */
static double wall_exit(const ga_t *g, const double p[3], const double d[3],
                        int *wall, int *axis)
{
	double lo[3], hi[3], best = 1e300;
	int i, bw = GA_XP, ba = 0;
	lo[0] = g->x0; lo[1] = g->y0; lo[2] = g->z0;
	hi[0] = g->x1; hi[1] = g->y1; hi[2] = g->z1;
	for (i = 0; i < 3; i++) {
		double t;
		if (d[i] > 1e-15) {
			t = (hi[i] - p[i]) / d[i];
			if (t < best) { best = t; ba = i; bw = GA_XP + 2 * i; }
		}
		else if (d[i] < -1e-15) {
			t = (lo[i] - p[i]) / d[i];
			if (t < best) { best = t; ba = i; bw = GA_XM + 2 * i; }
		}
	}
	*wall = bw;
	*axis = ba;
	return best;
}

/* 障害物への最短入射 (無ければ負を返す) */
static double obstacle_entry(const ga_t *g, const double p[3], const double d[3],
                             double tlimit, int *axis)
{
	double best = -1.0;
	int n, ba = -1;
	for (n = 0; n < g->ngeom; n++) {
		const ga_geom_t *o = &g->geom[n];
		double tnear = -1e300, tfar = 1e300;
		int i, na = -1, hit = 1;
		if (!o->ok) continue;
		for (i = 0; i < 3; i++) {
			if (fabs(d[i]) < 1e-15) {
				if (p[i] < o->lo[i] || p[i] > o->hi[i]) { hit = 0; break; }
			}
			else {
				double t1 = (o->lo[i] - p[i]) / d[i];
				double t2 = (o->hi[i] - p[i]) / d[i];
				if (t1 > t2) { double t = t1; t1 = t2; t2 = t; }
				if (t1 > tnear) { tnear = t1; na = i; }
				if (t2 < tfar) tfar = t2;
				if (tnear > tfar) { hit = 0; break; }
			}
		}
		if (!hit || na < 0) continue;
		if (tnear <= 1e-9 || tnear >= tfar) continue;
		if (tnear >= tlimit) continue;
		if (best < 0.0 || tnear < best) { best = tnear; ba = na; }
	}
	*axis = ba;
	return best;
}

int ga_rays(ga_t *g)
{
	const double golden = GA_PI * (3.0 - sqrt(5.0));
	const double smax = GA_C0 * g->duration;
	const double einit = 1.0 / (double)g->nrays;
	double lo[3], hi[3];
	int i, prog = 0;

	lo[0] = g->x0; lo[1] = g->y0; lo[2] = g->z0;
	hi[0] = g->x1; hi[1] = g->y1; hi[2] = g->z1;

	printf("solve: %d rays, image order %d, %d samples at %d Hz, %d receivers "
	       "(進捗は %d 分割)\n",
	       g->nrays, g->order, g->nsamples, GA_FS, g->nrecv, GA_PROG_TOTAL);
	fflush(stdout);

	for (i = 0; i < g->nrays; i++) {
		double p[3], d[3], w[GA_NBAND];
		double z, rxy, phi, s = 0.0;
		int nref = 0, hitobs = 0, scattered = 0, bounce, b;

		z   = 1.0 - (2.0 * i + 1.0) / (double)g->nrays;
		rxy = sqrt((z * z < 1.0) ? (1.0 - z * z) : 0.0);
		phi = golden * (double)i;
		d[0] = rxy * cos(phi);
		d[1] = rxy * sin(phi);
		d[2] = z;
		p[0] = g->srcx; p[1] = g->srcy; p[2] = g->srcz;
		for (b = 0; b < GA_NBAND; b++) w[b] = einit;

		for (bounce = 0; bounce < GA_MAXBOUNCE; bounce++) {
			int wall, waxis, oaxis;
			double tw, tobs, thit;
			int is_obs = 0, last = 0, r;

			tw = wall_exit(g, p, d, &wall, &waxis);
			if (!(tw > 0.0) || tw > 1e299) break;   /* 退化 : 打ち切る */
			tobs = obstacle_entry(g, p, d, tw, &oaxis);
			thit = tw;
			if (tobs > 0.0) { thit = tobs; is_obs = 1; }
			if (s + thit >= smax) { thit = smax - s; last = 1; }
			if (thit <= 0.0) break;

			/* 受音球の通過。鏡像法が厳密に受け持つ経路 — 次数 order 以下で
			 * 障害物にも当たらず一度も拡散していない鏡面経路 — は数えない
			 * (二重計上の回避)。拡散した時点で鏡像法の対象外になる。 */
			if (!(nref <= g->order && !hitobs && !scattered)) {
				for (r = 0; r < g->nrecv; r++) {
					const ga_recv_t *rv = &g->recv[r];
					double u[3], tc, dist2, rad = rv->radius;
					u[0] = rv->x - p[0]; u[1] = rv->y - p[1]; u[2] = rv->z - p[2];
					tc = u[0] * d[0] + u[1] * d[1] + u[2] * d[2];
					if (tc <= 0.0 || tc >= thit) continue;
					dist2 = (u[0] - tc * d[0]) * (u[0] - tc * d[0])
					      + (u[1] - tc * d[1]) * (u[1] - tc * d[1])
					      + (u[2] - tc * d[2]) * (u[2] - tc * d[2]);
					if (dist2 >= rad * rad) continue;
					{
						double sd = s + tc;
						int bin = (int)((sd / GA_C0) / GA_BIN_S);
						double *e;
						if (bin < 0 || bin >= g->nbins) continue;
						e = g->echo + (((size_t)r * g->nbins) + bin) * GA_NBAND;
						for (b = 0; b < GA_NBAND; b++)
							e[b] += w[b] * exp(-g->air_e[b] * tc)
							      / (4.0 * GA_PI * GA_PI * rad * rad);
						g->recv[r].nhit++;
						g->nray_detect++;
					}
				}
			}

			/* 前進 + 空気吸収 */
			for (b = 0; b < GA_NBAND; b++) w[b] *= exp(-g->air_e[b] * thit);
			p[0] += thit * d[0];
			p[1] += thit * d[1];
			p[2] += thit * d[2];
			s += thit;
			g->nbounce++;
			if (last) break;

			/* 反射 : 割合 s で拡散 (Lambert)、残りは鏡面。抽選は Halton 列 */
			{
				int    ax  = is_obs ? oaxis : waxis;
				double sgn = (d[ax] > 0.0) ? -1.0 : 1.0;   /* 面の内向き法線 */
				uint32_t q = (uint32_t)g->qidx * 3u;
				int    diffuse = 0;
				if (g->scatter > 0.0)
					diffuse = (hash01(q + 3u) < g->scatter);
				if (is_obs) {
					hitobs = 1;      /* 障害物は剛体 (alpha = 0) */
				}
				else {
					for (b = 0; b < GA_NBAND; b++) w[b] *= 1.0 - g->alpha[wall][b];
				}
				if (diffuse) {
					lambert_dir(d, ax, sgn, hash01(q + 1u), hash01(q + 2u));
					scattered = 1;
				}
				else {
					d[ax] = -d[ax];
				}
				g->qidx++;
				/* 面から僅かに離して自己交差を避ける */
				p[ax] += sgn * GA_EPS;
				if (!is_obs) {
					if (p[waxis] < lo[waxis]) p[waxis] = lo[waxis];
					if (p[waxis] > hi[waxis]) p[waxis] = hi[waxis];
				}
			}
			nref++;
			{
				double wmax = 0.0;
				for (b = 0; b < GA_NBAND; b++) if (w[b] > wmax) wmax = w[b];
				if (wmax <= GA_ENERGY_FLOOR * einit) break;
			}
		}

		{
			int k = (int)(((long long)(i + 1) * GA_PROG_RAYS) / g->nrays);
			while (prog < k) {
				++prog;
				printf("progress %d/%d\n", prog, GA_PROG_TOTAL);
				fflush(stdout);
			}
		}
	}
	while (prog < GA_PROG_RAYS) {
		++prog;
		printf("progress %d/%d\n", prog, GA_PROG_TOTAL);
		fflush(stdout);
	}
	ga_log(g, "ray tracing: %d rays, %ld reflections, %ld receiver-sphere "
	       "detections", g->nrays, g->nbounce, g->nray_detect);
	for (i = 0; i < g->nrecv; i++)
		ga_log(g, "point #%d: %ld ray detections (sphere r = %.4g m)",
		       i + 1, g->recv[i].nhit, g->recv[i].radius);
	return 0;
}
