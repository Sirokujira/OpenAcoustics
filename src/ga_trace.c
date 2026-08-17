/* ga_trace.c — 早期反射 (鏡像法) と後期残響 (光線追跡)
 *
 * ── 鏡像法 (ga_images) ────────────────────────────────────────────
 * 反射面は「軸平行の有限矩形」の集合 surf[] : 室の 6 面 + 障害物 (AABB)
 * 1 個あたり 6 面。この集合に対する**一般化された鏡像法**なので、室壁だけの
 * 経路も、障害物を含む経路 (壁→衝立→受音点 など) も同じ手順で扱える。
 * 面が軸平行なので鏡映は 1 座標の反転で済む。
 *
 *   1. 音源から始めて、面 s の音場側にある像だけを s で鏡映して深さを進める
 *      (Borish の枝刈り : 面の裏側にある像は、その面では反射できない)。
 *   2. 各ノードで受音点から逆追跡し、反射点が**その面の矩形の内側**に
 *      落ちるかを確かめる (可視性判定の前半 = 有限面の判定)。
 *   3. 復元した折れ線の全区間を障害物 AABB と交差判定する
 *      (可視性判定の後半 = 遮蔽)。
 * 直方体だけの室ではこれは Allen–Berkley の閉形式と厳密に一致する
 * (次数 3 で像は 63 個)。
 *
 * 振幅 (バンド別) :
 *   A_b = 1/(4 pi r) * prod_k R_k,b * prod_k sqrt(1 - s_k) * 10^(-a_b r / 20)
 *   - 1/(4 pi r) : 自由音場の直接音を 1/(4 pi r) に正規化する規約
 *     (FDTD 側と合成するときの共通規約)。像でも経路長 r で同じ形になる。
 *   - R_k,b : 反射 k の圧力反射係数。角度依存吸音が off なら sqrt(1-alpha)、
 *     on なら局所反応の R(theta) = (zeta cos - 1)/(zeta cos + 1)。
 *     展開空間では経路が 1 本の直線なので、軸 a の面での cos(theta) は
 *     どの反射でも |u_a| (u = 展開直線の単位方向) で与えられる。
 *   - sqrt(1 - s_k,b) : 面ごと (バンド別) の散乱係数で鏡面成分から抜ける分。
 *   - 10^(-a_b r/20) : ISO 9613-1 の空気吸収 (音圧レベル a_b [dB/m])。
 *   到達時刻は t = r/c。t = 0 は音源発火時刻なので遅延は入れない。
 *
 * ── 光線追跡 (ga_rays) ────────────────────────────────────────────
 * 方向は球面フィボナッチ格子 (決定的な準一様分布 — 乱数は使わない)。
 *   z_i = 1 - (2i+1)/N,  phi_i = i * pi (3 - sqrt(5))
 * 各レイはバンド別エネルギー w_b を運び、面で R^2 を、空気中で exp(-m_b s) を
 * 掛ける。反射する面は鏡像法と同じ surf[] なので、吸音・散乱係数・角度依存の
 * 扱いが両者で完全に一致する。
 *
 * 検出は半径 R の受音球の通過本数で行う。自由音場で N 本のレイのうち
 * 球を通過するのは N R^2/(4 r^2) 本、1 本あたりのエネルギーは 1/N なので
 * 通過エネルギー和 = I * pi R^2 (I = 1/(4 pi r^2))。本ソルバーの振幅規約
 * A = 1/(4 pi r) では A^2 = I/(4 pi) なので、通過 1 本あたりの寄与は
 *   dE = (1/N) * w_b / (4 pi^2 R^2)
 * となる (自由音場で dE の総和が厳密に 1/(4 pi r)^2 になることで検算できる)。
 *
 * 二重計上の回避 : 鏡像法が厳密に受け持つ経路 — 反射回数が order 以下で
 * 一度も拡散反射していない鏡面経路 — はレイ側で検出しない。障害物の面も
 * 鏡像法の対象になったので、条件は「次数」と「拡散したか」だけで決まる。
 * 散乱係数 s_b のぶんは鏡像の振幅から面ごと・バンドごとに sqrt(1-s_b) で
 * 抜いてあり、レイ側は重み付き抽選 (下記) で期待値どおり s_b / (1-s_b) を
 * 運ぶので、「鏡面成分 = 鏡像法、拡散成分 = レイ」で過不足なく 1 回ずつ
 * 数えられる。
 *
 * 拡散反射 (Lambert) : 散乱係数は面ごと・バンド別 (.ofdx の absorption 行の
 * scattering、省略時は acoustic.ga.scattering — どちらも数値か配列)。
 * 抽選は基準確率 sref (バンド平均) で 1 回、バンド別はエネルギーの
 * 重み付けで実現する (ga.h と反射処理のコメント参照)。
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ga.h"

#define GA_EPS        1e-9
#define GA_SEG_EPS    1e-6     /* 遮蔽判定で線分の両端を縮める割合 */
#define GA_MAXBOUNCE  200000   /* 退化した幾何での無限ループ止め */
#define GA_PROG_RAYS  45       /* "progress k/50" のうちレイ追跡に割く分 */
#define GA_PROG_TOTAL 50

/* ── 線分 vs AABB (交差すれば 1) ────────────────────────────────
 * 接触 (面上をかすめる・面で反射する) を遮蔽と誤判定しないよう、箱をごく
 * 僅かに縮め、線分の両端も縮めてから判定する。 */
static int seg_hits_box(const double p0[3], const double p1[3],
                        const double lo[3], const double hi[3])
{
	double tmin = GA_SEG_EPS, tmax = 1.0 - GA_SEG_EPS;
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

/* ── 反射面の基本演算 ───────────────────────────────────────────── */

/* 面 s による鏡映 (軸平行なので 1 座標の反転) */
static void surf_mirror(const ga_surf_t *s, const double p[3], double q[3])
{
	q[0] = p[0]; q[1] = p[1]; q[2] = p[2];
	q[s->axis] = 2.0 * s->coord - p[s->axis];
}

/* 面 s の音場側から測った符号付き距離 (> 0 なら反射できる側にいる) */
static double surf_front(const ga_surf_t *s, const double p[3])
{
	return (p[s->axis] - s->coord) * s->nrm;
}

/* 入射角 cos(theta) での圧力反射係数 (符号つき — 位相反転も表す) */
static double surf_refl(const ga_t *g, const ga_surf_t *s, int band, double ct)
{
	double zc;
	if (!g->angle_dep) return s->refl[band];
	zc = s->zeta[band] * ct;
	return (zc - 1.0) / (zc + 1.0);
}

/* ── 鏡像法 ─────────────────────────────────────────────────────── */

typedef struct {
	ga_t   *g;
	int     ri;
	double *band;
	double  recv[3];
	double  src[3];
	double  sgain;     /* 音源のゲイン (振幅に掛かる — ADR-0010 Decision 7) */
	double  sdelay;    /* 音源の発火時刻 [s] (到達は sdelay + r/c) */
	long    nodes;
	int     truncated;
	double *seen;      /* 重複除去 : 1 像あたり x, y, z, depth の 4 要素 */
	int     nseen, cseen;
	int     oom;
} imgctx_t;

/* 同じ位置・同じ次数の像に複数の面の並びから到達することがある。
 * 音源と受音点が室の対称面上に厳密に乗っているときの「稜線に落ちる経路」が
 * 典型で、例えば (y+ のあと z+) と (z+ のあと y+) が同じ 1 本の経路を表す。
 * 経路としては 1 本なので、最初に見つけたものだけを採用する。
 * (この重複除去が無いと、対称配置で像が二重計上されるか、逆に稜線判定で
 *  両方とも落ちて像が欠ける。直方体では像の位置と次数が経路を一意に決める。) */
static int img_seen(imgctx_t *c, const double I[3], int depth)
{
	int i;
	for (i = 0; i < c->nseen; i++) {
		const double *e = c->seen + 4 * i;
		if ((int)e[3] == depth &&
		    fabs(e[0] - I[0]) < 1e-9 &&
		    fabs(e[1] - I[1]) < 1e-9 &&
		    fabs(e[2] - I[2]) < 1e-9)
			return 1;
	}
	if (c->nseen >= c->cseen) {
		int cap = c->cseen ? c->cseen * 2 : 256;
		double *p = (double *)realloc(c->seen, (size_t)cap * 4 * sizeof(double));
		if (!p) { c->oom = 1; return 1; }   /* 確保できなければ以降は数えない */
		c->seen = p;
		c->cseen = cap;
	}
	c->seen[4 * c->nseen + 0] = I[0];
	c->seen[4 * c->nseen + 1] = I[1];
	c->seen[4 * c->nseen + 2] = I[2];
	c->seen[4 * c->nseen + 3] = (double)depth;
	c->nseen++;
	return 0;
}

/* 深さ depth の像 (chain[depth]) が受音点から見えるかを判定し、見えるなら
 * band[] へ置く。chain[0] = 実音源、chain[k+1] = chain[k] を surf[seq[k]] で
 * 鏡映したもの。 */
static void img_visit(imgctx_t *c, double chain[][3], const int *seq, int depth)
{
	ga_t *g = c->g;
	ga_recv_t *rv = &g->recv[c->ri];
	const double *I = chain[depth];
	double u[3], r, amp[GA_NBAND], gain, ct[GA_ORDER_MAX];
	double pts[GA_ORDER_MAX + 2][3];
	int i, k, b;

	u[0] = c->recv[0] - I[0];
	u[1] = c->recv[1] - I[1];
	u[2] = c->recv[2] - I[2];
	r = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
	if (r < 1e-9) return;
	if (c->sdelay + r / GA_C0 >= g->duration) return;
	u[0] /= r; u[1] /= r; u[2] /= r;

	/* 逆追跡 : 受音点 -> 最後の反射面 -> ... -> 音源。反射点が有限矩形の
	 * 内側に落ちなければ、その像は幾何学的に存在しない。 */
	pts[0][0] = c->recv[0]; pts[0][1] = c->recv[1]; pts[0][2] = c->recv[2];
	for (i = 1; i <= depth; i++) {
		const ga_surf_t *s = &g->surf[seq[depth - i]];
		const double *P = pts[i - 1];
		const double *T = chain[depth - i + 1];
		double den = T[s->axis] - P[s->axis];
		double t, X[3];
		int a = s->axis, uu, vv;
		if (fabs(den) < 1e-15) return;
		t = (s->coord - P[s->axis]) / den;
		/* t = 0 は「前の反射点と同じ場所」= 稜線に落ちる縮退経路。物理的には
		 * 角での 2 回反射で、像としては実在するので許容する (重複は img_seen
		 * が除く)。t >= 1 は反射点が像そのものになる = 実在しない。 */
		if (!(t > -1e-9 && t < 1.0 - 1e-12)) return;
		if (t < 0.0) t = 0.0;
		X[0] = P[0] + t * (T[0] - P[0]);
		X[1] = P[1] + t * (T[1] - P[1]);
		X[2] = P[2] + t * (T[2] - P[2]);
		X[a] = s->coord;
		uu = (a + 1) % 3; vv = (a + 2) % 3;
		if (X[uu] < s->lo[0] || X[uu] > s->hi[0] ||
		    X[vv] < s->lo[1] || X[vv] > s->hi[1]) return;
		pts[i][0] = X[0]; pts[i][1] = X[1]; pts[i][2] = X[2];
	}
	pts[depth + 1][0] = c->src[0];
	pts[depth + 1][1] = c->src[1];
	pts[depth + 1][2] = c->src[2];

	if (img_seen(c, I, depth)) return;   /* 別の面の並びで到達済みの同一経路 */

	/* 遮蔽 : 復元した折れ線の全区間 */
	for (i = 0; i <= depth; i++) {
		if (any_obstacle_hits(g, pts[i], pts[i + 1])) {
			rv->nblocked++;
			return;
		}
	}

	/* 振幅 : 音源ゲイン・距離減衰・面ごとの反射係数と散乱・空気吸収。
	 * 散乱の sqrt(1-s_b) はバンド一様な面では gain へ (従来とビット等価)、
	 * バンド別の面ではバンドループ内で掛ける。 */
	gain = c->sgain / (4.0 * GA_PI * r);
	for (k = 0; k < depth; k++) {
		const ga_surf_t *s = &g->surf[seq[k]];
		ct[k] = fabs(u[s->axis]);
		if (s->suni) gain *= sqrt(1.0 - s->scatter[0]);
	}
	for (b = 0; b < GA_NBAND; b++) {
		double v = gain * pow(10.0, -g->air_db_m[b] * r / 20.0);
		for (k = 0; k < depth; k++) {
			const ga_surf_t *s = &g->surf[seq[k]];
			v *= surf_refl(g, s, b, ct[k]);
			if (!s->suni) v *= sqrt(1.0 - s->scatter[b]);
		}
		amp[b] = v;
	}
	/* 全バンドで消えている像は置かない (alpha = 1 の壁など) */
	for (b = 0; b < GA_NBAND; b++)
		if (fabs(amp[b]) > 1e-300) break;
	if (b >= GA_NBAND) return;
	ga_deposit(g, c->band, c->sdelay + r / GA_C0, amp);
	rv->nimage++;
}

static void img_recurse(imgctx_t *c, double chain[][3], int *seq, int depth)
{
	ga_t *g = c->g;
	int si;

	img_visit(c, chain, seq, depth);
	if (depth >= g->order) return;
	for (si = 0; si < g->nsurf; si++) {
		const ga_surf_t *s = &g->surf[si];
		if (depth > 0 && seq[depth - 1] == si) continue;
		if (surf_front(s, chain[depth]) <= 1e-12) continue;
		if (++c->nodes > GA_IMAGE_NODE_MAX) { c->truncated = 1; return; }
		surf_mirror(s, chain[depth], chain[depth + 1]);
		seq[depth] = si;
		img_recurse(c, chain, seq, depth + 1);
		if (c->truncated) return;
	}
}

int ga_images(ga_t *g, int ri, int si, double *band)
{
	ga_recv_t *rv = &g->recv[ri];
	imgctx_t c;
	double chain[GA_ORDER_MAX + 1][3];
	int seq[GA_ORDER_MAX];

	memset(&c, 0, sizeof(c));
	memset(seq, 0, sizeof(seq));
	c.g = g;
	c.ri = ri;
	c.band = band;
	c.recv[0] = rv->x; c.recv[1] = rv->y; c.recv[2] = rv->z;
	c.src[0] = g->feedpos[3 * si + 0];
	c.src[1] = g->feedpos[3 * si + 1];
	c.src[2] = g->feedpos[3 * si + 2];
	c.sgain  = g->srcgain[si];
	c.sdelay = g->srcdelay[si];
	chain[0][0] = c.src[0]; chain[0][1] = c.src[1]; chain[0][2] = c.src[2];

	/* nimage/nblocked は加算のみ (リセットとまとめログは ga_synth —
	 * 複数音源では全音源分の合計が受音点の統計になる) */
	img_recurse(&c, chain, seq, 0);
	free(c.seen);

	if (c.oom)
		ga_log(g, "warning: point #%d source #%d: out of memory while "
		       "de-duplicating image sources — some reflections may be missing",
		       ri + 1, si + 1);
	if (c.truncated)
		ga_log(g, "warning: point #%d source #%d: the image-source search hit "
		       "the node limit (%d) — reflections beyond that were not "
		       "enumerated. Lower acoustic.ga.image_order or use fewer "
		       "geometry entries (%d reflecting surfaces).",
		       ri + 1, si + 1, GA_IMAGE_NODE_MAX, g->nsurf);
	return 0;
}

/* ── 拡散反射 (Lambert) ─────────────────────────────────────────
 * 実在の壁は完全な鏡ではなく、凹凸・什器で入射エネルギーの一部を拡散する。
 * 拡散の割合 (散乱係数 s) は面ごとに与える。これは見た目の細部ではなく
 * **残響時間そのもの**に効く : 鏡面反射だけの直方体ではレイの方向余弦が
 * 保存されるため反射回数が方向ごとに固定され、場が混合しない。
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
 * 数 % 狂う)。ハッシュ列では 6.655 m と理論値に一致する。 */
static double hash01(uint32_t x)
{
	x ^= x >> 16; x *= 0x7feb352du;
	x ^= x >> 15; x *= 0x846ca68bu;
	x ^= x >> 16;
	return (double)(x >> 8) * (1.0 / 16777216.0);   /* [0,1) — 24 bit */
}

/* 法線 (軸 ax、向き sgn) の半球へ cos 重みで散乱させる */
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

/* 室内から見た壁までの距離 (最短) と surf[] 添字 (室壁は 0..GA_NWALL-1) */
static double wall_exit(const ga_t *g, const double p[3], const double d[3],
                        int *sidx, int *axis)
{
	double lo[3], hi[3], best = 1e300;
	int i, bs = GA_XP, ba = 0;
	lo[0] = g->x0; lo[1] = g->y0; lo[2] = g->z0;
	hi[0] = g->x1; hi[1] = g->y1; hi[2] = g->z1;
	for (i = 0; i < 3; i++) {
		double t;
		if (d[i] > 1e-15) {
			t = (hi[i] - p[i]) / d[i];
			if (t < best) { best = t; ba = i; bs = 2 * i + 1; }
		}
		else if (d[i] < -1e-15) {
			t = (lo[i] - p[i]) / d[i];
			if (t < best) { best = t; ba = i; bs = 2 * i; }
		}
	}
	*sidx = bs;
	*axis = ba;
	return best;
}

/* 障害物への最短入射 (無ければ負を返す)。当たった面の surf[] 添字も返す */
static double obstacle_entry(const ga_t *g, const double p[3], const double d[3],
                             double tlimit, int *sidx, int *axis)
{
	double best = -1.0;
	int n, ba = -1, bs = -1;
	for (n = 0; n < g->ngeom; n++) {
		const ga_geom_t *o = &g->geom[n];
		double tnear = -1e300, tfar = 1e300;
		int i, na = -1, nside = 0, hit = 1;
		if (!o->ok || o->surf0 < 0) continue;
		for (i = 0; i < 3; i++) {
			if (fabs(d[i]) < 1e-15) {
				if (p[i] < o->lo[i] || p[i] > o->hi[i]) { hit = 0; break; }
			}
			else {
				double t1 = (o->lo[i] - p[i]) / d[i];
				double t2 = (o->hi[i] - p[i]) / d[i];
				if (t1 > t2) { double t = t1; t1 = t2; t2 = t; }
				if (t1 > tnear) {
					tnear = t1; na = i;
					nside = (d[i] > 0.0) ? 0 : 1;   /* lo 面 / hi 面 */
				}
				if (t2 < tfar) tfar = t2;
				if (tnear > tfar) { hit = 0; break; }
			}
		}
		if (!hit || na < 0) continue;
		if (tnear <= 1e-9 || tnear >= tfar) continue;
		if (tnear >= tlimit) continue;
		if (best < 0.0 || tnear < best) {
			best = tnear;
			ba = na;
			bs = o->surf0 + 2 * na + nside;
		}
	}
	*axis = ba;
	*sidx = bs;
	return best;
}

int ga_rays(ga_t *g)
{
	const double golden = GA_PI * (3.0 - sqrt(5.0));
	const long total = (long)g->nsrc * g->nrays;
	double lo[3], hi[3];
	int i, sidx, prog = 0;

	lo[0] = g->x0; lo[1] = g->y0; lo[2] = g->z0;
	hi[0] = g->x1; hi[1] = g->y1; hi[2] = g->z1;

	printf("solve: %d rays x %d source(s), image order %d over %d surfaces, "
	       "%d samples at %d Hz, %d receivers (進捗は %d 分割)\n",
	       g->nrays, g->nsrc, g->order, g->nsurf, g->nsamples, GA_FS, g->nrecv,
	       GA_PROG_TOTAL);
	fflush(stdout);

	for (sidx = 0; sidx < g->nsrc; sidx++) {
	/* 音源ごとのゲイン・遅延 (ADR-0010 Decision 7) :
	 *   - エネルギーは gain^2 重み (後期残響はエネルギー加算なので)
	 *   - レイの時間予算は duration - delay (検出時刻は delay + 経路長/c)
	 * 既定 (gain = 1, delay = 0) は従来とビット等価。 */
	const double sdelay = g->srcdelay[sidx];
	const double sgain2 = g->srcgain[sidx] * g->srcgain[sidx];
	const double smax   = GA_C0 * (g->duration - sdelay);
	const double einit  = sgain2 / (double)g->nrays;
	/* 音源ごとに決定的ハッシュ列を先頭から使う : 単一音源の既定動作と
	 * 一致し、multi_source の各音源も単独実行と同じ列になる */
	g->qidx = 0;
	if (smax <= 0.0 || einit <= 0.0) continue;   /* 窓外 / gain = 0 は寄与なし */
	for (i = 0; i < g->nrays; i++) {
		double p[3], d[3], w[GA_NBAND];
		double z, rxy, phi, s = 0.0;
		int nref = 0, scattered = 0, bounce, b;

		z   = 1.0 - (2.0 * i + 1.0) / (double)g->nrays;
		rxy = sqrt((z * z < 1.0) ? (1.0 - z * z) : 0.0);
		phi = golden * (double)i;
		d[0] = rxy * cos(phi);
		d[1] = rxy * sin(phi);
		d[2] = z;
		p[0] = g->feedpos[3 * sidx + 0];
		p[1] = g->feedpos[3 * sidx + 1];
		p[2] = g->feedpos[3 * sidx + 2];
		for (b = 0; b < GA_NBAND; b++) w[b] = einit;

		for (bounce = 0; bounce < GA_MAXBOUNCE; bounce++) {
			int wsidx, waxis, osidx, oaxis, sidx, hitaxis;
			double tw, tobs, thit;
			int last = 0, r;

			tw = wall_exit(g, p, d, &wsidx, &waxis);
			if (!(tw > 0.0) || tw > 1e299) break;   /* 退化 : 打ち切る */
			tobs = obstacle_entry(g, p, d, tw, &osidx, &oaxis);
			thit = tw;
			sidx = wsidx;
			hitaxis = waxis;
			if (tobs > 0.0 && osidx >= 0) {
				thit = tobs; sidx = osidx; hitaxis = oaxis;
			}
			if (s + thit >= smax) { thit = smax - s; last = 1; }
			if (thit <= 0.0) break;

			/* 受音球の通過。鏡像法が厳密に受け持つ経路 — 次数 order 以下で
			 * 一度も拡散していない鏡面経路 — は数えない (二重計上の回避)。
			 * 障害物の面も鏡像法の対象なので、条件に障害物は入らない。 */
			if (!(nref <= g->order && !scattered)) {
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
						int bin = (int)((sdelay + sd / GA_C0) / GA_BIN_S);
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

			/* 反射 : 面の吸音を掛け、割合 s で拡散 (Lambert)、残りは鏡面 */
			{
				const ga_surf_t *sf = &g->surf[sidx];
				int    ax  = hitaxis;
				double sgn = (d[ax] > 0.0) ? -1.0 : 1.0;   /* 音場側の法線 */
				double ct  = fabs(d[ax]);
				uint32_t q = (uint32_t)g->qidx * 3u;
				int diffuse = 0;

				for (b = 0; b < GA_NBAND; b++) {
					double rr = surf_refl(g, sf, b, ct);
					w[b] *= rr * rr;                     /* エネルギー反射率 */
				}
				/* 拡散/鏡面の抽選は基準確率 sref で 1 回 (1 本のレイが
				 * 全バンドを運ぶため)。バンド別の散乱係数は重み付けで実現 :
				 * 拡散枝 w_b *= s_b/sref、鏡面枝 w_b *= (1-s_b)/(1-sref)。
				 * 期待値は各バンドで厳密に s_b / (1-s_b) となり、鏡像側の
				 * sqrt(1-s_b) と過不足なく対応する (二重計上なし)。
				 * バンド一様 (suni) なら重みは厳密に 1 なので掛けない
				 * (従来動作とビット等価)。 */
				if (sf->sref > 0.0)
					diffuse = (hash01(q + 3u) < sf->sref);
				if (diffuse) {
					if (!sf->suni)
						for (b = 0; b < GA_NBAND; b++)
							w[b] *= sf->scatter[b] / sf->sref;
					lambert_dir(d, ax, sgn, hash01(q + 1u), hash01(q + 2u));
					scattered = 1;
				}
				else {
					if (!sf->suni)
						for (b = 0; b < GA_NBAND; b++)
							w[b] *= (1.0 - sf->scatter[b]) / (1.0 - sf->sref);
					d[ax] = -d[ax];
				}
				g->qidx++;
				/* 面から僅かに離して自己交差を避ける */
				p[ax] += sgn * GA_EPS;
				if (sf->wall >= 0) {
					if (p[ax] < lo[ax]) p[ax] = lo[ax];
					if (p[ax] > hi[ax]) p[ax] = hi[ax];
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
			int k = (int)(((long long)((long)sidx * g->nrays + i + 1)
			               * GA_PROG_RAYS) / total);
			while (prog < k) {
				++prog;
				printf("progress %d/%d\n", prog, GA_PROG_TOTAL);
				fflush(stdout);
			}
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
