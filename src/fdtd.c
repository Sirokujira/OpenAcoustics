/* fdtd.c — 3 次元線形音響スタガード格子 FDTD カーネル (v1)
 *
 * 場の配置 : p はセル中心 (nx*ny*nz)、vx/vy/vz は面 ((nx+1)*ny*nz など)。
 * 更新式 (leapfrog) :
 *   v^{n+1/2} = v^{n-1/2} - (Δt/ρ) ∇p^n
 *   p^{n+1}   = p^n - ρc²Δt ∇·v^{n+1/2}
 * c = 343.0 m/s, ρ = 1.204 kg/m³ (20 ℃)。単一の一様格子 dx = .ofd の最小刻み。
 *
 * 時間刻み : Δt = 1/fs、fs = ceil(c·√3/(0.99·dx)) [整数 Hz]。
 *   3 次元の安定条件は c·Δt·√3/dx ≤ 1 なので、この fs は CFL 数 0.99 を
 *   満たす最小の整数サンプリング周波数 (WAV の fs と時間軸が厳密に一致する)。
 *
 * 音源 (ソフト音源、ガウシアン微分パルス) :
 *   s(t) = -u·exp(1/2 - u²/2),  u = (t - t0)/σ  (ピーク振幅 1 に正規化)
 *   設計根拠 : 格子分解能から使える帯域は fmax = c/(10·dx) (10 セル/波長 —
 *   1 次元の位相速度誤差 (kΔx)²(1-S²)/24 が fmax で約 1% に収まる)。
 *   ガウシアン微分のスペクトルは |S(ω)| ∝ ω·exp(-ω²σ²/2) で ω_p = 1/σ が
 *   ピーク。σ = 2/(π·fmax) と選ぶと ω_max = 2π·fmax = 4/σ となり、
 *   |S(4/σ)|/|S(1/σ)| = 4·e^{-8+1/2} ≈ 2.2e-3 (-53 dB) — 帯域上限で十分減衰。
 *   微分形なので DC 成分は 0 (剛壁閉領域でも圧力オフセットが溜まらない)。
 *   遅延 t0 = 5σ (t=0 での振幅は 5·e^{-12} ≈ 3e-5 — 立ち上がり誤差は無視できる)。
 *
 * 剛体 (geometry) : セル中心が形状内にあるセルを剛体とし、剛体に接する面の
 *   法線速度を 0 に固定する (面マスク)。本家 OpenFDTD の全 shape (直方体 1、
 *   楕円体 2、円柱 11-13、三角柱 31-33、角錐台 41-43、円錐台 51-53) を
 *   それぞれの形状の内外判定でボクセル化する (外接箱で近似しない)。
 *   パラメータの並びは本家 sol/ingeometry.c (OpenFDTD-X core/Geometry.h が
 *   正) と同じ。未知の shape は warning を出して無視する。
 *   .ofdx の acoustic.obstacles[] で材質 (吸音率) を持つ障害物は、流体との
 *   界面を外壁と同じ半陰的インピーダンス更新にする (係数は bnd_coeffs を
 *   外壁と共有。面と最寄りの p 節点の間は外壁と同じ半セル)。材質なしの
 *   障害物の面はマスクで v ≡ 0 (剛体、厳密に無損失)。
 *
 * 外壁境界 (局所反応インピーダンス) :
 *   Z = ρc·(1+√(1-α))/(1-√(1-α))。これは平面波の圧力反射係数
 *   R = (Z-ρc)/(Z+ρc) = √(1-α) の逆変換そのもの (エネルギー反射率 R² = 1-α)。
 *   α→1 で Z→ρc (垂直入射無反射)、α→0 で Z→∞ (剛壁) に漸近する。
 *   離散化 : 壁と最寄り p 節点の間の半セルの運動方程式
 *     ρ(dx/2)·dv/dt = -(p - p_wall),  p_wall = Z·v_n  (v_n は壁へ向かう速度)
 *   を Z 項について時間中心の半陰解法で解く :
 *     v+ = wA·v- ∓ wB·p,  wA = (a - Z/2)/(a + Z/2),  wB = 1/(a + Z/2),
 *     a = ρ·dx/(2Δt)   (符号は外向き法線が -軸 なら -、+軸 なら +)
 *   |wA| < 1 なので無条件安定。Z→∞ で wA→-1, wB→0 (v≡0 = 剛壁)。
 *   α = 0 (Z = ∞) は厳密に更新をスキップして剛壁にする (エネルギー保存)。
 *
 * 決定性 : 乱数不使用。OpenMP 並列は p 更新・v 更新のセル/面ループのみで、
 * いずれも要素ごとに独立でリダクションを持たない (OpenPEEC と同じ原則) —
 * スレッド数によらずビット単位で一致する。MSVC OpenMP 2.0 のため
 * ループ変数は事前宣言の int (OpenPEEC portability.md)。
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "acoustic.h"

#define PROGRESS_TOTAL 50   /* stdout の "progress a/b" の分母 */

/* ── 幾何 : AABB (shape 1 は厳密に同じ箱)。戻り値 0 = 対応しない shape ── */
static int geom_aabb(const ac_geom_t *g, double lo[3], double hi[3])
{
	const double *p = g->g;
	int i;
	double a, b;
	switch (g->shape) {
		case 1: case 2: case 11: case 12: case 13:
			/* 直方体 / 楕円体 / 円柱 : パラメータ自体が境界箱 */
			for (i = 0; i < 3; i++) {
				a = p[2*i]; b = p[2*i + 1];
				lo[i] = (a < b) ? a : b;
				hi[i] = (a < b) ? b : a;
			}
			return 1;
		case 31: case 32: case 33: {
			/* 三角柱 : 軸範囲 g[0..1] + 断面頂点 (g[i+2], g[i+5]) i=0..2 */
			int ax = g->shape - 31;             /* 0=x, 1=y, 2=z */
			int u = (ax + 1) % 3, w = (ax + 2) % 3;  /* 断面の 2 軸 (本家の順) */
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
			/* 角錐 / 円錐 : 軸範囲 g[0..1]、断面中心 (g[2],g[3])、
			 * 全幅 max(g[4],g[6]) / max(g[5],g[7]) (本家 ingeometry.c の並び) */
			int ax = (g->shape % 10) - 1;       /* 0=x, 1=y, 2=z */
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

/* 点 q = (x, y, z) が形状の内部 (境界を含む) か。外接箱 (geom_aabb) の
 * 内側であることは呼び出し側が保証する。軸付き形状の断面 2 軸は
 * u = (軸+1)%3, w = (軸+2)%3 (本家 ingeometry.c の (u,v) 順)。
 *   2      楕円体   : 外接箱に内接 — sum ((q-c)/r)^2 <= 1
 *   11-13  円柱     : 軸 = shape-11、断面は外接箱に内接する楕円
 *   31-33  三角柱   : 軸 = shape-31、軸範囲 g[0..1]、断面の頂点
 *                     (g[2+k], g[5+k]) k = 0..2 (向きは問わない)
 *   41-43  角錐台   : 軸 = shape%10-1、軸範囲 g[0..1]、断面中心 (g[2], g[3])、
 *   51-53  円錐台     g[0] 側の全幅 (g[4], g[5])、g[1] 側の全幅 (g[6], g[7])
 *                     — 間は線形補間。角錐は矩形、円錐は楕円の断面
 * 退化した形状 (半径 0・軸長 0・面積 0) は内部を持たない。 */
static int geom_inside(const ac_geom_t *g, double x, double y, double z)
{
	const double *p = g->g;
	double q[3];
	q[0] = x; q[1] = y; q[2] = z;
	switch (g->shape) {
		case 1:
			return 1;                      /* 外接箱そのもの */
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
			int ax = g->shape - 11, a;
			double sum = 0.0;
			for (a = 0; a < 3; a++) {
				double c, r;
				if (a == ax) continue;         /* 軸方向は外接箱の範囲だけ */
				c = 0.5 * (p[2*a] + p[2*a + 1]);
				r = 0.5 * fabs(p[2*a + 1] - p[2*a]);
				if (r <= 0.0) return 0;
				sum += ((q[a] - c) / r) * ((q[a] - c) / r);
			}
			return sum <= 1.0;
		}
		case 31: case 32: case 33: {
			int ax = g->shape - 31;
			int u = (ax + 1) % 3, w = (ax + 2) % 3, k;
			double area2 = (p[3] - p[2]) * (p[7] - p[5])
			             - (p[6] - p[5]) * (p[4] - p[2]);
			int neg = 0, pos = 0;
			if (area2 == 0.0) return 0;
			for (k = 0; k < 3; k++) {          /* 3 辺の外積の符号がそろうか */
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
			int ax = (g->shape % 10) - 1;
			int u = (ax + 1) % 3, w = (ax + 2) % 3;
			double len = p[1] - p[0], t, hu, hw, du, dw;
			if (len == 0.0) return 0;
			t = (q[ax] - p[0]) / len;
			if (t < 0.0 || t > 1.0) return 0;
			hu = 0.5 * (p[4] + (p[6] - p[4]) * t);
			hw = 0.5 * (p[5] + (p[7] - p[5]) * t);
			du = q[u] - p[2];
			dw = q[w] - p[3];
			if (g->shape < 50)                 /* 角錐台 : 矩形断面 */
				return fabs(du) <= hu && fabs(dw) <= hw;
			if (hu <= 0.0 || hw <= 0.0) return 0;   /* 円錐台 : 楕円断面 */
			return (du / hu) * (du / hu) + (dw / hw) * (dw / hw) <= 1.0;
		}
		default:
			return 0;
	}
}

/* 形状の外接箱と重なるセル番号の範囲 (セル中心が箱の中)。戻り値 0 = 未知の
 * shape か、格子と重ならない。shape 1 はこの範囲がそのまま剛体になる
 * (従来と同じ ceil/floor — 既存入力のボクセル化はビット単位で不変)。 */
static int geom_irange(const ac_t *ac, const ac_geom_t *g, int r[6])
{
	double lo[3], hi[3];
	if (!geom_aabb(g, lo, hi)) return 0;
	r[0] = (int)ceil((lo[0] - ac->x0) / ac->dx - 0.5);
	r[1] = (int)floor((hi[0] - ac->x0) / ac->dx - 0.5);
	r[2] = (int)ceil((lo[1] - ac->y0) / ac->dx - 0.5);
	r[3] = (int)floor((hi[1] - ac->y0) / ac->dx - 0.5);
	r[4] = (int)ceil((lo[2] - ac->z0) / ac->dx - 0.5);
	r[5] = (int)floor((hi[2] - ac->z0) / ac->dx - 0.5);
	if (r[0] < 0) r[0] = 0;
	if (r[2] < 0) r[2] = 0;
	if (r[4] < 0) r[4] = 0;
	if (r[1] > ac->nx - 1) r[1] = ac->nx - 1;
	if (r[3] > ac->ny - 1) r[3] = ac->ny - 1;
	if (r[5] > ac->nz - 1) r[5] = ac->nz - 1;
	return 1;
}

/* 形状 g がセル (i, j, k) を剛体にするか (ボクセル化の唯一の規則) */
static int geom_marks(const ac_t *ac, const ac_geom_t *g, const int r[6],
                      int i, int j, int k)
{
	if (i < r[0] || i > r[1] || j < r[2] || j > r[3] || k < r[4] || k > r[5])
		return 0;
	if (g->shape == 1) return 1;
	return geom_inside(g, ac->x0 + (i + 0.5) * ac->dx,
	                      ac->y0 + (j + 0.5) * ac->dx,
	                      ac->z0 + (k + 0.5) * ac->dx);
}

/* セルを剛体にした形状の番号 (0 起点。無ければ -1)。エラーメッセージ用 —
 * ボクセル化と同じ規則 (geom_marks) を使うので、剛体判定とずれない。 */
static int geom_containing(const ac_t *ac, int i, int j, int k)
{
	int n, r[6];
	for (n = 0; n < ac->ngeom; n++) {
		if (!geom_irange(ac, &ac->geom[n], r)) continue;
		if (geom_marks(ac, &ac->geom[n], r, i, j, k)) return n;
	}
	return -1;
}

static int snap_cell(double x, double x0, double dx, int n)
{
	int i = (int)floor((x - x0) / dx);
	if (i < 0) i = 0;
	if (i > n - 1) i = n - 1;
	return i;
}

/* 名前 → ファイル名に使える形 (ファイル名に使えない ASCII だけを '_' へ)。
 * 非 ASCII (UTF-8 の日本語など) はそのまま残す — GUI の受音点名は既定で
 * 「マイク N」なので、ここで潰すと GUI 側の自動割当 (rir_<受音点名>.wav の
 * 名前照合) が一致しなくなる。Windows の UTF-8 名は ac_fopen (_wfopen) が扱う。 */
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

/* 受音点の WAV ファイル名 : #1 は契約どおり rir.wav、以降は rir_<名前>.wav
 * (名前が空なら rir_2.wav 等 — GUI の rir_<受音点名>.wav 自動割当と噛み合う)。
 * #1 に名前があるときは別名 rir_<名前>.wav も出す (alias) — 受音点が複数ある
 * とき、GUI の自動割当は名前でしか照合できず rir.wav だけでは #1 が
 * 割り当てられないため。 */
static void recv_filename(ac_recv_t *r, int index)
{
	char safe[AC_NAME_MAX];
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

/* ── 局所反応境界の統計入射吸音率と、その逆問題 (角度依存吸音 opt-in 用) ──
 * 実数の規格化インピーダンス zeta の面に平面波が角度 theta で入るとき
 *   R(theta) = (zeta cos(theta) - 1) / (zeta cos(theta) + 1)
 * FDTD の局所反応境界 (v_n = p/Z の半陰的更新) はこの角度依存性を離散化から
 * **自動的に**持つので、決めるべきは zeta ただ 1 つ。吸音表の alpha は
 * **ランダム入射** (残響室法) の値なので、Paris の式 (cos 重み平均)
 *   alpha_stat = 8 zeta int_0^1 x^2/(zeta x + 1)^2 dx
 *              = (8/zeta^2) [ zeta + 1 - 2 ln(1+zeta) - 1/(1+zeta) ]
 * を逆に解いて zeta を得る (表の値を垂直入射として使うと斜入射で吸音が
 * 二重に効く)。alpha_stat は zeta = 1.567 付近で最大 0.951 の山型なので、
 * 逆問題は「山より右 (zeta が大きい = 空気より硬い、通常の吸音材)」の枝で
 * 二分法を使い、最大値を超える alpha はクランプして警告する。
 * 幾何音響側 (ga_setup.c の ga_alpha_stat / ga_zeta_from_alpha) と同じ式・
 * 同じ枝 — 両ソルバーが同じ壁を見るための対称性なので、片方だけ変えない
 * (別バイナリなので実装は共有しない)。 */
static double ac_alpha_stat(double zeta)
{
	if (zeta <= 0.0) return 0.0;
	return (8.0 / (zeta * zeta))
	     * (zeta + 1.0 - 2.0 * log(1.0 + zeta) - 1.0 / (1.0 + zeta));
}

static double ac_zeta_from_alpha(double alpha_stat)
{
	double lo, hi = 1.0e9, peak;
	int i;
	{                                    /* 山の位置を三分探索 (決定的) */
		double a = 0.1, b = 10.0;
		for (i = 0; i < 200; i++) {
			double m1 = a + (b - a) / 3.0, m2 = b - (b - a) / 3.0;
			if (ac_alpha_stat(m1) < ac_alpha_stat(m2)) a = m1; else b = m2;
		}
		peak = 0.5 * (a + b);
	}
	if (alpha_stat >= ac_alpha_stat(peak)) return peak;
	if (alpha_stat <= 0.0) return hi;
	lo = peak;
	for (i = 0; i < 200; i++) {          /* alpha_stat は [peak, inf) で単調減少 */
		double mid = 0.5 * (lo + hi);
		if (ac_alpha_stat(mid) > alpha_stat) lo = mid; else hi = mid;
	}
	return 0.5 * (lo + hi);
}

/* 実効吸音率 a から境界係数を作る (外壁と障害物の面で共通 — 同じ式で
 * なければ「吸音率 a の障害物の面 = 吸音率 a の外壁」が崩れる)。
 * 戻り値 1 = 剛 (更新しない — 厳密に無損失)。 */
static int bnd_coeffs(const ac_t *ac, double a, double *Zo, double *A,
                      double *B, double *zeta)
{
	const double c = AC_C0, rho = AC_RHO0;
	double Z, ma;
	if (a < 1e-6) {
		*Zo = *A = *B = *zeta = 0.0;
		return 1;
	}
	ma = rho * ac->dx / (2.0 * ac->dt);
	if (ac->angle_dep) {
		double amax = ac_alpha_stat(ac_zeta_from_alpha(2.0));
		*zeta = ac_zeta_from_alpha((a > amax) ? amax : a);
		Z = rho * c * *zeta;
	}
	else {
		double s = sqrt(1.0 - ((a > 1.0) ? 1.0 : a));
		Z = rho * c * (1.0 + s) / (1.0 - s + 1e-300);
		*zeta = Z / (rho * c);
	}
	*Zo = Z;
	*A = (ma - Z / 2.0) / (ma + Z / 2.0);
	*B = 1.0 / (ma + Z / 2.0);
	return 0;
}

int ac_setup(ac_t *ac)
{
	const double c = AC_C0, rho = AC_RHO0;
	double lx, ly, lz, vol, area, tsolve;
	long long cells;
	size_t np, nvx, nvy, nvz;
	unsigned short *gown = NULL;   /* 剛体セル → 形状番号+1 (吸音性の障害物があるときだけ) */
	int i, j, k, n, w;
	int nx, ny, nz;

	/* ── 格子 : dx = 最小刻み、セル数は領域を丸めて決める ── */
	ac->dx = ac->dxmin;
	if (ac->dx <= 0.0) {
		ac_err(ac, "invalid mesh spacing (dx = %g)", ac->dx);
		return 1;
	}
	ac->nx = (int)floor((ac->x1 - ac->x0) / ac->dx + 0.5);
	ac->ny = (int)floor((ac->y1 - ac->y0) / ac->dx + 0.5);
	ac->nz = (int)floor((ac->z1 - ac->z0) / ac->dx + 0.5);
	if (ac->nx < 1) ac->nx = 1;
	if (ac->ny < 1) ac->ny = 1;
	if (ac->nz < 1) ac->nz = 1;
	nx = ac->nx; ny = ac->ny; nz = ac->nz;

	cells = (long long)nx * ny * nz;
	if (cells > AC_MAX_CELLS) {
		ac_err(ac, "grid too large: %d x %d x %d = %lld cells exceeds the "
		       "limit of %d — coarsen the .ofd mesh (larger minimum spacing) "
		       "or reduce the domain size", nx, ny, nz, cells, AC_MAX_CELLS);
		return 1;
	}

	/* ── 時間刻み : CFL を満たす最小の整数 fs ── */
	ac->fs = (int)ceil(c * sqrt(3.0) / (AC_CFL * ac->dx));
	ac->dt = 1.0 / ac->fs;

	/* ── 音源パルス設計 (冒頭コメント参照) ── */
	ac->fmax  = c / (10.0 * ac->dx);
	ac->sigma = 2.0 / (AC_PI * ac->fmax);
	ac->t0    = 5.0 * ac->sigma;

	/* ── 実効吸音率 : 有効帯域 [0, fmax] と重なるオクターブバンドだけを平均 ──
	 * .ofdx の alpha[] は 125 Hz .. 4 kHz の 6 バンド。本ソルバーは低域担当で
	 * 上限が fmax = c/(10 dx) なので、そこより上のバンド (ホール規模なら
	 * 1 kHz 以上) の吸音率まで混ぜて壁を決めるのは誤り — 幾何音響側
	 * (バンド別に使う) とのクロスオーバー整合も崩れる。バンド b の範囲は
	 * [fc_b/sqrt(2), fc_b*sqrt(2)] なので、fc_b <= fmax*sqrt(2) のバンドが
	 * 帯域と重なる。1 つも重ならない (fmax < 88.4 Hz) 場合は最低バンド
	 * 125 Hz の値を使う (壁を既定値に落とさない)。 */
	{
		double fc, hi = AC_BAND_F0 * AC_SQRT2;
		int nb = 0, b;
		for (b = 0, fc = AC_BAND_F0; b < AC_NBAND_OFDX; b++, fc *= 2.0) {
			if (fc > ac->fmax * AC_SQRT2) break;
			nb++;
			hi = fc * AC_SQRT2;
		}
		ac->alpha_nband   = (nb > 0) ? nb : 1;
		ac->alpha_band_hi = hi;
		for (w = 0; w < AC_NWALL; w++) {
			double sum = 0.0;
			for (b = 0; b < ac->alpha_nband; b++) sum += ac->alpha_bands[w][b];
			ac->alpha[w] = sum / ac->alpha_nband;
		}
		if (nb > 0)
			ac_log(ac, "absorption: using %d of %d octave bands "
			       "(125 Hz .. %g Hz, band edges up to %.4g Hz <= fmax*sqrt(2))",
			       ac->alpha_nband, AC_NBAND_OFDX,
			       AC_BAND_F0 * pow(2.0, (double)(ac->alpha_nband - 1)), hi);
		else
			ac_log(ac, "absorption: warning: fmax = %.4g Hz is below the "
			       "lowest octave band (125 Hz, lower edge %.4g Hz) — using the "
			       "125 Hz alpha for all walls", ac->fmax, AC_BAND_F0 / AC_SQRT2);
	}

	/* ── 障害物の材質 (.ofdx acoustic.obstacles[]) : 外壁と同じ規則で実効値
	 * (有効帯域と重なるバンドの平均) を決める。指定なしは剛体 (alpha = 0)。 ── */
	if (ac->ngeom > 0) {
		ac->gvox = (ac_gvox_t *)calloc((size_t)ac->ngeom, sizeof(ac_gvox_t));
		if (!ac->gvox) {
			ac_err(ac, "out of memory (%d geometries)", ac->ngeom);
			return 1;
		}
	}
	for (n = 0; n < ac->ngeom; n++) {
		double sum = 0.0;
		int b;
		if (!ac->obs_has || !ac->obs_has[n]) continue;
		for (b = 0; b < ac->alpha_nband; b++)
			sum += ac->obs_bands[(size_t)n * AC_NBAND_OFDX + b];
		ac->gvox[n].alpha = sum / ac->alpha_nband;
	}

	/* ── 計算時間 : T = clamp(1.5·T_Sabine, 0.5, 3.0)
	 * T_Sabine = 0.161·V/A。A は吸音表 (壁 6 面 × 実効 α) から。
	 * 表が無い場合も既定 α=0.1 の 6 面で同じ式 (input_ofdx.c の既定値)。 ── */
	lx = nx * ac->dx; ly = ny * ac->dx; lz = nz * ac->dx;
	vol = lx * ly * lz;
	area = ac->alpha[AC_XM] * ly * lz + ac->alpha[AC_XP] * ly * lz
	     + ac->alpha[AC_YM] * lx * lz + ac->alpha[AC_YP] * lx * lz
	     + ac->alpha[AC_ZM] * lx * ly + ac->alpha[AC_ZP] * lx * ly;
	ac->tsab = (area > 1e-12) ? AC_SABINE_COEF * vol / area : -1.0;
	tsolve = (ac->tsab > 0.0) ? 1.5 * ac->tsab : AC_TMAX;
	if (tsolve < AC_TMIN) tsolve = AC_TMIN;
	if (tsolve > AC_TMAX) tsolve = AC_TMAX;

	/* ── 音源ごとのゲイン・遅延 (ADR-0010 Decision 7) ──
	 * .ofdx acoustic.sources[] が無ければ既定値 (gain = 1, delay = 0) で
	 * 確保する — 以降は常に配列越しに参照でき、既定値なら従来動作と
	 * ビット単位で一致する (1.0 * s(t - 0.0) は s(t) と IEEE で等価)。
	 * 発火する音源集合の max(delay) だけ計算時間を延長する (遅延で残響が
	 * 窓の外にこぼれない)。 */
	if (!ac->srcgain) {
		ac->srcgain  = (double *)malloc((size_t)(ac->nfeed > 0 ? ac->nfeed : 1)
		                                * sizeof(double));
		ac->srcdelay = (double *)calloc((size_t)(ac->nfeed > 0 ? ac->nfeed : 1),
		                                sizeof(double));
		if (!ac->srcgain || !ac->srcdelay) {
			ac_err(ac, "out of memory (%d source gains)", ac->nfeed);
			return 1;
		}
		for (n = 0; n < ac->nfeed; n++) ac->srcgain[n] = 1.0;
	}
	{
		double dmax = 0.0;
		int nact = (ac->multi_source && ac->nfeed > 0) ? ac->nfeed : 1;
		for (n = 0; n < nact && n < ac->nfeed; n++)
			if (ac->srcdelay[n] > dmax) dmax = ac->srcdelay[n];
		if (dmax > 0.0)
			ac_log(ac, "duration: extended by max source delay %.4g s "
			       "(.ofdx acoustic.sources)", dmax);
		tsolve += dmax;
	}
	ac->nsteps = (int)ceil(tsolve * ac->fs);
	ac->duration = ac->nsteps * ac->dt;

	ac_log(ac, "grid: dx = %.6g m, %d x %d x %d = %lld cells "
	       "(effective domain %.6g x %.6g x %.6g m)",
	       ac->dx, nx, ny, nz, cells, lx, ly, lz);
	ac_log(ac, "time: fs = %d Hz (CFL %.4f), dt = %.6e s",
	       ac->fs, c * ac->dt * sqrt(3.0) / ac->dx, ac->dt);
	ac_log(ac, "source: gaussian-derivative pulse, fmax = %.6g Hz "
	       "(c/(10 dx)), sigma = %.6e s, t0 = %.6e s",
	       ac->fmax, ac->sigma, ac->t0);
	if (ac->tsab > 0.0)
		ac_log(ac, "duration: T_Sabine = %.4g s (V = %.4g m^3, A = %.4g m^2) "
		       "-> T = %.4g s (%d steps)", ac->tsab, vol, area, tsolve, ac->nsteps);
	else
		ac_log(ac, "duration: all walls rigid (A = 0, T_Sabine = inf) "
		       "-> T = %.4g s (%d steps)", tsolve, ac->nsteps);

	/* ── 境界係数 ──
	 * 既定 : 吸音表の alpha を**垂直入射**として読み、R = sqrt(1-alpha) から
	 *   Z = rho c (1+R)/(1-R) を作る (従来動作)。
	 * angle_dependent_absorption = true : alpha を**ランダム入射**として読み、
	 *   Paris の式を逆に解いた zeta で Z = rho c zeta にする。境界の角度依存性
	 *   R(theta) = (zeta cos - 1)/(zeta cos + 1) は離散化から自動的に出る。 */
	{
		double amax = ac_alpha_stat(ac_zeta_from_alpha(2.0));   /* 山の高さ */
		int clamped = 0;
		for (w = 0; w < AC_NWALL; w++)
			if (ac->angle_dep && ac->alpha[w] > amax) clamped = 1;
		for (n = 0; n < ac->ngeom; n++)
			if (ac->angle_dep && ac->gvox[n].alpha > amax) clamped = 1;
		if (ac->angle_dep) {
			ac_log(ac, "angle-dependent absorption: locally reacting boundary "
			       "Z = rho c zeta; zeta solved from the random-incidence alpha "
			       "via Paris' formula (max representable alpha_stat = %.4f). "
			       "R(theta) = (zeta cos - 1)/(zeta cos + 1) follows from the "
			       "update itself", amax);
			if (clamped)
				ac_log(ac, "warning: some alpha exceed %.4f, the largest "
				       "random-incidence absorption a locally reacting real "
				       "impedance can have — clamped (the wall is treated as "
				       "the most absorbing such boundary, not as alpha = 1)",
				       amax);
		}
	}
	for (w = 0; w < AC_NWALL; w++) {
		/* 剛壁 (alpha = 0) は更新なし (v ≡ 0) — 厳密に無損失 */
		ac->wrigid[w] = bnd_coeffs(ac, ac->alpha[w], &ac->wZ[w], &ac->wA[w],
		                           &ac->wB[w], &ac->wzeta[w]);
		ac_log(ac, "wall %d (%s): alpha = %.4g (mean of %d band(s), "
		       "<= %.4g Hz), zeta = %.4g -> %s", w,
		       (w == AC_XM) ? "x-" : (w == AC_XP) ? "x+" :
		       (w == AC_YM) ? "y-" : (w == AC_YP) ? "y+" :
		       (w == AC_ZM) ? "z-" : "z+",
		       ac->alpha[w], ac->alpha_nband, ac->alpha_band_hi,
		       ac->wzeta[w], ac->wrigid[w] ? "rigid" : "impedance");
	}

	/* ── 配列確保 (calloc = 零初期化。決定性のため未初期化領域を作らない) ── */
	np  = (size_t)nx * ny * nz;
	nvx = (size_t)(nx + 1) * ny * nz;
	nvy = (size_t)nx * (ny + 1) * nz;
	nvz = (size_t)nx * ny * (nz + 1);
	ac->p     = (double *)calloc(np,  sizeof(double));
	ac->vx    = (double *)calloc(nvx, sizeof(double));
	ac->vy    = (double *)calloc(nvy, sizeof(double));
	ac->vz    = (double *)calloc(nvz, sizeof(double));
	ac->solid = (unsigned char *)calloc(np,  1);
	ac->mvx   = (unsigned char *)calloc(nvx, 1);
	ac->mvy   = (unsigned char *)calloc(nvy, 1);
	ac->mvz   = (unsigned char *)calloc(nvz, 1);
	ac->rec   = (double *)calloc((size_t)ac->nrecv * ac->nsteps, sizeof(double));
	if (!ac->p || !ac->vx || !ac->vy || !ac->vz || !ac->solid ||
	    !ac->mvx || !ac->mvy || !ac->mvz || !ac->rec) {
		ac_err(ac, "out of memory (%lld cells — coarsen the mesh)", cells);
		return 1;
	}

	/* ── 剛体ボクセル化 (セル中心が形状内 → 剛体) ──
	 * 外接箱と重なる範囲だけを走査し、各セル中心を形状の内外判定にかける。
	 * 形状ごとのセル数・重心・外接範囲を記録する (solver.log / metadata)。 */
	{
		/* 吸音性の障害物があるときだけ、各剛体セルを「どの形状が最後に剛体に
		 * したか」で引けるようにする (重なりは後の行が優先 — GUI と同じ)。
		 * 界面リストを作ったら捨てる一時配列。 */
		int any = 0;
		for (n = 0; n < ac->ngeom; n++)
			if (ac->gvox[n].alpha >= 1e-6) any = 1;
		if (any) {
			if (ac->ngeom > 65535) {
				ac_err(ac, "too many geometries with acoustic.obstacles "
				       "materials (%d > 65535)", ac->ngeom);
				return 1;
			}
			gown = (unsigned short *)calloc(np, sizeof(unsigned short));
			if (!gown) {
				ac_err(ac, "out of memory (obstacle material map, %lld cells)",
				       cells);
				return 1;
			}
		}
	}
	for (n = 0; n < ac->ngeom; n++) {
		const ac_geom_t *g = &ac->geom[n];
		ac_gvox_t *gv = &ac->gvox[n];
		int r[6];
		double sx = 0.0, sy = 0.0, sz = 0.0;
		long ilo = nx, ihi = -1, jlo = ny, jhi = -1, klo = nz, khi = -1;
		if (!geom_irange(ac, g, r)) {
			if (!geom_aabb(g, gv->lo, gv->hi))
				ac_log(ac, "warning: geometry #%d has unknown shape %d — ignored",
				       n + 1, g->shape);
			continue;
		}
		for (k = r[4]; k <= r[5]; k++)
			for (j = r[2]; j <= r[3]; j++)
				for (i = r[0]; i <= r[1]; i++) {
					if (!geom_marks(ac, g, r, i, j, k)) continue;
					ac->solid[((size_t)k * ny + j) * nx + i] = 1;
					if (gown)
						gown[((size_t)k * ny + j) * nx + i] = (unsigned short)(n + 1);
					gv->ncell++;
					sx += i; sy += j; sz += k;
					if (i < ilo) ilo = i;
					if (i > ihi) ihi = i;
					if (j < jlo) jlo = j;
					if (j > jhi) jhi = j;
					if (k < klo) klo = k;
					if (k > khi) khi = k;
				}
		if (gv->ncell > 0) {
			gv->cen[0] = ac->x0 + (sx / gv->ncell + 0.5) * ac->dx;
			gv->cen[1] = ac->y0 + (sy / gv->ncell + 0.5) * ac->dx;
			gv->cen[2] = ac->z0 + (sz / gv->ncell + 0.5) * ac->dx;
			gv->lo[0] = ac->x0 + ilo * ac->dx;  gv->hi[0] = ac->x0 + (ihi + 1) * ac->dx;
			gv->lo[1] = ac->y0 + jlo * ac->dx;  gv->hi[1] = ac->y0 + (jhi + 1) * ac->dx;
			gv->lo[2] = ac->z0 + klo * ac->dx;  gv->hi[2] = ac->z0 + (khi + 1) * ac->dx;
		}
		else {
			int a;
			for (a = 0; a < 3; a++) gv->cen[a] = gv->lo[a] = gv->hi[a] = 0.0;
		}
		ac_log(ac, "voxelizer: geometry #%d shape %d -> %ld cells "
		       "(%.6g m^3), centroid (%.4g, %.4g, %.4g) m, extent "
		       "[%.4g, %.4g] x [%.4g, %.4g] x [%.4g, %.4g] m", n + 1, g->shape,
		       gv->ncell, gv->ncell * ac->dx * ac->dx * ac->dx,
		       gv->cen[0], gv->cen[1], gv->cen[2],
		       gv->lo[0], gv->hi[0], gv->lo[1], gv->hi[1], gv->lo[2], gv->hi[2]);
		if (gv->ncell == 0)
			ac_log(ac, "warning: geometry #%d (shape %d) covers no cell center "
			       "— it has no effect at dx = %g m (thinner than the grid or "
			       "outside the domain)", n + 1, g->shape, ac->dx);
	}
	{
		size_t c2, nsolid = 0;
		for (c2 = 0; c2 < np; c2++) nsolid += ac->solid[c2];
		if (ac->ngeom > 0)
			ac_log(ac, "voxelizer: %d geometries (all treated as rigid) -> "
			       "%zu solid cells", ac->ngeom, nsolid);
	}

	/* ── 面マスク : 両側とも流体の内部面のみ更新する ── */
	for (k = 0; k < nz; k++)
		for (j = 0; j < ny; j++)
			for (i = 1; i < nx; i++)
				ac->mvx[((size_t)k * ny + j) * (nx + 1) + i] =
					!ac->solid[((size_t)k * ny + j) * nx + i - 1] &&
					!ac->solid[((size_t)k * ny + j) * nx + i];
	for (k = 0; k < nz; k++)
		for (j = 1; j < ny; j++)
			for (i = 0; i < nx; i++)
				ac->mvy[((size_t)k * (ny + 1) + j) * nx + i] =
					!ac->solid[((size_t)k * ny + j - 1) * nx + i] &&
					!ac->solid[((size_t)k * ny + j) * nx + i];
	for (k = 1; k < nz; k++)
		for (j = 0; j < ny; j++)
			for (i = 0; i < nx; i++)
				ac->mvz[((size_t)k * ny + j) * nx + i] =
					!ac->solid[((size_t)(k - 1) * ny + j) * nx + i] &&
					!ac->solid[((size_t)k * ny + j) * nx + i];

	/* ── 吸音性の障害物の界面 (流体セルと、材質を持つ剛体セルの間の面) ──
	 * 外壁と同じ半陰的インピーダンス更新をする (面と最寄りの p 節点の間は
	 * 外壁と同じ半セル)。剛体が面の負側なら外壁 x- と同じ式、正側なら x+ と
	 * 同じ式。数え上げ → 確保 → 記入の 2 パス。材質なしの剛体の面は従来どおり
	 * マスクで v ≡ 0 のまま (厳密に無損失)。 */
	if (gown) {
		long cnt = 0, pass;
		for (pass = 0; pass < 2; pass++) {
			long m = 0;
			int ax;
			for (ax = 0; ax < 3; ax++) {
				const int ex = (ax == 0), ey = (ax == 1), ez = (ax == 2);
				for (k = ez; k < nz; k++)
					for (j = ey; j < ny; j++)
						for (i = ex; i < nx; i++) {
							size_t ca = ((size_t)(k - ez) * ny + (j - ey)) * nx + (i - ex);
							size_t cb = ((size_t)k * ny + j) * nx + i;
							size_t cs, cf;
							int g, sg;
							if (ac->solid[ca] == ac->solid[cb]) continue;
							cs = ac->solid[ca] ? ca : cb;      /* 剛体側 */
							cf = ac->solid[ca] ? cb : ca;      /* 流体側 */
							sg = ac->solid[ca] ? -1 : +1;      /* 剛体が負側なら -1 */
							g = (int)gown[cs] - 1;
							if (g < 0 || ac->gvox[g].alpha < 1e-6) continue;
							if (pass == 1) {
								ac_oface_t *o = &ac->oface[m];
								double Zd, zd;
								o->axis = ax;
								o->sgn = sg;
								o->pi = cf;
								o->vi = (ax == 0) ? ((size_t)k * ny + j) * (nx + 1) + i
								      : (ax == 1) ? ((size_t)k * (ny + 1) + j) * nx + i
								                  : ((size_t)k * ny + j) * nx + i;
								bnd_coeffs(ac, ac->gvox[g].alpha, &Zd, &o->A, &o->B, &zd);
								ac->gvox[g].nface++;
							}
							m++;
						}
			}
			if (pass == 0) {
				cnt = m;
				ac->oface = (ac_oface_t *)calloc((size_t)(cnt > 0 ? cnt : 1),
				                                 sizeof(ac_oface_t));
				if (!ac->oface) {
					free(gown);
					ac_err(ac, "out of memory (%ld obstacle faces)", cnt);
					return 1;
				}
			}
		}
		ac->noface = cnt;
		free(gown);
		gown = NULL;
	}
	for (n = 0; n < ac->ngeom; n++) {
		ac_gvox_t *gv = &ac->gvox[n];
		double Zd, Ad, Bd;
		if (!ac->obs_has || !ac->obs_has[n]) continue;
		bnd_coeffs(ac, gv->alpha, &Zd, &Ad, &Bd, &gv->zeta);
		ac_log(ac, "obstacle #%d (shape %d): alpha = %.4g (mean of %d band(s), "
		       "<= %.4g Hz), zeta = %.4g -> %s, %ld fluid faces", n + 1,
		       ac->geom[n].shape, gv->alpha, ac->alpha_nband, ac->alpha_band_hi,
		       gv->zeta, (gv->alpha < 1e-6) ? "rigid" : "impedance", gv->nface);
	}

	/* ── 使用する音源集合 (複数音源の契約 : ADR-0010) ──
	 * 既定は feed #1 のみ (従来動作と完全一致)。multi_source なら全 feed に
	 * 同一パルス (共通 t0) を注入し、RIR は重ね合わせになる。 */
	ac->nsrc = (ac->multi_source && ac->nfeed > 0) ? ac->nfeed : 1;
	if (ac->nfeed > 1 && !ac->multi_source)
		ac_log(ac, "warning: %d feeds found — using feed #1 only "
		       "(set acoustic.multi_source in the .ofdx to sum all sources)",
		       ac->nfeed);
	if (ac->multi_source)
		ac_log(ac, "multi_source: %d feed(s) fire the same pulse "
		       "simultaneously — rir.wav is the superposition "
		       "(unit strength each, no 1/N)", ac->nsrc);
	ac->srccell = (int *)calloc((size_t)ac->nsrc * 3, sizeof(int));
	if (!ac->srccell) {
		ac_err(ac, "out of memory (%d sources)", ac->nsrc);
		return 1;
	}

	/* ── 音源・受音点のスナップ (剛体セル内は誤りなので正直に失敗する) ── */
	for (n = 0; n < ac->nsrc; n++) {
		double sx = ac->feedpos[3 * n + 0];
		double sy = ac->feedpos[3 * n + 1];
		double sz = ac->feedpos[3 * n + 2];
		int is = snap_cell(sx, ac->x0, ac->dx, nx);
		int js = snap_cell(sy, ac->y0, ac->dx, ny);
		int ks = snap_cell(sz, ac->z0, ac->dx, nz);
		if (ac->solid[((size_t)ks * ny + js) * nx + is]) {
			int gi = geom_containing(ac, is, js, ks);
			if (gi >= 0)
				ac_err(ac, "feed #%d position (%.4g, %.4g, %.4g) is inside "
				       "rigid geometry #%d (shape %d) — move the source out "
				       "of the object (e.g. 1.5 m above the stage floor)",
				       n + 1, sx, sy, sz, gi + 1, ac->geom[gi].shape);
			else
				ac_err(ac, "feed #%d position (%.4g, %.4g, %.4g) is inside "
				       "rigid geometry — move the source out of the object "
				       "(e.g. 1.5 m above the stage floor)", n + 1, sx, sy, sz);
			return 1;
		}
		ac->srccell[3 * n + 0] = is;
		ac->srccell[3 * n + 1] = js;
		ac->srccell[3 * n + 2] = ks;
		if (n == 0) { ac->isrc = is; ac->jsrc = js; ac->ksrc = ks; }
		ac_log(ac, "source #%d: (%.6g, %.6g, %.6g) m -> cell (%d, %d, %d), "
		       "center (%.6g, %.6g, %.6g) m, gain = %.4g, delay = %.4g s",
		       n + 1, sx, sy, sz, is, js, ks,
		       ac->x0 + (is + 0.5) * ac->dx,
		       ac->y0 + (js + 0.5) * ac->dx,
		       ac->z0 + (ks + 0.5) * ac->dx,
		       ac->srcgain[n], ac->srcdelay[n]);
	}
	for (n = 0; n < ac->nrecv; n++) {
		ac_recv_t *r = &ac->recv[n];
		r->ic = snap_cell(r->x, ac->x0, ac->dx, nx);
		r->jc = snap_cell(r->y, ac->y0, ac->dx, ny);
		r->kc = snap_cell(r->z, ac->z0, ac->dx, nz);
		if (ac->solid[((size_t)r->kc * ny + r->jc) * nx + r->ic]) {
			/* 剛体内の受音点は測れない。どの形状に入っているかと、
			 * 典型的な直し方 (耳の高さへ上げる) を添えて失敗する。 */
			int gi = geom_containing(ac, r->ic, r->jc, r->kc);
			if (gi >= 0)
				ac_err(ac, "point #%d (%.4g, %.4g, %.4g) is inside rigid "
				       "geometry #%d (shape %d) — move it out of the object, "
				       "e.g. to ear height above the floor/audience block "
				       "(ISO 3382-1: 1.2 m seated)",
				       n + 1, r->x, r->y, r->z, gi + 1, ac->geom[gi].shape);
			else
				ac_err(ac, "point #%d (%.4g, %.4g, %.4g) is inside rigid "
				       "geometry — move it out of the object, e.g. to ear "
				       "height above the floor (ISO 3382-1: 1.2 m seated)",
				       n + 1, r->x, r->y, r->z);
			return 1;
		}
		recv_filename(r, n);
		ac_log(ac, "point #%d%s%s: (%.6g, %.6g, %.6g) m -> cell (%d, %d, %d) "
		       "-> %s", n + 1, r->name[0] ? " " : "", r->name,
		       r->x, r->y, r->z, r->ic, r->jc, r->kc, r->file);
	}
	return 0;
}

/* 音源波形 (冒頭コメントの設計) */
static double src_pulse(const ac_t *ac, double t)
{
	const double u = (t - ac->t0) / ac->sigma;
	return -u * exp(0.5 - 0.5 * u * u);
}

int ac_run(ac_t *ac)
{
	const int nx = ac->nx, ny = ac->ny, nz = ac->nz;
	const double cv = ac->dt / (AC_RHO0 * ac->dx);          /* v 更新係数 */
	const double cp = AC_RHO0 * AC_C0 * AC_C0 * ac->dt / ac->dx; /* p 更新係数 */
	double *p = ac->p, *vx = ac->vx, *vy = ac->vy, *vz = ac->vz;
	const unsigned char *mvx = ac->mvx, *mvy = ac->mvy, *mvz = ac->mvz;
	const unsigned char *solid = ac->solid;
	int n, i, j, k, r, prog_done = 0;

	/* 実行条件を stdout にも 1 行出す。GUI のログ枠に見えるのは stdout だけで
	 * (solver.log は作業ディレクトリのファイル)、進捗行の分母は
	 * PROGRESS_TOTAL 固定なので、実際のステップ数がここに無いと
	 * 「progress は 50 固定？」が分からない。 */
	printf("solve: %d steps, %d cells, fs = %d Hz, %d receivers "
	       "(進捗は %d 分割)\n",
	       ac->nsteps, ac->nx * ac->ny * ac->nz, ac->fs, ac->nrecv,
	       PROGRESS_TOTAL);
	fflush(stdout);

	for (n = 0; n < ac->nsteps; n++) {
		/* ── v 更新 (内部面) : 面ごとに独立 — リダクションなし ── */
#ifdef _OPENMP
#pragma omp parallel for private(i, j)
#endif
		for (k = 0; k < nz; k++) {
			for (j = 0; j < ny; j++) {
				const double *pp = p + ((size_t)k * ny + j) * nx;
				double *v = vx + ((size_t)k * ny + j) * (nx + 1);
				const unsigned char *m = mvx + ((size_t)k * ny + j) * (nx + 1);
				for (i = 1; i < nx; i++)
					v[i] -= cv * m[i] * (pp[i] - pp[i - 1]);
			}
		}
#ifdef _OPENMP
#pragma omp parallel for private(i, j)
#endif
		for (k = 0; k < nz; k++) {
			for (j = 1; j < ny; j++) {
				const double *pa = p + ((size_t)k * ny + j - 1) * nx;
				const double *pb = p + ((size_t)k * ny + j) * nx;
				double *v = vy + ((size_t)k * (ny + 1) + j) * nx;
				const unsigned char *m = mvy + ((size_t)k * (ny + 1) + j) * nx;
				for (i = 0; i < nx; i++)
					v[i] -= cv * m[i] * (pb[i] - pa[i]);
			}
		}
#ifdef _OPENMP
#pragma omp parallel for private(i, j)
#endif
		for (k = 1; k < nz; k++) {
			for (j = 0; j < ny; j++) {
				const double *pa = p + ((size_t)(k - 1) * ny + j) * nx;
				const double *pb = p + ((size_t)k * ny + j) * nx;
				double *v = vz + ((size_t)k * ny + j) * nx;
				const unsigned char *m = mvz + ((size_t)k * ny + j) * nx;
				for (i = 0; i < nx; i++)
					v[i] -= cv * m[i] * (pb[i] - pa[i]);
			}
		}

		/* ── 外壁 (インピーダンス境界)。剛壁はスキップ (v ≡ 0) ── */
		if (!ac->wrigid[AC_XM]) {
			const double A = ac->wA[AC_XM], B = ac->wB[AC_XM];
			for (k = 0; k < nz; k++)
				for (j = 0; j < ny; j++)
					if (!solid[((size_t)k * ny + j) * nx]) {
						double *v = vx + ((size_t)k * ny + j) * (nx + 1);
						*v = A * (*v) - B * p[((size_t)k * ny + j) * nx];
					}
		}
		if (!ac->wrigid[AC_XP]) {
			const double A = ac->wA[AC_XP], B = ac->wB[AC_XP];
			for (k = 0; k < nz; k++)
				for (j = 0; j < ny; j++)
					if (!solid[((size_t)k * ny + j) * nx + nx - 1]) {
						double *v = vx + ((size_t)k * ny + j) * (nx + 1) + nx;
						*v = A * (*v) + B * p[((size_t)k * ny + j) * nx + nx - 1];
					}
		}
		if (!ac->wrigid[AC_YM]) {
			const double A = ac->wA[AC_YM], B = ac->wB[AC_YM];
			for (k = 0; k < nz; k++)
				for (i = 0; i < nx; i++)
					if (!solid[((size_t)k * ny) * nx + i]) {
						double *v = vy + ((size_t)k * (ny + 1)) * nx + i;
						*v = A * (*v) - B * p[((size_t)k * ny) * nx + i];
					}
		}
		if (!ac->wrigid[AC_YP]) {
			const double A = ac->wA[AC_YP], B = ac->wB[AC_YP];
			for (k = 0; k < nz; k++)
				for (i = 0; i < nx; i++)
					if (!solid[((size_t)k * ny + ny - 1) * nx + i]) {
						double *v = vy + ((size_t)k * (ny + 1) + ny) * nx + i;
						*v = A * (*v) + B * p[((size_t)k * ny + ny - 1) * nx + i];
					}
		}
		if (!ac->wrigid[AC_ZM]) {
			const double A = ac->wA[AC_ZM], B = ac->wB[AC_ZM];
			for (j = 0; j < ny; j++)
				for (i = 0; i < nx; i++)
					if (!solid[(size_t)j * nx + i]) {
						double *v = vz + (size_t)j * nx + i;
						*v = A * (*v) - B * p[(size_t)j * nx + i];
					}
		}
		if (!ac->wrigid[AC_ZP]) {
			const double A = ac->wA[AC_ZP], B = ac->wB[AC_ZP];
			for (j = 0; j < ny; j++)
				for (i = 0; i < nx; i++)
					if (!solid[((size_t)(nz - 1) * ny + j) * nx + i]) {
						double *v = vz + ((size_t)nz * ny + j) * nx + i;
						*v = A * (*v) + B * p[((size_t)(nz - 1) * ny + j) * nx + i];
					}
		}

		/* ── 吸音性の障害物の面 : 外壁と同じ式 (面ごとに独立 — リダクションなし) ── */
		{
			long f;
			for (f = 0; f < ac->noface; f++) {
				const ac_oface_t *o = &ac->oface[f];
				double *v = ((o->axis == 0) ? vx : (o->axis == 1) ? vy : vz) + o->vi;
				if (o->sgn < 0)
					*v = o->A * (*v) - o->B * p[o->pi];
				else
					*v = o->A * (*v) + o->B * p[o->pi];
			}
		}

		/* ── p 更新 : セルごとに独立 — リダクションなし。剛体セルは周囲の
		 * v が 0 に固定されているため発散が 0 で自然に p = 0 のまま ── */
#ifdef _OPENMP
#pragma omp parallel for private(i, j)
#endif
		for (k = 0; k < nz; k++) {
			for (j = 0; j < ny; j++) {
				double *pp = p + ((size_t)k * ny + j) * nx;
				const double *fx = vx + ((size_t)k * ny + j) * (nx + 1);
				const double *fy0 = vy + ((size_t)k * (ny + 1) + j) * nx;
				const double *fy1 = vy + ((size_t)k * (ny + 1) + j + 1) * nx;
				const double *fz0 = vz + ((size_t)k * ny + j) * nx;
				const double *fz1 = vz + ((size_t)(k + 1) * ny + j) * nx;
				for (i = 0; i < nx; i++)
					pp[i] -= cp * (fx[i + 1] - fx[i] + fy1[i] - fy0[i]
					               + fz1[i] - fz0[i]);
			}
		}

		/* ── ソフト音源 (直列)。multi_source なら全音源セルへ同一波形の
		 * パルス — 離散更新は線形なので、結果は各音源単独の RIR の厳密な和。
		 * 音源ごとのゲイン・遅延 (ADR-0010 Decision 7) は gain * s(t - delay)
		 * (既定 gain = 1 / delay = 0 は従来とビット等価)。 ── */
		{
			const double t = (n + 1) * ac->dt;
			int si2;
			for (si2 = 0; si2 < ac->nsrc; si2++)
				p[((size_t)ac->srccell[3 * si2 + 2] * ny
				   + ac->srccell[3 * si2 + 1]) * nx
				  + ac->srccell[3 * si2 + 0]] +=
					ac->srcgain[si2]
					* src_pulse(ac, t - ac->srcdelay[si2]);
		}

		/* ── 受音点記録 (サンプル n の時刻は (n+1)·dt — 1 サンプル未満の
		 * オフセットで、パルス幅 σ ≫ dt に対して無視できる) ── */
		for (r = 0; r < ac->nrecv; r++)
			ac->rec[(size_t)r * ac->nsteps + n] =
				p[((size_t)ac->recv[r].kc * ny + ac->recv[r].jc) * nx
				  + ac->recv[r].ic];

		/* ── 進捗 (mock_acoustic_solver.c と同じ "progress a/b") ── */
		{
			int step = (int)(((long long)(n + 1) * PROGRESS_TOTAL) / ac->nsteps);
			while (prog_done < step) {
				++prog_done;
				printf("progress %d/%d\n", prog_done, PROGRESS_TOTAL);
				fflush(stdout);
			}
		}
	}
	ac_log(ac, "solve: %d steps done", ac->nsteps);
	return 0;
}

void ac_free(ac_t *ac)
{
	free(ac->geom);  ac->geom = NULL;
	free(ac->gvox);  ac->gvox = NULL;
	free(ac->obs_has); ac->obs_has = NULL;
	free(ac->obs_bands); ac->obs_bands = NULL;
	free(ac->oface); ac->oface = NULL;
	free(ac->feedpos); ac->feedpos = NULL;
	free(ac->srccell); ac->srccell = NULL;
	free(ac->srcgain); ac->srcgain = NULL;
	free(ac->srcdelay); ac->srcdelay = NULL;
	free(ac->recv);  ac->recv = NULL;
	free(ac->freq1); ac->freq1 = NULL;
	free(ac->p);     ac->p = NULL;
	free(ac->vx);    ac->vx = NULL;
	free(ac->vy);    ac->vy = NULL;
	free(ac->vz);    ac->vz = NULL;
	free(ac->solid); ac->solid = NULL;
	free(ac->mvx);   ac->mvx = NULL;
	free(ac->mvy);   ac->mvy = NULL;
	free(ac->mvz);   ac->mvz = NULL;
	free(ac->rec);   ac->rec = NULL;
}
