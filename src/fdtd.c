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
 *   法線速度を 0 に固定する (面マスク)。shape 1 (直方体) は厳密、他の shape は
 *   AABB (軸平行境界箱) 近似で warning を出す。
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

/* 点を含む剛体形状の番号 (0 起点。無ければ -1)。エラーメッセージ用 —
 * ボクセル化と同じ AABB 判定を使う (形状の近似も含めて実際に剛体になった
 * 領域と一致させるため)。 */
static int geom_containing(const ac_t *ac, double x, double y, double z)
{
	double lo[3], hi[3];
	int n;
	for (n = 0; n < ac->ngeom; n++) {
		if (!geom_aabb(&ac->geom[n], lo, hi)) continue;
		if (x >= lo[0] && x <= hi[0] && y >= lo[1] && y <= hi[1] &&
		    z >= lo[2] && z <= hi[2])
			return n;
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

int ac_setup(ac_t *ac)
{
	const double c = AC_C0, rho = AC_RHO0;
	double lx, ly, lz, vol, area, tsolve;
	long long cells;
	size_t np, nvx, nvy, nvz;
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

	/* ── 計算時間 : T = clamp(1.5·T_Sabine, 0.5, 3.0)
	 * T_Sabine = 0.161·V/A。A は吸音表 (壁 6 面 × 帯域平均 α) から。
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

	/* ── 境界係数 ── */
	for (w = 0; w < AC_NWALL; w++) {
		double a = ac->alpha[w];
		if (a < 1e-6) {
			ac->wrigid[w] = 1;   /* 剛壁 : 更新なし (v ≡ 0) — 厳密に無損失 */
			ac->wZ[w] = 0.0;
			ac->wA[w] = ac->wB[w] = 0.0;
		}
		else {
			double s = sqrt(1.0 - ((a > 1.0) ? 1.0 : a));
			double Z = rho * c * (1.0 + s) / (1.0 - s + 1e-300);
			double ma = rho * ac->dx / (2.0 * ac->dt);
			ac->wrigid[w] = 0;
			ac->wZ[w] = Z;
			ac->wA[w] = (ma - Z / 2.0) / (ma + Z / 2.0);
			ac->wB[w] = 1.0 / (ma + Z / 2.0);
		}
		ac_log(ac, "wall %d (%s): alpha = %.4g -> %s", w,
		       (w == AC_XM) ? "x-" : (w == AC_XP) ? "x+" :
		       (w == AC_YM) ? "y-" : (w == AC_YP) ? "y+" :
		       (w == AC_ZM) ? "z-" : "z+",
		       ac->alpha[w],
		       ac->wrigid[w] ? "rigid" : "impedance");
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

	/* ── 剛体ボクセル化 (セル中心が形状内 → 剛体) ── */
	for (n = 0; n < ac->ngeom; n++) {
		double lo[3], hi[3];
		int i0, i1, j0, j1, k0, k1;
		if (!geom_aabb(&ac->geom[n], lo, hi)) {
			ac_log(ac, "warning: geometry #%d has unknown shape %d — ignored",
			       n + 1, ac->geom[n].shape);
			continue;
		}
		if (ac->geom[n].shape != 1)
			ac_log(ac, "warning: geometry #%d shape %d approximated by its "
			       "AABB (only shape 1 = box is voxelized exactly in v1)",
			       n + 1, ac->geom[n].shape);
		/* AABB とセル中心格子の重なり範囲だけ走査 */
		i0 = (int)ceil((lo[0] - ac->x0) / ac->dx - 0.5);
		i1 = (int)floor((hi[0] - ac->x0) / ac->dx - 0.5);
		j0 = (int)ceil((lo[1] - ac->y0) / ac->dx - 0.5);
		j1 = (int)floor((hi[1] - ac->y0) / ac->dx - 0.5);
		k0 = (int)ceil((lo[2] - ac->z0) / ac->dx - 0.5);
		k1 = (int)floor((hi[2] - ac->z0) / ac->dx - 0.5);
		if (i0 < 0) i0 = 0;
		if (j0 < 0) j0 = 0;
		if (k0 < 0) k0 = 0;
		if (i1 > nx - 1) i1 = nx - 1;
		if (j1 > ny - 1) j1 = ny - 1;
		if (k1 > nz - 1) k1 = nz - 1;
		for (k = k0; k <= k1; k++)
			for (j = j0; j <= j1; j++)
				for (i = i0; i <= i1; i++)
					ac->solid[((size_t)k * ny + j) * nx + i] = 1;
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
			int gi = geom_containing(ac, sx, sy, sz);
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
			int gi = geom_containing(ac, r->x, r->y, r->z);
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
