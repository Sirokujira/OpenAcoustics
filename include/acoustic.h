/* acoustic.h — ofdx_acoustic_fdtd 共有定義 (コンテキスト構造体と各段の入口)
 *
 * OpenFDTD-X の室内音響用外部 FDTD ソルバー (ADR-0007 の出力契約に準拠)。
 * OpenPEEC と同じ設計規約: グローバル変数禁止、状態は ac_t コンテキスト
 * 構造体 1 個を main で確保して各関数に渡す。C11、外部ライブラリ依存なし。
 */
#ifndef ACOUSTIC_H
#define ACOUSTIC_H

#include <stdio.h>

/* ── 識別情報 (metadata.json に出力) ────────────────────────────── */
#define AC_SOLVER_NAME       "ofdx_acoustic_fdtd"
#define AC_SOLVER_VERSION    "1.0.0"
#define AC_CONTRACT_VERSION  1

/* ── 物理定数 (20 ℃ 空気) ──────────────────────────────────────── */
#define AC_C0    343.0        /* 音速 [m/s] */
#define AC_RHO0  1.204        /* 密度 [kg/m^3] */

/* ── 数値仕様 (v1) ─────────────────────────────────────────────── */
#define AC_CFL           0.99        /* Courant 数の上限 (c·dt·sqrt(3)/dx <= 0.99) */
#define AC_MAX_CELLS     30000000    /* セル総数の上限 (超えたら非零終了) */
#define AC_ALPHA_DEFAULT 0.1        /* 吸音表が無い/該当 role が無い壁の吸音率 */
#define AC_TMIN          0.5        /* 計算時間 T = clamp(1.5 T_Sabine, 0.5, 3.0) [s] */
#define AC_TMAX          3.0
#define AC_SABINE_COEF   0.161      /* T_Sabine = 0.161 V/A */

#define AC_PATH_MAX  4096
#define AC_NAME_MAX  64
#define AC_PI        3.14159265358979323846

/* 外壁インデックス (alpha[] / 境界係数の並び) */
enum { AC_XM = 0, AC_XP, AC_YM, AC_YP, AC_ZM, AC_ZP, AC_NWALL };

/* ── 入力要素 ──────────────────────────────────────────────────── */

typedef struct {
	int    shape;      /* .ofd の shape (1=直方体は厳密、他は AABB 近似) */
	double g[8];       /* shape パラメータ (本家 sol/ingeometry.c と同一) */
} ac_geom_t;

typedef struct {
	char   name[AC_NAME_MAX];   /* .ofd point 行の "# 名前" (空可) */
	char   file[AC_NAME_MAX + 16]; /* 出力 WAV ファイル名 (setup で決定) */
	double x, y, z;             /* 指定座標 [m] */
	int    ic, jc, kc;          /* スナップ後のセル (setup で決定) */
} ac_recv_t;

/* ── コンテキスト (唯一の状態。main で 1 個確保) ────────────────── */
typedef struct {
	/* パスとログ */
	char  workdir[AC_PATH_MAX];
	char  ofd_path[AC_PATH_MAX];
	char  ofd_name[AC_NAME_MAX * 4];  /* metadata 用のファイル名 */
	char  ofdx_path[AC_PATH_MAX + 8];
	FILE *logfp;                       /* solver.log (main が開閉) */

	/* .ofd 入力 */
	char       title[256];
	double     x0, x1, y0, y1, z0, z1; /* 計算領域 [m] */
	double     dxmin;                  /* メッシュ最小刻み [m] */
	int        nonuniform;             /* 非一様メッシュだったか (警告用) */
	ac_geom_t *geom;
	int        ngeom;
	double     srcx, srcy, srcz;       /* feed #1 の位置 */
	int        nfeed;
	ac_recv_t *recv;
	int        nrecv;
	double    *freq1;                  /* frequency1 (記録のみ、物理には不使用) */
	int        nfreq1;

	/* .ofdx (吸音率、壁ごと。帯域平均済み) */
	int    have_ofdx;
	double alpha[AC_NWALL];

	/* 格子と時間 */
	double dx;
	int    nx, ny, nz;
	int    fs;            /* サンプリング周波数 [Hz] (= 1/dt, CFL を満たす最小整数) */
	double dt;
	int    nsteps;
	double duration;      /* = nsteps * dt */
	double tsab;          /* Sabine 残響時間 (A=0 のとき -1) */
	double fmax;          /* 音源設計帯域 = c/(10 dx) */
	double sigma, t0;     /* ガウシアン微分パルスの幅と遅延 */
	int    isrc, jsrc, ksrc;

	/* 境界係数 (半陰的局所インピーダンス更新 v+ = wA v- ± wB p) */
	int    wrigid[AC_NWALL];  /* 1 = 剛壁 (更新しない) */
	double wZ[AC_NWALL];      /* インピーダンス [Pa s/m] (剛壁は 0) */
	double wA[AC_NWALL], wB[AC_NWALL];

	/* 場 (フラット配列。VLA 禁止 — OpenPEEC portability.md) */
	double        *p;              /* nx*ny*nz (セル中心) */
	double        *vx;             /* (nx+1)*ny*nz (x 面) */
	double        *vy;             /* nx*(ny+1)*nz (y 面) */
	double        *vz;             /* nx*ny*(nz+1) (z 面) */
	unsigned char *solid;          /* nx*ny*nz : 1 = 剛体セル */
	unsigned char *mvx, *mvy, *mvz; /* 面マスク : 1 = 更新する (両側とも流体) */

	/* 受音点の記録 (nrecv * nsteps) */
	double *rec;
} ac_t;

/* ── 各段の入口 ────────────────────────────────────────────────── */

/* main.c : ログ (solver.log とエラー時は stderr へ) */
void ac_log(ac_t *ac, const char *fmt, ...);
void ac_err(ac_t *ac, const char *fmt, ...);

/* input_ofd.c */
int ac_find_input(ac_t *ac);              /* workdir から唯一の .ofd を探す */
int ac_read_ofd(ac_t *ac);

/* input_ofdx.c */
int ac_read_ofdx(ac_t *ac);               /* 無ければ既定値のまま 0 を返す */

/* fdtd.c */
int  ac_setup(ac_t *ac);
int  ac_run(ac_t *ac);                    /* stdout に "progress a/b" を出す */
void ac_free(ac_t *ac);

/* wav.c : float32 モノラル WAV (44 byte ヘッダ、WAVE_FORMAT_IEEE_FLOAT) */
int ac_write_wav_f32(const char *path, int fs, const double *x, int n);

/* jsonout.c */
int ac_write_metadata(ac_t *ac);
int ac_write_metrics(ac_t *ac);

#endif /* ACOUSTIC_H */
