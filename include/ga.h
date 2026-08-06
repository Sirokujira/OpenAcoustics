/* ga.h — ofdx_acoustic_ga 共有定義 (幾何音響ソルバーのコンテキストと各段の入口)
 *
 * OpenFDTD-X の室内音響用「高域担当」ソルバー (ADR-0007 の出力契約に準拠)。
 * ofdx_acoustic_fdtd が低域 (fmax = c/(10 dx)) を、こちらが高域を受け持ち、
 * 帯域分割・合成は GUI (OpenFDTD-X) 側が行う。したがって本ソルバーは
 * 「高域の RIR を規約どおりに出す」ことだけに責任を持つ。
 *
 * 設計規約は ofdx_acoustic_fdtd と共有 (CLAUDE.md / AGENTS.md) :
 *   - C11、外部ライブラリ依存なし。WAV は src/wav.c を共用する。
 *   - グローバル変数禁止。状態は ga_t コンテキスト 1 個を main で確保して渡す。
 *   - C99 VLA 禁止、場・信号はすべて double、WAV 書き出し時のみ float32。
 *   - 乱数不使用。レイ方向は球面フィボナッチ格子 (決定的準一様分布)。
 *
 * FDTD 側のコードには一切手を入れない (別バイナリとして独立)。共有するのは
 * src/wav.c だけで、そのために必要な ac_fopen は ga_main.c が定義する。
 */
#ifndef GA_H
#define GA_H

#include <stdio.h>

/* ── 識別情報 (metadata.json に出力) ────────────────────────────── */
#define GA_SOLVER_NAME       "ofdx_acoustic_ga"
#define GA_SOLVER_VERSION    "1.0.0"
#define GA_CONTRACT_VERSION  1

/* ── 物理定数 (20 ℃ 空気 — FDTD 側と同一) ──────────────────────── */
#define GA_C0    343.0        /* 音速 [m/s] */
#define GA_RHO0  1.204        /* 密度 [kg/m^3] */

/* ── 出力サンプリング周波数 (固定) ─────────────────────────────── */
#define GA_FS    48000

/* ── オクターブバンド (125 Hz 〜 8 kHz の 7 バンド) ──────────────
 * .ofdx の alpha[] は 6 バンド (125〜4000 Hz) なので、8 kHz バンドは
 * 4 kHz の値を外挿して使う (solver.log に明示する)。 */
#define GA_NBAND       7
#define GA_NBAND_OFDX  6
#define GA_BAND_F0     125.0    /* 最低バンドの中心周波数 */

/* ── 幾何音響のパラメータ既定値 (.ofdx acoustic.ga で上書き可) ──── */
#define GA_ORDER_DEFAULT 2      /* 鏡像法の次数 */
#define GA_ORDER_MIN     1
#define GA_ORDER_MAX     3
#define GA_RAYS_DEFAULT  30000  /* 光線追跡の本数 */
#define GA_RAYS_MIN      100
#define GA_RAYS_MAX      2000000
#define GA_ALPHA_DEFAULT 0.1    /* 吸音表が無い/該当 role が無い壁の吸音率 */
#define GA_TMIN          0.5    /* 計算時間 T = clamp(1.5 T_Eyring, 0.5, 3.0) [s] */
#define GA_TMAX          3.0
#define GA_BIN_S         0.001  /* エコーグラムのビン幅 [s] */
#define GA_ENERGY_FLOOR  1e-8   /* レイの打ち切り (-80 dB) */
#define GA_RSPHERE_MIN   0.3    /* 受音検出球の半径の下限 / 上限 [m] */
#define GA_RSPHERE_MAX   1.5
#define GA_SCATTER_DEFAULT 0.1  /* 拡散反射 (Lambert) の割合 s : 0 = 鏡面のみ */

/* 空気吸収 (ISO 9613-1) の既定条件 */
#define GA_TEMP_C_DEFAULT     20.0
#define GA_HUMID_DEFAULT      50.0
#define GA_PRESS_KPA_DEFAULT  101.325

#define GA_PATH_MAX  4096
#define GA_NAME_MAX  64
#define GA_PI        3.14159265358979323846

/* 外壁インデックス (alpha[][] の並び — FDTD 側 AC_XM.. と同じ順序) */
enum { GA_XM = 0, GA_XP, GA_YM, GA_YP, GA_ZM, GA_ZP, GA_NWALL };

/* ── 入力要素 ──────────────────────────────────────────────────── */

/* 室内の剛体障害物。可視性判定 (遮蔽) と光線追跡の反射に使う。
 * shape 1 (直方体) は lo/hi が厳密に形状と一致し、他の shape は AABB 近似。 */
typedef struct {
	int    shape;
	double g[8];
	double lo[3], hi[3];
	int    ok;          /* 0 = 未知 shape (無視した) */
	int    exact;       /* 1 = shape 1 (厳密)、0 = AABB 近似 */
} ga_geom_t;

typedef struct {
	char   name[GA_NAME_MAX];        /* .ofd point 行末の "# 名前" (空可) */
	char   file[GA_NAME_MAX + 16];   /* 出力 WAV ファイル名 */
	char   alias[GA_NAME_MAX + 16];  /* 受音点 #1 の別名 (GUI の名前照合用) */
	double x, y, z;
	double radius;                   /* 検出球半径 [m] */
	/* 統計 (ログ用) */
	int    nimage;                   /* 可視だった鏡像音源の数 */
	int    nblocked;                 /* 遮蔽で棄却した鏡像音源の数 */
	long   nhit;                     /* 検出球を通過したレイの数 */
} ga_recv_t;

/* ── コンテキスト (唯一の状態。main で 1 個確保) ────────────────── */
typedef struct {
	/* パスとログ */
	char  workdir[GA_PATH_MAX];
	char  ofd_path[GA_PATH_MAX];
	char  ofd_name[GA_NAME_MAX * 4];
	char  ofdx_path[GA_PATH_MAX + 8];
	FILE *logfp;

	/* .ofd 入力 */
	char       title[256];
	double     x0, x1, y0, y1, z0, z1;  /* 室 (直方体) */
	double     dxmin;                   /* メッシュ最小刻み (受音点名の照合許容) */
	ga_geom_t *geom;
	int        ngeom;
	double     srcx, srcy, srcz;
	int        nfeed;
	ga_recv_t *recv;
	int        nrecv;

	/* .ofdx */
	int    have_ofdx;
	int    have_band_alpha;              /* 吸音表を読めたか */
	double alpha[GA_NWALL][GA_NBAND];    /* バンド別吸音率 */
	double refl[GA_NWALL][GA_NBAND];     /* 圧力反射係数 R = sqrt(1-alpha) */

	/* 幾何音響パラメータ */
	int    order;                        /* 鏡像法の次数 (1..3) */
	int    nrays;
	double scatter;                      /* 拡散反射の割合 (0..1) */
	long   qidx;                         /* 準乱数列 (Halton) の位置 — 決定的 */
	double temp_c, humid, press_kpa;
	int    air_on;
	double duration_user;                /* > 0 なら計算時間を強制 */
	double rsphere_user;                 /* > 0 なら検出球半径を強制 */

	/* 空気吸収 (ISO 9613-1) */
	double air_db_m[GA_NBAND];           /* 減衰 [dB/m] (音圧レベル) */
	double air_e[GA_NBAND];              /* エネルギー減衰係数 m [1/m] */

	/* 時間軸 */
	int    nsamples;
	double duration;
	double teyring[GA_NBAND];            /* Eyring 残響時間 (無損失は -1) */
	double tsabine[GA_NBAND];            /* Sabine 残響時間 (参考値、同上) */
	double t60max;                       /* Eyring のバンド最大 (無損失は -1) */
	double tsabmax;                      /* Sabine のバンド最大 (無損失は -1) */
	double flo, fhi;                     /* 有効帯域 [Hz] */
	double vol, surf;                    /* 室容積 / 表面積 */

	/* エコーグラム (nrecv * nbins * GA_NBAND) */
	int     nbins;
	double *echo;

	/* 出力 RIR (nrecv * nsamples) */
	double *rir;

	/* 合成の作業領域 */
	double *band;        /* GA_NBAND * nsamples : バンド別インパルス列 */
	int     fftn;
	double *fr, *fi;     /* FFT 作業 (fftn) */
	double *yr, *yi;     /* スペクトル累算 (fftn) */
	double *twr, *twi;   /* 回転因子 (fftn/2) */
	int    *brev;        /* ビット反転表 (fftn) */

	/* 統計 */
	long   nbounce;
	long   nray_detect;
} ga_t;

/* ── 各段の入口 ────────────────────────────────────────────────── */

/* ga_main.c : ログ (solver.log と、エラー時は stderr にも) */
void ga_log(ga_t *g, const char *fmt, ...);
void ga_err(ga_t *g, const char *fmt, ...);

/* UTF-8 パスで開く fopen。宣言は acoustic.h と同一 —
 * src/wav.c を FDTD 側と共用するため、名前も acoustic.h に合わせてある
 * (定義は ga_main.c。FDTD の main.c とは別バイナリなので衝突しない)。 */
FILE *ac_fopen(const char *path, const char *mode);

/* src/wav.c (FDTD 側と共用) : float32 モノラル WAV */
int ac_write_wav_f32(const char *path, int fs, const double *x, int n);

/* ga_input.c */
int ga_find_input(ga_t *g);
int ga_read_ofd(ga_t *g);

/* ga_ofdx.c */
int ga_read_ofdx(ga_t *g);

/* ga_setup.c */
int    ga_setup(ga_t *g);
double ga_air_alpha_db_m(double f_hz, double temp_c, double humid_pct,
                         double press_kpa);   /* ISO 9613-1 [dB/m] */

/* ga_trace.c : 鏡像法 (早期反射) と光線追跡 (後期残響)
 *   ga_images : 受音点 r の可視な鏡像音源を band[] へ置く
 *   ga_rays   : 全受音点のエコーグラム echo[] を 1 パスで埋める (進捗はここ) */
int ga_images(ga_t *g, int r, double *band);
int ga_rays(ga_t *g);

/* ga_synth.c : バンド別インパルス列 + エコーグラム -> RIR
 *   ga_deposit : 反射 1 発 (時刻 t [s]、バンド別振幅 amp[GA_NBAND]) を
 *                band[] へ 3 次ラグランジュ補間で置く (DC 利得・1 次モーメント厳密) */
void ga_deposit(ga_t *g, double *band, double t, const double *amp);
int  ga_synth(ga_t *g);

/* ga_json.c */
int ga_write_metadata(ga_t *g);
int ga_write_metrics(ga_t *g);

/* ga_main.c */
void ga_free(ga_t *g);

#endif /* GA_H */
