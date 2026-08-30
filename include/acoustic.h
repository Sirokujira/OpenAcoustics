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
#define AC_NBAND_OFDX    6          /* .ofdx alpha[] のバンド数 (125 Hz .. 4 kHz) */
#define AC_BAND_F0       125.0      /* 最低バンドの中心周波数 [Hz] */
#define AC_TMIN          0.5        /* 計算時間 T = clamp(1.5 T_Sabine, 0.5, 3.0) [s] */
#define AC_TMAX          3.0
#define AC_SABINE_COEF   0.161      /* T_Sabine = 0.161 V/A */

/* 音源ごとのゲイン・遅延 (.ofdx acoustic.sources[] — ADR-0010 Decision 7)
 * の値域。範囲外は既定値に落とさず非零終了する (数値を捏造しない)。
 * 幾何音響側 (ga.h) と同じ値 — 両ソルバーで対称。 */
#define AC_GAIN_MAX   1000.0        /* |gain| の上限 (負は極性反転) */
#define AC_DELAY_MAX  1.0           /* delay_s の上限 [s] (下限は 0) */

#define AC_PATH_MAX  4096
#define AC_NAME_MAX  64
#define AC_PI        3.14159265358979323846
#define AC_SQRT2     1.41421356237309504880   /* オクターブバンドの端 (fc/√2, fc·√2) */

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
	/* 受音点 #1 用の別名 rir_<名前>.wav (契約の rir.wav に加えて出す。
	 * GUI の「フォルダから自動割当」は名前照合なので、#1 も名前付きの
	 * ファイルが無いと割り当てられない)。空なら出力しない。 */
	char   alias[AC_NAME_MAX + 16];
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
	double     srcx, srcy, srcz;       /* feed #1 の位置 (metadata の source 用) */
	double    *feedpos;                /* 全 feed の座標 (3 * nfeed) */
	int        nfeed;
	ac_recv_t *recv;
	int        nrecv;
	double    *freq1;                  /* frequency1 (記録のみ、物理には不使用) */
	int        nfreq1;

	/* .ofdx (吸音率、壁ごと) */
	int    have_ofdx;
	/* バンド別吸音率 (125 Hz .. 4 kHz)。読んだ値をそのまま保持する。 */
	double alpha_bands[AC_NWALL][AC_NBAND_OFDX];
	/* 実効吸音率 (境界インピーダンスに使う 1 値)。ac_setup が
	 * **有効帯域 [0, fmax] と重なるオクターブバンドだけ**を平均して決める
	 * (fmax = c/(10 dx) は格子で決まる上限)。低域担当のソルバーが 4 kHz の
	 * 吸音率まで混ぜて壁を決めるのは誤りで、幾何音響側 (バンド別に使う) との
	 * クロスオーバー整合も崩れるため。重なるバンドが無い (fmax < 88.4 Hz) 場合は
	 * 最低バンド 125 Hz の値を使う。 */
	double alpha[AC_NWALL];
	int    alpha_nband;                /* 実効値に使ったバンド数 */
	double alpha_band_hi;              /* 使ったバンドの上端 [Hz] (ログ/metadata) */
	/* 複数音源 (契約は OpenFDTD-X の ADR-0010) : .ofdx acoustic.multi_source
	 * (既定 false = feed #1 のみ = 従来動作)。true で全 feed に同一の
	 * ガウシアン微分パルスを注入し (共通 t0)、rir.wav は重ね合わせになる。
	 * 幾何音響側 (ofdx_acoustic_ga) と**対称に**実装する — ハイブリッド合成
	 * (ADR-0008) のバンドエネルギー整合は両ソルバーが同じ音源集合を使う前提。 */
	int    multi_source;

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
	int    isrc, jsrc, ksrc;           /* feed #1 のセル (metadata 用) */
	int    nsrc;                       /* 使用する音源数 (= multi ? nfeed : 1) */
	int   *srccell;                    /* 使用音源のセル (i,j,k の 3 * nsrc) */
	/* 音源ごとのゲイン・遅延 (ADR-0010 Decision 7) : .ofdx acoustic.sources[]。
	 * feed 順に nfeed 要素 (省略は gain = 1 / delay = 0 = 従来動作)。
	 * feed i には gain_i * s(t - delay_i) を注入する (s は共通パルス)。
	 * NULL のままなら全既定値 (setup が既定値で確保する)。 */
	double *srcgain;                   /* nfeed 要素 (既定 1) */
	double *srcdelay;                  /* nfeed 要素 [s] (既定 0) */

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

/* UTF-8 パスで開く fopen。Windows の fopen は ANSI コードページ解釈なので
 * 日本語を含むパス・ファイル名 (GUI の受音点名は既定で「マイク N」) が
 * 開けない。Windows だけ UTF-16 へ変換して _wfopen を使う。
 * 他 OS は fopen そのまま。 */
FILE *ac_fopen(const char *path, const char *mode);

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
