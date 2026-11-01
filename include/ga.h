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

/* 音源ごとのゲイン・遅延 (.ofdx acoustic.sources[] — ADR-0010 Decision 7)
 * の値域。範囲外は既定値に落とさず非零終了する (数値を捏造しない)。
 * FDTD 側 (acoustic.h) と同じ値 — 両ソルバーで対称。 */
#define GA_GAIN_MAX   1000.0   /* |gain| の上限 (負は極性反転) */
#define GA_DELAY_MAX  1.0      /* delay_s の上限 [s] (下限は 0) */

/* 1 次指向性の名前 -> a (D = a + (1-a) cos)。FDTD 側 (acoustic.h) と同じ値 */
#define GA_DIR_SUBCARDIOID   0.7
#define GA_DIR_CARDIOID      0.5
#define GA_DIR_SUPERCARDIOID 0.36602540378443865   /* (sqrt(3)-1)/2 : 前後比最大 */
#define GA_DIR_HYPERCARDIOID 0.25                  /* 指向係数最大 */

#define GA_PATH_MAX  4096
#define GA_NAME_MAX  64
#define GA_PI        3.14159265358979323846

/* 外壁インデックス (alpha[][] の並び — FDTD 側 AC_XM.. と同じ順序) */
enum { GA_XM = 0, GA_XP, GA_YM, GA_YP, GA_ZM, GA_ZP, GA_NWALL };

/* 鏡像法の探索コスト上限 (障害物が多いと面数^次数 で増えるため) */
#define GA_IMAGE_NODE_MAX 4000000

/* ── 入力要素 ──────────────────────────────────────────────────── */

/* 室内の障害物 (既定は剛体 alpha = 0)。遮蔽・鏡像法の反射面・光線追跡の
 * 反射に使う。本家 OpenFDTD の全 shape を**形状どおり**に扱う (ga_shape.c —
 * 内外判定は FDTD 側 fdtd.c の geom_inside と同じ規則):
 *   GA_SOLID_BOX    shape 1         軸平行の 6 面 (従来どおりの経路)
 *   GA_SOLID_POLY   shape 31〜33 (三角柱) / 41〜43 (角錐台)
 *                   平面の凸多面体。面 (5 / 6 枚) は一般の平面多角形として
 *                   鏡像法の面集合に入る
 *   GA_SOLID_CURVED shape 2 (楕円体) / 11〜13 (楕円柱) / 51〜53 (円錐台)
 *                   曲面は鏡面像を作れないので鏡像法の面集合に入らず、
 *                   反射は光線追跡が受け持つ (曲面に当たった経路はレイ側で
 *                   検出する — 二重計上の回避と取りこぼしの防止)
 * lo/hi は外接箱 (交差判定の前段の絞り込みと受音点の距離の下限に使う)。
 * 材質は .ofdx の acoustic.obstacles[] / acoustic.ga.obstacles[]。 */
#define GA_SOLID_BOX    1
#define GA_SOLID_POLY   2
#define GA_SOLID_CURVED 3
#define GA_POLY_FMAX    6       /* 多面体の面の最大数 (角錐台 = 6) */
typedef struct {
	int    shape;
	double g[8];
	double lo[3], hi[3];
	int    ok;          /* 0 = 未知 shape / 体積 0 (無視した) */
	int    kind;        /* GA_SOLID_BOX / GA_SOLID_POLY / GA_SOLID_CURVED */
	int    surf0;       /* この障害物の面の先頭 surf[] 添字 (無効なら -1)。
	                     * 箱は 2*axis + side の 6 面、多面体は面 f が surf0 + f、
	                     * 曲面は材質だけを持つ 1 面 (axis = GA_SURF_CURVED) */
	/* 多面体の面 (GA_SOLID_POLY のみ) : 外向き単位法線 fn と平面 fn.x = fd。
	 * 頂点 fv は fn のまわりに反時計回り。fnv = 0 は面積 0 の面 (角錐台の
	 * 尖った端など — 半空間としては使うが鏡像法には入れない)。 */
	int    nface;
	double fn[GA_POLY_FMAX][3];
	double fd[GA_POLY_FMAX];
	int    fnv[GA_POLY_FMAX];
	double fv[GA_POLY_FMAX][4][3];
	int    has_mat;     /* 1 = .ofdx で材質指定あり (下の mat_src のどちらか) */
	/* 材質の出所 : 0 = なし、1 = acoustic.obstacles (FDTD と共有するキー)、
	 * 2 = acoustic.ga.obstacles (幾何音響専用 — 1 より優先。既存入力の
	 * 後方互換のため。JSON の並び順に依らない)。 */
	int    mat_src;
	double mat_alpha[GA_NBAND];  /* バンド別吸音率 (既定 0 = 剛体) */
	double mat_scatter[GA_NBAND]; /* 散乱係数 (バンド別。[0] < 0 = 室の既定) */
} ga_geom_t;

/* 反射面。鏡像法 (一般化された面集合による鏡像) と光線追跡が共通で使う。
 *   axis >= 0 : 軸平行の有限矩形 (室の 6 面と箱の障害物の 6 面)。鏡映は
 *               1 座標の反転、面内判定は残り 2 座標の範囲判定、入射角は
 *               cos(theta) = |d[axis]| で済む (従来どおりの経路 — 箱だけの
 *               入力の結果はビット単位で不変)。
 *   axis = GA_SURF_PLANE  : 一般の平面凸多角形 (多面体の面)。法線 n・平面
 *               n.x = dpl・頂点 v[nv] で表し、鏡映・面内判定・入射角は
 *               ベクトル演算で行う。
 *   axis = GA_SURF_CURVED : 曲面の障害物の材質だけを持つ面 (鏡像法には
 *               入らない。法線は光線追跡が当たった点ごとに求める)。 */
#define GA_SURF_PLANE  (-1)
#define GA_SURF_CURVED (-2)
typedef struct {
	int    axis;             /* 法線の軸 : 0=x, 1=y, 2=z / GA_SURF_PLANE / _CURVED */
	double coord;            /* 平面の位置 (axis >= 0) */
	double nrm;              /* 音場側 (反射する側) を向く法線の符号 +1/-1 (axis >= 0) */
	double lo[2], hi[2];     /* 面内の矩形 (u = (axis+1)%3, v = (axis+2)%3) */
	double n[3];             /* 音場側を向く単位法線 (GA_SURF_PLANE) */
	double dpl;              /* 平面 n.x = dpl (GA_SURF_PLANE) */
	int    nv;               /* 多角形の頂点数 (GA_SURF_PLANE) */
	double v[4][3];          /* 頂点 (n のまわりに反時計回り) */
	int    img;              /* 1 = 鏡像法の面集合に入る */
	double alpha[GA_NBAND];
	double refl[GA_NBAND];   /* 垂直入射の圧力反射係数 sqrt(1-alpha) */
	double zeta[GA_NBAND];   /* 規格化インピーダンス (角度依存吸音のとき使う) */
	/* 散乱係数 (面ごと・バンド別)。レイの拡散/鏡面の抽選は 1 本のレイが
	 * 全バンドを運ぶため 1 回しかできない — 基準確率 sref で抽選し、
	 * 拡散枝はバンドエネルギーに s_b/sref、鏡面枝に (1-s_b)/(1-sref) を
	 * 掛ける (重み付き抽選)。期待値は各バンドで厳密に s_b / (1-s_b) になり、
	 * 鏡像側の sqrt(1-s_b) と過不足なく対応する (二重計上なし)。
	 * 全バンド同値 (suni = 1) なら重みは厳密に 1 で従来とビット等価。 */
	double scatter[GA_NBAND];
	double sref;             /* 抽選の基準確率 (一様なら scatter[0]、他は平均) */
	int    suni;             /* 1 = 全バンド同値 */
	int    wall;             /* 室壁なら GA_XM.. / 障害物面なら -1 */
	int    geom;             /* 障害物番号 (0 起点) / 室壁なら -1 */
} ga_surf_t;

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
	double     srcx, srcy, srcz;        /* feed #1 (metadata の source 用) */
	double    *feedpos;                 /* 全 feed の座標 (3 * nfeed) */
	int        nfeed;
	ga_recv_t *recv;
	int        nrecv;

	/* .ofdx */
	int    have_ofdx;
	int    have_band_alpha;              /* 吸音表を読めたか */
	double alpha[GA_NWALL][GA_NBAND];    /* バンド別吸音率 */
	double refl[GA_NWALL][GA_NBAND];     /* 圧力反射係数 R = sqrt(1-alpha) */
	double wall_scatter[GA_NWALL][GA_NBAND]; /* 面ごとの散乱係数
	                                     * (バンド別。[w][0] < 0 = 既定を使う) */

	/* 反射面リスト (室 6 面 + 障害物 6 面/個) */
	ga_surf_t *surf;
	int        nsurf;

	/* 幾何音響パラメータ */
	int    order;                        /* 鏡像法の次数 (1..3) */
	int    nrays;
	/* 複数音源 (契約は OpenFDTD-X の ADR-0010) : .ofdx acoustic.multi_source
	 * (既定 false = feed #1 のみ)。true で全 feed を強度 1 で t = 0 に
	 * 同時発火し、rir.wav は重ね合わせになる。1/N 正規化はしない。
	 * FDTD 側と**対称に**実装する (ハイブリッド合成のバンドエネルギー
	 * 整合 (ADR-0008) は両ソルバーが同じ音源集合を使う前提)。 */
	int    multi_source;
	int    nsrc;                         /* 使用する音源数 (= multi ? nfeed : 1) */
	/* 音源ごとのゲイン・遅延 (ADR-0010 Decision 7) : .ofdx acoustic.sources[]。
	 * feed 順に nfeed 要素 (省略は gain = 1 / delay = 0 = 従来動作)。
	 * 音源 i は t = srcdelay[i] に強度 srcgain[i] で発火する。
	 * NULL のままなら全既定値 (setup が既定値で確保する)。 */
	double *srcgain;                     /* nfeed 要素 (既定 1) */
	double *srcdelay;                    /* nfeed 要素 [s] (既定 0) */
	/* 音源の指向性 (.ofdx acoustic.feeds[].directivity / aim — OpenFDTD-X
	 * ADR-0014)。1 次の指向性 D_b(theta) = a_b + (1-a_b) cos(theta) (theta は
	 * aim からの放射方向の角度)、正面で D = 1 (無指向と同じ振幅) に規格化。
	 * 直接音・鏡像は放射方向の D_b を振幅に、レイは発射方向の D_b^2 を
	 * エネルギーに掛ける。NULL = 全 feed 無指向 (従来動作)。FDTD 側と同じ
	 * キー・同じ意味 (FDTD は有効帯域内のバンド平均の a を使う)。 */
	double *srcdir;                      /* nfeed * GA_NBAND : バンド別の a (既定 1) */
	double *srcaim;                      /* nfeed * 3 : 向き (setup で正規化。0 = 未指定) */
	double  delay_max;                   /* 発火する音源集合の max(delay) */
	double scatter[GA_NBAND];            /* 拡散反射の割合の室既定 (バンド別) */
	int    angle_dep;                    /* 1 = 角度依存吸音 (局所反応) を使う */
	/* acoustic.ga.angle_dependent_absorption を読んだか。acoustic 直下の
	 * 同名キー (FDTD 側と共有する対称キー) より ga 側を優先するため。 */
	int    angle_dep_ga;
	long   qidx;                         /* 決定的ハッシュ列の位置 */
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
	double vol, area;                    /* 室容積 / 表面積 */

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

/* 局所反応境界 (実数の規格化インピーダンス zeta) の統計 (ランダム) 入射吸音率。
 * Paris の式 alpha_stat = 2 int_0^{pi/2} alpha(theta) cos sin dtheta を
 * R(theta) = (zeta cos - 1)/(zeta cos + 1) について積分した閉形式 :
 *   alpha_stat(zeta) = (8/zeta^2) [ zeta + 1 - 2 ln(1+zeta) - 1/(1+zeta) ]
 * 角度依存吸音では、吸音表の値 (= ランダム入射) からこの式を逆に解いて
 * zeta を決める。逆に解かずに表の値を垂直入射として使うと二重に効いてしまう。 */
double ga_alpha_stat(double zeta);
double ga_zeta_from_alpha(double alpha_stat);   /* 逆問題 (zeta >= peak の枝) */
#define GA_ALPHA_STAT_MAX 0.9514   /* 局所反応・実インピーダンスの上限 (概数) */

/* ga_trace.c : 鏡像法 (早期反射) と光線追跡 (後期残響)
 *   ga_images : 音源 si (feedpos の添字) から受音点 r への可視な鏡像を
 *               band[] へ置く。受音点統計 (nimage/nblocked) は加算するだけ
 *               なので、リセットとまとめログは呼び出し側 (ga_synth)。
 *   ga_rays   : 全音源 x 全受音点のエコーグラム echo[] を埋める (進捗はここ) */
/* ── 障害物の厳密形状 (ga_shape.c) ──
 *   ga_solid_prepare : shape から kind・外接箱・多面体の面を作る (0 = 未知 /
 *                      体積 0)
 *   ga_solid_inside  : 点が形状の内側 (境界を含む) か — FDTD 側の
 *                      geom_inside と同じ規則
 *   ga_solid_entry   : 直線 p + t d が t in (tmin, tmax) で形状に入る最初の t。
 *                      入れば 1 を返し、外向き単位法線 nrm と面番号 face
 *                      (多面体の面 0..nface-1 / 箱は 2*axis + side / 曲面は -1)
 *   ga_solid_distance: 点から形状までの距離の推定 (受音球の半径の上限用) */
int    ga_solid_prepare(ga_geom_t *o);
int    ga_solid_inside(const ga_geom_t *o, double x, double y, double z);
int    ga_solid_entry(const ga_geom_t *o, const double p[3], const double d[3],
                      double tmin, double tmax, double *t, double nrm[3],
                      int *face);
double ga_solid_distance(const ga_geom_t *o, double x, double y, double z);

int ga_images(ga_t *g, int r, int si, double *band);
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
