# OpenAcoustics — ofdx_acoustic_fdtd / ofdx_acoustic_ga

OpenFDTD-X (GUI) の室内音響用外部ソルバー (C11)。部屋のインパルス応答
(RIR) を計算し、ADR-0007 のファイル契約 (rir.wav / metadata.json /
metrics.json / solver.log) で GUI に返す。**バイナリは 2 本**ある。

| バイナリ | 担当 | 手法 |
|---|---|---|
| `ofdx_acoustic_fdtd` | 低域 | 3 次元線形音響のスタガード格子 FDTD |
| `ofdx_acoustic_ga` | 高域 | 幾何音響 (鏡像法 + 光線追跡) |

FDTD の上限は fmax = c/(10·dx) で、ホール規模で 4 kHz を得るのは原理的に
不可能なため、**低域 = FDTD / 高域 = 幾何音響**のハイブリッドにしている。
**帯域分割・両者の合成は GUI (OpenFDTD-X) 側の担当**で、ここではやらない。
2 本は独立したバイナリで、共有するソースは `src/wav.c` だけ。
**FDTD 側のコードに幾何音響の都合で手を入れないこと** (逆も同じ)。

OpenPEEC の姉妹プロジェクトで、ビルド規約・移植性規則を共有する。
外部ライブラリに依存しない (C11 + CMake、OpenMP のみ任意。
WAV / JSON / FFT は自前実装 — 追加しないこと)。

**`AGENTS.md`** に同じ規約を単独で読める形でまとめてある (Codex 等、
`CLAUDE.md` を読まないエージェント向け)。**規約を変えたら両方直すこと。**

## ビルド / テスト (変更後は必ず実行 — 全判定 PASS が完了条件)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"        # bin/ に 2 本できる

# 検証 : 解析解との比較 130 判定 (FDTD 50 + 幾何音響 80)。1 本のスクリプトが
# 両ソルバーを判定する (第 3 引数を省略すると実行ファイル名の fdtd -> ga
# 置換で幾何音響側を探すので、CI の呼び出しは変えなくてよい)
sh data/sample/acoustic_check.sh "$PWD/bin/ofdx_acoustic_fdtd" /tmp/ac-check

# E2E (ホール規模 20x15x10 m。同じ .ofd を両ソルバーが読める)
mkdir -p /tmp/ac-hall && cp data/sample/hall.ofd data/sample/hall.ofdx /tmp/ac-hall/
./bin/ofdx_acoustic_fdtd /tmp/ac-hall && grep "normal end" /tmp/ac-hall/solver.log
mkdir -p /tmp/ga-hall && cp data/sample/hall.ofd data/sample/hall.ofdx /tmp/ga-hall/
./bin/ofdx_acoustic_ga /tmp/ga-hall && grep "normal end" /tmp/ga-hall/solver.log
```

「ビルドが通る」は合格条件ではない。**130 判定すべて OK でなければ完了ではない**
(この検証群が物理と契約の番人になっている)。バイナリは `bin/` に出る
(CMakeLists が固定)。sanitize ビルド (CI と同じ ASan+UBSan) も `bin/` を
上書きするので、検証後は Release を作り直すこと。

## ソース構成

**ofdx_acoustic_fdtd (低域)**

| ファイル | 役割 |
|---|---|
| `include/acoustic.h` | コンテキスト構造体 `ac_t` と全定数 (c, ρ, CFL, 上限) |
| `src/main.c` | エントリ。引数 → 入力 → setup → run → 出力。ログ (`ac_log`/`ac_err`) |
| `src/input_ofd.c` | `.ofd` パーサ (OpenFDTD 互換、未知キーは無視) と入力探索 |
| `src/input_ofdx.c` | `.ofdx` (JSON サイドカー) の tolerant scanner — 必要キーのみ読む |
| `src/fdtd.c` | setup (格子/音源設計/境界係数/ボクセル化) と leapfrog 本体 |
| `src/jsonout.c` | metadata.json / metrics.json |

**ofdx_acoustic_ga (高域)**

| ファイル | 役割 |
|---|---|
| `include/ga.h` | コンテキスト構造体 `ga_t` と全定数 (バンド, fs, 既定値) |
| `src/ga_main.c` | エントリ。ログ (`ga_log`/`ga_err`) と `ac_fopen` の定義 |
| `src/ga_input.c` | `.ofd` パーサ (幾何音響が使う部分集合) と入力探索 |
| `src/ga_ofdx.c` | `.ofdx` の tolerant scanner (バンド別 α と `acoustic.ga`) |
| `src/ga_setup.c` | 室/反射面集合/境界 (Paris の逆解)/ISO 9613-1/Eyring/確保 |
| `src/ga_trace.c` | 一般化鏡像法 (面集合・可視性判定) と光線追跡 (エコーグラム) |
| `src/ga_synth.c` | 分数遅延の反射配置、バンド重み、自前 radix-2 FFT、RIR 合成 |
| `src/ga_json.c` | metadata.json / metrics.json |

**共有**

| ファイル | 役割 |
|---|---|
| `src/wav.c` | float32 モノラル WAV 書き出し (44 byte ヘッダ)。両ターゲットに入る。必要な `ac_fopen` は `main.c` と `ga_main.c` がそれぞれ定義する (別バイナリなので衝突しない) |

## FDTD (ofdx_acoustic_fdtd) の物理不変条件 (壊すと結果が静かに狂う)

1. **CFL**: fs = ceil(c·√3/(0.99·dx)) — CFL を満たす**最小の整数** Hz。
   変えると box_modes の fs 判定 (独立に同式で計算) が落ちる。
2. **剛壁 (α=0) は厳密に無損失**: 音源はガウシアン微分 (DC 成分 0) なので
   離散系はエネルギー保存する。境界更新に散逸を持ち込むと (d) が落ちる。
3. **インピーダンス境界 R = √(1−α)**: Z = ρc(1+√(1−α))/(1−√(1−α))。
   α→0 は剛壁に厳密一致、α→1 は Z→ρc (無反射)。tube ケース (±3%) が番人。
4. **決定性**: 乱数不使用。OpenMP 並列 (p 更新・v 更新のセルループ) は
   リダクションを持たず、スレッド数によらず**ビット単位で一致** (OpenPEEC と
   同じ原則)。リダクションを持つ並列化を足すと (e) が落ちる。その場合は
   ReadMe の「ビット一致」の主張ごと見直すこと。
5. **契約 (ADR-0007)**: 出力ファイル名・WAV ヘッダ・`progress a/b` 書式・
   metadata キーは OpenFDTD-X の
   `tests/acoustics/mock_acoustic_solver.c` (参照実装) と
   `docs/adr/0007-acoustic-solver-contract.md` が正。
   **契約を変えるときは OpenFDTD-X 側の mock と ADR-0007 に追随する**
   (こちらだけ先に変えない)。metadata.json はキー追加のみ可 (削除・改名禁止)。
6. **数値を捏造しない**: 入力が読めない・feed/point が無い・セル総数 > 3000 万
   は非零終了 + stderr に理由。合成 RIR を出して正常終了しない。
7. **複数音源 (multi_source) は両ソルバーで対称**: 契約は OpenFDTD-X の
   ADR-0010。`.ofdx` の `acoustic.multi_source` (既定 false = feed #1 のみ =
   従来動作) を両ソルバーが同じ意味で読む。true では全 feed に同一パルスを
   注入した重ね合わせ (1/N 正規化なし)。音源ごとのゲイン・遅延は
   `acoustic.sources[]` (同 Decision 7 : gain, delay_s — 既定 1 / 0、
   範囲外は非零終了) を両ソルバーが同じ意味で読み、音源 i を t = delay_i に
   強度 gain_i で発火させる (計算時間は max(delay) だけ自動延長)。
   ハイブリッド合成 (ADR-0008) のバンドエネルギー整合は両ソルバーが同じ
   音源集合を使う前提なので、**片方だけ変えてはいけない**。
   番人は (h)(h2) 線形性判定と (K)(K2)。

## 幾何音響 (ofdx_acoustic_ga) の不変条件 — 番人は acoustic_check.sh の (A)〜(G)

1. **振幅規約**: 自由音場の直接音が **1/(4πr)**。離散 RIR ではこれを
   「到達の標本和 (単位標本利得)」として実現する。反射は 3 次ラグランジュ
   分数遅延 (4 タップ) で置くので標本和 = 振幅・1 次モーメント = 到達時刻が
   厳密。**FDTD 側と合成するときの共通規約なので変えてはいけない** ((A)(B))。
2. **時間原点**: **t = 0 が音源発火時刻**。直接音は t = r/c。metadata の
   `t0_s` は常に 0。遅延を入れると GUI 側の合成が壊れる ((A)(B)(G))。
3. **出力 fs = 48000 Hz 固定**。metadata の `sample_rate` と WAV ヘッダの
   両方に出る ((G))。
4. **バンド重みは単位分割**: Σ_b W_b(f) = 1 かつ W_b(fc_b') = δ。前者が
   無損失ケースのエネルギー保存 ((E)) を、後者がバンド中心の DFT 読み取り
   ((D)) を支えている。零位相 (実数重み) なので群遅延 0 — ここに位相を
   持ち込むと (2) が壊れる。
5. **決定性**: 乱数生成器を使わない。レイ方向は球面フィボナッチ格子、
   拡散反射の抽選は固定の 32 bit 整数ハッシュ、後期残響の符号列は固定
   初期値の LFSR。**OpenMP は使わない** (エコーグラム加算がリダクションに
   なりビット一致を壊すため)。並列化するならビット一致の主張ごと見直すこと ((F))。
   準乱数列 (Halton) を拡散反射に使ってはいけない — 直方体の反射列と共鳴して
   平均自由行程が 4V/S から +4.8% ずれ、残響時間がそのまま狂う (実測済み)。
6. **二重計上の禁止**: 鏡像法が受け持つのは「次数 order 以下・一度も拡散して
   いない鏡面経路」だけで、レイ側はその条件を外れた経路のみ検出する。
   障害物の面も鏡像法の反射面集合に入っているので、条件に障害物は入らない。
   散乱係数 s_b (面ごと・バンド別) のぶんは鏡像の振幅からバンドごとに
   sqrt(1-s_b) で抜き、抜けた分をレイが受け持つ。レイの拡散/鏡面の抽選は
   面ごとの基準確率 s* (バンド平均) で 1 回だけ行い、拡散枝はバンド
   エネルギーに s_b/s*、鏡面枝に (1-s_b)/(1-s*) を掛ける (重み付き抽選 —
   期待値は各バンドで厳密、全バンド同値なら重みは 1 でビット等価)。
   どちらかを変えるなら両方直すこと。番人は (L) の閉形式 DFT 判定。
7. **鏡像法は面集合 surf[] (軸平行の有限矩形) に対する一般形**: 室 6 面 +
   障害物 AABB の各 6 面。障害物の無い直方体では Allen–Berkley の閉形式と
   厳密に一致しなければならない (次数 3 で像 63 個、次数 2 で 25 個)。
   音源と受音点が対称面上に厳密に乗ると稜線に落ちる経路が複数の面の並びから
   到達するので、像の位置と次数で重複除去している (これを外すと二重計上、
   逆に稜線を弾くと像が 12 個欠ける)。
8. **角度依存吸音は opt-in** (`.ofdx` の angle_dependent_absorption、既定
   false = 従来どおり R = sqrt(1-alpha))。on のときは局所反応境界
   R(theta) = (zeta cos - 1)/(zeta cos + 1) で、zeta は吸音表の**ランダム
   入射** alpha から Paris の式 alpha_stat(zeta) =
   (8/zeta^2)[zeta+1-2ln(1+zeta)-1/(1+zeta)] を逆に解いて決める。
   表の値を垂直入射として使うと吸音が二重に効くので、この逆解を外さないこと。
   alpha_stat の上限は 0.951 (zeta ~ 1.55) — それを超える alpha は
   クランプして warning を出す (完全吸音と偽らない)。
9. **空気吸収は ISO 9613-1**。実装は ISO 9613-2:1996 Table 2 の
   10 ℃/70 %RH 行 (0.4/1.0/1.9/3.7/9.7/32.8/117 dB/km) と 20 ℃/70 %RH 行
   (0.3/1.1/2.8/5.0/9.0/22.9/76.6) を再現することで検証してある。
   温度・湿度は**空気吸収にのみ**使い、音速は FDTD 側と揃えて 343.0 m/s 固定。
10. **残響時間は Eyring と一致しない前提で扱う**: 吸音が一様な直方体では
   初期減衰率が一致するだけで、後半は場が混合せず緩くなる。吸音が偏った室
   (客席床だけ吸音が大きいホール等) では Eyring より 10〜40% 長く出る
   (散乱係数で縮む)。これは実装の誤差ではなく幾何音響の物理なので、
   Eyring に合わせるために減衰を細工しないこと。検証 (C) が初期減衰率で
   比較しているのはこのため。
11. **有効帯域を偽らない**: metadata の `valid_band_hz` は
   [Schroeder 周波数, min(fs/2, 8 kHz バンド上端 = 11313.7 Hz)]。
   下限より下・上限より上の値を「使える」と書かないこと。
12. **数値を捏造しない**: 室外/剛体内の音源・受音点、`.ofdx` の値域外
   (`image_order` = 9 等) は非零終了 + stderr に理由。既定値に黙って
   落とさない ((G))。

## 移植性の絶対規則 (OpenPEEC portability.md を踏襲 — MSVC で実際に踏んだもの)

1. **C99 VLA 禁止** (MSVC C2057/C2466)。`malloc` + 明示インデックスの
   フラット配列 (`p[(k*ny + j)*nx + i]` の形)。
2. **OpenMP for のインデックスは事前宣言の `int`** (MSVC は OpenMP 2.0、C3015)。
3. `<complex.h>` は使わない (音響では不要)。
4. float\*/double\* の取り違え禁止 (Windows で 0xC0000005)。場は全て `double`、
   WAV 書き出し時のみ float32 に変換する。
5. libm リンクは CMake の `MATH_LIB` 変数経由 (Windows には m.lib が無い)。
6. MSVC フラグは CMakeLists の既存ブロックに従う (`/utf-8`,
   `_USE_MATH_DEFINES`, `_CRT_SECURE_NO_WARNINGS`)。円周率は `AC_PI`。
7. 外部ライブラリを追加しない。OpenMP のみ任意で `#ifdef _OPENMP` ガード。
8. ディレクトリ走査は POSIX `opendir` / Windows `FindFirstFileA` の 2 経路
   (input_ofd.c)。パス結合は `/` で良い (Windows も受け付ける)。

## 設計の規則

- グローバル変数は使わない。状態は `ac_t` コンテキスト構造体 1 個を
  main で確保して関数に渡す。確保したメモリは `ac_free()` で必ず解放する
  (CI の sanitize ジョブが LeakSanitizer で検査する。`ac_t` に動的メンバを
  足したらここにも足す)。
- `.ofd` の書式を変えない・既存キーの意味を変えない (本家 OpenFDTD 互換)。
  未知キーは無視 (前方互換)。拡張は `.ofdx` 側に足す (読むキーは
  `acoustic.absorption[]` と `acoustic.receivers[]` — 増やすときも未知キー
  無視を保つ)。receivers は GUI が持つ受音点名の引き当てにだけ使い、
  **座標一致 (許容 = dxmin)** で `.ofd` の point 行と対応させる
  (GUI はマイク配置時に point と receivers の両方へ同じ座標の行を作る)。
- 新機能には data/sample/ の**解析解付き検証ケース**を追加し、
  `acoustic_check.sh` に判定を足す (CI 3 OS + sanitize で自動実行される)。
  期待値は**コードと独立な出所** (教科書の公式・解析解) にすること。
  `.ofd` の先頭コメントに導出・許容誤差の根拠を書く。
- 既定値は「キー省略時に従来動作と完全一致」(後方互換)。
- 幾何音響側も同じ規約 : グローバル変数禁止、状態は `ga_t` 1 個、確保した
  メモリは `ga_free()` で必ず解放する (`ga_t` に動的メンバを足したらここにも
  足す)。エコーグラム・バンド別信号・FFT 作業配列はサンプル数に比例して
  大きくなるので、LeakSanitizer の検査が実用上も効く。
- 幾何音響の拡張は `.ofdx` の `acoustic.ga` に足す (未知キー無視を保つ)。
  FDTD 側は `acoustic.ga` を未知キーとして読み飛ばすので、1 つの `.ofdx` を
  両ソルバーで共有できる。**この性質を壊さないこと**。
- 幾何音響で未対応のもの (非直方体の室・回折・音源の指向性) は ReadMe の
  「幾何音響の制約」に書いてある。音源ごとのゲイン・遅延は
  `.ofdx` の `acoustic.sources[]` で対応済み (ADR-0010 Decision 7 —
  両ソルバー対称)。バンド別の散乱係数は `scattering` の配列指定で対応済み
  (基準確率での抽選 + バンド別の重み付け — 不変条件 6)。障害物の材質は
  `.ofdx` の `acoustic.ga.obstacles[]` で対応済み (既定は剛体 — FDTD 側は
  常に剛体なので、指定すると高域のみ材質が効く点は ReadMe に明記)。
  実装したら制約表からも消すこと。逆に、原理的に無理なもの (回折、
  Schroeder 周波数より下) を「対応した」と書かないこと。
- 未実装を実装済みと偽らない: 対応できない入力 (AABB 近似・複数 feed・
  非一様メッシュ) は solver.log に warning を明示する。

## CI

`.github/workflows/ci.yml`: Linux / macOS / Windows (MSVC + Ninja) の 3 OS +
`sanitize` (Linux, ASan + UBSan + LeakSanitizer)。検証スクリプトは全ジョブとも
同一の `data/sample/acoustic_check.sh` を `shell: bash` (Windows は Git Bash)
で実行し、**2 本のソルバーを 1 回で判定する**。Linux はホール規模の E2E
スモークを両ソルバーで実行。タグ `v*` push で Release に**両方の**バイナリを
添付する (`ofdx_acoustic-<os>-<arch>` の 1 アーカイブにまとめる)。
