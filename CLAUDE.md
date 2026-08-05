# OpenAcoustics — ofdx_acoustic_fdtd

OpenFDTD-X (GUI) の室内音響用外部 FDTD ソルバー (C11)。
3 次元線形音響のスタガード格子 FDTD で部屋のインパルス応答 (RIR) を計算し、
ADR-0007 のファイル契約 (rir.wav / metadata.json / metrics.json / solver.log)
で GUI に返す。OpenPEEC の姉妹プロジェクトで、ビルド規約・移植性規則を
共有する。外部ライブラリに依存しない (C11 + CMake、OpenMP のみ任意。
WAV / JSON は自前実装 — 追加しないこと)。

**`AGENTS.md`** に同じ規約を単独で読める形でまとめてある (Codex 等、
`CLAUDE.md` を読まないエージェント向け)。**規約を変えたら両方直すこと。**

## ビルド / テスト (変更後は必ず実行 — 全判定 PASS が完了条件)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 検証 : 解析解との比較 40 判定 (モード周波数 / 反射係数 / 到達時刻 /
# エネルギー保存 / T_Sabine / 決定性 / 契約 / 異常系)
sh data/sample/acoustic_check.sh "$PWD/bin/ofdx_acoustic_fdtd" /tmp/ac-check

# E2E (ホール規模 20x15x10 m、約 2 秒)
mkdir -p /tmp/ac-hall && cp data/sample/hall.ofd data/sample/hall.ofdx /tmp/ac-hall/
./bin/ofdx_acoustic_fdtd /tmp/ac-hall && grep "normal end" /tmp/ac-hall/solver.log
```

「ビルドが通る」は合格条件ではない。**40 判定すべて OK でなければ完了ではない**
(この検証群が物理と契約の番人になっている)。バイナリは `bin/` に出る
(CMakeLists が固定)。sanitize ビルド (CI と同じ ASan+UBSan) も `bin/` を
上書きするので、検証後は Release を作り直すこと。

## ソース構成

| ファイル | 役割 |
|---|---|
| `include/acoustic.h` | コンテキスト構造体 `ac_t` と全定数 (c, ρ, CFL, 上限) |
| `src/main.c` | エントリ。引数 → 入力 → setup → run → 出力。ログ (`ac_log`/`ac_err`) |
| `src/input_ofd.c` | `.ofd` パーサ (OpenFDTD 互換、未知キーは無視) と入力探索 |
| `src/input_ofdx.c` | `.ofdx` (JSON サイドカー) の tolerant scanner — 必要キーのみ読む |
| `src/fdtd.c` | setup (格子/音源設計/境界係数/ボクセル化) と leapfrog 本体 |
| `src/wav.c` | float32 モノラル WAV 書き出し (44 byte ヘッダ) |
| `src/jsonout.c` | metadata.json / metrics.json |

## 物理不変条件 (壊すと結果が静かに狂う — 番人は acoustic_check.sh)

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
- 未実装を実装済みと偽らない: 対応できない入力 (AABB 近似・複数 feed・
  非一様メッシュ) は solver.log に warning を明示する。

## CI

`.github/workflows/ci.yml`: Linux / macOS / Windows (MSVC + Ninja) の 3 OS +
`sanitize` (Linux, ASan + UBSan + LeakSanitizer)。検証スクリプトは全ジョブとも
同一の `data/sample/acoustic_check.sh` を `shell: bash` (Windows は Git Bash)
で実行する。Linux はホール規模の E2E スモークも実行。
タグ `v*` push で Release にバイナリ添付。
