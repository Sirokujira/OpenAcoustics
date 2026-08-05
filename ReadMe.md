# OpenAcoustics — ofdx_acoustic_fdtd

[OpenFDTD-X](../OpenFDTD-X) (Qt6 GUI) の室内音響用外部 FDTD ソルバー。
GUI から QProcess で起動され、部屋のインパルス応答 (RIR) を計算して
ファイル契約 (ADR-0007) で返す処理カーネル。

- 3 次元線形音響のスタガード格子 FDTD (leapfrog):
  p (セル中心) / vx, vy, vz (面)。
  p<sup>n+1</sup> = p<sup>n</sup> − ρc²Δt ∇·v、
  v<sup>n+1/2</sup> = v<sup>n−1/2</sup> − (Δt/ρ)∇p。
  c = 343.0 m/s、ρ = 1.204 kg/m³ (20 ℃ 空気)。
- 入力は OpenFDTD 互換 `.ofd` + `.ofdx` (JSON サイドカー、吸音率)。
- 外部ライブラリに依存しない (C11 + CMake、OpenMP のみ任意)。
  WAV / JSON は自前実装。[OpenPEEC](../OpenPEEC) と同じビルド規約・
  移植性規則を共有する姉妹プロジェクト。
- 決定性: 乱数不使用。OpenMP 並列はリダクションを持たず、
  スレッド数によらず結果は**ビット単位で一致**する。

## ビルド

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 検証 (解析解との比較 — 40 判定)
sh data/sample/acoustic_check.sh "$PWD/bin/ofdx_acoustic_fdtd" /tmp/ac-check
```

必要環境: CMake 3.18+、C11 コンパイラ (gcc / clang / MSVC)。
OpenMP は見つかれば使う (無ければシリアル — macOS AppleClang 等)。
Windows は MSVC + Ninja でビルドできる (CI と同じ)。

## 使い方

```
ofdx_acoustic_fdtd <working_dir> [<input_file.ofd>]
```

- `input_file` を省略すると `working_dir` 直下の**唯一の** `.ofd` を探す
  (無い / 複数あるときは非零終了 + stderr に理由)。
- 同じ basename の `.ofdx` があれば吸音率を読む (無ければ全壁 α = 0.1)。
- `mpiexec` 経由の起動も許容 (MPI 非対応 — 単プロセスとして動く)。
- 進捗は stdout に `progress a/b` 行 (正規表現 `^progress\s+(\d+)\s*/\s*(\d+)$`、
  50 分割)。
- 異常時 (入力なし・格子過大など) は**非零終了コード + stderr へ明確な理由**。
  合成 RIR を捏造して正常終了することはない。

### 出力 (working_dir 直下 — ADR-0007 契約)

| ファイル | 必須 | 内容 |
|---|---|---|
| `rir.wav` | ○ | 受音点 #1 の RIR (float32 モノラル WAV、fs = サンプリング周波数) |
| `rir_<名前>.wav` | — | 受音点ごとの RIR。名前は **`.ofd` point 行末の `# 名前` → `.ofdx` の `acoustic.receivers[]` (座標一致で引き当て)** の順に決まり、どちらも無ければ `rir_2.wav` 等の連番。受音点 #1 は契約どおり `rir.wav` を出すが、名前があるときは別名 `rir_<名前>.wav` も併せて出す (GUI の自動割当は名前照合のため) |
| `metadata.json` | ○ | `contract_version: 1`、ソルバー名/バージョン、格子 (Δx・セル数)、fs、音源/受音点座標、音速、実行条件 |
| `metrics.json` | — | T_Sabine 等の参考値 (GUI は自前計算を正とし突合表示のみ) |
| `solver.log` | ○ | 実行ログ (`normal end` で正常終了) |

### OpenFDTD-X (GUI) からの接続

GUI 側の `AcousticRunner` は次の順でソルバーバイナリを解決する:

1. 環境変数 `OFDX_ACOUSTIC_SOLVER` (絶対パス直接指定、最優先)
2. `OPENFDTD_ACOUSTICS_HOME` 配下
3. アプリ実行ディレクトリの `kernel/`
4. `PATH`

```bash
export OFDX_ACOUSTIC_SOLVER=/path/to/OpenAcoustics/bin/ofdx_acoustic_fdtd
```

または OpenFDTD-X のカーネルパス設定ダイアログで直接指定する。
受音点が複数あるときの `rir_<受音点名>.wav` は GUI の自動割当パターンと
噛み合う。

## 対応キーワード (.ofd)

書式は本家 OpenFDTD と互換。未知キーは無視する (前方互換)。

| キーワード | 対応 | 備考 |
|---|---|---|
| `title` | ○ | metadata / ログに記録 |
| `xmesh` / `ymesh` / `zmesh` | ○ | 計算領域と刻み。**単一の一様格子 dx = 最小刻み**に丸める (非一様は warning) |
| `geometry` | ○ | 剛体 (法線速度 0) としてボクセル化。shape 1 (直方体) は厳密、**その他の shape は AABB 近似** (warning) |
| `feed` | ○ | 音源位置。**#1 のみ使用** (複数あれば warning)。励振波形はガウシアン微分のソフト音源 (下記) |
| `point` | ○ | 受音点 (全て)。行末 `# 名前` が WAV ファイル名になる |
| `frequency1` / `frequency2` | 記録のみ | 音響ソルバーでは物理に使用しない |
| `material` / `abc` / `pml` / `rfeed` ほか | 無視 | 電磁界用 (前方互換で読み飛ばす) |

### .ofdx (JSON サイドカー)

`acoustic.absorption[]` の `enabled` 行から壁の吸音率 α を読む。
6 バンド `alpha` 配列は**帯域平均**を使う (v1)。role の対応:

| role | 壁 |
|---|---|
| 4 (Floor) | z− |
| 1 (Ceiling) | z+ |
| 2 (SideWall) | y± |
| 3 (RearWall) | x± |
| その他 / 欠落 | 既定 α = 0.1 |

外壁は局所反応インピーダンス Z = ρc(1+√(1−α))/(1−√(1−α)) の半陰的更新。
α → 1 で Z → ρc (無反射)、α → 0 で剛壁に厳密に一致する。

## 数値仕様 (v1)

- dx = `.ofd` メッシュの最小刻み。Δt = 1/fs、
  **fs = ceil(c·√3/(0.99·dx))** (CFL を満たす最小の整数 Hz — WAV の
  サンプリング周波数がそのまま整数になる)。
- 音源: ガウシアン微分パルス (DC 成分 0)。帯域は格子分解能
  **fmax = c/(10·dx)** に合わせて設計 (σ = 2/(π·fmax)、遅延 t0 = 5σ)。
- 計算時間 T = clamp(1.5·T_Sabine, 0.5 s, 3.0 s)。
  T_Sabine = 0.161·V/A、A は吸音表 × 壁面積から (表が無ければ α=0.1 の 6 面)。
- セル総数 > 3000 万は非零終了 (メッシュを粗くする案内を stderr に出す)。

## 検証ケース (data/sample/ — 期待値はすべて実装から独立な解析解)

`acoustic_check.sh` が CI (3 OS + sanitize) で全判定を実行する。
各 `.ofd` の先頭コメントに期待値の導出が書いてある。

| ケース | 検証内容 | 解析解 | 許容 |
|---|---|---|---|
| `box_modes` | 剛体閉箱 4×3×2.5 m の軸モード周波数 | f<sub>lmn</sub> = (c/2)√((l/Lx)²+…) | ±3% |
| `box_modes` | 剛壁 + DC なし音源のエネルギー保存 | 離散系は無損失 | 窓間減衰 ≤ 0.5 dB |
| `tube_a03` / `tube_a09` | 準 1D 管の端面反射係数 (α = 0.3 / 0.9) | R = √(1−α) | ±3% |
| `anechoic` | 直接音到達時刻 (全壁 α=1、r = 3 m) | t = t0 + r/c | σ + 2/fs |
| `box_pillar` | 剛体障害物 + 既定吸音の T_Sabine | 0.161·V/A = 0.409 s | ±1% |
| `box_pillar` | 決定性 (再実行 / スレッド数 1 vs 4) | — | **ビット一致** |
| `box_modes` ほか | 契約 (WAV ヘッダ・metadata・progress 書式・異常系の非零終了) | ADR-0007 | 完全一致 |
| `hall` | E2E: ホール規模 20×15×10 m (dx = 0.25 m、192,000 セル) | 直接音 t0 + r/c、Sabine 減衰 | スモーク |

## v1 の制約 (正直な現状)

- **周波数上限 fmax = c/(10·dx)**: dx = 0.1 m で 343 Hz、dx = 0.25 m で
  137 Hz。可聴帯域全体の RIR ではない (低域のモード・初期反射の検討用)。
- **吸音率は帯域平均**: 6 バンド α 配列の平均 1 値を全帯域に使う。
  周波数依存インピーダンス境界は将来課題。
- **単一音源**: `feed` が複数あっても #1 のみ (warning を出す)。
- **shape 1 (直方体) 以外の geometry は AABB 近似** (warning を出す)。
- **壁は局所反応**: 斜入射の角度依存吸音は Z 一定の範囲でのみ表現される。
- **空気吸収なし** (`air_a` は読み飛ばす)。長距離・高域では過大評価になる。
- **MPI 非対応**: `mpiexec` 起動は許容するが単プロセスで動く。

## ライセンス / 関連リポジトリ

- OpenFDTD-X (GUI): `.ofd`/`.ofdx` の書き手。契約は
  `docs/adr/0007-acoustic-solver-contract.md` と
  `tests/acoustics/mock_acoustic_solver.c` (参照実装) が正。
- OpenPEEC: ビルド規約・移植性規則の出典 (姉妹プロジェクト)。
