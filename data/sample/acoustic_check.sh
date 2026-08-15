#!/bin/sh
# acoustic_check.sh — ソルバー 2 本の検証 (CI 用、OpenPEEC peec_check.sh の流儀)
#
# data/sample/ の各ケースを実行し、rir.wav を解析解と比較する。
# 期待値の導出は各 .ofd ファイルのコメント参照 (すべて実装から独立な解析解)。
#
# ── ofdx_acoustic_fdtd (低域担当、FDTD) の判定 :
#  (a) 剛体閉箱 4x3x2.5 m の軸モード f_100/f_010/f_001 が ±3% (走査 DFT)
#  (b) 準 1D 管の端面反射 R = sqrt(1-alpha) が ±3% (alpha=0.3, 0.9)
#  (c) 直接音到達 : |p| ピークが t = t0 + r/c ± (sigma + 2/fs)
#  (d) 剛壁のみ (alpha=0) でエネルギー保存 (1.4 s 窓 2 つの減衰 < 0.5 dB)
#  (e) 決定性 : 同一入力 2 回 / OMP_NUM_THREADS=1 と 4 で rir.wav ビット一致
#  (f) 契約 : 出力 4 ファイル、WAV ヘッダ、progress 行書式、異常系の非零終了
#  (h) 複数音源 (multi_source) : 2 feed の同時発火が各 feed 単独の和に一致
#      (離散更新の線形性、L2 相対 1e-5)
#
# ── ofdx_acoustic_ga (高域担当、幾何音響) の判定 :
#  (A) 自由音場 : 直接音の時刻 t = r/c と振幅 1/(4 pi r) が ±1%
#  (B) 剛体床 1 枚 : 直接音 + 1 次反射の 2 発が鏡像法の閉形式と ±1%
#  (C) 直方体箱 : バンド別 T60 (500 Hz / 4 kHz) が Eyring 式と ±5%
#  (D) 空気吸収 : 100 m・4 kHz の超過減衰が ISO 9613-2 Table 2 と ±0.5 dB
#  (E) 全 alpha=0 + 空気吸収 off : エコーグラムが減衰しない (< 0.5 dB)
#  (F) 決定性 : 再実行 / OMP_NUM_THREADS=1 と 4 で rir.wav ビット一致
#  (G) 契約 : 出力ファイル、WAV ヘッダ (48 kHz)、progress 行書式、
#      metadata の valid_band_hz、遮蔽・近似の warning、異常系の非零終了
#  (H) 剛体反射板の 1 次鏡面反射が閉形式と ±1% + 有限面からはみ出す像の棄却
#  (B2) 面ごとの散乱係数が室全体の既定を上書きする (床 s = 0 / 1)
#  (I) 局所反応境界の角度依存反射 R(theta) が Paris の式の逆解と ±1%
#  (J) 障害物の材質 (acoustic.ga.obstacles) : alpha で鏡面反射が
#      sqrt(1-alpha) 倍、scattering = 1 で鏡面像が消える、不正 index は非零終了
#  (K) 複数音源 (multi_source) : 2 音源の直接音がそれぞれ 1/(4 pi r_i) ±1%、
#      既定は feed #1 のみ + warning (後方互換)、室外の音源は非零終了
#
# WAV の読みは od -t f4 (float32 リトルエンディアン)。CI の 3 OS
# (Linux / macOS / Windows Git Bash) はいずれもリトルエンディアンかつ
# coreutils/BSD od が -j/-t f4 に対応している。
#
# 使い方 :
#   acoustic_check.sh <ofdx_acoustic_fdtd(絶対パス)> [作業ディレクトリ] [<ofdx_acoustic_ga>]
# 第 3 引数を省略すると第 1 引数のファイル名の "fdtd" を "ga" に置き換えた
# パスを使う (.exe も含めて解決できる)。

set -e

SOLVER="$1"
WORK="${2:-.}"
SRC="$(cd "$(dirname "$0")" && pwd)"
PI=3.14159265358979

if [ -z "$SOLVER" ]; then
	echo "Usage: acoustic_check.sh <ofdx_acoustic_fdtd> [workdir] [ofdx_acoustic_ga]" >&2
	exit 2
fi

GA_SOLVER="${3:-}"
if [ -z "$GA_SOLVER" ]; then
	GA_SOLVER=$(printf '%s' "$SOLVER" | sed 's/ofdx_acoustic_fdtd/ofdx_acoustic_ga/')
fi

mkdir -p "$WORK"
status=0

# run_case <case名> <ファイル...> : 作業サブディレクトリに copy して実行。
# 入力ファイル名は渡さない (working_dir 直下の唯一の .ofd 自動探索も検証対象)
run_case() {
	_name="$1"; shift
	_dir="$WORK/$_name"
	rm -rf "$_dir"
	mkdir -p "$_dir"
	for _f in "$@"; do cp "$SRC/$_f" "$_dir/"; done
	if "$SOLVER" "$_dir" > "$_dir/stdout.log" 2> "$_dir/stderr.log"; then :; else
		echo "*** $_name: solver failed (see $_dir/stderr.log)" >&2
		status=1
		return 1
	fi
}

# chk <label> <actual> <expected> <tol(相対)>
chk() {
	awk -v a="$2" -v e="$3" -v tol="$4" -v lb="$1" 'BEGIN {
		d = (a - e) / e; ad = (d < 0) ? -d : d;
		printf "%-26s actual=%.6g expected=%.6g -> %s %+.4f%%\n", lb, a, e, (ad <= tol) ? "OK" : "NG", d * 100;
		exit (ad <= tol) ? 0 : 1
	}' || status=1
}

# chkabs <label> <actual> <expected> <tol(絶対)> <単位>
chkabs() {
	awk -v a="$2" -v e="$3" -v tol="$4" -v u="$5" -v lb="$1" 'BEGIN {
		d = a - e; ad = (d < 0) ? -d : d;
		printf "%-26s actual=%.6g expected=%.6g -> %s %+.4f %s\n", lb, a, e, (ad <= tol) ? "OK" : "NG", d, u;
		exit (ad <= tol) ? 0 : 1
	}' || status=1
}

# ok <label> <条件の真偽 (0=真)> : grep 等の結果をそのまま判定にする
say_ok() { printf "%-26s -> OK%s\n" "$1" "${2:-}"; }
say_ng() { printf "%-26s -> NG%s\n" "$1" "${2:-}" >&2; status=1; }

# WAV の float32 サンプル列を 1 行 1 値で dump (44 byte 標準ヘッダ想定)
dump() {
	od -An -v -t f4 -j 44 "$1" | awk '{for (i = 1; i <= NF; i++) print $i}'
}

# WAV ヘッダのリトルエンディアン整数 (u16/u32) を取り出す
u16() { od -An -v -t u2 -j "$2" -N 2 "$1" | tr -d ' \t'; }
u32() { od -An -v -t u4 -j "$2" -N 4 "$1" | tr -d ' \t'; }
tag4() { dd if="$1" bs=1 skip="$2" count=4 2>/dev/null; }

echo "--- (f) contract: rigid box (explicit input file name)"
# この 1 ケースだけ入力ファイル名を明示して起動する経路も検証する
d="$WORK/box_modes"
rm -rf "$d"; mkdir -p "$d"
cp "$SRC/box_modes.ofd" "$SRC/box_modes.ofdx" "$d/"
if "$SOLVER" "$d" box_modes.ofd > "$d/stdout.log" 2> "$d/stderr.log"; then
	printf "%-26s -> OK\n" "box_modes run"
else
	printf "%-26s -> NG (solver failed)\n" "box_modes run" >&2
	status=1
fi
for f in rir.wav rir_mid.wav metadata.json metrics.json solver.log; do
	if [ -s "$d/$f" ]; then
		printf "%-26s -> OK\n" "exists: $f"
	else
		printf "%-26s -> NG (missing)\n" "exists: $f" >&2
		status=1
	fi
done
grep -q "normal end" "$d/solver.log" || { echo "*** no 'normal end' in solver.log" >&2; status=1; }
# fs は仕様式 ceil(c*sqrt(3)/(0.99*dx)) から独立に計算する (dx=0.1 -> 6001)
fs_box=$(awk 'BEGIN{x=343*sqrt(3)/(0.99*0.1); f=int(x); if (f<x) f++; print f}')
chk "metadata sample_rate" "$(grep '"sample_rate"' "$d/metadata.json" | tr -dc '0-9')" "$fs_box" 0
chk "metadata contract_ver" "$(grep '"contract_version"' "$d/metadata.json" | tr -dc '0-9')" 1 0
# WAV ヘッダ : RIFF/WAVE、fmt tag 3 (IEEE float)、mono、32 bit、fs、サイズ整合
[ "$(tag4 "$d/rir.wav" 0)" = "RIFF" ] && [ "$(tag4 "$d/rir.wav" 8)" = "WAVE" ] \
	&& printf "%-26s -> OK\n" "wav RIFF/WAVE tags" \
	|| { printf "%-26s -> NG\n" "wav RIFF/WAVE tags" >&2; status=1; }
chk "wav format tag (float)" "$(u16 "$d/rir.wav" 20)" 3 0
chk "wav channels" "$(u16 "$d/rir.wav" 22)" 1 0
chk "wav bits" "$(u16 "$d/rir.wav" 34)" 32 0
chk "wav sample rate" "$(u32 "$d/rir.wav" 24)" "$fs_box" 0
# T = clamp(1.5*T_Sabine, 0.5, 3.0)。剛壁 (A=0) なので T=3.0 s -> steps = 3*fs
size=$(wc -c < "$d/rir.wav" | tr -d ' \t')
chk "wav file size" "$size" "$(awk -v f="$fs_box" 'BEGIN{print 44 + 4*3*f}')" 0
chk "wav data size field" "$(u32 "$d/rir.wav" 40)" "$(awk -v s="$size" 'BEGIN{print s-44}')" 0
chk "wav riff size field" "$(u32 "$d/rir.wav" 4)" "$(awk -v s="$size" 'BEGIN{print s-8}')" 0
# progress 行 : "progress" を含む行はすべて ^progress <a>/<b>$ 書式で 10 行以上
nprog=$(grep -c "progress" "$d/stdout.log" || true)
ngood=$(grep -Ec "^progress [0-9]+/[0-9]+$" "$d/stdout.log" || true)
if [ "$nprog" -ge 10 ] && [ "$nprog" = "$ngood" ]; then
	printf "%-26s -> OK (%d lines)\n" "progress format" "$nprog"
else
	printf "%-26s -> NG (progress=%s well-formed=%s)\n" "progress format" "$nprog" "$ngood" >&2
	status=1
fi

echo "--- (a) axial modes of the rigid box (analytic: f = c*l/(2L), +-3%)"
# f_100 = 343/8 = 42.875, f_010 = 343/6 = 57.1667, f_001 = 343/5 = 68.6 Hz。
# 走査 DFT (刻み 0.1 Hz、窓 ±3.5% — 隣接モード 71.46 Hz を含まない)。
# 導出・格子分散の評価は box_modes.ofd のコメント参照。
dump "$d/rir.wav" > "$d/rir.txt"
awk -v fs="$fs_box" -v pi="$PI" '
	{ s[NR] = $1 }
	END {
		n = NR; nt = 3;
		T[1] = 42.875; T[2] = 57.1666667; T[3] = 68.6;
		bad = 0;
		for (m = 1; m <= nt; m++) {
			fe = T[m]; best = -1; bf = 0;
			for (f = fe * 0.965; f <= fe * 1.035; f += 0.1) {
				cr = 0; ci = 0; w = 2 * pi * f / fs;
				for (k = 1; k <= n; k++) { cr += s[k] * cos(w * k); ci += s[k] * sin(w * k) }
				a = cr * cr + ci * ci;
				if (a > best) { best = a; bf = f }
			}
			dv = (bf - fe) / fe; ad = (dv < 0) ? -dv : dv;
			printf "%-26s actual=%.6g expected=%.6g -> %s %+.4f%%\n",
				"mode " fe " Hz", bf, fe, (ad <= 0.03) ? "OK" : "NG", dv * 100;
			if (ad > 0.03) bad = 1;
		}
		exit bad
	}' "$d/rir.txt" || status=1

echo "--- (d) energy conservation with rigid walls (< 0.5 dB)"
# 全壁剛体 + DC の無い音源 -> 離散系は無損失。受音点 2 点の平均自乗圧を
# 1.4 s 窓 2 つで比較する (根拠は box_modes.ofd のコメント)。
dump "$d/rir_mid.wav" > "$d/rir_mid.txt"
awk -v fs="$fs_box" '
	{ t = FNR / fs; s = $1;
	  if (t > 0.2 && t <= 1.6) e1 += s * s;
	  else if (t > 1.6 && t <= 3.0) e2 += s * s }
	END {
		db = 10 * log(e1 / e2) / log(10); ad = (db < 0) ? -db : db;
		printf "%-26s decay=%.4f dB -> %s (<= 0.5 dB)\n", "rigid-box energy", db, (ad <= 0.5) ? "OK" : "NG";
		exit (ad <= 0.5) ? 0 : 1
	}' "$d/rir.txt" "$d/rir_mid.txt" || status=1

echo "--- (b) tube end reflection R = sqrt(1-alpha) (+-3%)"
# 測定は窓内エネルギー比 sqrt(E_refl/E_inc) — 無損失 1D 格子では窓内
# エネルギーは伝搬不変 (Parseval) で、ピーク比の格子分散バイアスを受けない。
# 窓 [2,14] / [15.5,28] / [30,42] ms の導出は tube_a03.ofd のコメント参照。
tube_R() {   # tube_R <case dir> : E 比 2 つ (左端・右端) を出力
	dump "$1/rir.wav" | awk -v fs="$2" '
		{ t = NR / fs; s = $1;
		  if (t >= 0.002 && t <= 0.014) e1 += s * s;
		  else if (t >= 0.0155 && t <= 0.028) e2 += s * s;
		  else if (t >= 0.030 && t <= 0.042) e3 += s * s }
		END { printf "%.9e %.9e", sqrt(e2 / e1), sqrt(e3 / e1) }'
}
fs_tube=$(awk 'BEGIN{x=343*sqrt(3)/(0.99*0.05); f=int(x); if (f<x) f++; print f}')
run_case tube_a03 tube_a03.ofd tube_a03.ofdx && {
	res=$(tube_R "$WORK/tube_a03" "$fs_tube")
	chk "tube a=0.3 R (left)"  "${res% *}" 0.836660027 0.03   # sqrt(0.7)
	chk "tube a=0.3 R (right)" "${res#* }" 0.836660027 0.03
}
run_case tube_a09 tube_a09.ofd tube_a09.ofdx && {
	res=$(tube_R "$WORK/tube_a09" "$fs_tube")
	chk "tube a=0.9 R (left)"  "${res% *}" 0.316227766 0.03   # sqrt(0.1)
	chk "tube a=0.9 R (right)" "${res#* }" 0.316227766 0.03
}

echo "--- (c) direct sound arrival t = t0 + r/c"
# 全壁 alpha=1 の部屋で r=3.0 m。受音波形は s(t) の時間微分 (Ricker 型、
# 対称) なので |p| ピークはパルス中心に一致する。許容 sigma + 2/fs の
# 根拠は anechoic.ofd のコメント参照。
run_case anechoic anechoic.ofd anechoic.ofdx && {
	dump "$WORK/anechoic/rir.wav" | awk -v fs="$fs_box" -v pi="$PI" '
		BEGIN { sg = 20 * 0.1 / (pi * 343); t0 = 5 * sg; te = t0 + 3.0 / 343;
		        w = 3 * sg; tol = sg + 2 / fs }
		{ t = NR / fs; a = ($1 < 0) ? -$1 : $1;
		  if (t >= te - w && t <= te + w && a > pk) { pk = a; tp = t } }
		END {
			d = tp - te; ad = (d < 0) ? -d : d;
			printf "%-26s actual=%.6g expected=%.6g -> %s (tol %.4g s)\n",
				"direct arrival [s]", tp, te, (ad <= tol) ? "OK" : "NG", tol;
			exit (ad <= tol) ? 0 : 1
		}' || status=1
}

echo "--- smoke: rigid obstacles + default absorption (no .ofdx)"
run_case box_pillar box_pillar.ofd && {
	d="$WORK/box_pillar"
	grep -q "normal end" "$d/solver.log" \
		&& printf "%-26s -> OK\n" "box_pillar normal end" \
		|| { printf "%-26s -> NG\n" "box_pillar normal end" >&2; status=1; }
	# 複数 feed の warning と AABB 近似の warning が明示されること (v1 仕様)
	grep -q "warning: 2 feeds" "$d/solver.log" \
		&& printf "%-26s -> OK\n" "multi-feed warning" \
		|| { printf "%-26s -> NG\n" "multi-feed warning" >&2; status=1; }
	grep -q "AABB" "$d/solver.log" \
		&& printf "%-26s -> OK\n" "AABB approx warning" \
		|| { printf "%-26s -> NG\n" "AABB approx warning" >&2; status=1; }
	# 吸音表なし -> 既定 alpha=0.1 の 6 面から T_Sabine = 0.409 s (導出は .ofd)
	chk "default-alpha T_Sabine" "$(awk -F: '/t_sabine_s/ {gsub(/[ ,\r]/, "", $2); print $2}' "$d/metrics.json")" 0.4093 0.01
}

echo "--- (e) determinism (bit-identical reruns and thread invariance)"
# 乱数不使用。OpenMP 並列 (p 更新・v 更新のセルループ) はリダクションを
# 持たないので、スレッド数によらずビット単位で一致する (OpenPEEC と同じ原則)。
d="$WORK/box_pillar"
cp "$d/rir.wav" "$d/rir_run1.wav"
if "$SOLVER" "$d" > /dev/null 2>&1 && cmp -s "$d/rir_run1.wav" "$d/rir.wav"; then
	printf "%-26s -> OK (2 回実行がビット一致)\n" "rerun determinism"
else
	printf "%-26s -> NG\n" "rerun determinism" >&2
	status=1
fi
OMP_NUM_THREADS=1 "$SOLVER" "$d" > /dev/null 2>&1
cp "$d/rir.wav" "$d/rir_n1.wav"
OMP_NUM_THREADS=4 "$SOLVER" "$d" > /dev/null 2>&1
if cmp -s "$d/rir_n1.wav" "$d/rir.wav"; then
	printf "%-26s -> OK (1 と 4 スレッドがビット一致)\n" "thread invariance"
else
	printf "%-26s -> NG\n" "thread invariance" >&2
	status=1
fi

echo "--- (g) receiver names from the .ofdx sidecar (GUI 連携)"
# GUI (OpenFDTD-X) は受音点名を .ofdx acoustic.receivers[] に持ち、.ofd の
# point 行には名前を書かない。座標一致で名前を引き当て、rir_<受音点名>.wav
# を出すことで GUI の「フォルダから自動割当」(名前照合) と噛み合う。
# 期待するファイル名の根拠は named_recv.ofd の先頭コメント。
run_case named_recv named_recv.ofd named_recv.ofdx && {
	d="$WORK/named_recv"
	for f in "rir.wav" "rir_マイク_1.wav" "rir_mic-B.wav" "rir_3.wav"; do
		if [ -f "$d/$f" ]; then
			printf "%-26s -> OK\n" "named recv: $f"
		else
			printf "%-26s -> NG (missing)\n" "named recv: $f" >&2; status=1
		fi
	done
	# 無効行 (enabled=false) と座標の合わない行は名前を配らない
	if [ -f "$d/rir_disabled.wav" ] || [ -f "$d/rir_elsewhere.wav" ]; then
		printf "%-26s -> NG (disabled/unmatched row was used)\n" "named recv: filtering" >&2
		status=1
	else
		printf "%-26s -> OK (disabled と不一致行は無視)\n" "named recv: filtering"
	fi
}

echo "--- (f) error paths (honest non-zero exit, no fabricated output)"
# 入力なし
rm -rf "$WORK/empty"; mkdir -p "$WORK/empty"
if "$SOLVER" "$WORK/empty" > /dev/null 2> "$WORK/empty/stderr.log"; then
	printf "%-26s -> NG (should fail)\n" "no-input exit" >&2; status=1
elif [ -s "$WORK/empty/stderr.log" ] && [ ! -f "$WORK/empty/rir.wav" ]; then
	printf "%-26s -> OK (non-zero, stderr reason, no rir.wav)\n" "no-input exit"
else
	printf "%-26s -> NG (missing stderr reason or fabricated rir)\n" "no-input exit" >&2; status=1
fi
# セル総数 > 3000 万
rm -rf "$WORK/toolarge"; mkdir -p "$WORK/toolarge"
cp "$SRC/toolarge.ofd" "$WORK/toolarge/"
if "$SOLVER" "$WORK/toolarge" > /dev/null 2> "$WORK/toolarge/stderr.log"; then
	printf "%-26s -> NG (should fail)\n" "30M-cell limit" >&2; status=1
elif grep -qi "coarsen" "$WORK/toolarge/stderr.log" && [ ! -f "$WORK/toolarge/rir.wav" ]; then
	printf "%-26s -> OK (refused with mesh advice)\n" "30M-cell limit"
else
	printf "%-26s -> NG (no advice in stderr)\n" "30M-cell limit" >&2; status=1
fi
# 引数なし (usage)
if "$SOLVER" > /dev/null 2>&1; then
	printf "%-26s -> NG (should fail)\n" "usage exit" >&2; status=1
else
	printf "%-26s -> OK\n" "usage exit"
fi


echo "--- (h) multi-source superposition (linearity, .ofdx acoustic.multi_source)"
# super2.ofd の 2 feed を multi_source で同時発火した RIR は、各 feed 単独の
# RIR の標本ごとの和に一致する (離散更新の線形性 — 導出と許容誤差の根拠は
# super2.ofd の先頭コメント。float32 量子化 ~6e-8 が支配、許容 1e-5)。
for v in ab a b; do rm -rf "$WORK/super_$v"; mkdir -p "$WORK/super_$v"; done
cp "$SRC/super2.ofd" "$SRC/super2.ofdx" "$WORK/super_ab/"
grep -v "^feed = z 3.15" "$SRC/super2.ofd" > "$WORK/super_a/super2.ofd"
cp "$SRC/super2.ofdx" "$WORK/super_a/"
grep -v "^feed = z 0.55" "$SRC/super2.ofd" > "$WORK/super_b/super2.ofd"
cp "$SRC/super2.ofdx" "$WORK/super_b/"
sup_ok=1
for v in ab a b; do
	"$SOLVER" "$WORK/super_$v" > /dev/null 2>&1 || sup_ok=0
done
if [ "$sup_ok" = "1" ]; then
	for v in ab a b; do dump "$WORK/super_$v/rir.wav" > "$WORK/super_$v/rir.txt"; done
	paste "$WORK/super_ab/rir.txt" "$WORK/super_a/rir.txt" "$WORK/super_b/rir.txt" | awk \
		'{ dd = $1 - ($2 + $3); e2 += dd * dd; r2 += $1 * $1 }
		END {
			rel = sqrt(e2 / r2);
			printf "%-26s L2 rel err=%.3e -> %s (<= 1e-5)\n", "fdtd superposition", rel, (rel <= 1e-5) ? "OK" : "NG";
			exit (rel <= 1e-5) ? 0 : 1
		}' || status=1
	grep -q '"multi_source": true' "$WORK/super_ab/metadata.json" \
		&& say_ok "fdtd multi_source metadata" || say_ng "fdtd multi_source metadata"
	grep -Fq '{ "pos_m": [3.15, 2.35, 1.95] }' "$WORK/super_ab/metadata.json" \
		&& say_ok "fdtd sources listed" || say_ng "fdtd sources listed"
else
	say_ng "fdtd superposition" " (solver failed)"
fi

###############################################################################
# ofdx_acoustic_ga (幾何音響、高域担当) の検証
#
# 振幅規約 : 自由音場の直接音が 1/(4 pi r)。離散 RIR ではこれを
# 「到達近傍の標本和 (= 単位標本利得)」として測る — 3 次ラグランジュ分数
# 遅延で置いているので標本和と 1 次モーメントが厳密に一致する。
# 時間原点 : t = 0 が音源発火時刻。標本 n の時刻は n/fs (fs = 48000 固定)。
###############################################################################
echo ""
echo "=== ofdx_acoustic_ga (geometric acoustics, high band) ==="

if [ ! -x "$GA_SOLVER" ] && [ ! -f "$GA_SOLVER" ]; then
	echo "*** ga solver not found: $GA_SOLVER" >&2
	status=1
fi

ga_run() {
	_name="$1"; shift
	_dir="$WORK/$_name"
	rm -rf "$_dir"
	mkdir -p "$_dir"
	for _f in "$@"; do cp "$SRC/$_f" "$_dir/"; done
	if "$GA_SOLVER" "$_dir" > "$_dir/stdout.log" 2> "$_dir/stderr.log"; then :; else
		echo "*** $_name: ga solver failed (see $_dir/stderr.log)" >&2
		status=1
		return 1
	fi
}

# ga_sum <wavのdumpファイル> <中心時刻[s]> <半幅[s]> : 窓内の標本和 (= 振幅)
ga_sum() {
	awk -v fs=48000 -v tc="$2" -v hw="$3" '
		{ t = (NR - 1) / fs; if (t > tc - hw && t < tc + hw) s += $1 }
		END { printf "%.10e", s }' "$1"
}
# ga_peak <dump> <中心時刻> <半幅> : 窓内で |p| が最大になる時刻 [s]
ga_peak() {
	awk -v fs=48000 -v tc="$2" -v hw="$3" '
		{ t = (NR - 1) / fs; a = ($1 < 0) ? -$1 : $1;
		  if (t > tc - hw && t < tc + hw && a > pk) { pk = a; tp = t } }
		END { printf "%.10e", tp }' "$1"
}
# ga_mag <dump> <周波数[Hz]> : その周波数の DFT 振幅
#   バンド重みはオクターブ中心を節点とする単位分割ハットなので、
#   f = バンド中心では他バンドが 0 になりそのバンドの利得だけが読める。
ga_mag() {
	awk -v fs=48000 -v f0="$2" -v pi="$PI" '
		{ s[NR-1] = $1; n = NR }
		END {
			w = 2 * pi * f0 / fs;
			for (i = 0; i < n; i++) { cr += s[i] * cos(w * i); ci -= s[i] * sin(w * i) }
			printf "%.10e", sqrt(cr * cr + ci * ci)
		}' "$1"
}
# ga_t60 <dump> <中心周波数> : バンドパス -> Schroeder 逆積分 ->
#   -5 .. -15 dB の最小二乗直線の傾きから T60 [s]。
#   バンドパスは RBJ の 2 次 (定 0 dB ピーク利得) を 2 段。Q = 2.5。
#   -25 dB まで使わない理由は ga_t60.ofd の先頭コメント参照
#   (鏡面直方体の減衰は指数関数の重ね合わせで対数軸に凸 = 後半が緩い)。
ga_t60() {
	awk -v fs=48000 -v f0="$2" -v Q=2.5 -v pi="$PI" '
		{ x[NR-1] = $1; n = NR }
		END {
			w0 = 2 * pi * f0 / fs; al = sin(w0) / (2 * Q); cw = cos(w0);
			b0 = al / (1 + al); b2 = -al / (1 + al);
			a1 = (-2 * cw) / (1 + al); a2 = (1 - al) / (1 + al);
			for (p = 0; p < 2; p++) {
				z1 = 0; z2 = 0; y1 = 0; y2 = 0;
				for (i = 0; i < n; i++) {
					v = x[i];
					y = b0 * v + b2 * z2 - a1 * y1 - a2 * y2;
					z2 = z1; z1 = v; y2 = y1; y1 = y; x[i] = y;
				}
			}
			s = 0;
			for (i = n - 1; i >= 0; i--) { s += x[i] * x[i]; S[i] = s }
			if (S[0] <= 0) { printf "-1"; exit }
			for (i = 0; i < n; i++) {
				db = 10 * log(S[i] / S[0]) / log(10);
				if (db <= -5 && db >= -15) {
					t = i / fs; sx += t; sy += db; sxx += t * t; sxy += t * db; m++
				}
			}
			if (m < 100) { printf "-1"; exit }
			slope = (m * sxy - sx * sy) / (m * sxx - sx * sx);
			printf "%.10e", -60 / slope
		}' "$1"
}

echo "--- (G) ga contract: free-field case (explicit input file name)"
d="$WORK/ga_freefield"
rm -rf "$d"; mkdir -p "$d"
cp "$SRC/ga_freefield.ofd" "$SRC/ga_freefield.ofdx" "$d/"
if "$GA_SOLVER" "$d" ga_freefield.ofd > "$d/stdout.log" 2> "$d/stderr.log"; then
	say_ok "ga_freefield run"
else
	say_ng "ga_freefield run" " (solver failed)"
fi
for f in rir.wav rir_near.wav rir_far.wav metadata.json metrics.json solver.log; do
	if [ -s "$d/$f" ]; then say_ok "ga exists: $f"; else say_ng "ga exists: $f" " (missing)"; fi
done
grep -q "normal end" "$d/solver.log" || { echo "*** no 'normal end' in ga solver.log" >&2; status=1; }
chk "ga metadata sample_rate" "$(grep '"sample_rate"' "$d/metadata.json" | tr -dc '0-9')" 48000 0
chk "ga metadata contract_ver" "$(grep '"contract_version"' "$d/metadata.json" | tr -dc '0-9')" 1 0
# 有効帯域の上端は 8 kHz バンドの上端 8000*sqrt(2) = 11313.7 Hz (ナイキスト未満)
chk "ga valid_band upper" \
	"$(sed -n 's/.*"valid_band_hz": \[[^,]*, \([0-9.e+-]*\)\].*/\1/p' "$d/metadata.json")" \
	"$(awk 'BEGIN{print 8000*sqrt(2)}')" 0.001
# 時間原点の規約 : t0_s = 0 (音源発火時刻)
grep -q '"t0_s": 0' "$d/metadata.json" && say_ok "ga time origin t0_s=0" \
	|| say_ng "ga time origin t0_s=0"
[ "$(tag4 "$d/rir.wav" 0)" = "RIFF" ] && [ "$(tag4 "$d/rir.wav" 8)" = "WAVE" ] \
	&& say_ok "ga wav RIFF/WAVE tags" || say_ng "ga wav RIFF/WAVE tags"
chk "ga wav format tag (float)" "$(u16 "$d/rir.wav" 20)" 3 0
chk "ga wav channels" "$(u16 "$d/rir.wav" 22)" 1 0
chk "ga wav bits" "$(u16 "$d/rir.wav" 34)" 32 0
chk "ga wav sample rate" "$(u32 "$d/rir.wav" 24)" 48000 0
size=$(wc -c < "$d/rir.wav" | tr -d ' \t')
chk "ga wav data size field" "$(u32 "$d/rir.wav" 40)" "$(awk -v s="$size" 'BEGIN{print s-44}')" 0
chk "ga wav riff size field" "$(u32 "$d/rir.wav" 4)" "$(awk -v s="$size" 'BEGIN{print s-8}')" 0
nprog=$(grep -c "progress" "$d/stdout.log" || true)
ngood=$(grep -Ec "^progress [0-9]+/[0-9]+$" "$d/stdout.log" || true)
if [ "$nprog" -ge 10 ] && [ "$nprog" = "$ngood" ]; then
	printf "%-26s -> OK (%d lines)\n" "ga progress format" "$nprog"
else
	printf "%-26s -> NG (progress=%s well-formed=%s)\n" "ga progress format" "$nprog" "$ngood" >&2
	status=1
fi

echo "--- (A) free field: direct sound t = r/c and amplitude 1/(4 pi r) (+-1%)"
# 導出は ga_freefield.ofd の先頭コメント。r = 10 m、c = 343 m/s。
dump "$d/rir.wav" > "$d/near.txt"
dump "$d/rir_far.wav" > "$d/far.txt"
t_dir=$(awk 'BEGIN{printf "%.10f", 10.0/343.0}')
chk "ga direct amplitude" "$(ga_sum "$d/near.txt" "$t_dir" 0.001)" \
	"$(awk -v pi="$PI" 'BEGIN{printf "%.10e", 1/(4*pi*10.0)}')" 0.01
chk "ga direct arrival [s]" "$(ga_peak "$d/near.txt" "$t_dir" 0.001)" "$t_dir" 0.01

echo "--- (D) air absorption at 4 kHz over 100 m (ISO 9613-2 Table 2, +-0.5 dB)"
# 10 degC / 70 %RH / 101.325 kPa の 4 kHz は 32.8 dB/km = 0.0328 dB/m。
# 距離減衰 1/(4 pi r) を外した超過減衰 = a (r2 - r1) = 3.28 dB。
m1=$(ga_mag "$d/near.txt" 4000)
m2=$(ga_mag "$d/far.txt" 4000)
chkabs "ga air 4 kHz / 100 m" \
	"$(awk -v a="$m1" -v b="$m2" 'BEGIN{printf "%.6f", 20*log(a*10.0/(b*110.0))/log(10)}')" \
	3.28 0.5 dB

echo "--- (B) rigid floor: direct + 1st-order image, per-surface scattering (+-1%)"
# r0 = sqrt(32) = 5.65685425 m、r1 = sqrt(52) = 7.21110255 m (ga_floor.ofd 参照)。
# 床だけ面ごとの散乱係数 s = 0 を与えてある (室の既定は 0.5)。面ごとの指定が
# 効いていれば鏡面成分は減らず 1/(4 pi r) のまま — 下の (B2) と対になる判定。
ga_run ga_floor ga_floor.ofd ga_floor.ofdx && {
	d="$WORK/ga_floor"
	dump "$d/rir.wav" > "$d/rir.txt"
	for k in 0 1; do
		r=$(awk -v k=$k 'BEGIN{printf "%.10f", (k==0)?sqrt(32):sqrt(52)}')
		te=$(awk -v r="$r" 'BEGIN{printf "%.10f", r/343.0}')
		chk "ga floor #$k amplitude" "$(ga_sum "$d/rir.txt" "$te" 0.0015)" \
			"$(awk -v r="$r" -v pi="$PI" 'BEGIN{printf "%.10e", 1/(4*pi*r)}')" 0.01
		chk "ga floor #$k arrival [s]" "$(ga_peak "$d/rir.txt" "$te" 0.0015)" "$te" 0.01
	done
}

echo "--- (C) band T60 of a 10 m cube vs the Eyring formula (+-5%)"
# 期待値は ga_t60.ofd の先頭コメントの手計算 (Eyring + 空気吸収 4mV、
# 空気吸収係数は ISO 9613-2 Table 2 の 10 degC / 70 %RH の行)。
ga_run ga_t60 ga_t60.ofd ga_t60.ofdx && {
	d="$WORK/ga_t60"
	dump "$d/rir.wav" > "$d/rir.txt"
	chk "ga T60 500 Hz [s]"  "$(ga_t60 "$d/rir.txt" 500)"  2.4782 0.05
	chk "ga T60 4 kHz [s]"   "$(ga_t60 "$d/rir.txt" 4000)" 1.1942 0.05
}

echo "--- (E) lossless room (alpha = 0, air absorption off): no decay (< 0.5 dB)"
ga_run ga_lossless ga_lossless.ofd ga_lossless.ofdx && {
	d="$WORK/ga_lossless"
	dump "$d/rir.wav" > "$d/rir.txt"
	awk -v fs=48000 '
		{ t = (NR - 1) / fs; s = $1;
		  if (t >= 0.30 && t < 0.60) e1 += s * s;
		  else if (t >= 0.60 && t < 0.90) e2 += s * s }
		END {
			db = 10 * log(e1 / e2) / log(10); ad = (db < 0) ? -db : db;
			printf "%-26s decay=%.4f dB -> %s (<= 0.5 dB)\n", "ga lossless energy", db, (ad <= 0.5) ? "OK" : "NG";
			exit (ad <= 0.5) ? 0 : 1
		}' "$d/rir.txt" || status=1
}

echo "--- (G) occlusion by rigid obstacles + honest warnings"
ga_run ga_pillar ga_pillar.ofd ga_pillar.ofdx && {
	d="$WORK/ga_pillar"
	grep -q "normal end" "$d/solver.log" && say_ok "ga_pillar normal end" \
		|| say_ng "ga_pillar normal end"
	grep -q "warning: 2 feeds" "$d/solver.log" && say_ok "ga multi-feed warning" \
		|| say_ng "ga multi-feed warning"
	grep -q "AABB" "$d/solver.log" && say_ok "ga AABB approx warning" \
		|| say_ng "ga AABB approx warning"
	# 受音点 #1 は衝立で直接音が遮られ、#2 は遮られない (導出は ga_pillar.ofd)
	dump "$d/rir.wav" > "$d/blocked.txt"
	dump "$d/rir_clear.wav" > "$d/clear.txt"
	tb=$(awk 'BEGIN{printf "%.10f", 8.0/343.0}')
	tc=$(awk 'BEGIN{printf "%.10f", sqrt(73.0)/343.0}')
	awk -v a="$(ga_sum "$d/blocked.txt" "$tb" 0.0008)" -v pi="$PI" 'BEGIN {
		ref = 1/(4*pi*8.0); r = (a < 0 ? -a : a) / ref;
		printf "%-26s |sum|/free-field=%.4f -> %s (< 0.05)\n", "ga occluded direct", r, (r < 0.05) ? "OK" : "NG";
		exit (r < 0.05) ? 0 : 1 }' || status=1
	chk "ga clear direct amplitude" "$(ga_sum "$d/clear.txt" "$tc" 0.0008)" \
		"$(awk -v pi="$PI" 'BEGIN{printf "%.10e", 1/(4*pi*sqrt(73.0))}')" 0.01
	# 遮蔽で棄却した像の数が受音点 #1 でゼロでないこと (可視性判定が働いた証拠)
	grep -qE "point #1: image sources up to order 3 over [0-9]+ surfaces -> [0-9]+ visible, [1-9][0-9]* blocked" \
		"$d/solver.log" \
		&& say_ok "ga blocked-image count" \
		|| say_ng "ga blocked-image count" " (visibility test did not reject anything)"
}

echo "--- (H) early specular reflection off a rigid obstacle panel (+-1%)"
# 障害物 (geometry) 自身の低次鏡面反射が鏡像法に入っていること、および
# 反射点が有限面からはみ出す像は棄却されること。導出は ga_panel.ofd 参照。
ga_run ga_panel ga_panel.ofd ga_panel.ofdx && {
	d="$WORK/ga_panel"
	dump "$d/rir.wav" > "$d/under.txt"
	dump "$d/rir_beyond.wav" > "$d/beyond.txt"
	for k in 0 1; do
		r=$(awk -v k=$k 'BEGIN{printf "%.10f", (k==0)?sqrt(17.17):sqrt(77.65)}')
		te=$(awk -v r="$r" 'BEGIN{printf "%.10f", r/343.0}')
		chk "ga panel #$k amplitude" "$(ga_sum "$d/under.txt" "$te" 0.001)" \
			"$(awk -v r="$r" -v pi="$PI" 'BEGIN{printf "%.10e", 1/(4*pi*r)}')" 0.01
		chk "ga panel #$k arrival [s]" "$(ga_peak "$d/under.txt" "$te" 0.001)" "$te" 0.01
	done
	# 板の縁の外の受音点 : 直接音だけで、板反射の時刻には何も来ない
	rb=$(awk 'BEGIN{printf "%.10f", sqrt(145.17)}')
	tb=$(awk -v r="$rb" 'BEGIN{printf "%.10f", r/343.0}')
	chk "ga beyond direct" "$(ga_sum "$d/beyond.txt" "$tb" 0.001)" \
		"$(awk -v r="$rb" -v pi="$PI" 'BEGIN{printf "%.10e", 1/(4*pi*r)}')" 0.01
	tp=$(awk 'BEGIN{printf "%.10f", sqrt(144+0.81+60.84)/343.0}')
	awk -v a="$(ga_sum "$d/beyond.txt" "$tp" 0.001)" -v pi="$PI" 'BEGIN {
		ref = 1/(4*pi*sqrt(145.17)); r = (a < 0 ? -a : a) / ref;
		printf "%-26s |sum|/direct=%.4f -> %s (< 0.05, 反射点が板の外)\n",
			"ga off-panel rejected", r, (r < 0.05) ? "OK" : "NG";
		exit (r < 0.05) ? 0 : 1 }' || status=1
}

echo "--- (B2) per-surface scattering overrides the room default"
# 床だけ s = 0 (室の既定は 0.5)。既定値が使われてしまえば鏡面成分が
# sqrt(0.5) = 0.707 倍になり、上の (B) の振幅判定が ±1% を大きく外れる。
# もう一方の端 (床 s = 1 = 完全拡散) では床の鏡面像そのものが消える。
grep -q '"image_sources": 2' "$WORK/ga_floor/metadata.json" \
	&& say_ok "ga per-surface s=0" || say_ng "ga per-surface s=0" " (image_sources != 2)"
rm -rf "$WORK/ga_floor_s1"; mkdir -p "$WORK/ga_floor_s1"
cp "$SRC/ga_floor.ofd" "$WORK/ga_floor_s1/"
sed 's/"scattering": 0.0,/"scattering": 1.0,/' "$SRC/ga_floor.ofdx" \
	> "$WORK/ga_floor_s1/ga_floor.ofdx"
if "$GA_SOLVER" "$WORK/ga_floor_s1" > /dev/null 2>&1 && \
   grep -q '"image_sources": 1' "$WORK/ga_floor_s1/metadata.json"; then
	say_ok "ga per-surface s=1" " (完全拡散で床の鏡面像が消える)"
else
	say_ng "ga per-surface s=1" " (image_sources != 1)"
fi

echo "--- (I) angle-dependent absorption of a locally reacting floor (+-1%)"
# zeta は Paris の式 alpha_stat(zeta) = (8/zeta^2)[zeta+1-2ln(1+zeta)-1/(1+zeta)]
# を二分法で逆に解く (published な閉形式の独立実装)。導出は ga_angle.ofd 参照。
ga_zeta() {   # ga_zeta <alpha_stat>
	awk -v a="$1" 'function f(z) { return (8/(z*z))*(z + 1 - 2*log(1+z) - 1/(1+z)) }
		BEGIN { lo = 1.5537; hi = 1e9;
			for (i = 0; i < 200; i++) { m = 0.5*(lo+hi); if (f(m) > a) lo = m; else hi = m }
			printf "%.10f", 0.5*(lo+hi) }'
}
ga_run ga_angle ga_angle.ofd ga_angle.ofdx && {
	d="$WORK/ga_angle"
	dump "$d/rir.wav" > "$d/rir.txt"
	r1=$(awk 'BEGIN{printf "%.10f", sqrt(457.94)}')
	t1=$(awk -v r="$r1" 'BEGIN{printf "%.10f", r/343.0}')
	ct=$(awk -v r="$r1" 'BEGIN{printf "%.10f", 7.5/r}')
	zeta=$(ga_zeta 0.30)
	chk "ga angle-dep R(theta)" "$(ga_sum "$d/rir.txt" "$t1" 0.001)" \
		"$(awk -v r="$r1" -v z="$zeta" -v ct="$ct" -v pi="$PI" \
		   'BEGIN{printf "%.10e", (z*ct - 1)/(z*ct + 1)/(4*pi*r)}')" 0.01
	# 同じ入力で角度依存を off にすると R = sqrt(1-alpha) (既定 = 従来動作)
	rm -rf "$WORK/ga_angle_off"; mkdir -p "$WORK/ga_angle_off"
	cp "$SRC/ga_angle.ofd" "$WORK/ga_angle_off/"
	sed 's/"angle_dependent_absorption": true/"angle_dependent_absorption": false/' \
		"$SRC/ga_angle.ofdx" > "$WORK/ga_angle_off/ga_angle.ofdx"
	if "$GA_SOLVER" "$WORK/ga_angle_off" > /dev/null 2>&1; then
		dump "$WORK/ga_angle_off/rir.wav" > "$WORK/ga_angle_off/rir.txt"
		chk "ga angle-indep sqrt(1-a)" \
			"$(ga_sum "$WORK/ga_angle_off/rir.txt" "$t1" 0.001)" \
			"$(awk -v r="$r1" -v pi="$PI" 'BEGIN{printf "%.10e", sqrt(0.7)/(4*pi*r)}')" 0.01
	else
		say_ng "ga angle-indep sqrt(1-a)" " (solver failed)"
	fi
}

echo "--- (J) per-obstacle material via acoustic.ga.obstacles (+-1%)"
# 同じ ga_panel の板に alpha = 0.36 を与えると 1 次板反射は
# sqrt(1-0.36) = 0.8 倍になる。吸音で失われたエネルギーは (散乱と違い)
# 戻らないので判定窓は無雑音のまま。導出は ga_panel.ofd の先頭コメント。
rm -rf "$WORK/ga_panel_mat"; mkdir -p "$WORK/ga_panel_mat"
cp "$SRC/ga_panel.ofd" "$WORK/ga_panel_mat/"
sed 's/"scattering": 0.0/"scattering": 0.0, "obstacles": [ { "geometry": 1, "alpha": [0.36, 0.36, 0.36, 0.36, 0.36, 0.36] } ]/' \
	"$SRC/ga_panel.ofdx" > "$WORK/ga_panel_mat/ga_panel.ofdx"
if "$GA_SOLVER" "$WORK/ga_panel_mat" > /dev/null 2>&1; then
	d="$WORK/ga_panel_mat"
	dump "$d/rir.wav" > "$d/under.txt"
	r1=$(awk 'BEGIN{printf "%.10f", sqrt(77.65)}')
	t1=$(awk -v r="$r1" 'BEGIN{printf "%.10f", r/343.0}')
	chk "ga obstacle alpha=0.36" "$(ga_sum "$d/under.txt" "$t1" 0.001)" \
		"$(awk -v r="$r1" -v pi="$PI" 'BEGIN{printf "%.10e", sqrt(1-0.36)/(4*pi*r)}')" 0.01
	r0=$(awk 'BEGIN{printf "%.10f", sqrt(17.17)}')
	t0d=$(awk -v r="$r0" 'BEGIN{printf "%.10f", r/343.0}')
	chk "ga obstacle direct intact" "$(ga_sum "$d/under.txt" "$t0d" 0.001)" \
		"$(awk -v r="$r0" -v pi="$PI" 'BEGIN{printf "%.10e", 1/(4*pi*r)}')" 0.01
else
	say_ng "ga obstacle material run" " (solver failed)"
fi
# 板の散乱係数 = 1 (alpha 省略 = 剛体のまま) -> 板の鏡面像が拡散に回って
# 消える (受音点 #1 の image_sources が 2 -> 1。#2 はもともと 1)。
rm -rf "$WORK/ga_panel_s1"; mkdir -p "$WORK/ga_panel_s1"
cp "$SRC/ga_panel.ofd" "$WORK/ga_panel_s1/"
sed 's/"scattering": 0.0/"scattering": 0.0, "obstacles": [ { "geometry": 1, "scattering": 1.0 } ]/' \
	"$SRC/ga_panel.ofdx" > "$WORK/ga_panel_s1/ga_panel.ofdx"
if "$GA_SOLVER" "$WORK/ga_panel_s1" > /dev/null 2>&1 && \
   ! grep -q '"image_sources": 2,' "$WORK/ga_panel_s1/metadata.json"; then
	say_ok "ga obstacle scattering=1" " (板の鏡面像が消える)"
else
	say_ng "ga obstacle scattering=1" " (panel image should vanish)"
fi
# 存在しない geometry を指す行は非零終了 (黙って無視しない)
rm -rf "$WORK/ga_badobs"; mkdir -p "$WORK/ga_badobs"
cp "$SRC/ga_panel.ofd" "$WORK/ga_badobs/"
sed 's/"scattering": 0.0/"scattering": 0.0, "obstacles": [ { "geometry": 9 } ]/' \
	"$SRC/ga_panel.ofdx" > "$WORK/ga_badobs/ga_panel.ofdx"
if "$GA_SOLVER" "$WORK/ga_badobs" > /dev/null 2> "$WORK/ga_badobs/stderr.log"; then
	say_ng "ga bad obstacle index" " (should fail)"
elif grep -qi "geometry" "$WORK/ga_badobs/stderr.log"; then
	printf "%-26s -> OK (refused with the reason)\n" "ga bad obstacle index"
else
	say_ng "ga bad obstacle index" " (no reason in stderr)"
fi

echo "--- (K) multi-source: two feeds, two closed-form arrivals (+-1%)"
# 契約 (OpenFDTD-X ADR-0010) : multi_source で全 feed が強度 1・t = 0 で
# 同時発火し、rir.wav は重ね合わせ。導出は ga_two.ofd の先頭コメント。
ga_run ga_two ga_two.ofd ga_two.ofdx && {
	d="$WORK/ga_two"
	dump "$d/rir.wav" > "$d/rir.txt"
	for k in 1 2; do
		r=$(awk -v k=$k 'BEGIN{print (k==1)?10.0:20.0}')
		te=$(awk -v r="$r" 'BEGIN{printf "%.10f", r/343.0}')
		chk "ga two-src #$k amplitude" "$(ga_sum "$d/rir.txt" "$te" 0.001)" \
			"$(awk -v r="$r" -v pi="$PI" 'BEGIN{printf "%.10e", 1/(4*pi*r)}')" 0.01
		chk "ga two-src #$k arrival" "$(ga_peak "$d/rir.txt" "$te" 0.001)" "$te" 0.01
	done
	grep -q '"multi_source": true' "$d/metadata.json" \
		&& say_ok "ga multi_source metadata" || say_ng "ga multi_source metadata"
	grep -Fq '{ "pos_m": [35, 10, 10] }' "$d/metadata.json" \
		&& say_ok "ga sources listed" || say_ng "ga sources listed"
}
# 既定 (multi_source なし) は feed #1 のみ + warning (後方互換)
rm -rf "$WORK/ga_two_def"; mkdir -p "$WORK/ga_two_def"
cp "$SRC/ga_two.ofd" "$WORK/ga_two_def/"
sed '/"multi_source"/d' "$SRC/ga_two.ofdx" > "$WORK/ga_two_def/ga_two.ofdx"
if "$GA_SOLVER" "$WORK/ga_two_def" > /dev/null 2>&1; then
	grep -q "warning: 2 feeds" "$WORK/ga_two_def/solver.log" \
		&& say_ok "ga default single-source" " (warning)" \
		|| say_ng "ga default single-source" " (no warning)"
	dump "$WORK/ga_two_def/rir.wav" > "$WORK/ga_two_def/rir.txt"
	t2=$(awk 'BEGIN{printf "%.10f", 20.0/343.0}')
	awk -v a="$(ga_sum "$WORK/ga_two_def/rir.txt" "$t2" 0.001)" -v pi="$PI" 'BEGIN {
		ref = 1/(4*pi*20.0); rr = (a < 0 ? -a : a) / ref;
		printf "%-26s |sum|/would-be=%.4f -> %s (< 0.05, feed #2 不使用)\n", "ga default no 2nd arrival", rr, (rr < 0.05) ? "OK" : "NG";
		exit (rr < 0.05) ? 0 : 1 }' || status=1
else
	say_ng "ga default single-source" " (solver failed)"
fi
# multi_source では全音源が室内になければならない (室外は非零終了)
rm -rf "$WORK/ga_two_out"; mkdir -p "$WORK/ga_two_out"
sed 's/^feed = z 35.0/feed = z 150.0/' "$SRC/ga_two.ofd" > "$WORK/ga_two_out/ga_two.ofd"
cp "$SRC/ga_two.ofdx" "$WORK/ga_two_out/"
if "$GA_SOLVER" "$WORK/ga_two_out" > /dev/null 2> "$WORK/ga_two_out/stderr.log"; then
	say_ng "ga multi-src outside exit" " (should fail)"
elif grep -qi "feed #2" "$WORK/ga_two_out/stderr.log"; then
	printf "%-26s -> OK (feed #2 を指して拒否)\n" "ga multi-src outside exit"
else
	say_ng "ga multi-src outside exit" " (no feed # in stderr)"
fi

echo "--- (F) ga determinism (bit-identical reruns and thread invariance)"
# 乱数生成器を使わない (レイ方向は球面フィボナッチ格子、拡散反射の抽選は
# 決定的な整数ハッシュ、後期残響の符号列は固定初期値の LFSR)。
# ソルバー自体は OpenMP を使わないので、スレッド数にも依存しない。
d="$WORK/ga_lossless"
cp "$d/rir.wav" "$d/rir_run1.wav"
if "$GA_SOLVER" "$d" > /dev/null 2>&1 && cmp -s "$d/rir_run1.wav" "$d/rir.wav"; then
	printf "%-26s -> OK (2 回実行がビット一致)\n" "ga rerun determinism"
else
	say_ng "ga rerun determinism"
fi
OMP_NUM_THREADS=1 "$GA_SOLVER" "$d" > /dev/null 2>&1
cp "$d/rir.wav" "$d/rir_n1.wav"
OMP_NUM_THREADS=4 "$GA_SOLVER" "$d" > /dev/null 2>&1
if cmp -s "$d/rir_n1.wav" "$d/rir.wav"; then
	printf "%-26s -> OK (1 と 4 スレッドがビット一致)\n" "ga thread invariance"
else
	say_ng "ga thread invariance"
fi

echo "--- (G) ga error paths (honest non-zero exit, no fabricated output)"
rm -rf "$WORK/ga_empty"; mkdir -p "$WORK/ga_empty"
if "$GA_SOLVER" "$WORK/ga_empty" > /dev/null 2> "$WORK/ga_empty/stderr.log"; then
	say_ng "ga no-input exit" " (should fail)"
elif [ -s "$WORK/ga_empty/stderr.log" ] && [ ! -f "$WORK/ga_empty/rir.wav" ]; then
	printf "%-26s -> OK (non-zero, stderr reason, no rir.wav)\n" "ga no-input exit"
else
	say_ng "ga no-input exit" " (missing stderr reason or fabricated rir)"
fi
# .ofdx の値域外 (image_order = 9) は黙って既定値に落とさず非零終了する
rm -rf "$WORK/ga_badorder"; mkdir -p "$WORK/ga_badorder"
cp "$SRC/ga_floor.ofd" "$WORK/ga_badorder/"
sed 's/"image_order": 2/"image_order": 9/' "$SRC/ga_floor.ofdx" \
	> "$WORK/ga_badorder/ga_floor.ofdx"
if "$GA_SOLVER" "$WORK/ga_badorder" > /dev/null 2> "$WORK/ga_badorder/stderr.log"; then
	say_ng "ga bad image_order exit" " (should fail)"
elif grep -qi "image_order" "$WORK/ga_badorder/stderr.log" && \
     [ ! -f "$WORK/ga_badorder/rir.wav" ]; then
	printf "%-26s -> OK (refused with the offending key)\n" "ga bad image_order exit"
else
	say_ng "ga bad image_order exit" " (no reason in stderr)"
fi
# 音源が室外 (数値を捏造せず失敗する)
rm -rf "$WORK/ga_outside"; mkdir -p "$WORK/ga_outside"
sed 's/^feed = z 3.0 5.0 5.0/feed = z 30.0 5.0 5.0/' "$SRC/ga_floor.ofd" \
	> "$WORK/ga_outside/ga_floor.ofd"
if "$GA_SOLVER" "$WORK/ga_outside" > /dev/null 2> "$WORK/ga_outside/stderr.log"; then
	say_ng "ga source-outside exit" " (should fail)"
elif grep -qi "outside the room" "$WORK/ga_outside/stderr.log"; then
	printf "%-26s -> OK (refused with the reason)\n" "ga source-outside exit"
else
	say_ng "ga source-outside exit" " (no reason in stderr)"
fi
if "$GA_SOLVER" > /dev/null 2>&1; then
	say_ng "ga usage exit" " (should fail)"
else
	say_ok "ga usage exit"
fi

if [ "$status" -ne 0 ]; then
	echo "*** acoustic validation FAILED" >&2
else
	echo "acoustic validation passed"
fi
exit $status
