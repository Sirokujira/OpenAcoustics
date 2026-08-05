#!/bin/sh
# acoustic_check.sh — ofdx_acoustic_fdtd 検証 (CI 用、OpenPEEC peec_check.sh の流儀)
#
# data/sample/ の各ケースを実行し、rir.wav を解析解と比較する。
# 期待値の導出は各 .ofd ファイルのコメント参照 (すべて実装から独立な解析解)。
#
# 判定 :
#  (a) 剛体閉箱 4x3x2.5 m の軸モード f_100/f_010/f_001 が ±3% (走査 DFT)
#  (b) 準 1D 管の端面反射 R = sqrt(1-alpha) が ±3% (alpha=0.3, 0.9)
#  (c) 直接音到達 : |p| ピークが t = t0 + r/c ± (sigma + 2/fs)
#  (d) 剛壁のみ (alpha=0) でエネルギー保存 (1.4 s 窓 2 つの減衰 < 0.5 dB)
#  (e) 決定性 : 同一入力 2 回 / OMP_NUM_THREADS=1 と 4 で rir.wav ビット一致
#  (f) 契約 : 出力 4 ファイル、WAV ヘッダ、progress 行書式、異常系の非零終了
#
# WAV の読みは od -t f4 (float32 リトルエンディアン)。CI の 3 OS
# (Linux / macOS / Windows Git Bash) はいずれもリトルエンディアンかつ
# coreutils/BSD od が -j/-t f4 に対応している。
#
# 使い方 : acoustic_check.sh <ofdx_acoustic_fdtd 実行ファイル(絶対パス)> [作業ディレクトリ]

set -e

SOLVER="$1"
WORK="${2:-.}"
SRC="$(cd "$(dirname "$0")" && pwd)"
PI=3.14159265358979

if [ -z "$SOLVER" ]; then
	echo "Usage: acoustic_check.sh <ofdx_acoustic_fdtd> [workdir]" >&2
	exit 2
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

if [ "$status" -ne 0 ]; then
	echo "*** acoustic validation FAILED" >&2
else
	echo "acoustic validation passed"
fi
exit $status
