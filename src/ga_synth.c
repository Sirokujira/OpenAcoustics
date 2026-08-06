/* ga_synth.c — バンド別インパルス列 + エコーグラム -> RIR (48 kHz float32)
 *
 * ── 反射 1 発の置き方 (ga_deposit) ───────────────────────────────
 * 幾何音響の反射は「時刻 t、振幅 A の理想インパルス」なので、離散時系列へは
 * 3 次ラグランジュ補間 (4 タップ) の分数遅延で置く。ラグランジュ補間は定数と
 * 1 次式を厳密に再現するため
 *   sum_n h[n] = 1  (DC 利得 = 1)、sum_n n h[n] = t*fs  (1 次モーメント厳密)
 * が成り立つ。したがって RIR 上でこの到達の**標本和が振幅 A に等しく**、
 * 重心が正確に t に来る。離散 RIR の「振幅」は畳み込みで再現される量、
 * すなわち単位標本利得 (= 標本和) のことなので、これが振幅規約
 * A = 1/(4 pi r) の厳密な離散表現になる (検証 (a)(b) はこの和を測る)。
 *
 * ── バンド合成 ────────────────────────────────────────────────────
 * バンド重み W_b(f) はオクターブ中心 fc_b = 125*2^b を節点とする
 * 「対数周波数上の raised-cosine ハット」:
 *   W_b(f) = (1 + cos(pi |log2(f/fc_b)|)) / 2   (|log2(f/fc_b)| <= 1)
 *   W_0(f) = 1 (f <= 125 Hz)、W_6(f) = 1 (f >= 8 kHz)
 * この族は
 *   (1) sum_b W_b(f) = 1  (単位分割 — 全バンド同利得なら入力を厳密に復元)
 *   (2) W_b(fc_b') = delta_bb'  (バンド中心では他バンドが漏れない)
 * を厳密に満たす。(1) は「吸音・空気吸収が周波数によらない場合に理想
 * インパルスがそのまま残る」ことを保証し、(2) は検証 (d) の 4 kHz 成分の
 * 読み取りを厳密にする。零位相 (実数の W_b) なので群遅延は 0 で、到達時刻に
 * 遅れを持ち込まない (t = 0 は音源発火時刻という規約を壊さない)。
 * 実装は自前の radix-2 FFT (外部ライブラリ禁止) で
 *   Y(f) = sum_b W_b(f) X_b(f)
 * を作り、逆変換して RIR にする。零詰め長は nsamples + 4096 以上の 2 冪なので
 * 巡回畳み込みの回り込みは出力範囲に入らない。
 *
 * ── 後期残響 (エコーグラム -> 雑音) ──────────────────────────────
 * ビン bin のバンド別エネルギー E_b (= 到達振幅の 2 乗和) を、そのビン内の
 * 標本へ等エネルギーの雑音として配る :
 *   x_b[n] = s[n] * sqrt(D_b(t_n) / fs),  D_b = E_b / bin幅  [エネルギー/秒]
 * とすると sum_n x_b[n]^2 = E_b が厳密に成り立つ。D_b はビン中心間で線形補間
 * するので包絡は滑らかになる。s[n] は **決定的** な ±1 系列 (32 bit 最大長
 * LFSR、固定初期値) で、乱数ではない — 再実行・スレッド数によらずビット一致
 * する。全バンドで同一の s[n] を使うのが要点で、これにより W_b の単位分割
 * (1) がそのまま効き、バンド境界でも振幅が滑らかに補間される。
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ga.h"

#define GA_PROG_RAYS  45
#define GA_PROG_TOTAL 50

/* ── 反射 1 発を置く (3 次ラグランジュ分数遅延) ──────────────────── */
void ga_deposit(ga_t *g, double *band, double t, const double *amp)
{
	double x = t * (double)GA_FS;
	double n0d = floor(x);
	double frac = x - n0d;
	double D = 1.0 + frac;          /* タップ n0-1 を原点とした遅延 (1 <= D < 2) */
	double h[4];
	int n0 = (int)n0d, k, b;

	if (!(x >= 0.0) || x >= (double)g->nsamples) return;
	h[0] = (D - 1.0) * (D - 2.0) * (D - 3.0) / (-6.0);
	h[1] = (D - 0.0) * (D - 2.0) * (D - 3.0) / (2.0);
	h[2] = (D - 0.0) * (D - 1.0) * (D - 3.0) / (-2.0);
	h[3] = (D - 0.0) * (D - 1.0) * (D - 2.0) / (6.0);
	for (k = 0; k < 4; k++) {
		int idx = n0 - 1 + k;
		if (idx < 0 || idx >= g->nsamples) continue;
		for (b = 0; b < GA_NBAND; b++)
			band[(size_t)b * g->nsamples + idx] += amp[b] * h[k];
	}
}

/* ── バンド重み (単位分割の raised-cosine ハット) ────────────────── */
static double band_weight(int b, double f)
{
	double fc = GA_BAND_F0 * (double)(1 << b);
	double u;
	if (b == 0 && f <= fc) return 1.0;
	if (b == GA_NBAND - 1 && f >= fc) return 1.0;
	if (f <= 0.5 * fc || f >= 2.0 * fc) return 0.0;
	u = log(f / fc) / log(2.0);
	if (u < 0.0) u = -u;
	return 0.5 * (1.0 + cos(GA_PI * u));
}

/* ── radix-2 FFT (自前実装。テーブルは setup 済みの配列を使う) ───── */
static void fft_tables(ga_t *g)
{
	int n = g->fftn, i, k, bits = 0;
	while ((1 << bits) < n) bits++;
	for (i = 0; i < n; i++) {
		int rev = 0;
		for (k = 0; k < bits; k++)
			if (i & (1 << k)) rev |= 1 << (bits - 1 - k);
		g->brev[i] = rev;
	}
	for (i = 0; i < n / 2; i++) {
		double a = -2.0 * GA_PI * (double)i / (double)n;
		g->twr[i] = cos(a);
		g->twi[i] = sin(a);
	}
}

static void fft_run(ga_t *g, double *re, double *im, int inverse)
{
	int n = g->fftn, i, j, len, half, step;
	for (i = 0; i < n; i++) {
		int r = g->brev[i];
		if (r > i) {
			double t;
			t = re[i]; re[i] = re[r]; re[r] = t;
			t = im[i]; im[i] = im[r]; im[r] = t;
		}
	}
	for (len = 2; len <= n; len <<= 1) {
		half = len >> 1;
		step = n / len;
		for (i = 0; i < n; i += len) {
			for (j = 0; j < half; j++) {
				double wr = g->twr[j * step];
				double wi = inverse ? -g->twi[j * step] : g->twi[j * step];
				double ar = re[i + j], ai = im[i + j];
				double br = re[i + j + half], bi = im[i + j + half];
				double vr = br * wr - bi * wi;
				double vi = br * wi + bi * wr;
				re[i + j] = ar + vr; im[i + j] = ai + vi;
				re[i + j + half] = ar - vr; im[i + j + half] = ai - vi;
			}
		}
	}
}

/* ── 決定的 ±1 系列 (32 bit 最大長 LFSR : x^32 + x^22 + x^2 + x + 1) ──
 * 乱数ではない。初期値は固定で、受音点ごとに決まった回数だけ進めるので
 * 再実行・スレッド数・OS によらず同じ系列になる。 */
static unsigned int lfsr_step(unsigned int *s)
{
	unsigned int bit = ((*s >> 0) ^ (*s >> 10) ^ (*s >> 30) ^ (*s >> 31)) & 1u;
	*s = (*s >> 1) | (bit << 31);
	return *s;
}

/* ── 後期残響をバンド別インパルス列へ ───────────────────────────── */
static void synth_tail(ga_t *g, int ri, double *band)
{
	const double *echo = g->echo + (size_t)ri * g->nbins * GA_NBAND;
	unsigned int st = 0x13579BDFu + (unsigned int)ri * 0x9E3779B9u;
	int n, b, i;

	if (st == 0u) st = 0x1u;
	for (i = 0; i < 64 + ri * 977; i++) lfsr_step(&st);   /* 受音点ごとに位相をずらす */

	for (n = 0; n < g->nsamples; n++) {
		double t = (double)n / GA_FS;
		double u = t / GA_BIN_S - 0.5;
		double fu = floor(u);
		int i0 = (int)fu;
		double fr = u - fu;
		double sgn = (lfsr_step(&st) & 1u) ? 1.0 : -1.0;
		int i1;
		if (i0 < 0) { i0 = 0; fr = 0.0; }
		if (i0 > g->nbins - 1) { i0 = g->nbins - 1; fr = 0.0; }
		i1 = (i0 + 1 < g->nbins) ? i0 + 1 : i0;
		for (b = 0; b < GA_NBAND; b++) {
			double d = (1.0 - fr) * echo[(size_t)i0 * GA_NBAND + b]
			         + fr * echo[(size_t)i1 * GA_NBAND + b];
			if (d <= 0.0) continue;
			band[(size_t)b * g->nsamples + n] +=
				sgn * sqrt(d / (GA_BIN_S * GA_FS));
		}
	}
}

int ga_synth(ga_t *g)
{
	int r, b, n, prog = GA_PROG_RAYS;

	fft_tables(g);

	for (r = 0; r < g->nrecv; r++) {
		double *out = g->rir + (size_t)r * g->nsamples;
		double peak = 0.0;
		int k;

		memset(g->band, 0, (size_t)GA_NBAND * g->nsamples * sizeof(double));
		if (ga_images(g, r, g->band) != 0) return 1;
		synth_tail(g, r, g->band);

		memset(g->yr, 0, (size_t)g->fftn * sizeof(double));
		memset(g->yi, 0, (size_t)g->fftn * sizeof(double));
		for (b = 0; b < GA_NBAND; b++) {
			const double *x = g->band + (size_t)b * g->nsamples;
			double sum = 0.0;
			for (n = 0; n < g->nsamples; n++) sum += fabs(x[n]);
			if (sum <= 0.0) continue;      /* 空のバンドは FFT を省く */
			memcpy(g->fr, x, (size_t)g->nsamples * sizeof(double));
			memset(g->fr + g->nsamples, 0,
			       (size_t)(g->fftn - g->nsamples) * sizeof(double));
			memset(g->fi, 0, (size_t)g->fftn * sizeof(double));
			fft_run(g, g->fr, g->fi, 0);
			for (k = 0; k < g->fftn; k++) {
				/* 負の周波数は折り返して同じ重みを掛ける (出力を実数に保つ) */
				int kk = (k <= g->fftn / 2) ? k : (g->fftn - k);
				double f = (double)kk * GA_FS / (double)g->fftn;
				double W = band_weight(b, f);
				if (W == 0.0) continue;
				g->yr[k] += W * g->fr[k];
				g->yi[k] += W * g->fi[k];
			}
		}
		fft_run(g, g->yr, g->yi, 1);
		for (n = 0; n < g->nsamples; n++) {
			double v = g->yr[n] / (double)g->fftn;
			out[n] = v;
			if (fabs(v) > peak) peak = fabs(v);
		}
		if (peak <= 0.0)
			ga_log(g, "warning: point #%d received no energy (all paths blocked "
			       "or fully absorbed) — rir is silent", r + 1);
		ga_log(g, "point #%d: rir peak |p| = %.6g (%d samples)",
		       r + 1, peak, g->nsamples);

		{
			int q = GA_PROG_RAYS
			      + (int)(((long long)(r + 1) * (GA_PROG_TOTAL - GA_PROG_RAYS))
			              / g->nrecv);
			while (prog < q) {
				++prog;
				printf("progress %d/%d\n", prog, GA_PROG_TOTAL);
				fflush(stdout);
			}
		}
	}
	while (prog < GA_PROG_TOTAL) {
		++prog;
		printf("progress %d/%d\n", prog, GA_PROG_TOTAL);
		fflush(stdout);
	}
	return 0;
}
