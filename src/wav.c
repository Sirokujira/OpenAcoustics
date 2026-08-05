/* wav.c — float32 モノラル WAV 書き出し (自前実装 — 外部ライブラリ禁止)
 *
 * OpenFDTD-X/tests/acoustics/mock_acoustic_solver.c の write_wav_header と
 * 同型の 44 byte 標準ヘッダ (RIFF + fmt(16) + data)。フォーマットは
 * WAVE_FORMAT_IEEE_FLOAT (= 3)、32 bit。GUI 側 WavReader は format tag 3 /
 * 32 bit float に対応している (src/acoustics/io/WavReader.cpp)。
 * バイト順は WAV 仕様どおりリトルエンディアンに明示シリアライズする
 * (ホストのエンディアンに依存しない)。
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "acoustic.h"

static void put_u32le(unsigned char *p, uint32_t v)
{
	p[0] = (unsigned char)(v & 0xffu);
	p[1] = (unsigned char)((v >> 8) & 0xffu);
	p[2] = (unsigned char)((v >> 16) & 0xffu);
	p[3] = (unsigned char)((v >> 24) & 0xffu);
}

static void put_u16le(unsigned char *p, uint16_t v)
{
	p[0] = (unsigned char)(v & 0xffu);
	p[1] = (unsigned char)((v >> 8) & 0xffu);
}

int ac_write_wav_f32(const char *path, int fs, const double *x, int n)
{
	unsigned char h[44];
	unsigned char buf[4096 * 4];
	const uint16_t channels = 1, bits = 32;
	const uint16_t block_align = (uint16_t)(channels * bits / 8);
	const uint32_t byte_rate = (uint32_t)fs * block_align;
	const uint32_t data_size = (uint32_t)n * block_align;
	FILE *fp;
	int i, m;

	memcpy(h + 0, "RIFF", 4);
	put_u32le(h + 4, 36u + data_size);
	memcpy(h + 8, "WAVE", 4);
	memcpy(h + 12, "fmt ", 4);
	put_u32le(h + 16, 16u);          /* fmt チャンクサイズ */
	put_u16le(h + 20, 3u);           /* WAVE_FORMAT_IEEE_FLOAT */
	put_u16le(h + 22, channels);
	put_u32le(h + 24, (uint32_t)fs);
	put_u32le(h + 28, byte_rate);
	put_u16le(h + 32, block_align);
	put_u16le(h + 34, bits);
	memcpy(h + 36, "data", 4);
	put_u32le(h + 40, data_size);

	fp = ac_fopen(path, "wb");
	if (!fp) return -1;
	if (fwrite(h, 1, sizeof(h), fp) != sizeof(h)) {
		fclose(fp);
		return -1;
	}
	for (i = 0; i < n; i += m) {
		int b;
		m = n - i;
		if (m > 4096) m = 4096;
		for (b = 0; b < m; b++) {
			float f = (float)x[i + b];
			uint32_t bitsval;
			memcpy(&bitsval, &f, 4);          /* 型 punning は memcpy で */
			put_u32le(buf + 4 * b, bitsval);
		}
		if (fwrite(buf, 4, (size_t)m, fp) != (size_t)m) {
			fclose(fp);
			return -1;
		}
	}
	return (fclose(fp) == 0) ? 0 : -1;
}
