/* ga_main.c — ofdx_acoustic_ga エントリ (ADR-0007 出力契約の I/O と進捗)
 *
 * 起動 : ofdx_acoustic_ga <working_dir> [<input_file.ofd>]
 *   - input_file 省略時は working_dir 直下の唯一の .ofd を探す。
 * 出力 (working_dir 直下) :
 *   rir.wav (+ rir_<名前>.wav) / metadata.json / metrics.json / solver.log
 * 進捗 : stdout に "progress a/b" 行 (AcousticRunner が解析する書式 —
 *        ofdx_acoustic_fdtd および OpenFDTD-X の mock_acoustic_solver.c と同型)。
 * 異常時 : 非零終了コード + stderr へ明確な理由 (合成データは出力しない)。
 *
 * ofdx_acoustic_fdtd と同じ CLI・同じ契約だが、中身は完全に独立した
 * 別バイナリ (FDTD 側のソースには手を入れない)。共有するのは src/wav.c のみ。
 */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ga.h"

#ifdef _WIN32
#include <windows.h>   /* MultiByteToWideChar (ac_fopen の UTF-8 → UTF-16) */
#endif

/* UTF-8 パス対応の fopen。Windows の fopen は ANSI コードページで解釈するため
 * 日本語を含むパス・ファイル名が開けない (GUI の受音点名は既定で「マイク N」)。
 * 名前が ac_ 接頭辞なのは src/wav.c を FDTD 側と共用しているため (ga.h の
 * 宣言参照)。FDTD の main.c とは別バイナリなのでシンボルは衝突しない。 */
FILE *ac_fopen(const char *path, const char *mode)
{
#ifdef _WIN32
	wchar_t wpath[GA_PATH_MAX + 160];
	wchar_t wmode[8];
	int n = MultiByteToWideChar(CP_UTF8, 0, path, -1, wpath,
	                            (int)(sizeof(wpath) / sizeof(wpath[0])));
	if (n <= 0) return fopen(path, mode);   /* 変換できなければ従来経路 */
	n = MultiByteToWideChar(CP_UTF8, 0, mode, -1, wmode,
	                        (int)(sizeof(wmode) / sizeof(wmode[0])));
	if (n <= 0) return fopen(path, mode);
	return _wfopen(wpath, wmode);
#else
	return fopen(path, mode);
#endif
}

void ga_log(ga_t *g, const char *fmt, ...)
{
	va_list ap;
	FILE *fp = (g && g->logfp) ? g->logfp : stderr;
	va_start(ap, fmt);
	vfprintf(fp, fmt, ap);
	va_end(ap);
	fputc('\n', fp);
	fflush(fp);
}

void ga_err(ga_t *g, const char *fmt, ...)
{
	va_list ap;
	fprintf(stderr, "%s: error: ", GA_SOLVER_NAME);
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fputc('\n', stderr);
	fflush(stderr);
	if (g && g->logfp) {
		fprintf(g->logfp, "error: ");
		va_start(ap, fmt);
		vfprintf(g->logfp, fmt, ap);
		va_end(ap);
		fputc('\n', g->logfp);
		fflush(g->logfp);
	}
}

void ga_free(ga_t *g)
{
	free(g->geom); g->geom = NULL;
	free(g->surf); g->surf = NULL;
	free(g->feedpos); g->feedpos = NULL;
	free(g->srcgain); g->srcgain = NULL;
	free(g->srcdelay); g->srcdelay = NULL;
	free(g->recv); g->recv = NULL;
	free(g->echo); g->echo = NULL;
	free(g->rir);  g->rir  = NULL;
	free(g->band); g->band = NULL;
	free(g->fr);   g->fr   = NULL;
	free(g->fi);   g->fi   = NULL;
	free(g->yr);   g->yr   = NULL;
	free(g->yi);   g->yi   = NULL;
	free(g->twr);  g->twr  = NULL;
	free(g->twi);  g->twi  = NULL;
	free(g->brev); g->brev = NULL;
}

static int join_path(char *out, size_t cap, const char *dir, const char *name)
{
	int r = snprintf(out, cap, "%s/%s", dir, name);
	return (r > 0 && (size_t)r < cap) ? 0 : -1;
}

/* 受音点ごとの WAV 書き出し (契約 : #1 は rir.wav) */
static int write_rirs(ga_t *g)
{
	char path[GA_PATH_MAX + 128];
	int r;
	for (r = 0; r < g->nrecv; r++) {
		const double *x = g->rir + (size_t)r * g->nsamples;
		if (join_path(path, sizeof(path), g->workdir, g->recv[r].file) != 0) {
			ga_err(g, "output path too long for receiver #%d", r + 1);
			return 1;
		}
		if (ac_write_wav_f32(path, GA_FS, x, g->nsamples) != 0) {
			ga_err(g, "cannot write %s", path);
			return 1;
		}
		ga_log(g, "output: %s (%d samples, float32 %d Hz)",
		       g->recv[r].file, g->nsamples, GA_FS);
		if (g->recv[r].alias[0] != '\0') {
			if (join_path(path, sizeof(path), g->workdir,
			              g->recv[r].alias) != 0) {
				ga_err(g, "output path too long for receiver #%d alias", r + 1);
				return 1;
			}
			if (ac_write_wav_f32(path, GA_FS, x, g->nsamples) != 0) {
				ga_err(g, "cannot write %s", path);
				return 1;
			}
			ga_log(g, "output: %s (alias of %s — 名前照合の自動割当用)",
			       g->recv[r].alias, g->recv[r].file);
		}
	}
	return 0;
}

int main(int argc, char **argv)
{
	static ga_t ga;               /* 状態はこの 1 個 (BSS — スタックに置かない) */
	char logpath[GA_PATH_MAX + 16];
	const char *args[2] = { NULL, NULL };
	int i, w, b, nargs = 0;

	memset(&ga, 0, sizeof(ga));
	for (w = 0; w < GA_NWALL; w++) {
		for (b = 0; b < GA_NBAND; b++)
			ga.alpha[w][b] = GA_ALPHA_DEFAULT;
		ga.wall_scatter[w] = -1.0;    /* < 0 = acoustic.ga.scattering を使う */
	}
	ga.order    = GA_ORDER_DEFAULT;
	ga.nrays    = GA_RAYS_DEFAULT;
	ga.scatter  = GA_SCATTER_DEFAULT;
	ga.temp_c   = GA_TEMP_C_DEFAULT;
	ga.humid    = GA_HUMID_DEFAULT;
	ga.press_kpa = GA_PRESS_KPA_DEFAULT;
	ga.air_on   = 1;

	/* 引数 : "-" で始まるものは無視 (ランチャ由来のオプションを許容) */
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			fprintf(stderr, "%s: warning: ignoring option '%s'\n",
			        GA_SOLVER_NAME, argv[i]);
			continue;
		}
		if (nargs < 2) args[nargs] = argv[i];
		nargs++;
	}
	if (nargs < 1 || nargs > 2) {
		fprintf(stderr, "usage: %s <working_dir> [<input_file.ofd>]\n",
		        GA_SOLVER_NAME);
		return 2;
	}
	if (strlen(args[0]) >= sizeof(ga.workdir)) {
		fprintf(stderr, "%s: error: working_dir path too long\n", GA_SOLVER_NAME);
		return 1;
	}
	strcpy(ga.workdir, args[0]);

	if (join_path(logpath, sizeof(logpath), ga.workdir, "solver.log") != 0 ||
	    (ga.logfp = ac_fopen(logpath, "w")) == NULL) {
		fprintf(stderr, "%s: error: cannot write solver.log in '%s' "
		        "(directory missing or not writable)\n",
		        GA_SOLVER_NAME, ga.workdir);
		return 1;
	}
	ga_log(&ga, "%s %s (contract_version %d)",
	       GA_SOLVER_NAME, GA_SOLVER_VERSION, GA_CONTRACT_VERSION);
	ga_log(&ga, "working dir: %s", ga.workdir);
	printf("%s %s\n", GA_SOLVER_NAME, GA_SOLVER_VERSION);
	fflush(stdout);

	if (nargs == 2) {
		FILE *fp = ac_fopen(args[1], "rb");
		if (fp) {
			fclose(fp);
			if (strlen(args[1]) >= sizeof(ga.ofd_path)) {
				ga_err(&ga, "input path too long");
				goto failed;
			}
			strcpy(ga.ofd_path, args[1]);
		}
		else if (join_path(ga.ofd_path, sizeof(ga.ofd_path),
		                   ga.workdir, args[1]) != 0 ||
		         (fp = ac_fopen(ga.ofd_path, "rb")) == NULL) {
			ga_err(&ga, "input file '%s' not found (neither as given nor in %s)",
			       args[1], ga.workdir);
			goto failed;
		}
		else {
			fclose(fp);
		}
	}
	else if (ga_find_input(&ga) != 0) {
		goto failed;
	}

	/* metadata 用ファイル名と .ofdx パス (同 basename の JSON サイドカー) */
	{
		const char *base = strrchr(ga.ofd_path, '/');
		const char *b2 = strrchr(ga.ofd_path, '\\');
		size_t len;
		if (b2 && (!base || b2 > base)) base = b2;
		base = base ? base + 1 : ga.ofd_path;
		snprintf(ga.ofd_name, sizeof(ga.ofd_name), "%.*s",
		         (int)sizeof(ga.ofd_name) - 1, base);
		len = strlen(ga.ofd_path);
		if (len >= 4 && strcmp(ga.ofd_path + len - 4, ".ofd") == 0)
			snprintf(ga.ofdx_path, sizeof(ga.ofdx_path), "%.*sofdx",
			         (int)(len - 3), ga.ofd_path);
		else
			snprintf(ga.ofdx_path, sizeof(ga.ofdx_path), "%s.ofdx", ga.ofd_path);
	}
	ga_log(&ga, "input: %s", ga.ofd_path);

	if (ga_read_ofd(&ga) != 0)  goto failed;
	if (ga_read_ofdx(&ga) != 0) goto failed;
	if (ga_setup(&ga) != 0)     goto failed;
	if (ga_rays(&ga) != 0)      goto failed;
	if (ga_synth(&ga) != 0)     goto failed;

	if (write_rirs(&ga) != 0)      goto failed;
	if (ga_write_metadata(&ga) != 0) goto failed;
	if (ga_write_metrics(&ga) != 0)  goto failed;

	ga_log(&ga, "normal end");
	printf("%s: done\n", GA_SOLVER_NAME);
	fclose(ga.logfp);
	ga.logfp = NULL;
	ga_free(&ga);
	return 0;

failed:
	if (ga.logfp) { fclose(ga.logfp); ga.logfp = NULL; }
	ga_free(&ga);
	return 1;
}
