/* main.c — ofdx_acoustic_fdtd エントリ (ADR-0007 出力契約の I/O と進捗)
 *
 * 起動 : ofdx_acoustic_fdtd <working_dir> [<input_file.ofd>]
 *   - input_file 省略時は working_dir 直下の唯一の .ofd を探す。
 *   - mpiexec 経由の起動も許容 (MPI 非対応 — 単プロセスとして動く)。
 * 出力 (working_dir 直下) :
 *   rir.wav (+ rir_<名前>.wav) / metadata.json / metrics.json / solver.log
 * 進捗 : stdout に "progress a/b" 行 (正規表現 ^progress ¥s+(¥d+)¥s*[/]¥s*(¥d+)$ の
 *        ¥ を \ に読み替え。AcousticRunner が解析する —
 *        OpenFDTD-X/tests/acoustics/mock_acoustic_solver.c と同型)。
 * 異常時 : 非零終了コード + stderr へ明確な理由 (合成データは出力しない)。
 */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "acoustic.h"

/* ログ : solver.log へ (開けていなければ黙って捨てない — stderr へ) */
void ac_log(ac_t *ac, const char *fmt, ...)
{
	va_list ap;
	FILE *fp = (ac && ac->logfp) ? ac->logfp : stderr;
	va_start(ap, fmt);
	vfprintf(fp, fmt, ap);
	va_end(ap);
	fputc('\n', fp);
	fflush(fp);
}

/* エラー : stderr と solver.log の両方へ */
void ac_err(ac_t *ac, const char *fmt, ...)
{
	va_list ap;
	fprintf(stderr, "%s: error: ", AC_SOLVER_NAME);
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fputc('\n', stderr);
	fflush(stderr);
	if (ac && ac->logfp) {
		fprintf(ac->logfp, "error: ");
		va_start(ap, fmt);
		vfprintf(ac->logfp, fmt, ap);
		va_end(ap);
		fputc('\n', ac->logfp);
		fflush(ac->logfp);
	}
}

static int join_path(char *out, size_t cap, const char *dir, const char *name)
{
	int r = snprintf(out, cap, "%s/%s", dir, name);
	return (r > 0 && (size_t)r < cap) ? 0 : -1;
}

/* 受音点ごとの WAV 書き出し */
static int write_rirs(ac_t *ac)
{
	char path[AC_PATH_MAX + 128];
	int r;
	for (r = 0; r < ac->nrecv; r++) {
		if (join_path(path, sizeof(path), ac->workdir, ac->recv[r].file) != 0) {
			ac_err(ac, "output path too long for receiver #%d", r + 1);
			return 1;
		}
		if (ac_write_wav_f32(path, ac->fs,
		                     ac->rec + (size_t)r * ac->nsteps, ac->nsteps) != 0) {
			ac_err(ac, "cannot write %s", path);
			return 1;
		}
		ac_log(ac, "output: %s (%d samples, float32 %d Hz)",
		       ac->recv[r].file, ac->nsteps, ac->fs);
	}
	return 0;
}

int main(int argc, char **argv)
{
	static ac_t ac;               /* 状態はこの 1 個 (BSS — スタックに置かない) */
	char logpath[AC_PATH_MAX + 16];
	const char *args[2] = { NULL, NULL };
	int i, nargs = 0;

	memset(&ac, 0, sizeof(ac));
	for (i = 0; i < AC_NWALL; i++) ac.alpha[i] = AC_ALPHA_DEFAULT;

	/* 引数 : "-" で始まるものは無視 (ランチャ由来のオプションを許容) */
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			fprintf(stderr, "%s: warning: ignoring option '%s'\n",
			        AC_SOLVER_NAME, argv[i]);
			continue;
		}
		if (nargs < 2) args[nargs] = argv[i];
		nargs++;
	}
	if (nargs < 1 || nargs > 2) {
		fprintf(stderr, "usage: %s <working_dir> [<input_file.ofd>]\n",
		        AC_SOLVER_NAME);
		return 2;
	}
	if (strlen(args[0]) >= sizeof(ac.workdir)) {
		fprintf(stderr, "%s: error: working_dir path too long\n", AC_SOLVER_NAME);
		return 1;
	}
	strcpy(ac.workdir, args[0]);

	/* solver.log (契約の必須出力) を最初に開く — 書けない working_dir は
	 * その時点で異常。以降のログはファイルにも残る。 */
	if (join_path(logpath, sizeof(logpath), ac.workdir, "solver.log") != 0 ||
	    (ac.logfp = fopen(logpath, "w")) == NULL) {
		fprintf(stderr, "%s: error: cannot write solver.log in '%s' "
		        "(directory missing or not writable)\n",
		        AC_SOLVER_NAME, ac.workdir);
		return 1;
	}
	ac_log(&ac, "%s %s (contract_version %d)",
	       AC_SOLVER_NAME, AC_SOLVER_VERSION, AC_CONTRACT_VERSION);
	ac_log(&ac, "working dir: %s", ac.workdir);
	printf("%s %s\n", AC_SOLVER_NAME, AC_SOLVER_VERSION);
	fflush(stdout);

	/* 入力ファイルの決定 */
	if (nargs == 2) {
		/* 指定名 : まずそのまま、開けなければ working_dir 直下として解決 */
		FILE *fp = fopen(args[1], "rb");
		if (fp) {
			fclose(fp);
			if (strlen(args[1]) >= sizeof(ac.ofd_path)) {
				ac_err(&ac, "input path too long");
				goto failed;
			}
			strcpy(ac.ofd_path, args[1]);
		}
		else if (join_path(ac.ofd_path, sizeof(ac.ofd_path),
		                   ac.workdir, args[1]) != 0 ||
		         (fp = fopen(ac.ofd_path, "rb")) == NULL) {
			ac_err(&ac, "input file '%s' not found (neither as given nor in %s)",
			       args[1], ac.workdir);
			goto failed;
		}
		else {
			fclose(fp);
		}
	}
	else if (ac_find_input(&ac) != 0) {
		goto failed;   /* 理由は ac_find_input が stderr へ出力済み */
	}

	/* metadata 用ファイル名と .ofdx パス (同 basename の JSON サイドカー) */
	{
		const char *base = strrchr(ac.ofd_path, '/');
		const char *b2 = strrchr(ac.ofd_path, '\\');
		size_t len;
		if (b2 && (!base || b2 > base)) base = b2;
		base = base ? base + 1 : ac.ofd_path;
		snprintf(ac.ofd_name, sizeof(ac.ofd_name), "%.*s",
		         (int)sizeof(ac.ofd_name) - 1, base);
		len = strlen(ac.ofd_path);
		/* foo.ofd -> foo.ofdx (末尾が .ofd でなければ単に .ofdx を付加) */
		if (len >= 4 && strcmp(ac.ofd_path + len - 4, ".ofd") == 0)
			snprintf(ac.ofdx_path, sizeof(ac.ofdx_path), "%.*sofdx",
			         (int)(len - 3), ac.ofd_path);
		else
			snprintf(ac.ofdx_path, sizeof(ac.ofdx_path), "%s.ofdx", ac.ofd_path);
	}
	ac_log(&ac, "input: %s", ac.ofd_path);

	if (ac_read_ofd(&ac) != 0) goto failed;
	if (ac_read_ofdx(&ac) != 0) goto failed;
	if (ac_setup(&ac) != 0) goto failed;
	if (ac_run(&ac) != 0) goto failed;

	if (write_rirs(&ac) != 0) goto failed;
	if (ac_write_metadata(&ac) != 0) goto failed;
	if (ac_write_metrics(&ac) != 0) goto failed;

	ac_log(&ac, "normal end");
	printf("%s: done\n", AC_SOLVER_NAME);
	fclose(ac.logfp);
	ac.logfp = NULL;
	ac_free(&ac);
	return 0;

failed:
	if (ac.logfp) { fclose(ac.logfp); ac.logfp = NULL; }
	ac_free(&ac);
	return 1;
}
