/* jsonout.c — metadata.json / metrics.json の書き出し (ADR-0007 契約)
 *
 * metadata.json (必須) : contract_version、ソルバー名/バージョン、格子
 * (Δx・セル数)、fs、音源/受音点座標、音速、実行条件。
 * 互換規則は「追加キーのみ・未知キー無視」 — キーの改名・削除・型変更は禁止。
 *
 * metrics.json (任意) : GUI は自前計算 (RirAnalyzer) を正とするため、
 * v1 では突合用の情報のみ (指標の数値は出さない — 校正なしの絶対値や
 * 動的レンジ不足の T30 を「それらしく」出さない方針と同根)。
 */
#include <stdio.h>
#include <string.h>
#include "acoustic.h"

/* JSON 文字列エスケープ (必要最小限 : \ " と制御文字) */
static void jput_str(FILE *fp, const char *s)
{
	fputc('"', fp);
	for (; *s; s++) {
		unsigned char c = (unsigned char)*s;
		if (c == '"' || c == '\\') {
			fputc('\\', fp);
			fputc(c, fp);
		}
		else if (c < 0x20) {
			fprintf(fp, "\\u%04x", c);
		}
		else {
			fputc(c, fp);
		}
	}
	fputc('"', fp);
}

int ac_write_metadata(ac_t *ac)
{
	char path[AC_PATH_MAX + 16];
	FILE *fp;
	int r;

	if (snprintf(path, sizeof(path), "%s/metadata.json", ac->workdir)
	    >= (int)sizeof(path))
		return 1;
	fp = ac_fopen(path, "w");
	if (!fp) {
		ac_err(ac, "cannot write metadata.json in %s", ac->workdir);
		return 1;
	}
	fprintf(fp, "{\n");
	fprintf(fp, "  \"contract_version\": %d,\n", AC_CONTRACT_VERSION);
	fprintf(fp, "  \"solver\": \"%s\",\n", AC_SOLVER_NAME);
	fprintf(fp, "  \"solver_version\": \"%s\",\n", AC_SOLVER_VERSION);
	fprintf(fp, "  \"input_file\": ");
	jput_str(fp, ac->ofd_name);
	fprintf(fp, ",\n");
	fprintf(fp, "  \"grid\": {\n");
	fprintf(fp, "    \"dx_m\": %.10g,\n", ac->dx);
	fprintf(fp, "    \"cells\": [%d, %d, %d],\n", ac->nx, ac->ny, ac->nz);
	fprintf(fp, "    \"origin_m\": [%.10g, %.10g, %.10g],\n",
	        ac->x0, ac->y0, ac->z0);
	fprintf(fp, "    \"size_m\": [%.10g, %.10g, %.10g]\n",
	        ac->nx * ac->dx, ac->ny * ac->dx, ac->nz * ac->dx);
	fprintf(fp, "  },\n");
	fprintf(fp, "  \"sample_rate\": %d,\n", ac->fs);
	fprintf(fp, "  \"duration_s\": %.10g,\n", ac->duration);
	fprintf(fp, "  \"steps\": %d,\n", ac->nsteps);
	fprintf(fp, "  \"speed_of_sound_mps\": %.10g,\n", AC_C0);
	fprintf(fp, "  \"density_kgm3\": %.10g,\n", AC_RHO0);
	fprintf(fp, "  \"source\": {\n");
	fprintf(fp, "    \"pos_m\": [%.10g, %.10g, %.10g],\n",
	        ac->x0 + (ac->isrc + 0.5) * ac->dx,
	        ac->y0 + (ac->jsrc + 0.5) * ac->dx,
	        ac->z0 + (ac->ksrc + 0.5) * ac->dx);
	fprintf(fp, "    \"type\": \"gaussian_derivative_soft\",\n");
	fprintf(fp, "    \"fmax_hz\": %.10g,\n", ac->fmax);
	fprintf(fp, "    \"sigma_s\": %.10g,\n", ac->sigma);
	fprintf(fp, "    \"t0_s\": %.10g\n", ac->t0);
	fprintf(fp, "  },\n");
	/* 複数音源の契約 (ADR-0010) : source は feed #1 のまま (互換)、使用した
	 * 全音源はスナップ後のセル中心で sources に列挙する。 */
	fprintf(fp, "  \"multi_source\": %s,\n", ac->multi_source ? "true" : "false");
	fprintf(fp, "  \"sources\": [\n");
	for (r = 0; r < ac->nsrc; r++)
		fprintf(fp, "    { \"pos_m\": [%.10g, %.10g, %.10g] }%s\n",
		        ac->x0 + (ac->srccell[3 * r + 0] + 0.5) * ac->dx,
		        ac->y0 + (ac->srccell[3 * r + 1] + 0.5) * ac->dx,
		        ac->z0 + (ac->srccell[3 * r + 2] + 0.5) * ac->dx,
		        (r + 1 < ac->nsrc) ? "," : "");
	fprintf(fp, "  ],\n");
	fprintf(fp, "  \"receivers\": [\n");
	for (r = 0; r < ac->nrecv; r++) {
		const ac_recv_t *rv = &ac->recv[r];
		fprintf(fp, "    { \"name\": ");
		jput_str(fp, rv->name);
		fprintf(fp, ", \"pos_m\": [%.10g, %.10g, %.10g], \"file\": ",
		        ac->x0 + (rv->ic + 0.5) * ac->dx,
		        ac->y0 + (rv->jc + 0.5) * ac->dx,
		        ac->z0 + (rv->kc + 0.5) * ac->dx);
		jput_str(fp, rv->file);
		fprintf(fp, " }%s\n", (r + 1 < ac->nrecv) ? "," : "");
	}
	fprintf(fp, "  ],\n");
	fprintf(fp, "  \"boundary_alpha\": [%.10g, %.10g, %.10g, %.10g, %.10g, %.10g],\n",
	        ac->alpha[AC_XM], ac->alpha[AC_XP], ac->alpha[AC_YM],
	        ac->alpha[AC_YP], ac->alpha[AC_ZM], ac->alpha[AC_ZP]);
	fprintf(fp, "  \"t_sabine_s\": %.10g,\n", ac->tsab);   /* -1 = A=0 (無限大) */
	fprintf(fp, "  \"rigid_geometries\": %d,\n", ac->ngeom);
	fprintf(fp, "  \"ofdx_sidecar\": %s\n", ac->have_ofdx ? "true" : "false");
	fprintf(fp, "}\n");
	return (fclose(fp) == 0) ? 0 : 1;
}

int ac_write_metrics(ac_t *ac)
{
	char path[AC_PATH_MAX + 16];
	FILE *fp;

	if (snprintf(path, sizeof(path), "%s/metrics.json", ac->workdir)
	    >= (int)sizeof(path))
		return 1;
	fp = ac_fopen(path, "w");
	if (!fp) {
		ac_err(ac, "cannot write metrics.json in %s", ac->workdir);
		return 1;
	}
	fprintf(fp, "{\n");
	fprintf(fp, "  \"note\": \"acoustic metrics are computed by the GUI "
	            "(RirAnalyzer); this file is informational\",\n");
	fprintf(fp, "  \"solver\": \"%s %s\",\n", AC_SOLVER_NAME, AC_SOLVER_VERSION);
	if (ac->tsab > 0.0)
		fprintf(fp, "  \"t_sabine_s\": %.10g\n", ac->tsab);
	else
		fprintf(fp, "  \"t_sabine_s\": -1\n");
	fprintf(fp, "}\n");
	return (fclose(fp) == 0) ? 0 : 1;
}
