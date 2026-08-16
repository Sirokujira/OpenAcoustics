/* ga_json.c — metadata.json / metrics.json の書き出し (ADR-0007 契約)
 *
 * キーの並びと意味は ofdx_acoustic_fdtd の src/jsonout.c を基準にし、
 * 幾何音響で意味を持たないものは「該当なし」を表す値にしてある
 * (キーの削除・改名はしない — 互換規則は「追加のみ」)。
 *   grid : 幾何音響に格子は無いので dx_m = 0 / cells = [0,0,0]。
 *          origin_m / size_m は室の直方体で意味を持つ。gridless = true を併記。
 *   source.type / sigma_s / t0_s : 音源は理想インパルス。**t0_s = 0** が
 *          「t = 0 は音源発火時刻」という本ソルバーの時間原点の規約。
 *   boundary_alpha : 従来キーは互換のため帯域平均を出し、バンド別の実値は
 *          追加キー boundary_alpha_bands に出す。
 * 追加キー : valid_band_hz / bands_hz / boundary_alpha_bands / image_order /
 *          rays / air / t_eyring_s / amplitude_convention / time_origin。
 */
#include <stdio.h>
#include <string.h>
#include "ga.h"

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

static void jput_darray(FILE *fp, const double *v, int n)
{
	int i;
	fputc('[', fp);
	for (i = 0; i < n; i++)
		fprintf(fp, "%s%.10g", i ? ", " : "", v[i]);
	fputc(']', fp);
}

int ga_write_metadata(ga_t *g)
{
	char path[GA_PATH_MAX + 16];
	FILE *fp;
	double fc[GA_NBAND], amean[GA_NWALL];
	int r, w, b;

	for (b = 0; b < GA_NBAND; b++) fc[b] = GA_BAND_F0 * (double)(1 << b);
	for (w = 0; w < GA_NWALL; w++) {
		double s = 0.0;
		for (b = 0; b < GA_NBAND; b++) s += g->alpha[w][b];
		amean[w] = s / GA_NBAND;
	}

	if (snprintf(path, sizeof(path), "%s/metadata.json", g->workdir)
	    >= (int)sizeof(path))
		return 1;
	fp = ac_fopen(path, "w");
	if (!fp) {
		ga_err(g, "cannot write metadata.json in %s", g->workdir);
		return 1;
	}
	fprintf(fp, "{\n");
	fprintf(fp, "  \"contract_version\": %d,\n", GA_CONTRACT_VERSION);
	fprintf(fp, "  \"solver\": \"%s\",\n", GA_SOLVER_NAME);
	fprintf(fp, "  \"solver_version\": \"%s\",\n", GA_SOLVER_VERSION);
	fprintf(fp, "  \"method\": \"geometric (image source + ray tracing)\",\n");
	fprintf(fp, "  \"input_file\": ");
	jput_str(fp, g->ofd_name);
	fprintf(fp, ",\n");
	fprintf(fp, "  \"gridless\": true,\n");
	fprintf(fp, "  \"grid\": {\n");
	fprintf(fp, "    \"dx_m\": 0,\n");
	fprintf(fp, "    \"cells\": [0, 0, 0],\n");
	fprintf(fp, "    \"origin_m\": [%.10g, %.10g, %.10g],\n",
	        g->x0, g->y0, g->z0);
	fprintf(fp, "    \"size_m\": [%.10g, %.10g, %.10g]\n",
	        g->x1 - g->x0, g->y1 - g->y0, g->z1 - g->z0);
	fprintf(fp, "  },\n");
	fprintf(fp, "  \"room\": { \"volume_m3\": %.10g, \"surface_m2\": %.10g },\n",
	        g->vol, g->area);
	fprintf(fp, "  \"sample_rate\": %d,\n", GA_FS);
	fprintf(fp, "  \"duration_s\": %.10g,\n", g->duration);
	fprintf(fp, "  \"steps\": %d,\n", g->nsamples);
	fprintf(fp, "  \"speed_of_sound_mps\": %.10g,\n", GA_C0);
	fprintf(fp, "  \"density_kgm3\": %.10g,\n", GA_RHO0);
	fprintf(fp, "  \"valid_band_hz\": [%.10g, %.10g],\n", g->flo, g->fhi);
	fprintf(fp, "  \"bands_hz\": ");
	jput_darray(fp, fc, GA_NBAND);
	fprintf(fp, ",\n");
	fprintf(fp, "  \"amplitude_convention\": \"free-field direct sound = "
	            "1/(4*pi*r); the sample sum of an arrival is its amplitude\",\n");
	fprintf(fp, "  \"time_origin\": \"source firing (direct sound at t = r/c)\",\n");
	/* 複数音源の契約 (ADR-0010) : source は feed #1 のまま (互換)、使用した
	 * 全音源は sources に列挙する。multi_source が false なら 1 個。
	 * gain / delay_s (Decision 7) はキー追加のみ (既定 1 / 0)。 */
	fprintf(fp, "  \"multi_source\": %s,\n", g->multi_source ? "true" : "false");
	fprintf(fp, "  \"sources\": [\n");
	{
		int n2;
		for (n2 = 0; n2 < g->nsrc; n2++)
			fprintf(fp, "    { \"pos_m\": [%.10g, %.10g, %.10g], "
			        "\"gain\": %.10g, \"delay_s\": %.10g }%s\n",
			        g->feedpos[3 * n2 + 0], g->feedpos[3 * n2 + 1],
			        g->feedpos[3 * n2 + 2], g->srcgain[n2], g->srcdelay[n2],
			        (n2 + 1 < g->nsrc) ? "," : "");
	}
	fprintf(fp, "  ],\n");
	fprintf(fp, "  \"source\": {\n");
	fprintf(fp, "    \"pos_m\": [%.10g, %.10g, %.10g],\n",
	        g->srcx, g->srcy, g->srcz);
	fprintf(fp, "    \"type\": \"ideal_impulse_omni\",\n");
	fprintf(fp, "    \"fmax_hz\": %.10g,\n", g->fhi);
	fprintf(fp, "    \"sigma_s\": 0,\n");
	fprintf(fp, "    \"t0_s\": 0\n");
	fprintf(fp, "  },\n");
	fprintf(fp, "  \"receivers\": [\n");
	for (r = 0; r < g->nrecv; r++) {
		const ga_recv_t *rv = &g->recv[r];
		fprintf(fp, "    { \"name\": ");
		jput_str(fp, rv->name);
		fprintf(fp, ", \"pos_m\": [%.10g, %.10g, %.10g], \"file\": ",
		        rv->x, rv->y, rv->z);
		jput_str(fp, rv->file);
		fprintf(fp, ", \"sphere_radius_m\": %.10g, \"image_sources\": %d, "
		        "\"image_sources_blocked\": %d, \"ray_detections\": %ld }%s\n",
		        rv->radius, rv->nimage, rv->nblocked, rv->nhit,
		        (r + 1 < g->nrecv) ? "," : "");
	}
	fprintf(fp, "  ],\n");
	fprintf(fp, "  \"boundary_alpha\": ");
	jput_darray(fp, amean, GA_NWALL);
	fprintf(fp, ",\n");
	fprintf(fp, "  \"boundary_alpha_bands\": [\n");
	for (w = 0; w < GA_NWALL; w++) {
		fprintf(fp, "    ");
		jput_darray(fp, g->alpha[w], GA_NBAND);
		fprintf(fp, "%s\n", (w + 1 < GA_NWALL) ? "," : "");
	}
	fprintf(fp, "  ],\n");
	fprintf(fp, "  \"image_order\": %d,\n", g->order);
	fprintf(fp, "  \"rays\": %d,\n", g->nrays);
	fprintf(fp, "  \"scattering\": %.10g,\n", g->scatter);
	fprintf(fp, "  \"scattering_walls\": ");
	{
		double sw[GA_NWALL];
		for (w = 0; w < GA_NWALL; w++) sw[w] = g->surf[w].scatter;
		jput_darray(fp, sw, GA_NWALL);
	}
	fprintf(fp, ",\n");
	fprintf(fp, "  \"angle_dependent_absorption\": %s,\n",
	        g->angle_dep ? "true" : "false");
	fprintf(fp, "  \"reflecting_surfaces\": %d,\n", g->nsurf);
	fprintf(fp, "  \"obstacle_materials\": [");
	{
		int first = 1, n;
		for (n = 0; n < g->ngeom; n++) {
			const ga_geom_t *o = &g->geom[n];
			if (!o->ok || !o->has_mat) continue;
			fprintf(fp, "%s\n    { \"geometry\": %d, \"alpha\": ",
			        first ? "" : ",", n + 1);
			jput_darray(fp, o->mat_alpha, GA_NBAND);
			fprintf(fp, ", \"scattering\": %.10g }",
			        (o->mat_scatter >= 0.0) ? o->mat_scatter : g->scatter);
			first = 0;
		}
		fprintf(fp, "%s],\n", first ? "" : "\n  ");
	}
	fprintf(fp, "  \"air\": { \"enabled\": %s, \"temperature_c\": %.10g, "
	            "\"humidity_percent\": %.10g, \"pressure_kpa\": %.10g, "
	            "\"attenuation_db_per_m\": ",
	        g->air_on ? "true" : "false", g->temp_c, g->humid, g->press_kpa);
	jput_darray(fp, g->air_db_m, GA_NBAND);
	fprintf(fp, " },\n");
	fprintf(fp, "  \"t_eyring_s\": ");
	jput_darray(fp, g->teyring, GA_NBAND);
	fprintf(fp, ",\n");
	fprintf(fp, "  \"t_sabine_s\": %.10g,\n", g->tsabmax);  /* -1 = 無損失 */
	fprintf(fp, "  \"t_sabine_bands_s\": ");
	jput_darray(fp, g->tsabine, GA_NBAND);
	fprintf(fp, ",\n");
	fprintf(fp, "  \"t_eyring_max_s\": %.10g,\n", g->t60max);
	fprintf(fp, "  \"rigid_geometries\": %d,\n", g->ngeom);
	fprintf(fp, "  \"ofdx_sidecar\": %s\n", g->have_ofdx ? "true" : "false");
	fprintf(fp, "}\n");
	return (fclose(fp) == 0) ? 0 : 1;
}

int ga_write_metrics(ga_t *g)
{
	char path[GA_PATH_MAX + 16];
	FILE *fp;

	if (snprintf(path, sizeof(path), "%s/metrics.json", g->workdir)
	    >= (int)sizeof(path))
		return 1;
	fp = ac_fopen(path, "w");
	if (!fp) {
		ga_err(g, "cannot write metrics.json in %s", g->workdir);
		return 1;
	}
	fprintf(fp, "{\n");
	fprintf(fp, "  \"note\": \"acoustic metrics are computed by the GUI "
	            "(RirAnalyzer); this file is informational\",\n");
	fprintf(fp, "  \"solver\": \"%s %s\",\n", GA_SOLVER_NAME, GA_SOLVER_VERSION);
	fprintf(fp, "  \"t_eyring_s\": ");
	jput_darray(fp, g->teyring, GA_NBAND);
	fprintf(fp, ",\n");
	fprintf(fp, "  \"t_sabine_s\": %.10g,\n", g->tsabmax);
	fprintf(fp, "  \"t_eyring_max_s\": %.10g,\n", g->t60max);
	fprintf(fp, "  \"valid_band_hz\": [%.10g, %.10g]\n", g->flo, g->fhi);
	fprintf(fp, "}\n");
	return (fclose(fp) == 0) ? 0 : 1;
}
