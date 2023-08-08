#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


params.gtf="/home/ezotova/230227_regulotyping/data/rnaseq/eqtl/combined.annotation.collapsed.gtf"
params.chr_list="/home/ezotova/230227_regulotyping/data/rnaseq/eqtl/vcf_chr_list"
params.genotype_covs="/home/ezotova/230227_regulotyping/data/rnaseq/eqtl/chroms1-22.phaseI+II+stim.annotated.ancestral.top5_PCA_covs.tsv"
params.other_covs="/home/ezotova/230227_regulotyping/data/rnaseq/eqtl/covs.tsv"
params.genotyping_panel="/net/seq/data/projects/regulotyping-phaseI+II+stim/genotyping_panel/chroms1-22.phaseI+II+stim.annotated.ancestral.vcf.gz"


process prepare_expression {
        container "/home/ezotova/230227_regulotyping/gtex-pipeline/qtl/apptainer_sandbox"

	input: tuple val(sample_id), path(tpm_gct), path(raw_counts_gct), path(participant_lookup)

	output: tuple val("${sample_id}"), path("${sample_id}.expression.bed.gz") 
	
	script:
	"""
	/home/ezotova/230227_regulotyping/gtex-pipeline/qtl/src/eqtl_prepare_expression.py \
	${tpm_gct} \
	${raw_counts_gct} \
	${params.gtf} \
	${participant_lookup} \
	${sample_id} \
	--chr ${params.chr_list} \
	--tpm_threshold 0.1 \
	--count_threshold 6 \
	--sample_frac_threshold 0.2 \
	--normalization_method tmm
	"""
}

process estimate_hidden_covs {
	container "docker://broadinstitute/gtex_eqtl:V8"

	input: tuple val(sample_id), path(expression_bed)

	output: tuple val("${sample_id}"), path("${sample_id}.PEER_covariates.txt")

	script:
	"""
	Rscript /home/ezotova/230227_regulotyping/gtex-pipeline/qtl/src/run_PEER.R ${expression_bed} ${sample_id} 15
	"""
}

process merge_covs {
	container "/home/ezotova/230227_regulotyping/gtex-pipeline/qtl/apptainer_sandbox"

	input: tuple val(sample_id), path(peer_covs)

	output: tuple val("${sample_id}"), path("${sample_id}.combined_covariates.txt")

	script:
	"""
	python /home/ezotova/230227_regulotyping/gtex-pipeline/qtl/src/combine_covariates.py \
	${peer_covs} \
	${sample_id} \
	    --genotype_pcs ${params.genotype_covs} \
	    --add_covariates ${params.other_covs}
	"""
}

process run_eqtl {
	container "docker://broadinstitute/gtex_eqtl:V8"
	cpus 36

	input: tuple path(sample_id), path(expression_bed), path(covs)
		
	script:
	"""
	/opt/fastqtl/python/run_FastQTL_threaded.py \
		${params.genotyping_panel} \
		${expression_bed} \
		${sample_id} \
		--covariates ${covs} \
		--window 1e6 --chunks 100 --threads ${cpus}


	/opt/fastqtl/python/run_FastQTL_threaded.py \
		${params.genotyping_panel} \
		${expression_bed} \
		${sample_id} \
		--covariates ${covs} \
		--window 1e6 --chunks 100 --threads ${cpus} \
		--permute 1000 10000
	"""
}

workflow {
	expression = Channel.of( "phaseI+II+stim_h.CD19+", "/home/ezotova/230227_regulotyping/data/rnaseq/eqtl/inputs/phaseI+II+stim_tpm_h.CD19+.gct", "/home/ezotova/230227_regulotyping/data/rnaseq/eqtl/inputs/phaseI+II+stim_raw_counts_h.CD19+.gct", "/home/ezotova/230227_regulotyping/data/rnaseq/eqtl/inputs/phaseI+II+stim_participant_lookup_h.CD19+.tsv" ).collate(4).view()
	| prepare_expression
	covs = estimate_hidden_covs(expression)
	| merge_covs
	expression.join(covs).view()
//	| run_eqtl


}
