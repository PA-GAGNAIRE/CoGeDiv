'''
VariantCalling.snakefile
Pierre Barry
Map paired-end reads to reference genome and perform variant calling
-----------------------------------------------------------------------

Requirements:

Usage: 

snakemake \
	--snakefile VariantCalling.snakefile \

'''


### FUNCTION

def get_input_gvcf(species):
	SAMPLES_sp=[SAMPLES[i] for i in range(len(SAMPLES)) if SAMPLES[i][0:5]=="species"]
	L=[]
	for i in range(len(SAMPLES_sp)):
		L+=["/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/"+{species}+"/VariantCalling/"+SAMPLES_sp[i]+"_gcvf_first.g.vcf"]
	L=np.asmatrix(L)	
	with open("/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/"+species+"/VariantCalling/name.list", "wb") as f:
         for line in L:
             np.savetxt(f, line, fmt='%s', delimiter=";")
			 
def get_input_recal_gvcf(species):
	SAMPLES_sp=[SAMPLES[i] for i in range(len(SAMPLES)) if SAMPLES[i][0:5]=="species"]
	L=[]
	for i in range(len(SAMPLES_sp)):
		L+=["/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/"+{species}+"/VariantCalling/"+SAMPLES_sp[i]+"_gcvf_recal.g.vcf"]
	L=np.asmatrix(L)	
	with open("/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/"+species+"/VariantCalling/name.list", "wb") as f:
         for line in L:
             np.savetxt(f, line, fmt='%s', delimiter=";")
			 
###### BWA - MEM
rule index_reference:
	input:
		ref=".../ref_{species}.fa"
	output: 
		index_name="{species}_cogediv"		
	message:
		"Index reference : {wildcards.species}"
	shell:
		"bwa index "
		"-p {output.index_name} "
		"-a bwtsw "
		"{input.ref}"

rule reference_mapping:
	input:
		ref=".../ref_{species}.fa"
		fastp_R1="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/fastp/{sample}_R1_fastp.fastq.gz",
		fastp_R2="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/fastp/{sample}_R2_fastp.fastq.gz"
	output:
		align_sam="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_align.sam"
	message:
		"Reference mapping : {wildcards.sample}"
	shell : 
		"bwa mem "
		"-M "
		"-t 16 "
		"{input.ref} "
		"{input.fastp_R1} "
		"{input.fastp_R2} "
		"> {output.align_sam}"
		
rule samtools_convert_sam_bam:
	input:
		align_sam="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_align.sam"
	output:
		align_bam="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_align.bam"
	message:
		"Convert sam to bam : {wildcards.sample}"
	shell : 
		"samtools " 
		"view "
		"{input.align_sam} "
		"-bS "
		"-o {output.align_bam}"

rule samtools_sort:
	input:
		align_bam="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_align.bam"
	output:
		align_bam_sorted="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_align_sorted.bam"
	message:
		"Sort bam : {wildcards.sample}"
	shell : 
		"samtools " 
		"sort "
		"-bS "
		"{input.align_bam} "
		"-o {output.align_bam_sorted}"		

rule samtools_fixmate:
	input:
		align_bam_sorted="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_align_sorted.bam"
	output:
		touch("samtools_fixmate.done")
	message:
		"Fixmate bam : {wildcards.sample}"
	shell : 
		"samtools " 
		"fixmate "
		"-m "
		"{input.align_bam_sorted}"	

rule samtools_markduplicate:
	input:
		align_bam_sorted="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_align_sorted.bam"
	output:
		align_bam_sorted_dedup="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_align_sorted_dedup.bam"
	message:
		"Mark duplicates : {wildcards.sample}"
	shell : 
		"samtools " 
		"markdup "
		"-r " # remove duplicate reads
		"-S " # mark supplementary reads of duplicates as duplicates
		"-s " # print some summary stats
		"{input.align_bam_sorted} "
		"{output.align_bam_sorted_dedup}"	
		
rule samtools_stats:
	input:
		align_bam_sorted_dedup="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_align_sorted_dedup.bam"
	output:
		samtools_stats="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_samtools_stats.txt"
	message:
		"Samtools stats : {wildcards.sample}"
	shell : 
		"samtools " 
		"stats "
		"{input.align_bam_sorted_dedup} > "
		"{output.samtools_stats}"
		
rule plot_bamstats:
	input:
		samtools_stats="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_samtools_stats.txt"
	output:
		plot_samtools="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}_plotsamtools"
	message:
		"Plot bamstats : {wildcards.sample}"
	shell:
		"plot-bamstats "
		"-p "
		"{output.plot_samtools} "
		"{input.samtools_stats} "

### GATK
	
rule generate_gvcf_first:
	input:
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa",
		markdup="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}/{sample}_markdup.bam"
	output: 
		gvcf_first="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{sample}_gcvf_first.g.vcf"		
	message:
		"Generate GVCF first : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		-nct 16 "
		-T HaplotypeCaller "
		-R {input.reference_genome} "
		-I {input.markdup} "
		-o {output.gcvf_first} "
		--heterozygosity 0.005 "
		--genotyping_mode DISCOVERY "
		--emitRefConfidence GVCF "
		--variant_index_type LINEAR "
		--variant_index_parameter 128000"

		
rule joint_genotyping_first:
	input:
		variant_list=**get_input_gvcf({species}),
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		joint_gvcf_first="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_joint_gcvf_first.vcf"		
	message:
		"Joint genotyping first : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-nt 16 "
		"-T GenotypeGVCFs "
		"-R {input.reference_genome} "
		"-V {input.variant_list} "
		"--heterozygosity	0.005 "
		"-o {output.joint_gvcf_first}"
		
## Genotype refinement

# GenotypePosteriors

rule genotype_posteriors_first:
	input:
		joint_gvcf_first="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_joint_gcvf_first.vcf",
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		gvcf_genotype_posterior_first="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_genotype_posterior_first.vcf"		
	message:
		"Compute genotype posteriors first : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-T CalculateGenotypePosteriors "
		"-R {input.reference_genome} "
		"-V {input.joint_gvcf_first} "
		"-o {output.gvcf_genotype_posterior_first} "
		"--skipPopulationPriors "
		"--skipFamilyPriors"

# Variant filtering
##SNP

rule selectvariants_snp_first:
	input:
		gvcf_genotype_posterior_first="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_genotype_posterior_first.vcf"		
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		gvcf_refineSNP="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_refineSNP_first.vcf"		
	message:
		"Select Variant to SNP first : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-nt 16 "
		"-T SelectVariants "
		"-R {input.reference_genome} "
		"--variant {input.gvcf_genotype_posterior_first} "
		"-o {output.gvcf_refineSNP} "
		"-selectType SNP"

rule filtervariants_snp_first:
	input:
		gvcf_refineSNP="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_refineSNP_first.vcf"
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		gvcf_refineSNP_HQ="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_refineSNP_HQ_first.vcf"		
	message:
		"Variant Filtration to SNP first : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-T VariantFiltration "
		"-R {input.reference_genome} "
		"--variant {input.gvcf_refineSNP} "
		"-o {output.gvcf_refineSNP_HQ} "
		"--filterExpression "QD < 10.0 || FS > 60.0 || MQ < 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" "
		"--filterName "lowQ" "
		"--genotypeFilterExpression "GQ < 30.0" " 
		"--genotypeFilterName "lowGQ""

##Indel

rule selectvariants_indel_first:
	input:
		gvcf_genotype_posterior_first="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_genotype_posterior_first.vcf"		
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		gvcf_refineINDEL="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_refineINDEL_first.vcf"		
	message:
		"Select Variant to INDEL first : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-nt 16 "
		"-T SelectVariants "
		"-R {input.reference_genome} "
		"--variant {input.gvcf_genotype_posterior_first} "
		"-o {output.gvcf_refineINDEL} "
		"-selectType INDEL"

rule filtervariants_indel_first:
	input:
		gvcf_refineINDEL="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_refineINDEL_first.vcf"
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		gvcf_refineINDEL_HQ="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_refineINDEL_HQ_first.vcf"		
	message:
		"Variant Filtration to INDEL first : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-T VariantFiltration "
		"-R {input.reference_genome} "
		"--variant {input.gvcf_refineINDEL} "
		"-o {output.gvcf_refineINDEL_HQ} "
		"--filterExpression "QD < 10.0 || FS > 200.0 || ReadPosRankSum < -20.0" "
		"--filterName "lowQ" "
		"--genotypeFilterExpression "GQ < 30.0" " 
		"--genotypeFilterName "lowGQ""
		
## Combine variants

rule combine_variants_first:
	input:
		gvcf_refineSNP_HQ="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_refineSNP_HQ_first.vcf"	
		gvcf_refineINDEL_HQ="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_refineINDEL_HQ_first.vcf"		
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		gvcf_refine_HQ="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_refine_HQ_first.vcf"		
	message:
		"Combine variants first : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-nt 16 "
		"-R {input.reference_genome} "
		"-T CombineVariants "
		"-o {output.gvcf_refine_HQ} "
		"-V:snp {input.gvcf_refineSNP_HQ} "
		"-V:indel {input.gvcf_refineINDEL_HQ} "
		"--filteredrecordsmergetype KEEP_IF_ALL_UNFILTERED "
		"--genotypemergeoption PRIORITIZE "
		"-priority snp,indel"	

## Finalize

rule database_vcf_first:
	input:
		gvcf_refine_HQ="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_refine_HQ_first.vcf"		
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		dbase_first="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_dbase_first.vcf"		
	message:
		"Finalize variants first : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-nt 16 "
		"-R {input.reference_genome} "
		"-T SelectVariants "
		"--variant {input.gvcf_refine_HQ} "
		"-o {output.dbase_first} "
		"--excludeFiltered"	
		
## BASE RECALIBRATION

# Recalibrate quality scores

# first round
rule base_recalibration_firstround:
	input:
		markdup="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}/{sample}_markdup.bam",
		dbase_first="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_dbase_first.vcf",		
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
	    recal_table="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_recal.table"
	message:
		"Base recalibration - first round : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-nct 16 "
		"-R {input.reference_genome} "
		"-T BaseRecalibrator "
		"-I {input.markdup} "
		"-knownSites {input.dbase_first} "
		"-o {output.recal_table}"
		
# second round
rule base_recalibration_secondround:
	input:
		markdup="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}/{sample}_markdup.bam",
		dbase_first="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_dbase_first.vcf",	
	    recal_table="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_recal.table",		
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
	    recal_post_table="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_recal_post.table"
	message:
		"Base recalibration - second round : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-nct 16 "
		"-R {input.reference_genome} "
		"-T BaseRecalibrator "
		"-I {input.markdup} "
		"-knownSites {input.dbase_first} "
		"-BQSR {input.recal_table} "
		"-o {output.recal_post_table}"
		
# plot
rule plot_BQSR:
	input:
		markdup="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}/{sample}_markdup.bam",
	    recal_table="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_recal.table",	
		recal_post_table="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_recal_post.table",
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
	    recal_plots="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_recal_plots.pdf"
	message:
		"Plot base recalibration : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-R {input.reference_genome} "
		"-T AnalyzeCovariates "
		"-before {input.recal_table} "
		"-after {input.recal_post_table} "
		"-plots {output.recal_plots}"				
		
# print reads		
rule recalibrate_bam:
	input:
		markdup="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Mapping/{sample}/{sample}_markdup.bam",
	    recal_table="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_recal.table",	
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
	    recal_bam="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_markdup_recal.bam"
	message:
		"Regenerate recalibrate bam : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-nct 16 "
		"-R {input.reference_genome} "
		"-T PrintReads "
		"-I {input.markdup} "
		"-BQSR {input.recal_table} "
		"-o {output.recal_bam}" 				

# Haplotype caller after recalibration

rule generate_gvcf_second:
	input:
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa",
		recal_bam="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_markdup_recal.bam"
	output: 
		gvcf_recal="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{sample}_gcvf_recal.g.vcf"		
	message:
		"Generate GVCF after recalibration : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		-nct 16 "
		-T HaplotypeCaller "
		-R {input.reference_genome} "
		-I {input.recal_bam} "
		-o {output.gvcf_recal} "
		--heterozygosity 0.005 "
		--genotyping_mode DISCOVERY "
		--emitRefConfidence GVCF "
		--variant_index_type LINEAR "
		--variant_index_parameter 128000"		
		
rule joint_genotyping_recal:
	input:
		variant_list=**get_input_recal_gvcf({species}),
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		joint_gvcf_recal="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_joint_gcvf_recal.vcf"		
	message:
		"Joint genotyping after recalibration : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-nt 16 "
		"-T GenotypeGVCFs "
		"-R {input.reference_genome} "
		"-V {input.variant_list} "
		"--heterozygosity	0.005 "
		"-o {output.joint_gvcf_recal}"		
		
## VQSR

rule vqsr:
	input:
		joint_gvcf_recal="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_joint_gcvf_recal.vcf",
		dbase_first="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_dbase_first.vcf",
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		joint_gvcf_recal="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_joint_gcvf_recal.vcf"		
	message:
		"VQSR : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-nt 16 "
		"-T VariantRecalibrator "
		"-R {input.reference_genome} "
		"-input {input.joint_gvcf_recal} "
		"-resource:DB_HQ,known=false,training=true,truth=true,prior=10.0 {input.dbase_first} "
		"-an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum "
		"-mode SNP "
		"-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "
		"-recalFile ${infolder}/joint_bwa_mem_mdup_IR_recal_SNP.recal " <- ? 
		"-tranchesFile ${infolder}/joint_bwa_mem_mdup_IR_recal_SNP.tranches " <- ?
		"-rscriptFile ${infolder}/joint_bwa_mem_mdup_IR_recal_SNP_plots.R" 
		
## Final Genotype Refinement

# GenotypePosteriors

rule genotype_posteriors_recal:
	input:
		joint_gvcf_recal="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_joint_gcvf_recal.vcf",
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		gvcf_genotype_posterior_recal="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_genotype_posterior_recal.vcf"		
	message:
		"Compute genotype posteriors after recalibration : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-T CalculateGenotypePosteriors "
		"-R {input.reference_genome} "
		"-V {input.joint_gvcf_recal} "
		"-o {output.gvcf_genotype_posterior_recal} "
		"--skipPopulationPriors "
		"--skipFamilyPriors"		

rule filtervariants_snp_recal:
	input:
		gvcf_genotype_posterior_recal="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_genotype_posterior_recal.vcf"
		reference_genome="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/Reference_Genome/referencegenome_{species}.fa"
	output: 
		gvcf_refineVCF_HQ="/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/{species}/VariantCalling/{species}_gvcf_refineVCF_HQ_recal.vcf"		
	message:
		"Variant Filtration to SNP after recalibration : {wildcards.sample}"
	shell:
		"java -Xmx4g -jar /usr/share/java/GenomeAnalysisTK.jar "
		"-T VariantFiltration "
		"-R {input.reference_genome} "
		"--variant {input.gvcf_genotype_posterior_recal} "
		"-o {output.gvcf_refineVCF_HQ} "
		"--genotypeFilterExpression "GQ < 30.0" " 
		"--genotypeFilterName "lowGQ""		
		
		
		
		
		
		
		
		

