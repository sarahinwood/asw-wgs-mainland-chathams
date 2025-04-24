#!/usr/bin/env python3
import peppy

##ASW CHATHAMS VS NZ MAINLAND

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

## req. for admixture
def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

#containers
bbduk_container = 'https://github.com/deardenlab/container-bbmap/releases/download/0.0.3/container-bbmap.bbmap_38.90.sif'
fastqc_container = 'docker://biocontainers/fastqc:v0.11.9_cv8'
multiqc_container = 'docker://ewels/multiqc:1.9'
bwa_container = 'docker://staphb/bwa:0.7.17'
gatk_container = 'docker://broadinstitute/gatk:4.3.0.0'
bcftools_container = 'docker://staphb/bcftools:1.16'
plink2_container = 'docker://pgscatalog/plink2:2.00a5.10'
glnexus_container = 'docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1'
admixture_container = 'docker://evolbioinfo/admixture:v1.3.0'
vcftools_container = 'docker://biocontainers/vcftools:v0.1.16-1-deb_cv1'


#########
# RULES #
#########

rule target:
    input:
        ## trim, fastQC --> multiQC
        expand('output/01_mapping/multiqc/multiqc_report.html'),
        ## mapped --> indexed bams & QC
        expand('output/01_mapping/bwa_GCA_014170235_1/{sample}_sorted.bam.bai', sample=all_samples),
        ## mapping QC
        'output/01_mapping/mapping_qc_analysis/mapping_QC_res.csv',
        expand('output/01_mapping/samtools_GCA_014170235_1/coverage/{sample}_coverage.out', sample=all_samples),
        expand('output/01_mapping/samtools_GCA_014170235_1/depth/{sample}_depth.out', sample=all_samples),
        ## generate gvcf
        'output/02_variants/ASW_agresearch_combined.gvcf.bcf',
        ## filtering
        expand('output/03_filtering/vcf_stats/ASW_agresearch_{filtering}_stats.txt', filtering=["combined", "combined_filtered"]),
        ## plink - fst, pca, admixture
        'output/04_plink/fst/population_fst_results.fst.summary',
        'output/04_plink/no_ldpruning/filtered_snps_plink_pca.eigenvec',
        'output/04_plink/ld_pruned/admixture/admixture_cv_res.out',
        ## mashtree
        'output/05_mashtree/asw_mashtree_rooted_bootstrap.dnd'
        ## other pop stats
        'output/06_pop_stats/heterozygosity_rates.het',
        expand('output/06_pop_stats/nucleotide_diversity_{population}_100kb.windowed.pi', population=["Chatham", "Fortrose", "Lincoln", "Stewart"]),
        expand('output/06_pop_stats/tajimas_d_{population}_100kb.Tajima.D', population=["Chatham", "Fortrose", "Lincoln", "Stewart"])


#####################
## other pop stats ##
#####################

rule tajimas_d:
    input:
        vcf = 'output/03_filtering/ASW_agresearch_combined_filtered.vcf.gz',
        pops = 'data/{population}_pop.txt'
    output:
        'output/06_pop_stats/tajimas_d_{population}_100kb.Tajima.D'
    params:
        wd = 'output/06_pop_stats/tajimas_d_{population}_100kb'
    log:
        'output/logs/tajimas_d_{population}.log'
    singularity:
        vcftools_container
    shell:
        'vcftools '
        '--gzvcf {input.vcf} '
        '--keep {input.pops} '
        '--TajimaD 100000 '
        '--out {params.wd} '
        '&> {log}'



rule nucleotide_diversity:
    input:
        vcf = 'output/03_filtering/ASW_agresearch_combined_filtered.vcf.gz',
        pops = 'data/{population}_pop.txt'
    output:
        'output/06_pop_stats/nucleotide_diversity_{population}_100kb.windowed.pi'
    params:
        wd = 'output/06_pop_stats/nucleotide_diversity_{population}_100kb'
    log:
        'output/logs/nucleotide_diversity_{population}.log'
    singularity:
        vcftools_container
    shell:
        'vcftools '
        '--gzvcf {input.vcf} '
        '--keep {input.pops} '
        '--window-pi 100000 '
        '--out {params.wd} '
        '&> {log}'



rule heterozygosity_rates:
    input:
        'output/03_filtering/ASW_agresearch_combined_filtered.vcf.gz'
    output:
        'output/06_pop_stats/heterozygosity_rates.het'
    params:
        wd = 'output/06_pop_stats/heterozygosity_rates'
    log:
        'output/logs/heterozygosity_rates.log'
    singularity:
        vcftools_container
    shell:
        'vcftools --gzvcf {input} --het --out {params.wd} 2> {log}'



############
# mashtree #
############


rule mashtree_boot:
    input:
        samples = expand('output/05_mashtree/consensus/{sample}.fasta', sample=all_samples)
    output:
        'output/05_mashtree/asw_mashtree_rooted_bootstrap.dnd'
    singularity:
        mashtree_container
    log:
        'output/logs/mashtree.log'
    shell:
        'nice mashtree_bootstrap.pl '
        '{input.samples} '
        '--reps 100 '
        '--numcpus 12 '
        '-- --min-depth 0 '
        '2> {log} '
        '--outtree {output} '
        '> mashtree_execution.err'



# convert vcf to fastas=
rule bcf_consensus:
    input:
        vcf = '/output/03_filtering/ASW_agresearch_combined_filtered.vcf.gz',
        asw_ref = 'data/GCA_014170235_1.fa'
    output:
        'output/05_mashtree/consensus/{sample}.fasta'
    log:
        'output/logs/bcf_consensus/{sample}.log'
    params:
        sample = "{sample}"
    singularity:
        bcftools_container
    shell:
        'nice bcftools consensus '
        '--haplotype I '
        '--fasta-ref {input.asw_ref} '
        '--output {output} '
        '-s {params.sample} '
        '2> {log} '
        '{input.vcf}'



rule bcftools_index:
    input:
        'outout/03_filtering/ASW_agresearch_combined_filtered.vcf.gz'
    output:
        'output/03_filtering/ASW_agresearch_combined_filtered.vcf.gz.tbi'
    singularity:
       bcftools_container
    shell:
        'bcftools index {input} -t'




###########################
## ADMIXTURE on WGS only ##
###########################



# determine best K value from CV
rule grep_admixture_cvs_WGS_only:
    input:
        expand('output/04_plink/ld_pruned/admixture/admixture.{k}.log', k=["2", "3", "4", "5"])
    output:
        'output/04_plink/ld_pruned/admixture/admixture_cv_res.out'
    shell:
        'grep CV {input} > {output}'



# run admixture on pruned bed file
rule admixture_WGS_only:
    input:
        bed = 'output/04_plink/ld_pruned/admixture/admixture.bed',
        bim = 'output/04_plink/ld_pruned/admixture/admixture.bim',
        fam =  'output/04_plink/ld_pruned/admixture/admixture.fam'
    output:
        q = 'output/04_plink/ld_pruned/admixture/admixture.{k}.Q',
        log = 'output/04_plink/ld_pruned/admixture/admixture.{k}.log'
    params:
        wd = 'output/04_plink/ld_pruned/admixture/',
        bed = lambda wildcards, input: resolve_path(input.bed),
        k = '{k}'
    log:
        str(pathlib2.Path(resolve_path('output/04_plink/ld_pruned/admixture/'),
                            'admixture.{k}.log'))
    threads:
        10
    singularity:
        admixture_container
    shell:
        'cd {params.wd} || exit 1 ;'
        'admixture '
        '--cv ' ## outputs cross-validation scores (lowest = best fit) - cv=5 is default
        '--seed=43 ' ## set seed for reproducibility
        '-j{threads} '
        '{params.bed} '
        '{params.k} ' # no. populations
        '> {log}'



## copy fam to have renamed in filename also
rule copy_fam_WGS_only:
    input:
        'output/04_plink/ld_pruned/ld_pruned_snps_plink_pca.fam'
    output:
        'output/04_plink/ld_pruned/admixture/admixture.fam'
    shell:
        'cp {input} {output}'



## copy bed to have renamed in filename also
rule copy_bed_WGS_only:
    input:
        'output/04_plink/ld_pruned/ld_pruned_snps_plink_pca.bed'
    output:
        'output/04_plink/ld_pruned/admixture/admixture.bed'
    shell:
        'cp {input} {output}'



## convert Chr names in bim to only numbers
rule rename_bim_for_admixture2_WGS_only:
    input:
        'output/04_plink/ld_pruned/admixture/admixture_renamed1.bim'
    output:
        'output/04_plink/ld_pruned/admixture/admixture.bim'
    shell:
        "sed 's/\.//g' {input}>{output} "



## convert Chr names in bim to only numbers
rule rename_bim_for_admixture1_WGS_only:
    input:
        'output/04_plink/ld_pruned/ld_pruned_snps_plink_pca.bim'
    output:
        temp('output/04_plink/ld_pruned/admixture/admixture_renamed1.bim')
    shell:
        "sed 's/JACEGR*//g' {input}>{output} "



################
## plink pcas ##
################

# pca on pruned
rule plink2_pca_ld:
    input:
        vcf = 'output/03_filtering/ASW_agresearch_combined_filtered.vcf.gz',
        pruned = 'output/04_plink/ld_pruned/ld_pruned_snps.prune.in',
        freq = 'output/04_plink/no_ldpruning/filtered_snps_freq.afreq'
    output:
        eigenvec = 'output/04_plink/ld_pruned/ld_pruned_snps_plink_pca.eigenvec',
        bim = 'output/04_plink/ld_pruned/ld_pruned_snps_plink_pca.bim',
        bed = 'output/04_plink/ld_pruned/ld_pruned_snps_plink_pca.bed',
        fam = 'output/04_plink/ld_pruned/ld_pruned_snps_plink_pca.fam'
    params:
        out = 'output/04_plink/ld_pruned/ld_pruned_snps_plink_pca',
    log:
        'logs/plink_pca_LD.log'
    singularity:
        plink2_container
    shell:
        'plink2 '
        '--vcf {input.vcf} '
        '--read-freq {input.freq} '
        '--set-all-var-ids @:# '
        '--allow-extra-chr '
        '--vcf-half-call missing '
        '--extract {input.pruned} '
        '--pca 20 '
        '--make-bed '
        '--out {params.out} '
        '&> {log}'



# prune dataset of variants in linkage - PCA relies on independent variables
rule plink2_ld_prune:
    input:
        'output/03_filtering/ASW_agresearch_combined_filtered.vcf.gz'
    output:
        'output/04_plink/ld_pruned/ld_pruned_snps.prune.in'
    params:
        indep = '10 10 0.1',     # 10 kb window, 10 bp step size, phased r2 threshold < 0.1
        out = 'output/04_plink/ld_pruned/ld_pruned_snps'
    log:
        'output/logs/plink_prune_linkage.log'
    singularity:
        plink2_container
    shell:
        'plink2 '
        '--vcf {input} '
        '--set-all-var-ids @:# '
        '--allow-extra-chr '
        '--vcf-half-call missing '
        '--bad-ld '
        '--indep-pairwise {params.indep} '
        '--out {params.out} '
        '&> {log}'



rule plink2_pca_no_ldpruning:
    input:
        vcf = 'output/03_filtering/ASW_agresearch_combined_filtered.vcf.gz',
        freq = 'output/04_plink/no_ldpruning/filtered_snps_freq.afreq'
    output:
        'output/04_plink/no_ldpruning/filtered_snps_plink_pca.eigenvec'
    params:
        out = 'output/04_plink/no_ldpruning/filtered_snps_plink_pca',
    log:
        'output/logs/plink2_pca_no_ldpruning.log'
    singularity:
        plink2_container
    shell:
        'plink2 '
        '--vcf {input.vcf} '
        '--read-freq {input.freq} '
        '--set-all-var-ids @:# '
        '--allow-extra-chr '
        '--vcf-half-call missing '
        '--pca 20 '
        '--make-bed '
        '--out {params.out} '
        '&> {log}'



rule plink2_allele_freqs:
    input:
        'output/03_filtering/ASW_agresearch_combined_filtered.vcf.gz'
    output:
        'output/04_plink/no_ldpruning/filtered_snps_freq.afreq'
    params:
        wd = 'output/04_plink/no_ldpruning/filtered_snps_freq'
    log:
        'output/logs/plink2_allele_freqs.log'
    singularity:
        plink2_container
    shell:
        'plink2 '
        '--vcf {input} '
        '--freq '
        '--set-all-var-ids @:# '
        '--allow-extra-chr '
        '--vcf-half-call missing '
        '--out {params.wd} '
        '&> {log}'



###############
## plink fst ##
###############

rule plink2_fstmatrix:
    input:
        pgen = 'output/04_plink/fst/fst.pgen',
    output:
        'output/04_plink/fst/population_fst_results.fst.summary'
    params:
        in_path = 'output/04_plink/fst/fst',
        out_path = 'output/04_plink/fst/population_fst_results'
    log:
        'output/logs/plink_fstmatrix_ploidy_2_mitoout_population.log'
    singularity:
        plink2_container
    shell:
        'plink2 '
        '--pfile {params.in_path} '
        '--set-all-var-ids @:# '
        '--allow-extra-chr '
        '--fst region '
        '--out {params.out_path} '
        '&> {log}'



rule plink2_make_pgen_for_fst:
    input:
        vcf = 'output/03_filtering/ASW_agresearch_combined_filtered.vcf.gz',
        populations_psam = 'data/plink2_populations.psam'
    output:
        pgen = 'output/04_plink/fst/fst.pgen',
        psam  = 'output/04_plink/fst/fst.psam'
    params:
        wd = 'output/04_plink/fst/fst'
    log:
        'output/logs/plink2_make_pgen_for_fst.log'
    singularity:
        plink2_container
    shell:
        """
        plink2 \
        --vcf {input.vcf} \
        --set-all-var-ids @:# \
        --allow-extra-chr \
        --vcf-half-call missing \
        --make-pgen \
        --out {params.wd} \
        &> {log}

        cp {input.populations_psam} {output.psam}
        """


#############
# filtering #
#############

# filtration stats
rule variant_stats:
    input:
        'output/03_filtering/ASW_agresearch_{filtering}.vcf.gz'
    output:
        'output/03_filtering/vcf_stats/ASW_agresearch_{filtering}_stats.txt'
    log:
        'logs/{filtering}_bcfvariant_stats.log'
    singularity:
        bcftools_container
    shell:
        'bcftools stats {input} > {output} 2> {log}'



rule bcf_filtering:
    input:
        'output/03_filtering/ASW_agresearch_combined.vcf.gz'
    output:
        'output/03_filtering/ASW_agresearch_combined_filtered.vcf.gz'
    log:
        'logs/bcf_filtering.log'
    singularity:
        bcftools_container
    shell:
        'bcftools view '
        '-m2 -M2 ' 
        '-i "F_MISSING<0.2 & QUAL>40 & AVG(GQ)>20" ' 
        '-v snps ' 
        '-o {output} ' 
        '-q 0.05 ' 
        '&> {log} '
        '{input}'



rule bcftools_index:
    input:
        'output/03_filtering/ASW_agresearch_combined.vcf.gz'
    output:
        'output/03_filtering/ASW_agresearch_combined.vcf.gz.tbi'
    singularity:
        bcftools_container
    shell:
        'bcftools index {input} -t'



rule bcftools_convert:
    input: 
        gvcf = 'output/02_variants/ASW_agresearch_combined.gvcf.bcf',
        fa = 'data/GCA_014170235_1.fa'
    output:
        'output/03_filtering/ASW_agresearch_combined.vcf.gz'
    singularity:
        bcftools_container
    log:
        'output/logs/bcftools_convert.log'
    shell:
        'bcftools convert --gvcf2vcf -f {input.fa} {input.gvcf} -O v -o {output} &> {log}'

#####################
## variant calling ##
#####################


rule glnexus_cli:
    input:
        expand('output/02_variants/gatkhap/{sample}_gatk.g.vcf.gz', sample=all_samples)
    output:
        'output/02_variants/ASW_agresearch_combined.gvcf.bcf'
    singularity:
        glnexus_container
    log:
        'output/logs/glnexus.log'
    threads:
        16
    shell:
        'glnexus_cli -m 16 -t {threads} --config gatk {input} > {output} 2> {log}'



rule gatk_haplotypecaller:
    input:
        genome = 'data/GCA_014170235_1.fa',
        genome_dict = 'data/GCA_014170235_1.dict',
        genome_index = 'data/GCA_014170235_1.fa.fai',
        bam = 'output/02_variants/markdup/{sample}_markdup.cram',
        index = 'output/02_variants/markdup/{sample}_markdup.cram.crai'
    output:
        gvcf = 'output/02_variants/gatkhap/{sample}_gatk.g.vcf.gz'
    threads:
        3
    log:
        'output/logs/gatkhap/{sample}.log'
    singularity:
        gatk_container
    shell:
        'nice gatk --java-options "-Xmx12g" HaplotypeCaller '
        '-R {input.genome} '
        '-I {input.bam} '
        '--ploidy 2 '
        '-ERC GVCF ' # output gVCF files
        '-O {output.gvcf} '
        '2> {log}'



rule samtools_index_gatk:
    input:
        'output/02_variants/markdup/{sample}_markdup.cram'
    output:
        'output/02_variants/markdup/{sample}_markdup.cram.crai'
    log:
        'output/logs/samtools_index_gatk/{sample}.log'
    shell:
        'samtools index '
        '{input} '
        '2> {log}'




#################################
## transform bams and markdups ##
#################################

rule samtools_markdup:
    input:
        cram = 'output/02_variants/bwa_GCA_014170235_1_crams/{sample}_fixed_sorted.cram',
        fa = 'data/GCA_014170235_1.fa'
    output:
        'output/02_variants/markdup/{sample}_markdup.cram'
    threads:
        2
    log:
        'output/logs/markdup/{sample}.log'
    shell:
        'samtools markdup --threads 2 --reference {input.fa} -O CRAM {input.cram} {output}'


rule sort_bams_position:
    input:
        cram = 'output/02_variants/bwa_GCA_014170235_1_crams/{sample}_fixed.cram',
        fa = 'data/GCA_014170235_1.fa'
    output:
        temp('output/02_variants/bwa_GCA_014170235_1_crams/{sample}_fixed_sorted.cram')
    threads:
        2
    shell:
        'samtools sort --threads 2 --reference {input.fa} {input.cram} -o {output}'


rule fix_read_mates:
    input:
        cram = 'output/02_variants/bwa_GCA_014170235_1_crams/{sample}_name_sorted_cr.cram',
        fa = 'data/GCA_014170235_1.fa'
    output:
        temp('output/02_variants/bwa_GCA_014170235_1_crams/{sample}_fixed.cram')
    threads:
        2
    shell:
        'samtools fixmate -r -c -m --threads 2 --reference {input.fa} {input.cram} {output}'


rule sort_bams_name:
    input:
        'output/02_variants/bwa_GCA_014170235_1_crams/{sample}_sorted_cr.cram'
    output:
        temp('output/02_variants/bwa_GCA_014170235_1_crams/{sample}_name_sorted_cr.cram')
    threads:
        2
    shell:
        'samtools sort -n -O CRAM {input} -o {output}'


rule samtools_convert_bamcram:
    input:
        bam = 'output/02_variants/bwa_GCA_014170235_1_crams/{sample}_sorted_rgs.bam',
        fa = 'data/GCA_014170235_1.fa',
        fai = 'data/GCA_014170235_1.fa.fai'
    output:
        temp('output/02_variants/bwa_GCA_014170235_1_crams/{sample}_sorted_cr.cram')
    log:
        'output/logs/convert_bamcram/{sample}.log'
    threads:
        2
    shell:
        'samtools view -C --reference {input.fa} -t {input.fai} --threads 2 {input.bam} -o {output} 2> {log}'



rule gatk_readgroups:
    input:
        'output/01_mapping/bwa_GCA_014170235_1/{sample}_sorted.bam'
    output:
        temp('output/02_variants/bwa_GCA_014170235_1_crams/{sample}_sorted_rgs.bam')
    params:
        sample_id = '{sample}'
    log:
        'output/logs/gatk_readgroups/{sample}.log'
    singularity:
        gatk_container
    threads:
        10
    shell:
        'gatk --java-options "-Xmx8g" AddOrReplaceReadGroups '
        'I={input} O={output} '
        'RGID={params.sample_id} RGLB={params.sample_id} RGPU={params.sample_id} RGSM={params.sample_id} RGPL=ILLUMINA '
        '2> {log}'



#############################
## prep ref fasta for GATK ##
#############################

rule gatk_seq_dict:
    input:
        'data/GCA_014170235_1.fa'
    output:
        'data/GCA_014170235_1.dict'
    log:
        'output/logs/gatk_seq_dict.log'
    singularity:
        gatk_container
    shell:
        'gatk CreateSequenceDictionary -R {input} -O {output} 2> {log}'

rule faidx:
    input:
        fasta = 'data/GCA_014170235_1.fa'
    output:
        fai = 'data/GCA_014170235_1.fa.fai'
    log:
        'output/logs/faidx.log'
    shell:
        'samtools faidx {input.fasta} -o {output.fai} 2> {log} '


##########
# map qc #
##########


rule samtools_coverage:
    input:
        bam = 'output/01_mapping/bwa_GCA_014170235_1/{sample}_sorted.bam'
    output:
        coverage_out = 'output/01_mapping/samtools_GCA_014170235_1/coverage/{sample}_coverage.out'
    log:
        'output/logs/samtools_GCA_014170235_1_coverage/samtools_GCA_014170235_1_coverage_{sample}.log'
    shell:
        'nice samtools coverage '
        '{input.bam} '
        '-o {output.coverage_out} '
        '2> {log}'


rule samtools_depth:
    input:
        sorted_bam = 'output/01_mapping/bwa_GCA_014170235_1/{sample}_sorted.bam'
    output:
        depth_out = 'output/01_mapping/samtools_GCA_014170235_1/depth/{sample}_depth.out'
    log:
        'output/logs/samtools_depth/samtools_depth_GCA_014170235_1_{sample}.log'
    shell:
        'nice samtools depth '
        '{input.sorted_bam} '
        '-aa ' # write positions/contigs even when depth = 0
        '> {output.depth_out} '
        '2> {log}'

#######
# map #
#######

rule samtools_flagstat:
    input:
        'output/01_mapping/bwa_GCA_014170235_1/{sample}_sorted.bam'
    output:
        temp('output/01_mapping/bwa_GCA_014170235_1/samtools_flagstat/{sample}.out')
    log:
        'output/logs/samtools_flagstat/samtools_flagstat_{sample}_GCA_014170235_1.log'
    shell:
        'nice samtools flagstat '
        '{input} > {output} '
        '2> {log}'

rule samtools_index:
    input:
        bam = 'output/01_mapping/bwa_GCA_014170235_1/{sample}_sorted.bam'
    output:
        index = 'output/01_mapping/bwa_GCA_014170235_1/{sample}_sorted.bam.bai'
    log:
        'output/logs/samtools_index/samtools_index_{sample}.log'
    shell:
        'nice samtools index '
        '{input.bam} '
        '2> {log}'

rule samtools_sort:
    input:
        sam = 'output/01_mapping/bwa_GCA_014170235_1/{sample}_bwa_mem.sam'
    output:
        sorted_bam = 'output/01_mapping/bwa_GCA_014170235_1/{sample}_sorted.bam'
    log:
        'output/logs/samtools_sort/samtools_sort_GCA_014170235_1_{sample}.log'
    shell:
        'nice samtools sort '
        '{input.sam} '
        '-o {output.sorted_bam} '
        '2> {log}'

rule bwa_mem:
    input:
        index = 'output/01_mapping/bwa_GCA_014170235_1/index/index.bwt',
        r1 = 'output/01_mapping/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/01_mapping/bbduk_trim/{sample}_r2.fq.gz',
    output:
        sam = temp('output/01_mapping/bwa_GCA_014170235_1/{sample}_bwa_mem.sam')
    params:
        index_dir = 'output/01_mapping/bwa_GCA_014170235_1/index/index'
    threads:
        20
    log:
        'output/logs/bwa_mem/{sample}_GCA_014170235_1.log'
    singularity:
        bwa_container
    shell:
        'nice bwa mem '
        '-t {threads} '
        '{params.index_dir} '
        '{input.r1} {input.r2} '
        '> {output.sam} '
        '2> {log}'


rule bwa_index:
    input:
        genome = 'data/GCA_014170235_1.fa'
    output:
        index = 'output/01_mapping/bwa_GCA_014170235_1/index/index.bwt'
    params:
        outdir = 'output/01_mapping/bwa_GCA_014170235_1/index/index'
    threads:
        20
    log:
        'output/logs/bwa_index_GCA_014170235_1.log'
    singularity:
        bwa_container
    shell:
        'nice bwa index '
        '{input.genome} '
        '-p {params.outdir} '
        '2> {log} '



#######################
# trim and initial QC #
#######################


rule multiqc:
    input:
        fastqc = expand('output/01_mapping/fastqc/{sample}_r{n}_fastqc.zip', sample=all_samples, n=[1,2])
    output:
        'output/01_mapping/multiqc/multiqc_report.html'
    params:
        outdir = 'output/01_mapping/multiqc',
        indirs = ['output/01_mapping/fastqc']
    log:
        'output/logs/multiqc.log'
    container:
        multiqc_container
    shell:
        'multiqc '
        '-f '
        '-o {params.outdir} '
        '{params.indirs} '
        '2> {log}'


rule fastqc:
    input:
        expand('output/01_mapping/bbduk_trim/{sample}_r{n}.fq.gz',
            sample=all_samples, n=[1,2])
    output:
        temp(expand('output/01_mapping/fastqc/{sample}_r{n}_fastqc.zip', sample=all_samples, n=[1,2]))
    params:
        outdir = directory('output/01_mapping/fastqc')
    log:
        'output/logs/fastqc.log'
    singularity:
        fastqc_container
    shell:
        'mkdir -p {params.outdir} ; '
        'fastqc --outdir {params.outdir} {input} 2> {log}'


rule bbduk_trim:
    input:
        r1 = 'output/joined/{sample}_r1.fq.gz',
        r2 = 'output/joined/{sample}_r2.fq.gz'
    output:
        r1 = 'output/01_mapping/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/01_mapping/bbduk_trim/{sample}_r2.fq.gz',
        stats = 'output/01_mapping/bbduk_trim/stats/{sample}_stats.txt'
    params:
        adapters = '/adapters.fa'
    log:
        'output/logs/bbduk_trim/{sample}.log'
    singularity:
        bbduk_container
    threads:
        80
    shell:
        'nice bbduk.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo forcetrimmod=5 qtrim=r trimq=15 '
        'stats={output.stats} '
        'threads={threads} '
        '&> {log}'





