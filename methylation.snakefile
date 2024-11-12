from os.path import join
from collections import defaultdict
import pandas as pd
import glob
outdir = config["output_dir"] 
scratch = config["scratch_dir"]

sample_table = pd.read_table(config["sample_table"])
samples = sample_table['Sample_name'].to_list()
tissue_dict=defaultdict(list)

sample_tissue_map = { row.Sample_name : row.Tissue for i,row in sample_table.iterrows() }
technology_map = { row.Sample_name : row.Technology for i,row in sample_table.iterrows() }
input_map = { row.Sample_name : row.Methylation_input for i,row in sample_table.iterrows()}
for sample,tissue in sample_tissue_map.items():
  tissue_dict[tissue].append(sample)

for tissue in tissue_dict.keys():
  print("tissue: %s , N = %d" % (tissue, len(tissue_dict[tissue])))

unique_tissues = list(tissue_dict.keys())
print(unique_tissues)

autosomes = ["chr%d" % i for i in range(1,23)]
sex_chroms = ["chrX", "chrY"]

#set up default params
if not "params" in config:
  config["params"] = {}
MAX_DEPTH = config["params"]["MAX_DEPTH"] if "MAX_DEPTH" in config["params"] else 30
MIN_SEG_SIZE = config["params"]["MIN_SEG_SIZE"] if "MIN_SEG_SIZE" in config["params"] else 20
MIN_ABS_ZSCORE = config["params"]["MIN_ABS_ZSCORE"] if "MIN_ABS_ZSCORE" in config["params"] else 3
MIN_ABS_DELTA = config["params"]["MIN_ABS_DELTA"] if "MIN_ABS_DELTA" in config["params"] else 0.3

rule all:
    input:
        expand(join(outdir,"{sample}/methylation/{sample}.GRCh38.tissue_{tissue}.METAFORA.outlier_regions.ALL_AUTOSOMES.bed"),
               tissue=unique_tissues, sample=samples)
        #expand(join(outdir, "methylation_results/UDN_Cohort.tissue_{tissue}.Methylation_outliers.combined.tsv"), tissue=[unique_tissues]),

rule modBam2Bed:
    threads: 4
    resources:
        mem=32,
        time=12
    input:
        bam = lambda w: input_map[w.sample],
        ref = lambda w: config["reference_params"]["GRCh38"]["fasta"]
    params:
        tmp_dir = join(scratch, "{sample}_GRCh38_methylation"),
        tmp_bam = join(scratch, "{sample}_GRCh38_methylation/primary_alignment.bam"),
        tmp_bed = join(scratch, "{sample}_GRCh38_methylation/out.bed") 
    output:
        meth_bed = join(outdir, "{sample}/methylation/{sample}.GRCh38.ONT.cpg_methylation.bed.gz"),
        tbi = join(outdir, "{sample}/methylation/{sample}.GRCh38.ONT.cpg_methylation.bed.gz.tbi")
    conda: '/oak/stanford/groups/smontgom/tannerj/RUSH_AD/envs/methylation.yaml'
    shell: """
        # in process of updating scripts
        mkdir -p {params.tmp_dir}
        #ml samtools #TODO: add to conda yaml in future
        #samtools view -b -F 256 {input.bam} > {params.tmp_bam} #filter out non-primary alignments
        #samtools index {params.tmp_bam}
        modbam2bed -t {threads} --cpg --combine -m 5mC -d 40 {input.ref} {input.bam} > {params.tmp_bed}
        bgzip -c {params.tmp_bed} > {output.meth_bed}
        tabix {output.meth_bed}
        rm -rf {params.tmp_dir}
    """

rule format_modBam2Bed:
    threads: 16
    resources:
      time=4,
      mem=128
    input:
      join(outdir, "{sample}/methylation/{sample}.GRCh38.ONT.cpg_methylation.bed.gz")
    params:
      script = "format_modBam2Bed.R",
      meth_bed = join(outdir, "{sample}/methylation/{sample}.GRCh38.ONT.METAFORA_formatted.cpg_methylation.bed")
    output:
      meth_bed = join(outdir, "{sample}/methylation/{sample}.GRCh38.ONT.METAFORA_formatted.cpg_methylation.bed.gz"),
      tbi = join(outdir, "{sample}/methylation/{sample}.GRCh38.ONT.METAFORA_formatted.cpg_methylation.bed.gz.tbi")
    conda: 'r'
    shell: """
      Rscript {params.script} --input {input} --output {params.meth_bed}
      bgzip {params.meth_bed}
      tabix -p bed -S 1 {output.meth_bed}
    """

rule pacbio_cpg_tools:
  threads: 16
  resources:
    time=4,
    mem=32
  input:
    bam = lambda w: bam_map[w.sample],
    ref = lambda w: config["reference_params"]["GRCh38"]["fasta"]
  params:
    cpg_tools_command = config["pb_cpg_tools"]["bin"],
    model = config["pb_cpg_tools"]["model"],
    prefix = join(outdir, "{sample}/methylation/{sample}.GRCh38.PacBio.cpg_methylation"),
  output:
    meth_bed = join(outdir, "{sample}/methylation/{sample}.GRCh38.PacBio.cpg_methylation.combined.bed")
  shell: """
    {params.cpg_tools_command} \
        --bam {input.bam} \
        --output-prefix {params.prefix} \
        --model {params.model} \
        --pileup-mode model \
        --modsites-mode reference \
        --ref {input.ref} \
        --threads {threads}
  """

rule format_pacbio:
  threads: 1 
  resources:
    time=4,
    mem=24
  input:
    join(outdir, "{sample}/methylation/{sample}.GRCh38.PacBio.cpg_methylation.combined.bed")
  params:
    tmp_bed = join(outdir, "{sample}/methylation/{sample}.GRCh38.PacBio.METAFORA_formatted.cpg_methylation.bed")
  output:
    meth_bed = join(outdir, "{sample}/methylation/{sample}.GRCh38.PacBio.METAFORA_formatted.cpg_methylation.bed.gz"),
    tbi = join(outdir, "{sample}/methylation/{sample}.GRCh38.PacBio.METAFORA_formatted.cpg_methylation.bed.gz.tbi")
  conda: '/oak/stanford/groups/smontgom/tannerj/RUSH_AD/envs/methylation.yaml'
  shell: """ 
    echo -e "chromosome\tstart\tend\tdepth\tbeta" > {params.tmp_bed}
    awk '{{print $1,$2,$3,$6,($9/100)}}' {input} | sed 's/ /\t/g' >> {params.tmp_bed} 

    bgzip {params.tmp_bed} 
    tabix -p bed -S 1 {output.meth_bed}
  """
  
rule create_tissue_sample_reference:
  threads: 1
  resources:
    time=24,
    mem=128
  input:
    lambda wildcards: expand(join(outdir, "{sample}/methylation/{sample}.GRCh38.{tech}.METAFORA_formatted.cpg_methylation.bed.gz"),
                             zip, sample=tissue_dict[wildcards.tissue], tech=[ technology_map[s] for s in tissue_dict[wildcards.tissue]]),
  params:
    bed = config["reference_cpg_bed"],
    filelist=join(outdir,"tmp.{tissue}.GRCh38.file_list.chrom_{chr}.txt"),
    script = "calculate_population_mean_betas.R",
    tmp_beta = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.betas.mat"),
    tmp_depth = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.coverage.mat"),
    tmp_mean = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv"),
    min_segment_cpgs = MIN_SEG_SIZE
  output:
    beta_mat = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.betas.mat.gz"),
    beta_tbi = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.betas.mat.gz.tbi"),
    depth_mat = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.coverage.mat.gz"),
    depth_tbi = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.coverage.mat.gz.tbi"),
    mean = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv.gz"),
    mean_tbi = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv.gz.tbi"),
    seg_beta = join(outdir,"methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Meth_segments.GRCh38.tissue_{tissue}.segment_betas.chrom_{chr}.bed"),
    seg_depth = join(outdir,"methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Meth_segments.GRCh38.tissue_{tissue}.segment_coverage.chrom_{chr}.mat")
  conda: "r"
  shell: """
    ls {input} > {params.filelist}
    Rscript {params.script} --filelist {params.filelist} \
        --chrom {wildcards.chr} \
        --cpgs {params.bed} \
        --min_segment_cpgs {params.min_segment_cpgs} \
        --beta_mat {params.tmp_beta} \
        --depth_mat {params.tmp_depth} \
        --pop_mean {params.tmp_mean} \
        --segment_beta {output.seg_beta} \
        --segment_depth {output.seg_depth}
    rm {params.filelist}

    bgzip {params.tmp_beta}
    bgzip {params.tmp_depth}
    bgzip {params.tmp_mean}

    tabix -p bed -S 1 {output.beta_mat}
    tabix -p bed -S 1 {output.depth_mat}
    tabix -p bed -S 1 {output.mean}
  """

rule compute_hidden_factors:
  threads: 16
  resources:
    time=12,
    mem=128
  input:
    seg_beta = expand(join(outdir,"methylation_results/Population_methylation.GRCh38.tissue_{{tissue}}/Meth_segments.GRCh38.tissue_{{tissue}}.segment_betas.chrom_{chr}.bed"),chr=autosomes + sex_chroms),
    seg_depth = expand(join(outdir,"methylation_results/Population_methylation.GRCh38.tissue_{{tissue}}/Meth_segments.GRCh38.tissue_{{tissue}}.segment_coverage.chrom_{chr}.mat"),chr=autosomes + sex_chroms)
  params:
    script = "compute_hidden_factors.R",
    seg_beta = ','.join(list(expand(join(outdir,"methylation_results/Population_methylation.GRCh38.tissue_{{tissue}}/Meth_segments.GRCh38.tissue_{{tissue}}.segment_betas.chrom_{chr}.bed"),chr=autosomes))),
    seg_depth = ','.join(list(expand(join(outdir,"methylation_results/Population_methylation.GRCh38.tissue_{{tissue}}/Meth_segments.GRCh38.tissue_{{tissue}}.segment_coverage.chrom_{chr}.mat"),chr=autosomes))),
    sex_seg_beta = ','.join(list(expand(join(outdir,"methylation_results/Population_methylation.GRCh38.tissue_{{tissue}}/Meth_segments.GRCh38.tissue_{{tissue}}.segment_betas.chrom_{chr}.bed"),chr=sex_chroms))),
    sex_seg_depth = ','.join(list(expand(join(outdir,"methylation_results/Population_methylation.GRCh38.tissue_{{tissue}}/Meth_segments.GRCh38.tissue_{{tissue}}.segment_coverage.chrom_{chr}.mat"),chr=sex_chroms))),
    covariates_arg = ("--covariates %s" % config["covariates"]) if "covariates" in config else ""
  output:
    global_pcs = join(outdir, "methylation_results/Global_Methylation_PCA_GRCh38_tissue_{tissue}/PCA_covariates.txt"),
    correlation_plot = join(outdir, "methylation_results/Global_Methylation_PCA_GRCh38_tissue_{tissue}/Pairwise_sample_correlations.pdf"),
    pc_plot = join(outdir, "methylation_results/Global_Methylation_PCA_GRCh38_tissue_{tissue}/PCA_biplot.pdf"),
    sex_plots = join(outdir, "methylation_results/Global_Methylation_PCA_GRCh38_tissue_{tissue}/Sex_chromosome_estimate_plots.pdf"),
    sex_summary_table = join(outdir, "methylation_results/Global_Methylation_PCA_GRCh38_tissue_{tissue}/Sex_chromosome_estimates_summary.txt")
  conda: "r"
  shell: """
    Rscript {params.script} \
        --seg_beta {params.seg_beta} \
        --seg_depth {params.seg_depth} \
        --sex_seg_beta {params.sex_seg_beta} \
        --sex_seg_depth {params.sex_seg_depth} \
        --correlation_plot_out {output.correlation_plot} \
        --pc_biplot_out {output.pc_plot} \
        --sex_plot_out {output.sex_plots} \
        --sex_chrom_summary_out {output.sex_summary_table} \
        --global_meth_pcs_out {output.global_pcs} \
        {params.covariates_arg}
  """

rule call_outliers:
  threads: 4
  resources:
    time=16,
    mem=128
  input:
    beta_mat = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.betas.mat.gz"),
    depth_mat = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.coverage.mat.gz"),
    global_pcs = join(outdir, "methylation_results/Global_Methylation_PCA_GRCh38_tissue_{tissue}/PCA_covariates.txt")
    #global_pcs = join(outdir, "methylation_results/Global_Methylation_PCA_GRCh38_tissue_{tissue}/Global_meth_pcs.GRCh38.tissue_{tissue}.txt")
  params:
    script = "call_methylation_outliers.R",
    MAX_DEPTH = MAX_DEPTH, 
    MIN_SEG_SIZE = MIN_SEG_SIZE, 
    MIN_ABS_ZSCORE = MIN_ABS_ZSCORE, 
    MIN_ABS_DELTA = MIN_ABS_DELTA 
  output:
    outlier_bed = temp(join(outdir,"{sample}/methylation/{sample}.GRCh38.tissue_{tissue}.METAFORA.outlier_regions.chrom_{chr}.bed")),
    outlier_z_mat = temp(join(outdir, "{sample}/methylation/{sample}.GRCh38.tissue_{tissue}.METAFORA.outlier_region.zscore.chrom_{chr}.mat"))
  conda: 'r'
  shell: """ 
    outlier_plot_dir="$(dirname {output.outlier_bed})/outlier_plots"
    mkdir -p $outlier_plot_dir
    Rscript {params.script} \
        --sample {wildcards.sample} \
        --chrom {wildcards.chr} \
        --beta_mat {input.beta_mat} \
        --depth_mat {input.depth_mat} \
        --global_meth_pcs {input.global_pcs} \
        --outlier_bed {output.outlier_bed} \
        --outlier_z_mat {output.outlier_z_mat} \
        --plot_dir $outlier_plot_dir \
        --min_seg_size {params.MIN_SEG_SIZE} \
        --min_abs_zscore {params.MIN_ABS_ZSCORE} \
        --min_abs_delta {params.MIN_ABS_DELTA} \
        --max_depth {params.MAX_DEPTH}
  """

rule call_outliers_sex_chroms:
  threads: 4
  resources:
    time=16,
    mem=128
  input:
    beta_mat = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.betas.mat.gz"),
    depth_mat = join(outdir, "methylation_results/Population_methylation.GRCh38.tissue_{tissue}/Population_methylation.GRCh38.tissue_{tissue}.chrom_{chr}.coverage.mat.gz"),
    global_pcs = join(outdir, "methylation_results/Global_Methylation_PCA_GRCh38_tissue_{tissue}/PCA_covariates.txt"),
    sex_chrom_df = join(outdir, "methylation_results/Global_Methylation_PCA_GRCh38_tissue_{tissue}/Sex_chromosome_estimates.tsv")
    #global_pcs = join(outdir, "methylation_results/Global_Methylation_PCA_GRCh38_tissue_{tissue}/Global_meth_pcs.GRCh38.tissue_{tissue}.txt")
  params:
    script = "call_methylation_outliers.sex_chromosomes.R"
  output:
    outlier_bed = temp(join(outdir,"{sample}/methylation/{sample}.GRCh38.tissue_{tissue}.METAFORA.outlier_regions.sex_chroms.chrom_{chr}.bed")),
    outlier_z_mat = temp(join(outdir, "{sample}/methylation/{sample}.GRCh38.tissue_{tissue}.METAFORA.outlier_region.zscore.sex_chroms.chrom_{chr}.mat"))
  conda: 'r'
  shell: """ 
    outlier_plot_dir="$(dirname {output.outlier_bed})/outlier_plots"
    mkdir -p $outlier_plot_dir
    Rscript {params.script} \
        --sample {wildcards.sample} \
        --chrom {wildcards.chr} \
        --beta_mat {input.beta_mat} \
        --depth_mat {input.depth_mat} \
        --global_meth_pcs {input.global_pcs} \
        --sex_chromosome_tsv {input.sex_chrom_df} \
        --outlier_bed {output.outlier_bed} \
        --outlier_z_mat {output.outlier_z_mat} \
        --plot_dir $outlier_plot_dir
  """

rule combine_chrom_outliers:
  threads: 1 
  resources:
    time=4,
    mem=36
  input:
    beds = expand(join(outdir, "{{sample}}/methylation/{{sample}}.GRCh38.tissue_{{tissue}}.METAFORA.outlier_regions.chrom_{chr}.bed"), chr=autosomes),
    mats = expand(join(outdir, "{{sample}}/methylation/{{sample}}.GRCh38.tissue_{{tissue}}.METAFORA.outlier_region.zscore.chrom_{chr}.mat"), chr=autosomes) 
  output:
    outlier_bed = join(outdir,"{sample}/methylation/{sample}.GRCh38.tissue_{tissue}.METAFORA.outlier_regions.ALL_AUTOSOMES.bed"),
    outlier_z_mat = join(outdir, "{sample}/methylation/{sample}.GRCh38.tissue_{tissue}.METAFORA.outlier_regions.ALL_AUTOSOMES.zscore.mat")
  shell: """ 
    header=$(for i in $(ls {input.beds}); do head -1 $i; done | sort | tail -1)
    echo -e "$header" > {output.outlier_bed}
    cat {input.beds} | grep -v "$header" | grep -v "^$" >> {output.outlier_bed} || true

    header=$(for i in $(ls {input.mats}); do head -1 $i; done | sort | tail -1) 
    echo -e "$header" > {output.outlier_z_mat}
    cat {input.mats} | grep -v "$header" | grep -v "^$" >> {output.outlier_z_mat} || true
  """

rule combine_sex_chrom_outliers:
  threads: 1 
  resources:
    time=4,
    mem=36
  input:
    beds = expand(join(outdir, "{{sample}}/methylation/{{sample}}.GRCh38.tissue_{{tissue}}.METAFORA.outlier_regions.sex_chroms.chrom_{chr}.bed"), chr=sex_chroms),
    mats = expand(join(outdir, "{{sample}}/methylation/{{sample}}.GRCh38.tissue_{{tissue}}.METAFORA.outlier_region.zscore.sex_chroms.chrom_{chr}.mat"), chr=sex_chroms)
  output:
    outlier_bed = join(outdir,"{sample}/methylation/{sample}.GRCh38.tissue_{tissue}.METAFORA.outlier_regions.ALL_SEX_CHROM.bed"),
    outlier_z_mat = join(outdir, "{sample}/methylation/{sample}.GRCh38.tissue_{tissue}.METAFORA.outlier_regions.ALL_SEX_CHROM.zscore.mat")
  shell: """ 
    header=$(for i in $(ls {input.beds}); do head -1 $i; done | sort | tail -1)
    echo -e "$header" > {output.outlier_bed}
    cat {input.beds} | grep -v "$header" | grep -v "^$" >> {output.outlier_bed} || true

    header=$(for i in $(ls {input.mats}); do head -1 $i; done | sort | tail -1) 
    echo -e "$header" > {output.outlier_z_mat}
    cat {input.mats} | grep -v "$header" | grep -v "^$" >> {output.outlier_z_mat} || true
  """

rule combine_outliers:
  threads: 1 
  resources:
     time=1,
     mem=25
  input:
    lambda w: expand(join(outdir,"{sample}/methylation/{sample}.GRCh38.tissue_{{tissue}}.METAFORA.outlier_regions.ALL_AUTOSOMES.bed"),
                             sample=tissue_dict[w.tissue]),
    lambda w: expand(join(outdir,"{sample}/methylation/{sample}.GRCh38.tissue_{{tissue}}.METAFORA.outlier_regions.ALL_SEX_CHROM.bed"),
                             sample=tissue_dict[w.tissue]),
  output:
     join(outdir, "methylation_results/UDN_Cohort.tissue_{tissue}.Methylation_outliers.combined.tsv")
  shell: """ 
      header=$(for i in $(ls {input}); do head -1 $i; done | sort | tail -1)
      echo "$header" > {output}
      cat {input} | grep -v "$header" | grep -v "^$" >> {output} || true
  """
