from os.path import join
from collections import defaultdict
import pandas as pd
import glob
outdir = "/oak/stanford/groups/smontgom/tannerj/ADRC_LRS/output"
scratch="/tmp/tannerj"
metafora_dir = "/oak/stanford/groups/smontgom/tannerj/Metafora"
sample_table = pd.read_table("./PBMC_ADRC_sample_table.txt")
samples = sample_table['Sample_name'].to_list()
tissues = sample_table['Tissue'].to_list()

valid_chroms=["chr%d" % i for i in range(1,23)]
tissue_dict= defaultdict(list)
for sample,tiss in zip(samples,tissues):
  tissue_dict[tiss].append(sample)

rule all:
    input:
        expand(join(outdir,"{sample}/methylation_calls/{sample}.{ref_name}.tissue_{tissue}.METAFORA.outlier_regions.ALL_AUTOSOMES.bed"),
                tissue="PBMC",
               ref_name="GRCh38",
               sample=samples)
        #expand(join(outdir, "methylation_results/Meth_segments.{ref_name}.tissue_{tissue}.segment_betas.ALL_AUTOSOMES.bed"),
        #       tissue="PBMC",
        #       ref_name="GRCh38")

rule nanopolish_cpg_betas:
    threads: 1
    resources: 
      time=8, 
      mem=32
    input:
      join(outdir, "{sample}/methylation_calls/{sample}.{ref_name}.f5c.methylation_calls.tsv.gz"),
    params:
      tmp_dir = join(scratch, "{sample}"),
      tmp_calls = join(scratch, "{sample}/{sample}.{ref_name}.tmp.methylation.tsv"),
      tmp_betas = join(scratch, "{sample}/{sample}.{ref_name}.tmp.cpg_sites.betas.tsv"),
      script = join(metafora_dir,"calculate_methylation_frequency.py"),
      LLR_thresh = 1.5
    output:
      bed = join(outdir, "{sample}/methylation_calls/{sample}.{ref_name}.cpg_sites.betas.tsv.gz"),
      tbi = join(outdir, "{sample}/methylation_calls/{sample}.{ref_name}.cpg_sites.betas.tsv.gz.tbi")
    conda: "r"
    shell: """
       mkdir -p {params.tmp_dir}
       zcat {input} > {params.tmp_calls}
       python {params.script} -c {params.LLR_thresh} -s {params.tmp_calls} > {params.tmp_betas} 
       bgzip -c {params.tmp_betas} > {output.bed}
       tabix -p bed -S 1 {output.bed}
       rm -rf {params.tmp_dir}
    """

rule create_tissue_sample_reference:
  threads: 1
  resources:
    time=48,
    mem=124
  input:
    lambda wildcards: expand(join(outdir, "{sample}/methylation_calls/{sample}.{{ref_name}}.cpg_sites.betas.tsv.gz"),
                             sample=list(tissue_dict[wildcards.tissue])),
  params:
    bed = "/oak/stanford/groups/euan/projects/tannerj/UDN/methylation_results/GRCh38.CpG_sites.bed",
    filelist=join(outdir,"tmp.{tissue}.{ref_name}.file_list.txt"),
    script = join(metafora_dir,"calculate_population_mean_betas.R"),
    tmp_dir = join(scratch, "create_pop_mean_{tissue}_{chr}"),
    tmp_beta = join(scratch, "create_pop_mean_{tissue}_{chr}/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.betas.mat"),
    tmp_depth = join(scratch, "create_pop_mean_{tissue}_{chr}/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.coverage.mat"),
    tmp_mean = join(scratch, "create_pop_mean_{tissue}_{chr}/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv"),
  output:
    beta_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.betas.mat.gz"),
    beta_tbi = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.betas.mat.gz.tbi"),
    depth_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.coverage.mat.gz"),
    depth_tbi = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.coverage.mat.gz.tbi"),
    mean = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv.gz"),
    mean_tbi = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv.gz.tbi"),
  conda: "r"
  shell: """
    ls {input} > {params.filelist}
    mkdir -p {params.tmp_dir}
    Rscript {params.script} --filelist {params.filelist} \
        --chrom {wildcards.chr} \
        --cpgs {params.bed} \
        --beta_mat {params.tmp_beta} \
        --depth_mat {params.tmp_depth} \
        --pop_mean {params.tmp_mean} 

    bgzip -c {params.tmp_beta} > {output.beta_mat}
    bgzip -c {params.tmp_depth} > {output.depth_mat}
    bgzip -c {params.tmp_mean} > {output.mean}

    tabix -p bed -S 1 {output.beta_mat}
    tabix -p bed -S 1 {output.depth_mat}
    tabix -p bed -S 1 {output.mean}
    rm -rf {params.tmp_dir}
  """

rule segment_mean_beta:
  threads: 1
  resources:
    time=8,
    mem=124
  input:
    pop_mean = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv.gz"),
    beta_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.betas.mat.gz"),
    depth_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.coverage.mat.gz")
  params:
    script = join(metafora_dir,"segment_population_mean_betas.R")
  output:
    seg_bed = join(outdir, "methylation_results/Meth_segments.{ref_name}.tissue_{tissue}.chrom_{chr}.segment_betas.bed")
  conda: "r"
  shell: """
    Rscript {params.script} \
        --chrom {wildcards.chr} \
        --population_beta {input.pop_mean} \
        --beta_mat {input.beta_mat} \
        --depth_mat {input.depth_mat} \
        --segment_bed_out {output.seg_bed} \
  """

rule combine_chrom_segments:
  threads: 1 
  resources:
    time=4,
    mem=36
  input:
    expand(join(outdir, "methylation_results/Meth_segments.{{ref_name}}.tissue_{{tissue}}.chrom_{chr}.segment_betas.bed"),
           chr=valid_chroms)
  output:
    seg_bed = join(outdir, "methylation_results/Meth_segments.{ref_name}.tissue_{tissue}.segment_betas.ALL_AUTOSOMES.bed"),
  shell: """ 
    head -1 {input[0]} > {output}
    cat {input} | grep -v "^chrom" >> {output}
  """

rule find_variable_segments:
  threads: 1
  resources:
    time=8,
    mem=124
  input:
    seg_bed = join(outdir, "methylation_results/Meth_segments.{ref_name}.tissue_{tissue}.segment_betas.ALL_AUTOSOMES.bed"),
  params:
    script = join(metafora_dir,"find_variable_segments.R")
  output:
    varseg_bed = join(outdir, "methylation_results/Meth_segments.{ref_name}.tissue_{tissue}.segment_betas.ALL_AUTOSOMES.variable_segments.bed"),
    varseg_plot_file = join(outdir, "methylation_results/Variable_segments.{ref_name}.tissue_{tissue}.all_autosomes.segment_CV2.pdf"),
    pca_out = join(outdir, "methylation_results/Gloabal_Methylation_PCA_{ref_name}_tissue_{tissue}/Global_meth_pcs.{ref_name}.tissue_{tissue}.txt")
  conda: "r"
  shell: """
    pca_plot_dir=$(dirname {output.pca_out})
    Rscript {params.script} \
        --segment_betas {input.seg_bed} \
        --varsegment_bed_out {output.varseg_bed} \
        --varseg_plot_file {output.varseg_plot_file} \
        --pca_plot_dir $pca_plot_dir \
        --global_meth_pcs_out {output.pca_out} 
  """

rule call_outliers:
  threads: 1
  resources:
    time=16,
    mem=96
  input:
    pop_mean = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv.gz"),
    beta_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.betas.mat.gz"),
    depth_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.chrom_{chr}.coverage.mat.gz"),
    global_pcs = join(outdir, "methylation_results/Gloabal_Methylation_PCA_{ref_name}_tissue_{tissue}/Global_meth_pcs.{ref_name}.tissue_{tissue}.txt")
  params:
    script = join(metafora_dir,"call_methylation_outliers.R")
  output:
    outlier_bed = temp(join(outdir,"{sample}/methylation_calls/{sample}.{ref_name}.tissue_{tissue}.METAFORA.chrom_{chr}.outlier_regions.bed")),
    outlier_z_mat = temp(join(outdir, "{sample}/methylation_calls/{sample}.{ref_name}.tissue_{tissue}.METAFORA.chrom_{chr}.outlier_regions.zscore.mat"))
  conda: 'r'
  shell: """ 
    outlier_plot_dir="$(dirname {output.outlier_bed})/outlier_plots"
    mkdir -p $outlier_plot_dir
    Rscript {params.script} \
        --chrom {wildcards.chr} \
        --sample {wildcards.sample} \
        --population_beta {input.pop_mean} \
        --beta_mat {input.beta_mat} \
        --depth_mat {input.depth_mat} \
        --global_meth_pcs {input.global_pcs} \
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
    beds = expand(join(outdir, "{{sample}}/methylation_calls/{{sample}}.{{ref_name}}.tissue_{{tissue}}.METAFORA.chrom_{chr}.outlier_regions.bed"),
           chr=valid_chroms),
    mats = expand(join(outdir, "{{sample}}/methylation_calls/{{sample}}.{{ref_name}}.tissue_{{tissue}}.METAFORA.chrom_{chr}.outlier_regions.zscore.mat"),
           chr=valid_chroms)
  output:
    outlier_bed = join(outdir,"{sample}/methylation_calls/{sample}.{ref_name}.tissue_{tissue}.METAFORA.outlier_regions.ALL_AUTOSOMES.bed"),
    outlier_z_mat = join(outdir, "{sample}/methylation_calls/{sample}.{ref_name}.tissue_{tissue}.METAFORA.outlier_regions.ALL_AUTOSOMES.zscore.mat")
  shell: """ 
    header=$(for i in $(ls {input.beds}); do head -1 $i; done | sort | tail -1)
    echo -e "$header" > {output.outlier_bed}
    cat {input.beds} | grep -v "$header" | grep -v "^$" >> {output.outlier_bed}

    header=$(for i in $(ls {input.mats}); do head -1 $i; done | sort | tail -1) 
    echo -e "$header" > {output.outlier_z_mat}
    cat {input.mats} | grep -v "$header" | grep -v "^$" >> {output.outlier_z_mat}
  """
#
#rule combine_outliers:
#  threads: 1 
#  resources:
#     time=1,
#     mem=25
#  input:
#     expand(join(outdir,"{sample}/methylation_calls/{sample}.GRCh38.tissue_{tissue}.METAFORA.all_autosomes.outlier_regions.bed"),zip,sample=samples,tissue=tissues)
#  output:
#     join(outdir, "methylation_results/ADRC.Methylation_outliers.combined.bed")
#  shell: """ 
#      echo -e "chrom\tstart\tend\twidth\tstrand\tID\tnum.mark\tseg.mean\tstartRow\tendRow\tseg_id\tzscore " > {output}
#      cat {input} | grep -v "^seqnames" >> {output}
#  """
