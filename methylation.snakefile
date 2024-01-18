from os.path import join
from collections import defaultdict
import pandas as pd
import glob
outdir = "./output"
scratch="/tmp/tannerj"
sample_table = pd.read_table("./DORADO_input.bams.txt")
samples = sample_table['Sample_name'].to_list()

tissue_dict=defaultdict(list)
for sample in samples: 
  _,tissue=sample.split('_')
  tissue_dict[tissue].append(sample)

tissues = [sample.split('_')[1] for sample in samples]

rule all:
    input:
        expand(join(outdir, "{sample}/methylation/{sample}.{ref_name}.METAFORA_formatted.cpg_methylation.bed"),
          sample="UDN496308_Blood",
          ref_name="GRCh38")
        #join(outdir, "methylation_results/UDN_Cohort.Methylation_outliers.combined.bed")

rule phased_modBam2Bed:
    threads: 4
    resources:
        mem=32,
        time=12
    input:
        bam = join(outdir, "{sample}/{sample}.{ref_name}.methylated.sorted.whastshap_PHASED.bam"),
        bai = join(outdir, "{sample}/{sample}.{ref_name}.methylated.sorted.whastshap_PHASED.bam.bai"),
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"]
    params:
        tmp_dir = join(scratch, "{sample}_{ref_name}_methylation"),
        tmp_bam = join(scratch, "{sample}_{ref_name}_methylation/primary_alignment.bam"),
        tmp_bed = join(scratch, "{sample}_{ref_name}_methylation/out.bed") 
    output:
        meth_bed = join(outdir, "{sample}/methylation/{sample}.{ref_name}.cpg_methylation.bed.gz"),
        tbi = join(outdir, "{sample}/methylation/{sample}.{ref_name}.cpg_methylation.bed.gz.tbi")
    conda: '/oak/stanford/groups/smontgom/tannerj/RUSH_AD/envs/methylation.yaml'
    shell: """
        # in process of updating scripts
        mkdir -p {params.tmp_dir}
        ml samtools #TODO: add to conda yaml in future
        samtools view -b -F 256 {input.bam} > {params.tmp_bam} #filter out non-primary alignments
        samtools index {params.tmp_bam}
        modbam2bed -t {threads} --cpg --combine -m 5mC -d 40 {input.ref} {params.tmp_bam} > {params.tmp_bed}
        bgzip -c {params.tmp_bed} > {output.meth_bed}
        tabix {output.meth_bed}
        rm -rf {params.tmp_dir}
    """

rule modBam2Bed:
    threads: 4
    resources:
        mem=32,
        time=12
    input:
        bam = join(outdir, "{sample}/{sample}.{ref_name}.methylated.sorted.bam"),
        bai = join(outdir, "{sample}/{sample}.{ref_name}.methylated.sorted.bam.bai"),
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"]
    params:
        tmp_dir = join(scratch, "{sample}_{ref_name}_methylation"),
        tmp_bam = join(scratch, "{sample}_{ref_name}_methylation/primary_alignment.bam"),
        tmp_bed = join(scratch, "{sample}_{ref_name}_methylation/out.bed") 
    output:
        meth_bed = join(outdir, "{sample}/methylation/{sample}.{ref_name}.cpg_methylation.bed.gz"),
        tbi = join(outdir, "{sample}/methylation/{sample}.{ref_name}.cpg_methylation.bed.gz.tbi")
    conda: '/oak/stanford/groups/smontgom/tannerj/RUSH_AD/envs/methylation.yaml'
    shell: """
        # in process of updating scripts
        mkdir -p {params.tmp_dir}
        ml samtools #TODO: add to conda yaml in future
        samtools view -F 256 {input.bam} > {params.tmp_bam} #filter out non-primary alignments
        modbam2bed -t {threads} --cpg --combine -m 5mC -d 40 {input.ref} {params.tmp_bam} > {params.tmp_bed}
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
      join(outdir, "{sample}/methylation/{sample}.{ref_name}.cpg_methylation.bed.gz")
    params:
      script = "format_modBam2Bed.R"
    output:
      meth_bed = temp(join(outdir, "{sample}/methylation/{sample}.{ref_name}.METAFORA_formatted.cpg_methylation.bed"))
    conda: 'r'
    shell: """
      Rscript {params.script} --input {input} --output {output}
    """

rule create_tissue_sample_reference:
  threads: 1
  resources:
    time=16,
    mem=128
  input:
    lambda wildcards: expand(join(outdir, "{sample}/methylation/{sample}.{{ref_name}}.METAFORA_formatted.cpg_methylation.bed.gz"),
                             sample=tissue_dict[wildcards.tissue]),
  params:
    bed = "/oak/stanford/groups/euan/projects/tannerj/UDN/methylation_results/GRCh38.CpG_sites.bed",
    filelist=join(outdir,"tmp.{tissue}.{ref_name}.file_list.txt"),
    script = "calculate_population_mean_betas.R",
    tmp_beta = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.betas.mat"),
    tmp_depth = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.coverage.mat"),
    tmp_mean = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.population_mean_betas.tsv"),
  output:
    beta_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.betas.mat.gz"),
    beta_tbi = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.betas.mat.gz.tbi"),
    depth_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.coverage.mat.gz"),
    depth_tbi = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.coverage.mat.gz.tbi"),
    mean = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.population_mean_betas.tsv.gz"),
    mean_tbi = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.population_mean_betas.tsv.gz.tbi"),
  conda: "r"
  shell: """
    ls {input} > {params.filelist}
    Rscript {params.script} --filelist {params.filelist} \
        --cpgs {params.bed} \
        --beta_mat {params.tmp_beta} \
        --depth_mat {params.tmp_depth} \
        --pop_mean {params.tmp_mean}  
    rm {params.filelist}

    module load tabix
    bgzip {params.tmp_beta}
    bgzip {params.tmp_depth}
    bgzip {params.tmp_mean}

    tabix -p bed -S 1 {output.beta_mat}
    tabix -p bed -S 1 {output.depth_mat}
    tabix -p bed -S 1 {output.mean}
  """

rule segment_mean_beta:
  threads: 16
  resources:
    time=16,
    mem=128
  input:
    pop_mean = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.population_mean_betas.tsv.gz"),
    beta_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.betas.mat.gz"),
    depth_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.coverage.mat.gz")
  params:
    script = "segment_population_mean_betas.R"
  output:
    seg_bed = join(outdir, "methylation_results/Meth_segments.{ref_name}.tissue_{tissue}.segment_betas.bed"),
    varseg_bed = join(outdir, "methylation_results/Meth_segments.{ref_name}.tissue_{tissue}.segment_betas.variable_segments.bed"),
    varseg_plot_file = join(outdir, "methylation_results/Variable_segments.{ref_name}.tissue_{tissue}.segment_CV2.pdf"),
    pca_out = join(outdir, "methylation_results/Gloabal_Methylation_PCA_{ref_name}_tissue_{tissue}/Global_meth_pcs.{ref_name}.tissue_{tissue}.txt")
  conda: "r"
  shell: """
    pca_plot_dir=$(dirname {output.pca_out})
    Rscript {params.script} \
        --population_beta {input.pop_mean} \
        --beta_mat {input.beta_mat} \
        --depth_mat {input.depth_mat} \
        --segment_bed_out {output.seg_bed} \
        --varsegment_bed_out {output.varseg_bed} \
        --varseg_plot_file {output.varseg_plot_file} \
        --pca_plot_dir $pca_plot_dir \
        --global_meth_pcs_out {output.pca_out} 
  """

rule call_outliers:
  threads: 4
  resources:
    time=16,
    mem=146
  input:
    pop_mean = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.population_mean_betas.tsv.gz"),
    beta_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.betas.mat.gz"),
    depth_mat = join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.coverage.mat.gz"),
    global_pcs = join(outdir, "methylation_results/Gloabal_Methylation_PCA_{ref_name}_tissue_{tissue}/Global_meth_pcs.{ref_name}.tissue_{tissue}.txt")
  params:
    script = "call_methylation_outliers.R"
  output:
    outlier_bed = join(outdir,"{sample}/methylation/{sample}.{ref_name}.tissue_{tissue}.METAFORA.outlier_regions.bed"),
    outlier_z_mat = join(outdir, "{sample}/methylation/{sample}.{ref_name}.tissue_{tissue}.METAFORA.outlier_region.zscore.mat")
  conda: 'r'
  shell: """ 
    outlier_plot_dir="$(dirname {output.outlier_bed})/outlier_plots"
    mkdir -p $outlier_plot_dir
    Rscript {params.script} \
        --sample {wildcards.sample} \
        --population_beta {input.pop_mean} \
        --beta_mat {input.beta_mat} \
        --depth_mat {input.depth_mat} \
        --global_meth_pcs {input.global_pcs} \
        --outlier_bed {output.outlier_bed} \
        --outlier_z_mat {output.outlier_z_mat} \
        --plot_dir $outlier_plot_dir
  """

rule combine_outliers:
  threads: 1 
  resources:
     time=1,
     mem=25
  input:
     expand(join(outdir,"{sample}/methylation/{sample}.GRCh38.tissue_{tissue}.METAFORA.outlier_regions.bed"),zip,sample=samples,tissue=tissues)
  output:
     join(outdir, "methylation_results/UDN_Cohort.Methylation_outliers.combined.bed")
  shell: """ 
      echo -e "chrom\tstart\tend\twidth\tstrand\tID\tnum.mark\tseg.mean\tstartRow\tendRow\tseg_id\tzscore " > {output}
      cat {input} | grep -v "^seqnames" >> {output}
  """
