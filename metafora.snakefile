from os.path import join
from collections import defaultdict
import pandas as pd
import glob
outdir = config["output_dir"] 
scratch = config["scratch_dir"]

sample_table = pd.read_table(config["sample_table"])
samples = sample_table['Sample_name'].to_list()

sample_tissues = sample_table['Tissue'].to_list()
tissue_dict=defaultdict(list)

phased_bam = defaultdict(lambda: False)
sample_tissue_map = { row.Sample_name : row.Tissue for i,row in sample_table.iterrows() }
technology_map = { row.Sample_name : row.Technology for i,row in sample_table.iterrows() }
input_map = { row.Sample_name : row.Methylation_input for i,row in sample_table.iterrows()}
if "Phased" in sample_table.columns:
  for i,row in sample_table.iterrows():
    phased_bam[row.Sample_name] = (str(row.Phased).upper() in ["TRUE", "T"])

for sample,tissue in sample_tissue_map.items():
  tissue_dict[tissue].append(sample)

for tissue in tissue_dict.keys():
  print("tissue: %s , N = %d" % (tissue, len(tissue_dict[tissue])))

unique_tissues = list(tissue_dict.keys())

def get_input_samples(wildcards):
  inputs = []
  for s in tissue_dict[wildcards.tissue]:
    if technology_map[s] == "Metafora":
      inputs.append(input_map[s])
    elif technology_map[s] == "PacBio":
      inputs.append(join(outdir,"sample_level_data/"+s+"/"+s+".tech_PacBio.METAFORA_formatted.cpg_methylation.bed.gz"))
    elif technology_map[s] == "ONT":
      inputs.append(join(outdir,"sample_level_data/"+s+"/"+s+".tech_ONT.METAFORA_formatted.cpg_methylation.bed.gz"))
  return inputs

#set up default params
if not "params" in config:
  config["params"] = {}
MAX_DEPTH = config["params"]["MAX_DEPTH"] if "MAX_DEPTH" in config["params"] else 30
MIN_SEG_SIZE = config["params"]["MIN_SEG_SIZE"] if "MIN_SEG_SIZE" in config["params"] else 20
MIN_ABS_ZSCORE = config["params"]["MIN_ABS_ZSCORE"] if "MIN_ABS_ZSCORE" in config["params"] else 3
MIN_ABS_DELTA = config["params"]["MIN_ABS_DELTA"] if "MIN_ABS_DELTA" in config["params"] else 0.25
SKIP_SEX_CHROMOSOME_ESTIMATION = config["params"]["SKIP_SEX_CHROMOSOME_ESTIMATION"] if "SKIP_SEX_CHROMOSOME_ESTIMATION" in config["params"] else "FALSE"

reference_seqnames_table = config["reference_seqnames_table"] if "reference_seqnames_table" in config else "./GRCh38_ref_seqnames_table.txt"
autosomes = []
sex_chroms = []
chrX_seqname = "SKIP"
chrY_seqname = "SKIP"

for line in open(reference_seqnames_table, 'r'):
  seqname,chrom_type = line.strip().split()
  if chrom_type == "autosome": autosomes.append(seqname)
  elif chrom_type == "sex_chromosome_X": 
    sex_chroms.append(seqname)
    chrX_seqname = seqname
  elif chrom_type == "sex_chromosome_Y":
    sex_chroms.append(seqname)
    chrY_seqname = seqname
  else:
    print("unregnozied chromosome type (%s) in reference_seqnames_table, must be one of [autosome, sex_chromosome_X, sex_chromosome_Y]")
    raise SystemExit

if len(sex_chroms) != 2 or chrX_seqname == "SKIP" or chrY_seqname == "SKIP":
  print("User did not correctly specify X and Y chromosome seqnames, so sex chromosome copy number estimation and outlier calling will be skipped")
  SKIP_SEX_CHROMOSOME_ESTIMATION = "TRUE"
  sex_chroms = []

if SKIP_SEX_CHROMOSOME_ESTIMATION in ["TRUE","T","True","true",True]:
  print("Skipping sex chromosome copy number estimation and outlier calling")
  SKIP_SEX_CHROMOSOME_ESTIMATION = "TRUE"
  sex_chroms=[]

rule all:
    input:
      #expand(join(outdir, "sample_level_data/{sample}/{sample}.chrX_inactivation_skew.summary_dat.txt"), sample=samples),
      expand(join(outdir, "sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_report.html"), zip, sample=samples, tissue=sample_tissues),
      join(outdir, "summary_figures/METAFORA.outlier_count_per_sample_tissue.tsv")

rule create_cpg_reference:
    threads: 16
    resources: 
      mem=24,
      time=12
    input:
      ref = config["reference_fasta"]
    params:
      script = "scripts/create_reference_cpg_bed.R",
      valid_chroms = ','.join(autosomes + sex_chroms)
    output:
      cpg_bed = join(outdir, "cpg_reference.bed")
    conda: 'envs/metafora.yaml'
    shell: """
      Rscript {params.script} \
          --reference {input.ref} \
          --valid_chroms {params.valid_chroms} \
          --cpg_bed_out {output.cpg_bed}
    """


rule modBam2Bed:
    threads: 4
    resources:
        mem=32,
        time=12
    input:
        bam = lambda w: input_map[w.sample],
        ref = config["reference_fasta"]
    output:
        meth_bed = join(outdir, "sample_level_data/{sample}/{sample}.tech_ONT.cpg_methylation.bed.gz"),
        tbi = join(outdir, "sample_level_data/{sample}/{sample}.tech_ONT.cpg_methylation.bed.gz.tbi")
    conda: 'envs/methylation.yaml'
    shell: """
        modbam2bed -t {threads} --cpg --combine -m 5mC -d 40 {input.ref} {input.bam} | bgzip > {output.meth_bed}
        tabix {output.meth_bed}
    """

rule format_modBam2Bed:
    threads: 16
    resources:
      time=4,
      mem=128
    input:
      join(outdir, "sample_level_data/{sample}/{sample}.tech_ONT.cpg_methylation.bed.gz")
    params:
      script = "scripts/format_modBam2Bed.R",
      meth_bed = join(outdir, "sample_level_data/{sample}/{sample}.tech_ONT.METAFORA_formatted.cpg_methylation.bed")
    output:
      meth_bed = join(outdir, "sample_level_data/{sample}/{sample}.tech_ONT.METAFORA_formatted.cpg_methylation.bed.gz"),
      tbi = join(outdir, "sample_level_data/{sample}/{sample}.tech_ONT.METAFORA_formatted.cpg_methylation.bed.gz.tbi")
    conda: 'envs/metafora.yaml'
    shell: """
      Rscript {params.script} --input {input} --output {params.meth_bed}
      bgzip {params.meth_bed}
      tabix -p bed -S 1 {output.meth_bed}
    """

rule modBam2Bed_HP:
    threads: 4
    resources:
        mem=32,
        time=12
    input:
        bam = lambda w: input_map[w.sample],
        ref = config["reference_fasta"]
    output:
        meth_bed = join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_{hp}.tech_ONT.cpg_methylation.bed.gz"),
        tbi = join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_{hp}.tech_ONT.cpg_methylation.bed.gz.tbi")
    conda: 'envs/methylation.yaml'
    shell: """
        modbam2bed -t {threads} --cpg --combine -m 5mC -d 40 --haplotype={wildcards.hp} {input.ref} {input.bam} | bgzip > {output.meth_bed}
        tabix {output.meth_bed}
    """


rule format_modBam2Bed_HP:
    threads: 16
    resources:
      time=4,
      mem=128
    input:
      join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_{hp}.tech_ONT.cpg_methylation.bed.gz")
    params:
      script = "scripts/format_modBam2Bed.R",
      meth_bed = join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_{hp}.tech_ONT.METAFORA_formatted.cpg_methylation.bed")
    output:
      meth_bed = join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_{hp}.tech_ONT.METAFORA_formatted.cpg_methylation.bed.gz"),
      tbi = join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_{hp}.tech_ONT.METAFORA_formatted.cpg_methylation.bed.gz.tbi")
    conda: 'envs/metafora.yaml'
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
    bam = lambda w: input_map[w.sample],
    ref = lambda w: config["reference_fasta"]
  params:
    prefix = join(outdir, "sample_level_data/{sample}/{sample}.tech_PacBio.cpg_methylation"),
  output:
    join(outdir, "sample_level_data/{sample}/{sample}.tech_PacBio.cpg_methylation.combined.bed.gz"),
    #lambda w:  "" if not phased_bam[w.sample] else ["sample_level_data/{sample}/{sample}.tech_PacBio.cpg_methylation.hap1.bed.gz", "sample_level_data/{sample}/{sample}.tech_PacBio.cpg_methylation.hap2.bed.gz"]
  conda: "envs/pb_cpg_tools.yaml"
  shell: """
    aligned_bam_to_cpg_scores \
        --bam {input.bam} \
        --output-prefix {params.prefix} \
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
    join(outdir, "sample_level_data/{sample}/{sample}.tech_PacBio.cpg_methylation.combined.bed.gz")
  params:
    tmp_bed = join(outdir, "sample_level_data/{sample}/{sample}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed")
  output:
    meth_bed = join(outdir, "sample_level_data/{sample}/{sample}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed.gz"),
    tbi = join(outdir, "sample_level_data/{sample}/{sample}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed.gz.tbi")
  conda: 'envs/methylation.yaml'
  shell: """ 
    echo -e "chromosome\tstart\tend\tdepth\tbeta" > {params.tmp_bed}
    zcat {input} | grep -v "^#" | awk '{{print $1,$2,$3,$6,($9/100)}}' | sed 's/ /\t/g' >> {params.tmp_bed} 

    bgzip {params.tmp_bed} 
    tabix -p bed -S 1 {output.meth_bed}
  """

rule format_pacbio_hap:
  threads: 1 
  resources:
    time=4,
    mem=24
  input:
    join(outdir, "sample_level_data/{sample}/{sample}.tech_PacBio.cpg_methylation.combined.bed.gz")
  params:
    hp_bed = join(outdir, "sample_level_data/{sample}/{sample}.tech_PacBio.cpg_methylation.hap{hp}.bed.gz"),
    tmp_bed = join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_{hp}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed")
  output:
    meth_bed = join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_{hp}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed.gz"),
    tbi = join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_{hp}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed.gz.tbi")
  conda: 'envs/methylation.yaml'
  shell: """ 
    echo -e "chromosome\tstart\tend\tdepth\tbeta" > {params.tmp_bed}
    zcat {params.hp_bed} | grep -v "^#" | awk '{{print $1,$2,$3,$6,($9/100)}}' | sed 's/ /\t/g' >> {params.tmp_bed} 

    bgzip {params.tmp_bed} 
    tabix -p bed -S 1 {output.meth_bed}
  """
  
rule create_tissue_sample_reference:
  threads: 1
  resources:
    time=24,
    mem=128
  input:
    meth_beds = get_input_samples,
    cpg_bed = join(outdir, "cpg_reference.bed")
  params:
    filelist=join(outdir,"tmp.{tissue}.file_list.chrom_{chr}.txt"),
    script = "scripts/calculate_population_mean_betas.R",
    tmp_beta = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.betas.mat"),
    tmp_depth = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.coverage.mat"),
    tmp_mean = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv"),
    min_segment_cpgs = MIN_SEG_SIZE
  output:
    beta_mat = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.betas.mat.gz"),
    beta_tbi = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.betas.mat.gz.tbi"),
    depth_mat = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.coverage.mat.gz"),
    depth_tbi = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.coverage.mat.gz.tbi"),
    mean = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv.gz"),
    mean_tbi = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv.gz.tbi"),
    seg_beta = join(outdir,"Population_methylation.tissue_{tissue}/Meth_segments.tissue_{tissue}.segment_betas.chrom_{chr}.bed"),
    seg_depth = join(outdir,"Population_methylation.tissue_{tissue}/Meth_segments.tissue_{tissue}.segment_coverage.chrom_{chr}.mat")
  conda: "envs/metafora.yaml"
  shell: """
    ls {input.meth_beds} > {params.filelist}
    Rscript {params.script} \
        --filelist {params.filelist} \
        --chrom {wildcards.chr} \
        --cpgs {input.cpg_bed} \
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

#rule append_samples_to_reference:
#  threads: 1
#  resources:
#    time=24,
#    mem=128 
#  input: 
#    meth_beds = get_input_samples,
#    reference_beta_mat = config["sample_reference"]["beta_mat"],
#    reference_depth_mat = config["sample_reference"]["depth_mat"]
#  output:
#    beta_mat = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.betas.mat.gz"),
#    beta_tbi = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.betas.mat.gz.tbi"),
#    depth_mat = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.coverage.mat.gz"),
#    depth_tbi = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.coverage.mat.gz.tbi"),
#    mean = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv.gz"),
#    mean_tbi = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.population_mean_betas.tsv.gz.tbi"),
#    seg_beta = join(outdir,"Population_methylation.tissue_{tissue}/Meth_segments.tissue_{tissue}.segment_betas.chrom_{chr}.bed"),
#    seg_depth = join(outdir,"Population_methylation.tissue_{tissue}/Meth_segments.tissue_{tissue}.segment_coverage.chrom_{chr}.mat")

rule compute_hidden_factors:
  threads: 16
  resources:
    time=12,
    mem=128
  input:
    seg_beta = expand(join(outdir,"Population_methylation.tissue_{{tissue}}/Meth_segments.tissue_{{tissue}}.segment_betas.chrom_{chr}.bed"),chr=autosomes + sex_chroms),
    seg_depth = expand(join(outdir,"Population_methylation.tissue_{{tissue}}/Meth_segments.tissue_{{tissue}}.segment_coverage.chrom_{chr}.mat"),chr=autosomes + sex_chroms)
  params:
    script = "scripts/compute_hidden_factors.R",
    seg_beta = ','.join(list(expand(join(outdir,"Population_methylation.tissue_{{tissue}}/Meth_segments.tissue_{{tissue}}.segment_betas.chrom_{chr}.bed"),chr=autosomes))),
    seg_depth = ','.join(list(expand(join(outdir,"Population_methylation.tissue_{{tissue}}/Meth_segments.tissue_{{tissue}}.segment_coverage.chrom_{chr}.mat"),chr=autosomes))),
    sex_seg_beta = "SKIP" if SKIP_SEX_CHROMOSOME_ESTIMATION == "TRUE" else ','.join(list(expand(join(outdir,"Population_methylation.tissue_{{tissue}}/Meth_segments.tissue_{{tissue}}.segment_betas.chrom_{chr}.bed"),chr=sex_chroms))),
    sex_seg_depth = "SKIP" if SKIP_SEX_CHROMOSOME_ESTIMATION == "TRUE" else ','.join(list(expand(join(outdir,"Population_methylation.tissue_{{tissue}}/Meth_segments.tissue_{{tissue}}.segment_coverage.chrom_{chr}.mat"),chr=sex_chroms))),
    chrX_seqname =  chrX_seqname,
    chrY_seqname =  chrY_seqname,
    covariates = config["covariates"] if "covariates" in config else "SKIP",
    plot_out_dir = join(outdir, "Global_Methylation_PCA_tissue_{tissue}/")
  output:
    global_pcs = join(outdir, "Global_Methylation_PCA_tissue_{tissue}/PCA_covariates.txt"),
    cor_summary_table = join(outdir, "Global_Methylation_PCA_tissue_{tissue}/Mean_pairwise_correlation_summary.txt"),
    sex_summary_table = join(outdir, "Global_Methylation_PCA_tissue_{tissue}/Sex_chromosome_estimates_summary.txt")
  conda: "envs/metafora.yaml"
  shell: """
    Rscript {params.script} \
        --seg_beta {params.seg_beta} \
        --seg_depth {params.seg_depth} \
        --sex_seg_beta {params.sex_seg_beta} \
        --sex_seg_depth {params.sex_seg_depth} \
        --correlation_summary_out {output.cor_summary_table} \
        --sex_chrom_summary_out {output.sex_summary_table} \
        --global_meth_pcs_out {output.global_pcs} \
        --plot_out_dir {params.plot_out_dir} \
        --chrX_seqname {params.chrX_seqname} \
        --chrY_seqname {params.chrY_seqname} \
        --covariates {params.covariates} 
  """

rule call_outliers:
  threads: 4
  resources:
    time=16,
    mem=64
  input:
    beta_mat = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.betas.mat.gz"),
    depth_mat = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.coverage.mat.gz"),
    global_pcs = join(outdir, "Global_Methylation_PCA_tissue_{tissue}/PCA_covariates.txt")
  params:
    script = "scripts/call_methylation_outliers.R",
    MAX_DEPTH = MAX_DEPTH, 
    MIN_SEG_SIZE = MIN_SEG_SIZE, 
    MIN_ABS_ZSCORE = MIN_ABS_ZSCORE, 
    MIN_ABS_DELTA = MIN_ABS_DELTA 
  output:
    outlier_bed = temp(join(outdir,"sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_regions.chrom_{chr}.bed")),
    outlier_z_mat = temp(join(outdir, "sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_regions.zscore.chrom_{chr}.mat"))
  conda: 'envs/metafora.yaml'
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
        --max_depth {params.MAX_DEPTH} \
        --tissue {wildcards.tissue}
  """

rule call_outliers_sex_chroms:
  threads: 4
  resources:
    time=16,
    mem=128
  input:
    beta_mat = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.betas.mat.gz"),
    depth_mat = join(outdir, "Population_methylation.tissue_{tissue}/Population_methylation.tissue_{tissue}.chrom_{chr}.coverage.mat.gz"),
    global_pcs = join(outdir, "Global_Methylation_PCA_tissue_{tissue}/PCA_covariates.txt"),
    sex_chrom_df = join(outdir, "Global_Methylation_PCA_tissue_{tissue}/Sex_chromosome_estimates_summary.txt")
  params:
    script = "scripts/call_methylation_outliers.sex_chromosomes.R",
    chrY_seqname = chrY_seqname,
    MAX_DEPTH = MAX_DEPTH, 
    MIN_SEG_SIZE = MIN_SEG_SIZE, 
    MIN_ABS_ZSCORE = MIN_ABS_ZSCORE, 
    MIN_ABS_DELTA = MIN_ABS_DELTA 
  output:
    outlier_bed = temp(join(outdir,"sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_regions.sex_chroms.chrom_{chr}.bed")),
    outlier_z_mat = temp(join(outdir, "sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_regions.sex_chroms.zscore.chrom_{chr}.mat"))
  conda: 'envs/metafora.yaml'
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
        --chrY_seqname {params.chrY_seqname} \
        --outlier_bed {output.outlier_bed} \
        --outlier_z_mat {output.outlier_z_mat} \
        --plot_dir $outlier_plot_dir \
        --min_seg_size {params.MIN_SEG_SIZE} \
        --min_abs_zscore {params.MIN_ABS_ZSCORE} \
        --min_abs_delta {params.MIN_ABS_DELTA} \
        --max_depth {params.MAX_DEPTH} \
        --tissue {wildcards.tissue} 
  """

rule combine_chrom_outliers:
  threads: 1 
  resources:
    time=4,
    mem=12
  input:
    beds = expand(join(outdir, "sample_level_data/{{sample}}/{{sample}}.tissue_{{tissue}}.METAFORA.outlier_regions.chrom_{chr}.bed"), chr=autosomes) +
           expand(join(outdir, "sample_level_data/{{sample}}/{{sample}}.tissue_{{tissue}}.METAFORA.outlier_regions.sex_chroms.chrom_{chr}.bed"), chr=sex_chroms),
    mats = expand(join(outdir, "sample_level_data/{{sample}}/{{sample}}.tissue_{{tissue}}.METAFORA.outlier_regions.zscore.chrom_{chr}.mat"), chr=autosomes) +
           expand(join(outdir, "sample_level_data/{{sample}}/{{sample}}.tissue_{{tissue}}.METAFORA.outlier_regions.sex_chroms.zscore.chrom_{chr}.mat"), chr=sex_chroms)
  output:
    outlier_bed = join(outdir,"sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_regions.bed"),
    outlier_z_mat = join(outdir, "sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_regions.zscore.mat")
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
    lambda w: expand(join(outdir,"sample_level_data/{sample}/{sample}.tissue_{{tissue}}.METAFORA.outlier_regions.bed"),
                     sample=tissue_dict[w.tissue]),
  output:
    join(outdir, "METAFORA.tissue_{tissue}.methylation_outliers.combined.tsv")
  shell: """ 
      header=$(for i in $(ls {input}); do head -1 $i; done | sort | tail -1)
      echo "$header" > {output}
      cat {input} | grep -v "$header" | grep -v "^$" >> {output} || true
  """

rule summary_plots:
  threads: 1 
  resources:
    time=2,
    mem=24
  input:
    outlier_bed = expand(join(outdir, "METAFORA.tissue_{tissue}.methylation_outliers.combined.tsv"), tissue=unique_tissues),
    covariates = expand(join(outdir, "Global_Methylation_PCA_tissue_{tissue}/PCA_covariates.txt"),tissue=unique_tissues)
  params:
    script = "scripts/plot_summary_figures.R",
    plot_dir = join(outdir, "summary_figures"),
    outlier_files = ','.join(expand(join(outdir, "METAFORA.tissue_{tissue}.methylation_outliers.combined.tsv"), tissue=unique_tissues)),
    covariates_files = ','.join(expand(join(outdir, "Global_Methylation_PCA_tissue_{tissue}/PCA_covariates.txt"),tissue=unique_tissues)),
    MIN_ABS_ZSCORE = MIN_ABS_ZSCORE, 
    MIN_ABS_DELTA = MIN_ABS_DELTA 

  output:
    join(outdir, "summary_figures/METAFORA.outlier_count_per_sample_tissue.tsv")
  conda: "envs/metafora.yaml"
  shell: """
    Rscript {params.script} \
      --outlier_files {params.outlier_files} \
      --covariates_files {params.covariates_files} \
      --plot_dir_out {params.plot_dir} \
      --summary_out {output} \
      --min_abs_delta {params.MIN_ABS_DELTA} \
      --min_abs_zscore {params.MIN_ABS_ZSCORE} 
  """

rule annotate_haplotype_delta:
  threads: 1 
  resources:
    time=4,
    mem=48
  input:
    outlier_bed = join(outdir,"sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_regions.bed"),
    combined_bed = lambda w: join(outdir, "sample_level_data/{sample}/{sample}.tech_" + technology_map[w.sample] + ".METAFORA_formatted.cpg_methylation.bed.gz"),
    hap1_bed = lambda w: join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_1.tech_" + technology_map[w.sample] + ".METAFORA_formatted.cpg_methylation.bed.gz"),
    hap2_bed = lambda w: join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_2.tech_" + technology_map[w.sample] + ".METAFORA_formatted.cpg_methylation.bed.gz")
  output:
    annotated_bed = join(outdir, "sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_regions.haplotype_annotated.bed")
  conda: 'envs/metafora.yaml'
  shell: """
    Rscript scripts/annotate_haplotype_methylation.R \
        --outlier_bed {input.outlier_bed} \
        --combined {input.combined_bed} \
        --hap1 {input.hap1_bed} \
        --hap2 {input.hap2_bed} \
        --annotated_out {output.annotated_bed}
  """

rule calculate_X_skew:
  threads: 1 
  resources:
    time=4,
    mem=48
  input:
    combined_bed = lambda w: join(outdir, "sample_level_data/{sample}/{sample}.tech_" + technology_map[w.sample] + ".METAFORA_formatted.cpg_methylation.bed.gz"),
    hap1_bed = lambda w: join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_1.tech_" + technology_map[w.sample] + ".METAFORA_formatted.cpg_methylation.bed.gz"),
    hap2_bed = lambda w: join(outdir, "sample_level_data/{sample}/{sample}.Haplotype_2.tech_" + technology_map[w.sample] + ".METAFORA_formatted.cpg_methylation.bed.gz")
  params:
    X_bed = config["chrX_inactivation_regions"]
  output:
    join(outdir, "sample_level_data/{sample}/{sample}.chrX_inactivation_skew.summary_dat.txt")
  conda: 'envs/metafora.yaml'
  shell: """
    Rscript scripts/calculate_X_skew.R \
        --sample {wildcards.sample} \
        --X_region_bed {params.X_bed} \
        --combined {input.combined_bed} \
        --hap1 {input.hap1_bed} \
        --hap2 {input.hap2_bed} \
        --out_tsv {output}
  """

rule make_outlier_report:
  threads: 1 
  resources:
    time=4,
    mem=12
  input:
    outlier_bed = lambda w: join(outdir,"sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_regions" + (".haplotype_annotated.bed" if phased_bam[w.sample] else ".bed")),
    outlier_z_mat = join(outdir, "sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_regions.zscore.mat"),
    covariates = join(outdir, "Global_Methylation_PCA_tissue_{tissue}/PCA_covariates.txt")
  params:
    annos = config["annotation_track_tsv"],
    sample_table = config['sample_table'],
    tmp_bed_out = join(outdir, "sample_level_data/{sample}/tmp.outliers.bed"),
    tmp_tsv_out = join(outdir, "sample_level_data/{sample}/tmp.outliers.tsv"),
    tmp_config_json = join(outdir, "sample_level_data/{sample}/tmp.track_config.json")
  output:
    report = join(outdir, "sample_level_data/{sample}/{sample}.tissue_{tissue}.METAFORA.outlier_report.html")
  conda: 'envs/igv_report.yaml'
  shell: """
        Rscript scripts/create_report_json.R \
            --sample {wildcards.sample} \
            --outlier_bed {input.outlier_bed} \
            --outlier_z_mat {input.outlier_z_mat} \
            --covariates {input.covariates} \
            --sample_table {params.sample_table} \
            --annos {params.annos} \
            --output_bed {params.tmp_bed_out} \
            --output_tsv {params.tmp_tsv_out} \
            --output_json {params.tmp_config_json}
        
        create_report {params.tmp_tsv_out} \
            --genome hg38 \
            --flanking 2000 \
            --track-config {params.tmp_config_json} \
            --info-columns seg_id coordinates num.mark pop_median delta zscore combined_depth haplotype_coverage_bias hap_delta Tissue \
            --sequence 1 --begin 2 --end 3 \
            --output {output.report}

        rm -f {params.tmp_bed_out} {params.tmp_tsv_out} {params.tmp_config_json}
  """
