from os.path import join
import pandas as pd
from collections import defaultdict
import glob
outdir = config["output_dir"]
scratch_dir = config["scratch_dir"]
scratch_dir="/tmp/tannerj"
project_dir=config["project_dir"]
samples = list(config["samples"].keys())
tissue_dict = defaultdict(list)
for sample in samples: 
  if sample in ['UDN631726','UDN890454']: continue
  tissue_dict[config["samples"][sample]["ONT"]["tissue"]].append(sample)
for tissue in tissue_dict:
  print(tissue, len(tissue_dict[tissue]))
process_tissues = [ tissue for tissue in tissue_dict.keys() if len(tissue_dict[tissue]) > 30 ]

rerun_samples=list()
for tissue in process_tissues:
  for sample in tissue_dict[tissue]:
    if sample.startswith('UDN'): rerun_samples.append(sample)


rule all:
  input:
#       expand(join(outdir, "{sample}/{ref_name}/methylation_calls/{sample}.{ref_name}.f5c.methylation_calls.tsv.gz"),
#           sample=rerun_samples,
#           ref_name="GRCh38"),
    expand(join(outdir, "methylation_results/Population_methylation.{ref_name}.tissue_{tissue}.population_mean_betas.tsv.gz"),
           ref_name="GRCh38", 
          tissue='Fibroblast')

rule remora_cpg_betas:
  threads: 1
    

rule nanopolish_cpg_betas:
    threads: 1
    resources: 
      time=8, 
      mem=32
    input:
      join(outdir, "{sample}/{ref_name}/methylation_calls/{sample}.{ref_name}.f5c.methylation_calls.tsv.gz"),
    params:
      tmp_dir = join(scratch_dir, "{sample}"),
      tmp_calls = join(scratch_dir, "{sample}/{sample}.{ref_name}.tmp.methylation.tsv"),
      tmp_betas = join(scratch_dir, "{sample}/{sample}.{ref_name}.tmp.cpg_sites.betas.tsv"),
      script = join(project_dir,"scripts/calculate_methylation_frequency.py"),
      LLR_thresh = 1.5
    output:
      join(outdir, "{sample}/{ref_name}/methylation_calls/{sample}.{ref_name}.cpg_sites.betas.tsv")
    shell: """
       mkdir -p {params.tmp_dir}
       zcat {input} > {params.tmp_calls}
       python {params.script} -c {params.LLR_thresh} -s {params.tmp_calls} > {params.tmp_betas} 
       mv {params.tmp_betas} {output}
       rm {params.tmp_calls}
    """

rule create_tissue_sample_reference:
  threads: 1
  resources:
    time=16,
    mem=128
  input:
    lambda wildcards: expand(join(outdir, "{sample}/{{ref_name}}/methylation_calls/{sample}.{{ref_name}}.cpg_sites.betas.tsv"),
                             sample=tissue_dict[wildcards.tissue]),
  params:
    bed = join(project_dir,"methylation_results/GRCh38.CpG_sites.bed"),
    filelist=join(outdir,"tmp.{tissue}.{ref_name}.file_list.txt"),
    script = join(project_dir,"scripts/calculate_population_mean_betas.R"),
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
    sort -k1,1 -k2,2n {params.tmp_beta} | bgzip -c > {output.beta_mat}
    sort -k1,1 -k2,2n {params.tmp_depth} | bgzip -c > {output.depth_mat}
    sort -k1,1 -k2,2n {params.tmp_mean} | bgzip -c > {output.mean}

    tabix -p bed -S 1 {output.beta_mat}
    tabix -p bed -S 1 {output.depth_mat}
    tabix -p bed -S 1 {output.mean}
  """

chrom_list = ["chr{}".format(i) for i in range(1,23)]
rule find_candidate_outlier_regions:
  threads: 1
  resources:
    time=24,
    mem=48
  input:
    betas = join(outdir, "{sample}/{ref_name}/nanopolish_methylation/{sample}.{ref_name}.cpg_sites.betas.tsv"),
    tissue_pop_mean = lambda wildcards: expand(join(outdir,"methylation_results/Population_mean_betas.{{ref_name}}.tissue_{tissue}"), tissue=sample_tissues[wildcards.sample])
  params:
    script = join(project_dir, "scripts/find_cnadidate_outlier_regions.Rscript")
  output:
    join(outdir, "{sample}/{ref_name}/nanopolish_methylation/{sample}.{ref_name}.candidate_methyation_outliers.regions.bed")
  shell: """
    Rscript {params.script} {input.betas} {input.tissue_pop_mean} {output}
  """

rule call_methylation_outliers:
  threads: 1 
  resources:
    time=24,
    mem=48
  input:
    cands = join(outdir, "{sample}/{ref_name}/nanopolish_methylation/{sample}.{ref_name}.candidate_methyation_outliers.regions.bed"),
    beta_matrix= lambda wildcards: expand(join(outdir,"methylation_results/Population_methylation.{{ref_name}}.tissue_{tissue}.betas.mat.gz"), tissue=sample_tissues[wildcards.sample]),
    depth_matrix= lambda wildcards: expand(join(outdir,"methylation_results/Population_methylation.{{ref_name}}.tissue_{tissue}.coverage.mat.gz"), tissue=sample_tissues[wildcards.sample])
  params:
    script = join(project_dir, "scripts/methylation_outlier_calling.Rscript")
  output:
    join(outdir, "{sample}/{ref_name}/nanopolish_methylation/{sample}.{ref_name}.methyation_outliers.regions.tsv")
  shell: """
    Rscript {params.script} {wildcards.sample} {input.cands} {input.beta_matrix} {input.depth_matrix} {output}
  """
