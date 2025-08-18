version 1.0

workflow metafora {
    input {
        File sample_table
        File ref_fa
        File cpg_reference_script
        File create_tissue_reference_script
        File compute_hidden_factors_script
        File call_outliers_script
        File call_sex_outliers_script
        File annotate_haplotype_delta_script
        String tissue
        File covariates
    }

    call CreateCpgReference {
        input:
            ref_fa = ref_fa,
            script = cpg_reference_script
    }

    call ParseSampleTable {
        input:
            sample_table = sample_table
    }

    Array[String] autosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
    Array[String] sex_chromosomes = ["chrX","chrY"]

    scatter (sample in ParseSampleTable.samples) {
        call PacbioCpgToolsHap {
            input:
                sample=sample.name,
                bam=sample.input_file,
                index=sample.input_file_index,
                ref_fa=ref_fa
            }
        call FormatPacBioHap as FormatPacBioHapH1 {
            input:
                haplotype_bed = PacbioCpgToolsHap.pacbio_sample_meth_bed_H1,
                sample = sample.name,
                hp = "1"
            }
        call FormatPacBioHap as FormatPacBioHapH2 {
            input:
                haplotype_bed = PacbioCpgToolsHap.pacbio_sample_meth_bed_H2,
                sample = sample.name,
                hp = "2"
            }
        call FormatPacBio as FormatPacBioCombined{
            input:
                sample = sample.name,
                bed = PacbioCpgToolsHap.pacbio_sample_meth_bed
            }
        SampleHaplotypes pacbio_haplotype_data = {
            "name" : sample.name,
            "tissue": sample.tissue,
            "haplotype_bed_combined": FormatPacBioCombined.meth_bed,
            "haplotype_bed_combined_index": FormatPacBioCombined.meth_bed_tbi,
            "haplotype_bed_H1" : FormatPacBioHapH1.meth_bed,
            "haplotype_bed_H1_index" : FormatPacBioHapH1.meth_bed_tbi,
            "haplotype_bed_H2" : FormatPacBioHapH2.meth_bed,
            "haplotype_bed_H2_index" : FormatPacBioHapH2.meth_bed_tbi
            }

        }

    output {
        Array[SampleHaplotypes] beds = pacbio_haplotype_data
    }
}

struct Sample {
        String name
        String tissue
        String technology
        Boolean phased
        String input_file
        String input_file_index
    }

struct SampleHaplotypes {
    String name
    String tissue
    File haplotype_bed_combined
    File haplotype_bed_combined_index
    File haplotype_bed_H1
    File haplotype_bed_H1_index
    File haplotype_bed_H2
    File haplotype_bed_H2_index
}

struct ChromosomeReference {
    String chromosome
    File beta_matrix
    File beta_matrix_tbi
    File depth_matrix
    File depth_matrix_tbi
}

#struct SampleHaplotypeAndOutliers

task ParseSampleTable {
    input {
        File sample_table
    }

    command <<<
    awk 'BEGIN {
    FS="\t";
    print "["
}
{
    if (NR > 1) {
        printf "  {\n"
        printf "    \"name\": \"%s\",\n", $1
        printf "    \"tissue\": \"%s\",\n", $2
        printf "    \"technology\": \"%s\",\n", $3
        printf "    \"phased\": %s,\n", ($4 == "TRUE" ? "true" : "false")
        printf "    \"input_file\": \"%s\",\n", $5
        printf "    \"input_file_index\": \"%s\"\n", $6
        printf "  }"
        if (NR < total_lines) {
            printf ",\n"
        } else {
            printf "\n"
        }
    }
}
END {
    print "]"
}' total_lines=$(wc -l < ~{sample_table}) ~{sample_table} > samples.json
    cat samples.json
    >>>

    output {
        Array[Sample] samples = read_json("samples.json")
    }
}

task CreateCpgReference {
    input {
        File ref_fa
        File script
    }

    command <<<
      which Rscript
      Rscript ~{script} \
          --reference ~{ref_fa} \
          --valid_chroms chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
          --cpg_bed_out cpg_reference.bed
    >>>

    runtime {
        docker: "quay.io/jonnguye/metafora:1.4"
        memory: "24G"
        cpu: 16
        disks: "local-disk 128 SSD"
    }

    output {
        File cpg_bed = "cpg_reference.bed"
    }
}

task ModBamToBed {
        input {
                String sample
                File bam
                File index
                File ref_fa
            }
        
        command <<<
            modbam2bed -t 4 --cpg --combine -m 5mC -d 40 ~{ref_fa} ~{bam} | bgzip > ~{sample}.tech_ONT.cpg_methylation.bed.gz 
        tabix ~{sample}.tech_ONT.cpg_methylation.bed.gz
        >>>
        
        runtime {
                #docker: "quay.io/jonnguye/metafora_methylation:1.0"
                memory: "32G"
                cpu: 4
                disks: "local-disk 128 SSD"
            }

        output {
                File meth_bed = "~{sample}.tech_ONT.cpg_methylation.bed.gz"
                File meth_bed_tbi = "~{sample}.tech_ONT.cpg_methylation.bed.gz.tbi"
            }
    }

task PacbioCpgTools {
        input {
                String sample
                File bam
                File index
                File ref_fa
            }

        command <<<
        aligned_bam_to_cpg_scores \
        --bam ~{bam} \
        --output-prefix ~{sample}.tech_PacBio.cpg_methylation \
        --modsites-mode reference \
        --ref ~{ref_fa} \
        --threads 16
        >>>
        
        runtime {
                docker: "quay.io/pacbio/pb-cpg-tools:3.0.0_build1"
                memory: "32G"
                cpu: 4
                disks: "local-disk 128 SSD"
            }

        output {
            File pacbio_sample_meth_bed = "~{sample}.tech_PacBio.cpg_methylation.combined.bed.gz"
            }
} 

task PacbioCpgToolsHap {
        input {
                String sample
                File bam
                File index
                File ref_fa
            }

        command <<<
        aligned_bam_to_cpg_scores \
        --bam ~{bam} \
        --output-prefix ~{sample}.tech_PacBio.cpg_methylation \
        --modsites-mode reference \
        --ref ~{ref_fa} \
        --threads 16
        >>>
        
        runtime {
                docker: "quay.io/pacbio/pb-cpg-tools:3.0.0_build1"
                memory: "32G"
                cpu: 4
                disks: "local-disk 128 SSD"
            }

        output {
            File pacbio_sample_meth_bed = "~{sample}.tech_PacBio.cpg_methylation.combined.bed.gz"
            File pacbio_sample_meth_bed_H1 = "~{sample}.tech_PacBio.cpg_methylation.hap1.bed.gz"
            File pacbio_sample_meth_bed_H2 = "~{sample}.tech_PacBio.cpg_methylation.hap2.bed.gz"
            }
} 

task FormatPacBio{
        input {
                File bed
                String sample
            }
        command <<<
        echo -e "chromosome\tstart\tend\tdepth\tbeta" > "~{sample}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed"
        zcat ~{bed} | grep -v "^#" | awk '{{print $1,$2,$3,$6,($9/100)}}' | sed 's/ /\t/g' >> "~{sample}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed"

        bgzip "~{sample}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed"
        tabix -p bed -S 1 "~{sample}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed.gz"
        >>>

        runtime {
            docker: "quay.io/jonnguye/metafora_methylation:1.0"
            memory: "24G"
            cpu: 1
            disks: "local-disk 128 SSD"
            }

        output {
            File meth_bed = "~{sample}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed.gz"
            File meth_bed_tbi = "~{sample}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed.gz.tbi"
        }
    }

task FormatPacBioHap{
        input {
                File haplotype_bed
                String sample
                String hp
            }
        command <<<
            echo -e "chromosome\tstart\tend\tdepth\tbeta" > "~{sample}.Haplotype_~{hp}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed"
            zcat ~{haplotype_bed} | grep -v "^#" | awk '{{print $1,$2,$3,$6,($9/100)}}' | sed 's/ /\t/g' >> "~{sample}.Haplotype_~{hp}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed"

            bgzip "~{sample}.Haplotype_~{hp}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed"
            tabix -p bed -S 1 "~{sample}.Haplotype_~{hp}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed.gz"
        >>>
        
        runtime {
            docker: "quay.io/jonnguye/metafora_methylation:1.0"
            memory: "24G"
            cpu: 1
            disks: "local-disk 128 SSD"
            }

        output {
            File meth_bed = "~{sample}.Haplotype_~{hp}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed.gz"
            File meth_bed_tbi = "~{sample}.Haplotype_~{hp}.tech_PacBio.METAFORA_formatted.cpg_methylation.bed.gz.tbi"
        }
    }

task CreateTissueSampleReference {
        input {
                Array[File] meth_beds
                Array[File] meth_beds_index
                File cpg_bed
                String chr
                String tissue
                File rscript
                Int min_segment_size = 20
            }

        command <<<

        for file in ~{sep=' ' meth_beds}; do
            echo $file >> tmp.~{tissue}.file_list.chrom_~{chr}.txt
        done

        Rscript ~{rscript} \
            --filelist tmp.~{tissue}.file_list.chrom_~{chr}.txt \
            --chrom ~{chr} \
            --cpgs ~{cpg_bed} \
            --min_segment_cpgs ~{min_segment_size} \
            --beta_mat Population_methylation.tissue_~{tissue}.chrom_~{chr}.betas.mat \
            --depth_mat Population_methylation.tissue_~{tissue}.chrom_~{chr}.coverage.mat \
            --pop_mean Population_methylation.tissue_~{tissue}.chrom_~{chr}.population_mean_betas.tsv \
            --segment_beta Meth_segments.tissue_~{tissue}.segment_betas.chrom_~{chr}.bed \
            --segment_depth Meth_segments.tissue_~{tissue}.segment_coverage.chrom_~{chr}.mat

        bgzip Population_methylation.tissue_~{tissue}.chrom_~{chr}.betas.mat
        bgzip Population_methylation.tissue_~{tissue}.chrom_~{chr}.coverage.mat
        bgzip Population_methylation.tissue_~{tissue}.chrom_~{chr}.population_mean_betas.tsv

        tabix -p bed -S 1 Population_methylation.tissue_~{tissue}.chrom_~{chr}.betas.mat.gz
        tabix -p bed -S 1 Population_methylation.tissue_~{tissue}.chrom_~{chr}.coverage.mat.gz
        tabix -p bed -S 1 Population_methylation.tissue_~{tissue}.chrom_~{chr}.population_mean_betas.tsv.gz

        >>>
        
        runtime {
                docker: "quay.io/jonnguye/metafora:1.4"
                memory: "128G"
                cpu: 1
                disks: "local-disk 128 SSD"
            }

        output {
                String chromosome = chr
                File beta_matrix = "Population_methylation.tissue_${tissue}.chrom_${chr}.betas.mat.gz"
                File beta_matrix_tbi = "Population_methylation.tissue_${tissue}.chrom_${chr}.betas.mat.gz.tbi"
                File depth_matrix = "Population_methylation.tissue_${tissue}.chrom_${chr}.coverage.mat.gz" 
                File depth_matrix_tbi = "Population_methylation.tissue_${tissue}.chrom_${chr}.coverage.mat.gz.tbi"
                File mean = "Population_methylation.tissue_${tissue}.chrom_${chr}.population_mean_betas.tsv.gz"
                File mean_tbi = "Population_methylation.tissue_${tissue}.chrom_${chr}.population_mean_betas.tsv.gz.tbi"
                File seg_beta = "Meth_segments.tissue_${tissue}.segment_betas.chrom_${chr}.bed"
                File seg_depth = "Meth_segments.tissue_${tissue}.segment_coverage.chrom_${chr}.mat"
            }
    }

task ComputeHiddenFactors{
        input {
            File rscript
            Array[File] segment_betas
            Array[File] segment_depths
            Array[File] sex_segment_betas
            Array[File] sex_segment_depths
            File covariates
            String tissue
            }

        command <<<
        mkdir Global_Methylation_PCA_tissue_~{tissue}
            Rscript ~{rscript} \
        --seg_beta ~{sep=',' segment_betas} \
        --seg_depth ~{sep=',' segment_depths} \
        --sex_seg_beta ~{sep=',' sex_segment_betas} \
        --sex_seg_depth ~{sep=',' sex_segment_depths} \
        --correlation_summary_out Global_Methylation_PCA_tissue_~{tissue}/Mean_pairwise_correlation_summary.txt \
        --sex_chrom_summary_out Global_Methylation_PCA_tissue_~{tissue}/Sex_chromosome_estimates_summary.txt \
        --global_meth_pcs_out Global_Methylation_PCA_tissue_~{tissue}/PCA_covariates.txt \
        --plot_out_dir Global_Methylation_PCA_tissue_~{tissue}/ \
        --chrX_seqname chrX \
        --chrY_seqname chrY \
        --covariates ~{covariates}
        >>>

        runtime {
            docker: "quay.io/jonnguye/metafora:1.4"
            memory: "128GB"
            cpu: 16
            disks: "local-disk 128 SSD"
            }
        output {
            File global_pcs = "Global_Methylation_PCA_tissue_~{tissue}/PCA_covariates.txt"
            File cor_summary_table = "Global_Methylation_PCA_tissue_~{tissue}/Mean_pairwise_correlation_summary.txt"
            File sex_summary_table = "Global_Methylation_PCA_tissue_~{tissue}/Sex_chromosome_estimates_summary.txt"
            }
    }

task CallOutliers {
        input {
                File rscript
                String sample
                String tissue
                String chrom
                File beta_matrix
                File beta_matrix_tbi
                File depth_matrix
                File depth_matrix_tbi
                File global_pcs
                Int max_depth
                Int min_segment_size
                Int min_abs_zscore
                Float min_abs_delta
            }

        String outlier_bed_param = "~{sample}.tissue_~{tissue}.METAFORA.outlier_regions.chrom_~{chrom}.bed"
        String outlier_z_mat_param = "~{sample}.tissue_~{tissue}.METAFORA.outlier_regions.zscore.chrom_~{chrom}.mat"

        command <<<
        mkdir -p outlier_plots
        Rscript ~{rscript} \
            --sample ~{sample} \
            --chrom ~{chrom} \
            --beta_mat ~{beta_matrix} \
            --depth_mat ~{depth_matrix} \
            --global_meth_pcs ~{global_pcs} \
            --outlier_bed ~{outlier_bed_param} \
            --outlier_z_mat ~{outlier_z_mat_param} \
            --plot_dir outlier_plots \
            --min_seg_size ~{min_segment_size} \
            --min_abs_zscore ~{min_abs_zscore} \
            --min_abs_delta ~{min_abs_delta} \
            --max_depth ~{max_depth} \
            --tissue ~{tissue}
        >>>
        
        runtime {
            docker: "quay.io/jonnguye/metafora:1.4"
            cpu: 4
            memory: "64GB"
            disks: "local-disk 128 SSD"
            }

        output {
            File outlier_bed = "~{outlier_bed_param}"
            File outlier_z_mat = "~{outlier_z_mat_param}"
            }
    }


task CallOutliersSexChromosomes {
        input {
                File rscript
                String sample
                String tissue
                String chrom
                File beta_matrix
                File beta_matrix_tbi
                File depth_matrix
                File depth_matrix_tbi
                File global_pcs
                File sex_chrom_df
                Int max_depth
                Int min_segment_size
                Int min_abs_zscore
                Float min_abs_delta
            }

        String outlier_bed_param = "~{sample}.tissue_~{tissue}.METAFORA.outlier_regions.chrom_~{chrom}.bed"
        String outlier_z_mat_param = "~{sample}.tissue_~{tissue}.METAFORA.outlier_regions.zscore.chrom_~{chrom}.mat"

        command <<<
        mkdir -p outlier_plots
        Rscript ~{rscript} \
            --sample ~{sample} \
            --chrom ~{chrom} \
            --beta_mat ~{beta_matrix} \
            --depth_mat ~{depth_matrix} \
            --global_meth_pcs ~{global_pcs} \
            --sex_chromosome_tsv ~{sex_chrom_df} \
            --chrY_seqname chrY \
            --outlier_bed ~{outlier_bed_param} \
            --outlier_z_mat ~{outlier_z_mat_param} \
            --plot_dir outlier_plots \
            --min_seg_size ~{min_segment_size} \
            --min_abs_zscore ~{min_abs_zscore} \
            --min_abs_delta ~{min_abs_delta} \
            --max_depth ~{max_depth} \
            --tissue ~{tissue}
        >>>
        
        runtime {
            docker: "quay.io/jonnguye/metafora:1.4"
            cpu: 4
            memory: "64GB"
            disks: "local-disk 128 SSD"
            }

        output {
            File outlier_bed = "~{outlier_bed_param}"
            File outlier_z_mat = "~{outlier_z_mat_param}"
            }
    }

task CombineChromOutliers {
        input {
                String sample
                String tissue
                Array[File] beds
                Array[File] mats
            }

        String outlier_bed_param = "~{sample}.tissue_~{tissue}.METAFORA.outlier_regions.bed"
        String outlier_z_mat_param = "~{sample}.tissue_~{tissue}.METAFORA.outlier_regions.zscore.mat"

        command <<<
            header=$(for i in $(ls ~{sep=' ' beds}); do head -1 $i; done | sort | tail -1)
            echo -e "$header" > ~{outlier_bed_param}
            cat ~{sep=' ' beds} | grep -v "$header" | grep -v "^$" >> ~{outlier_bed_param} || true

            header=$(for i in $(ls ~{sep=' ' mats}); do head -1 $i; done | sort | tail -1) 
            echo -e "$header" > ~{outlier_z_mat_param}
            cat ~{sep=' ' mats} | grep -v "$header" | grep -v "^$" >> ~{outlier_z_mat_param} || true
        >>>
        
        runtime {
            docker: "quay.io/jonnguye/metafora:1.4"
            cpu: 1
            memory: "12GB"
            disks: "local-disk 128 SSD"
            }

        output {
            File outlier_bed = "~{outlier_bed_param}"
            File outlier_z_mat = "~{outlier_z_mat_param}"
            }
    }

task CombineOutliers {
        input {
                String tissue
                Array[File] beds
            }

        String output_file = "METAFORA.tissue_~{tissue}.methylation_outliers.combined.tsv"

        command <<<
            header=$(for i in $(ls ~{sep=' ' beds}); do head -1 $i; done | sort | tail -1)
            echo "$header" > ~{output_file}
            cat {input} | grep -v "$header" | grep -v "^$" >> ~{output_file} || true
        >>>
        
        runtime {
            docker: "quay.io/jonnguye/metafora:1.4"
            cpu: 1
            memory: "24GB"
            disks: "local-disk 128 SSD"
            }

        output {
                File combined_outliers = "~{output_file}"
            }
    }

task SummaryPlots {
        input {
                File rscript
                File outlier_bed
                File covariates
                Float min_abs_zscore
                Float min_abs_delta
            }

        command <<<
        mkdir summary_figures
                Rscript ~{rscript} \
      --outlier_files ~{outlier_bed} \
      --covariates_files ~{covariates} \
      --plot_dir_out summary_figures \
      --summary_out summary_figures/METAFORA.outlier_count_per_sample_tissue.tsv \
      --min_abs_delta ~{min_abs_delta} \
      --min_abs_zscore ~{min_abs_zscore} 
        >>>
        
        runtime {
            docker: "quay.io/jonnguye/metafora:1.4"
            cpu: 1
            memory: "24GB"
            disks: "local-disk 128 SSD"
            }

        output {
                Array[File] summary_figures = glob("summary_figures/*")
            }
    }

task AnnotateHaplotypeDelta {
        input {
                File rscript
                String sample
                String tissue
                File outlier_bed
                File combined_bed
                File combined_bed_tbi
                File hap1_bed
                File hap1_bed_tbi
                File hap2_bed
                File hap2_bed_tbi
            }

        String annotated_bed_output = "~{sample}.tissue_~{tissue}.METAFORA.outlier_regions.haplotype_annotated.bed"

        command <<<
            Rscript ~{rscript} \
        --outlier_bed ~{outlier_bed} \
        --combined ~{combined_bed} \
        --hap1 ~{hap1_bed} \
        --hap2 ~{hap2_bed} \
        --annotated_out ~{annotated_bed_output}
        >>>
        
        runtime {
            docker: "quay.io/jonnguye/metafora:1.4"
            cpu: 1
            memory: "48GB"
            disks: "local-disk 128 SSD"
            }

        output {
                File annotated_bed = "~{annotated_bed_output}"
            }
    }
