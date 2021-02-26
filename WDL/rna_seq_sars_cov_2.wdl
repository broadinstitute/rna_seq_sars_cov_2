version 1.0

#import "https://api.firecloud.org/ga4gh/v1/tools/CTAT:ctat_mutations/versions/1/plain-WDL/descriptor" as ctat_mutations
#import "https://api.firecloud.org/ga4gh/v1/tools/CTAT:trinityrnaseq/versions/1/plain-WDL/descriptor" as trinity_assembly
#import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cumulus/versions/34/plain-WDL/descriptor" as cumulus

import "./ctat_mutations.wdl" as ctat_mutations
import "./trinityrnaseq.wdl" as trinityrnaseq
import "./cumulus.wdl" as cumulus

workflow rna_seq_sars_cov_2 {
    input {
        String sample_id
        Array[File] read1
        Array[File] read2

        File viral_ref_fasta

        File reference # tar.gz file of RSEM or cellranger reference
        Int star_cpu = 12
        String star_memory = "43G"
        Int rsem_cpu = 4
        String rsem_memory = "32G"
        Int? genome_guided_max_intron = 29000

        String minimap2_memory = "2G"
        String minimap2_flags = "asm20"

        # resources
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        File? gtf

        File? db_snp_vcf
        File? db_snp_vcf_index

        File? gnomad_vcf
        File? gnomad_vcf_index

        File? rna_editing_vcf
        File? rna_editing_vcf_index

        File? repeat_mask_bed
        File? ref_splice_adj_regions_bed

        Boolean call_host_variants = true
        Boolean apply_bqsr = true


        String type # one of bulk, 10x
        # single cell
        Int downsample_max_coverage = 100
        String downsample_bam_memory = "30G"
        Int cellranger_cpu = 32
        String cellranger_memory = "120G"
        Int cellranger_expect_cells = 3000
        String cellranger_chemistry = "auto"
    }


    if(type=='10x') {
        call cellranger_count {
            input:
                read1=read1,
                read2=read2,
                cpu=cellranger_cpu,
                memory=cellranger_memory,
                id=sample_id,
                reference=reference,
                expect_cells=cellranger_expect_cells,
                chemistry=cellranger_chemistry
        }

        call cumulus.cumulus as cumulus {
            input:
                input_file=cellranger_count.filtered_feature_bc_matrix,
                output_name=sample_id,
                output_directory=""
        }

        call downsample_bam {
            input:
                input_bam=select_first([cellranger_count.bam]),
                input_bam_index=select_first([cellranger_count.bam_index]),
                output_name=sample_id + '.downsampled',
                memory=downsample_bam_memory,
                max_coverage=downsample_max_coverage
        }
    }

    if(type=='bulk') {
        call star {
            input:
                read1=read1,
                read2=read2,
                cpu=star_cpu,
                memory=star_memory,
                base_name=sample_id,
                reference=reference,
                output_unmapped_reads=false
        }

        call rsem {
            input:
                sample_name=sample_id,
                reference=reference,
                is_paired=defined(read2),
                cpu=rsem_cpu,
                memory=rsem_memory,
                bam=star.transcript_coords_bam
        }
    }

    if(type=='smartseq2') {
        scatter(i in range(length(read1))) {
            String cell_name = basename(basename(basename(read1[i], ".gz"), ".fq"), ".fastq")
            Boolean is_paired = i < length(read2)
            call star as star_ss2{
                input:
                    read1=[read1[i]],
                    read2=if is_paired then [read2[i]] else [],
                    cpu=star_cpu,
                    memory=star_memory,
                    base_name=cell_name,
                    reference=reference,
                    output_unmapped_reads=false
            }

            call rsem as rsem_ss2 {
                input:
                    sample_name=cell_name,
                    reference=reference,
                    is_paired=is_paired,
                    cpu=rsem_cpu,
                    memory=rsem_memory,
                    bam=star_ss2.transcript_coords_bam
            }

            call ctat_mutations.ctat_mutations as ctat_mutations_ss2 {
                input:
                    bam=star_ss2.sorted_bam,
                    bai=star_ss2.sorted_bam_index,
                    merge_extra_fasta=false,
                    sample_id=cell_name,
                    extra_fasta=viral_ref_fasta,
                    filter_cancer_variants=false,
                    call_variants=call_host_variants,
                    ref_dict=ref_dict,
                    ref_fasta=ref_fasta,
                    ref_fasta_index=ref_fasta_index,
                    gtf=gtf,
                    db_snp_vcf=db_snp_vcf,
                    db_snp_vcf_index=db_snp_vcf_index,
                    gnomad_vcf=gnomad_vcf,
                    gnomad_vcf_index=gnomad_vcf_index,
                    rna_editing_vcf=rna_editing_vcf,
                    rna_editing_vcf_index=rna_editing_vcf_index,
                    repeat_mask_bed=repeat_mask_bed,
                    ref_splice_adj_regions_bed=ref_splice_adj_regions_bed,
                    haplotype_caller_args_for_extra_reads="-dont-use-soft-clipped-bases --stand-call-conf 20 --recover-dangling-heads true --sample-ploidy 1"
            }

            Int? extra_bam_number_of_reads = ctat_mutations_ss2.extra_bam_number_of_reads

            if(extra_bam_number_of_reads>0) {
                call trinityrnaseq.trinityrnaseq as trinityrnaseq_ss2 {
                    input:
                        output_base_name=cell_name,
                        genome_guided_bam=ctat_mutations_ss2.extra_bam,
                        genome_guided_max_intron=genome_guided_max_intron,
                }

                call minimap2 as minimap2_ss2 {
                    input:
                        ref_fasta=viral_ref_fasta,
                        memory=minimap2_memory,
                        flags=minimap2_flags,
                        assembled_fasta=select_first([trinityrnaseq_ss2.guided_fasta]),
                        output_base_name=cell_name
                }

                call create_report as create_report_ss2 {
                    input:
                        vcf=ctat_mutations.extra_vcf,
                        ref_fasta=viral_ref_fasta,
                        assembled_bam=minimap2_ss2.bam,
                        assembled_bam_index=minimap2_ss2.bam_index,
                        viral_bam=ctat_mutations_ss2.extra_bam,
                        viral_bam_index=ctat_mutations_ss2.extra_bam_index
                }
            }
        }
        call generate_rsem_matrix {
            input:
                gene_results=rsem_ss2.rsem_gene,
                count_results=rsem_ss2.rsem_cnt,
                output_name=sample_id,
                memory="2G"
        }
        call cumulus.cumulus as cumulus_ss2 {
            input:
                input_file=generate_rsem_matrix.count_matrix,
                output_name=sample_id,
                output_directory=""
        }
    }

    if(type=='10x' || type=='bulk') {
        call ctat_mutations.ctat_mutations as ctat_mutations {
            input:
                bam=select_first([star.sorted_bam, downsample_bam.bam]),
                bai=select_first([star.sorted_bam_index, downsample_bam.bam_index]),
                merge_extra_fasta=false,
                sample_id=sample_id,
                extra_fasta=viral_ref_fasta,
                filter_cancer_variants=false,
                call_variants=call_host_variants,
                ref_dict=ref_dict,
                ref_fasta=ref_fasta,
                ref_fasta_index=ref_fasta_index,
                gtf=gtf,
                db_snp_vcf=db_snp_vcf,
                db_snp_vcf_index=db_snp_vcf_index,
                gnomad_vcf=gnomad_vcf,
                gnomad_vcf_index=gnomad_vcf_index,
                rna_editing_vcf=rna_editing_vcf,
                rna_editing_vcf_index=rna_editing_vcf_index,
                repeat_mask_bed=repeat_mask_bed,
                ref_splice_adj_regions_bed=ref_splice_adj_regions_bed,
                haplotype_caller_args_for_extra_reads="-dont-use-soft-clipped-bases --stand-call-conf 20 --recover-dangling-heads true --sample-ploidy 1"
        }

        Int? extra_bam_number_of_reads = ctat_mutations.extra_bam_number_of_reads

        if(extra_bam_number_of_reads>0) {
            call trinityrnaseq.trinityrnaseq as trinityrnaseq {
                input:
                    output_base_name=sample_id,
                    genome_guided_bam=ctat_mutations.extra_bam,
                    genome_guided_max_intron=genome_guided_max_intron,
            }

            call minimap2 {
                input:
                    ref_fasta=viral_ref_fasta,
                    memory=minimap2_memory,
                    flags=minimap2_flags,
                    assembled_fasta=select_first([trinityrnaseq.guided_fasta]),
                    output_base_name=sample_id
            }

            call create_report {
                input:
                    vcf=ctat_mutations.extra_vcf,
                    ref_fasta=viral_ref_fasta,
                    assembled_bam=minimap2.bam,
                    assembled_bam_index=minimap2.bam_index,
                    viral_bam=ctat_mutations.extra_bam,
                    viral_bam_index=ctat_mutations.extra_bam_index
            }
        }
    }

    output {
        File? star_output_log_final =star.output_log_final
        File? star_output_SJ = star.output_SJ
        File? transcript_coords_bam = star.transcript_coords_bam
        File? aligned_sorted_bam = star.sorted_bam
        File? aligned_sorted_bam_index = star.sorted_bam_index

        File? rsem_gene = rsem.rsem_gene
        File? rsem_isoform = rsem.rsem_isoform
        File? rsem_cnt = rsem.rsem_cnt
        File? rsem_model = rsem.rsem_model
        File? rsem_theta = rsem.rsem_theta


        Array[File]? star_output_log_final_ss2 =star_ss2.output_log_final
        Array[File]? star_output_SJ_ss2 = star_ss2.output_SJ
        Array[File]? transcript_coords_bam_ss2 = star_ss2.transcript_coords_bam
        Array[File]? aligned_sorted_bam_ss2 = star_ss2.sorted_bam
        Array[File]? aligned_sorted_bam_index_ss2 = star_ss2.sorted_bam_index

        Array[File]? rsem_gene_ss2 = rsem_ss2.rsem_gene
        Array[File]? rsem_isoform_ss2 = rsem_ss2.rsem_isoform
        Array[File]? rsem_cnt_ss2 = rsem_ss2.rsem_cnt
        Array[File]? rsem_model_ss2 = rsem_ss2.rsem_model
        Array[File]? rsem_theta_ss2 = rsem_ss2.rsem_theta

        File? viral_bam=ctat_mutations.extra_bam
        File? viral_bam_index=ctat_mutations.extra_bam_index
        File? minimap2_bam = minimap2.bam
        File? minimap2_bam_index = minimap2.bam_index
        File? haplotype_caller_vcf = ctat_mutations.haplotype_caller_vcf
        File? trinity_assembly_guided_fasta = trinityrnaseq.fasta
        File? trinity_assembly_guided_gene_trans_map = trinityrnaseq.gene_trans_map
        File? igvjs_viewer = create_report.igvjs_viewer

        Array[File?]? viral_bam_ss2=ctat_mutations_ss2.extra_bam
        Array[File?]? viral_bam_index_ss2=ctat_mutations_ss2.extra_bam_index
        Array[File?]? minimap2_bam_ss2 = minimap2_ss2.bam
        Array[File?]? minimap2_bam_index_ss2 = minimap2_ss2.bam_index
        Array[File?]? haplotype_caller_vcf_ss2 = ctat_mutations_ss2.haplotype_caller_vcf
        Array[File?]? trinity_assembly_guided_fasta_ss2 = trinityrnaseq_ss2.fasta
        Array[File?]? trinity_assembly_guided_gene_trans_map_ss2 = trinityrnaseq_ss2.gene_trans_map
        Array[File?]? igvjs_viewer_ss2 = create_report_ss2.igvjs_viewer
        File? count_matrix_ss2 = generate_rsem_matrix.count_matrix
        File? qc_report_ss2 = generate_rsem_matrix.qc_report

        File? cellranger_web_summary = cellranger_count.web_summary
        File? cellranger_metrics_summary = cellranger_count.metrics_summary
        File? cellranger_filtered_feature_bc_matrix = cellranger_count.filtered_feature_bc_matrix
        File? cellranger_raw_feature_bc_matrix = cellranger_count.raw_feature_bc_matrix
        File? cellranger_molecule_info = cellranger_count.molecule_info
        File? cellranger_bam = cellranger_count.bam
        File? cellranger_bam_index = cellranger_count.bam_index
        File? downsampled_bam = downsample_bam.bam
        File? downsampled_bam_index = downsample_bam.bam_index



        Array[File]? filt_xlsx = cumulus.output_filt_xlsx
        Array[File]? filt_plot = cumulus.output_filt_plot
        Array[File]? hvf_plot = cumulus.output_hvf_plot
        Array[File]? dbl_plot = cumulus.output_dbl_plot
        Array[File?]? de_h5ad = cumulus.output_de_h5ad
        Array[File?]? de_xlsx = cumulus.output_de_xlsx
        Array[Array[File]?]? pdfs = cumulus.output_pdfs


        Array[File]? filt_xlsx_ss2 = cumulus_ss2.output_filt_xlsx
        Array[File]? filt_plot_ss2 = cumulus_ss2.output_filt_plot
        Array[File]? hvf_plot_ss2 = cumulus_ss2.output_hvf_plot
        Array[File]? dbl_plot_ss2 = cumulus_ss2.output_dbl_plot
        Array[File?]? de_h5ad_ss2 = cumulus_ss2.output_de_h5ad
        Array[File?]? de_xlsx_ss2 = cumulus_ss2.output_de_xlsx
        Array[Array[File]?]? pdfs_ss2 = cumulus_ss2.output_pdfs


    }
}

task minimap2 {
    input {
        File ref_fasta
        File assembled_fasta
        String output_base_name
        String flags
        String memory
    }

    command <<<
        set -e

        minimap2 -a -cx ~{flags} --cs ~{ref_fasta} ~{assembled_fasta} > ~{output_base_name}.sam
        samtools sort ~{output_base_name}.sam > ~{output_base_name}.bam
        samtools index ~{output_base_name}.bam
    >>>

    output {
        File bam = "~{output_base_name}.bam"
        File bam_index = "~{output_base_name}.bam.bai"
    }

    runtime {
        docker: "trinityctat/rna_seq_sars_cov_2:1.0.0"
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(ref_fasta, "GB") + size(assembled_fasta, "GB") + 2) + " HDD"
        cpu: 1
        preemptible: 2
    }
}

task create_report {
    input {
        File? vcf
        File? viral_bam
        File? viral_bam_index
        File? assembled_bam
        File? assembled_bam_index
        File ref_fasta

    }

    output {
        File igvjs_viewer = "igvjs_viewer.html"
    }

    runtime {
        preemptible: 2
        disks: "local-disk " + ceil(size(ref_fasta, "GB")+size(viral_bam, "GB")+size(assembled_bam, "GB")+10) + " HDD"
        docker: "trinityctat/rna_seq_sars_cov_2:1.0.0"
        cpu: 1
        memory: "2G"
    }
    command <<<
        samtools faidx ~{ref_fasta}
        samtools index ~{vcf}
        create_report ~{vcf} ~{ref_fasta} --tracks ~{assembled_bam} ~{viral_bam}
    >>>
}

task star_solo {
    input {
        File reference
        Array[File] read1
        Array[File]? read2
        String base_name
        Int cpu
        String memory
        File? whitelist
        Boolean output_unmapped_reads
    }
    Float extra_disk_space = 30
    Boolean is_paired = defined(read2)
    Float fastq_disk_space_multiplier = if(is_paired) then 20 else 10
    Boolean is_gzip = sub(read1[0], "^.+\\.(gz)$", "GZ") == "GZ"
    Float fastq_size = size(read1, "GB")*fastq_disk_space_multiplier
    Boolean no_whitelist = !defined(whitelist)

    command <<<
        set -e

        genomeDir="~{reference}"


        if [ -f "${genomeDir}" ] ; then
            mkdir genome_dir
            tar xf ~{reference} -C genome_dir --strip-components 1
            genomeDir="genome_dir"
        fi

        STAR \
        --genomeDir $genomeDir \
        --runThreadN ~{cpu} \
        --readFilesIn ~{sep=',' read1} ~{sep=',' read2} \
        ~{true='--readFilesCommand "gunzip -c"' false='' is_gzip} \
        --outSAMtype BAM SortedByCoordinate \
        --soloType Droplet \
        ~{"--soloCBwhitelist " + whitelist} \
        ~{true='--soloCBwhitelist None' false='' no_whitelist} \
        --outFileNamePrefix ~{base_name}. \
        ~{true='--outReadsUnmapped Fastx' false='' output_unmapped_reads}

        if [ "~{output_unmapped_reads}" == "true" ] ; then
            mv ~{base_name}.Unmapped.out.mate1 ~{base_name}.Unmapped.out.mate1.fastq
            mv ~{base_name}.Unmapped.out.mate2 ~{base_name}.Unmapped.out.mate2.fastq
        fi
        samtools index "~{base_name}.Aligned.sortedByCoord.out.bam"

    >>>

    output {
        #        File bam = "~{base_name}.Aligned.out.bam"
        File sorted_bam = "~{base_name}.Aligned.sortedByCoord.out.bam"
        File sorted_bam_index = "~{base_name}.Aligned.sortedByCoord.out.bam.bai"

        File output_log_final = "~{base_name}.Log.final.out"
        File output_SJ = "~{base_name}.SJ.out.tab"
        File barcodes = "~{base_name}.Solo.out/barcodes.tsv"
        File gene_stats = "~{base_name}.Solo.out/Gene.stats"
        File genes = "~{base_name}.Solo.out/genes.tsv"
        File matrix = "~{base_name}.Solo.out/matrix.mtx"
        Array[File] unmapped_reads = glob("~{base_name}.Unmapped.out.*")
    }

    runtime {
        preemptible: 2
        disks: "local-disk " + ceil(fastq_size + size(reference, "GB")*8 + extra_disk_space) + " HDD"
        docker: "trinityctat/rna_seq_sars_cov_2:1.0.0"
        cpu: cpu
        memory: memory
    }

}

task star {
    input {
        File reference
        Array[File] read1
        Array[File]? read2
        String base_name
        Boolean output_unmapped_reads
        Int cpu
        String memory
    }
    Float extra_disk_space = 30
    Boolean is_paired = defined(read2)
    Float fastq_disk_space_multiplier = if(is_paired) then 20 else 10
    Boolean is_gzip = sub(read1[0], "^.+\\.(gz)$", "GZ") == "GZ"
    Float fastq_size = size(read1, "GB")*fastq_disk_space_multiplier
    command <<<
        set -e

        genomeDir="~{reference}"

        if [ -f "${genomeDir}" ] ; then
            mkdir genome_dir
            tar xf ~{reference} -C genome_dir --strip-components 1
            genomeDir="genome_dir"
        fi

        STAR \
        --genomeDir $genomeDir \
        --runThreadN ~{cpu} \
        --readFilesIn ~{sep=',' read1} ~{sep=',' read2} \
        ~{true='--readFilesCommand "gunzip -c"' false='' is_gzip} \
        --outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM \
        --outSAMheaderHD \@HD VN:1.4 SO:unsorted \
        --twopassMode Basic \
        --limitBAMsortRAM 30000000000 \
        --outSAMmapqUnique 60 \
        --outFileNamePrefix ~{base_name}. \
        ~{true='--outReadsUnmapped Fastx' false='' output_unmapped_reads}


        if [ "~{output_unmapped_reads}" == "true" ] ; then
            mv ~{base_name}.Unmapped.out.mate1 ~{base_name}.Unmapped.out.mate1.fastq
            mv ~{base_name}.Unmapped.out.mate2 ~{base_name}.Unmapped.out.mate2.fastq
        fi

        samtools sort ~{base_name}.Aligned.out.bam -o ~{base_name}.Aligned.sortedByCoord.out.bam
        samtools index "~{base_name}.Aligned.sortedByCoord.out.bam"

    >>>

    output {
        File transcript_coords_bam = "~{base_name}.Aligned.toTranscriptome.out.bam"
        #        File bam = "~{base_name}.Aligned.out.bam"
        File sorted_bam = "~{base_name}.Aligned.sortedByCoord.out.bam"
        File sorted_bam_index = "~{base_name}.Aligned.sortedByCoord.out.bam.bai"

        File output_log_final = "~{base_name}.Log.final.out"
        File output_SJ = "~{base_name}.SJ.out.tab"
        Array[File] unmapped_reads = glob("~{base_name}.Unmapped.out.*")
    }

    runtime {
        preemptible: 2
        disks: "local-disk " + ceil(fastq_size + size(reference, "GB")*8 + extra_disk_space) + " HDD"
        docker: "trinityctat/rna_seq_sars_cov_2:1.0.0"
        cpu: cpu
        memory: memory
    }

}

task rsem {
    input {
        File reference
        File bam
        String sample_name
        Int cpu
        String memory
        Boolean is_paired

    }
    Float extra_disk_space = 10
    Float disk_space_multiplier = 4

    command <<<
        set -e

        genomeDir="~{reference}"

        if [ -f "${genomeDir}" ] ; then
            mkdir genome_dir
            tar xf ~{reference} -C genome_dir --strip-components 1
            genomeDir="genome_dir"
        fi

        REFERENCE_NAME="$(basename `ls $genomeDir/*.grp` .grp)"

        rsem-calculate-expression ~{true="--paired-end" false="" is_paired} \
        --bam \
        -p ~{cpu} \
        ~{bam} \
        $genomeDir/$REFERENCE_NAME ~{sample_name}
    >>>

    output {
        File rsem_gene = "~{sample_name}.genes.results"
        File rsem_isoform = "~{sample_name}.isoforms.results"
        File rsem_cnt = "~{sample_name}.stat/~{sample_name}.cnt"
        File rsem_model = "~{sample_name}.stat/~{sample_name}.model"
        File rsem_theta = "~{sample_name}.stat/~{sample_name}.theta"
    }

    runtime {
        docker: "trinityctat/rna_seq_sars_cov_2:1.0.0"
        memory:memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(reference, "GB")*5 + disk_space_multiplier * size(bam, "GB")  + extra_disk_space)+ " HDD"
        cpu: cpu
        preemptible: 2
    }
}

task downsample_bam {
    input {
        File input_bam
        File input_bam_index
        String output_name
        Int max_coverage
        String memory
    }

    command <<<
        bamsifter -c ~{max_coverage} -o ~{output_name}.bam ~{input_bam}
        samtools index ~{output_name}.bam
    >>>

    runtime {
        docker: "trinityctat/rna_seq_sars_cov_2:1.0.0"
        bootDiskSizeGb: 12
        memory: memory
        disks: "local-disk " + ceil(1 + size(input_bam, "GB")*2) + " HDD"
        preemptible: 2
        cpu: 1
    }

    output {
        File bam = "~{output_name}.bam"
        File bam_index = "~{output_name}.bam.bai"
    }
}

task cellranger_count {
    input {
        File reference
        Array[File] read1
        Array[File]? read2
        String id
        Int cpu
        String memory
        Int expect_cells
        String chemistry
    }

    command <<<
        set -e

        python <<CODE
        import os
        import subprocess
        reference="~{reference}"
        if not os.path.isdir(reference):
            os.mkdir("reference")
            subprocess.check_call(["tar", "xf", reference, "-C", "reference", "--strip-components", "1"])
            reference = "reference"

        id="~{id}"
        chemistry="~{chemistry}"
        expect_cells="~{expect_cells}"
        read1_list="~{sep=',' read1}".split(",")
        read2_list="~{sep=',' read2}".split(",")

        os.mkdir("fastqs")
        for i in range(len(read1_list)):
            read1=read1_list[i]
            read2=read2_list[i]
            renamed_read1 = "{}_S1_L00{}_R1_001.fastq".format(id, i+1)
            renamed_read2 = "{}_S1_L00{}_R2_001.fastq".format(id, i+1)
            if read1.endswith('.gz'):
                renamed_read1 += '.gz'
                renamed_read2 += '.gz'
            os.rename(read1, os.path.join("fastqs", renamed_read1))
            os.rename(read2, os.path.join("fastqs", renamed_read2))

        subprocess.check_call(["cellranger", "count", "--id=" + id, "--transcriptome=" + reference, "--fastqs=fastqs", "--expect-cells=" + expect_cells, "--chemistry=" + chemistry, "--nosecondary"])
        CODE
    >>>

    runtime {
        docker: "trinityctat/rna_seq_sars_cov_2:1.0.0"
        bootDiskSizeGb: 12
        memory: memory
        disks: "local-disk " + ceil(100 + size(read1, "GB")*10) + " HDD"
        preemptible: 2
        cpu: cpu
    }

    #- Run summary HTML:                         /opt/sample345/outs/web_summary.html
    #- Run summary CSV:                          /opt/sample345/outs/metrics_summary.csv
    #- BAM:                                      /opt/sample345/outs/possorted_genome_bam.bam
    #- BAM index:                                /opt/sample345/outs/possorted_genome_bam.bam.bai
    #- Filtered feature-barcode matrices MEX:    /opt/sample345/outs/filtered_feature_bc_matrix
    #- Filtered feature-barcode matrices HDF5:   /opt/sample345/outs/filtered_feature_bc_matrix.h5
    #- Unfiltered feature-barcode matrices MEX:  /opt/sample345/outs/raw_feature_bc_matrix
    #- Unfiltered feature-barcode matrices HDF5: /opt/sample345/outs/raw_feature_bc_matrix.h5
    #- Secondary analysis output CSV:            /opt/sample345/outs/analysis
    #- Per-molecule read information:            /opt/sample345/outs/molecule_info.h5
    output {
        File web_summary = "~{id}/outs/web_summary.html"
        File metrics_summary = "~{id}/outs/metrics_summary.csv"
        File bam = "~{id}/outs/possorted_genome_bam.bam"
        File bam_index = "~{id}/outs/possorted_genome_bam.bam.bai"
        File filtered_feature_bc_matrix = "~{id}/outs/filtered_feature_bc_matrix.h5"
        File raw_feature_bc_matrix = "~{id}/outs/raw_feature_bc_matrix.h5"
        File molecule_info = "~{id}/outs/molecule_info.h5"
    }
}

task filter_bam {
    input {
        String region
        File input_bam
        File input_bam_index
        String output_name
    }

    command <<<
        samtools view -h -b ~{input_bam} ~{region} > ~{output_name}.bam
        samtools view -c -F 260 ~{output_name}.bam > "nreads.txt"
        samtools index "~{output_name}.bam"
    >>>

    runtime {
        docker: "trinityctat/rna_seq_sars_cov_2:1.0.0"
        bootDiskSizeGb: 12
        memory: "2G"
        disks: "local-disk " + ceil(1 + size(input_bam, "GB")*2) + " HDD"
        preemptible: 2
        cpu: 1
    }

    output {
        File bam = "~{output_name}.bam"
        File bam_index = "~{output_name}.bam.bai"
        Int nreads = read_int('nreads.txt')
    }
}


task generate_rsem_matrix {
    input {
        Array[File] gene_results
        Array[File] count_results
        String output_name
        String memory
    }
    Float disk_space = size(gene_results, "GB")+2+size(count_results, "GB")+2 + 10

    command <<<
        set -e

        python <<CODE
        import os
        import numpy as np
        import pandas as pd
        gene_names = None
        barcodes = []
        cntmat = []
        normalize_tpm_by_sequencing_depth = True
        for result_file in "~{sep=',' gene_results}".split(','):
            barcodes.append(os.path.basename(result_file)[:-len('.genes.results')])
            df = pd.read_table(result_file, header = 0, index_col = 0)
            if gene_names is None:
                gene_names = df.index
            if normalize_tpm_by_sequencing_depth:
                tot_counts = df['expected_count'].sum()
                counts = df['TPM'].values / 10.0 # convert TPMs into TP100Ks
                denom = counts.sum()
                if denom > 0:
                    counts = (counts / denom * tot_counts + 0.5).astype(int)
            else:
                counts = df['TPM'].values
            cntmat.append(counts)
        df_idx = pd.Index(gene_names, name = 'GENE')
        df_out = pd.DataFrame(data = np.stack(cntmat, axis = 1), index = df_idx, columns = barcodes)
        df_out.to_csv('~{output_name}.dge.txt.gz', sep = '\t', compression = 'gzip')

        arr = []
        barcodes = []
        for result_file in "~{sep=',' count_results}".split(','):
            barcodes.append(os.path.basename(result_file)[:-len('.cnt')])
            with open(result_file) as fin:
                Ns = [int(x) for x in next(fin).strip().split(' ')]
                align_values = [int(x) for x in next(fin).strip().split(' ')]
                res = [str(Ns[3]), str(round(Ns[1] * 100.0 / Ns[3], 2)) + "%", str(round(align_values[0] * 100.0 / Ns[3], 2)) + "%"] if Ns[3] > 0 else ['0%', '0']
                arr.append(res)
        df = pd.DataFrame(data = np.array(arr), index = barcodes, columns = ["Total reads", "Alignment rate", "Unique rate"])
        df.index.name = "Cell"
        df.to_csv('~{output_name}.qc.stat.tsv', sep = '\t')
        CODE
    >>>

    output {
        File count_matrix = "~{output_name}.dge.txt.gz"
        File qc_report = "~{output_name}.qc.stat.tsv"
    }

    runtime {
        docker: "trinityctat/rna_seq_sars_cov_2:1.0.0"
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: 2
    }
}
