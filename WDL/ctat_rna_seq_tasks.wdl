version 1.0

workflow ctat_rna_seq_tasks {}


task minimap2 {
    input {
        File? ref_fasta
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
        docker: "trinityctat/ctat_rna_seq:1.0.0"
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
        File? ref_fasta

    }

    output {
        File igvjs_viewer = "igvjs_viewer.html"
    }

    runtime {
        preemptible: 2
        disks: "local-disk " + ceil(size(ref_fasta, "GB")+size(viral_bam, "GB")+size(assembled_bam, "GB")+10) + " HDD"
        docker: "trinityctat/ctat_rna_seq:1.0.0"
        cpu: 1
        memory: "2G"
    }
    command <<<
        samtools faidx ~{ref_fasta}
        samtools index ~{vcf}
        create_report ~{vcf} ~{ref_fasta} --tracks ~{assembled_bam} ~{viral_bam}
    >>>
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
            compress="pigz"
            if [[ $genomeDir == *.bz2 ]] ; then
                compress="pbzip2"
            fi
            tar -I $compress -xf $genomeDir -C genome_dir --strip-components 1
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
        docker: "trinityctat/ctat_rna_seq:1.0.0"
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
            compress="pigz"
            if [[ $genomeDir == *.bz2 ]] ; then
                compress="pbzip2"
            fi
            tar -I $compress -xf $genomeDir -C genome_dir --strip-components 1
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
        docker: "trinityctat/ctat_rna_seq:1.0.0"
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
        docker: "trinityctat/ctat_rna_seq:1.0.0"
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
        docker: "trinityctat/ctat_rna_seq:1.0.0"
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
        docker: "trinityctat/ctat_rna_seq:1.0.0"
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
        docker: "trinityctat/ctat_rna_seq:1.0.0"
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: 2
    }
}
