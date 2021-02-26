version 1.0

workflow rna_seq_sars_cov_2_create_reference {
    input {
        Array[File] fasta
        Array[File] gtf
        String reference_name
        Int cpu = 16
        Int memory_gb = 43
        Int disk_space =  120
        String type
        String docker = "trinityctat/rna_seq_sars_cov_2:1.0.0"
        Boolean tar_gz_output = false
    }

    if(type=='smartseq2' || type=='bulk') {
        call rsem_prepare_reference {
            input:
                fasta=fasta,
                gtf=gtf,
                reference_name = reference_name,
                cpu=cpu,
                memory=memory_gb,
                disk_space=disk_space,
                docker=docker,
                tar_gz=tar_gz_output
        }
    }
    if(type=='10x') {
        call cellranger_mkref {
            input:
                fasta=fasta,
                gtf=gtf,
                reference_name = reference_name,
                cpu=cpu,
                memory=memory_gb,
                disk_space=disk_space,
                docker=docker,
                tar_gz=tar_gz_output
        }
    }
    output {
        Directory? reference = rsem_prepare_reference.reference
        File? tarred_reference = rsem_prepare_reference.tarred_reference
        File? fasta_index = rsem_prepare_reference.fasta_index
        File? fasta_output =  rsem_prepare_reference.fasta_output
        File? dict = rsem_prepare_reference.dict

        Directory? reference_10x = cellranger_mkref.reference
        File? tarred_reference_10x = cellranger_mkref.tarred_reference
    }
}

task rsem_prepare_reference {
    input {
        Array[File] fasta
        Array[File] gtf

        String reference_name
        Int cpu
        Int memory
        Int disk_space
        String docker
        Boolean tar_gz
    }

    command <<<
        set -e
        mkdir ~{reference_name}

        python <<CODE
        import os
        from subprocess import check_call

        def extract_files(file_list):
            extracted_list = []
            for i in range(len(file_list)):
                f = file_list[i]
                name, ext = os.path.splitext(f)
                if ext == '.gz':
                    check_call(['gunzip', '-f', f])
                    extracted_list.append(name)
                else:
                    extracted_list.append(f)
            return extracted_list


        fasta_list = extract_files('~{sep="," fasta}'.split(','))
        gtf_list = extract_files('~{sep="," gtf}'.split(','))

        gtf = gtf_list[0]
        fasta = fasta_list[0]
        if len(fasta_list) > 1:
            reference = '~{reference_name}'
            gtf = reference + '.gtf'
            fasta = reference + '.fa'
            for i in range(len(fasta_list)):
                check_call('cat {} >> {}'.format(fasta_list[i], reference + '.fa'), shell=True)
                check_call('cat {} >> {}'.format(gtf_list[i], reference + '.gtf'), shell=True)
        else:
            check_call(['mv', fasta, reference + '.fa'])
            fasta = reference + '.fa'
        args = ['rsem-prepare-reference', '--gtf', gtf, '--star', '-p', '~{cpu}', fasta, '~{reference_name}/rsem_ref']
        check_call(args)
        if '~{tar_gz}' == 'true':
            check_call(['tar', 'czf', '~{reference_name}.tar.gz', '~{reference_name}'])
            check_call(['mv', "~{reference_name}", 'tmp'])
        check_call(['samtools', 'faidx', fasta])
        check_call(['gatk', 'CreateSequenceDictionary', '-R', fasta])
        CODE
    >>>

    output {
        File? tarred_reference = "~{reference_name}.tar.gz"
        File fasta_index = "~{reference_name}.fa.fai"
        File fasta_output = "~{reference_name}.fa"
        File dict = "~{reference_name}.dict"
        Directory? reference = "~{reference_name}"
    }

    runtime {
        disks: "local-disk " + disk_space + " HDD"
        docker: docker
        preemptible: 2
        cpu: cpu
        memory: memory
    }
}

task cellranger_mkref {
    input {
        Array[File] fasta
        Array[File] gtf
        String reference_name
        Int disk_space
        Int cpu
        Int memory
        String docker
        Boolean tar_gz
    }

    command <<<
        set -e

        python <<CODE
        import os
        from subprocess import check_call

        def extract(file_list):
            extracted_list = []
            for i in range(len(file_list)):
                f = file_list[i]
                name, ext = os.path.splitext(f)
                if ext == '.gz':
                    check_call(['gunzip', '-f', f])
                    extracted_list.append(name)
                else:
                    extracted_list.append(f)
            return extracted_list

        fasta_list = extract('~{sep="," fasta}'.split(','))
        gtf_list = extract('~{sep="," gtf}'.split(','))
        genome_list = []
        for f in fasta_list:
            name, ext = os.path.splitext(f)
            genome_list.append(os.path.basename(name))
        call_args = ['cellranger', 'mkref']
        for genome, fasta, gtf in zip(genome_list, fasta_list, gtf_list):
            call_args.extend(['--genome=' + genome, '--fasta=' + fasta, '--genes=' + gtf])
        call_args.extend(['--nthreads=~{cpu}', '--memgb=~{memory}'])
        check_call(call_args)
        output_name = '_and_'.join(genome_list)
        check_call(['mv', output_name, '~{reference_name}'])
        if '~{tar_gz}' == 'true':
            check_call(['tar', 'czf', '~{reference_name}.tar.gz', '~{reference_name}'])
            check_call(['mv', "~{reference_name}", 'tmp'])
        CODE
    >>>

    output {
        Directory? reference = "~{reference_name}"
        File? tarred_reference = "~{reference_name}.tar.gz"
    }

    runtime {
        docker: docker
        memory: "~{memory}G"
        disks: "local-disk ~{disk_space} HDD"
        cpu: cpu
        preemptible: 2
    }
}
