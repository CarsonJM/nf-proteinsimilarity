#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CarsonJM/nf-protsim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/CarsonJM/nf-protsim
----------------------------------------------------------------------------------------
    Overview:
        1. Download latest ICTV VMR (Nextflow)
        2. Download ICTV genomes (process - VMR_to_fasta.py)
        3. Create DIAMOND database of ICTV genomes (process - DIAMOND)
        4. Split query viruses into chunks (process - seqkit)
        5. Align query virus genomes to ICTV database (process - DIAMOND)
        6. Perform self alignment of query genomes (process - DIAMOND)
        7. Calculate self score and normalized protein similarity (process - python)
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process VMRTOFASTA {
    label "process_single"

    input:
    path(vmr)

    output:
    path("${vmr.getBaseName()}.fna")            , emit: ictv_fasta
    path("processed_accessions_b.fa_names.tsv") , emit: processed_acc
    path("bad_accessions_b.tsv")                , emit: bad_acc
    path(".command.log")                        , emit: log
    path(".command.sh")                         , emit: log

    script:
    """
    # process VMR accessions
    VMR_to_fasta.py \\
        -mode VMR \\
        -ea B \\
        -VMR_file_name ${vmr} \\
        -v

    # download fasta file using current vmr
    VMR_to_fasta.py \\
        -email ${params.email} \\
        -mode fasta \\
        -ea b \\
        -fasta_dir ./ictv_fastas \\
        -VMR_file_name ${vmr} \\
        -v

    cat ictv_fastas/*/*.fa > ${vmr.getBaseName()}.fna
    """
}

process DIAMOND_MAKEDB {
    label "process_super_high"

    input:
    path(ictv_fasta)

    output:
    path("${ictv_fasta.getBaseName()}.dmnd")                , emit: dmnd
    path("${ictv_fasta.getBaseName()}.pyrodigalgv.faa.gz")  , emit: faa

    script:
    """
    # predict genes from FastA
    pyrodigal-gv \\
        -i ${ictv_fasta} \\
        -a ${ictv_fasta.getBaseName()}.pyrodigalgv.faa \\
        --jobs ${task.cpus}

    diamond \\
        makedb \\
        --threads ${task.cpus} \\
        --in ${ictv_fasta.getBaseName()}.pyrodigalgv.faa \\
        -d ${ictv_fasta.getBaseName()}

    gzip ${ictv_fasta.getBaseName()}.pyrodigalgv.faa
    """
}

process SEQKIT_SPLIT2 {
    label 'process_high'

    input:
    path(fasta)

    output:
    path("split_fastas/*")  , emit: split_fastas

    script:
    """
    seqkit \\
        split2 \\
            ${fasta} \\
            --threads ${task.cpus} \\
            --by-size ${params.chunk_size} \\
            --out-dir split_fastas
    """
}

process DIAMOND_BLASTP {
    label 'process_super_high'

    input:
    tuple val(meta), path(fasta)
    path(dmnd_db)

    output:
    tuple val(meta), path("${meta.id}.diamond_blastp.tsv.gz")   , emit: tsv
    tuple val(meta), path("${meta.id}.pyrodigalgv.faa.gz")      , emit: faa

    script:
    """
    # predict genes from FastA
    pyrodigal-gv \\
        -i ${fasta} \\
        -a ${meta.id}.pyrodigalgv.faa \\
        --jobs ${task.cpus}

    # align genes to DIAMOND reference db
    diamond \\
        blastp \\
        ${params.diamond_args} \\
        --query ${meta.id}.pyrodigalgv.faa \\
        --db ${dmnd_db} \\
        --threads ${task.cpus} \\
        --outfmt 6 \\
        --out ${meta.id}.diamond_blastp.tsv

    gzip ${meta.id}.diamond_blastp.tsv ${meta.id}.pyrodigalgv.faa
    """
}

process DIAMOND_SELF {
    label 'process_super_high'

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("${meta.id}.diamond_blastp.tsv.gz")   , emit: tsv

    script:
    """
    # make DIAMOND db for self alignment
    diamond \\
        makedb \\
        --threads ${task.cpus} \\
        --in ${faa} \\
        -d ${meta.id}

    # align genes to DIAMOND self db
    diamond \\
        blastp \\
        --masking none \\
        -k 1000 \\
        -e 1e-3 \\
        --faster \\
        --query ${faa} \\
        --db ${meta.id}.dmnd \\
        --threads ${task.cpus} \\
        --outfmt 6 \\
        --out ${meta.id}.diamond_blastp.tsv

    gzip ${meta.id}.diamond_blastp.tsv
    """
}

process SELFSCORE {
    label 'process_single'

    input:
    tuple val(meta), path(self_tsv)

    output:
    tuple val(meta), path("${meta.id}.selfscore.tsv")   , emit: tsv

    script:
    """
    self_score.py \\
        --input ${self_tsv} \\
        --output ${meta.id}.selfscore.tsv
    """
}

process NORMSCORE {
    label 'process_high'

    input:
    tuple val(meta), path(self_tsv), path(ref_tsv)

    output:
    tuple val(meta), path("${meta.id}.normscore.tsv")   , emit: tsv

    script:
    """
    norm_score.py \\
        --input ${ref_tsv} \\
        --self_score ${self_tsv} \\
        --min_score 0 \\
        --threads ${task.cpus} \\
        --output ${meta.id}.normscore.tsv
    """
}

process COMBINESCORES {
    label 'process_single'
    storeDir "."

    input:
    path(norm_scores)

    output:
    path("${params.output}")    , emit: tsv

    script:
    """
    # iterate over scores
    for table in ${norm_scores[0]}; do
       head -n 1 \${table} > ${params.output}
    done

    for table in ${norm_scores}; do
        tail -n +2 \${table} >> ${params.output}
    done
    """
}

// Run entry workflow
workflow {

    main:
    // Check if output file already exists
    def output_file = file("${params.output}")
    println workflow.configFiles
    if (!output_file.exists()) {

        ch_query_fasta = channel.fromPath(params.query_fasta).collect()

        // Prepare reference fasta file
        if (file("${params.ref_fasta}").exists()) {
            ch_ref_fna  = channel.fromPath(params.ref_fasta)
        } else if (params.vmr_url) {
            // 1. Download VMR file if necessary (Nextflow)
            ch_ictv_vmr = channel.fromPath(params.vmr_url)

            // 2. Download fasta (process - VMR_to_fasta.py)
            VMRTOFASTA(
                ch_ictv_vmr
            )
            ch_ref_fna = VMRTOFASTA.out.ictv_fasta
        }

        if (!file("${params.ref_db}").exists()) {
            // 3. Create DIAMOND database (process - pyrodiga-gv + DIAMOND)
            DIAMOND_MAKEDB(
                ch_ref_fna
            )
            ch_dmnd_db = DIAMOND_MAKEDB.out.dmnd
        } else {
            ch_dmnd_db  = channel.fromPath(params.ref_db)
        }

        // 4. Split query fasta file (process - seqkit)
        SEQKIT_SPLIT2(
            ch_query_fasta
        )

        ch_split_fastas = SEQKIT_SPLIT2.out.split_fastas
            .map { file -> file }
            .flatten()
            .map { file ->
                [ [ id: file.getBaseName() ], file ]
            }

        // 5. Run DIAMOND against ref db (process - pyrodigal-gv + DIAMOND)
        DIAMOND_BLASTP(
            ch_split_fastas,
            ch_dmnd_db.collect()
        )

        // 6. Run DIAMOND self alignment (process - DIAMOND)
        DIAMOND_SELF(
            DIAMOND_BLASTP.out.faa
        )

        // 7. Calculate self score (process - self_score.py)
        SELFSCORE(
            DIAMOND_SELF.out.tsv
        )

        // 8. Calculate normalized bitscore (process - norm_score.py)
        NORMSCORE(
            SELFSCORE.out.tsv.combine(DIAMOND_BLASTP.out.tsv, by:0)
        )

        // 9. Combine results (process - CAT)
        COMBINESCORES(
            NORMSCORE.out.tsv.map { _meta, tsvs -> [ tsvs ] }.collect()
        )

    } else {
        println "Output file [${params.output}] already exists! Skipping nf-proteinsimilarity."
    }

    // Delete intermediate and Nextflow-specific files
    def remove_tmp = params.remove_tmp
    workflow.onComplete {
        if (output_file.exists()) {
            def work_dir = new File("./work/")
            def nextflow_dir = new File("./.nextflow/")
            def launch_dir = new File(".")

            work_dir.deleteDir()
            nextflow_dir.deleteDir()
            launch_dir.eachFileRecurse { file ->
                if (file.name ==~ /\.nextflow\.log.*/) {
                    file.delete()
                }
            }

            if (remove_tmp) {
                def tmp_dir = new File("./tmp/")
                tmp_dir.deleteDir()
            }
        }
    }
}