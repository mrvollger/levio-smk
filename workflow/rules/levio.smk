rule fastq_input:
    input:
        reads=lambda wc: config["inputs"][wc.sm],
    output:
        fastq=temp("temp/{sm}/fastq_input/{sm}.fastq.gz"),
    threads: MAX_THREADS // 2
    resources:
        mem_mb=32 * 1024,
    conda:
        DEFAULT_ENV
    shell:
        """
        # check if input is bam
        if [[ {input.reads} =~ .*\\.(bam|sam|cram) ]]; then
            samtools collate -@ {threads} -u -O {input.reads} \
                | samtools fastq -@ {threads} -0 /dev/null \
                | bgzip -@ {threads} \
                > {output.fastq}
        else
            ln -s {input.reads} {output.fastq}
        fi
        """


rule bwa_mem2:
    input:
        ref=DSA,
        fastq=rules.fastq_input.output.fastq,
    output:
        bam="results/DSA/{sm}-bwamem2.bam",
        csi="results/DSA/{sm}-bwamem2.bam.csi",
    threads: MAX_THREADS
    resources:
        mem_mb=MAX_THREADS * 4 * 1024,
    conda:
        DEFAULT_ENV
    params:
        reset_mapq=workflow.source_path("../scripts/reset-mapq.py"),
    shell:
        """
        bwa-mem2 mem -t {threads} -p {input.ref} {input.fastq} \
            | python {params.reset_mapq} \
            | samtools sort -@ {threads} -m 1G -o {output.bam} --write-index
        """


rule leviosam2:
    input:
        bam=rules.bwa_mem2.output.bam,
        levio_index=LEVIO_INDEX,
        ref=REF,
    output:
        bam=temp("temp/{sm}/leviosam2/{sm}-committed.bam"),
        deferred=temp("temp/{sm}/leviosam2/{sm}-deferred.bam"),
        unliftable=temp("temp/{sm}/leviosam2/{sm}-unliftable.bam"),
    threads: MAX_THREADS
    resources:
        mem_mb=MAX_THREADS * 4 * 1024,
    conda:
        DEFAULT_ENV
    params:
        # maximum number of CIGAR opts to change, also the max gap size that can be spanned
        G=config.get("levio_G", 100_000),
        S=config.get(
            "levio_S",
            "-S mapq:0 -S hdist:100000 -S isize:100000 -S aln_score:75 -S clipped_frac:0.50",
        ),
    shell:
        """
        PRE="temp/{wildcards.sm}/leviosam2/{wildcards.sm}"
        {LEVIO_EXE} lift \
            -t {threads} -T 100000 \
            -G {params.G} {params.S} \
            -C {input.levio_index} \
            -a {input.bam} \
            -p $PRE \
            -f {input.ref} -m \
            -O bam
        """


rule leviosam2_sorted:
    input:
        bam=rules.leviosam2.output.bam,
    output:
        bam="results/{sm}-leviosam2.bam",
        csi="results/{sm}-leviosam2.bam.csi",
    threads: MAX_THREADS
    resources:
        mem_mb=MAX_THREADS * 4 * 1024,
    conda:
        DEFAULT_ENV
    shell:
        """
        samtools sort \
            -@ {threads} -m 3G \
            -o {output.bam} --write-index \
            {input.bam}
        """
