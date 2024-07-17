rule fastq_input:
    input:
        reads=lambda wc: config["inputs"][wc.sm],
    output:
        fastq=temp("temp/{sm}/fastq_input/{sm}.fastq.gz"),
    threads: SORT_THREADS
    resources:
        mem_mb=SORT_THREADS * 4 * 1024,
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
        bwa_h=config.get("bwa_h", "0"),
    shell:
        """
        bwa-mem2 mem \
            -t {threads} -h {params.bwa_h} \
            -p {input.ref} {input.fastq} \
            | python {params.reset_mapq} \
            | samtools sort -@ {threads} -m 1G -o {output.bam} --write-index
        """


rule leviosam2_index:
    input:
        chain=LEVIO_CHAIN,
        fai=FAI,
    output:
        index="results/leviosam2-index/index.clft",
    conda:
        DEFAULT_ENV
    threads: 1
    resources:
        mem_mb=64*1024,
    shell:
        """
        {LEVIO_EXE} index \
            -p results/leviosam2-index/index.clft \
            -c {input.chain} \
            -F {input.fai}
        """

rule leviosam2:
    input:
        bam=rules.bwa_mem2.output.bam,
        csi=rules.bwa_mem2.output.csi,
        levio_index=rules.leviosam2_index.output.index,
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
        G=config.get("levio_G", 50_000),
        # "-S aln_score:100 -S clipped_frac:0.95",
        S=config.get(
            "levio_S",
            f"-S mapq:0 -S hdist:{50_000} -S isize:{50_000}",
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
    threads: SORT_THREADS
    resources:
        mem_mb=SORT_THREADS * 4 * 1024,
    conda:
        DEFAULT_ENV
    shell:
        """
        samtools sort \
            -@ {threads} -m 3G \
            -o {output.bam} --write-index \
            {input.bam}
        """
