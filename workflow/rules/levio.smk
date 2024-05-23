

rule fastq_input:
    input:
        reads=lambda wc: config["inputs"][wc.sm],
    output:
        fastq=temp("temp/{sm}/fastq_input/{chrom}.fastq.gz"),
    threads: 8
    resources:
        mem_mb=32 * 1024,
    conda:
        DEFAULT_ENV
    shell:
        """
        # check if input is bam
        if [[ {input.reads} =~ .*\.(bam|sam|cram) ]]; then
            samtools collate -@ {threads} -u -O {input.redads} -T tmp/collate \
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
        bam=temp("temp/{sm}/bwa_mem2/{chrom}.bam"),
        csi=temp("temp/{sm}/bwa_mem2/{chrom}.bam.csi"),
    threads: 16
    resources:
        mem_mb=32 * 1024,
    conda:
        DEFAULT_ENV
    shell:
        """
        bwa-mem2 mem -t {threads} -p {input.ref} {input.fastq} \
            | samtools sort -@ {threads} -m 4G -o {output.bam} --write-index
        """


rule leviosam2:
    input:
        bam=rules.bwa_mem2.output.bam,
        levio_index=LEVIO_INDEX,
        ref=REF,
    output:
        bam=temp("temp/{sm}/leviosam2/{sm}-committed.bam"),
    threads: 16
    resources:
        mem_mb=32 * 1024,
    conda:
        DEFAULT_ENV
    params:
        G=config.get("levio_G", 50),
        S=config.get("levio_S", "-S mapq:0 -S hdist:50 -S isize:10000"),
    shell:
        """
        PRE="temp/{wildcards.sm}/leviosam2/{wildcards.sm}"
        {LEVIO_EXE} lift \
            -t {threads} -T 100000 \
            -G {params.G} -S {params.S} \
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
    threads: 16
    resources:
        mem_mb=64 * 1024,
    conda:
        DEFAULT_ENV
    shell:
        """
        samtools sort \
            -@ {threads} -m 3G \ 
            -o {output.bam} --write-index \
            {input.bam}
        """
