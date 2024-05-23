rule fastq_input:
    input:
        reads=get_input,
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
