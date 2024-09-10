#
# Simple function that allows the pipeline to use fastqs or already aligned short paired reads
#
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
            samtools collate -f -@ {threads} -u -O {input.reads} -T $(dirname {output.fastq}) \
                | samtools fastq -@ {threads} -0 /dev/null \
                | bgzip -@ {threads} \
                > {output.fastq}
        else
            ln -s {input.reads} {output.fastq}
        fi
        """

#
# Alignment of fastqs to the 6Gbp donor specific assembly (DSA)
# These outputs will have many mapq 0 reads that will later be adjusted. 
# But this is only becuase of the diploid setup, not a lack of confidence in the alignment.
#
rule bwa_mem2:
    input:
        ref=DSA,
        fastq=rules.fastq_input.output.fastq,
    output:
        cram="results/DSA/{sm}-bwamem2.cram",
        csi="results/DSA/{sm}-bwamem2.cram.csi",
    threads: MAX_THREADS
    resources:
        mem_mb=MAX_THREADS * 4 * 1024,
    conda:
        DEFAULT_ENV
    params:
        bwa_h=config.get("bwa_h", "0"),
    shell:
        """
        bwa-mem2 mem \
            -t {threads} -h {params.bwa_h} \
            -p {input.ref} {input.fastq} \
            | samtools sort \
                -@ {threads} -m 1G \
                -O CRAM -T {input.ref} --output-fmt-option embed_ref=1 \
                -o {output.cram} --write-index
        """

#
# Index the chain file for leviosam2. 
#
# Input here is a chain file that defines the alignment between the DSA and the reference genome
# at a contig level (>100 kbp of alignment).
# 
# The output is a special leviosam2 index file that is used to lift over the alignments from the DSA to the reference genome.
#
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
            -p results/leviosam2-index/index \
            -c {input.chain} \
            -F {input.fai}
        """

#
# Lift over the alignments from the DSA to the reference genome using the chain file / leviosam2 index.
#
# This is not a realignment, but a lift over of the reads from the DSA to the reference genome.
#
rule leviosam2:
    input:
        cram=rules.bwa_mem2.output.cram,
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
            -a {input.cram} \
            -p $PRE \
            -f {input.ref} -m \
            -O bam
        """


#
# This step sorted the leviosam2 output and fixes some tags in the CRAM file.
#
# Specifically, the MAPQ is reset to 60 for all reads that were previously aligned to the DSA.
# And the XS tag is set to zero for all reads that were aligned to the DSA.
# This is a hueristic that we may need to return to in the future.
#
# Other tags and fields like CIGAR, bitflags, and MD are correctly updated by 
# leviosam2 during liftover.
#
rule leviosam2_sorted:
    input:
        bam=rules.leviosam2.output.bam,
        ref=REF,
    output:
        cram="results/{sm}-leviosam2.cram",
        csi="results/{sm}-leviosam2.cram.csi",
    threads: SORT_THREADS
    resources:
        mem_mb=SORT_THREADS * 4 * 1024,
    conda:
        DEFAULT_ENV
    params:
        reset_mapq=workflow.source_path("../scripts/reset-mapq.py"),
    shell:
        """
        python {params.reset_mapq} -t {threads} {input.bam} \
            | samtools sort \
                -@ {threads} -m 3G \
                -O CRAM -T {input.ref} --output-fmt-option embed_ref=1 \
                -o {output.cram} --write-index 
        """
