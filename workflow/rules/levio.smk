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
        runtime=4 * 60,
    conda:
        DEFAULT_ENV
    shell:
        """
        # check if input is bam
        if [[ {input.reads} =~ .*\\.(bam|sam|cram) ]]; then
            samtools view -u -F 2304 -@ {threads} {input.reads} \
                | samtools collate -@ {threads} -u -O -T $(dirname {output.fastq}) \
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
        crai="results/DSA/{sm}-bwamem2.cram.crai",
    threads: MAX_THREADS
    resources:
        mem_mb=MAX_THREADS * 8 * 1024,
        runtime=16 * 60,
    conda:
        DEFAULT_ENV
    params:
        bwa_h=config.get("bwa_h", "0"),
    shell:
        """
        bwa-mem2 mem \
            -t {threads} -h {params.bwa_h} \
            -p {input.ref} {input.fastq}\
            | samtools sort \
                -@ {threads} -m 1G \
                -O CRAM --reference {input.ref} --output-fmt-option embed_ref=1 \
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
        mem_mb=64 * 1024,
        runtime=16 * 60,
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
        crai=rules.bwa_mem2.output.crai,
        levio_index=rules.leviosam2_index.output.index,
        ref=REF,
    output:
        lifted=temp("temp/{sm}/leviosam2/{sm}-committed.bam"),
        deferred=temp("temp/{sm}/leviosam2/{sm}-deferred.bam"),
        unliftable=temp("temp/{sm}/leviosam2/{sm}-unliftable.bam"),
    threads: MAX_THREADS
    resources:
        mem_mb=MAX_THREADS * 4 * 1024,
        runtime=16 * 60,
    conda:
        DEFAULT_ENV
    params:
        # maximum number of CIGAR opts to change, also the max gap size that can be spanned
        G=config.get("levio_G", 100_000),
        # Using -S clipped_frac 0.05 means when a read has >5% clipped bases, it is deferred. A lower value is more stringent (by deferring more reads).
        # we have to include some amount of clipping to avoid a segfault in leviosam2? (I think)
        # aln_score is the minumum score before the alignment is lifted over
        S=config.get(
            "levio_S",
            f"-S mapq:0 -S hdist:{10_000} -S isize:{1_000} -S clipped_frac:0.25 -S aln_score:100",
        ),
        # number of reads per thread
        T=config.get("levio_T", 256),
    shell:
        """
        PRE="temp/{wildcards.sm}/leviosam2/{wildcards.sm}"
        {LEVIO_EXE} lift -t {threads} -a {input.cram} \
            -T {params.T} -G {params.G} {params.S} \
            -C {input.levio_index} -p $PRE -f {input.ref} -m -O bam
        """
        # ^ bam is the only option, no CRAM.
        #samtools view -@ {threads} -u {input.cram} \


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
        lifted=rules.leviosam2.output.lifted,
        ref=REF,
    output:
        cram="results/{sm}-leviosam2.cram",
        crai="results/{sm}-leviosam2.cram.crai",
    threads: SORT_THREADS
    resources:
        mem_mb=SORT_THREADS * 4 * 1024,
        runtime=16 * 60,
    conda:
        DEFAULT_ENV
    params:
        reset_mapq=workflow.source_path("../scripts/reset-mapq.py"),
    shell:
        """
        python {params.reset_mapq} -t {threads} {input.lifted} \
            | samtools sort \
                -@ {threads} -m 3G \
                -O CRAM --reference {input.ref} --output-fmt-option embed_ref=1 \
                -o {output.cram} --write-index 
        """
