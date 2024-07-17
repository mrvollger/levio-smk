
# Levio-smk

This is a Snakemake pipeline that aligns raw sequencing reads to a Donor-specific assembly (DSA) and then lifts those reads back to a reference genome using leviosam2 with custom parameters. 

Note that leviosam2 requires a chain file between the DSA and the reference genome and this chain file is not created by this pipeline.