#!/usr/bin/env python
from pathlib import Path
from typing import Optional
import sys

import defopt
import pysam
from tqdm import tqdm


def main(
    infile: Optional[Path] = None,
    *,
    outfile: Path = "-",
    threads: int = 16,
    verbose: int = 0,
    mapq: int = 60,
    XS: int = 0,
):
    """
    Author Mitchell R. Vollger
    """
    if infile is None:
        infile = sys.stdin
    if outfile == "-":
        outfile = sys.stdout

    bam = pysam.AlignmentFile(infile, threads=threads / 2)
    obam = pysam.AlignmentFile(outfile, "wbu", template=bam, threads=threads / 2)

    for rec in tqdm(bam.fetch(until_eof=True)):
        rec.mapping_quality = mapq
        if rec.has_tag("XS"):
            rec.set_tag("XS", XS)
        obam.write(rec)

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
