# bam2wig.py

Convert a BAM file into wiggle format. The BAM file must be sorted and indexed. SAM format is not supported.

## Usage

```bash
bam2wig.py -i input.bam -s chrom.sizes -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM format (sorted, indexed) | Required |
| `-s`, `--chromSize` | Chromosome size file (tab-separated: chrom, size) | Required |
| `-o`, `--out-prefix` | Prefix of output wiggle file(s) | Required |
| `-t`, `--wigsum` | Normalize to this wigsum (e.g., 1000000000 for 10M 100nt reads) | None |
| `-u`, `--skip-multi-hits` | Skip non-uniquely mapped reads | Off |
| `-d`, `--strand` | Strand rule: "1++,1--,2+-,2-+" or "1+-,1-+,2++,2--" or "none" | None |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

Generates wiggle (`.wig`) file(s). If strand-specific, produces separate files for each strand.
