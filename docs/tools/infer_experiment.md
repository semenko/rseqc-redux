# infer_experiment.py

Infer the strandedness of an RNA-seq experiment from a BAM file. Determines whether the experiment is:

1. Single-end or paired-end
2. If strand-specific, how reads were stranded (sense vs antisense)

## Usage

```bash
infer_experiment.py -i input.bam -r gene_model.bed
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in SAM or BAM format | Required |
| `-r`, `--refgene` | Reference gene model in BED format | Required |
| `-s`, `--sample-size` | Number of reads to sample | 200000 |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

Prints strandedness fractions to stdout:

```
This is PairEnd Data
Fraction of reads failed to determine: 0.0172
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4903
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4925
```

- Values near 0.5/0.5 indicate unstranded
- Dominant fraction indicates the strand rule to use with other tools
