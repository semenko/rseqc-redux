# inner_distance.py

Calculate the inner distance (insert size) between read pairs for paired-end RNA-seq data. The inner distance is the gap between the end of read 1 and the start of read 2 on the transcript.

## Usage

```bash
inner_distance.py -i input.bam -r gene_model.bed -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-r`, `--refgene` | Reference gene model in BED format | Required |
| `-k`, `--sample-size` | Number of read pairs to use | 1000000 |
| `-l`, `--lower-bound` | Lower bound of inner distance (bp) | -250 |
| `-u`, `--upper-bound` | Upper bound of inner distance (bp) | 250 |
| `-s`, `--step` | Step size (bp) for histogram | 5 |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

- `prefix.inner_distance.txt` — inner distances per read pair
- `prefix.inner_distance_freq.txt` — inner distance histogram
- `prefix.inner_distance_plot.r` — R script for histogram
- `prefix.inner_distance_plot.pdf` — inner distance histogram (if R is available)

!!! tip
    Negative inner distances indicate overlapping read pairs, which is common for short fragments.
