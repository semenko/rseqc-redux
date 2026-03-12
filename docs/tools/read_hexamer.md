# read_hexamer.py

Calculate hexamer (6-mer) frequency from read sequences and optionally compare to a reference genome or transcriptome.

## Usage

```bash
read_hexamer.py -i reads.fq
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input` | Read sequences in FASTA or FASTQ format. Multiple files separated by `,` | Required |
| `-r`, `--refgenome` | Reference genome sequence in FASTA format | Optional |
| `-g`, `--refgene` | Reference mRNA sequences in FASTA format | Optional |

## Output

Prints tab-separated hexamer frequencies to stdout. Columns include the hexamer sequence and observed frequency, plus reference frequencies if provided.
