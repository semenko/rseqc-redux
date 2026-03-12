# sc_seqLogo.py

Generate a DNA sequence logo from FASTA or FASTQ format files. Useful for visualizing nucleotide composition of sample barcodes, cell barcodes, and molecular barcodes (UMIs).

## Usage

```bash
sc_seqLogo.py -i barcodes.fq -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--infile` | Input file in FASTQ, FASTA, or plain sequence format (plain text or gzip/bzip2 compressed) | Required |
| `-o`, `--outfile` | Prefix of output files | Required |
| `--iformat` | Input format: `fq` or `fa` | fq |
| `--oformat` | Output format: `pdf`, `png`, or `svg` | pdf |
| `-n`, `--nseq-limit` | Maximum number of sequences to process | All |
| `--font-name` | Logo font: `sans`, `serif`, `mono`, etc. | sans |
| `--stack-order` | Stacking order: `big_on_top`, `small_on_top`, `fixed` | big_on_top |
| `--flip-below` | Flip characters below X-axis upside down | Off |
| `--shade-below` | Shading amount for below-axis characters (0.0–1.0) | 0.0 |
| `--fade-below` | Fading amount for below-axis characters (0.0–1.0) | 0.0 |
| `--excludeN` | Exclude sequences containing N | Off |
| `--highlight-start` | Start position for highlighted region | None |
| `--highlight-end` | End position for highlighted region | None |
| `--verbose` | Print detailed debugging information | Off |

## Output

Sequence logo file in the specified format.
