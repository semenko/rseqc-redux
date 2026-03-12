# sc_seqQual.py

Generate a heatmap from a FASTQ file to visualize sequencing quality across read positions.

## Usage

```bash
sc_seqQual.py -i reads.fq -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--infile` | Input file in FASTQ format | Required |
| `-o`, `--outfile` | Prefix of output files | Required |
| `-n`, `--nseq-limit` | Maximum number of sequences to process | All |
| `--cell-width` | Cell width (points) in heatmap | 12 |
| `--cell-height` | Cell height (points) in heatmap | 10 |
| `--font-size` | Font size (points) | 6 |
| `--angle` | Column label angle (0, 45, 90, 270, 315) | 45 |
| `--text-color` | Color of cell numbers | black |
| `--file-type` | Output format: `pdf`, `png`, `tiff`, `bmp`, `jpeg` | pdf |
| `--no-num` | Omit numerical values from cells | Off |
| `--verbose` | Print detailed debugging information | Off |

## Output

Quality heatmap showing Phred scores per read position.
