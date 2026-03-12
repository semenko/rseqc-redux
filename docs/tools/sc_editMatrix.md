# sc_editMatrix.py

Generate heatmaps to visualize error-corrected nucleotide changes in cell barcodes and UMIs. Shows the position (X-axis), type of edit (Y-axis, e.g., C→T), and frequency (color).

## Usage

```bash
sc_editMatrix.py -i input.bam -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--infile` | Input BAM file | Required |
| `-o`, `--outfile` | Prefix of output files | Required |
| `--limit` | Number of alignments to process | All |
| `--cr-tag` | BAM tag for raw cellular barcode | CR |
| `--cb-tag` | BAM tag for corrected cellular barcode | CB |
| `--ur-tag` | BAM tag for raw UMI | UR |
| `--ub-tag` | BAM tag for corrected UMI | UB |
| `--cell-width` | Cell width (points) in heatmap | 15 |
| `--cell-height` | Cell height (points) in heatmap | 10 |
| `--font-size` | Font size (points) | 8 |
| `--angle` | Column label angle (0, 45, 90, 270, 315) | 45 |
| `--text-color` | Color of cell numbers | black |
| `--file-type` | Output format: `pdf`, `png`, `tiff`, `bmp`, `jpeg` | pdf |
| `--verbose` | Print detailed running information | Off |
| `--no-num` | Omit numerical values from cells | Off |

## Output

Heatmap files showing barcode and UMI error correction patterns.
