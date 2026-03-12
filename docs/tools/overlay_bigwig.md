# overlay_bigwig.py

Perform pairwise operations on two BigWig files. Both files must use the same reference genome.

## Usage

```bash
overlay_bigwig.py -i file1.bw -j file2.bw -a mean -o output.wig
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--bwfile1` | First BigWig file | Required |
| `-j`, `--bwfile2` | Second BigWig file | Required |
| `-a`, `--action` | Operation: `mean`, `sum`, `diff`, `ratio`, `max`, `min`, `geometric_mean` | Required |
| `-o`, `--output` | Output wiggle file | Required |
| `-c`, `--chunk` | Chromosome chunk size (bp) for processing | 100000 |

## Output

Wiggle file with the result of the pairwise operation.
