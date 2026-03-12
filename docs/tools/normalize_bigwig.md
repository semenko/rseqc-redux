# normalize_bigwig.py

Normalize a BigWig signal to a fixed wigsum (equivalent to a target total read coverage). Outputs a wiggle file.

## Usage

```bash
normalize_bigwig.py -i input.bw -o output.wig
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--bwfile` | Input BigWig file | Required |
| `-o`, `--output` | Output wiggle file | Required |
| `-t`, `--wigsum` | Target wigsum for normalization | 100000000 |
| `-r`, `--refgene` | Reference gene model in BED format (restricts normalization to gene regions) | Optional |
| `-c`, `--chunk` | Chromosome chunk size (bp) for processing | 500000 |
| `-f`, `--format` | Output format: `wig` or `bgr` (bedGraph) | bgr |

## Output

Normalized signal in wiggle or bedGraph format.
