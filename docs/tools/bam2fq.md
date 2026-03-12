# bam2fq.py

Convert alignments in BAM or SAM format into FASTQ format.

## Usage

```bash
bam2fq.py -i input.bam -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-o`, `--out-prefix` | Prefix of output FASTQ file(s) | Required |
| `-s`, `--single-end` | Treat as single-end sequencing | Off |
| `-c`, `--compress` | Compress output FASTQ file(s) using gzip | Off |

## Output

- **Paired-end** (default): `prefix_1.fq` and `prefix_2.fq`
- **Single-end** (`-s`): `prefix.fq`
- With `-c`: files are gzip-compressed (`.fq.gz`)
