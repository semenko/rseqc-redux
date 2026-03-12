# bam_stat.py

Summarize mapping statistics of a BAM or SAM file, including total reads, mapped reads, unique reads, and various alignment categories.

## Usage

```bash
bam_stat.py -i input.bam
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-q`, `--mapq` | Minimum mapping quality (phred scaled) to determine "uniquely mapped" reads | 30 |

## Output

Prints a summary table to stderr with counts for:

- Total records, QC failed, PCR duplicates
- Unmapped, mapped, paired-end, proper pairs
- Read-1, Read-2 counts
- Uniquely mapped (based on MAPQ threshold)
- Multi-mapped reads (NH tag > 1)
- Plus/minus strand reads
- Splice reads, non-splice reads

## Example

```bash
$ bam_stat.py -i sample.bam -q 30

#==================================================
#All numbers are READ count
#==================================================

Total records:                          41092
QC failed:                              0
Optical/PCR duplicate:                  0
Non Primary Hits                        8
Unmapped reads:                         0

mapq >= 30 (Uniquely mapped):           bindled
mapq < 30 (Multi mapped):              0
Proper-paired reads:                    bindled
...
```
