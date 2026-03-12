# Tool Overview

rseqc-redux provides 33 command-line tools for RNA-seq quality control. All tools follow the same pattern:

```bash
tool_name.py [options]
```

## All tools

### BAM Statistics

| Tool | Description |
|------|-------------|
| [bam_stat.py](bam_stat.md) | Summarize mapping statistics of a BAM file |
| [bam2fq.py](bam2fq.md) | Convert BAM alignments to FASTQ format |
| [bam2wig.py](bam2wig.md) | Convert BAM file to wiggle format |
| [divide_bam.py](divide_bam.md) | Equally divide BAM file into n parts |
| [split_bam.py](split_bam.md) | Split BAM by gene list |
| [split_paired_bam.py](split_paired_bam.md) | Split paired-end BAM into two single-end BAMs |

### Read Quality

| Tool | Description |
|------|-------------|
| [read_quality.py](read_quality.md) | Phred quality score distribution across read positions |
| [read_GC.py](read_GC.md) | GC content distribution of reads |
| [read_NVC.py](read_NVC.md) | Nucleotide frequency per read position (NVC plot) |
| [read_hexamer.py](read_hexamer.md) | Hexamer frequency analysis |
| [read_duplication.py](read_duplication.md) | Sequence- and mapping-based duplication rates |
| [read_distribution.py](read_distribution.md) | Read distribution over genomic features |

### Alignment Profiles

| Tool | Description |
|------|-------------|
| [clipping_profile.py](clipping_profile.md) | Soft-clipping profile across read positions |
| [deletion_profile.py](deletion_profile.md) | Deletion distribution across read positions |
| [insertion_profile.py](insertion_profile.md) | Insertion distribution across read positions |
| [mismatch_profile.py](mismatch_profile.md) | Mismatch distribution across read positions |

### Gene Body & Expression

| Tool | Description |
|------|-------------|
| [geneBody_coverage.py](geneBody_coverage.md) | RNA-seq read coverage over gene body |
| [geneBody_coverage2.py](geneBody_coverage2.md) | Gene body coverage from BigWig input |
| [FPKM_count.py](FPKM_count.md) | Raw count, FPM, and FPKM per gene |
| [FPKM-UQ.py](FPKM_UQ.md) | Count, FPKM, and FPKM-UQ per gene |
| [RPKM_saturation.py](RPKM_saturation.md) | RPKM saturation analysis |
| [tin.py](tin.md) | Transcript integrity number |

### Junction Analysis

| Tool | Description |
|------|-------------|
| [junction_annotation.py](junction_annotation.md) | Annotate splice junctions against gene model |
| [junction_saturation.py](junction_saturation.md) | Junction discovery saturation analysis |

### Fragment & Strandedness

| Tool | Description |
|------|-------------|
| [infer_experiment.py](infer_experiment.md) | Infer strandedness of RNA-seq experiment |
| [inner_distance.py](inner_distance.md) | Inner distance (insert size) of RNA-seq fragments |
| [RNA_fragment_size.py](RNA_fragment_size.md) | Fragment size statistics per gene |

### BigWig Utilities

| Tool | Description |
|------|-------------|
| [normalize_bigwig.py](normalize_bigwig.md) | Normalize BigWig signal to fixed wigsum |
| [overlay_bigwig.py](overlay_bigwig.md) | Pairwise operations on two BigWig files |

### Single Cell

| Tool | Description |
|------|-------------|
| [sc_bamStat.py](sc_bamStat.md) | Single-cell RNA-seq mapping statistics |
| [sc_editMatrix.py](sc_editMatrix.md) | Barcode/UMI error correction heatmaps |
| [sc_seqLogo.py](sc_seqLogo.md) | DNA sequence logo from FASTQ/FASTA |
| [sc_seqQual.py](sc_seqQual.md) | Sequencing quality heatmap from FASTQ |
