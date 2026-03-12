# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Cover |   Missing |
|-------------------------------- | -------: | -------: | ------: | --------: |
| rseqc/BED.py                    |      239 |       18 |     92% |21, 24-25, 28-29, 39, 41, 43, 113, 121-123, 140, 142, 184-186, 196 |
| rseqc/FrameKmer.py              |       65 |        1 |     98% |        86 |
| rseqc/SAM.py                    |     1532 |      900 |     41% |44-45, 49-52, 77, 82-83, 217, 229, 240, 275, 285-286, 294-295, 301, 340-343, 351-353, 357, 367-379, 385, 387, 391-395, 398, 407, 410, 412-416, 418, 489, 502, 596-620, 632, 671, 711-712, 722, 728, 757-759, 774, 781, 862, 876, 878, 932-1032, 1045, 1059, 1061, 1073, 1078-1087, 1115-1217, 1231-1418, 1426-1605, 1620-1771, 1786-1985, 1992-2147, 2155-2231 |
| rseqc/\_\_init\_\_.py           |        1 |        0 |    100% |           |
| rseqc/annoGene.py               |      239 |       82 |     66% |32, 34, 51, 53, 55, 92, 94, 96, 105-107, 124, 126, 128, 140-142, 157, 159, 161, 163, 197, 199, 201, 203, 216, 218, 229-246, 251-307 |
| rseqc/bam\_cigar.py             |      136 |        0 |    100% |           |
| rseqc/changePoint.py            |        8 |        0 |    100% |           |
| rseqc/cigar.py                  |      107 |        4 |     96% |69, 93, 119, 149 |
| rseqc/cli\_common.py            |       35 |        4 |     89% |     48-51 |
| rseqc/fasta.py                  |      172 |       15 |     91% |129, 205-213, 245-246, 267-271, 275 |
| rseqc/fastq.py                  |      134 |       51 |     62% |69, 99-116, 145, 152, 154-155, 269-321 |
| rseqc/getBamFiles.py            |       44 |        3 |     93% | 57-58, 72 |
| rseqc/heatmap.py                |       24 |        1 |     96% |        57 |
| rseqc/ireader.py                |       19 |        1 |     95% |        21 |
| rseqc/mystat.py                 |      110 |        5 |     95% |59, 91-92, 105-106 |
| rseqc/orf.py                    |      158 |       27 |     83% |25-26, 30-31, 37-38, 42-43, 92-93, 97-100, 104-105, 109-112, 142, 159-160, 167-168, 181-182, 188-189 |
| rseqc/quantile.py               |       19 |        0 |    100% |           |
| rseqc/scbam.py                  |      379 |      127 |     66% |187, 247-352, 357-362, 367, 371-376, 382-387, 390, 393-398, 410-446, 503, 573, 583, 606, 622-625, 633-634, 640, 665, 671-672 |
| rseqc/twoList.py                |       32 |        0 |    100% |           |
| rseqc/wiggle.py                 |        1 |        0 |    100% |           |
| scripts/FPKM\_UQ.py             |      119 |      119 |      0% |    16-245 |
| scripts/FPKM\_count.py          |      218 |      188 |     14% |48-416, 420 |
| scripts/RNA\_fragment\_size.py  |       89 |       73 |     18% |33-88, 94-147, 151 |
| scripts/RPKM\_saturation.py     |       98 |       76 |     22% |32-78, 82-205, 209 |
| scripts/\_\_init\_\_.py         |        0 |        0 |    100% |           |
| scripts/bam2fq.py               |       40 |       40 |      0% |      6-75 |
| scripts/bam2wig.py              |       43 |       43 |      0% |     7-135 |
| scripts/bam\_stat.py            |       20 |       20 |      0% |      6-44 |
| scripts/clipping\_profile.py    |       29 |       29 |      0% |      8-65 |
| scripts/deletion\_profile.py    |       35 |       35 |      0% |      6-85 |
| scripts/divide\_bam.py          |       40 |       40 |      0% |      7-68 |
| scripts/geneBody\_coverage2.py  |       99 |       99 |      0% |     7-168 |
| scripts/geneBody\_coverage.py   |      213 |      177 |     17% |43-46, 64-103, 111-158, 163-249, 253-357, 361 |
| scripts/infer\_experiment.py    |       39 |       39 |      0% |      4-77 |
| scripts/inner\_distance.py      |       33 |       33 |      0% |    15-105 |
| scripts/insertion\_profile.py   |       29 |       29 |      0% |      7-64 |
| scripts/junction\_annotation.py |      106 |       35 |     67% |32, 48, 97, 109, 158-232, 236 |
| scripts/junction\_saturation.py |       41 |       41 |      0% |     8-132 |
| scripts/mismatch\_profile.py    |       35 |       35 |      0% |      6-84 |
| scripts/normalize\_bigwig.py    |      111 |      111 |      0% |     4-168 |
| scripts/overlay\_bigwig.py      |       56 |       56 |      0% |     4-104 |
| scripts/read\_GC.py             |       23 |       23 |      0% |      4-48 |
| scripts/read\_NVC.py            |       24 |       24 |      0% |      4-60 |
| scripts/read\_distribution.py   |      201 |      182 |      9% |36-123, 148-354, 358 |
| scripts/read\_duplication.py    |       24 |       24 |      0% |      8-61 |
| scripts/read\_hexamer.py        |       55 |       43 |     22% |20-97, 101 |
| scripts/read\_quality.py        |       24 |       24 |      0% |      4-69 |
| scripts/sc\_bamStat.py          |       31 |       31 |      0% |     7-117 |
| scripts/sc\_editMatrix.py       |       39 |       39 |      0% |     9-164 |
| scripts/sc\_seqLogo.py          |       47 |       47 |      0% |     8-169 |
| scripts/sc\_seqQual.py          |       33 |       33 |      0% |     6-136 |
| scripts/split\_bam.py           |       80 |       80 |      0% |     4-137 |
| scripts/split\_paired\_bam.py   |       62 |       62 |      0% |     4-100 |
| scripts/tin.py                  |      225 |      187 |     17% |47-55, 62-68, 75-93, 101-146, 154-176, 186-223, 241-370, 374 |
| **TOTAL**                       | **5815** | **3356** | **42%** |           |


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/semenko/rseqc-redux/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/semenko/rseqc-redux/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2Fsemenko%2Frseqc-redux%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.