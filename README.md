# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Cover |   Missing |
|-------------------------------- | -------: | -------: | ------: | --------: |
| rseqc/BED.py                    |      213 |       18 |     92% |17, 20-21, 24-25, 35, 37, 39, 109, 117-119, 136, 138, 180-182, 192 |
| rseqc/FrameKmer.py              |       45 |        0 |    100% |           |
| rseqc/SAM.py                    |     1431 |      669 |     53% |80, 82-84, 190, 199, 221, 229, 263, 266-270, 308-309, 315, 352-355, 366-371, 375, 450, 596-620, 630, 666, 719, 762, 769, 848, 914-915, 921, 931-1010, 1023, 1047, 1052-1061, 1089-1187, 1201-1380, 1387-1557, 1572-1715, 1731-1732, 1749, 1752, 1784-1793, 1820, 1833-1835, 1837-1838, 1848, 1853-1873, 1892-1897, 1903-1904, 1923-2048, 2056-2120 |
| rseqc/\_\_init\_\_.py           |        1 |        0 |    100% |           |
| rseqc/bam\_cigar.py             |       50 |        0 |    100% |           |
| rseqc/cli\_common.py            |       55 |       15 |     73% |19-21, 32, 43, 52, 62, 67-70, 108-111 |
| rseqc/fastq.py                  |      134 |       51 |     62% |69, 99-116, 145, 152, 154-155, 269-321 |
| rseqc/getBamFiles.py            |       42 |        2 |     95% |     57-58 |
| rseqc/heatmap.py                |       24 |        1 |     96% |        57 |
| rseqc/ireader.py                |       19 |        1 |     95% |        21 |
| rseqc/mystat.py                 |       17 |        0 |    100% |           |
| rseqc/scbam.py                  |      261 |      114 |     56% |187, 247-352, 357-362, 367, 371-376, 382-387, 390, 393-398, 410-446 |
| rseqc/twoList.py                |       29 |        0 |    100% |           |
| scripts/FPKM\_UQ.py             |      117 |       52 |     56% |43-80, 175-181, 190-234, 243 |
| scripts/FPKM\_count.py          |      217 |      187 |     14% |54-396, 400 |
| scripts/RNA\_fragment\_size.py  |       89 |       28 |     69% |67, 69, 71, 73, 75, 80, 82, 96-137, 141 |
| scripts/RPKM\_saturation.py     |       93 |       35 |     62% |52, 62-63, 88-182, 186 |
| scripts/\_\_init\_\_.py         |        0 |        0 |    100% |           |
| scripts/bam2fq.py               |       36 |       30 |     17% | 14-56, 60 |
| scripts/bam2wig.py              |       38 |       32 |     16% |15-111, 123 |
| scripts/bam\_stat.py            |       16 |       11 |     31% | 13-27, 31 |
| scripts/clipping\_profile.py    |       23 |       18 |     22% | 22-53, 57 |
| scripts/deletion\_profile.py    |       28 |       23 |     18% | 13-62, 66 |
| scripts/divide\_bam.py          |       36 |       29 |     19% | 17-60, 64 |
| scripts/geneBody\_coverage2.py  |       95 |       85 |     11% |22-111, 115-150, 154 |
| scripts/geneBody\_coverage.py   |      209 |       77 |     63% |64-65, 88-90, 120, 124-125, 132-133, 139, 141, 143, 145, 156, 231-242, 252-343, 347 |
| scripts/infer\_experiment.py    |       38 |       32 |     16% | 12-63, 67 |
| scripts/inner\_distance.py      |       26 |       21 |     19% | 29-86, 90 |
| scripts/insertion\_profile.py   |       23 |       18 |     22% | 21-52, 56 |
| scripts/junction\_annotation.py |       98 |       29 |     70% |37, 53, 102, 114, 163-216, 220 |
| scripts/junction\_saturation.py |       36 |       31 |     14% |22-109, 113 |
| scripts/mismatch\_profile.py    |       28 |       23 |     18% | 13-61, 65 |
| scripts/normalize\_bigwig.py    |      110 |      101 |      8% |16-158, 162 |
| scripts/overlay\_bigwig.py      |       55 |       48 |     13% |14-99, 103 |
| scripts/read\_GC.py             |       18 |       13 |     28% | 18-30, 34 |
| scripts/read\_NVC.py            |       19 |       14 |     26% | 18-37, 41 |
| scripts/read\_distribution.py   |      209 |      106 |     49% |185-380, 384 |
| scripts/read\_duplication.py    |       19 |       14 |     26% | 22-49, 53 |
| scripts/read\_hexamer.py        |       54 |       42 |     22% |20-96, 100 |
| scripts/read\_quality.py        |       19 |       14 |     26% | 18-46, 50 |
| scripts/sc\_bamStat.py          |       30 |       24 |     20% |15-112, 116 |
| scripts/sc\_editMatrix.py       |       38 |       32 |     16% |17-148, 163 |
| scripts/sc\_seqLogo.py          |       46 |       40 |     13% |16-152, 168 |
| scripts/sc\_seqQual.py          |       32 |       26 |     19% |14-121, 135 |
| scripts/split\_bam.py           |       72 |       65 |     10% |20-124, 128 |
| scripts/split\_paired\_bam.py   |       58 |       52 |     10% | 13-90, 94 |
| scripts/tin.py                  |      221 |       74 |     67% |79, 81, 86, 88, 92, 102-103, 124-127, 160, 162, 167, 169, 196-197, 201, 203, 205, 207, 220, 241-361, 365 |
| **TOTAL**                       | **4567** | **2297** | **50%** |           |


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