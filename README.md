# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Cover |   Missing |
|-------------------------------- | -------: | -------: | ------: | --------: |
| rseqc/BED.py                    |      213 |       18 |     92% |17, 20-21, 24-25, 35, 37, 39, 109, 117-119, 136, 138, 180-182, 192 |
| rseqc/FrameKmer.py              |       45 |        0 |    100% |           |
| rseqc/SAM.py                    |     1431 |      669 |     53% |80, 82-84, 190, 199, 221, 229, 263, 266-270, 308-309, 315, 352-355, 366-371, 375, 450, 596-620, 630, 666, 719, 762, 769, 848, 914-915, 921, 931-1010, 1023, 1047, 1052-1061, 1089-1187, 1201-1380, 1387-1557, 1572-1715, 1731-1732, 1749, 1752, 1784-1793, 1820, 1833-1835, 1837-1838, 1848, 1853-1873, 1892-1897, 1903-1904, 1923-2048, 2056-2120 |
| rseqc/\_\_init\_\_.py           |        1 |        0 |    100% |           |
| rseqc/bam\_cigar.py             |       50 |        0 |    100% |           |
| rseqc/cli\_common.py            |       35 |        4 |     89% |     48-51 |
| rseqc/fastq.py                  |      134 |       51 |     62% |69, 99-116, 145, 152, 154-155, 269-321 |
| rseqc/getBamFiles.py            |       42 |        2 |     95% |     57-58 |
| rseqc/heatmap.py                |       24 |        1 |     96% |        57 |
| rseqc/ireader.py                |       19 |        1 |     95% |        21 |
| rseqc/mystat.py                 |       17 |        0 |    100% |           |
| rseqc/scbam.py                  |      261 |      114 |     56% |187, 247-352, 357-362, 367, 371-376, 382-387, 390, 393-398, 410-446 |
| rseqc/twoList.py                |       29 |        0 |    100% |           |
| scripts/FPKM\_UQ.py             |      119 |       53 |     55% |44-81, 176-182, 191-236, 245 |
| scripts/FPKM\_count.py          |      218 |      188 |     14% |48-416, 420 |
| scripts/RNA\_fragment\_size.py  |       90 |       29 |     68% |67, 69, 71, 73, 75, 80, 82, 96-149, 153 |
| scripts/RPKM\_saturation.py     |       98 |       38 |     61% |46, 56-57, 82-205, 209 |
| scripts/\_\_init\_\_.py         |        0 |        0 |    100% |           |
| scripts/bam2fq.py               |       40 |       33 |     18% | 15-71, 75 |
| scripts/bam2wig.py              |       43 |       36 |     16% |16-123, 135 |
| scripts/bam\_stat.py            |       20 |       14 |     30% | 14-40, 44 |
| scripts/clipping\_profile.py    |       29 |       22 |     24% | 17-61, 65 |
| scripts/deletion\_profile.py    |       35 |       28 |     20% | 15-81, 85 |
| scripts/divide\_bam.py          |       40 |       32 |     20% | 18-64, 68 |
| scripts/geneBody\_coverage2.py  |       99 |       88 |     11% |23-112, 116-164, 168 |
| scripts/geneBody\_coverage.py   |      214 |       80 |     63% |66-67, 90-92, 122, 126-127, 134-135, 141, 143, 145, 147, 158, 233-244, 254-358, 362 |
| scripts/infer\_experiment.py    |       39 |       33 |     15% | 12-73, 77 |
| scripts/inner\_distance.py      |       33 |       26 |     21% |24-101, 105 |
| scripts/insertion\_profile.py   |       29 |       22 |     24% | 16-60, 64 |
| scripts/junction\_annotation.py |      106 |       35 |     67% |32, 48, 97, 109, 158-232, 236 |
| scripts/junction\_saturation.py |       41 |       34 |     17% |17-128, 132 |
| scripts/mismatch\_profile.py    |       35 |       28 |     20% | 15-80, 84 |
| scripts/normalize\_bigwig.py    |      111 |      102 |      8% |16-164, 168 |
| scripts/overlay\_bigwig.py      |       56 |       49 |     12% |14-100, 104 |
| scripts/read\_GC.py             |       23 |       16 |     30% | 13-44, 48 |
| scripts/read\_NVC.py            |       24 |       17 |     29% | 13-56, 60 |
| scripts/read\_distribution.py   |      217 |      112 |     48% |187-398, 402 |
| scripts/read\_duplication.py    |       24 |       17 |     29% | 17-57, 61 |
| scripts/read\_hexamer.py        |       55 |       43 |     22% |20-97, 101 |
| scripts/read\_quality.py        |       24 |       17 |     29% | 13-65, 69 |
| scripts/sc\_bamStat.py          |       31 |       25 |     19% |15-113, 117 |
| scripts/sc\_editMatrix.py       |       39 |       33 |     15% |17-149, 164 |
| scripts/sc\_seqLogo.py          |       47 |       41 |     13% |16-153, 169 |
| scripts/sc\_seqQual.py          |       33 |       27 |     18% |14-122, 136 |
| scripts/split\_bam.py           |       80 |       71 |     11% |16-133, 137 |
| scripts/split\_paired\_bam.py   |       62 |       55 |     11% |14-96, 100 |
| scripts/tin.py                  |      226 |       78 |     65% |80, 82, 87, 89, 93, 103-104, 125-128, 161, 163, 168, 170, 197-198, 202, 204, 206, 208, 221, 242-371, 375 |
| **TOTAL**                       | **4681** | **2382** | **49%** |           |


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