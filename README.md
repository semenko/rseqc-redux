# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Cover |   Missing |
|-------------------------------- | -------: | -------: | ------: | --------: |
| rseqc/BED.py                    |      213 |       18 |     92% |17, 20-21, 24-25, 35, 37, 39, 109, 117-119, 136, 138, 180-182, 192 |
| rseqc/FrameKmer.py              |       45 |        0 |    100% |           |
| rseqc/SAM.py                    |     1435 |      337 |     77% |90, 120, 205, 218, 229, 308, 360-361, 378-379, 411, 416, 422, 507, 531, 636, 672, 719, 725, 768, 775, 854, 921, 928-929, 931-1016, 1029, 1053, 1058-1067, 1095-1193, 1245, 1249, 1256, 1264-1265, 1273-1276, 1281, 1287-1292, 1299-1301, 1308-1310, 1319-1356, 1361-1362, 1396-1397, 1413, 1417-1418, 1422, 1439, 1453-1454, 1457, 1459, 1466-1467, 1516-1517, 1519-1520, 1526-1527, 1580-1581, 1591, 1594-1595, 1600, 1621, 1627-1628, 1638, 1660, 1678, 1737-1738, 1755, 1758, 1790-1799, 1826, 1839-1841, 1843-1844, 1854, 1859-1879, 1898-1903, 1909-1910, 1934, 1943-1944, 1952-1955, 1957, 1964, 1971, 1976-1990, 1999-2054, 2067, 2073-2074, 2084-2105, 2114-2115 |
| rseqc/\_\_init\_\_.py           |        1 |        0 |    100% |           |
| rseqc/bam\_cigar.py             |       50 |        0 |    100% |           |
| rseqc/cli\_common.py            |       55 |        2 |     96% |     69-70 |
| rseqc/fastq.py                  |      134 |       17 |     87% |69, 110, 112-113, 145, 152, 154-155, 282-283, 288-290, 303-304, 319-320 |
| rseqc/getBamFiles.py            |       42 |        2 |     95% |     57-58 |
| rseqc/heatmap.py                |       24 |        1 |     96% |        57 |
| rseqc/ireader.py                |       19 |        1 |     95% |        21 |
| rseqc/mystat.py                 |       17 |        0 |    100% |           |
| rseqc/scbam.py                  |      261 |       56 |     79% |189, 257-258, 311-352, 357, 360, 410-448 |
| rseqc/twoList.py                |       29 |        0 |    100% |           |
| scripts/FPKM\_UQ.py             |      117 |       52 |     56% |43-80, 175-181, 190-234, 243 |
| scripts/FPKM\_count.py          |      217 |      165 |     24% |113-115, 118-119, 124-156, 158-166, 168, 170-176, 178-396, 400 |
| scripts/RNA\_fragment\_size.py  |       89 |       15 |     83% |67, 69, 71, 73, 75, 80, 82, 116-117, 119-121, 125-126, 141 |
| scripts/RPKM\_saturation.py     |       93 |       14 |     85% |62-63, 155-156, 158-159, 161-162, 164-165, 168-182, 186 |
| scripts/\_\_init\_\_.py         |        0 |        0 |    100% |           |
| scripts/bam2fq.py               |       36 |       19 |     47% | 38-56, 60 |
| scripts/bam2wig.py              |       38 |       12 |     68% |90-91, 94, 102-108, 111, 123 |
| scripts/bam\_stat.py            |       16 |        2 |     88% |    27, 31 |
| scripts/clipping\_profile.py    |       23 |        9 |     61% |42-43, 47-53, 57 |
| scripts/deletion\_profile.py    |       28 |       11 |     61% |40-41, 45-47, 50-52, 55-62, 66 |
| scripts/divide\_bam.py          |       36 |        2 |     94% |    51, 64 |
| scripts/geneBody\_coverage2.py  |       95 |       18 |     81% |23-24, 39, 50, 56-58, 70-71, 73, 135-137, 140-141, 149-150, 154 |
| scripts/geneBody\_coverage.py   |      191 |       41 |     79% |64-65, 88-90, 119, 131-149, 225-235, 285-286, 306-314, 322-326, 332, 340 |
| scripts/infer\_experiment.py    |       38 |       18 |     53% |40-41, 43, 45-63, 67 |
| scripts/inner\_distance.py      |       26 |        7 |     73% |70-71, 74-75, 77-86, 90 |
| scripts/insertion\_profile.py   |       23 |        9 |     61% |41-42, 46-52, 56 |
| scripts/junction\_annotation.py |       98 |       17 |     83% |37, 53, 102, 114, 200-216, 220 |
| scripts/junction\_saturation.py |       36 |       13 |     64% |83-84, 86-87, 89-90, 92-93, 95-96, 99-109, 113 |
| scripts/mismatch\_profile.py    |       28 |       11 |     61% |40-41, 45-47, 50-52, 55-61, 65 |
| scripts/normalize\_bigwig.py    |      110 |       41 |     63% |51-52, 58-59, 69-85, 96, 100, 106-108, 118-130, 141, 157-158, 162 |
| scripts/overlay\_bigwig.py      |       55 |       10 |     82% |56-57, 80-81, 84-85, 87, 89, 91, 103 |
| scripts/read\_GC.py             |       18 |        3 |     83% | 29-30, 34 |
| scripts/read\_NVC.py            |       19 |        5 |     74% |32-33, 36-37, 41 |
| scripts/read\_distribution.py   |      209 |       81 |     61% |237-313, 317-380, 384 |
| scripts/read\_duplication.py    |       19 |        3 |     84% | 48-49, 53 |
| scripts/read\_hexamer.py        |       54 |       15 |     72% |55-56, 66-74, 77-85, 95-96, 100 |
| scripts/read\_quality.py        |       19 |        5 |     74% |41-42, 45-46, 50 |
| scripts/sc\_bamStat.py          |       30 |        5 |     83% |86, 99-100, 112, 116 |
| scripts/sc\_editMatrix.py       |       38 |       10 |     74% |108, 118-119, 131-148, 163 |
| scripts/sc\_seqLogo.py          |       46 |       12 |     74% |123, 133-134, 136-138, 140-141, 143-144, 149, 168 |
| scripts/sc\_seqQual.py          |       32 |        4 |     88% |87, 97-98, 135 |
| scripts/split\_bam.py           |       52 |        3 |     94% |93-94, 105 |
| scripts/split\_paired\_bam.py   |       58 |        2 |     97% |    75, 94 |
| scripts/tin.py                  |      210 |       38 |     82% |35, 39, 80, 82, 87, 89, 93, 103-104, 125-128, 161, 163, 168, 170, 190, 200-210, 287, 290-291, 293, 305-306, 338-340, 359 |
| **TOTAL**                       | **4522** | **1106** | **76%** |           |


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