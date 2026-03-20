# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Cover |   Missing |
|-------------------------------- | -------: | -------: | ------: | --------: |
| rseqc/BED.py                    |      201 |       15 |     93% |17, 20-21, 24-25, 35, 111-113, 132, 134, 167-169, 179 |
| rseqc/FrameKmer.py              |       39 |        5 |     87% |     28-32 |
| rseqc/SAM.py                    |     1476 |      112 |     92% |267, 464, 548, 755, 761, 806, 891, 965, 967, 1063, 1276, 1280, 1287, 1295, 1303-1304, 1309, 1315-1322, 1330-1332, 1339-1341, 1350-1381, 1386-1387, 1424-1425, 1441, 1445-1446, 1467, 1481-1482, 1485, 1494-1495, 1538-1539, 1548-1549, 1587-1588, 1598, 1601-1602, 1628, 1634-1635, 1645, 1668, 1741-1742, 1759, 1762, 1829, 1842-1844, 1846-1847, 1857, 1904-1905, 1929, 1938-1939, 1950, 1959, 1966, 1978, 1984, 2059, 2065-2066, 2078, 2085, 2087, 2089, 2096 |
| rseqc/\_\_init\_\_.py           |        1 |        0 |    100% |           |
| rseqc/bam\_cigar.py             |       33 |        0 |    100% |           |
| rseqc/cli\_common.py            |      144 |        4 |     97% |92-93, 144-145 |
| rseqc/fastq.py                  |      121 |       21 |     83% |40-41, 63-66, 97, 100, 132, 139, 141-142, 269-270, 275-277, 290-291, 306-307 |
| rseqc/heatmap.py                |       19 |        1 |     95% |        57 |
| rseqc/mystat.py                 |        8 |        0 |    100% |           |
| rseqc/scbam.py                  |      253 |        8 |     97% |162, 232-234, 285, 300, 334, 407 |
| scripts/FPKM\_UQ.py             |      112 |       47 |     58% |43-75, 170-176, 185-229, 238 |
| scripts/FPKM\_count.py          |      180 |      142 |     21% |103-124, 126-166, 171-175, 181-299, 307-308, 312 |
| scripts/RNA\_fragment\_size.py  |       79 |       10 |     87% |70, 72, 74, 76, 78, 83, 85, 119-120, 131 |
| scripts/RPKM\_saturation.py     |       97 |       14 |     86% |62-63, 159-160, 162-163, 165-166, 168-169, 172-186, 190 |
| scripts/\_\_init\_\_.py         |        0 |        0 |    100% |           |
| scripts/bam2fq.py               |       36 |       19 |     47% | 38-56, 60 |
| scripts/bam2wig.py              |       35 |       10 |     71% |91, 99-105, 108, 119 |
| scripts/bam\_stat.py            |       16 |        2 |     88% |    27, 31 |
| scripts/clipping\_profile.py    |       23 |        7 |     70% | 47-53, 57 |
| scripts/deletion\_profile.py    |       28 |        9 |     68% |45-47, 50-52, 55-62, 66 |
| scripts/divide\_bam.py          |       37 |        2 |     95% |    54, 67 |
| scripts/geneBody\_coverage2.py  |       80 |       12 |     85% |29-30, 47, 59-60, 62, 126-128, 131-132, 141 |
| scripts/geneBody\_coverage.py   |      193 |       47 |     76% |61-62, 115, 123, 135-159, 232-240, 289-290, 310-318, 326-330, 336, 344 |
| scripts/infer\_experiment.py    |       34 |       15 |     56% | 41-59, 63 |
| scripts/inner\_distance.py      |       26 |        7 |     73% |70-71, 74-75, 77-86, 90 |
| scripts/insertion\_profile.py   |       23 |        7 |     70% | 46-52, 56 |
| scripts/junction\_annotation.py |       98 |       17 |     83% |37, 53, 102, 114, 200-216, 220 |
| scripts/junction\_saturation.py |       36 |        5 |     86% |89-90, 99-109, 113 |
| scripts/mismatch\_profile.py    |       28 |        9 |     68% |45-47, 50-52, 55-61, 65 |
| scripts/normalize\_bigwig.py    |      112 |       15 |     87% |51-52, 58-59, 77, 83, 100, 106-108, 124, 139, 155-156, 162 |
| scripts/overlay\_bigwig.py      |       64 |       10 |     84% |78-79, 102-103, 106-107, 109, 111, 113, 127 |
| scripts/read\_GC.py             |       18 |        3 |     83% | 29-30, 34 |
| scripts/read\_NVC.py            |       19 |        5 |     74% |32-33, 36-37, 41 |
| scripts/read\_distribution.py   |      192 |       84 |     56% |206-307, 319, 323 |
| scripts/read\_duplication.py    |       19 |        3 |     84% | 48-49, 53 |
| scripts/read\_hexamer.py        |       48 |       27 |     44% |47-48, 53-88, 92 |
| scripts/read\_quality.py        |       19 |        5 |     74% |41-42, 45-46, 50 |
| scripts/sc\_bamStat.py          |       30 |        5 |     83% |87, 99-100, 112, 116 |
| scripts/sc\_editMatrix.py       |       38 |       10 |     74% |108, 118-119, 131-148, 163 |
| scripts/sc\_seqLogo.py          |       46 |       13 |     72% |133-134, 136-138, 140-141, 143-144, 149, 151-152, 168 |
| scripts/sc\_seqQual.py          |       32 |        5 |     84% |87, 97-98, 104-112 |
| scripts/split\_bam.py           |       51 |        3 |     94% |93-94, 105 |
| scripts/split\_paired\_bam.py   |       55 |        2 |     96% |    69, 88 |
| scripts/tin.py                  |      228 |       53 |     77% |34, 38, 64, 66, 71, 73, 77, 87-88, 94, 134, 136, 141, 143, 163, 173-183, 260, 263-264, 266, 278-279, 308-310, 330, 337, 345-361, 392 |
| **TOTAL**                       | **4427** |  **790** | **82%** |           |


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