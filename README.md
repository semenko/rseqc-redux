# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Cover |   Missing |
|-------------------------------- | -------: | -------: | ------: | --------: |
| rseqc/BED.py                    |      213 |       18 |     92% |17, 20-21, 24-25, 35, 37, 39, 109, 117-119, 136, 138, 180-182, 192 |
| rseqc/FrameKmer.py              |       45 |        0 |    100% |           |
| rseqc/SAM.py                    |     1431 |      329 |     77% |190, 199, 221, 229, 268, 270, 308, 315, 352, 369-370, 450, 630, 666, 719, 762, 769, 848, 914-915, 921, 931-1010, 1023, 1047, 1052-1061, 1089-1187, 1239, 1243, 1250, 1258-1259, 1267-1270, 1275, 1281-1286, 1293-1295, 1302-1304, 1313-1350, 1355-1356, 1390-1391, 1407, 1411-1412, 1416, 1433, 1447-1448, 1451, 1453, 1460-1461, 1510-1511, 1513-1514, 1520-1521, 1574-1575, 1585, 1588-1589, 1594, 1615, 1621-1622, 1632, 1654, 1672, 1731-1732, 1749, 1752, 1784-1793, 1820, 1833-1835, 1837-1838, 1848, 1853-1873, 1892-1897, 1903-1904, 1928, 1937-1938, 1946-1949, 1951, 1958, 1965, 1970-1984, 1993-2048, 2061, 2067-2068, 2078-2099, 2108-2109 |
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
| scripts/geneBody\_coverage.py   |      209 |       29 |     86% |64-65, 88-90, 120, 124-125, 132-133, 139, 141, 143, 145, 156, 231-242, 292-293, 311-312, 317, 347 |
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
| scripts/tin.py                  |      221 |       31 |     86% |79, 81, 86, 88, 92, 102-103, 124-127, 160, 162, 167, 169, 196-197, 203, 205, 207, 220, 293, 296-297, 299, 311-312, 344-346, 365 |
| **TOTAL**                       | **4547** | **1079** | **76%** |           |


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