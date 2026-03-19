# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Cover |   Missing |
|-------------------------------- | -------: | -------: | ------: | --------: |
| rseqc/BED.py                    |      213 |       18 |     92% |17, 20-21, 24-25, 35, 37, 39, 111, 119-121, 140, 142, 179-181, 191 |
| rseqc/FrameKmer.py              |       45 |        0 |    100% |           |
| rseqc/SAM.py                    |     1467 |      110 |     93% |269, 667, 763, 849, 1021, 1095, 1097, 1234, 1238, 1245, 1253, 1261-1262, 1267, 1273-1281, 1288-1290, 1297-1299, 1308-1339, 1344-1345, 1382-1383, 1399, 1403-1404, 1425, 1439-1440, 1443, 1452-1453, 1496-1497, 1506-1507, 1545-1546, 1556, 1559-1560, 1586, 1592-1593, 1603, 1626, 1699-1700, 1717, 1720, 1788, 1801-1803, 1805-1806, 1816, 1863-1864, 1888, 1897-1898, 1909, 1918, 1925, 1937, 1943, 2018, 2024-2025, 2037, 2044, 2046, 2048, 2055 |
| rseqc/\_\_init\_\_.py           |        1 |        0 |    100% |           |
| rseqc/bam\_cigar.py             |       50 |        0 |    100% |           |
| rseqc/cli\_common.py            |      139 |        4 |     97% |84-85, 136-137 |
| rseqc/fastq.py                  |      145 |       17 |     88% |85, 126, 128-129, 161, 168, 170-171, 298-299, 304-306, 319-320, 335-336 |
| rseqc/heatmap.py                |       24 |        1 |     96% |        57 |
| rseqc/mystat.py                 |       17 |        0 |    100% |           |
| rseqc/scbam.py                  |      277 |       81 |     71% |173, 243-245, 299-347, 350, 390, 393-398, 401-402, 404, 406, 409-446 |
| scripts/FPKM\_UQ.py             |      112 |       47 |     58% |43-75, 170-176, 185-229, 238 |
| scripts/FPKM\_count.py          |      180 |      142 |     21% |103-124, 126-166, 171-175, 181-299, 307-308, 312 |
| scripts/RNA\_fragment\_size.py  |       73 |       10 |     86% |64, 66, 68, 70, 72, 77, 79, 113-114, 125 |
| scripts/RPKM\_saturation.py     |       97 |       14 |     86% |62-63, 159-160, 162-163, 165-166, 168-169, 172-186, 190 |
| scripts/\_\_init\_\_.py         |        0 |        0 |    100% |           |
| scripts/bam2fq.py               |       36 |       19 |     47% | 38-56, 60 |
| scripts/bam2wig.py              |       35 |       10 |     71% |91, 99-105, 108, 120 |
| scripts/bam\_stat.py            |       16 |        2 |     88% |    27, 31 |
| scripts/clipping\_profile.py    |       23 |        7 |     70% | 47-53, 57 |
| scripts/deletion\_profile.py    |       28 |        9 |     68% |45-47, 50-52, 55-62, 66 |
| scripts/divide\_bam.py          |       37 |        2 |     95% |    54, 67 |
| scripts/geneBody\_coverage2.py  |       80 |       12 |     85% |29-30, 47, 59-60, 62, 126-128, 131-132, 141 |
| scripts/geneBody\_coverage.py   |      172 |       40 |     77% |61-62, 99, 111-129, 201-209, 258-259, 279-287, 295-299, 305, 313 |
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
| scripts/read\_hexamer.py        |       48 |       10 |     79% |47-48, 69-77, 87-88, 92 |
| scripts/read\_quality.py        |       19 |        5 |     74% |41-42, 45-46, 50 |
| scripts/sc\_bamStat.py          |       30 |        5 |     83% |87, 99-100, 112, 116 |
| scripts/sc\_editMatrix.py       |       38 |       10 |     74% |108, 118-119, 131-148, 163 |
| scripts/sc\_seqLogo.py          |       46 |       12 |     74% |123, 133-134, 136-138, 140-141, 143-144, 149, 168 |
| scripts/sc\_seqQual.py          |       32 |        4 |     88% |87, 97-98, 135 |
| scripts/split\_bam.py           |       51 |        3 |     94% |93-94, 105 |
| scripts/split\_paired\_bam.py   |       55 |        2 |     96% |    69, 88 |
| scripts/tin.py                  |      176 |       35 |     80% |33, 37, 63, 65, 70, 72, 76, 86-87, 93, 133, 135, 140, 142, 162, 172-182, 259, 262-263, 265, 277-278, 310-312, 330 |
| **TOTAL**                       | **4431** |  **811** | **82%** |           |


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