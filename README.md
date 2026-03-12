# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                  |    Stmts |     Miss |   Cover |   Missing |
|---------------------- | -------: | -------: | ------: | --------: |
| rseqc/BED.py          |     2099 |     1817 |     13% |33-259, 273-308, 315-361, 365-442, 452, 454, 456, 528, 530, 532, 540-542, 562, 564, 578-665, 675, 677, 679, 714-716, 730, 732, 734, 755-757, 765, 792-914, 918-958, 965-1026, 1032-1116, 1121-1160, 1177-1500, 1507-1592, 1596-1650, 1659-1664, 1672-1746, 1756-1938, 1945-2136, 2144-2248, 2256-2435, 2443-2535, 2542-2634, 2640-2733, 2743-2896 |
| rseqc/FrameKmer.py    |       71 |       44 |     38% |21-35, 48-64, 71-76, 87, 90-97 |
| rseqc/SAM.py          |     1936 |     1318 |     32% |42-45, 77, 154, 177, 194, 216, 228-230, 239, 259-263, 267-271, 277, 281, 286-290, 296, 301, 306-310, 316, 322, 326-376, 381-429, 462, 471, 505-881, 886-887, 905, 983-1007, 1021, 1060, 1100-1101, 1113, 1119, 1153-1155, 1173, 1180, 1261, 1273, 1275, 1282-1287, 1314-1396, 1411, 1423, 1425, 1432-1437, 1464-1548, 1554-1555, 1567, 1571, 1596, 1609-1611, 1619, 1647-1834, 1842-2023, 2030-2061, 2068-2218, 2233-2433, 2438-2619, 2634-2790, 2798-2875, 2879-2885 |
| rseqc/\_\_init\_\_.py |        1 |        0 |    100% |           |
| rseqc/annoGene.py     |      243 |       83 |     66% |30, 32, 48, 50, 52, 89, 91, 93, 104-106, 122, 124, 126, 139-141, 156, 158, 160, 162, 197, 199, 201, 203, 218, 220, 231-250, 255-310 |
| rseqc/bam\_cigar.py   |      146 |       34 |     77% |59, 74, 76, 80-84, 99, 101, 103, 108, 123, 127-132, 147, 149, 152-155, 170, 174-179, 201-208 |
| rseqc/changePoint.py  |        9 |        0 |    100% |           |
| rseqc/cigar.py        |      110 |        4 |     96% |69, 93, 119, 149 |
| rseqc/dotProduct.py   |       28 |        4 |     86% |     11-14 |
| rseqc/fasta.py        |      165 |      109 |     34% |43-44, 62-65, 73-76, 82-97, 101-114, 120-123, 131-141, 146-199, 207-224, 228-233, 237-241, 245 |
| rseqc/fastq.py        |      108 |       53 |     51% |65, 95-119, 145, 152, 154-155, 212-258 |
| rseqc/getBamFiles.py  |       45 |       13 |     71% |46, 53, 56-64, 67-68, 73 |
| rseqc/heatmap.py      |       26 |        4 |     85% |57, 133-135 |
| rseqc/ireader.py      |       19 |        4 |     79% |     18-21 |
| rseqc/mystat.py       |      104 |       16 |     85% |46-63, 90-91, 103-104 |
| rseqc/orf.py          |      158 |      103 |     35% |25-26, 30-31, 37-38, 42-43, 85-196 |
| rseqc/quantile.py     |       19 |        3 |     84% |31, 50, 52 |
| rseqc/scbam.py        |      376 |      334 |     11% |17-20, 92-192, 221-423, 441-484, 535-651 |
| rseqc/twoList.py      |       33 |        6 |     82% |35-36, 53-54, 59-60 |
| rseqc/wiggle.py       |      104 |       39 |     62% |91, 113-114, 119-120, 126-127, 131-132, 138-140, 144-146, 150-152, 156-158, 164-181 |
| **TOTAL**             | **5800** | **3988** | **31%** |           |


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