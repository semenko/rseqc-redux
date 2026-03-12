# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                  |    Stmts |     Miss |   Cover |   Missing |
|---------------------- | -------: | -------: | ------: | --------: |
| rseqc/BED.py          |      248 |       18 |     93% |21, 24-25, 28-29, 39, 41, 43, 118, 126-128, 148, 150, 193-195, 205 |
| rseqc/FrameKmer.py    |       72 |       45 |     38% |21-36, 49-65, 72-77, 88, 91-98 |
| rseqc/SAM.py          |     1501 |      984 |     34% |43-44, 49-51, 81-82, 167-168, 198, 206, 220, 240, 245-247, 275, 281, 285-390, 395-443, 471, 480, 513-514, 530, 609-633, 645, 684, 724-725, 735, 741, 770-772, 787, 794, 875, 887, 889, 896-901, 929-1019, 1032, 1044, 1046, 1053-1058, 1086-1178, 1192-1379, 1387-1566, 1581-1732, 1747-1946, 1953-2108, 2116-2192 |
| rseqc/\_\_init\_\_.py |        1 |        0 |    100% |           |
| rseqc/annoGene.py     |      243 |       82 |     66% |34, 36, 53, 55, 57, 94, 96, 98, 107-109, 126, 128, 130, 142-144, 159, 161, 163, 165, 199, 201, 203, 205, 220, 222, 233-250, 255-311 |
| rseqc/bam\_cigar.py   |      146 |       34 |     77% |59, 74, 76, 80-84, 99, 101, 103, 108, 123, 127-132, 147, 149, 152-155, 170, 174-179, 201-208 |
| rseqc/changePoint.py  |        9 |        0 |    100% |           |
| rseqc/cigar.py        |      110 |        4 |     96% |69, 93, 119, 149 |
| rseqc/dotProduct.py   |       28 |        4 |     86% |     11-14 |
| rseqc/fasta.py        |      172 |      112 |     35% |48-49, 67-70, 78-81, 87-102, 106-119, 125-129, 137-147, 152-213, 228-247, 251-263, 267-271, 275 |
| rseqc/fastq.py        |      134 |       80 |     40% |69, 99-116, 145, 152, 154-155, 185-221, 269-321 |
| rseqc/getBamFiles.py  |       46 |       13 |     72% |46, 54, 57-65, 68-69, 74 |
| rseqc/heatmap.py      |       24 |        1 |     96% |        57 |
| rseqc/ireader.py      |       19 |        4 |     79% |     18-21 |
| rseqc/mystat.py       |      104 |       16 |     85% |46-63, 90-91, 103-104 |
| rseqc/orf.py          |      158 |      103 |     35% |25-26, 30-31, 37-38, 42-43, 85-196 |
| rseqc/quantile.py     |       19 |        3 |     84% |31, 50, 52 |
| rseqc/scbam.py        |      380 |      126 |     67% |183, 243-352, 357-361, 366, 373-376, 383-387, 396, 406-443, 501, 571, 581, 604, 620-623, 631-632, 638, 663, 669-670 |
| rseqc/twoList.py      |       33 |        6 |     82% |35-36, 53-54, 59-60 |
| rseqc/wiggle.py       |      102 |       36 |     65% |93, 115-116, 121-122, 128-129, 133-134, 140-142, 146-148, 152-154, 158-160, 166-180 |
| **TOTAL**             | **3549** | **1671** | **53%** |           |


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