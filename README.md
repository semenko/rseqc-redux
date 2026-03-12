# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                  |    Stmts |     Miss |   Cover |   Missing |
|---------------------- | -------: | -------: | ------: | --------: |
| rseqc/BED.py          |      248 |       18 |     93% |21, 24-25, 28-29, 39, 41, 43, 118, 126-128, 148, 150, 193-195, 205 |
| rseqc/FrameKmer.py    |       72 |       45 |     38% |21-36, 49-65, 72-77, 88, 91-98 |
| rseqc/SAM.py          |     1534 |      897 |     42% |44-45, 49-52, 77, 82-83, 217, 229, 240, 275, 285-286, 294-295, 301, 340-343, 351-353, 357, 368-379, 387, 391-392, 394, 398, 405, 410, 412-413, 415, 418-420, 490, 598-622, 634, 673, 713-714, 724, 730, 759-761, 776, 783, 864, 878, 880, 934-1034, 1047, 1061, 1063, 1075, 1080-1089, 1117-1219, 1233-1420, 1428-1607, 1622-1773, 1788-1987, 1994-2149, 2157-2233 |
| rseqc/\_\_init\_\_.py |        1 |        0 |    100% |           |
| rseqc/annoGene.py     |      239 |       82 |     66% |32, 34, 51, 53, 55, 92, 94, 96, 105-107, 124, 126, 128, 140-142, 157, 159, 161, 163, 197, 199, 201, 203, 216, 218, 229-246, 251-307 |
| rseqc/bam\_cigar.py   |      136 |       34 |     75% |59, 74, 76, 80-84, 99, 101, 103, 108, 123, 127-132, 147, 149, 152-155, 170, 174-179, 201-208 |
| rseqc/changePoint.py  |        8 |        0 |    100% |           |
| rseqc/cigar.py        |      107 |        4 |     96% |69, 93, 119, 149 |
| rseqc/cli\_common.py  |       29 |        0 |    100% |           |
| rseqc/fasta.py        |      172 |      112 |     35% |48-49, 67-70, 78-81, 87-102, 106-119, 125-129, 137-147, 152-213, 228-247, 251-263, 267-271, 275 |
| rseqc/fastq.py        |      134 |       51 |     62% |69, 99-116, 145, 152, 154-155, 269-321 |
| rseqc/getBamFiles.py  |       46 |       13 |     72% |46, 54, 57-65, 68-69, 74 |
| rseqc/heatmap.py      |       24 |        1 |     96% |        57 |
| rseqc/ireader.py      |       19 |        4 |     79% |     18-21 |
| rseqc/mystat.py       |      110 |        5 |     95% |59, 91-92, 105-106 |
| rseqc/orf.py          |      158 |      103 |     35% |25-26, 30-31, 37-38, 42-43, 85-196 |
| rseqc/quantile.py     |       19 |        3 |     84% |31, 50, 52 |
| rseqc/scbam.py        |      380 |      124 |     67% |187, 247-352, 357-362, 367, 372-373, 383-387, 390, 393-395, 410-448, 477, 505, 575, 585, 608, 624-627, 635-636, 642, 667, 673-674 |
| rseqc/twoList.py      |       33 |        6 |     82% |35-36, 53-54, 59-60 |
| rseqc/wiggle.py       |        1 |        0 |    100% |           |
| **TOTAL**             | **3470** | **1502** | **57%** |           |


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