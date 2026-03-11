# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                  |    Stmts |     Miss |   Cover |   Missing |
|---------------------- | -------: | -------: | ------: | --------: |
| rseqc/BED.py          |     2109 |     1843 |     13% |57-283, 297-332, 339-385, 389-466, 476, 478, 480, 552, 554, 556, 564-566, 586, 588, 602-689, 699, 701, 703, 714, 734-736, 750, 752, 754, 775-777, 782-808, 812-934, 938-978, 985-1047, 1053-1137, 1142-1181, 1198-1521, 1528-1613, 1617-1671, 1680-1685, 1693-1768, 1778-1961, 1968-2159, 2167-2271, 2279-2458, 2466-2562, 2569-2665, 2670-2763, 2774-2929, 2984 |
| rseqc/FrameKmer.py    |       71 |       44 |     38% |21-35, 46-64, 71-76, 87, 90-97 |
| rseqc/SAM.py          |     3948 |     3860 |      2% |57-62, 66-210, 216-452, 456-491, 495-512, 516-591, 596-643, 648-736, 741-763, 767-794, 802-888, 893-914, 919-944, 948-1002, 1010-1043, 1052-1057, 1065-1262, 1267-1371, 1378-1649, 1658-1964, 1971-2048, 2053-2162, 2168-2307, 2314-2454, 2461-2645, 2652-2745, 2752-2915, 2919-2977, 2981-3025, 3036-3047, 3051-3141, 3151-3262, 3272-3389, 3394-3442, 3447-3503, 3522-3903, 3907-4022, 4028-4107, 4111-4165, 4170-4258, 4264-4403, 4408-4547, 4552-4644, 4649-4846, 4854-5030, 5043-5074, 5081-5231, 5246-5449, 5454-5638, 5643-5647, 5654-5808, 5816-5899, 5903-5909 |
| rseqc/\_\_init\_\_.py |        1 |        0 |    100% |           |
| rseqc/annoGene.py     |      251 |      161 |     36% |41, 43, 59, 61, 63, 100, 102, 104, 115-117, 127-156, 163-199, 204-235, 241-260, 265-322 |
| rseqc/bam\_cigar.py   |      154 |       34 |     78% |73, 88, 90, 94-98, 113, 115, 117, 122, 137, 141-146, 161, 163, 166-169, 184, 188-193, 215-222 |
| rseqc/changePoint.py  |       31 |       12 |     61% |     44-57 |
| rseqc/cigar.py        |      118 |        4 |     97% |81, 105, 131, 161 |
| rseqc/dotProduct.py   |       35 |        4 |     89% |     18-21 |
| rseqc/fasta.py        |      171 |      109 |     36% |60-61, 79-82, 90-93, 99-114, 118-131, 137-140, 148-158, 163-216, 224-241, 245-250, 255-259, 263 |
| rseqc/fastq.py        |      108 |       53 |     51% |67, 97-121, 147, 154, 156-157, 214-260 |
| rseqc/getBamFiles.py  |       53 |       13 |     75% |57, 64, 67-75, 78-79, 84 |
| rseqc/heatmap.py      |       25 |       22 |     12% |     49-92 |
| rseqc/ireader.py      |       19 |        4 |     79% |     20-23 |
| rseqc/mystat.py       |      114 |       35 |     69% |63-80, 85-96, 109-110, 115-124, 139-140 |
| rseqc/orf.py          |      161 |      100 |     38% |39-40, 44-45, 49-50, 54-55, 94-201 |
| rseqc/quantile.py     |       25 |        7 |     72% |31, 50, 52, 61-63, 67 |
| rseqc/scbam.py        |      387 |      373 |      4% |31-37, 42-57, 64-68, 86-191, 220-421, 439-487, 538-662 |
| rseqc/twoList.py      |       39 |        6 |     85% |40-41, 58-59, 64-65 |
| rseqc/wiggle.py       |      112 |       79 |     29% |39-51, 58-59, 66-67, 74-75, 82-84, 91-93, 106-122, 128-129, 134-135, 141-142, 147-148, 154-155, 159-160, 166-168, 172-174, 178-180, 184-186, 192-209 |
| **TOTAL**             | **7932** | **6763** | **15%** |           |


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