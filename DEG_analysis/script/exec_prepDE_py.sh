#!/bin/bash

SECONDS=0


# change working directory
cd /Users/hidaka/Desktop/bulk_kinoshita/DEG_analysis


# prepDE.py
# prepDE.pyをダウンロード
mkdir -p prepDE_py
curl -o ./prepDE_py/prepDE.py https://ccb.jhu.edu/software/stringtie/dl/prepDE.py

# prepDE.pyをpython3のコードにする
2to3 -w prepDE_py/prepDE.py

# prepDE.pyの実行
python prepDE_py/prepDE.py -i sample_name_path.txt

# fileの移動
mv gene_count_matrix.csv transcript_count_matrix.csv data/
