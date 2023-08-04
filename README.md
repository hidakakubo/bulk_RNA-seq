# bulk_RNA-seq

# Prepare the environment

### 仮想環境の構築
ターミナルで実行
```
conda create -n RNA-seq
```
```
conda activate RNA-seq
```

### packageのインストール
ターミナルで実行
```
conda install -c bioconda fastp==0.22.0
```
```
conda install -c bioconda hisat2==2.2.1
```
```
conda install -c bioconda samtools==1.3.1
```
```
conda install -c bioconda stringtie==2.1.1
```

### Rのpackage
RstudioのConsoleで実行
```
install.packages("tidyverse")
```
```
install.packages("BiocManager")
```
```
BiocManager::install("DESeq2")
```

# Mapping
やること
- raw dataのトリミング&クオリティチェック
- リファレンスゲノムへのmapping
- 発現量の定量化

### 準備が必要なディレクトリ
data  
解析に使うfastqファイル。dataディレクトリに格納。

HISAT2  
マッピングに必要なリファレンスゲノム（fasta）と定量化に必要なgtfファイルをダウンロードする。  
例：[ensembl](https://asia.ensembl.org/info/data/ftp/index.html)からファイルをダウンロード、gunzipコマンドで解凍。HISAT2ディレクトリに格納。

### scriptの変更点
最初に書いてあるchange working directoryのpathを書き換える。(mappingディレクトリまでcdできればよい)  
file_name.txtの作り方はサンプル名によって違う。  
HISAT2とStringTieで指定している、リファレンスゲノム(fasta)とアノテーションファイル(gtf)の名前を変更する  
スレッド数はお好みで変更するべし。

### ディレクトリ構造
解析開始時以下のようになっていれば良い
```
mapping
├── HISAT2
│   ├── Mus_musculus.gtf　(gtfファイル)
│   └── Mus_musculus.fa　(リファレンスゲノムのfastaファイル)
├── data
│   ├── SAMPLE_R1.fastq.gz
│   └── SAMPLE_R2.fastq.gz
└── script
    └── exec_mapping.sh
```

### 解析の実行
```
cd Desktop/bulk_RNA-seq/mapping/script
```
```
bash exec_mapping_pair.sh
```

# DEGanalysis
やること
- 各サンプルの発現量データをまとめた、gene_count_matrixの作成
- 発現変動解析

### 準備が必要なディレクトリ
data  
以下のようなサンプルのメタデータをmetadata.csvとして置いておく。(要手作り)
|sample_id|gene|
|:----|:----|
|TG_1|TG|
|TG_2|TG|
|TG_3|TG|
|WT_1|WT|
|WT_2|WT|
|WT_3|WT|

### scriptの変更点
- exec_prepDE_py.sh  
最初に書いてあるchange working directoryのpathを書き換える。(DEG_analysisディレクトリまでcdできればよい)  
- exec_DESeq2.R  
最初の方に書いてあるpathを書き換える。(DEG_analysisディレクトリまでのpath)  
39行目くらいのDESeqDataSetFromMatrix関数の引数を書き換える。  
```
# Step 2: construct a DESeqDataSet object ----------
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ #ここにmetadataにある実験の処理条件を入力(今回はgene))
```
53行目くらいのコードで基準となる処理条件を設定する。以下はWTと比べてTGがどれくらい変化しているかを知りたい場合。
```
dds$gene <- relevel(dds$gene, ref = "WT")
```

### ディレクトリ構造
解析開始時以下のようになっていれば良い
```
DEG_analysis
├── data
│   └── metadata.csv
├── sample_name_path.txt (mapping終了時に自動で作成される)
└── script
    ├── exec_DESeq2.R
    └── exec_prepDE_py.sh
```

### 解析の実行
gene_count_matrixの作成
```
cd Desktop/bulk_RNA-seq/DEG_analysis/script
```
```
bash exec_prepDE_py.sh
```
発現変動解析  
[こちらの動画](https://www.youtube.com/watch?v=OzNzO8qwwp0)を参考に1行ずつ実行するのをおすすめする。一気に実行する場合は以下のコード。
```
Rscript exec_DESeq2.R
```
