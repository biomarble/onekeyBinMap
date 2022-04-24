

## 分析环境准备

### 建立分析环境

```bash
mamba create -n linkageMap ncurses=6.3  python=2.7 numpy pillow openjdk mstmap r-ggplot2 r-qtl r-rlecuyer r-snow perl-getopt-long csvtk ngsep
```

### 激活环境
```bash
conda activate linkageMap
```

### 切换目录

```bash
cd /mnt/d/workshop/6.genetic_map/inbred
```

##  Bin的构建

Bin的构建使用[SNPBinner](https://github.com/solgenomics/SNPbinner)软件，软件适用于各种F2/RIL群体。

### 安装snpbinner
```bash
git  clone https://github.com/solgenomics/snpbinner
python SNPbinner/setup.py

```

### 一键式snpbinner

```bash
perl scripts/snpBinnerCommands.pl  -i data/demo.inbred.vcf  -F  BY804  -M  B73  -od BinResult 
```

### 结果文件

|文件夹||
|-|-|
|0.input|转换后的snpbinner输入文件|
|1.bin |SNPbinner的crosspoint分析与bin划分结果|
|2.plot|每个样本在每个连锁群的bin和crosspoint的可视化结果|
|3.final|最终的bin基因型文件(binGeno)和bin信息文件(binInfo)|


## 遗传图谱构建

### MSTmap

MSTmap是一款经典的遗传图快速排序软件。


#### 一键式分析

```bash
perl scripts/runMSTmap.pl  -g BinResult/3.final/all.binGeno.csv  -info BinResult/3.final/all.binInfo.csv -miss 0.5 -pcut 0.001  -pop RIL6 -od MapResult 
```

#### 结果文件

|文件|含义|
|-|-|
|MST.map.txt|最终的map文件|
|filter.all.pvalues.csv|每个标记的基因型计数，偏分离显著性，缺失率等|
|filter.valid|最终构建遗传图谱使用的标记名称列表|
|LG.\*.MSTmapIn.txt|每个连锁群的MSTmap输入文件|
|LG.\*.MSTmapOut.txt|每个连锁群的MSTmap输出结果文件|