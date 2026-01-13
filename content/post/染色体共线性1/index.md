---
title: 染色体共线性分析
description: 以两种跳蛛为例进行的染色体共线性分析初尝试
slug: synteny-analysis
date: 2026-01-12 15:52:00+0000
image: cover.jpg
categories:
    - 演化生物学
    - 生信
tags:
    - 共线性
    - 基因组 
weight: 1       # You can add weight to some posts to override the default sorting (date descending)
---

## 染色体共线性分析

以 *Siler cupreus* 和 *Portia taiwannica* (两种跳蛛) 为例

## 1. 物种间共线性分析

数据要求：  
1. 去冗余之后的蛋白质序列 *.faa
2. 基因组注释文件 *.gff

软件要求：MCScanX、mmseqs

### 1.1 blast
数据初步处理：只保留 gene ID，去掉 -RA 后的所有内容

    sed 's/-RA.*//' siler.faa > siler_filitered.fa
    sed 's/-RA.*//' portia.faa > portia_filitered.fa

#### 1.1.1 基于蛋白序列构建blast数据库

    mmseqs createdb siler_filitered.fa portia_filitered.fa all_mmseqs.db

#### 1.1.2 使用目标物种的 cds/pep/faa/fa 序列与此前构建的数据库进行比对，结果文件为*.blast

    mmseqs easy-search all.fa all_mmseqs.db all.blast tmp -s 7.5 --alignment-mode 3 --num-iterations 4 -e 1e-5 --max-accept 5 --threads 70

第三列序列相似性需要百分制，如果全都小于1的话，需扩大100倍

    awk -F'\t' 'BEGIN{OFS="\t"} {$3=sprintf("%.2f", $3*100); print}' all.blast > all_m8.blast

最终数据库格式如下：（gene ID 需要简化, 可参考后文替换命令）

    PT00001424	PT00001424	100	278	0	0	1	278	1	278	2.946E-177	548


### 1.2 *.bed 文件准备
要求：每个基因保留一个蛋白序列（通常为主转录本）

#### 1.2.1 将.gff文件转换成MCScanX所需格式（该命令需根据具体物种调整）

MCScanX 不需要原始*gff文件的所有内容，只需要提取其中的第 1, 9, 4, 5 列，其中第 9 列只提取 gene ID 即可

    awk -F "\t" '{a=substr($9,4,25)}$0~/gene/{print $1 "\t" a "\t" $4 "\t" $5}' siler_maker.gff > siler.gff
    
    awk -F "\t" '{a=substr($9,4,25)}$0~/gene/{print $1 "\t" a "\t" $4 "\t" $5}' portia_maker.gff > portia.gff

    cat siler.gff portia.gff > all.gff

串联得到的 all.gff 还需进一步调整格式，在这里我查看了文件前 5 行，如下：

    Scup_Un155	Siler_cupreus_00018918;Na	5062	30418
    Scup_Un155	Siler_cupreus_00018918-RA	5062	30418
    Scup_Un110	Siler_cupreus_00018926;Na	37738	48457
    Scup_Un110	Siler_cupreus_00018926-RA	37738	48457
    Scup_Chr6	Siler_cupreus_00012661;Na	359050	412595
    Scup_Chr6	Siler_cupreus_00012661-RA	359050	412595

接下来根据 all.gff 内容，需要作如下调整：  
1. 删除所有带有 “;” 的行
2. 删除所有未成功挂载的染色体（带有 "Un" 的行）
3. 简化第一列的命名
4. 简化第二列的命名
5. 去重复

#删除所有带有 “;” 的行

    grep -v ";" all.bed > all1.bed

#删除所有未成功挂载的染色体（带有 "Un" 的行）

    grep -v "Un" all1.bed > all2.bed

#简化第一列的命名

    sed "s/Scup_Chr/sc/g" all2.bed > all3.bed
    sed 's/Ptai_Chr/Pt/g' all3.bed > all4.bed

#简化第二列的命名

    #第二列中 -** 的去除
    awk 'BEGIN{OFS="\t"} {
    for (i=1; i<=NF; i++) {
        sub(/-.*/, "", $i);  
        }
    print
    }' all4.bed > all5.bed

    #第二列中 gene ID 前缀简化  
    sed 's/Siler_cupreus_/SC/g' all5.bed > all6.bed
    sed 's/Portia_taiwanica_/PT/g' all6.bed > all7.bed

#去重复

    sort all7.bed | uniq > all.gff

此时得到的all.gff文件如下（当然，在此之前需要删除最初的 all.gff )

    pt1	PT00000001	923338	1234673
    pt1	PT00000002	1283053	1329480
    pt1	PT00000003	1423863	1459005
    pt1	PT00000004	1510151	1568055
    pt1	PT00000005	1568392	1582644


### 1.3 运行 MCScanX
将 all.blast 和 all.gff 置于一个单独文件夹（这里使用 data/），运行 MCScanX

    MCScanX data/all -b 2 -s 5 -e 1e-5

**注意：**  
1. all.blast 和 all.gff 两个文件命名除了后缀必须完全相同
2. 两个文件中的gene ID 必须完全相同
3. MCScanX 运行命令在文件夹后必须跟上统一的文件名

### 1.4 可视化

使用上一步得到的 all.gff 和 all.collinearity 文件

在 https://www.chiplot.online/ 中绘图