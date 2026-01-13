---
title: 基于 CAFE 的基因家族分析
description: 本流程旨在从 CAFE5 输出结果中识别显著扩张或收缩的基因家族，通过树上标注和统计分析，直观展示显著扩张/收缩家族在系统发育树上的分布特征。
slug: cafe-genefamily1
date: 2025-11-21 00:00:00+0000
image: 1.jpg
categories:
    - 演化生物学
tags:
    - 比较基因组学
    - CAFE
    - 基因家族
    - 演化
---

## 1. 背景与目的

基因家族扩张与收缩是基因组进化中的重要现象。通过比较不同物种的基因家族大小，可以揭示基因在进化过程中经历的增减变化，进而理解物种如何适应不同生态环境以及其独特的进化机制。基因家族扩张往往与新功能获得、适应性特征的进化相关；而收缩可能反映某些功能的丧失或选择压力的变化。

---
## 2. CAFE5 简介

CAFE5（Computational Analysis of Gene Family Evolution）是一种常用的基因家族进化分析工具，它基于系统发育树和基因家族大小数据，利用最大似然法推测基因家族在各个分支上经历的扩张与收缩事件。该方法不仅可以量化基因家族变化，还可以判断这些变化是否显著，从而帮助研究者识别潜在的功能相关基因家族。

     Orthogroups.GeneCount.tsv          时间树（tree.txt）
             │                             │
             └──────────────┬──────────────┘
                            ▼
                    CAFE5 输入处理
                            │
                            ▼
                        CAFE5 分析
                            │
                            ▼
                    输出 Gamma_* 文件
                            │
                            ▼
                筛选显著扩张/收缩家族
                (Gamma_family_results.txt + Gamma_change.tab)
                            │
                            ▼
            ┌───────────────────────────┐
            │                           │
    sig_change_tsv.py           sig_change_map_to_tree.py
    输出 Gamma_change_sig.tsv       输出 cleaned_tree_sig_only.txt
            │                           │
            ▼                           ▼
    后续统计/功能分析           树上可视化显著家族扩张/收缩

---

## 3. 输入数据准备
在使用CAFE5时，至少需要准备两个输入文件：

- Orthogroups.GeneCount.tsv：基因家族的计数文件。
- tree.txt：系统发育树文件，包含物种分化时间。
### 3.1 从 OrthoFinder 的 Orthogroups.GeneCount.tsv 生成 CAFE5 输入文件

Orthogroups.GeneCount.tsv 文件用于记录每个基因家族（Orthogroup）在不同物种中的基因拷贝数，是 CAFE5 分析的核心输入之一。该文件通常来自 OrthoFinder 或 OrthoMCL 等软件，其中每一行对应一个基因家族，每一列对应一个物种的基因数量。为了让 CAFE5 正确读取，我们需要对原始文件进行一些格式检查和整理。

    cp  ../1_OrthoFinder/Results_Apr18/Orthogroups/Orthogroups.GeneCount.tsv 11.18.cafe
    sed 's/_//g' Orthogroups.GeneCount.tsv | awk 'BEGIN{OFS="\t"} {$NF=""; print}' | awk '{print "(null)\t"$$0}' | sed '1s/(null)/Desc/' > cafe.input.tsv

生成之后还需要剔除不同物种间拷贝数差异过大的基因家族，否则会报错，可以使用官方提供的脚本：https://github.com/hahnlab/cafe_tutorial/blob/main/python_scripts/cafetutorial_clade_and_size_filter.py

    python /home/salticidae/install/CAFE5-master/scripts/cafetutorial_clade_and_size_filter.py -i cafe.input.tsv -o gene_family_filter.txt -s
---

### 3.2 生成时间树
时间树的构建见博客xxx

---
## 4. 运行 cafe5

CAFE5的运行命令如下：

    cafe5 -i gene_family_filter.txt -t cafe.input.tree -o out_gamma_k1 -c 80 -k 1 -p    

一般运行 k=1~5，根据其输出的 Base/Gamma_results.txt 文件判断哪个拟合程度最好。
 - lnL（似然值）：越小越好。表征模型拟合度，-lnL 越小、模型越好。
 - Alpha：是否存在谱系异质性。Alpha 越大，表示更强的“家族进化速率的变异”；如果 Alpha ≈ 0，则说明 Gamma 模型没必要。
 - 失败家族（failure rates >20%）：太多失败说明模型可能不稳定，或某些家族数据异常。
---
## 5. 结果解读

    Gamma_asr.tre  # 每个基因家族的树文件
    Gamma_branch_probabilities.tab  # 每个分支计算的概率
    Gamma_category_likelihoods.txt
    Gamma_change.tab    # 每一个基因家族在每个节点的收缩与扩张数目
    Gamma_clade_results.txt # 每个节点基因家族的扩张/收缩数目
    Gamma_count.tab # 每一个基因家族在每个节点的数目
    Gamma_family_likelihoods.txt
    Gamma_family_results.txt    # 基因家族变化的p值和是否显著的结果
    Gamma_report.cafe
    Gamma_results.txt   # 模型，最终似然值，最终Lambda值等参数信息

### 5.1 每个节点显著收缩/扩张的基因家族数目可视化

将显著扩张/收缩的基因家族数目体现在树上，需要三个文件：

    Gamma_family_results.txt
    Gamma_clade_results.txt
    Gamma_asr.tre

从cafe5的输出文件 Gamma_asr.tre 中获得树文件，写入 id_tree.txt

去掉节点多余的内容：

    sed -E 's/(<[0-9]+>)[^:,;)]+/\1/g' id_tree.txt > cleaned_tree.txt

运行脚本 `sig0.05_change_map_to_tree.py` 
- 从 Gamma_family_results.txt 读入显著家族只保留 "y" 的基因家族；
- 从 Gamma_change.tab 选取显著家族对应的行；
- 每个节点分别统计：所有显著家族的扩张数，所有显著家族的收缩数；
- 最后将其 map 到树上：写入 cleaned_tree_sig0.05_only.txt
---
    #sig0.05_change_map_to_tree.py
    #!/usr/bin/env python3
    import re
    import pandas as pd

    # -----------------------------
    # 1. 读取显著家族列表
    # -----------------------------
    sig_fams = set()
    with open("Gamma_family_results.txt") as f:
        next(f)  # 跳过标题
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3 and parts[2].lower() == "y":
                sig_fams.add(parts[0])

    print(f"显著家族数: {len(sig_fams)}")

    # -----------------------------
    # 2. 读取 CAFE family × node 变化矩阵
    # -----------------------------
    df = pd.read_csv("Gamma_change.tab", sep="\t")

    node_cols = df.columns[1:]   # 第一列是 FamilyID

    # 仅显著家族
    df_sig = df[df["FamilyID"].isin(sig_fams)]
    print(f"显著家族矩阵形状: {df_sig.shape}")

    # -----------------------------
    # 3. 统计显著扩张/收缩的“家族数量”
    # -----------------------------
    node_change = {}

    for node in node_cols:
        changes = df_sig[node]

        inc = (changes > 0).sum()   # 扩张家族数量
        dec = (changes < 0).sum()   # 收缩家族数量

        node_change[node] = (int(inc), int(dec))

    print("每个节点显著扩张/收缩数量统计完毕。")

    # -----------------------------
    # 4. 读取树
    # -----------------------------
    with open("cleaned_tree.txt") as f:
        tree = f.read()

    # -----------------------------
    # 5. 替换树中的节点名称
    # -----------------------------
    for node, (inc, dec) in node_change.items():

        if inc == 0 and dec == 0:
            continue    # 两者都不显著则跳过

        new_label = node
        if inc > 0:
            new_label += f"+{inc}"
        if dec > 0:
            new_label += f"-{dec}"

        tree = re.sub(re.escape(node), new_label, tree)

    # -----------------------------
    # 6. 输出
    # -----------------------------
    with open("cleaned_tree_sig0.05_only.txt", "w") as f:
        f.write(tree)

    print("写入完成：cleaned_tree_sig0.05_only.txt")
---

### 5.2 过滤每个节点基因家族的扩张/收缩数目文件中不显著的基因家族
运行过滤脚本 `sig0.05_change_tsv.py`
- 以 Gamma_family_results.txt 里显著家族为准，从 Gamma_change.tab 中去掉不显著家族；
- 输出一个新的 Gamma_change_sig0.05.tsv。
---
    #sig0.05_change_tsv.py
    #!/usr/bin/env python3
    import pandas as pd

    # -----------------------------
    # 1. 读取显著家族列表
    # -----------------------------
    sig_fams = set()
    with open("Gamma_family_results.txt") as f:
        next(f)  # 跳过标题
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3 and parts[2].lower() == "y":
                sig_fams.add(parts[0])

    print(f"显著家族数: {len(sig_fams)}")

    # -----------------------------
    # 2. 读取 Gamma_change.tab
    # -----------------------------
    df = pd.read_csv("Gamma_change.tab", sep="\t")

    # -----------------------------
    # 3. 筛选显著家族
    # -----------------------------
    df_sig = df[df["FamilyID"].isin(sig_fams)]

    # -----------------------------
    # 4. 输出新的显著家族矩阵
    # -----------------------------
    output_file = "Gamma_change_sig0.05.tsv"
    df_sig.to_csv(output_file, sep="\t", index=False)

    print(f"完成，已生成 {output_file}，仅包含显著家族")
---
