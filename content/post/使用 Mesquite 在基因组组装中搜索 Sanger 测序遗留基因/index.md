---
title: 使用 Mesquite 在基因组组装中搜索 Sanger 测序遗留基因
description: 本流程基于 Mesquite 4.01
slug: cafe-genefamily
date: 2025-11-04 00:00:00+0000
image: 1.jpg
categories:
    - 演化生物学
tags:
    - 系统发育
    - bycatch
---

感谢 Wayne P. Maddison 对流程实现的帮助

---

## 前言

系统发育研究常受限于高通量测序（NGS）数据难以覆盖国外物种，导致取样范围不足，影响结果的代表性与说服力。然而，NCBI 等公共数据库中往往已有这些物种的 Sanger 测序数据。针对这一问题，可采用下文策略：从已有的 NGS 数据中“钓取”（bycatch）目标基因的同源序列，将其与 NCBI 下载的 Sanger 序列进行整合，从而构建取样更全面、更具代表性的系统发育树。

---

## 1. 前提条件
 - 确保已安装 BLAST
   - 本流程依赖 BLAST 工具。
   - BLAST 可能已包含在 Mesquite 安装目录下的 apps 文件夹中。
   - 若未包含，则需在您的计算机上单独安装 BLAST。
 - 准备基因组组装 FASTA 文件
   - 将包含目标基因可能存在的 contig 组装结果（FASTA 格式）放入一个专用目录中。
   - 这些 FASTA 文件是后续 BLAST 搜索的目标数据库。
 - 准备目标序列
   - 下载高质量的目标序列

## 2. 主要操作步骤

### 2.1 启动 Mesquite

### 2.2 将 FASTA 文件转换为 BLAST 可用数据库
- 在 Mesquite 的 Log 窗口中，选择菜单：`Utilities > Make BLASTable Files from FASTA`
- 在弹出的对话框中，选择包含组装 FASTA 文件的目录。
- Mesquite 将为该目录中的每个 FASTA 文件生成对应的 BLAST 数据库。

### 2.3 加载目标序列（Target Sequence）
- 在 Mesquite 中打开一个包含**目标基因序列**的文件。
- 该目标序列应为高质量的已知序列，来自您的研究物种或其近缘种（关系越近，bycatch 结果越好）。
- 后续将用此序列对上述组装数据库进行 BLAST 比对。
  
### 2.4 在矩阵编辑器中选中目标序列行
- 在 `Character Matrix Editor` 中，点击选中代表目标序列的那一行。

### 2.5 启动本地 BLAST 搜索
- 选择菜单：`Matrix > Search > Top BLAST Matches`，然后选择 `BLAST Local Server`。
- 点击 OK，将弹出 BLAST 配置对话框。
  
### 2.6 配置 BLAST 搜索参数
该对话框为通用 BLAST 设置界面，部分选项在此用途下不适用，具体设置如下：
- a. “Additional BLAST options”（附加 BLAST 选项）中，无需手动指定最大 E 值（eValue）或字长（word size），这些将在下一步对话框中设置。
- b. 取消勾选 `BLAST databases in default location`（BLAST 数据库位于默认位置）。
- c. 在 `Path to folder`（数据库文件夹路径）中，选择第 2.2 中生成 BLAST 数据库的同一目录（即**您的组装 FASTA 所在目录**）。
- d. `Databases to search`（要搜索的数据库）默认为 *，表示搜索该目录下所有数据库。如只需搜索部分数据库，可在此处输入具体的数据库名称列表。
- e. 关于序列标题（header）格式的注意事项：
  - 某些组装软件（如 CLC）在 FASTA 序列标题中包含物种名，例如：>MySpecies Specimen 2 | NODE 348830...
  - 而其他软件（如 SPAdes）则不包含，例如：>NODE 348830...
  - 如果您的组装 FASTA 标题中不包含物种名，建议勾选 `Prepend database name to hit names`（在命中序列名前添加数据库名）。因为数据库名源自 FASTA 文件名，而文件名通常包含物种信息，此举有助于后续识别来源。


### 2.7 设置比对与序列导入选项
- 建议勾选 `reverse complement if needed and align imported sequences`（如需则反向互补，并对导入序列进行比对），以便直观查看结果。
- 是否允许内部插入空位（gaps）？
  - 如果目标序列为蛋白质编码区（如 COI 基因），不要勾选 `allow new internal gaps`；
  - 如果目标序列包含非编码区，则建议勾选此项。


### 2.8 执行 BLAST 搜索并导入结果
- 点击 OK 后，Mesquite 将调用 BLAST 对指定数据库进行搜索；
- 自动提取满足条件的 contig 序列，并将其作为新行添加到当前矩阵中。

### 2.9 修剪命中序列
- 根据目标参考序列，对命中的序列进行修剪

### 2.10 处理同一物种的多个命中
- 若某物种有多个 contig 被命中（例如一个 contig 覆盖基因前端，另一个覆盖后端），Mesquite 会为每个 contig 添加独立行。
- 此时，您可能希望将这些行合并为一条完整序列。Mesquite 将此操作称为 `Merge Taxa`（合并分类单元）。
- 合并方式有两种入口：
  - 在 `Character Matrix Editor` 中：`Matrix > Taxon Utilities > Merge Taxa`
  - 或在 `List of Taxa` 窗口中：`List > Taxon Utilities > Merge Taxa`
- 手动合并选定行
  - 先在矩阵中选中需要合并的多行；
  - 选择 `Taxon Utilities > Merge Selected Taxa`；
  - 系统将提供多种合并策略供选择。
- 按名称自动匹配合并
  - 选择 `Taxon Utilities > Merge Taxa by Name Matching`；
  - 系统会弹出对话框，询问名称中哪些部分需一致才能视为同一物种（例如前缀、文件名主体等）。
- 选择合并策略
  - 无论采用手动还是自动匹配方式，下一步都会出现合并选项对话框。
  - 推荐选择：`blend`：在无冲突位置合并碱基，冲突处标记；或 `refuse`（拒绝合并重叠区）：保留重叠区域供人工检查。

### 2.11 导出数据