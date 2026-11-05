# 最新系统发育分析核实摘要

核实对象是 8 个基因家族的 9 个已完成树分析：Fhs 有 `high5` 与 `expanded52` 两棵完成树；另七个家族各一棵。所有输入比对 ID、IQ-TREE 报告序列数及最终树 tip ID 已逐一核对，完全一致。S3 收录 648 行（两棵 Fhs 树的序列分别计数），S4 收录 9 行。

## 查询和建树序列选择

| 家族/分析 | BLAST query 数 | 最终病毒 tips | 规则概要 |
|---|---:|---:|---|
| Fhs high5 | 5 | 5 | 只用 5 条高病毒来源置信度 Fhs 查询；同一检索招募的参考也用于另一棵 Fhs 树。 |
| Fhs expanded52 | 5 | 52 | 与 high5 共用检索；扩展入树集含 5 high、2 medium、45 low，新增 47 条没有作为 BLAST query。Fhs 病毒候选未做 CD-HIT。 |
| FolD | 86 | 18 | 125 条候选按 70% identity、双向覆盖各≥80%去冗余为 86 条 BLAST query；按已注释家族代表与未分类远缘代表选入树。 |
| NrdG | 386 | 22 | 612 条候选以相同 70%/80%条件去冗余为 386 条 BLAST query；入树集按家族代表与未分类远缘代表筛选。 |
| FolE | 49 | 15 | usable候选经 70%/80% CD-HIT 后查询；命名家族一条代表，未分类最多10条远缘代表。 |
| GlyA | 72 | 16 | 同上。 |
| IscS | 89 | 10 | 同上。 |
| IscU-like | 2 | 2 | 同上；本次候选数为2。 |
| MetK | 48 | 15 | 同上。 |

这里“高/中/低置信度”只描述病毒来源证据，不代表 ORF 或蛋白完整性。Fhs 当前流程并未把 52 条都作为数据库查询；旧文中“52条病毒序列全部用于检索”需替换。其余家族按序列可用性、病毒家族代表性和 ORF/结构域质量排序；不能概括为统一按高病毒置信度筛选。

## 实际检索参数

数据库为本地 NCBI **RefSeq protein**（非 nr），BLASTDB v5；记录的数据库元数据时间为 2026-08-27（445,458,864 sequences）。分类限制通过 `bacteria_archaea.taxids` 设置，为细菌和古菌。BLASTP 为 2.16.0+（FolD/NrdG 对应日志未单独留版本输出；其脚本/项目记录支持该版本，但论文表述宜注明版本记录范围）。所有搜索 `-max_target_seqs 20`，输出列为 `qseqid sacc staxids ssciname stitle pident length qlen slen qcovs evalue bitscore`，并用64线程。

Fhs 和 FolE/GlyA/IscS/IscU-like/MetK 搜索 E-value 为 0.001；FolD 和 NrdG 已保存检索的 E-value 为 1×10⁻⁵。命中后统一要求 bit score≥50、E-value≤0.001。没有额外设置 NCBI identity、query coverage 或 subject coverage 阈值。七个家族每个最终病毒树代表保留一个最佳 hit，排序为 bit score 降序、E-value 升序、query coverage 降序、accession 升序，并合并重复 accession；Fhs 则对5个 query 各取最多5个不同 accession（按相同质量排序）后合并。NCBI 参考不再 CD-HIT；选中 accession 以 blastdbcmd 抽取全长蛋白。BLAST task subtype及矩阵、gap、低复杂度等额外 flag 未在对应命令中记录，不能补写默认值。

## 猪肠道参考和其他参照

PiBAC关联的整合代表基因组集含 MAG 与分离株基因组，不是纯 MAG，也不是 UPGG；项目记录的集合规模为617个基因组。以 Prodigal 2.6.3 `-p meta` 预测蛋白，HMMER 3.4 搜索。各家族筛选证据如下：

| 家族 | 模型和猪肠道同源筛选要求 |
|---|---|
| Fhs | K01938；全序列 score≥913.80，最佳 domain HMM coverage≥0.70；partial=00、无内部终止。Fhs当前树沿用该批完成筛选并经CD-HIT的参考面板。 |
| FolD | PF00763.29 与 PF02882.26 均通过 gathering cutoff；每模型 HMM coverage≥0.80，N端模型先于C端模型，蛋白长230–350 aa；partial=00、无内部终止。 |
| NrdG | K04068 full-sequence score≥112.47，HMM coverage≥0.70，HMMER target length 150–300 aa；partial=00、无内部终止。身份为 putative NrdG-family。 |
| FolE | K01495 score≥66.33；PF01227.28通过GA，HMM及蛋白覆盖度均≥0.70；partial=00、无内部终止，并作类型区分。 |
| GlyA | K00600 score≥84.97，HMM coverage≥0.70；partial=00、无内部终止。 |
| IscS | K04487 score≥264.63，HMM coverage≥0.70；partial=00、无内部终止。 |
| IscU-like | PF01592.22通过GA，HMM coverage≥0.70、蛋白长≤250 aa；partial=00、无内部终止。 |
| MetK | K00789 score≥70.50，HMM coverage≥0.70；partial=00、无内部终止。 |

HMM coverage按命中 HMM 坐标区间/模型长度计算。其后按70% identity、双方序列覆盖各≥80%去冗余；最新五家族/Fhs运行使用 CD-HIT 4.8.1 `-n 5 -aS 0.80 -aL 0.80 -g 1 -d 0 -T 64 -M 16000`，FolD/NrdG已归档的病毒中间集用 `-n 4`；FolD/NrdG PiBAC CD-HIT运行的word size和完整命令未记录，不能套用 `-n 5`。Pig-gut最终面板最多40条，每个基因组一条、每属最多三条，优先覆盖病毒query前五命中；最新五家族使用BLASTP `-evalue 0.001 -max_target_seqs 1000`，命中需 bitscore≥50、qcovs≥50%，按未覆盖query数、支持query数、分数和ID贪心选取。FolD/NrdG另有单独面板日志：`-evalue 1e-5 -max_target_seqs 20`，病毒query覆盖≥50%。无top-five支持的候选可在名额剩余时补入。

参考中另有经过核对来源的 UniProt reviewed 功能序列及家族边界控制；FolD/NrdG 的 `CUR` 标签指 eggNOG v5来源的广义细胞同源序列，其中 accession/分类未核实者不能称作实验验证参照。PIBAC数据库来源、模型证据及序列分别见 `results/network5_rebuild_20261026/`、`results/folD_nrdG_pibac_refs_20261024/` 与 Fhs当前运行目录。

## 比对、修剪和树

MAFFT策略已统一为 `--auto`；Fhs与五个10月26日重建家族用64线程上限，FolD/NrdG用16线程。MAFFT 7.525在Fhs与五家族日志中可核实；FolD/NrdG run-specific版本未记录。trimAl均为 `-automated1`；Fhs版本文件记录v1.5.rev1，其余run-specific trimAl版本未记录。IQ-TREE均为3.1.2，`-m MFP`独立选模、BIC，`-B 1000 -alrt 1000 --bnni -seed 20261025`。Fhs及五家族用 `-T AUTO --threads-max 64`；FolD/NrdG用 `-T 16`。表S4给出各树模型及修剪后alignment长度。所有推断命令均未指定 `-o`；树按无根树报告。现有最终树流程未记录midpoint rooting，不能在Methods中称作midpoint-rooted。

## 旧方法描述需更新

- 删除 Fhs“5条高置信度只用于BLAST，随后无条件加入几十条低置信度序列”的旧流程版本描述；替换为当前5条query、top-five NCBI accessions，并分别报告 high5 和 expanded52 两棵完成树。
- 删除“NCBI使用nr数据库”及“所有家族使用同一个E-value”表述；实际为RefSeq protein，FolD/NrdG搜索E-value为1e-5，其余为0.001。
- 删除“所有家族NCBI命中均取单一best hit”；Fhs每条query取最多5条，其他家族才是一条best hit/selected viral query。
- 删除“所有序列以同一CD-HIT命令去冗余”；Fhs病毒query候选未聚类，FolD/NrdG病毒预筛word size为4，五家族/Fhs参考聚类word size为5；FolD/NrdG PiBAC聚类完整参数待确认。
- 删除“各家族 MAFFT/trimAl 版本和线程完全一致”“所有树按midpoint定根”以及旧的单一模型参数；目前策略均为MAFFT `--auto`，但线程上限不同；树未指定外群，模型按树独立选择，具体参数见表S4。
