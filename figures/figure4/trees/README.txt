AVG 三个关键基因的系统发育图（2026-11-05 修订）

本版只含三棵当前树，主图从左至右为 A Fhs、B FolD、C NrdG。
Fhs 使用用户提供的 Fhs_high5_NCBI_top5_54tips.treefile：5 条病毒候选及 49 条参考序列。
病毒候选为 vOTU_1650、2142、314、422、429，均保留。
使用这棵已有的 high5 分析树重新定根和绘图；未从 expanded52 树中直接删叶生成。
面板标题只写 Fhs；不加 high5/expanded52 字样；所有病毒使用相同样式，不加星号。
本版输入及图包只包含上述三棵树，不含 expanded52、历史 Fhs 或九面板版本。

快速查看
figures/AVG3_midpoint_combined.pdf / .svg / .png：一行三列圆形拓扑总览。
figures/current/：三棵树的独立圆形、矩形折叠、完整未折叠图和 midpoint Newick。
圆形图只表示拓扑，径向距离不代表进化距离；矩形图保留枝长及 substitutions/site 标尺。
总览中仅显示病毒、Curated 和折叠分支标签，其余参考末端保留彩色点。
完整参考标签、支持度见单树图及数据表。

运行
首次安装依赖：install.packages(c("ape", "digest"))
命令行：Rscript plot_AVG_trees.R
RStudio：打开脚本并点击 Source。脚本与 data/ 必须放在同一文件夹。
脚本用纯 R 读取原始树，重新 midpoint 定根、计算折叠并绘制全部图片。
无需 Python 或之前的 PNG。结果写入 figures/ 和 audit/。

数据
data/source/current_trees/：三棵原始 Newick 树。
data/source/tree_manifest.tsv：本版树清单与原始 SHA-256 校验值。
data/metadata/：末端 ID、来源、原始分类、统一门名、显示标签。
data/tree_order.tsv：面板标题、顺序及主图选择。
data/expected_collapse_summary.tsv、expected_collapsed_members.tsv：此前同一输入树的
折叠结果，用于核验 R 重新计算的结果，而不是替代计算。
这些是绘图所需的数据，不包含未随原始树提供的蛋白 FASTA 或比对矩阵。

定根与折叠
所有树均按原始枝长 midpoint 定根；支持值按无根分割对应到正确分支。
验证定根前后全部末端间距离、无根拓扑及支持值对应关系不变。
仅折叠至少两条序列、全部为 NCBI/PIBAC、同一已解析门、且同时满足
SH-aLRT >=80 和 UFBoot >=95 的最大后代分支。
所有病毒候选及 CUR 功能参考/边界对照都保留；空分类不折叠。
三角形 n 为所包含参考序列数量，不表示支持度。
名称别名统一，但 Bacillota、Bacillota_A/B/C 等类别继续独立。

配色与支持值
橙色：病毒候选；绿色：PIBAC；蓝色：NCBI；紫色：Curated。
Curated 表示输入树中人工选择的参考类别，不自动等同于实验验证。
深色内部节点：SH-aLRT >=80 且 UFBoot >=95；浅灰：未同时满足阈值。
单树数字顺序为 SH-aLRT/UFBoot。

可编辑参数
脚本顶部 COLORS、LABEL_PT、PANEL_LABEL_PT、PNG_DPI 控制配色、字号和分辨率。
OVERVIEW_ALL_LABELS <- TRUE 可恢复总览中所有参考文字。
FOCUS 为空，默认没有星号或单独强调的病毒。
总览画布为 27 x 9.5 英寸，PDF/SVG 可矢量缩放。

audit/ 提供定根验证、折叠统计、折叠成员和实际 R 运行版本信息。
