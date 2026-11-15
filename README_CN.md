# Scaffold Repair v2

读取已经生成的分子，检查结构问题，执行有限的修复与对照，再按同一标准复查。
不调用生成模型，不需要模型权重或 GPU。

## 1. 安装和运行

把压缩包解压到 `/home/ljl/projects/`，得到完整的新文件夹 `scaffold_repair_standalone_v2`。
旧项目、旧结果都保留，不要把新旧文件零散混合。

```bash
source /home/ljl/projects/flowr_root/.venv/bin/activate
cd /home/ljl/projects/scaffold_repair_standalone_v2
python -m unittest discover -s tests -v
```

继续使用原来已运行成功的 Python 环境；不需要重新安装依赖，也不用切换到 DiffSBDD 的 Conda 环境。
`requirements.txt` 只用于以后新建独立环境。

打开 `config.py`，默认已填好你这批 DiffSBDD 五个样本的路径：

```python
INPUT_MODE = "sdf"
INPUT_SDF = "/home/ljl/projects/DiffSBDD/outputs/official_example_001/generated.sdf"
PROTEIN_PDB = "/home/ljl/projects/DiffSBDD/example/5ndu.pdb"
REFERENCE_SDF = "examples/diffsbdd_example_001/reference_heavy.sdf"
ATOM_MAPPING_FILE = "examples/diffsbdd_example_001/atom_mapping.json"
INPUT_RECORD_INDICES = None
OUTPUT_DIR = "outputs/diffsbdd_example_001_v2"
```

参考配体和对应表已附在包中。**这份对应表只适用于之前上传的这五条记录。**
换批次/模型时，生成 SDF、蛋白、参考配体、固定原子对应表必须相互匹配，不能照搬原子序号。

然后运行：

```bash
# 读取原始分子，检查、定位和分类；不重新生成
python run.py --check-only

# 上一步成功后，针对失败样本执行有限干预和对照
python run.py --intervene-only
```

也可以用 `python run.py` 连续完成两步。重复实验只改 `OUTPUT_DIR` 为新目录。
没有 5、20 或 100 个的数量上限，也不扫描旧输入历史；只防止覆盖已有实验。
**v1 的 failed 目录不能直接进入 v2 干预，要先用 v2 重新检查原始输入。**

## 2. 输入要求

- 生成分子和参考配体用重原子三维 SDF；固定原子映射从 0 编号。
- 蛋白和配体已在同一坐标系，单位 Å；工具不猜测对齐或原子顺序。
- `INPUT_RECORD_INDICES = None` 检查全部记录；`[0, 1]` 只检查前两条。
- 不同蛋白/参考配体混合的一批数据，用 `INPUT_MODE = "manifest"`，格式见 `examples/manifest.json`。
- 旧 FLOWR 记录仍可选 `INPUT_MODE = "flowr_records"`，填写 `FLOWR_GENERATED_DIR`；`PROTEIN_PDB = ""` 表示从记录读取。此适配器不是跨模型必需依赖。
- 所有路径可以修改；相对路径以 config.py 所在文件夹为基准。

## 3. 这次改了什么

| 文件 | 作用 |
|---|---|
| `scaffold_repair/receptor.py` | 准备受体副本，排除正在被替换的原参考配体，不改原 PDB |
| `posebusters/tools/protein.py` | 修正蛋白、水、有机和无机辅因子类别交叉的问题 |
| `scaffold_repair/topology.py`、`chemistry.py` | 固定片段键/电荷恢复、局部化学价异常的键修复、有限跨片段连接候选 |
| `scaffold_repair/conformers.py` | 固定片段约束下的构象重建；不调用生成模型 |
| `scaffold_repair/coverage.py` | 文献问题分类、内部几何定位、附加任务/性质检查与批次指标 |
| `scaffold_repair/diagnosis.py` | 定位断连片段、固定键差异、候选连接，而不只报告碰撞 |
| `scaffold_repair/engine.py` | 分别执行图结构、几何及组合干预，同时保留对照 |
| `scaffold_repair/pipeline.py` | 区分未支持、无合适候选、干预后仍失败、运行错误和成功 |

规则依据原子、键、固定片段映射、距离和化学价，不依据模型名、样本编号或这五个分子的特殊位置。
这是扩展了结构问题覆盖范围的通用接口，不是所有模型、所有失败都能修好的保证。

## 4. 判定标准和有限搜索

成功仍需五项全部通过：RDKit 化学有效、整分子连通、固定子结构匹配、映射固定原子的元素/电荷/内部键保留、随包 dock PoseBusters 启用的检查全部通过。
数字阈值没有放宽。旧环境原子选择错误已修正，旧、新标签不能不加区分地混合统计。
本次干预前后使用同一版本评估代码和同一准备好的受体。
Murcko scaffold 是否改变另列；通过质量检查不自动等于实现 scaffold hopping，更不代表生物活性。

| 代号 | 内容 |
|---|---|
| A / B / AB / C | 原坐标修复：固定原子归位、生成部分优化、两者结合、全配体优化 |
| F | 根据明确映射恢复已知固定片段键/电荷，必要时有限删除冲突的固定—生成区连接 |
| L | 为相邻断连片段提出能通过化学价检查的单键连接 |
| FL | 固定片段修复与连接修复的组合 |
| E / EL | 对 RDKit 指出的异常区域做有限局部键修改；必要时组合连接修复 |
| G / F1G1 等 | 保留图结构、在固定原子约束下重建构象；可基于原图或已提出的图修复候选 |
| F1_B 等 | 在某个图修复候选上优化几何，固定原子坐标保持不变 |
| F1_C 等 | 图修复后允许全配体移动；放宽固定坐标的结果单列 |
| F1_UFF_control 等 | 同一图修复后的 UFF 对照，不加环境排斥，不作为主修复产物选择 |

每个图候选先单独评价“只改图”，再做几何组合；原图的几何分支用于检验“只移动坐标是否足够”。
组合才成功时，记录组合有效，不能把两个因素分别说成独立原因。
如果 UFF 对照也成功，不能声称环境排斥是修复必需因素。

固定预算：最多 4 个图候选、128 次图合法性尝试；每个连接方案最多增加 3 个键。
再加最多 2 个构象候选，使用固定种子；每个种子最多 30 次嵌入尝试、20 秒。
原有几何及对照共 19 个分支；每个图/构象候选增加 4 个分支，总计每样本最多 43 个评价分支。
最多 24 次坐标优化，每次最多 80 步、200 次目标函数调用，达到上限即结束，不自动追加轮次。
候选筛选不使用最终通过标签调参。记录所有通过的主干预；导出时优先不放宽固定坐标，再优先较少改键，再按预定顺序。
原本成功的分子保留原样；对照结果不冒充主干预修复。

## 5. 受体准备不是分子修复

默认只在某个 HETATM 残基的全部重原子，与参考配体元素和坐标唯一匹配时，排除该原参考配体。
不通过删除碰撞水或其他辅因子让分子通过。保留原 PDB，并记录移除残基、校验值；保留原子坐标不变。
参考配体本来不在 PDB 中时，不排除任何残基。其他未匹配的有机辅因子仍保留。

只有明确知道条件配体身份、但自动匹配不适用时，才手动指定。例如本示例可用：

```python
RECEPTOR_OPTIONS = {
    "auto_exclude_reference": False,
    "match_tolerance_A": 0.01,
    "exclude_residues": [{"chain": "A", "number": 201, "name": "8V2"}],
}
```

这组残基编号仅属于当前示例，不能照搬到其他蛋白。不要为通过检查排除实际需要保留的辅因子。
若配体与蛋白坐标不一致，应纠正输入，不能强行删除原子或扩大容差。

## 6. 输出在哪里

以下路径均位于 OUTPUT_DIR 下：

| 路径 | 内容 |
|---|---|
| `generated/success/`、`generated/failed/` | 原始分子检查后的分类 |
| `generated/samples/样本ID/diagnosis.json` | 结构问题、原子位置、候选解释；无轨迹时不猜生成时间步 |
| `generated/samples/样本ID/receptor_preparation.json` | 受体准备记录 |
| `generated/samples/样本ID/research_assessment.json` | 内部几何位置、局部立体检查、附着点/骨架改变和性质指标；扩展结果独立于原五项结构门槛 |
| `coverage_catalog.json`、`batch_metrics.json` | 问题覆盖范围、文献、批次重复/多样性和未评估项目；两个阶段均输出 |
| `intervention/arms/样本ID/planning.json` | 预定候选、具体改键、预算和分支 |
| `intervention/arms/样本ID/` | 全部分支分子、精确坐标、检查及优化记录 |
| `intervention/success/`、`intervention/failed/` | 最终通过的分子和未修复分子 |
| `intervention/summary.json` | 原失败数、修复数、仍失败数、实际执行数及错误数 |
| `intervention/results.json` | 全部干预和对照的详细结果 |

结束后提交这两个文件即可：

```text
outputs/diffsbdd_example_001_v2/generated/review.zip
outputs/diffsbdd_example_001_v2/intervention/review.zip
```

## 7. 已验证范围和边界

测试和实测结果见 `VALIDATION.json`、`validation/`。
文献依据、已覆盖和不能自动处理的类别，见 `LITERATURE_COVERAGE_CN.md`。没有假定存在一个能自动修复所有药物设计问题的算法。
每次运行记录输入、软件和代码版本。`generated/integrity.json` 只检查本次证据未中途改变，不排除旧输入。

不增加/删除原子，不替换元素，不任意重建生成区全部键级，不处理金属配位连接，也不补长距离缺失连接链。
这些情况保留为未解决；“无候选”与“已干预但失败”分别记录。分支运行报错单列为未评价。
最终 SDF 只能支持终态结构诊断与干预，不能定位网络层或生成时间步；可选轨迹格式见 `DATA_FORMAT.md`。
修复成功支持相应结构干预有效，不证明唯一的上游根因。

可选：运行原 FLOWR 三样本的 57 分支回归对照（比单元测试慢）：

```bash
python tests/regression.py --output outputs/flowr_regression_v2
```

该对照仅针对随包历史样本；其通过不意味着所有历史数据的标签在修正环境分类后都不变。
