# Scaffold Repair：独立的检查与干预工具

**模型负责生成；本工具读取已经生成的分子，检查、定位、干预、再检查。**

它可以放在 FLOWR.ROOT 项目外。不会运行 FLOWR 的生成程序，也不需要模型权重。
当前保留原来的判别条件、A/B/AB/C 干预、15 个对照及固定的优化预算。
原项目不需要覆盖或删除。

## 1. 第一次：用之前的 FLOWR 记录核对

将整个文件夹解压到：

```text
/home/ljl/projects/scaffold_repair_standalone
```

打开本文件夹里的 `config.py`，先只确认：

```python
INPUT_MODE = "flowr_records"
FLOWR_GENERATED_DIR = "/home/ljl/projects/flowr_root/outputs/experiment_002/generated"
OUTPUT_DIR = "outputs/experiment_001"
```

第一个路径改成你实际的旧 `generated` 目录。读取的是这个指定批次，不会扫描其他历史实验或排除旧输入。
输入目录需要保留 `group_*/quality`、`diagnostics`、`runs` 和完整蛋白 PDB。
只有 `failed` 中单独的 SDF 不足以严格复现原来的高精度干预。

在 WSL 终端运行：

```bash
source /home/ljl/projects/flowr_root/.venv/bin/activate
cd /home/ljl/projects/scaffold_repair_standalone

python run.py --check-only
```

这一步：**读取旧分子，重新检查、定位并分类；不重新生成。**

上一条成功结束后，运行：

```bash
python run.py --intervene-only
```

这一步：**读取本次失败目录，执行原有干预与对照，再分类。**

也可以一次执行两步：

```bash
python run.py
```

原来的 `--generate-only` 作为 `--check-only` 的别名保留，但在独立工具中只代表导入与检查。

## 2. 结果在哪里

| 路径（位于 OUTPUT_DIR 下） | 内容 |
|---|---|
| `generated/success/` | 原始分子中通过全部判别条件的样本 |
| `generated/failed/` | 原始分子中未通过的样本 |
| `generated/samples/样本ID/diagnosis.json` | 失败检查、涉及原子/蛋白残基、候选原因、可用的轨迹观察 |
| `generated/samples/样本ID/evaluation.json` | 逐项判别结果 |
| `generated/summary.json` | 检查总数、成功数、失败数、支持的修复范围 |
| `intervention/success/` | 保留的原始成功样本，以及选出的成功修复结果 |
| `intervention/failed/` | 未修好、暂不支持或未选择干预的样本，状态分别记录 |
| `intervention/arms/样本ID/` | 每个干预/对照的 SDF、精确坐标、判别结果、优化记录 |
| `intervention/summary.json` | 修好几个、未修好几个、C 修好几个、多少分支没有完成评价 |
| `intervention/samples.csv`、`arms.csv` | 可用于进一步制作论文表格的数据 |

运行后可提交：

```text
outputs/experiment_001/generated/review.zip
outputs/experiment_001/intervention/review.zip
```

再次实验只改 `OUTPUT_DIR` 为新目录，例如 `outputs/experiment_002`。
可以重复使用同一输入；没有 20 个或 100 个的数量上限。程序只防止覆盖已有结果和混用被修改的证据。

## 3. 换成其他模型的 SDF

先在该模型项目中生成分子。然后改本工具 `config.py`：

```python
INPUT_MODE = "sdf"
INPUT_SDF = "/你的路径/generated.sdf"
PROTEIN_PDB = "/你的路径/protein.pdb"
REFERENCE_SDF = "/你的路径/reference.sdf"
ATOM_MAPPING_FILE = "/你的路径/atom_mapping.json"
INPUT_RECORD_INDICES = None
OUTPUT_DIR = "outputs/another_model_001"
```

仍然运行第 1 节的两条命令。

**四种输入各自做什么：**

1. 生成 SDF：要检查的最终分子，使用不含显式氢的重原子表示。
2. 完整蛋白 PDB：用于原来的口袋环境与碰撞检查，保留原实验使用的水/辅因子等。
3. 参考 SDF：生成前作为条件的参考配体，同样使用重原子表示。
4. 对应表：说明生成分子中的哪些原子，对应参考配体中哪些应保留的原子。

坐标单位必须为 Å，配体与蛋白必须已经在同一坐标系。本工具不会自动对齐或移动蛋白。
对应表不能仅凭原子数量猜测；模型导出时应保留原子 ID 或明确映射。

`examples/atom_mapping.json` 是格式例子；其中的数字只是示例，需要改成实际映射：

```json
{
  "schema": "scaffold_repair_mapping_v1",
  "samples": [
    {
      "id": "sample_001",
      "input_id": "reference_ligand_001",
      "sdf_index": 0,
      "reference_index": 0,
      "fixed_generated_indices": [0, 1, 2, 3],
      "fixed_reference_indices": [24, 25, 26, 27]
    }
  ]
}
```

含义：第一条生成分子的第 0/1/2/3 号原子，分别对应第一条参考分子的第 24/25/26/27 号原子。
**所有序号从 0 开始。每条选中的生成记录都需要一项映射。**
同一蛋白、参考配体和固定条件的不同生成样本使用相同 `input_id`；`id` 则每个样本唯一。

想检查 5 个，就提供 5 条 SDF 记录；或设置 `INPUT_RECORD_INDICES = [0, 1, 2, 3, 4]`。
这里控制的是检查多少个已生成样本，不是控制外部模型生成多少个。

## 4. 多种输入或高精度记录

一批包含不同蛋白、参考配体时，用 `examples/manifest.json` 格式，每个样本分别指定文件：

```python
INPUT_MODE = "manifest"
MANIFEST_PATH = "/你的路径/manifest.json"
```

清单中的相对路径以清单文件所在目录为准。

每项还可添加：

| 可选字段 | 要求 |
|---|---|
| `coordinates_npy` | `(N,3)` 的 float32/float64 原始高精度坐标；原子顺序与该 SDF 完全一致 |
| `target_coordinates_npy` | 同形状、同精度的原始条件坐标；A 使用其中的固定原子坐标 |
| `trajectory_npz`、`trajectory_metadata` | 可选轨迹及其说明，两者一起提供，见 `DATA_FORMAT.md` |

**只有 SDF 时，可以做最终结构检查与干预；不能保证与原始高精度轨迹的优化结果逐位相同。**
SDF 通常只有小数点后四位，边界样本及优化过程可能受舍入影响。
FLOWR 记录模式会恢复原有高精度坐标，并先核对旧基线标签；不一致时停止。

## 5. 方法与停止条件

成功仍为下列五项全部通过：化学有效、分子连通、固定子结构匹配、固定原子化学身份保留、原版 dock PoseBusters 全部通过。
是否更换 Murcko scaffold 单独报告，不擅自加入这五项的成功定义。
因此“通过检查”不自动等于“实现了新的 scaffold hopping”。

| 干预 | 保持的原方法 |
|---|---|
| A | 将固定原子放回条件坐标 |
| B | 固定原子不动，优化生成部分的坐标，同时考虑环境排斥 |
| AB | 先 A，再 B |
| C | 整个配体一起优化；允许固定原子的坐标改变，化学身份保持 |

再加 3 个 UFF 对照和 12 个位移匹配的随机方向对照，共 **19 个分支**。
每个样本最多做 6 次有界优化；每次最多 80 步、200 次目标函数调用，不成功就记录，不继续套新的搜索。
同一输入的原始成功样本会作为对照；其输出始终保留原样。
多个主干预成功会全部记录；导出分子按固定顺序 A→B→AB→C 选择，随机对照不作为修复产物。

当前保留的修复范围是**坐标/几何问题**。不能用移动坐标来声称修好了断键、错误键型、缺失原子或化学身份变化；这类样本明确列为暂不支持。
UFF 缺参数等运行错误记为“该分支未评价”，不会冒充一次分子失败或成功。

有轨迹时报告“观察到问题的步骤”和“持续出现的窗口”。没有轨迹时明确写未知。
网络内部层/模块的原因不能从最终 SDF 推断。本版本没有新增神经网络干预。
按你的操作性定义，能让整分子通过的干预可记作有效修复因素；它不证明唯一的上游根因。

## 6. 环境和复核

第一次可直接使用原来的虚拟环境，保证 RDKit/NumPy/SciPy 等版本一致；独立工具不导入 FLOWR 代码。
需要完全独立的环境时：

```bash
cd /home/ljl/projects/scaffold_repair_standalone
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

项目自带与原项目逐字节一致的 PoseBusters 副本，防止换成另一发行版后判别标准变化。
Torch 仅在 FLOWR 适配器读取已存在的 `inputs.pt` 时可选使用；有可验证 NPZ 记录时不需要 Torch。
不需要 GPU，也不需要 checkpoint。

实际回归验证记录见 `VALIDATION.json`。测试使用随包保留的真实历史数据，不能据此宣称另一个模型已取得相同修复率。

```bash
# 可选：重新执行真实 3 样本、57 分支对照；耗时明显长于基础测试。
python tests/regression.py --output outputs/regression_001
```

本包的 `tests/data/regression3.zip` 是你此前提交的检查/干预结果，用作回归基准，不是本次新生成的数据。
`SOURCE_PROVENANCE.json` 和 `SOURCE_CHECK.json` 记录原算法来源及校验。

## 7. 后续在哪里拓展

新模型一般只需要适配 `scaffold_repair/adapters/` 的数据格式；判别和几何干预不依赖其网络结构。
`evaluation.py` 管判别，`diagnosis.py` 管位置与轨迹观察，`geometry_intervention.py` 保留原优化公式，`engine.py` 管固定干预与对照，`pipeline.py` 管输入快照和分类。
将来新增断键/键型修复或网络干预，应作为独立方法新增并验证，不能把换模型路径等同于已经完成跨模型验证。
