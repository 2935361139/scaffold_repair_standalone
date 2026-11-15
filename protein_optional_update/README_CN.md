# 无蛋白输入更新：转换 → 检查 → 干预

这是现有 v2 和 universal_input_update 的增量更新，不是新的生成模型。不需要重新生成 DiffShape 数据，不加入分子对接。

## 1. 安装

将本包解压到 `/home/ljl/projects/scaffold_repair_standalone_v2`，得到 `protein_optional_update` 子目录。

在 WSL 终端运行：

```bash
cd /home/ljl/projects/scaffold_repair_standalone_v2
source /home/ljl/projects/flowr_root/.venv/bin/activate
python protein_optional_update/install_optional_protein.py --repo .
```

安装程序先核对当前代码，自动备份后更新。不会修改你的 config.py、import_config.py、输入数据或旧实验结果。如果显示 `Different local version`，没有写入任何修改；把完整提示发回来，不要强行覆盖。

## 2. 修改同一个 import_config.py

本次 DiffShape 仍使用原来的 raw 输入。只在 `IMPORT = {...}` 中添加或修改下面四项；其他项保持不变。请修改原有项，不要保留重复键。

```python
'protein_mode': 'absent',
'protein_pdb': '',
'output': 'inputs/diffshape_5cb2_002',
'experiment_output': 'outputs/diffshape_5cb2_002',
```

若 002 已经存在，后续命令中的相关目录名也统一换成 003 等新名字。不要删除旧结果，也不要覆盖旧检查后直接继续干预。

完整参考配置（路径确实与原先相同时可以使用）：

```python
IMPORT = {
    'format': 'diffshape_capture',
    'source': '/home/ljl/projects/Diff-shape/outputs/5cb2_raw_001/raw',
    'output': 'inputs/diffshape_5cb2_002',
    'model': 'Diff-Shape',
    'protein_mode': 'absent',
    'protein_pdb': '',
    'protein_frame_confirmed': False,
    'frame': 'fixed_centroid_translation',
    'units': 'angstrom',
    'maximum_registration_residual_A': 0.01,
    'experiment_output': 'outputs/diffshape_5cb2_002',
}
```

`absent` 是明确声明本次没有蛋白。默认仍为 `required`，因此误填/漏填蛋白路径不会自动变成无蛋白实验。

以后有匹配蛋白的数据：将 `protein_mode` 改为 `required`，填写 `protein_pdb`，在确认蛋白与参考配体坐标系一致后设置 `protein_frame_confirmed=True`。不要仅根据文件名猜测蛋白。

## 3. 按顺序运行

```bash
# 转换已有数据，不生成新分子、不修复、不传入作者的成功标签
python convert_inputs.py
```

看到 `STATUS: ready` 后，再运行：

```bash
# 检查、定位、按本次可评价的标准分为 success / failed
python run.py --config inputs/diffshape_5cb2_002/run_config.py --check-only

# 对失败样本执行固定预算的干预，并运行同输入的原通过样本作为对照
python run.py --config inputs/diffshape_5cb2_002/run_config.py --intervene-only
```

入口仍是 `import_config.py`。`run_config.py` 是转换程序自动生成的配置，不需要手动填写；上述命令不会读取根目录旧的 `config.py`。

## 4. 结果

```text
outputs/diffshape_5cb2_002/
  generated/success/     分子自身检查通过
  generated/failed/      分子自身检查未通过
  generated/review.zip   检查和定位报告
  intervention/success/  保留原通过样本，以及干预后通过的样本
  intervention/failed/   仍未通过的样本
  intervention/review.zip
```

把 `generated/review.zip` 和 `intervention/review.zip` 两个包发回来即可。

## 本次方法边界

- 有蛋白：继续原有五项条件的联合判定，原 PoseBusters dock 配置、阈值和有蛋白干预保持不变。
- 无蛋白：保留同一 dock 配置中所有分子自身检查及其原阈值，另保留化学有效性、连通性和固定片段检查。仅移除需要蛋白的加载、距离、碰撞和体积重叠检查。
- 被移除的检查逐项记录为 `not_evaluated`、`passed: null`，不进入已评价检查，也不算通过。
- 无蛋白分组标为 `molecule_only`。`success`、`full_success` 等兼容字段都只表示通过该行 `evaluation_scope` 下的全部检查，不能据此声称蛋白环境合格、对接成功或有活性。
- 同一批可含不同上下文的样本，但同一 input_id 不能混用不同蛋白条件。统计按 `evaluation_scope` 分开，不用无蛋白通过率直接比较原来的全环境通过率。
- 无蛋白时，B/AB/C 使用原来的分子内 UFF 和位移约束优化，不加入蛋白排斥项；重复的 UFF-only 环境消融分支不再执行。原有随机位移对照保留。图修复、构象提议和所有数值预算不变；没有补充针对某个样本的专用规则。
- C 类和 `_C` 分支若移动固定原子，仍单列为放宽坐标约束的结果，不与保持固定坐标的结果混同。
- 本次仍只进行有限结构层面的干预。没有生成轨迹就不能给出故障发生步数；终态修复成功不等于证明了唯一的神经网络原因。
- 骨架变化、连接位点细节等原有补充研究指标仍在 research_assessment.json 中，不会因为这次更新偷偷变成新的通过门槛。

## 可选的软件测试

```bash
python -m unittest discover -s tests -v
```

该测试只使用临时数据，不修改你的实验。包内 validation/ 保存本次软件测试和实测摘要。不同 RDKit/数值库版本可能导致临界样本结果不同，因此最终统计以你的本地运行包为准。
