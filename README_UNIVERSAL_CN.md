# 统一输入转换入口

本更新只添加通用的格式转换和必要的结构属性读写。检查阈值、干预规则、候选数量和成功条件不变。转换不是修复，不使用作者成功/失败标签。不能凭缺失信息推断正确映射、坐标、键或蛋白。

## 安装与本批操作

把本包内文件复制到 `/home/ljl/projects/scaffold_repair_standalone_v2/`，合并同名目录；不要删除原目录。新增 `import_config.py`，原 `config.py` 不改。

```bash
source /home/ljl/projects/flowr_root/.venv/bin/activate
cd /home/ljl/projects/scaffold_repair_standalone_v2
python install_inputs.py --repo .
python -m unittest discover -s tests -v
python convert_inputs.py --config import_config.py
```

逐条运行，前一步报错则停止。安装会备份并修改 `scaffold_repair/evaluation.py` 与 `scaffold_repair/adapters/sdf.py` 的输入读写部分。只接受此次提供的原 v2 文件版本；发现其他本地修改会停止，不覆盖。新结构属性只保存元素、形式电荷、氢计数、芳香标志及对应键图，不允许传入成功标签或原因。没有新属性的旧输入维持原行为。

配置已指向 `Diff-shape/outputs/5cb2_raw_001/raw`。不重新采样，无需模型权重/GPU，不使用 Diff-Shape 环境。

当前这批没有匹配蛋白，所以正常结果是10条已转换，状态 `needs_context_or_representation`，唯一待补充项为 protein PDB。转换生成 `manifest.pending.json`，不生成可误运行的 `run_config.py`。

请提交 `inputs/diffshape_5cb2_001/conversion_review.zip`，并补充与参考配体对应的蛋白 PDB。必须核对晶体配体与参考结构坐标关系，不能将其他蛋白或仅文件名相同的结构直接配上。当前 v2 强制要求蛋白；本更新不把缺少蛋白检测标成通过。

核对完成后，填写 `protein_pdb`，只有明确确认参考/PDB同坐标系才设 `protein_frame_confirmed=True`。若需额外坐标变换，应在适配配置中明确记录，不能为了减少碰撞任意对齐。更改 `output` 为新输入目录，重新转换已有原始图，不重新生成分子。状态为 `ready` 后：

```bash
python run.py --config inputs/新的输入目录/run_config.py --check-only
python run.py --config inputs/新的输入目录/run_config.py --intervene-only
```

输出目录由 `experiment_output` 指定。旧检查和干预结果不能在读写代码更换后继续恢复；新输入应重新检查。

## 支持范围

| format | source | 还需提供 |
|---|---|---|
| sdf | 多记录SDF文件 | 每条原子映射、参考SDF、蛋白及坐标约定 |
| manifest | v2格式清单JSON | 清单中所列文件；有轨迹/特殊target的清单请直接用原v2读取，转换器明确拒绝丢弃这些数据 |
| npz_graphs | NumPy图数组NPZ | graph_profile指定数组字段、元素和键类别表；每条固定映射及参考 |
| diffshape_capture | capture入口的raw目录 | 内置格式预设，读原始NPZ和采样掩码，不读取作者标签 |
| flowr_records | 原FLOWR完整generated目录 | 生成运行配置，交给原有FLOWR读取器，不重复实现其轨迹变换 |

模型名字只用于记录。新模型若输出这些结构，只改配置；私有二进制、未知类别编码、无可靠原子对应的输出仍需扩展解码或补充信息。不能保证任意文件零配置。

支持原子整数类别或原子序数、形式电荷、对称稠密整数类别键矩阵、N×3坐标。尚不解码 logits/one-hot、稀疏边表、padding、PDB配体键推断、SMILES的3D生成；这些应明确导出到上述规范，不能自动猜测。

## 通用SDF配置示例

```python
IMPORT = {
 'format': 'sdf', 'source': '/data/generated.sdf',
 'mapping_file': '/data/atom_mapping.json',
 'reference_sdf': '/data/reference.sdf',
 'protein_pdb': '/data/protein.pdb', 'protein_frame_confirmed': True,
 'frame': 'reference', 'units': 'angstrom',
 'output': 'inputs/batch_002', 'experiment_output': 'outputs/batch_002',
 'model': '任意模型名称'
}
```

mapping_file沿用v2的samples字段。必须覆盖全部SDF记录（包括不可读记录）；不能把第一条参考模板当成生成样本。混合多个参考时各条可指定reference_sdf和reference_index。固定映射下标均为输入原始编号，从0开始。

```json
{"samples":[{"id":"sample_001","sdf_index":0,
"fixed_generated_indices":[0,1],"fixed_reference_indices":[4,5]}]}
```

## 通用NPZ配置示例

```python
IMPORT = {
 'format': 'npz_graphs', 'source': '/data/graphs.npz',
 'mapping_file': '/data/samples.json',
 'reference_sdf': '/data/reference.sdf',
 'protein_pdb': '/data/protein.pdb', 'protein_frame_confirmed': True,
 'frame': 'reference', 'output': 'inputs/batch_003',
 'experiment_output': 'outputs/batch_003',
 'graph_profile': {
   'atom_encoding': 'index', 'atom_decoder': ['H','C','N','O','F'],
   'bond_decoder': {'1':'SINGLE','2':'DOUBLE','3':'TRIPLE','4':'AROMATIC'},
   'keys': {'atoms':'{id}__atom_types','charges':'{id}__charges',
            'bonds':'{id}__bond_types','positions':'{id}__positions'}
 }
}
```

每条samples包含id和两组固定原子下标，也可有graph_index。字段模板支持 `{id}` 和 `{index}`，因此可只改字段名/类别表而不修改Python。键0固定表示无键，荷电字段必须是实际形式电荷。原子用原子序数时改atom_encoding为atomic_number，不需要atom_decoder。

## 保真与边界

1. 生成图不sanitize、不重构键、不选最大碎片、不优化坐标。读回验证结构，转换错误保留原始输入和原因，整批不会只拿成功转换部分继续实验。
2. 显式氢不能直接丢弃：普通中性单键氢转换为重原子上的氢计数。带电/同位素/桥接等特殊氢不能由当前重原子接口无损表达，会明确停止该样本转换。保留源文件及哈希。氢计数保留后，原始过价错误不会因去氢而消失。
3. `frame=reference`声明分子已与参考一致；`unverified`只做待核对输出；`explicit_rigid_transform`接受rotation(3×3正旋转矩阵)和translation_A(3维)，行坐标公式x@R.T+t。
4. `fixed_centroid_translation`用已知固定重原子质心做整体平移注册；不旋转、不把单个原子强行放回参考位置。它定义了配体放置方式，并非准确恢复未知历史变换。只在固定片段残差低于指定阈值时接受；本批实测最大残差0.001904 Å，阈值0.01 Å。全部内部距离保持不变。蛋白与参考仍需另行核对。
5. 所有10个样本参与输入计数；conversion_report不是生成质量报告。真正的判定、问题定位和干预由v2另行执行。不读rdkit_constructed或作者cause字段作为判别结果。
6. 通用转换不能替代任务定义。当前5cb2官方配置名为linker-5cb2，不能仅凭运行入口就认定每个结果都是成功scaffold hopping；还需要报告所更换的结构和骨架变化。

## 验证

原v2 36项测试和新增14项输入测试通过（共50项）。新测试覆盖SDF/NPZ配置解码、标准接口读取、氢价态错误保留、芳香[nH]保留、未知模型名和作者标签不影响转换、遗漏映射拒绝、坐标单位、缺少蛋白阻断、输入重读、键端点顺序和手性/多碎片保留。安装重复运行不会重复修改；遇到未知版本会在修改前停止。

这批10个真实图已全部转换，且全部通过v2读取接口验证；内部重原子距离保持不变。3个原始化学无效图转换后仍无效，另外7个化学状态保持不变。这是格式保真验证，不是完整质量判定或修复率。读取接口测试中的临时PDB只用于测试文件访问，不用于科学评价；未为本批执行完整检查或干预。

本包基于本次可读取的原v2版本验证。它不是整个v2项目，不能放在一个空目录里单独运行。请保留原项目目录。若安装提示 `Different local version`，先停止并提交提示中的文件，不要强行覆盖版本检查。
