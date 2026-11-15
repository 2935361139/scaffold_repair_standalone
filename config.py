"""只改这里的输入、输出即可。读取已生成分子，不调用生成模型。"""

# sdf：一个生成 SDF；manifest：一批不同蛋白/参考；flowr_records：旧 FLOWR 记录。
INPUT_MODE = "sdf"

# 已按你上传的 DiffSBDD 五个样本填写。换一批数据时，四项须彼此对应。
INPUT_SDF = "/home/ljl/projects/DiffSBDD/outputs/official_example_001/generated.sdf"
PROTEIN_PDB = "/home/ljl/projects/DiffSBDD/example/5ndu.pdb"
REFERENCE_SDF = "examples/diffsbdd_example_001/reference_heavy.sdf"
ATOM_MAPPING_FILE = "examples/diffsbdd_example_001/atom_mapping.json"
# 上面的对应表仅针对已上传的五条记录，换模型/批次须核实原子对应，不能照搬数字。
INPUT_RECORD_INDICES = None  # None = 全部；[0, 1, 2] = 只检查 SDF 的前三条。

# 以下两项只有选择对应 INPUT_MODE 时才读取。
MANIFEST_PATH = "/请填写/样本清单.json"
FLOWR_GENERATED_DIR = "/home/ljl/projects/flowr_root/outputs/experiment_002/generated"

# 从蛋白副本排除“正在被替换的完整原参考配体”，防止分子与自身碰撞误报。
# 只在原子元素、坐标与整个参考配体唯一对应时自动排除；不按碰撞结果删水/辅因子。
RECEPTOR_OPTIONS = {
    "auto_exclude_reference": True,
    "match_tolerance_A": 0.01,
    "exclude_residues": [],  # 通常留空；需手动指定时见 README_CN.md。
}

# 相对路径以本 config.py 所在目录为基准。再次实验时改一个新目录名。
OUTPUT_DIR = "outputs/diffsbdd_example_001_v2"

# None 自动读取本次 OUTPUT_DIR/generated/failed。
# 也可填同一 v2 代码/环境产生的其他 generated/failed 路径；v1 结果要重新检查。
INTERVENTION_INPUT_DIR = None
FAILURE_IDS = None  # None = 全部失败样本；也可填 ["sample_001", "sample_002"]。
