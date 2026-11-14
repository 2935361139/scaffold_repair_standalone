"""日常只改这个文件。此工具读取已生成分子，不负责调用模型生成。"""

# 第一次：读取以前完整的 FLOWR generated 目录，核对拆分前后的效果。
# 后续其他模型：改为 "sdf"；一批包含不同蛋白/参考配体时用 "manifest"。
#INPUT_MODE = "flowr_records"
FLOWR_GENERATED_DIR = "/home/ljl/projects/flowr_root/outputs/experiment_002/generated"

INPUT_MODE = "sdf"
# INPUT_MODE = "sdf" 时使用下面四个路径。
# 输入 SDF 包含多少个记录，就检查多少个，不会重新生成分子。                                    
INPUT_SDF = "/home/ljl/projects/DiffSBDD/outputs/official_example_001/generated.sdf"             #INPUT_SDF = "已生成分子的SDF路径"
PROTEIN_PDB = "/home/ljl/projects/DiffSBDD/example/5ndu.pdb"  # 与生成分子处于同一坐标系的完整蛋白 PDB；FLOWR 模式通常留空。    PROTEIN_PDB = "对应完整蛋白的PDB路径"
REFERENCE_SDF = "inputs/diffsbdd_official_example_001/reference_heavy.sdf"    # REFERENCE_SDF = "对应参考配体的SDF路径"
ATOM_MAPPING_FILE = "inputs/diffsbdd_official_example_001/atom_mapping.json"  # 格式见 examples/atom_mapping.json。   ATOM_MAPPING_FILE = "固定原子对应表的JSON路径"
INPUT_RECORD_INDICES = None  # None = 全部；[0, 1, 2] = 只检查 SDF 的前三条。

# 多种蛋白/参考配体/高精度坐标/可选轨迹：每个样本单独指定，格式见 examples/manifest.json。
MANIFEST_PATH = "/请填写/样本清单.json"

# 相对路径以本 config.py 所在目录为基准。再次实验时改一个新目录名。
OUTPUT_DIR = "outputs/diffsbdd_example_001"

# None 自动读取本次 OUTPUT_DIR/generated/failed。
# 也可填本独立工具以前产生的其他 generated/failed 路径。
INTERVENTION_INPUT_DIR = None
FAILURE_IDS = None  # None = 全部失败样本；也可填 ["sample_001", "sample_002"]。
