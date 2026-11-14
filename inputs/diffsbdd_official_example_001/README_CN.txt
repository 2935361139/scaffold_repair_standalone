这两份输入仅匹配本次 official_example_001/generated.sdf 中的 5 个记录。
reference_heavy.sdf：原始参考配体去氢副本，49 个重原子的顺序、坐标和重原子间键保持不变。
atom_mapping.json：逐个核对后的 17 个固定原子对应表。不能直接复用到其他批次。
input_audit.json：输入校验、原文件哈希和原子匹配证据；不是分子质量判别阈值。
本包没有修改生成分子、蛋白、模型或检验代码。
