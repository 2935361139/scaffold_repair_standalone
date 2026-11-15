# 日常只改本文件。相对路径以本文件所在目录为基准。
IMPORT = {
    'format': 'diffshape_capture',  # sdf / manifest / npz_graphs / diffshape_capture / flowr_records
    'source': '/home/ljl/projects/Diff-shape/outputs/5cb2_raw_001/raw',
    'output': 'inputs/diffshape_5cb2_002',  # 转换文件保存位置；重复运行请换新目录名
    'model': 'Diff-Shape',  # 仅记录来源，不改变检查或干预规则
    'protein_pdb': '',  # 当前尚缺与参考配体匹配的蛋白 PDB；不要随便配一个
    'protein_mode': 'absent',
    'protein_frame_confirmed': False,  # 核对蛋白与参考配体的坐标系后，才能改 True
    # 用已知固定原子的质心做整体平移，保持分子内部距离和原有畸变。
    # 这是明确声明的放置方式，不代表准确恢复了未知的历史坐标变换。
    'frame': 'fixed_centroid_translation',
    'units': 'angstrom',
    'maximum_registration_residual_A': 0.01,
    'experiment_output': 'outputs/diffshape_5cb2_002',  # 后续检查/干预结果保存位置
}
