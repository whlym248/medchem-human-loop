# 固定姿态复评分与 ADMET-AI 参考接口

`run_reference.py` 从此前实际执行过的研究复核脚本整理而来，调用已安装的 AutoDock Vina 对 **两份用户提供的姿态**执行 `vina` 和 `vinardo` 的 `--score_only`，随后用已有 ADMET-AI 模型预测两个输入。该接口不进行新一轮对接搜索、不生成新分子、不训练模型、不运行 MD，也不下载软件或权重。

本次公开适配仅检查 Python 语法和 CLI 帮助；没有重新运行复评分或推理。它是需要自行准备输入的参考工具，不属于根目录开箱可运行的合成 demo。公开包不包含此前研究的结构、结果、配置或模型权重。

## 调用

在已经安装对应软件的独立 Python 3.11 环境中，使用自己有权使用的输入：

```bash
python extras/fixed_pose_admet/run_reference.py --help
python extras/fixed_pose_admet/run_reference.py --root /path/to/prepared-evidence --config config.local.json
```

`--root` 必须指向外部准备目录。默认输出写入该目录下的 `runs/<timestamp>`，若已存在则拒绝覆盖。软件限制 Vina 一次使用 1 个 CPU worker，ADMET 在 CPU 上以 2 个线程执行；运行时间与环境有关。不要把科学输入或模型存入公开仓库。

## 输入目录约定

```text
prepared-evidence/
  config.local.json
  protocol.json
  input_manifest.json
  expected_models.json
  data/
    candidates.json
    <自己提供的 SDF、PDBQT 和受体结构>
```

`config.local.json`：

```json
{"vina_executable":"vina","vina_expected_sha256":null,"admet_models_dir":null}
```

生产复核建议设置可执行文件 SHA；`admet_models_dir=null` 使用所安装 ADMET-AI 的默认模型目录。指定模型目录时建议使用绝对路径；相对路径按进程的当前工作目录解释。这里不自动下载缺少的权重。

`protocol.json`：复制下面结构后，**用实际口袋定义替换演示框坐标和尺寸**，并记录自己的适用性审阅。`[0,0,0]` 不是建议的结合口袋。

```json
{
  "seed":260929,
  "cpu_per_scoring_job":1,
  "box":{"center":[0,0,0],"size":[30,30,30]},
  "scoring_modes":["vina","vinardo"],
  "expected_versions":{"admet-ai":"2.0.1","chemprop":"2.2.2","torch":"2.10.0","rdkit":"2026.3.6"}
}
```

这些版本来自此前执行环境，见 `requirements-reference.txt`；不是所有平台的安装锁文件。本脚本依赖 ADMET-AI 的两组各五个模型及其输出列，更换版本须核对接口和输出语义。

`data/candidates.json` 是 **恰好两个**对象组成的数组，每个对象需提供下列字段（示例仅说明字段；不附配套结构）：

```json
{
  "candidate_id":"compound_A",
  "target":"target_A",
  "record_group":"user_supplied",
  "smiles":"CCO",
  "canonical_isomeric_smiles":"CCO",
  "formal_charge":0,
  "pose_pdbqt":"data/compound_A.pdbqt",
  "structure_file":"data/compound_A.sdf",
  "receptor_pdbqt":"data/receptor.pdbqt",
  "receptor_structure_file":"data/receptor.pdb",
  "source_docking_seed":11,
  "source_pose_rank":1
}
```

`candidate_id` 必须唯一且仅含 ASCII 字母、数字、下划线或连字符。脚本检查 SMILES 的规范形式和形式电荷，但**不会自动证明 SMILES、SDF 与 PDBQT 描述同一化学实体或质子化状态**，也不会验证姿态科学合理性；这些属于输入准备审阅。

`input_manifest.json`：`{"files":[{"path":"data/candidates.json","sha256":"<64位SHA256>"}, ...]}`。必须覆盖 `data/candidates.json`、`protocol.json`、`expected_models.json` 和全部引用结构文件，路径相对 `--root`，不能越出该目录。哈希验证完整性，不验证科学内容。

`expected_models.json`：`{"files":[{"relative_path":"admet_classification/model_0.pt","sha256":"<可信来源权重的SHA256>"}, ...]}`。清单应覆盖实际模型目录内的所有 `.pt` 文件。使用官方/可信来源的模型，核对其许可；模型文件不会随本仓库发布。

## 产出与解释

- `fresh_scoring.json`、逐次评分日志：新执行的固定姿态评分与命令/日志哈希。
- `fresh_admet_full.csv`、`results.csv`：新执行的预训练模型输出与整理表。
- `runtime_metadata.json`、`status.json`：输入、模型、版本、设备、时间和执行结果。

这里的 `fresh` 只表示本次执行重新计算，不能理解为新药发现或新训练模型。Vina 的能量格式标签不是实验结合自由能；ADMET 分类值不是人体毒性发生率。工具不输出候选优胜排序或最终用药建议。

运行元数据可能包含本机绝对路径。分享结果前应另行检查权限、隐私和第三方数据条款。该工具的输出不会自动接入 `assemble_evidence.py`：目前应由使用者审阅后按其输入格式整理，完整协调器仍属后续工作。
