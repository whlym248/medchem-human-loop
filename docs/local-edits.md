# 受限局部化学图编辑 / Bounded local edits

本模块实现人工指定一个位点后的有限图变换。它不使用训练模型，不从对接数据自动选择位点，不生成三维构象，也不运行 MD。结构能被 RDKit 解析、价态校验通过，不代表可合成、稳定、无毒或有效。

## 已实现的范围

`scripts/enumerate_local_edits.py` 支持以下固定规则，均在同一个人工指定的芳香 C–H 位点执行，每个候选只包含一次变换：

| 规则名 | 化学图操作 |
| --- | --- |
| `aromatic_ch_to_n` | 芳香 C–H 换为无氢的芳香 N |
| `aromatic_ch_add_cn` | 芳香 C–H 的 H 换为 –C≡N |
| `aromatic_ch_add_f` | 芳香 C–H 的 H 换为 F |
| `aromatic_ch_add_cl` | 芳香 C–H 的 H 换为 Cl |

这不是通用官能团编辑引擎。它拒绝带形式电荷的原子（包括净电荷为零的两性离子）、多片段盐、同位素、自由基、已有原子映射、显式独立氢原子，以及已指定或潜在的立体化学。限制是有意的；它不会悄悄移除盐、改变质子化状态、丢失手性或替用户选择未指定的异构体。

## 运行

需要 Python 3.10+ 和独立可选依赖 RDKit，见 `requirements-chem.txt`。本模块用 RDKit 2026.03.6 实际验证过；其它允许版本需在自己的环境重新执行测试。现有 MD 环境无需为运行本模块被修改。以下命令只做小分子化学图操作：

```bash
python scripts/enumerate_local_edits.py \
  --input examples/local_edit_request.json \
  --out outputs/local_edits_first

python -m unittest discover -s tests -p 'test_local_edits.py' -v
```

输出目录必须不存在；即使是已有空目录也不覆盖。没有 RDKit 时，化学图测试会明确显示 skip，命令会提示依赖缺失。Skip 不算该环境完成了化学验证。

示例为甲苯的简单教学操作，不含真实研究结果，也不是药物优化建议。它产生四个不同图结构。它不会修改输入 JSON。

## 输入与位点索引

```json
{
  "scope": "synthetic_example",
  "smiles": "Cc1ccccc1",
  "site_index": 3,
  "rules": ["aromatic_ch_to_n", "aromatic_ch_add_cn", "aromatic_ch_add_f", "aromatic_ch_add_cl"],
  "fixed_atom_indices": [0, 1],
  "max_candidates": 4
}
```

`site_index` 和 `fixed_atom_indices` 使用 **RDKit 解析输入 SMILES 时的零起始原子顺序**，不是输出规范 SMILES 的文字顺序。该示例第 0 号为甲基碳，第 1 号为与甲基相连的芳香碳，第 3 号是被选中的芳香 C–H。所有原输入原子保留编号，输出 atom map 为 `index + 1`；新原子依次使用后续 map。可通过报告的 `parent.atoms` 表核对原子。

支持的字段仅为上述六个：`smiles`、`site_index`、`rules` 必填；`fixed_atom_indices` 默认空列表，`max_candidates` 默认 4，可设 1–20。可选 `scope` 仅接受 `synthetic_example` 或 `existing_evidence_only`，会原样保留到报告顶层和 `request`，不能当作科学真实性认证。本仓库示例明确声明 `synthetic_example`；兼容旧输入省略 scope 时保持未声明，不擅自推断。请求最多 64 条规则、256 个原子、4096 字符 SMILES。候选顺序仅遵循规则输入顺序，不代表优先级。禁止任意 SMARTS、Python 代码或命令输入。

`fixed_atom_indices` 表示保护固定原子及其直接邻域。固定原子本身、其氢数、邻接原子编号与元素等标签、邻接键都必须不变。因此在固定原子旁把 C 换 N 也会被拒绝；在一个未固定的邻居上增加 F，若固定原子的直接邻接及邻居标签不变，则允许。固定集合不是空间距离约束；欲保护完整区域，应明确列出该区域的原子索引。

## 输出与检查

`local_edits.json` 包含：

- `schema_version: "local-edits/v1"`。
- `parent`：原始与规范 SMILES、各自 SHA256、带映射 SMILES、原输入索引表。
- `request`：标准化后的位点、规则、固定区域与上限。
- `candidates`：候选 ID、规则、规范及带映射 SMILES、结构 SHA256、父结构来源、编辑记录。
- `rejected`：未生成候选的规则及原因，包括不支持位点、保护区冲突、重复和上限。
- `provenance` 与 `warnings`：RDKit 版本、操作范围及限制。

`candidates.smi` 是方便接续工具读取的规范 SMILES / 候选 ID 两列文件。未生成候选时允许空文件，必须查看 JSON 拒绝原因，不能当作成功完成分子优化。

每个候选重新消毒校验（RDKit sanitization），核对芳香性、零形式电荷、无自由基、无新立体化学。检查原有原子映射不丢失、未选中原子标签与氢数保持、原有骨架键完全不变、新片段只连接所选位点、固定邻域保持。按去除 atom map 后的规范 SMILES 去重。输入 SHA 指原始 SMILES 字符串的 UTF-8 SHA256；规范 SHA 则用于化学图表示去重，它不是实验样品身份认证。

API 用法：

```python
from scripts.enumerate_local_edits import enumerate_edits

report = enumerate_edits(request)
```

错误请求（未知规则、无效原子索引、未支持的电荷或立体输入等）抛 `ValueError`，缺少 RDKit 抛 `RuntimeError`。合法请求中的单条变换失败留在 `rejected`，不自动换位点或扩大变换范围。结果不包含时间戳或随机行为，相同环境及输入下可重复得到一致 JSON。

## 如何接入 workflow

先由人根据结构与证据选择编辑位置和固定区域，检查父结构索引表，再运行有限枚举。候选只是待检验假说；后续仍需独立进行构象与质子化处理、对接、ADMET 预测、合成可行性审阅和按需进行的 MD。当前输出不自动衔接这些耗费资源的步骤，也不依据一个图有效性标记批准计算或推荐药物。

局部 diffusion、自动策略选择、亲和力提升的因果判断，以及从案例中学得的生成策略均不属于此模块的已实现能力。

## 从候选提案创建审阅包

`scripts/proposals_to_evidence.py` 将上述 JSON 转成 `assemble_evidence.py` 和 `workflow_round.py` 可读取的输入。它只用 Python 标准库核对来源字段、字符串 SHA256、候选 ID、规则和位点对应关系；不会独立验证 SMILES 化学图或重新执行生成器。即使调用方伪造了一份哈希自洽的报告，导出结果也不代表其科学真实性得到认证。

以本仓库的合成教学示例为例，先执行上面的枚举命令，然后运行：

```bash
python scripts/proposals_to_evidence.py \
  --proposals outputs/local_edits_first/local_edits.json \
  --target toy_target_A --target toy_target_B \
  --scope synthetic_example \
  --out outputs/proposal_evidence_first

python scripts/workflow_round.py prepare \
  --input outputs/proposal_evidence_first/input.json \
  --parent parent \
  --out outputs/proposal_review_first

python scripts/workflow_round.py status \
  --round outputs/proposal_review_first
```

所有输出目录都必须是新目录。导出目录包含 `input.json`、`parent.smi` 以及每个候选的 `edit_<sha前12位>.smi`。本示例为母体加四个候选、两个目标，共十条候选–目标记录。`workflow_round.py` 生成 `report.html`、`review_packet.json`、被冻结的输入和文件清单，等待人审阅。这里不运行任何 MD、对接、AI 推理，也不自动给出批准决定。

导出时所有 `properties` 都为空、对接分数为 `null`、MD 状态为 `not_available`；报告中即使夹带分数或完成标签也不会被带入。审阅包的差值保持缺失，不能把它们当成零分或“没有风险”。上游生成报告的完整字节 SHA256 保存在 `proposal_provenance.generation_file_sha256`。父结构统一使用 ID `parent`，候选 ID 使用被核验的安全文件名，避免输入内容变成路径。

`--scope` 必填，可选 `synthetic_example` 或 `existing_evidence_only`。示例必须使用前者。若提案已明确声明 scope，导出不能改变它，尤其不能把 synthetic 升格为真实证据。旧生成报告未携带 scope 时，CLI 值只是调用方未认证的声明，输出会显式标记 `scope_basis: caller_declaration_only` 和 `scope_authenticity_verified: false`。`existing_evidence_only` 只表示整理实际已有的文件，并不表示候选已有亲和力、ADMET 或 MD 结果；数值仍全部缺失。

这条桥接命令适合建立新的提案审阅包。它固定母体 ID 为 `parent`，不会自动把上一轮某个获准候选重新命名并继承其批准。要连接 `--previous-round`，必须满足协调器关于已批准父体 ID、结构文件 SHA 和 scope 均保持一致的检查；不能改名或重建文件后把它们假装成同一个已批准输入。
