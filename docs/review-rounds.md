# 从证据到人工决策的一轮流程

`scripts/workflow_round.py` 是第一版轻量协调器：冻结已经提供的证据，可接入本地案例检索和透明路线规则，生成静态 HTML 审阅页，再记录人的决定。它不调度分子模拟或远程机器，不自动运行对接、ADMET、编辑器或扩散模型。

## 用合成数据走一轮

从仓库根目录执行，Python 3.10+，仅需标准库：

```bash
python scripts/workflow_round.py prepare --input examples/synthetic_evidence.json --parent toy_parent --out outputs/round1
```

打开 `outputs/round1/report.html` 查看各候选在各靶点的证据、缺失项和相对母体的数值差。所有示例数值均为人为构造；它们没有经过真实 docking 或 MD。

记录一个仅用于练习的决定：

```bash
python scripts/workflow_round.py decide --round outputs/round1 --candidate toy_variant --action approve_next_round --reviewer demo-reviewer --reason "Synthetic demonstration only; no real candidate approval" --acknowledge-limitations
python scripts/workflow_round.py status --round outputs/round1
```

可选择 `approve_next_round`、`reject`、`request_evidence` 或 `defer`。真实使用时应由实际审阅者提供决定、姓名/标识和具体理由。程序不会认证输入姓名的人是谁。`approve_next_round` 仅指研究方向选择，不代表技术 pilot 审批、药理安全验证或付费资源授权。

修改决定会生成下一条记录，不覆盖旧事件。例如新增 `defer` 会成为该候选的最新决定，同时保留先前意见。

## 案例与路线建议

```bash
python scripts/workflow_round.py prepare --input examples/synthetic_evidence.json --parent toy_parent --out outputs/round_with_context --cases examples/cases/synthetic_cases.jsonl --query aromatic --route-request examples/route_request.json
```

案例库按文本和标签匹配，只有满足记录审核条件的条目会进入结果。合成案例只允许用于合成输入，且永远不形成真实案例支持。真实证据包默认排除合成案例。

路线规则使用搜索到的真实文献案例支持状态；它不会采信调用者单独宣称“已有支持”的布尔值。文献记录中的审核状态仍是提供者声明，必须检查文献、实验条件和可迁移性。没有支持时，建议可能是先补证据；diffusion 没有实现后端时只能等待后续模块，不能假装调用成功。

## 下一轮与身份检查

```bash
python scripts/workflow_round.py prepare --input examples/synthetic_evidence.json --parent toy_variant --previous-round outputs/round1 --out outputs/round2
```

这只是用同一组合成数据演示轮次关联，没有产生新的科研结果。真实下一轮需要新候选和重新审阅的证据输入。

指定 `--previous-round` 后，新母体必须在上一轮的**最新**决定中获准继续，且 SMILES 字符串和结构文件 SHA256 均一致。程序拒绝把合成轮次改标签后接成真实证据轮次。任何规范化或质子化改变应先作为新的明确结构身份处理；精确文本/哈希匹配本身不等于化学同一性验证。

## 冻结与追溯范围

输出包括：

- `evidence/`：本轮输入与相关结构文件的快照。
- `review_packet.json`、`report.html`：审阅材料、来源、差异和限制。
- `round_manifest.json`：冻结文件 SHA256，以及准备脚本和原始输入的摘要。
- `decisions/decision_0001.json` 等：按顺序保存的决定，绑定本轮清单与前一事件哈希。

再次决策前会重新检查冻结文件。证据改变、事件中间缺失、身份或链不匹配会拒绝继续；已有轮次目录不能覆盖。文件哈希和链是本地完整性检查，**不是数字签名或不可篡改审计系统**：有写权限的人可以整体重写目录；没有外部保存的链头时，也不能检测末尾历史被整体删除。重要版本应另行备份并保留独立摘要。

## 如何解释风险与收益

本版展示来源、缺失状态、同靶点相对母体的算术差异，帮助人逐项检查取舍。它不把分数变化自动标成“收益”或“毒性增加”，也不按总分排名。ADMET 端点名称不同或缺失时不计算该差值；即使名称相同，也不自动证明模型、单位、协议和适用域相同。

`comparable_protocols_verified=false` 提醒审阅者核查可比性。对接差值不是实测亲和力差，单种子 RMSD 也不应用作不同体系的亲和力排名。真正的自动风险收益建议、校准与完整计算调度仍需后续实现和评估。
