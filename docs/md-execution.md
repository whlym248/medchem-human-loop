# MD 工具的输入、审阅与分析

本页说明工具接口；它不提供可直接开始生产计算的研究数据。请先运行 README 的合成证据示例。MD 执行器处理用户自备的兼容 OpenMM 系统，不负责参数化任意新分子。

## 输入布局

系统包通过 `--root` 指定，与代码仓库可以分开存放：

```text
your-bundle/
  systems/
    YOUR_CASE/
      system.xml
      complex.cif
      positions_nm.npy
      metadata.json
  package_manifest.json       # 可选；存在时校验其中的输入哈希
  tasks.json                  # 数组入口使用的清单，见 scripts/README.md
  inputs/
    charged_ligands/
      YOUR_CASE.json          # 分析所需的匹配配体化学图
  reviews/                    # 审阅操作写入
  runs_v2/                    # MD 执行器写入
```

`YOUR_CASE` 是单层案例标识，不是路径。`metadata.json` 必须与其一致，并包含 `system_sha256`。质量检查和轨迹分析还使用配体原子索引、靶点与配体化学身份等字段；完整要求请按公开源码和 [scripts/README.md](../scripts/README.md) 准备。

输入不是仅靠文件名就有效：XML 的粒子数、CIF 原子顺序、位置数组、配体结构映射和元数据必须匹配。准备系统时还需检查模型的参数来源、质子化、约束、周期性盒和力项是否与执行器假设兼容。

执行器具体要求一个 `MonteCarloBarostat`，并依赖位置约束的全局参数 `k_protein` 与 `k_ligand`。它不能直接运行任意 OpenMM XML；缺失这些约定时，应先准备兼容输入并重新审阅，而不是删去检查。

## 环境

使用隔离环境。`requirements-md.txt` 覆盖执行与 pilot 检查；`requirements-analysis.txt` 增加分析所需库。记录实际 Python、OpenMM、驱动和分析软件版本。依赖清单不是完整环境锁定，也不保证二进制检查点的跨硬件兼容。

执行器提供 CUDA 等平台选项。平台可用不代表输入或结果质量已通过验证。公开整理中的自动化测试不运行 MD。

## 从短测试到生产的顺序

1. 固定输入结构和参数，核查身份与起点。
2. 在预定计算环境运行 `pilot`。
3. 用 `check_pilot_quality.py` 读取真实日志和最终结构，检查数值、温度/密度、几何和速度。
4. 审阅报告及必要的原始数据，解决关键异常后，使用 `review` 明确记录检查结果和理由。
5. 由 `production` 执行器确认审阅对应同一输入和环境，再开始正式采样。
6. 完成后按有效片段清单分析，保存方法、证据和限制。

各命令的具体参数可查看脚本 `--help`。`review` 的四个检查项是声明已实际审阅的结果，不是为绕过检查而默认全部勾选的开关。

## 执行器默认协议

本版本保留原型协议：2 fs 步长、310 K、Langevin Middle 积分、预期 1 bar 压力，以及分阶段位置约束。执行器不会重新设定输入系统中 barostat 的温度和压力；准备系统时必须检查它们与协议一致。技术 pilot 的平衡与采样很短；生产默认包含 0.5 ns 分阶段平衡，生产长度由 `--ns` 指定。

这些参数是代码默认值，不是适用于所有体系的推荐标准。任何协议修改都应独立审阅并重新记录身份，不应静默套用已有审阅或检查点。

## 断点与有效片段

每次执行尝试写入独立目录。输入哈希、平台、种子、协议和目标长度共同构成运行身份。生产断点采用两个槽保存二进制与哈希元数据；恢复时需验证身份、哈希和步数。

同一目录的系统锁用于阻止并发执行。不要通过盲删锁、替换身份 JSON 或忽略哈希错误强行继续。先查明实际进程、有效检查点和原环境。

每次恢复都可能使旧尝试中的部分尾部失效。`segment_manifest.json` 明确指定有效帧前缀；CSV 中的绝对步数是重要依据。分析脚本核对清单及身份，而不是假设所有 DCD 都可完整拼接。

## 分析

`analyze_production_md.py` 读取已存在的数据，不创建 MD Simulation 或启动新的采样。它要求系统与配体化学结构能够映射，并默认拒绝未完成轨迹。

```bash
python scripts/analyze_production_md.py --bundle-root your-bundle --case YOUR_CASE --run-dir your-bundle/runs_v2/production/YOUR_CASE/seed_11 --out-dir analysis/YOUR_CASE_seed11
```

这是接口示例；仓库不含对应真实数据，不能直接复现研究轨迹。分析输出包括帧覆盖信息、RMSD、接触和氢键表、图以及分析规则与软件版本。结果解释见 [scientific-interpretation.md](scientific-interpretation.md)。

可重复提供 `--highlight-residue` 来指定需要关注的残基编号。新的研究体系应检查编号、原子选择和化学映射，不能仅替换案例名就认为分析已完成适配。
