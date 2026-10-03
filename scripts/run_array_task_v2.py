"""Run exactly one predeclared task through the audited v2 runner."""
import argparse
import json
from pathlib import Path
from run_md_v2 import ROOT, run


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--mode", choices=["pilot", "production"], required=True)
    parser.add_argument("--platform", choices=["CUDA", "OpenCL", "CPU", "Reference"], default="CUDA")
    parser.add_argument("--restart-from-inputs", action="store_true")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    tasks = json.loads((root / "tasks.json").read_text(encoding="utf-8"))[args.mode]
    if not 0 <= args.index < len(tasks):
        parser.error("Task index is outside the frozen manifest")
    item = tasks[args.index]
    prior = root / "runs_v2" / args.mode / item["case_id"] / f"seed_{item['seed']}" / "run_identity.json"
    # Completed tasks still pass through identity/length checks in run().
    result = run(item["case_id"], args.mode, item["seed"], args.platform, item.get("ns", 20),
        resume=prior.exists(), root=root, restart_from_inputs=args.restart_from_inputs)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return result


if __name__ == "__main__":
    main()
