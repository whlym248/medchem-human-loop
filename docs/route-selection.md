# Transparent route proposal / 透明方法路线建议

`scripts/route_proposal.py` now implements a small **deterministic decision
table**. It suggests whether to collect evidence, consider the separate limited
rule-edit module, or defer a fragment-rebuilding proposal. It does not choose a
clinical candidate, run generation or launch MD. It is not reinforcement
learning and does not contain a learned reward/utility model.

中文概述：我们先把“为何建议这种方法”做成可以检查的规则。当前局部约束 diffusion
没有实现，因此程序不能通过一个 `true` 参数假装它已经可用。规则输出仍需要人检查
证据、可改部位、保留区域和风险收益。

## Input and CLI

```console
python scripts/route_proposal.py --request examples/route_request.json
```

The bundled input truthfully declares no real case support and returns
`collect_evidence`. Output goes to standard output and includes the input file's
SHA256. JSON files are capped at 64 KiB; duplicate keys are rejected.

```json
{
  "edit_scope": "small_local",
  "fixed_region_required": true,
  "case_support": false
}
```

Required fields:

- `edit_scope`: `small_local`, `fragment_rebuild` or `unclear`.
- `fixed_region_required`: a JSON Boolean; if true, the proposal requires an
  explicit protected-atom/core mapping and downstream verification.
- `case_support`: a JSON Boolean supplied by the caller, based on relevant,
  reviewed **nonsynthetic** records. The route function does not inspect the
  underlying cases itself; the round coordinator can retrieve them first.

Optional `diffusion_backend_available` is accepted only as a Boolean and is
always ignored for capability purposes. This release has no validated diffusion
adapter. Unknown fields or truthy substitutes such as `"yes"`/`1` are rejected.

## Rules

| Condition | `recommendation` | Meaning |
| --- | --- | --- |
| Scope unclear | `collect_evidence` | Specify the hypothesis and editable/protected regions |
| Fragment rebuild | `defer_diffusion` | The proposed scale exceeds the limited editing module; validated diffusion support is absent |
| Small local edit, no real case support asserted | `collect_evidence` | Review relevant source evidence first |
| Small local edit, real case support asserted | `limited_rule_edits` | Consider only supported enumerated edits after human review |

If a fixed region is required, all routes add an explicit protected-region
check. This does **not** mean that this function receives atom maps, enforces a
constraint or validates stereochemistry. Those checks belong to the chemical
editing step. A suggested local edit is not guaranteed to be accepted by that
module: it may lack a supported operation or have invalid inputs.

All outputs include human-readable `reasons` and `required_next_steps`, plus:

```json
{
  "policy_version": "deterministic_rules_v1",
  "human_approval_required": true,
  "not_learned_policy": true,
  "executes_generation": false,
  "executes_MD": false,
  "case_support_verified_by_this_function": false,
  "blocked_methods": [
    {
      "method": "constrained_diffusion",
      "status": "not_implemented",
      "reason": "No installed and validated adapter is provided by this release."
    }
  ]
}
```

## Python integration

```python
from route_proposal import recommend
from case_library import search_cases

cases = search_cases("my-reviewed-cases.jsonl", query="aza polarity")
proposal = recommend({
    "edit_scope": "small_local",
    "fixed_region_required": True,
    "case_support": cases["real_case_support_available"],
})
```

These modules live in `scripts`; use the CLI or make that directory importable.
`recommend(request) -> dict` has no side effects and does not mutate the caller's
object. It returns the request with its reasoning so a review packet can retain
the decision context. Synthetic examples always remain teaching inputs, not a
shortcut to claiming real supporting literature.

No affinity, toxicity, ADMET or overall value is scored by these rules. A person
must still decide whether the proposed experiment is worth doing; calling the
function is not their approval. Diffusion integration, predictive comparison of
methods and any learned policy remain future work and require separate evidence
and testing.
