#!/usr/bin/env python3
"""Transparent, conservative route proposals. No learned policy and no execution."""
import argparse
import hashlib
import json
from pathlib import Path

REQUIRED = {"edit_scope", "fixed_region_required", "case_support"}
OPTIONAL = {"diffusion_backend_available"}
SCOPES = ("small_local", "fragment_rebuild", "unclear")


def recommend(request):
    """Suggest a review route; caller's case_support is an unverified evidence assertion."""
    if not isinstance(request, dict) or not REQUIRED.issubset(request) or set(request) - REQUIRED - OPTIONAL:
        raise ValueError("request requires edit_scope, fixed_region_required, case_support; "
                         "only optional field is diffusion_backend_available")
    if request["edit_scope"] not in SCOPES:
        raise ValueError("edit_scope must be small_local, fragment_rebuild, or unclear")
    for key in ("fixed_region_required", "case_support", "diffusion_backend_available"):
        if key in request and type(request[key]) is not bool:
            raise ValueError(f"{key} must be a Boolean")
    scope = request["edit_scope"]
    reasons, next_steps = [], []
    if scope == "unclear":
        route = "collect_evidence"
        reasons.append("The editable region and intended chemical change are not sufficiently specified.")
        next_steps.append("Record the binding hypothesis, editable atoms, protected atoms and desired endpoints.")
    elif scope == "fragment_rebuild":
        route = "defer_diffusion"
        reasons.append("Fragment rebuilding is outside the limited rule-edit capability.")
        reasons.append("A constrained diffusion backend has not been implemented or validated in this repository.")
        next_steps.append("Keep the proposal pending until a separately validated backend and review contract exist.")
    elif not request["case_support"]:
        route = "collect_evidence"
        reasons.append("No reviewed, relevant nonsynthetic case support has been asserted.")
        next_steps.append("Curate and review applicable case evidence before requesting a supported local edit.")
    else:
        route = "limited_rule_edits"
        reasons.append("A small local edit and reviewed case support were supplied to this deterministic rule.")
        reasons.append("Only enumerated edits supported by the separate rule-edit module may be proposed.")
        next_steps.append("Inspect original cases and choose an allowed edit; generation still requires its own input validation.")
    if request["fixed_region_required"]:
        reasons.append("Protected atoms/core must be explicitly mapped and checked by the editing module.")
        next_steps.append("Reject products whose protected region, stereochemistry or required bonds change unexpectedly.")
    if not request["case_support"] and scope != "small_local":
        reasons.append("Relevant reviewed nonsynthetic case support is still missing.")
    if "diffusion_backend_available" in request:
        reasons.append("The supplied diffusion availability flag is ignored; a Boolean does not install or validate a backend.")
    next_steps.append("A human must review evidence, chemical validity and tradeoffs before any downstream calculation.")
    return {
        "schema_version": 1, "policy_version": "deterministic_rules_v1",
        "recommendation": route, "reasons": reasons, "required_next_steps": next_steps,
        "blocked_methods": [{"method": "constrained_diffusion", "status": "not_implemented",
                             "reason": "No installed and validated adapter is provided by this release."}],
        "human_approval_required": True, "not_learned_policy": True,
        "executes_generation": False, "executes_MD": False,
        "case_support_verified_by_this_function": False,
        "limitations": "A route is not a predicted benefit, safety decision or scientific approval. "
                       "Case support must come from reviewed nonsynthetic evidence, not demo matches. "
                       "No reward training, RL or molecular-utility model is used.",
        "request": dict(request),
    }


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_constant(value):
    raise ValueError(f"Invalid JSON numeric constant: {value}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        with args.request.open("rb") as handle:
            raw = handle.read(65537)
        if len(raw) > 65536:
            raise ValueError("Request exceeds the 64 KiB limit")
        request = json.loads(raw.decode("utf-8-sig"), object_pairs_hook=_unique_object,
                             parse_constant=_reject_constant)
        result = recommend(request)
        result["request_sha256"] = hashlib.sha256(raw).hexdigest()
    except (ValueError, OSError, UnicodeError) as error:
        parser.error(str(error))
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    return result


if __name__ == "__main__":
    main()
