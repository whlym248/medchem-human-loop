#!/usr/bin/env python3
"""Local, provenance-preserving lexical retrieval; no chemistry model or web I/O."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import re
from urllib.parse import urlsplit

MAX_LIBRARY_BYTES = 10 * 1024 * 1024
MAX_RECORDS = 10000
RECORD_KEYS = {
    "id", "scope", "parent_smiles", "product_smiles", "edit_tags", "target",
    "context", "observations", "endpoints", "source", "review",
}
DISCLAIMER = (
    "Lexical/tag retrieval only, not molecular similarity or learned chemical intuition. "
    "Reviewed is a supplied attestation, not authenticated identity or verified scientific truth. "
    "Keep each endpoint's protocol and source; no pooled effect or optimization conclusion is computed."
)


def _text(value, label, *, nonempty=True, maximum=12000):
    if not isinstance(value, str) or len(value) > maximum or (nonempty and not value.strip()):
        raise ValueError(f"{label} must be a {'nonempty ' if nonempty else ''}string (max {maximum})")
    if any(ord(c) < 32 and c not in "\n\r\t" for c in value):
        raise ValueError(f"{label} contains control characters")
    return value


def _keys(value, expected, label):
    if not isinstance(value, dict) or set(value) != expected:
        raise ValueError(f"{label} requires exactly these fields: {', '.join(sorted(expected))}")


def validate_record(record):
    """Validate schema and provenance completeness, never biological truth or SMILES chemistry."""
    _keys(record, RECORD_KEYS, "record")
    identifier = _text(record["id"], "id", maximum=128)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", identifier):
        raise ValueError("id must be a simple identifier")
    if record["scope"] not in ("synthetic_example", "literature_curated"):
        raise ValueError("scope must be synthetic_example or literature_curated")
    for key in ("parent_smiles", "product_smiles"):
        _text(record[key], key, maximum=4096)
    for key in ("target", "context"):
        _text(record[key], key)
    tags = record["edit_tags"]
    if not isinstance(tags, list) or not tags or len(tags) > 50:
        raise ValueError("edit_tags must be a nonempty list of at most 50 strings")
    for tag in tags:
        _text(tag, "edit tag", maximum=128)
    if len({t.casefold() for t in tags}) != len(tags):
        raise ValueError("edit_tags contains a duplicate")
    observations = record["observations"]
    if not isinstance(observations, list) or len(observations) > 100:
        raise ValueError("observations must be a list of at most 100 strings")
    for observation in observations:
        _text(observation, "observation")
    endpoints = record["endpoints"]
    if not isinstance(endpoints, list) or len(endpoints) > 100:
        raise ValueError("endpoints must be a list of at most 100 endpoint objects")
    for endpoint in endpoints:
        _keys(endpoint, {"name", "value", "unit", "protocol"}, "endpoint")
        for key in ("name", "unit", "protocol"):
            _text(endpoint[key], f"endpoint.{key}")
        value = endpoint["value"]
        if isinstance(value, bool) or not isinstance(value, (str, int, float)):
            raise ValueError("endpoint.value must be text or a finite number; not a Boolean")
        if isinstance(value, str):
            _text(value, "endpoint.value")
        elif isinstance(value, float) and not math.isfinite(value):
            raise ValueError("endpoint.value must be finite")
    source = record["source"]
    _keys(source, {"url", "citation", "locator", "license_note"}, "source")
    for key in ("citation", "locator", "license_note"):
        _text(source[key], f"source.{key}")
    if source["url"] is None:
        if record["scope"] == "literature_curated":
            raise ValueError("literature_curated source.url must be traceable HTTP(S)")
    else:
        url = _text(source["url"], "source.url", maximum=2048)
        parsed = urlsplit(url)
        if parsed.scheme not in ("https", "http") or not parsed.netloc or parsed.username or parsed.password:
            raise ValueError("source.url must be an HTTP(S) URL without credentials")
    review = record["review"]
    _keys(review, {"status", "reviewer", "note"}, "review")
    if review["status"] not in ("draft", "reviewed", "rejected"):
        raise ValueError("review.status must be draft, reviewed, or rejected")
    for key in ("reviewer", "note"):
        _text(review[key], f"review.{key}", nonempty=review["status"] == "reviewed")
    # A name and note can be checked for presence, not identity or authenticity.
    return None


def _reject_constant(value):
    raise ValueError(f"Invalid JSON numeric constant: {value}")


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _load(path):
    path = Path(path)
    if not path.is_file():
        raise ValueError("Library path must identify an existing regular JSONL file")
    if path.stat().st_size > MAX_LIBRARY_BYTES:
        raise ValueError("Library exceeds the 10 MiB limit")
    # Bounded read, including protection against a file growing after stat().
    with path.open("rb") as handle:
        raw = handle.read(MAX_LIBRARY_BYTES + 1)
    if len(raw) > MAX_LIBRARY_BYTES:
        raise ValueError("Library exceeds the 10 MiB limit")
    records, ids = [], set()
    for line_number, line in enumerate(raw.decode("utf-8-sig").splitlines(), 1):
        if not line.strip():
            continue
        if len(records) >= MAX_RECORDS:
            raise ValueError("Library exceeds the 10000-record limit")
        try:
            record = json.loads(line, parse_constant=_reject_constant, object_pairs_hook=_unique_object)
            validate_record(record)
            if record["id"] in ids:
                raise ValueError(f"Duplicate case id: {record['id']}")
        except (ValueError, TypeError) as error:
            raise ValueError(f"Invalid case on JSONL line {line_number}: {error}") from error
        ids.add(record["id"])
        records.append(record)
    return records, hashlib.sha256(raw).hexdigest()


def load_library(path):
    """Read and validate a local JSONL file. Does not follow any record URL or path."""
    return _load(path)[0]


def validate_library(path):
    records, digest = _load(path)
    return {
        "schema_version": 1, "valid": True, "record_count": len(records),
        "library_sha256": digest,
        "scope_counts": {s: sum(r["scope"] == s for r in records)
                         for s in ("synthetic_example", "literature_curated")},
        "review_status_counts": {s: sum(r["review"]["status"] == s for r in records)
                                 for s in ("draft", "reviewed", "rejected")},
        "review_attestations_authenticated": False, "scientific_claims_verified": False,
        "limitations": DISCLAIMER,
    }


def _tokens(text):
    return set(re.findall(r"\w+", text.casefold(), flags=re.UNICODE))


def search_cases(path, query="", *, tags=(), target=None, allow_synthetic=False, limit=10):
    """Find reviewed records by exact lexical tokens, edit tags and optional exact target.

    All query tokens and all requested tags must match. No hit means no support;
    results are sorted by case ID, not by potency, relevance probability or efficacy.
    """
    _text(query, "query", nonempty=False, maximum=1000)
    if not isinstance(tags, (list, tuple)) or len(tags) > 50:
        raise ValueError("tags must be a list or tuple of at most 50 tags")
    for tag in tags:
        _text(tag, "query tag", maximum=128)
    if target is not None:
        _text(target, "target", maximum=1000)
    if type(allow_synthetic) is not bool:
        raise ValueError("allow_synthetic must be a Boolean")
    if type(limit) is not int or not 1 <= limit <= 100:
        raise ValueError("limit must be an integer from 1 to 100")
    query_tokens = _tokens(query)
    if not query_tokens and not tags and target is None:
        raise ValueError("Provide query tokens, at least one tag, or an exact target")
    records, digest = _load(path)
    matches = []
    wanted_tags = {t.casefold() for t in tags}
    for record in records:
        if record["review"]["status"] != "reviewed":
            continue
        if record["scope"] == "synthetic_example" and not allow_synthetic:
            continue
        if target is not None and record["target"].casefold() != target.casefold():
            continue
        if not wanted_tags.issubset({t.casefold() for t in record["edit_tags"]}):
            continue
        # Do not match URLs, reviewer names or endpoint numbers as scientific support.
        search_text = " ".join([
            record["id"], record["target"], record["context"],
            *record["edit_tags"], *record["observations"],
            *[e["name"] for e in record["endpoints"]],
        ])
        if query_tokens.issubset(_tokens(search_text)):
            matches.append(copy.deepcopy(record))
    matches.sort(key=lambda r: r["id"])
    selected = matches[:limit]
    return {
        "schema_version": 1, "library_sha256": digest,
        "query": query, "tags": list(tags), "target": target,
        "allow_synthetic": allow_synthetic, "match_count": len(matches),
        "returned_count": len(selected), "truncated": len(matches) > limit,
        "matches": selected,
        "real_case_support_available": any(r["scope"] == "literature_curated" for r in selected),
        "method": "lexical_all_tokens_and_exact_tags_and_exact_target",
        "sort": "case_id_not_scientific_rank",
        "review_attestations_authenticated": False, "scientific_claims_verified": False,
        "limitations": DISCLAIMER,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="command", required=True)
    check = actions.add_parser("validate", help="Validate schema/provenance, not scientific truth")
    check.add_argument("--library", required=True, type=Path)
    search = actions.add_parser("search", help="Retrieve reviewed local cases")
    search.add_argument("--library", required=True, type=Path)
    search.add_argument("--query", default="")
    search.add_argument("--tag", action="append", default=[])
    search.add_argument("--target")
    search.add_argument("--allow-synthetic", action="store_true")
    search.add_argument("--limit", default=10, type=int)
    args = parser.parse_args(argv)
    try:
        result = (validate_library(args.library) if args.command == "validate" else
                  search_cases(args.library, args.query, tags=args.tag, target=args.target,
                               allow_synthetic=args.allow_synthetic, limit=args.limit))
    except (ValueError, OSError, UnicodeError) as error:
        parser.error(str(error))
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    return result


if __name__ == "__main__":
    main()
