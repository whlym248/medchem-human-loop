# Local case library / 本地可追溯案例库

This module implements **local JSONL validation and lexical/tag retrieval**. It
does not crawl papers, learn molecular intuition, calculate chemical similarity,
validate an experiment or train a model. The real literature library shipped in
`examples/cases/literature_curated.jsonl` is deliberately empty: no classic drug
optimization cases have been claimed or fabricated to fill it.

中文概述：第一步是把案例存成结构明确、可以回看出处的记录，再按文字和标签查找。
这还不是“训练模型学会科学家直觉”。文献事实需要人阅读原文核对；程序只检查字段是否
齐全，不会因为写着 `reviewed` 就认证其真实性。

## Run without scientific dependencies

From the repository root, using Python 3.10+:

```console
python scripts/case_library.py validate --library examples/cases/literature_curated.jsonl
python scripts/case_library.py validate --library examples/cases/synthetic_cases.jsonl
python scripts/case_library.py search --library examples/cases/synthetic_cases.jsonl --query "aza polarity"
python scripts/case_library.py search --library examples/cases/synthetic_cases.jsonl --query "aza polarity" --allow-synthetic
```

The third command returns **zero matches**. Synthetic examples require an
explicit flag. The fourth returns TOY-002, still labeled synthetic and with
`real_case_support_available: false`. Draft or rejected records are excluded even
with that flag. No calculation or network request occurs. JSON is printed to
standard output; redirect it to your own new output file if desired.

## One JSON object per line

Every record has exactly these fields:

| Field | Contract |
| --- | --- |
| `id` | Unique simple identifier, letters/digits/underscore/dot/hyphen, max 128 characters |
| `scope` | `synthetic_example` or `literature_curated` |
| `parent_smiles`, `product_smiles` | Nonempty strings, max 4096 characters each; syntax, atom mapping, stereochemistry and structure identity are not chemically validated here |
| `edit_tags` | Nonempty array of up to 50 distinct tags, compared case-insensitively |
| `target` | Explicit target description; retrieval's target filter is an exact case-insensitive match |
| `context` | Conditions, hypothesis and applicability boundaries |
| `observations` | Array of source-grounded observations; keep assumptions and unknowns explicit |
| `endpoints` | Array of `{name, value, unit, protocol}`; value is finite numeric or descriptive text, not Boolean |
| `source` | `{url, citation, locator, license_note}`; citation/locator/license note are nonempty |
| `review` | `{status, reviewer, note}`; status is `draft`, `reviewed` or `rejected` |

For literature, `source.url` must be an HTTP(S) source without embedded
credentials; a DOI resolver URL is appropriate. `locator` should identify the
exact table/figure/page and parent/product IDs. A URL is checked for basic format,
not fetched or checked for accessibility. Synthetic records may use `null` URL.
Do not use a paper title alone as the locator or copy restricted full text into
this library. A license note records your assessment; it does not grant rights.

Use `protocol` to retain assay type, target construct, conditions and units
needed for interpretation. If those are unknown, explicitly state what is
unknown instead of inventing them. Preserve measured and predicted endpoints
separately. If one record cannot describe multiple differing protocols clearly,
split it into traceable records. No averaging or cross-protocol conclusion is
performed by the module.

`reviewed` requires a nonempty reviewer identifier and note. A human curator
should read the source, verify structures/endpoints/conditions and record the
scope of review before setting it. The module cannot authenticate the person or
their diligence: every output states `review_attestations_authenticated: false`
and `scientific_claims_verified: false`. A malicious or mistaken supplied record
can pass schema checks. Retrieval therefore supports human review, never replaces it.

## Retrieval and API

```python
from pathlib import Path
from case_library import load_library, validate_library, search_cases

records = load_library(Path("my-reviewed-cases.jsonl"))
validation = validate_library(Path("my-reviewed-cases.jsonl"))
result = search_cases(
    Path("my-reviewed-cases.jsonl"),
    query="aza polarity",
    tags=["small_local"],
    target="MY_TARGET",
    allow_synthetic=False,
    limit=10,
)
```

When imported from the repository root, add `scripts` to your Python module path
or use the CLI. `load_library` returns records and raises `ValueError` for schema
errors. `validate_library` returns counts and the library file's SHA256.
`search_cases` returns `matches`, query/filter parameters, SHA256, total and
returned counts, a truncation flag and explicit limitations.

All query tokens must occur in the combined ID/target/context/tags/observations/
endpoint-name text. Tokens use Unicode `\w+` and case folding: this is not
Chinese word segmentation, stemming, synonym expansion, SMILES similarity or
substructure search. All requested tags must match exactly; an optional target
must also match. A nonempty query, tag or target is required. URLs, reviewer
names, endpoint values and SMILES are not the lexical search corpus. Results
are sorted by case ID, **not scientific merit or confidence**. A maximum of 100
results can be returned per call.

Libraries are local regular files, at most 10 MiB and 10,000 records. The loader
rejects duplicate JSON keys, duplicate case IDs, invalid numeric constants and
incomplete schemas. It does not follow source URLs or execute record content.
The SHA identifies loaded bytes, not their scientific validity.

## Supplying evidence to the route proposal

Use `result["real_case_support_available"]`, not merely `bool(result["matches"])`,
when filling the route module's `case_support` field. The first only counts
returned reviewed records marked `literature_curated`; the second could be true
because of a tutorial fixture. Even the first is an indication of supplied case
records, not proof that their chemistry transfers to the current target.

The next useful curation task is a small, manually verified set of cases with
clear rights and comparable endpoints. Adding actual literature content is
still outstanding. Large-scale ingestion, embeddings, molecular similarity and
model training remain future work.
