# Source acknowledgements

This package retains selected functions and the PoseBusters fork from the supplied
FLOWR research repository: https://github.com/2935361139/flowr_root_analysis .
The repository's MIT license is retained in LICENSE. Extracted functions and
byte-identical files are recorded in SOURCE_PROVENANCE.json.

PoseBusters upstream: https://github.com/maabuu/posebusters .
Its MIT license is retained in posebusters/LICENSE (retrieved from upstream).
The executable Python/YAML files are based on the user's supplied fork; they
are not substituted with the latest upstream implementation. V2 modifies
`posebusters/tools/protein.py` to make environment classes mutually exclusive
and to exclude organic residues from other classes when requested. Numeric
thresholds and the other evaluation modules remain unchanged. This is a locally
modified evaluator, not a claim of upstream endorsement.

The independent input adapters, orchestration, documentation and tests were
prepared for this research workflow. RDKit, NumPy, SciPy, pandas, PyYAML and tqdm
are external dependencies governed by their respective distributions' licenses.
