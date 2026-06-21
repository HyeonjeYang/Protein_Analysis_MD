# UniProt Sequence Fetching

Use `pamd search-uniprot QUERY --reviewed --organism "Homo sapiens"` to inspect
candidates. Use `pamd fetch-sequence` with `--interactive` or explicit
`--accession` when a query is ambiguous. Cached FASTA and metadata files are
written under `data/sequences/`.
