from __future__ import annotations

import pytest

from idrptm.uniprot import AmbiguousUniProtQuery, fetch_sequence, search_uniprot


class FakeResponse:
    def __init__(self, payload):
        self.payload = payload

    def raise_for_status(self) -> None:
        return None

    def json(self):
        return self.payload


class FakeSession:
    def __init__(self, payload):
        self.payload = payload
        self.calls = []

    def get(self, url, **kwargs):
        self.calls.append((url, kwargs))
        return FakeResponse(self.payload)


def _payload():
    return {
        "results": [
            {
                "primaryAccession": "P00001",
                "uniProtkbId": "AAA_HUMAN",
                "entryType": "UniProtKB reviewed (Swiss-Prot)",
                "proteinDescription": {
                    "recommendedName": {"fullName": {"value": "Alpha protein"}}
                },
                "genes": [{"geneName": {"value": "AAA"}}],
                "organism": {"scientificName": "Homo sapiens"},
                "sequence": {"length": 4, "value": "AAAA"},
            },
            {
                "primaryAccession": "P00002",
                "uniProtkbId": "BBB_HUMAN",
                "entryType": "UniProtKB reviewed (Swiss-Prot)",
                "proteinDescription": {
                    "recommendedName": {"fullName": {"value": "Beta protein"}}
                },
                "genes": [{"geneName": {"value": "BBB"}}],
                "organism": {"scientificName": "Homo sapiens"},
                "sequence": {"length": 4, "value": "BBBB"},
            },
        ]
    }


def test_search_uniprot_uses_mock_response() -> None:
    candidates = search_uniprot("AAA", organism="Homo sapiens", session=FakeSession(_payload()))

    assert candidates[0].accession == "P00001"
    assert candidates[0].sequence == "AAAA"


def test_fetch_sequence_requires_selection_for_ambiguous_query(tmp_path) -> None:
    with pytest.raises(AmbiguousUniProtQuery):
        fetch_sequence(
            "protein",
            organism="Homo sapiens",
            cache_dir=tmp_path,
            session=FakeSession(_payload()),
        )


def test_fetch_sequence_writes_cache_with_accession(tmp_path) -> None:
    entry_payload = _payload()["results"][0]
    result = fetch_sequence(
        "AAA",
        accession="P00001",
        organism="Homo sapiens",
        cache_dir=tmp_path,
        session=FakeSession(entry_payload),
    )

    assert result.fasta_path.exists()
    assert result.metadata_path.exists()
    assert "AAAA" in result.fasta_path.read_text(encoding="utf-8")
