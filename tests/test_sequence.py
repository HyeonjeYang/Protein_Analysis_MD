from __future__ import annotations

from idrptm.sequence import parse_raw_sequence, read_fasta


def test_raw_sequence_parsing_normalizes_whitespace() -> None:
    record = parse_raw_sequence(name="raw seq", sequence="ms tg\n")

    assert record.name == "raw_seq"
    assert record.sequence == "MSTG"


def test_fasta_parsing_reads_single_record(tmp_path) -> None:
    fasta = tmp_path / "sequence.fasta"
    fasta.write_text(">idr fragment\nMS\nTG\n", encoding="utf-8")

    records = read_fasta(fasta)

    assert len(records) == 1
    assert records[0].name == "idr"
    assert records[0].description == "fragment"
    assert records[0].sequence == "MSTG"
