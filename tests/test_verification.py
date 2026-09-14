from protaine.verification import (
    check_off_target_homology,
    compute_cai,
    verify_translation,
)


def test_verify_translation_accepts_matching_sequence():
    # AUG AAA GUG UAA -> M K V (stop)
    assert verify_translation("AUGAAAGUGUAA", "MKV") is True


def test_verify_translation_rejects_frameshift():
    # Dropping one nucleotide shifts the reading frame.
    assert verify_translation("AUGAAGUGUAA", "MKV") is False


def test_verify_translation_rejects_premature_stop():
    # A premature stop after the first residue truncates the translation.
    assert verify_translation("AUGUAAGUGUAA", "MKV") is False


def test_verify_translation_ignores_trailing_stop_symbol_in_target():
    assert verify_translation("AUGAAAGUGUAA", "MKV*") is True


def test_compute_cai_of_most_frequent_codons_is_one():
    from protaine.codon_tables import load_codon_table

    table = load_codon_table("human")
    protein = "MKVLAT"
    best_codons = "".join(table.best_codon(aa).codon for aa in protein)
    assert compute_cai(best_codons, organism="human") == 1.0


def test_compute_cai_is_lower_for_rarer_codons():
    from protaine.codon_tables import load_codon_table

    table = load_codon_table("human")
    protein = "MKVLAT"
    best = "".join(table.best_codon(aa).codon for aa in protein)
    worst = "".join(table.synonymous_codons(aa)[-1].codon for aa in protein)
    assert compute_cai(worst, organism="human") <= compute_cai(best, organism="human")


def test_compute_cai_of_empty_sequence_is_zero():
    assert compute_cai("", organism="human") == 0.0


class _FakeHSP:
    def __init__(self, identities, align_length, expect):
        self.identities = identities
        self.align_length = align_length
        self.expect = expect


class _FakeAlignment:
    def __init__(self, accession, hit_def, hsp):
        self.accession = accession
        self.hit_def = hit_def
        self.hsps = [hsp]


class _FakeRecord:
    def __init__(self, alignments):
        self.alignments = alignments


class _FakeHandle:
    def close(self):
        pass


def test_check_off_target_homology_parses_hits(monkeypatch):
    fake_record = _FakeRecord(
        [_FakeAlignment("NM_000000", "Homo sapiens some gene mRNA", _FakeHSP(95, 100, 1e-40))]
    )

    monkeypatch.setattr("Bio.Blast.NCBIWWW.qblast", lambda *a, **k: _FakeHandle())
    monkeypatch.setattr("Bio.Blast.NCBIXML.read", lambda handle: fake_record)

    hits = check_off_target_homology("AUGAAAGUGUAA", organism="human")

    assert len(hits) == 1
    assert hits[0].accession == "NM_000000"
    assert hits[0].identity_pct == 95.0
    assert hits[0].e_value == 1e-40


def test_check_off_target_homology_returns_empty_when_no_hits(monkeypatch):
    monkeypatch.setattr("Bio.Blast.NCBIWWW.qblast", lambda *a, **k: _FakeHandle())
    monkeypatch.setattr("Bio.Blast.NCBIXML.read", lambda handle: _FakeRecord([]))

    assert check_off_target_homology("AUGAAAGUGUAA", organism="human") == []
