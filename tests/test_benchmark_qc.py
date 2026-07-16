"""QC benchmark tests with an independently-counted oracle.

The cleaned set is built by dropping a FIXED, RECORDED number of reads per source
class, then the benchmark is asserted to recover those exact hand-computed rates.
The oracle numbers come from the construction here, not from re-running the
benchmark's own classification.
"""

import gzip
from pathlib import Path

import pytest

from viroforge.benchmarking import benchmark_qc, read_labels, read_names
from viroforge.benchmarking.parsers import strip_mate


def _write_fastq(path, records, gz=False):
    """records: list of (name, source, is_dup). Writes a 4-line-per-read FASTQ."""
    lines = []
    for name, source, is_dup in records:
        tags = f"source={source}"
        if is_dup:
            tags += " pcr_duplicate=true"
        lines.append(f"@{name} {tags}")
        lines.append("ACGT" * 30)
        lines.append("+")
        lines.append("I" * 120)
    text = "\n".join(lines) + "\n"
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        Path(path).write_text(text)


def _make_raw():
    """Known composition: 100 viral, 40 host, 20 rrna, 10 phix, 10 low-complexity,
    plus 15 PCR duplicates (10 viral-derived, 5 host-derived)."""
    recs = []
    for i in range(100):
        recs.append((f"vir{i}_0/1", "viral", False))
    for i in range(40):
        recs.append((f"host{i}_0/1", "host_dna", False))
    for i in range(20):
        recs.append((f"rrna{i}_0/1", "rrna", False))
    for i in range(10):
        recs.append((f"phix{i}_0/1", "phix", False))
    for i in range(10):
        recs.append((f"lc{i}_0/1", "artifact_low_complexity", False))
    for i in range(10):
        recs.append((f"vir{i}_0/1_dup{i}", "viral", True))
    for i in range(5):
        recs.append((f"host{i}_0/1_dup{i}", "host_dna", True))
    return recs


def _cleaned_from(raw_recs, drop):
    """Drop the first `drop[class]` reads of each class; return surviving records.

    Classes: viral/host_dna/rrna/phix/artifact_low_complexity/dup. Records the
    survivors; drop counts are the oracle.
    """
    seen = {}
    kept = []
    for name, source, is_dup in raw_recs:
        cls = "dup" if is_dup else source
        n = seen.get(cls, 0)
        seen[cls] = n + 1
        if n < drop.get(cls, 0):
            continue  # dropped
        kept.append((name, source, is_dup))
    return kept


# Oracle: fixed drop counts.
DROP = {
    "host_dna": 38,               # keep 2  -> removal 0.95
    "rrna": 20,                   # keep 0  -> removal 1.0
    "phix": 10,                   # keep 0  -> removal 1.0
    "artifact_low_complexity": 8, # keep 2  -> removal 0.8
    "viral": 3,                   # keep 97 -> over-filter, retention 0.97
    "dup": 12,                    # keep 3  -> dedup 12/15
}


@pytest.fixture
def dataset(tmp_path):
    raw = _make_raw()
    cleaned = _cleaned_from(raw, DROP)
    rpath = tmp_path / "raw_R1.fastq"
    cpath = tmp_path / "cleaned_R1.fastq"
    _write_fastq(rpath, raw)
    _write_fastq(cpath, cleaned)
    return rpath, cpath


def test_oracle_metrics(dataset):
    rpath, cpath = dataset
    r = benchmark_qc(read_labels(rpath), read_names(cpath))

    assert r["match_rate"] == 1.0
    assert r["reliable"] is True

    rt = r["removal_rate_by_type"]
    assert rt["host_dna"]["removal_rate"] == pytest.approx(0.95)
    assert rt["rrna"]["removal_rate"] == pytest.approx(1.0)
    assert rt["phix"]["removal_rate"] == pytest.approx(1.0)
    assert rt["artifact_low_complexity"]["removal_rate"] == pytest.approx(0.8)

    assert r["viral_retention"] == pytest.approx(0.97)
    assert r["over_filtering_rate"] == pytest.approx(0.03)

    # dedup is orthogonal: duplicates of any source, not just viral
    assert r["dedup"]["n_duplicates"] == 15
    assert r["dedup"]["dedup_rate"] == pytest.approx(12 / 15)

    # aggregate confusion matrix over non-duplicate reads
    c = r["contamination"]
    assert (c["tp"], c["fp"], c["fn"], c["tn"]) == (76, 3, 4, 97)
    assert c["precision"] == pytest.approx(76 / 79)
    assert c["recall"] == pytest.approx(76 / 80)


def test_mate_suffix_stripped_still_matches(tmp_path):
    raw = _make_raw()
    cleaned = _cleaned_from(raw, DROP)
    # cleaned tool dropped the trailing /1 mate suffix from names
    cleaned = [(strip_mate(name), s, d) for name, s, d in cleaned]
    _write_fastq(tmp_path / "raw_R1.fastq", raw)
    _write_fastq(tmp_path / "cleaned_R1.fastq", cleaned)
    r = benchmark_qc(read_labels(tmp_path / "raw_R1.fastq"),
                     read_names(tmp_path / "cleaned_R1.fastq"))
    assert r["match_rate_after_mate_strip"] is True
    assert r["match_rate"] == pytest.approx(1.0)
    assert r["viral_retention"] == pytest.approx(0.97)


def test_low_match_rate_flags_unreliable(tmp_path):
    raw = _make_raw()
    # cleaned reads are renamed so nothing maps back
    cleaned = [(f"RENAMED_{i}", s, d) for i, (n, s, d) in enumerate(_cleaned_from(raw, DROP))]
    _write_fastq(tmp_path / "raw_R1.fastq", raw)
    _write_fastq(tmp_path / "cleaned_R1.fastq", cleaned)
    r = benchmark_qc(read_labels(tmp_path / "raw_R1.fastq"),
                     read_names(tmp_path / "cleaned_R1.fastq"))
    assert r["reliable"] is False
    assert r["match_rate"] < 0.5


def test_gzip_cleaned_input(tmp_path):
    raw = _make_raw()
    cleaned = _cleaned_from(raw, DROP)
    _write_fastq(tmp_path / "raw_R1.fastq", raw)
    _write_fastq(tmp_path / "cleaned_R1.fastq.gz", cleaned, gz=True)
    r = benchmark_qc(read_labels(tmp_path / "raw_R1.fastq"),
                     read_names(tmp_path / "cleaned_R1.fastq.gz"))
    assert r["match_rate"] == pytest.approx(1.0)
    assert r["contamination"]["tp"] == 76


def test_unknown_source_bucketed(tmp_path):
    raw = _make_raw() + [("weird0_0/1", "mystery", False)]
    cleaned = _cleaned_from(_make_raw(), DROP)  # weird read is removed
    _write_fastq(tmp_path / "raw_R1.fastq", raw)
    _write_fastq(tmp_path / "cleaned_R1.fastq", cleaned)
    r = benchmark_qc(read_labels(tmp_path / "raw_R1.fastq"),
                     read_names(tmp_path / "cleaned_R1.fastq"))
    assert r["unknown_source"]["n"] == 1
    assert r["unknown_source"]["removal_rate"] == pytest.approx(1.0)
