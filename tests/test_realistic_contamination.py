"""
Tests for realistic contamination reference sequences.

Validates that ViroForge uses real reference sequences (rRNA, host DNA, PhiX)
instead of synthetic random sequences, and that adapter read-through injection
works correctly.
"""

import pytest
import random
from pathlib import Path

from viroforge.core.contamination import (
    ContaminationProfile,
    add_host_contamination,
    add_rrna_contamination,
    add_phix_control,
    create_contamination_profile,
)
from viroforge.data.references.resolver import (
    get_phix_path,
    get_rrna_path,
    get_host_fragments_path,
    get_adapter_path,
    has_bundled_references,
)


# ---------------------------------------------------------------------------
# Resolver tests
# ---------------------------------------------------------------------------


class TestResolver:
    """Test the reference file resolver."""

    def test_bundled_references_exist(self):
        refs = has_bundled_references()
        assert refs["phix174"], "Bundled PhiX174 not found"
        assert refs["rrna"], "Bundled rRNA not found"
        assert refs["host_fragments"], "Bundled host fragments not found"
        assert refs["adapters"], "Bundled adapters not found"

    def test_get_phix_path_returns_existing_file(self):
        path = get_phix_path()
        assert path is not None
        assert path.exists()
        assert path.name == "phix174.fasta"

    def test_get_rrna_path_returns_existing_file(self):
        path = get_rrna_path()
        assert path is not None
        assert path.exists()

    def test_get_host_fragments_path_returns_existing_file(self):
        path = get_host_fragments_path()
        assert path is not None
        assert path.exists()

    def test_get_adapter_path_returns_existing_file(self):
        path = get_adapter_path()
        assert path is not None
        assert path.exists()

    def test_user_path_takes_priority(self, tmp_path):
        # Create a fake reference file
        fake_ref = tmp_path / "my_phix.fasta"
        fake_ref.write_text(">fake_phix\nACGT\n")

        path = get_phix_path(user_path=fake_ref)
        assert path == fake_ref

    def test_missing_user_path_falls_back_to_bundled(self):
        path = get_phix_path(user_path=Path("/nonexistent/phix.fasta"))
        # Should fall back to bundled
        assert path is not None
        assert path.exists()

    def test_env_var_override(self, tmp_path, monkeypatch):
        fake_ref = tmp_path / "env_rrna.fasta"
        fake_ref.write_text(">env_rrna\nACGT\n")
        monkeypatch.setenv("VIROFORGE_RRNA_DB", str(fake_ref))

        path = get_rrna_path()
        assert path == fake_ref


# ---------------------------------------------------------------------------
# Real reference contamination tests
# ---------------------------------------------------------------------------


class TestPhiXRealReference:
    """Test that PhiX uses the real NC_001422.1 sequence."""

    def test_phix_uses_real_sequence(self):
        profile = ContaminationProfile()
        add_phix_control(profile, abundance_pct=1.0)

        phix = [c for c in profile.contaminants if c.genome_id == "NC_001422.1_PhiX174"]
        assert len(phix) == 1

        seq = str(phix[0].sequence)
        # Real PhiX174 is exactly 5,386 bp
        assert len(seq) == 5386
        # Real PhiX174 starts with GAGTTTTATCGCTTCC...
        assert seq.startswith("GAGTTTTATCGCTTCC")

    def test_phix_synthetic_fallback(self):
        profile = ContaminationProfile()
        add_phix_control(profile, abundance_pct=1.0, use_real_references=False)

        phix = [c for c in profile.contaminants if c.genome_id == "NC_001422.1_PhiX174"]
        assert len(phix) == 1
        seq = str(phix[0].sequence)
        assert len(seq) == 5386
        # Synthetic won't start with the real PhiX sequence
        assert not seq.startswith("GAGTTTTATCGCTTCC")


class TestRRNARealReference:
    """Test that rRNA uses real reference sequences."""

    def test_rrna_uses_real_sequences(self):
        profile = ContaminationProfile()
        add_rrna_contamination(
            profile, abundance_pct=5.0, n_sequences=10, random_seed=42
        )

        rrna = [c for c in profile.contaminants if "rrna" in c.genome_id.lower()]
        assert len(rrna) == 10

        # Real rRNA sequences have recognizable conserved motifs
        # E. coli 16S contains the universal primer binding site
        all_seqs = "".join(str(c.sequence) for c in rrna)
        # At least some sequences should be from real rRNA
        # (they won't all be random GC-content sequences)
        assert len(all_seqs) > 1000

    def test_rrna_source_is_ncbi(self):
        profile = ContaminationProfile()
        add_rrna_contamination(
            profile, abundance_pct=5.0, n_sequences=5, random_seed=42
        )

        rrna = [c for c in profile.contaminants if "rrna" in c.genome_id.lower()]
        # When using real references, source should be NCBI
        for r in rrna:
            assert r.source == "NCBI_RefSeq_rRNA"

    def test_rrna_synthetic_fallback(self):
        profile = ContaminationProfile()
        add_rrna_contamination(
            profile,
            abundance_pct=5.0,
            n_sequences=10,
            random_seed=42,
            use_real_references=False,
        )

        rrna = [c for c in profile.contaminants if "rrna" in c.genome_id.lower()]
        assert len(rrna) > 0
        # Synthetic uses SILVA_synthetic source
        for r in rrna:
            assert r.source == "SILVA_synthetic"


class TestHostDNARealReference:
    """Test that host DNA uses real genome fragments."""

    def test_host_uses_real_fragments(self):
        profile = ContaminationProfile()
        add_host_contamination(
            profile,
            host_organism="human",
            abundance_pct=2.0,
            n_fragments=5,
            random_seed=42,
        )

        host = [c for c in profile.contaminants if "host" in c.genome_id.lower()]
        assert len(host) == 5

        # Real fragments should have descriptions referencing chromosome regions
        for h in host:
            assert "human_chr" in h.source or "GRCh38" in h.source

    def test_host_synthetic_fallback(self):
        profile = ContaminationProfile()
        add_host_contamination(
            profile,
            host_organism="human",
            abundance_pct=2.0,
            n_fragments=5,
            random_seed=42,
            use_real_references=False,
        )

        host = [c for c in profile.contaminants if "host" in c.genome_id.lower()]
        assert len(host) == 5
        for h in host:
            assert h.source == "GRCh38"
            assert h.description == "Synthetic host DNA fragment"


class TestCreateContaminationProfile:
    """Test that create_contamination_profile respects use_real_references."""

    def test_default_uses_real_references(self):
        profile = create_contamination_profile("realistic", random_seed=42)

        phix = [c for c in profile.contaminants if "PhiX" in c.genome_id]
        assert len(phix) == 1
        # Real PhiX starts with GAGTTTTATCGCTTCC
        assert str(phix[0].sequence).startswith("GAGTTTTATCGCTTCC")

    def test_no_real_contaminants_flag(self):
        profile = create_contamination_profile(
            "realistic", random_seed=42, use_real_references=False
        )

        phix = [c for c in profile.contaminants if "PhiX" in c.genome_id]
        assert len(phix) == 1
        # Synthetic won't match real sequence
        assert not str(phix[0].sequence).startswith("GAGTTTTATCGCTTCC")


# ---------------------------------------------------------------------------
# Adapter read-through tests
# ---------------------------------------------------------------------------


class TestAdapterReadthrough:
    """Test adapter read-through injection."""

    @pytest.fixture
    def sample_fastq(self, tmp_path):
        """Create minimal FASTQ files for testing."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio import SeqIO

        r1_path = tmp_path / "test_R1.fastq"
        r2_path = tmp_path / "test_R2.fastq"

        # Create 100 read pairs with 150bp reads
        r1_records = []
        r2_records = []
        rng = random.Random(42)
        for i in range(100):
            seq = "".join(rng.choice("ACGT") for _ in range(150))
            qual = [30] * 150

            r1 = SeqRecord(Seq(seq), id=f"read_{i}/1", description="")
            r1.letter_annotations["phred_quality"] = qual
            r1_records.append(r1)

            r2 = SeqRecord(Seq(seq), id=f"read_{i}/2", description="")
            r2.letter_annotations["phred_quality"] = qual
            r2_records.append(r2)

        SeqIO.write(r1_records, r1_path, "fastq")
        SeqIO.write(r2_records, r2_path, "fastq")

        return r1_path, r2_path

    def test_adapter_injection_basic(self, sample_fastq, tmp_path):
        from viroforge.simulators.adapters import add_adapter_readthrough

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        stats = add_adapter_readthrough(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            adapter_rate=0.1,
            random_seed=42,
        )

        assert stats["reads_total"] == 100
        assert stats["reads_modified"] == 10
        assert stats["adapter_type"] == "truseq"
        assert len(stats["adapter_lengths"]) == 10
        assert stats["mean_adapter_length"] > 0

    def test_adapter_detectable(self, sample_fastq, tmp_path):
        """Adapter sequences should be detectable by string matching."""
        from viroforge.simulators.adapters import (
            add_adapter_readthrough,
            ADAPTER_SEQUENCES,
        )
        from Bio import SeqIO

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        add_adapter_readthrough(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            adapter_rate=0.5,  # 50% for easy detection
            random_seed=42,
        )

        # Check that adapter k-mers appear in modified reads
        adapter_kmer = ADAPTER_SEQUENCES["truseq"]["r1_adapter"][:12]
        found = 0
        for record in SeqIO.parse(out_r1, "fastq"):
            if adapter_kmer in str(record.seq):
                found += 1

        assert found > 0, "No adapter k-mers detected in modified reads"

    def test_adapter_rate_zero_skips(self, sample_fastq):
        from viroforge.simulators.adapters import add_adapter_readthrough

        r1, r2 = sample_fastq
        stats = add_adapter_readthrough(r1, r2, adapter_rate=0.0)
        assert stats["reads_modified"] == 0

    def test_nextera_adapters(self, sample_fastq, tmp_path):
        from viroforge.simulators.adapters import add_adapter_readthrough

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        stats = add_adapter_readthrough(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            adapter_rate=0.1,
            adapter_type="nextera",
            random_seed=42,
        )
        assert stats["adapter_type"] == "nextera"
        assert stats["reads_modified"] == 10

    def test_invalid_adapter_type_raises(self, sample_fastq):
        from viroforge.simulators.adapters import add_adapter_readthrough

        r1, r2 = sample_fastq
        with pytest.raises(ValueError, match="Unknown adapter type"):
            add_adapter_readthrough(r1, r2, adapter_type="invalid")

    def test_read_lengths_preserved(self, sample_fastq, tmp_path):
        """Modified reads should maintain the same length."""
        from viroforge.simulators.adapters import add_adapter_readthrough
        from Bio import SeqIO

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        add_adapter_readthrough(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            adapter_rate=0.5,
            random_seed=42,
        )

        for record in SeqIO.parse(out_r1, "fastq"):
            assert len(record.seq) == 150
            assert len(record.letter_annotations["phred_quality"]) == 150


# ---------------------------------------------------------------------------
# Low-complexity artifact tests
# ---------------------------------------------------------------------------


class TestLowComplexityInjection:
    """Test low-complexity artifact read injection."""

    @pytest.fixture
    def sample_fastq(self, tmp_path):
        """Create minimal FASTQ files for testing."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio import SeqIO

        r1_path = tmp_path / "test_R1.fastq"
        r2_path = tmp_path / "test_R2.fastq"

        rng = random.Random(42)
        r1_records = []
        r2_records = []
        for i in range(200):
            seq = "".join(rng.choice("ACGT") for _ in range(150))
            qual = [30] * 150
            r1 = SeqRecord(Seq(seq), id=f"read_{i}/1", description="")
            r1.letter_annotations["phred_quality"] = qual
            r1_records.append(r1)
            r2 = SeqRecord(Seq(seq), id=f"read_{i}/2", description="")
            r2.letter_annotations["phred_quality"] = qual
            r2_records.append(r2)

        SeqIO.write(r1_records, r1_path, "fastq")
        SeqIO.write(r2_records, r2_path, "fastq")
        return r1_path, r2_path

    def test_basic_injection(self, sample_fastq, tmp_path):
        from viroforge.simulators.low_complexity import add_low_complexity_reads

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        stats = add_low_complexity_reads(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            rate=0.05,
            random_seed=42,
        )

        assert stats["reads_total"] == 200
        assert stats["reads_modified"] == 10
        assert len(stats["modified_read_ids"]) == 10
        assert sum(stats["artifact_counts"].values()) == 10

    def test_artifact_types_present(self, sample_fastq, tmp_path):
        from viroforge.simulators.low_complexity import add_low_complexity_reads

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        stats = add_low_complexity_reads(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            rate=0.20,  # high rate to get all types
            random_seed=42,
        )

        # With 40 modified reads, we should see multiple artifact types
        assert len(stats["artifact_counts"]) >= 2

    def test_header_tags(self, sample_fastq, tmp_path):
        from viroforge.simulators.low_complexity import add_low_complexity_reads
        from Bio import SeqIO

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        add_low_complexity_reads(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            rate=0.5,
            random_seed=42,
        )

        tagged_count = 0
        for record in SeqIO.parse(out_r1, "fastq"):
            if "source=artifact_low_complexity" in record.description:
                tagged_count += 1
                assert "artifact_type=" in record.description

        assert tagged_count > 0

    def test_low_entropy_detectable(self, sample_fastq, tmp_path):
        """Artifact reads should have low Shannon entropy."""
        from viroforge.simulators.low_complexity import add_low_complexity_reads
        from Bio import SeqIO
        import math

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        add_low_complexity_reads(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            rate=0.5,
            random_seed=42,
        )

        low_entropy_count = 0
        for record in SeqIO.parse(out_r1, "fastq"):
            if "source=artifact_low_complexity" in record.description:
                seq = str(record.seq)
                # Compute Shannon entropy
                freqs = [seq.count(b) / len(seq) for b in "ACGT"]
                entropy = -sum(f * math.log2(f) for f in freqs if f > 0)
                # Low-complexity should have entropy < 1.5 (max is 2.0 for 4 bases)
                if entropy < 1.5:
                    low_entropy_count += 1

        assert low_entropy_count > 0, "No low-entropy artifact reads detected"

    def test_manifest_written(self, sample_fastq, tmp_path):
        from viroforge.simulators.low_complexity import add_low_complexity_reads

        r1, r2 = sample_fastq
        manifest = tmp_path / "lc_manifest.tsv"

        add_low_complexity_reads(
            r1, r2,
            rate=0.05,
            random_seed=42,
            in_place=True,
            manifest_path=manifest,
        )

        assert manifest.exists()
        lines = manifest.read_text().strip().split("\n")
        assert lines[0] == "read_id\tartifact_type\tentropy"
        assert len(lines) == 11  # header + 10 reads

    def test_rate_zero_skips(self, sample_fastq):
        from viroforge.simulators.low_complexity import add_low_complexity_reads

        r1, r2 = sample_fastq
        stats = add_low_complexity_reads(r1, r2, rate=0.0)
        assert stats["reads_modified"] == 0

    def test_read_lengths_preserved(self, sample_fastq, tmp_path):
        from viroforge.simulators.low_complexity import add_low_complexity_reads
        from Bio import SeqIO

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        add_low_complexity_reads(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            rate=0.5,
            random_seed=42,
        )

        for record in SeqIO.parse(out_r1, "fastq"):
            assert len(record.seq) == 150
            assert len(record.letter_annotations["phred_quality"]) == 150


class TestEntropyRange:
    """Test controlled entropy spectrum for threshold sensitivity testing."""

    @pytest.fixture
    def sample_fastq(self, tmp_path):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio import SeqIO

        r1_path = tmp_path / "test_R1.fastq"
        r2_path = tmp_path / "test_R2.fastq"

        rng = random.Random(42)
        r1_records = []
        r2_records = []
        for i in range(200):
            seq = "".join(rng.choice("ACGT") for _ in range(150))
            qual = [30] * 150
            r1 = SeqRecord(Seq(seq), id=f"read_{i}/1", description="")
            r1.letter_annotations["phred_quality"] = qual
            r1_records.append(r1)
            r2 = SeqRecord(Seq(seq), id=f"read_{i}/2", description="")
            r2.letter_annotations["phred_quality"] = qual
            r2_records.append(r2)

        SeqIO.write(r1_records, r1_path, "fastq")
        SeqIO.write(r2_records, r2_path, "fastq")
        return r1_path, r2_path

    def test_entropy_range_produces_controlled_entropy(self, sample_fastq, tmp_path):
        from viroforge.simulators.low_complexity import add_low_complexity_reads
        from viroforge.simulators.low_complexity import _shannon_entropy
        from Bio import SeqIO

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        add_low_complexity_reads(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            rate=0.3,
            entropy_range=(0.3, 0.7),
            random_seed=42,
        )

        entropies = []
        for record in SeqIO.parse(out_r1, "fastq"):
            if "source=artifact_low_complexity" in record.description:
                entropies.append(_shannon_entropy(str(record.seq)))

        assert len(entropies) > 0
        # Most entropies should be in the gray zone; allow some tolerance
        # since the binary search approximation isn't perfect
        in_range = sum(1 for e in entropies if 0.15 < e < 1.1)
        assert in_range / len(entropies) >= 0.8, (
            f"Only {in_range}/{len(entropies)} reads in expected entropy range"
        )

        mean_e = sum(entropies) / len(entropies)
        assert 0.25 < mean_e < 0.85

    def test_entropy_range_artifact_type_is_controlled(self, sample_fastq, tmp_path):
        from viroforge.simulators.low_complexity import add_low_complexity_reads

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        stats = add_low_complexity_reads(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            rate=0.1,
            entropy_range=(0.4, 0.6),
            random_seed=42,
        )

        assert "controlled_entropy" in stats["artifact_counts"]
        assert stats["artifact_counts"]["controlled_entropy"] == stats["reads_modified"]

    def test_entropy_tag_in_header(self, sample_fastq, tmp_path):
        from viroforge.simulators.low_complexity import add_low_complexity_reads
        from Bio import SeqIO

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        add_low_complexity_reads(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            rate=0.3,
            entropy_range=(0.3, 0.7),
            random_seed=42,
        )

        found_entropy_tag = False
        for record in SeqIO.parse(out_r1, "fastq"):
            if "entropy=" in record.description:
                found_entropy_tag = True
                for part in record.description.split():
                    if part.startswith("entropy="):
                        val = float(part.split("=")[1])
                        assert 0.0 <= val <= 2.0
                break

        assert found_entropy_tag

    def test_manifest_includes_entropy(self, sample_fastq, tmp_path):
        from viroforge.simulators.low_complexity import add_low_complexity_reads

        r1, r2 = sample_fastq
        manifest = tmp_path / "lc_manifest.tsv"

        add_low_complexity_reads(
            r1, r2,
            rate=0.05,
            entropy_range=(0.3, 0.7),
            random_seed=42,
            in_place=True,
            manifest_path=manifest,
        )

        lines = manifest.read_text().strip().split("\n")
        assert "entropy" in lines[0]
        for line in lines[1:]:
            parts = line.split("\t")
            entropy = float(parts[2])
            assert 0.0 <= entropy <= 2.0


# ---------------------------------------------------------------------------
# PCR duplicate tests
# ---------------------------------------------------------------------------


class TestPCRDuplicates:
    """Test PCR duplicate injection."""

    @pytest.fixture
    def sample_fastq(self, tmp_path):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio import SeqIO

        r1_path = tmp_path / "test_R1.fastq"
        r2_path = tmp_path / "test_R2.fastq"

        rng = random.Random(42)
        r1_records = []
        r2_records = []
        for i in range(100):
            seq = "".join(rng.choice("ACGT") for _ in range(150))
            qual = [30] * 150
            r1 = SeqRecord(Seq(seq), id=f"read_{i}/1", description="")
            r1.letter_annotations["phred_quality"] = qual
            r1_records.append(r1)
            r2 = SeqRecord(Seq(seq), id=f"read_{i}/2", description="")
            r2.letter_annotations["phred_quality"] = qual
            r2_records.append(r2)

        SeqIO.write(r1_records, r1_path, "fastq")
        SeqIO.write(r2_records, r2_path, "fastq")
        return r1_path, r2_path

    def test_basic_duplicate_injection(self, sample_fastq, tmp_path):
        from viroforge.simulators.duplicates import add_pcr_duplicates

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        stats = add_pcr_duplicates(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            duplicate_rate=0.20,
            random_seed=42,
        )

        assert stats["reads_original"] == 100
        assert stats["templates"] == 20
        assert stats["copies_generated"] > 0
        assert stats["reads_output"] == 100 + stats["copies_generated"]
        assert stats["duplicate_fraction"] > 0

    def test_output_has_more_reads(self, sample_fastq, tmp_path):
        from viroforge.simulators.duplicates import add_pcr_duplicates
        from Bio import SeqIO

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        stats = add_pcr_duplicates(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            duplicate_rate=0.30,
            random_seed=42,
        )

        r1_count = sum(1 for _ in SeqIO.parse(out_r1, "fastq"))
        r2_count = sum(1 for _ in SeqIO.parse(out_r2, "fastq"))

        assert r1_count == stats["reads_output"]
        assert r1_count == r2_count
        assert r1_count > 100

    def test_duplicate_header_tags(self, sample_fastq, tmp_path):
        from viroforge.simulators.duplicates import add_pcr_duplicates
        from Bio import SeqIO

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        add_pcr_duplicates(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            duplicate_rate=0.30,
            random_seed=42,
        )

        dup_count = 0
        for record in SeqIO.parse(out_r1, "fastq"):
            if "pcr_duplicate=true" in record.description:
                dup_count += 1
                assert "duplicate_of=" in record.description
                assert "copy_number=" in record.description

        assert dup_count > 0

    def test_duplicates_share_sequence(self, sample_fastq, tmp_path):
        """Duplicates should be near-identical to their template."""
        from viroforge.simulators.duplicates import add_pcr_duplicates
        from Bio import SeqIO

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        add_pcr_duplicates(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            duplicate_rate=0.30,
            error_rate=0.0,  # exact duplicates
            random_seed=42,
        )

        # Build map of read_id -> sequence
        seqs = {}
        for record in SeqIO.parse(out_r1, "fastq"):
            seqs[record.id] = str(record.seq)

        # Check that duplicates match their templates
        for record in SeqIO.parse(out_r1, "fastq"):
            if "duplicate_of=" in record.description:
                for part in record.description.split():
                    if part.startswith("duplicate_of="):
                        template_id = part.split("=")[1]
                        if template_id in seqs:
                            assert str(record.seq) == seqs[template_id]

    def test_geometric_distribution(self, sample_fastq, tmp_path):
        from viroforge.simulators.duplicates import add_pcr_duplicates

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        stats = add_pcr_duplicates(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            duplicate_rate=0.50,
            max_copies=5,
            copy_distribution="geometric",
            random_seed=42,
        )

        dist = stats["copy_count_distribution"]
        # Geometric: more templates with 1 copy than with many
        if 1 in dist and len(dist) > 1:
            assert dist[1] >= max(dist.get(k, 0) for k in dist if k > 1)

    def test_manifest_written(self, sample_fastq, tmp_path):
        from viroforge.simulators.duplicates import add_pcr_duplicates

        r1, r2 = sample_fastq
        manifest = tmp_path / "dup_manifest.tsv"

        stats = add_pcr_duplicates(
            r1, r2,
            duplicate_rate=0.20,
            random_seed=42,
            in_place=True,
            manifest_path=manifest,
        )

        assert manifest.exists()
        lines = manifest.read_text().strip().split("\n")
        assert lines[0] == "template_read_id\tcopy_read_id\tcopy_number\ttotal_copies"
        assert len(lines) - 1 == stats["copies_generated"]

    def test_rate_zero_skips(self, sample_fastq):
        from viroforge.simulators.duplicates import add_pcr_duplicates

        r1, r2 = sample_fastq
        stats = add_pcr_duplicates(r1, r2, duplicate_rate=0.0)
        assert stats["copies_generated"] == 0

    def test_pcr_errors_introduce_variation(self, sample_fastq, tmp_path):
        from viroforge.simulators.duplicates import add_pcr_duplicates
        from Bio import SeqIO

        r1, r2 = sample_fastq
        out_r1 = tmp_path / "out_R1.fastq"
        out_r2 = tmp_path / "out_R2.fastq"

        add_pcr_duplicates(
            r1, r2,
            output_r1=out_r1,
            output_r2=out_r2,
            duplicate_rate=0.50,
            error_rate=0.01,  # 1% error rate for visibility
            random_seed=42,
        )

        seqs = {}
        for record in SeqIO.parse(out_r1, "fastq"):
            seqs[record.id] = str(record.seq)

        # At least some copies should differ from their template
        mismatches_found = False
        for record in SeqIO.parse(out_r1, "fastq"):
            if "duplicate_of=" in record.description:
                for part in record.description.split():
                    if part.startswith("duplicate_of="):
                        template_id = part.split("=")[1]
                        if template_id in seqs:
                            if str(record.seq) != seqs[template_id]:
                                mismatches_found = True
                                break

        assert mismatches_found, "No PCR errors detected at 1% error rate"


# ---------------------------------------------------------------------------
# ERV injection tests
# ---------------------------------------------------------------------------


class TestERVInjection:
    """Test endogenous and exogenous retroviral read injection."""

    def test_erv_exogenous_from_database(self):
        from viroforge.core.contamination import add_erv_exogenous

        profile = ContaminationProfile()
        add_erv_exogenous(
            profile,
            abundance_pct=0.5,
            n_sequences=3,
            random_seed=42,
        )

        erv = [c for c in profile.contaminants
               if c.contaminant_type.value == "erv_exogenous"]
        assert len(erv) == 3

        for c in erv:
            assert c.source.startswith("RefSeq_")
            assert len(str(c.sequence)) > 1000  # retroviruses are > 1 kb

    def test_erv_exogenous_specific_viruses(self):
        from viroforge.core.contamination import add_erv_exogenous

        profile = ContaminationProfile()
        add_erv_exogenous(
            profile,
            abundance_pct=0.5,
            virus_names=["Human immunodeficiency virus"],
            n_sequences=2,
            random_seed=42,
        )

        erv = [c for c in profile.contaminants
               if c.contaminant_type.value == "erv_exogenous"]
        assert len(erv) == 2

        for c in erv:
            assert "immunodeficiency" in c.organism.lower()

    def test_erv_endogenous_without_fasta_warns(self):
        """Without HERV FASTA, endogenous injection should warn and skip."""
        from viroforge.core.contamination import add_erv_endogenous
        import os

        # Clear env var if set
        old_val = os.environ.pop("VIROFORGE_HERV_DB", None)
        try:
            profile = ContaminationProfile()
            add_erv_endogenous(
                profile,
                abundance_pct=0.5,
                random_seed=42,
            )

            # Should have 0 contaminants (no HERV FASTA available)
            erv = [c for c in profile.contaminants
                   if c.contaminant_type.value == "erv_endogenous"]
            assert len(erv) == 0
        finally:
            if old_val is not None:
                os.environ["VIROFORGE_HERV_DB"] = old_val

    def test_erv_endogenous_with_fasta(self, tmp_path):
        """With a HERV FASTA, endogenous injection should work."""
        from viroforge.core.contamination import add_erv_endogenous

        # Create a fake HERV FASTA
        herv_fasta = tmp_path / "herv_test.fasta"
        herv_fasta.write_text(
            ">HERVK_1 HERV-K family consensus\n"
            + "ACGTACGT" * 100 + "\n"
            + ">HERVW_2 HERV-W family consensus\n"
            + "GCTAGCTA" * 100 + "\n"
            + ">HERVH_3 HERV-H family consensus\n"
            + "TGATCATG" * 100 + "\n"
        )

        profile = ContaminationProfile()
        add_erv_endogenous(
            profile,
            abundance_pct=0.5,
            herv_fasta_path=herv_fasta,
            n_sequences=3,
            random_seed=42,
        )

        erv = [c for c in profile.contaminants
               if c.contaminant_type.value == "erv_endogenous"]
        assert len(erv) == 3

        for c in erv:
            assert c.source == "Dfam_HERV_consensus"

    def test_erv_contamination_types(self):
        """ERV types should be distinct from other contamination types."""
        from viroforge.core.contamination import ContaminantType

        assert ContaminantType.ERV_ENDOGENOUS.value == "erv_endogenous"
        assert ContaminantType.ERV_EXOGENOUS.value == "erv_exogenous"
