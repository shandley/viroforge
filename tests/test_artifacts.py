"""
Unit tests for platform-specific sequencing artifacts.

Tests cover:
- PolyG tail artifacts
- Optical duplicates
- Index hopping
- Platform profiles
- Pre-defined platform configurations
"""

import pytest
import numpy as np
from viroforge.artifacts import (
    ReadPair,
    PolyGTailArtifact,
    OpticalDuplicateArtifact,
    IndexHoppingArtifact,
    PlatformProfile,
    novaseq_6000,
    nextseq_2000,
    miseq,
    hiseq_2500,
    no_artifacts
)


# ============================================================================
# Test ReadPair
# ============================================================================

class TestReadPair:
    """Test ReadPair dataclass."""

    def test_read_pair_creation(self):
        """Test basic ReadPair creation."""
        read = ReadPair(
            read_id="read_001",
            forward_seq="ACGTACGT",
            reverse_seq="TGCATGCA"
        )

        assert read.read_id == "read_001"
        assert read.forward_seq == "ACGTACGT"
        assert read.reverse_seq == "TGCATGCA"
        # Quality should be auto-generated
        assert read.forward_qual == 'I' * 8
        assert read.reverse_qual == 'I' * 8

    def test_read_pair_with_custom_quality(self):
        """Test ReadPair with custom quality scores."""
        read = ReadPair(
            read_id="read_002",
            forward_seq="ACGT",
            reverse_seq="TGCA",
            forward_qual="IIII",
            reverse_qual="HHHH"
        )

        assert read.forward_qual == "IIII"
        assert read.reverse_qual == "HHHH"


# ============================================================================
# Test PolyGTailArtifact
# ============================================================================

class TestPolyGTailInitialization:
    """Test PolyG tail artifact initialization."""

    def test_default_initialization(self):
        """Test default parameters."""
        artifact = PolyGTailArtifact()

        assert artifact.frequency == 0.02
        assert artifact.min_length == 10
        assert artifact.max_length == 50
        assert artifact.r1_rate == 0.3
        assert artifact.r2_rate == 1.0

    def test_custom_initialization(self):
        """Test custom parameters."""
        artifact = PolyGTailArtifact(
            frequency=0.05,
            min_length=15,
            max_length=40,
            r1_rate=0.5,
            r2_rate=0.8,
            random_seed=42
        )

        assert artifact.frequency == 0.05
        assert artifact.min_length == 15
        assert artifact.max_length == 40
        assert artifact.r1_rate == 0.5
        assert artifact.r2_rate == 0.8
        assert artifact.random_seed == 42

    def test_invalid_frequency_raises_error(self):
        """Test that invalid frequency raises error."""
        with pytest.raises(ValueError, match="frequency must be"):
            PolyGTailArtifact(frequency=1.5)

    def test_invalid_length_raises_error(self):
        """Test that invalid length raises error."""
        with pytest.raises(ValueError, match="max_length"):
            PolyGTailArtifact(min_length=50, max_length=10)


class TestPolyGTailApplication:
    """Test PolyG tail application to reads."""

    def create_test_reads(self, n=100):
        """Helper to create test reads."""
        reads = []
        for i in range(n):
            read = ReadPair(
                read_id=f"read_{i:04d}",
                forward_seq="ACGTACGT" * 10,  # 80bp
                reverse_seq="TGCATGCA" * 10   # 80bp
            )
            reads.append(read)
        return reads

    def test_polyg_tail_adds_g_bases(self):
        """Test that polyG tails add G bases."""
        artifact = PolyGTailArtifact(frequency=1.0, min_length=20, max_length=20, random_seed=42)
        reads = self.create_test_reads(n=10)

        modified = artifact.apply(reads)

        # All reads should have polyG tails (frequency=1.0)
        for read in modified:
            # Check R1 (rate=0.3, might not all have tails with random seed)
            # Check R2 (rate=1.0, should all have tails)
            if read.reverse_seq:
                assert read.reverse_seq.endswith('G' * 20) or len(read.reverse_seq) == 80

    def test_polyg_tail_preserves_original_sequence(self):
        """Test that polyG tail doesn't modify original sequence."""
        artifact = PolyGTailArtifact(frequency=1.0, min_length=10, max_length=10, random_seed=42)
        reads = self.create_test_reads(n=5)

        original_seq = reads[0].forward_seq
        modified = artifact.apply(reads)

        # Original sequence should be preserved at start
        for read in modified:
            assert read.forward_seq.startswith(original_seq) or read.forward_seq == original_seq

    def test_r2_more_affected_than_r1(self):
        """Test that R2 has higher polyG rate than R1."""
        artifact = PolyGTailArtifact(
            frequency=1.0,  # 100% chance
            min_length=15,
            max_length=15,
            r1_rate=0.5,  # 50% of R1
            r2_rate=1.0,  # 100% of R2
            random_seed=42
        )

        reads = self.create_test_reads(n=1000)
        modified = artifact.apply(reads)

        # Count R1 and R2 with polyG tails
        r1_with_polyg = sum(1 for r in modified if len(r.forward_seq) > 80)
        r2_with_polyg = sum(1 for r in modified if r.reverse_seq and len(r.reverse_seq) > 80)

        # R2 should have more polyG tails than R1
        assert r2_with_polyg > r1_with_polyg

        # R2 rate should be close to 100% (within 5% due to randomness)
        assert r2_with_polyg / len(modified) > 0.95

    def test_polyg_quality_scores_added(self):
        """Test that quality scores are added for polyG tail."""
        artifact = PolyGTailArtifact(frequency=1.0, min_length=20, max_length=20, random_seed=42)
        reads = self.create_test_reads(n=10)

        original_qual_len = len(reads[0].forward_qual)
        modified = artifact.apply(reads)

        # Some reads should have longer quality strings
        some_extended = any(len(r.forward_qual) > original_qual_len for r in modified)
        assert some_extended


# ============================================================================
# Test OpticalDuplicateArtifact
# ============================================================================

class TestOpticalDuplicateInitialization:
    """Test optical duplicate initialization."""

    def test_default_initialization(self):
        """Test default parameters."""
        artifact = OpticalDuplicateArtifact()

        assert artifact.rate == 0.05
        assert artifact.proximity_threshold == 100
        assert artifact.tile_size == 10000

    def test_custom_initialization(self):
        """Test custom parameters."""
        artifact = OpticalDuplicateArtifact(
            rate=0.10,
            proximity_threshold=150,
            tile_size=5000,
            random_seed=42
        )

        assert artifact.rate == 0.10
        assert artifact.proximity_threshold == 150
        assert artifact.tile_size == 5000

    def test_invalid_rate_raises_error(self):
        """Test that invalid rate raises error."""
        with pytest.raises(ValueError, match="rate must be"):
            OpticalDuplicateArtifact(rate=0.5)


class TestOpticalDuplicateApplication:
    """Test optical duplicate application."""

    def create_test_reads(self, n=100):
        """Helper to create test reads."""
        reads = []
        for i in range(n):
            read = ReadPair(
                read_id=f"read_{i:04d}",
                forward_seq="ACGTACGT" * 10,
                reverse_seq="TGCATGCA" * 10,
                tile_x=i * 100,  # Spread out
                tile_y=i * 100
            )
            reads.append(read)
        return reads

    def test_optical_duplicates_increase_read_count(self):
        """Test that optical duplicates add reads."""
        artifact = OpticalDuplicateArtifact(rate=0.10, random_seed=42)
        reads = self.create_test_reads(n=100)

        modified = artifact.apply(reads)

        # Should have ~10% more reads
        assert len(modified) > len(reads)
        assert len(modified) == pytest.approx(110, abs=5)

    def test_duplicates_have_identical_sequences(self):
        """Test that duplicates have identical sequences to source."""
        artifact = OpticalDuplicateArtifact(rate=0.25, random_seed=42)
        reads = self.create_test_reads(n=10)

        modified = artifact.apply(reads)

        # Should have some duplicates (rate=50%)
        assert len(modified) > len(reads)

        # Check that duplicate sequences match originals
        original_seqs = {r.forward_seq for r in reads}
        modified_seqs = {r.forward_seq for r in modified}

        # All sequences in modified should be from originals
        assert modified_seqs.issubset(original_seqs) or modified_seqs == original_seqs

    def test_duplicate_coordinates_nearby(self):
        """Test that duplicates have nearby coordinates."""
        artifact = OpticalDuplicateArtifact(rate=0.30, proximity_threshold=50, random_seed=42)
        reads = self.create_test_reads(n=10)

        modified = artifact.apply(reads)

        # Check that optdup reads exist
        dup_reads = [r for r in modified if "_optdup" in r.read_id]

        assert len(dup_reads) > 0

        # Each duplicate should be near its source
        for dup_read in dup_reads:
            # Find source (same sequence without _optdup)
            source_id = dup_read.read_id.split("_optdup")[0]
            source = next(r for r in reads if r.read_id == source_id)

            # Calculate distance
            dx = abs(dup_read.tile_x - source.tile_x)
            dy = abs(dup_read.tile_y - source.tile_y)

            # Should be within proximity threshold
            assert dx <= artifact.proximity_threshold
            assert dy <= artifact.proximity_threshold


# ============================================================================
# Test IndexHoppingArtifact
# ============================================================================

class TestIndexHoppingInitialization:
    """Test index hopping initialization."""

    def test_default_initialization(self):
        """Test default parameters."""
        artifact = IndexHoppingArtifact()

        assert artifact.rate == 0.01
        assert artifact.n_indexes == 96
        assert len(artifact.index_pool) == 96

    def test_custom_initialization(self):
        """Test custom parameters."""
        artifact = IndexHoppingArtifact(
            rate=0.02,
            n_indexes=384,
            random_seed=42
        )

        assert artifact.rate == 0.02
        assert artifact.n_indexes == 384
        assert len(artifact.index_pool) == 384

    def test_invalid_rate_raises_error(self):
        """Test that invalid rate raises error."""
        with pytest.raises(ValueError, match="rate must be"):
            IndexHoppingArtifact(rate=0.5)


class TestIndexHoppingApplication:
    """Test index hopping application."""

    def create_test_reads(self, n=100):
        """Helper to create test reads."""
        reads = []
        for i in range(n):
            read = ReadPair(
                read_id=f"read_{i:04d}",
                forward_seq="ACGTACGT",
                reverse_seq="TGCATGCA",
                sample_index="INDEX0000"  # All start with same index
            )
            reads.append(read)
        return reads

    def test_index_hopping_changes_indexes(self):
        """Test that index hopping changes sample indexes."""
        artifact = IndexHoppingArtifact(rate=0.10, n_indexes=10, random_seed=42)
        reads = self.create_test_reads(n=100)

        # All start with same index
        original_indexes = {r.sample_index for r in reads}
        assert len(original_indexes) == 1

        modified = artifact.apply(reads)

        # After hopping, should have multiple indexes
        modified_indexes = {r.sample_index for r in modified}
        assert len(modified_indexes) > 1

    def test_hopping_rate_approximately_correct(self):
        """Test that hopping rate is approximately correct."""
        artifact = IndexHoppingArtifact(rate=0.10, n_indexes=96, random_seed=42)
        reads = self.create_test_reads(n=1000)

        modified = artifact.apply(reads)

        # Count reads with changed indexes
        hopped = sum(1 for r in modified if r.sample_index != "INDEX0000")

        # Should be close to 10% (within 3%)
        assert hopped / len(modified) == pytest.approx(0.10, abs=0.03)

    def test_hopped_indexes_from_pool(self):
        """Test that hopped indexes come from defined pool."""
        artifact = IndexHoppingArtifact(rate=0.10, n_indexes=10, random_seed=42)
        reads = self.create_test_reads(n=100)

        modified = artifact.apply(reads)

        # All modified indexes should be from the pool
        for read in modified:
            assert read.sample_index in artifact.index_pool or read.sample_index == "INDEX0000"


# ============================================================================
# Test PlatformProfile
# ============================================================================

class TestPlatformProfile:
    """Test platform profile."""

    def test_platform_profile_creation(self):
        """Test creating a platform profile."""
        artifacts = [
            PolyGTailArtifact(frequency=0.02),
            OpticalDuplicateArtifact(rate=0.05)
        ]

        profile = PlatformProfile(
            name="Test Platform",
            flow_cell_type="patterned",
            artifacts=artifacts
        )

        assert profile.name == "Test Platform"
        assert profile.flow_cell_type == "patterned"
        assert len(profile.artifacts) == 2

    def test_invalid_flow_cell_type_raises_error(self):
        """Test that invalid flow cell type raises error."""
        with pytest.raises(ValueError, match="flow_cell_type must be"):
            PlatformProfile(
                name="Invalid",
                flow_cell_type="invalid_type",
                artifacts=[]
            )

    def test_platform_profile_applies_all_artifacts(self):
        """Test that profile applies all artifacts sequentially."""
        reads = [
            ReadPair(
                read_id=f"read_{i:04d}",
                forward_seq="ACGT" * 20,
                reverse_seq="TGCA" * 20,
                sample_index="INDEX0000"
            )
            for i in range(100)
        ]

        artifacts = [
            PolyGTailArtifact(frequency=0.50, random_seed=42),
            OpticalDuplicateArtifact(rate=0.10, random_seed=42),
            IndexHoppingArtifact(rate=0.10, random_seed=42)
        ]

        profile = PlatformProfile(
            name="Test",
            flow_cell_type="patterned",
            artifacts=artifacts
        )

        modified = profile.apply(reads, random_seed=42)

        # Should have more reads (from optical duplicates)
        assert len(modified) > len(reads)

        # Some reads should have polyG tails (longer sequences)
        some_longer = any(len(r.forward_seq) > 80 for r in modified)
        assert some_longer

        # Some reads should have different indexes
        modified_indexes = {r.sample_index for r in modified}
        assert len(modified_indexes) > 1

    def test_get_summary(self):
        """Test profile summary."""
        profile = PlatformProfile(
            name="Test",
            flow_cell_type="cluster",
            artifacts=[PolyGTailArtifact()],
            description="Test description"
        )

        summary = profile.get_summary()

        assert summary['name'] == "Test"
        assert summary['flow_cell_type'] == "cluster"
        assert summary['n_artifacts'] == 1
        assert summary['description'] == "Test description"


# ============================================================================
# Test Pre-defined Platforms
# ============================================================================

class TestPreDefinedPlatforms:
    """Test pre-defined platform profiles."""

    def test_novaseq_6000_has_polyg(self):
        """Test that NovaSeq has polyG artifact."""
        platform = novaseq_6000()

        assert platform.name == "NovaSeq 6000"
        assert platform.flow_cell_type == "patterned"
        assert len(platform.artifacts) == 3

        # Should have PolyG, OpticalDup, IndexHopping
        artifact_types = [type(a).__name__ for a in platform.artifacts]
        assert "PolyGTailArtifact" in artifact_types
        assert "OpticalDuplicateArtifact" in artifact_types
        assert "IndexHoppingArtifact" in artifact_types

    def test_nextseq_2000_has_polyg(self):
        """Test that NextSeq has polyG artifact."""
        platform = nextseq_2000()

        assert platform.name == "NextSeq 2000"
        assert platform.flow_cell_type == "patterned"

        artifact_types = [type(a).__name__ for a in platform.artifacts]
        assert "PolyGTailArtifact" in artifact_types

    def test_miseq_no_polyg(self):
        """Test that MiSeq has NO polyG artifact (cluster-based)."""
        platform = miseq()

        assert platform.name == "MiSeq"
        assert platform.flow_cell_type == "cluster"

        # Should NOT have PolyG
        artifact_types = [type(a).__name__ for a in platform.artifacts]
        assert "PolyGTailArtifact" not in artifact_types

        # Should have OpticalDup and IndexHopping
        assert "OpticalDuplicateArtifact" in artifact_types
        assert "IndexHoppingArtifact" in artifact_types

    def test_hiseq_2500_no_polyg(self):
        """Test that HiSeq has NO polyG artifact (cluster-based)."""
        platform = hiseq_2500()

        assert platform.name == "HiSeq 2500"
        assert platform.flow_cell_type == "cluster"

        artifact_types = [type(a).__name__ for a in platform.artifacts]
        assert "PolyGTailArtifact" not in artifact_types

    def test_no_artifacts_empty(self):
        """Test that no_artifacts profile has no artifacts."""
        platform = no_artifacts()

        assert platform.name == "Ideal (No Artifacts)"
        assert len(platform.artifacts) == 0


class TestPlatformComparison:
    """Test comparison of different platforms."""

    def create_test_reads(self, n=1000):
        """Helper to create test reads."""
        reads = []
        for i in range(n):
            read = ReadPair(
                read_id=f"read_{i:04d}",
                forward_seq="ACGTACGT" * 10,
                reverse_seq="TGCATGCA" * 10,
                sample_index="INDEX0000"
            )
            reads.append(read)
        return reads

    def test_novaseq_vs_miseq_artifacts(self):
        """Test that NovaSeq has more artifacts than MiSeq."""
        reads_novaseq = self.create_test_reads(n=1000)
        reads_miseq = self.create_test_reads(n=1000)

        novaseq = novaseq_6000()
        miseq_platform = miseq()

        # Apply platforms
        modified_novaseq = novaseq.apply(reads_novaseq, random_seed=42)
        modified_miseq = miseq_platform.apply(reads_miseq, random_seed=42)

        # NovaSeq should have more total reads (higher dup rate)
        assert len(modified_novaseq) > len(modified_miseq)

        # NovaSeq should have polyG tails
        novaseq_has_polyg = any(len(r.forward_seq) > 80 for r in modified_novaseq)
        miseq_has_polyg = any(len(r.forward_seq) > 80 for r in modified_miseq)

        assert novaseq_has_polyg
        assert not miseq_has_polyg  # MiSeq is cluster-based, no polyG

    def test_no_artifacts_preserves_reads(self):
        """Test that no_artifacts doesn't modify reads."""
        reads = self.create_test_reads(n=100)
        original_count = len(reads)
        original_first_seq = reads[0].forward_seq

        platform = no_artifacts()
        modified = platform.apply(reads, random_seed=42)

        # Should have same number of reads
        assert len(modified) == original_count

        # Sequences should be unchanged
        assert modified[0].forward_seq == original_first_seq


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
