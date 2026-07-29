"""Tests for bacterial background, the sample's own microbiome.

Bulk metagenomes (``--no-vlp``) were unrealistic without this: measured on
v0.16.0, an IBD gut run at ``--contamination-level heavy`` came out 81.6%
viral-origin where real bulk stool is 1-5% viral and 60-80% bacterial. No
setting could fix it, because maxing every existing contamination knob on that
collection tops out at 34.1% and none of it is bacterial. See issue #37.

These tests use an independent oracle: expected values are computed here from
the stated model rather than read back from the implementation.
"""

import pytest

from viroforge.benchmarking.qc import DEFAULT_KEEP_REMOVE
from viroforge.core.contamination import (
    BACTERIAL_COMMUNITY_PROFILES,
    ContaminantType,
    ContaminationProfile,
    add_bacterial_background,
    add_reagent_contamination,
)
from viroforge.enrichment.vlp import VLPEnrichment, VLPProtocol


def _profile(**kwargs):
    p = ContaminationProfile(name="test")
    defaults = dict(abundance_pct=70.0, n_fragments=20, fragment_length=1000,
                    random_seed=42)
    defaults.update(kwargs)
    add_bacterial_background(p, **defaults)
    return p


class TestOptIn:
    """Nothing changes for existing callers until they ask for it."""

    def test_default_adds_nothing(self):
        p = ContaminationProfile(name="t")
        add_bacterial_background(p)
        assert p.contaminants == []

    def test_zero_pct_adds_nothing(self):
        p = ContaminationProfile(name="t")
        add_bacterial_background(p, abundance_pct=0.0, n_fragments=10)
        assert p.contaminants == []

    def test_rejects_zero_fragments(self):
        p = ContaminationProfile(name="t")
        with pytest.raises(ValueError, match="n_fragments"):
            add_bacterial_background(p, abundance_pct=50.0, n_fragments=0)


class TestAbundance:
    def test_total_abundance_matches_requested_pct(self):
        p = _profile(abundance_pct=70.0)
        assert sum(c.abundance for c in p.contaminants) == pytest.approx(0.70)

    def test_abundance_split_evenly(self):
        n = 20
        p = _profile(abundance_pct=50.0, n_fragments=n)
        assert len(p.contaminants) == n
        for c in p.contaminants:
            assert c.abundance == pytest.approx(0.50 / n)


class TestTypeIsDistinctFromKitome:
    """The whole point of a separate type: the kitome baseline is cited."""

    def test_tagged_as_bacterial_background(self):
        p = _profile()
        assert {c.contaminant_type for c in p.contaminants} == {
            ContaminantType.BACTERIAL_BACKGROUND}

    def test_not_conflated_with_reagent_bacteria(self):
        assert (ContaminantType.BACTERIAL_BACKGROUND
                is not ContaminantType.REAGENT_BACTERIA)
        assert (ContaminantType.BACTERIAL_BACKGROUND.value
                != ContaminantType.REAGENT_BACTERIA.value)

    def test_both_can_coexist_and_stay_separable(self):
        p = ContaminationProfile(name="t")
        add_reagent_contamination(p, abundance_pct=0.5, n_genomes=2,
                                  random_seed=1)
        add_bacterial_background(p, abundance_pct=70.0, n_fragments=10,
                                 fragment_length=500, random_seed=1)
        by_type = {}
        for c in p.contaminants:
            by_type.setdefault(c.contaminant_type, []).append(c)
        assert len(by_type[ContaminantType.REAGENT_BACTERIA]) == 2
        assert len(by_type[ContaminantType.BACTERIAL_BACKGROUND]) == 10


class TestReproducibility:
    def _seqs(self, seed):
        return [str(c.sequence) for c in _profile(random_seed=seed).contaminants]

    def test_same_seed_reproducible(self):
        assert self._seqs(42) == self._seqs(42)

    def test_different_seed_differs(self):
        assert self._seqs(42) != self._seqs(7)


class TestCommunityProfiles:
    def test_gc_tracks_the_community(self):
        """Soil bacteria are GC-rich; gut and marine are not."""
        means = {}
        for community in ("gut", "soil", "marine"):
            p = _profile(community_type=community, n_fragments=40)
            gcs = [c.gc_content for c in p.contaminants]
            means[community] = sum(gcs) / len(gcs)

        for community, mean in means.items():
            expected = BACTERIAL_COMMUNITY_PROFILES[community]["gc_content"]
            assert mean == pytest.approx(expected, abs=3.0)

        assert means["soil"] > means["gut"] > means["marine"]

    def test_organisms_come_from_the_community(self):
        p = _profile(community_type="soil")
        expected = BACTERIAL_COMMUNITY_PROFILES["soil"]["genera"]
        for c in p.contaminants:
            assert any(g in c.organism for g in expected), c.organism

    def test_unknown_community_falls_back_to_gut(self):
        p = _profile(community_type="not_a_real_site")
        gut = BACTERIAL_COMMUNITY_PROFILES["gut"]["genera"]
        assert all(any(g in c.organism for g in gut) for c in p.contaminants)


class TestVLPRemoval:
    """VLP's main job in reality is removing this. It must not take the
    generic 'OTHER or unknown' path, which removes only 50-99%."""

    def test_filtration_removes_about_98_percent(self):
        p = _profile(abundance_pct=70.0, n_fragments=30)
        before = p.get_total_abundance()

        vlp = VLPEnrichment(protocol=VLPProtocol.tangential_flow_standard(),
                            random_seed=42)
        after_profile, _ = vlp.apply_contamination_reduction(p)
        removal = 1.0 - after_profile.get_total_abundance() / before

        # 0.2 um filtration against 1-5 um cells, same as reagent bacteria
        assert removal > 0.95, f"only {removal:.1%} removed; generic fallthrough?"

    def test_removed_like_reagent_bacteria_not_like_other(self):
        """Both bacterial types are whole cells, so they filter the same."""
        def removal_for(adder):
            p = ContaminationProfile(name="t")
            adder(p)
            before = p.get_total_abundance()
            vlp = VLPEnrichment(
                protocol=VLPProtocol.tangential_flow_standard(), random_seed=7)
            after, _ = vlp.apply_contamination_reduction(p)
            return 1.0 - after.get_total_abundance() / before

        bg = removal_for(lambda p: add_bacterial_background(
            p, abundance_pct=70.0, n_fragments=20, fragment_length=500,
            random_seed=1))
        kit = removal_for(lambda p: add_reagent_contamination(
            p, abundance_pct=0.5, n_genomes=20, random_seed=1))

        assert bg == pytest.approx(kit, abs=0.05)

    def test_no_vlp_leaves_it_in_place(self):
        p = _profile(abundance_pct=70.0, n_fragments=20)
        before = p.get_total_abundance()
        bulk = VLPEnrichment(protocol=VLPProtocol.no_vlp(), random_seed=42)
        after, _ = bulk.apply_contamination_reduction(p)
        assert after.get_total_abundance() == pytest.approx(before, rel=0.01)


class TestQCPolicy:
    def test_removed_by_default_for_virome_qc(self):
        assert DEFAULT_KEEP_REMOVE["bacterial_background"] == "remove"

    def test_distinct_policy_entry_from_reagent(self):
        assert "reagent_bacteria" in DEFAULT_KEEP_REMOVE
        assert "bacterial_background" in DEFAULT_KEEP_REMOVE
