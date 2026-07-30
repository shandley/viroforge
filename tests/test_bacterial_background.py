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


class TestCurationScriptAgreesWithModule:
    """The reference set and the runtime model must name the same genera.

    add_bacterial_background() matches fragments to a community by genus
    substring against the FASTA description. If the curation script fetches
    genera the module does not list, those fragments never match their own
    community and silently fall back to sampling across all taxa. If the module
    lists genera the script never fetches, that community has no real reference.
    Neither failure is visible at runtime, so it is pinned here.
    """

    @staticmethod
    def _script_taxa():
        import importlib.util
        from pathlib import Path

        path = Path(__file__).parent.parent / "scripts" / "curate_bacterial_background.py"
        spec = importlib.util.spec_from_file_location("cbb", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod.COMMUNITY_TAXA

    def test_same_communities(self):
        assert set(self._script_taxa()) == set(BACTERIAL_COMMUNITY_PROFILES)

    def test_genera_match_exactly(self):
        taxa = self._script_taxa()
        for community, names in taxa.items():
            fetched = {n.split()[0] for n in names}
            modelled = set(BACTERIAL_COMMUNITY_PROFILES[community]["genera"])
            assert fetched == modelled, (
                f"{community}: curation fetches {sorted(fetched - modelled)} "
                f"that the model does not list, and the model lists "
                f"{sorted(modelled - fetched)} with no reference sequence"
            )

    def test_every_taxon_is_binomial(self):
        """Genus-only entries would match too broadly on substring."""
        for community, names in self._script_taxa().items():
            for name in names:
                assert len(name.split()) >= 2, f"{community}: {name!r}"


class TestNestedAbundanceSemantics:
    """Contamination fractions are fractions of the FINAL read set.

    The viral community is scaled to fill 1 - contamination. Before bacterial
    background existed both blocks were concatenated and the total normalised,
    so a requested fraction c arrived as c/(1+c): the more contamination asked
    for, the further short the result fell. Small fractions hid it (8.6%
    arrived as 7.9%); 70% would have arrived as 41%.
    """

    @staticmethod
    def _combine(viral_weights, profile):
        import numpy as np
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        from viroforge.generator import FASTQGenerator

        gen = FASTQGenerator.__new__(FASTQGenerator)  # pure method, no setup
        viral = [SeqRecord(Seq("ACGT" * 50), id=f"v{i}")
                 for i in range(len(viral_weights))]
        _, ab = gen._combine_viral_and_contamination(
            viral, np.array(viral_weights, dtype=float), profile)
        ab = np.asarray(ab)
        n = len(viral)
        return ab[:n], ab[n:]

    def _bg_profile(self, pct):
        p = ContaminationProfile(name="t")
        add_bacterial_background(p, abundance_pct=pct, n_fragments=10,
                                 fragment_length=400, random_seed=1)
        return p

    def test_requested_fraction_is_delivered(self):
        for pct in (5.0, 25.0, 70.0, 90.0):
            _, contam = self._combine([0.5, 0.3, 0.2], self._bg_profile(pct))
            assert contam.sum() == pytest.approx(pct / 100.0), pct

    def test_old_behaviour_would_have_undershot(self):
        """Guards the specific arithmetic that was wrong."""
        _, contam = self._combine([0.5, 0.3, 0.2], self._bg_profile(70.0))
        naive = 0.70 / 1.70  # concatenate-then-normalise
        assert contam.sum() == pytest.approx(0.70)
        assert abs(contam.sum() - naive) > 0.25

    def test_total_is_one(self):
        viral, contam = self._combine([0.5, 0.3, 0.2], self._bg_profile(70.0))
        assert viral.sum() + contam.sum() == pytest.approx(1.0)

    def test_viral_structure_is_preserved(self):
        """Scaling must not distort relative abundance within the community."""
        weights = [0.5, 0.3, 0.15, 0.05]
        viral, _ = self._combine(weights, self._bg_profile(70.0))
        expected = [w / sum(weights) for w in weights]
        actual = list(viral / viral.sum())
        assert actual == pytest.approx(expected)

    def test_viral_fills_the_remainder(self):
        viral, _ = self._combine([1.0, 1.0], self._bg_profile(70.0))
        assert viral.sum() == pytest.approx(0.30)

    def test_no_contamination_leaves_viral_at_one(self):
        viral, contam = self._combine([0.6, 0.4], ContaminationProfile(name="t"))
        assert contam.size == 0
        assert viral.sum() == pytest.approx(1.0)

    def test_input_scale_does_not_matter(self):
        """VLP-reduced viral abundances must give the same answer."""
        a, _ = self._combine([0.5, 0.3, 0.2], self._bg_profile(70.0))
        b, _ = self._combine([5.0, 3.0, 2.0], self._bg_profile(70.0))
        c, _ = self._combine([0.05, 0.03, 0.02], self._bg_profile(70.0))
        assert list(b) == pytest.approx(list(a))
        assert list(c) == pytest.approx(list(a))

    def test_over_requested_contamination_is_clamped(self):
        """A profile summing past 1 must not erase the viral community."""
        viral, contam = self._combine([0.5, 0.5], self._bg_profile(140.0))
        assert contam.sum() < 1.0
        assert viral.sum() > 0
        assert viral.sum() + contam.sum() == pytest.approx(1.0)

    def test_zero_viral_abundance_is_an_error(self):
        with pytest.raises(ValueError, match="Viral abundances"):
            self._combine([0.0, 0.0], self._bg_profile(50.0))


class TestBundledReferenceAndFallback:
    """Since bacterial_fragments.fasta is bundled, the real path is the default
    and the synthetic path only runs when a reference cannot be resolved. Both
    need coverage, or the fallback rots unnoticed."""

    def test_bundled_reference_is_used_by_default(self):
        p = _profile(community_type="soil", n_fragments=15)
        assert not any(c.source == "synthetic" for c in p.contaminants)

    def test_bundled_reference_is_discoverable(self):
        from viroforge.data.references.resolver import get_bacterial_fragments_path

        path = get_bacterial_fragments_path()
        assert path is not None and path.exists()

    def test_synthetic_fallback_still_works(self):
        p = ContaminationProfile(name="t")
        add_bacterial_background(p, abundance_pct=50.0, n_fragments=10,
                                 fragment_length=500, random_seed=42,
                                 use_real_references=False)
        assert len(p.contaminants) == 10
        assert all(c.source == "synthetic" for c in p.contaminants)
        assert sum(c.abundance for c in p.contaminants) == pytest.approx(0.50)

    def test_modelled_gc_matches_the_bundled_reference(self):
        """gc_content shapes the synthetic fallback, so it must track the real
        set. Otherwise a run without the reference would misrepresent the
        community relative to a run with it."""
        import collections
        import re

        from Bio import SeqIO

        from viroforge.data.references.resolver import get_bacterial_fragments_path

        by = collections.defaultdict(list)
        for rec in SeqIO.parse(get_bacterial_fragments_path(), "fasta"):
            community = re.search(r"\[(\w+)\]", rec.description).group(1)
            seq = str(rec.seq).upper()
            by[community].append((seq.count("G") + seq.count("C")) / len(seq) * 100)

        for community, values in by.items():
            actual = sum(values) / len(values)
            modelled = BACTERIAL_COMMUNITY_PROFILES[community]["gc_content"]
            assert actual == pytest.approx(modelled, abs=1.0), (
                f"{community}: reference is {actual:.1f}% GC but the model says "
                f"{modelled:.1f}%. Remeasure after changing the taxa."
            )

    def test_every_community_has_reference_fragments(self):
        import re

        from Bio import SeqIO

        from viroforge.data.references.resolver import get_bacterial_fragments_path

        present = {re.search(r"\[(\w+)\]", r.description).group(1)
                   for r in SeqIO.parse(get_bacterial_fragments_path(), "fasta")}
        assert present == set(BACTERIAL_COMMUNITY_PROFILES)


class TestPerCollectionBaselines:
    """Each collection carries its own pre-VLP microbiome level.

    The baseline is what was in the SAMPLE, before any VLP step, so one number
    serves both workflows: VLP removes ~98% of it and leaves a realistic
    residual, while --no-vlp keeps the full bulk-metagenome background. It is
    therefore not scaled by --contamination-level, which models prep quality.
    """

    DB = "viroforge/data/viral_genomes.db"

    @staticmethod
    def _resolve(flag, meta):
        """Mirror of the resolution in run_generation()."""
        if flag is None:
            base = meta.get("default_bacterial_pct")
            return (base / 100.0) if base else 0.0, "collection baseline"
        return flag, "flag"

    def _meta(self, cid):
        from viroforge.core.collection import CollectionLoader

        meta, _ = CollectionLoader(self.DB).load_collection(cid)
        return meta

    def test_columns_reach_the_generator(self):
        meta = self._meta(1)
        assert meta.get("default_bacterial_pct") is not None
        assert meta.get("bacterial_community") is not None

    def test_baseline_used_when_flag_absent(self):
        frac, src = self._resolve(None, self._meta(1))
        assert src == "collection baseline"
        assert frac == pytest.approx(0.70)

    def test_explicit_flag_overrides_baseline(self):
        frac, src = self._resolve(0.25, self._meta(1))
        assert (frac, src) == (0.25, "flag")

    def test_explicit_zero_switches_it_off(self):
        """0.0 must be distinguishable from 'not given', or there is no way to
        disable a nonzero baseline."""
        frac, src = self._resolve(0.0, self._meta(1))
        assert frac == 0.0 and src == "flag"

    def test_parser_default_is_none_not_zero(self):
        from viroforge.generator import build_parser

        assert build_parser().parse_args([]).bacterial_fraction is None
        assert build_parser().parse_args(
            ["--bacterial-fraction", "0"]).bacterial_fraction == 0.0

    def test_biomass_ordering_across_sites(self):
        """Soil and gut are bacteria-dominated; blood is near-sterile."""
        pct = {cid: self._meta(cid).get("default_bacterial_pct")
               for cid in (1, 5, 6, 17, 18)}
        assert pct[6] > pct[5] > pct[1] > pct[18] > pct[17], pct
        assert pct[17] <= 5.0, "blood should be near-sterile"

    def test_every_collection_has_a_baseline(self):
        import sqlite3

        con = sqlite3.connect(self.DB)
        rows = con.execute(
            "SELECT collection_id, default_bacterial_pct, bacterial_community "
            "FROM body_site_collections ORDER BY collection_id").fetchall()
        con.close()
        assert len(rows) == 20
        for cid, pct, community in rows:
            assert pct is not None, f"collection {cid} has no bacterial baseline"
            assert community in BACTERIAL_COMMUNITY_PROFILES, (cid, community)

    def test_rna_collections_are_zero(self):
        """Bacterial background is wired into the DNA path only; an RNA
        library's bacterial signal is ribosomal and covered by rrna_pct."""
        for cid in (13, 14, 15):
            assert self._meta(cid).get("default_bacterial_pct") == 0.0


class TestContaminationLevelDoesNotScaleBacterial:
    """--contamination-level must leave bacterial background alone.

    The level models how well the VLP prep went; the bacterial baseline models
    how much microbiome was in the sample. They are independent, so a 'heavy'
    prep does not conjure extra bacteria into the tube. Host, rRNA and reagent
    all scale; this one must not.
    """

    @staticmethod
    def _bacterial_total(level):
        from viroforge.core.contamination import create_contamination_profile

        profile = create_contamination_profile(
            level,
            random_seed=42,
            use_real_references=False,
            collection_defaults={
                "default_host_pct": 5.0,
                "default_rrna_pct": 3.0,
                "default_reagent_pct": 0.5,
                "default_phix_pct": 0.1,
                "host_organism": "human",
            },
            bacterial_pct=70.0,
        )
        by_type = {}
        for c in profile.contaminants:
            by_type[c.contaminant_type] = by_type.get(c.contaminant_type, 0.0) + c.abundance
        return by_type

    def test_bacterial_is_identical_across_levels(self):
        clean = self._bacterial_total("clean")[ContaminantType.BACTERIAL_BACKGROUND]
        heavy = self._bacterial_total("heavy")[ContaminantType.BACTERIAL_BACKGROUND]
        assert clean == pytest.approx(0.70)
        assert heavy == pytest.approx(0.70)

    def test_host_still_scales_with_level(self):
        """Guards against 'fixing' this by exempting everything."""
        clean = self._bacterial_total("clean")[ContaminantType.HOST_DNA]
        heavy = self._bacterial_total("heavy")[ContaminantType.HOST_DNA]
        assert heavy > clean * 5, (clean, heavy)
