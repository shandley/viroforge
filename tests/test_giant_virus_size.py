"""Tests for giant virus particle size estimation."""

import numpy as np
import pytest

from viroforge.enrichment.vlp import VirionSizeEstimator


class TestGiantVirusSizeEstimation:
    """Verify that giant viruses get realistic diameter estimates."""

    def setup_method(self):
        self.estimator = VirionSizeEstimator()
        self.rng = np.random.default_rng(42)

    def test_mimivirus_diameter(self):
        est = self.estimator.estimate_size(
            1181549, 'dsDNA', rng=self.rng,
            genome_name='Acanthamoeba polyphaga mimivirus, complete genome'
        )
        assert 350 < est.estimated_diameter_nm < 550

    def test_pandoravirus_diameter(self):
        est = self.estimator.estimate_size(
            2473870, 'dsDNA', rng=self.rng,
            genome_name='Pandoravirus salinus, complete genome'
        )
        assert 800 < est.estimated_diameter_nm < 1200

    def test_pithovirus_diameter(self):
        est = self.estimator.estimate_size(
            610033, 'dsDNA', rng=self.rng,
            genome_name='Pithovirus sibericum isolate P1084-T, complete genome'
        )
        assert 1200 < est.estimated_diameter_nm < 1800

    def test_cafeteria_roenbergensis_virus_diameter(self):
        est = self.estimator.estimate_size(
            617453, 'dsDNA', rng=self.rng,
            genome_name='Cafeteria roenbergensis virus BV-PW1, complete genome'
        )
        assert 200 < est.estimated_diameter_nm < 400

    def test_normal_phage_unchanged(self):
        """Regular phages should not be affected by the giant virus lookup."""
        est = self.estimator.estimate_size(
            50000, 'dsDNA', rng=self.rng,
            genome_name='Escherichia phage T4, complete genome'
        )
        assert 40 < est.estimated_diameter_nm < 120

    def test_no_genome_name_falls_back_to_formula(self):
        """Without genome_name, even large genomes use the log-linear formula."""
        est = self.estimator.estimate_size(
            1181549, 'dsDNA', rng=self.rng,
            genome_name=None
        )
        assert 100 < est.estimated_diameter_nm < 200

    def test_giant_virus_confidence_is_high(self):
        est = self.estimator.estimate_size(
            2473870, 'dsDNA', rng=self.rng,
            genome_name='Pandoravirus salinus, complete genome'
        )
        assert est.confidence == 'high'

    def test_backward_compatible_without_genome_name(self):
        """Calling without genome_name kwarg should still work."""
        est = self.estimator.estimate_size(50000, 'dsDNA', rng=self.rng)
        assert est.estimated_diameter_nm > 0
