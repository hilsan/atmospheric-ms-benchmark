"""Tests for the similarity metrics in compare_spectra.py.

Unlike bin_and_scale()'s rounding rule, these are standard, well-defined
similarity formulas (cosine similarity, Tanimoto/Jaccard index, set-overlap
fractions) with no domain-specific ambiguity to resolve, so expected
values are computed directly from the documented formula rather than
needing separate ground truth. The one genuine domain fact confirmed by
the project owner was weighted_dot()'s exponents: mass exponent = 3.0,
intensity exponent = 0.6 (the standard Stein & Scott / NIST MS Search
convention) — the code previously had these swapped (m_exp=0.6, i_exp=3.0),
fixed as part of this test.

Toy spectra used throughout (asymmetric peak counts/overlap, to
distinguish fraction_sim_in_ref from fraction_ref_in_sim):
    spec1 (ref) = {10: 100.0, 20: 50.0, 30: 25.0}
    spec2 (sim) = {10: 80.0, 20: 60.0, 40: 10.0, 50: 5.0}
Overlap: {10, 20}. spec1 has 3 peaks, spec2 has 4.

entropy_similarity() is not tested for its actual entropy computation —
the ms_entropy package isn't available in this environment (or CI); the
real computation happens separately via patch_entropy.py on a different
Python 3.8 interpreter. Only the graceful-degradation path (returns NaN
when the library or inputs are unavailable) is tested here.
"""

import math

import pytest

from src.analysis.compare_spectra import (
    cosine_similarity,
    fraction_ref_in_sim,
    fraction_sim_in_ref,
    entropy_similarity,
    read_spectrum,
    tanimoto_index,
    weighted_dot,
)

SPEC1 = {10: 100.0, 20: 50.0, 30: 25.0}
SPEC2 = {10: 80.0, 20: 60.0, 40: 10.0, 50: 5.0}


class TestCosineSimilarity:
    def test_matches_formula(self):
        assert cosine_similarity(SPEC1, SPEC2) == pytest.approx(954.213404660788)

    def test_identical_spectra_score_1000(self):
        assert cosine_similarity(SPEC1, SPEC1) == pytest.approx(1000.0)

    def test_none_input_returns_nan(self):
        assert math.isnan(cosine_similarity(None, SPEC2))
        assert math.isnan(cosine_similarity(SPEC1, None))

    def test_disjoint_spectra_score_zero(self):
        assert cosine_similarity({10: 1.0}, {20: 1.0}) == pytest.approx(0.0)


class TestWeightedDot:
    def test_matches_formula_with_corrected_exponents(self):
        assert weighted_dot(SPEC1, SPEC2) == pytest.approx(91.96211102177416)

    def test_none_input_returns_nan(self):
        assert math.isnan(weighted_dot(None, SPEC2))
        assert math.isnan(weighted_dot(SPEC1, None))

    def test_no_matched_peaks_returns_nan(self):
        assert math.isnan(weighted_dot({10: 1.0}, {20: 1.0}))


class TestTanimotoIndex:
    def test_matches_formula(self):
        assert tanimoto_index(SPEC1, SPEC2) == pytest.approx(0.4)

    def test_identical_spectra_score_one(self):
        assert tanimoto_index(SPEC1, SPEC1) == pytest.approx(1.0)

    def test_none_input_returns_nan(self):
        assert math.isnan(tanimoto_index(None, SPEC2))

    def test_disjoint_spectra_score_zero(self):
        assert tanimoto_index({10: 1.0}, {20: 1.0}) == pytest.approx(0.0)


class TestFractionOverlap:
    def test_fraction_sim_in_ref_uses_sim_peak_count_as_denominator(self):
        # overlap {10,20} = 2 peaks; sim (SPEC2) has 4 peaks -> 2/4 = 50%
        assert fraction_sim_in_ref(SPEC1, SPEC2) == pytest.approx(50.0)

    def test_fraction_ref_in_sim_uses_ref_peak_count_as_denominator(self):
        # overlap {10,20} = 2 peaks; ref (SPEC1) has 3 peaks -> 2/3 = 66.67%
        assert fraction_ref_in_sim(SPEC1, SPEC2) == pytest.approx(66.66666666666667)

    def test_these_are_not_symmetric(self):
        # Different denominators (sim's peak count vs ref's peak count) mean
        # these two functions genuinely differ for asymmetric spectra —
        # this asymmetry is the whole point of having both.
        assert fraction_sim_in_ref(SPEC1, SPEC2) != fraction_ref_in_sim(SPEC1, SPEC2)

    def test_none_input_returns_nan(self):
        assert math.isnan(fraction_sim_in_ref(None, SPEC2))
        assert math.isnan(fraction_ref_in_sim(SPEC1, None))

    def test_empty_sim_returns_nan_for_sim_in_ref(self):
        assert math.isnan(fraction_sim_in_ref(SPEC1, {}))

    def test_empty_ref_returns_nan_for_ref_in_sim(self):
        assert math.isnan(fraction_ref_in_sim({}, SPEC2))


class TestEntropySimilarityGracefulDegradation:
    """Only the fallback path is tested — real entropy computation needs
    ms_entropy, unavailable in this environment (see module docstring)."""

    def test_none_input_returns_nan(self):
        assert math.isnan(entropy_similarity(None, SPEC2))
        assert math.isnan(entropy_similarity(SPEC1, None))

    def test_returns_nan_when_library_unavailable(self):
        # In this environment _HAS_MS_ENTROPY is False (package not
        # installed), so even valid spectra should yield NaN, not an error.
        result = entropy_similarity(SPEC1, SPEC2)
        assert math.isnan(result)


class TestReadSpectrum:
    """read_spectrum() is only ever called on spectra_all/10pct/top20.csv —
    files already written by bin_scale_intensities.py's bin_and_scale(),
    which already binned m/z to integers. So these tests use integer m/z
    input matching real usage, not synthetic fractional/tie-case m/z that
    the function never actually receives in production."""

    def test_missing_file_returns_none(self, tmp_path):
        assert read_spectrum(str(tmp_path / "does_not_exist.csv")) is None

    def test_reads_mz_intensity_pairs(self, tmp_path):
        f = tmp_path / "spectra_all.csv"
        f.write_text("10,100.0000\n20,50.0000\n30,25.0000\n")

        result = read_spectrum(str(f))

        assert result == {10: 100.0, 20: 50.0, 30: 25.0}

    def test_basepeak_rescales_to_given_max(self, tmp_path):
        f = tmp_path / "spectra_all.csv"
        f.write_text("10,100.0000\n20,50.0000\n")

        result = read_spectrum(str(f), basepeak=999)

        assert result[10] == pytest.approx(999.0)
        assert result[20] == pytest.approx(499.5)

    def test_skips_malformed_lines(self, tmp_path):
        f = tmp_path / "spectra_all.csv"
        f.write_text("10,100.0000\nnot,a,valid,line\n20,50.0000\n")

        result = read_spectrum(str(f))

        assert result == {10: 100.0, 20: 50.0}
