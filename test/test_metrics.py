#!/usr/bin/env python3
"""Tests for the FADVI metrics modules.

Covers :mod:`fadvi._disentanglement_metrics` (quantifying how well the latent
subspaces are separated).

The tests use small synthetic data with a *known* structure, so the expected
behaviour of each metric is unambiguous: a latent subspace built to encode a
label must be predictive of it, and two independent subspaces must show low
cross-correlation.
"""

import numpy as np
import pytest


from fadvi._disentanglement_metrics import (
    compute_cross_correlation,
    compute_orthogonality_score,
    evaluate_subspace_predictability,
)


SEED = 0


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def synthetic_latents():
    """Two informative latent spaces plus an independent one.

    ``z_batch`` encodes ``batch``, ``z_bio`` encodes ``label``, and the two are
    generated independently of one another, mimicking well-disentangled
    subspaces.
    """
    rng = np.random.default_rng(SEED)
    n, n_batch, n_label = 300, 3, 4

    batch = rng.integers(0, n_batch, n)
    label = rng.integers(0, n_label, n)

    # each subspace has a strong signal for its own factor and noise elsewhere
    z_batch = rng.normal(0, 0.3, (n, 8))
    z_batch[:, 0] += batch * 3.0
    z_bio = rng.normal(0, 0.3, (n, 8))
    z_bio[:, 0] += label * 3.0

    return dict(z_batch=z_batch, z_bio=z_bio, batch=batch, label=label)


# ---------------------------------------------------------------------------
# Disentanglement metrics
# ---------------------------------------------------------------------------

def test_predictability_recovers_encoded_factor(synthetic_latents):
    """A subspace that encodes a factor should predict it far above chance."""
    d = synthetic_latents
    res = evaluate_subspace_predictability(d["z_batch"], d["batch"])
    assert "accuracy" in res
    # 3 balanced classes -> chance is ~0.33; the signal is strong so expect >0.8
    assert res["accuracy"] > 0.8, f"expected high accuracy, got {res['accuracy']}"
    assert 0.0 <= res["accuracy"] <= 1.0


def test_predictability_is_near_chance_for_unrelated_factor(synthetic_latents):
    """The batch subspace should NOT predict the independent label factor."""
    d = synthetic_latents
    res = evaluate_subspace_predictability(d["z_batch"], d["label"])
    # 4 classes -> chance ~0.25; allow slack but it must be far below the
    # accuracy achieved on the factor the subspace actually encodes.
    assert res["accuracy"] < 0.6, (
        f"batch subspace unexpectedly predicts labels ({res['accuracy']})"
    )


def test_cross_correlation_low_for_independent_subspaces(synthetic_latents):
    """Independently generated subspaces should have low cross-correlation."""
    d = synthetic_latents
    mean_abs, matrix = compute_cross_correlation(d["z_batch"], d["z_bio"])
    assert 0.0 <= mean_abs <= 1.0
    assert matrix.shape == (d["z_batch"].shape[1], d["z_bio"].shape[1])
    assert mean_abs < 0.3, f"expected low cross-correlation, got {mean_abs}"


def test_cross_correlation_higher_for_duplicated_subspace(synthetic_latents):
    """Self-correlation must exceed correlation with an independent subspace.

    Note the score averages over ALL dimension pairs, and only the matching
    diagonal pairs are strongly correlated, so even a perfect duplicate scores
    well below 1. The meaningful check is the ordering, not an absolute value.
    """
    d = synthetic_latents
    self_corr, _ = compute_cross_correlation(d["z_batch"], d["z_batch"].copy())
    indep_corr, _ = compute_cross_correlation(d["z_batch"], d["z_bio"])
    assert self_corr > indep_corr, (
        f"self-correlation ({self_corr}) should exceed independent ({indep_corr})"
    )


def test_orthogonality_score_bounds_and_ordering(synthetic_latents):
    """Independent subspaces should be scored as more orthogonal than identical ones."""
    d = synthetic_latents
    indep = compute_orthogonality_score(d["z_batch"], d["z_bio"])
    same = compute_orthogonality_score(d["z_batch"], d["z_batch"].copy())

    for res in (indep, same):
        assert "orthogonality_score" in res
        assert 0.0 <= res["orthogonality_score"] <= 1.0

    assert indep["orthogonality_score"] > same["orthogonality_score"], (
        "independent subspaces should be more orthogonal than identical ones"
    )
