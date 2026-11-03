"""Numerical checks against a full OLS solve, independent of residual projection."""
import unittest

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from api.association import AssociationDataError, CovariateModel


class NominalTest(unittest.TestCase):
    def test_matches_full_ols_with_covariates_and_missing_hardcalls(self):
        rng = np.random.default_rng(7301)
        n, k = 83, 5
        covariates = rng.normal(size=(n, k))
        g = rng.integers(0, 3, size=n).astype(float)
        g[[3, 13, 70]] = np.nan
        imputed = np.where(np.isnan(g), np.nanmean(g), g)
        y = .7 * imputed + covariates @ rng.normal(size=k) + rng.normal(size=n)
        model = CovariateModel.from_frame(pd.DataFrame(covariates))
        actual = model.test([y, -y], g)
        design = np.column_stack([np.ones(n), covariates, imputed])
        coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
        errors = y - design @ coefficients
        se = np.sqrt((errors @ errors) / (n - k - 2) * np.linalg.inv(design.T @ design)[-1, -1])
        expected_p = 2 * student_t.sf(abs(coefficients[-1] / se), n - k - 2)
        self.assertAlmostEqual(actual[0]["slope"], coefficients[-1], places=12)
        self.assertAlmostEqual(actual[0]["slopeSE"], se, places=12)
        np.testing.assert_allclose(actual[0]["pValue"], expected_p, rtol=1e-11)
        self.assertEqual(actual[0]["pValue"], actual[1]["pValue"])
        self.assertEqual(actual[0]["slope"], -actual[1]["slope"])
        self.assertEqual(actual[0]["imputedGenotypeCount"], 3)

    def test_untestable_vectors_are_explicit(self):
        model = CovariateModel.from_frame(pd.DataFrame(index=range(8)))
        for y, g, reason in [
            (np.arange(8), [np.nan] * 8, "no_called_genotypes"),
            (np.arange(8), [1] * 8, "no_residual_genotype_variation"),
            ([1] * 8, np.arange(8) % 3, "no_residual_phenotype_variation"),
            ([1, 2, 3, 4, 5, 6, 7, np.nan], np.arange(8) % 3, "nonfinite_phenotype"),
        ]:
            with self.subTest(reason=reason):
                result = model.test(y, g)[0]
                self.assertEqual(result["associationStatus"], reason)
                self.assertIsNone(result["pValue"])
        with self.assertRaises(AssociationDataError):
            model.test(np.arange(8), [3] * 8)

    def test_rank_deficiency_and_insufficient_dof_rejected(self):
        for frame in [pd.DataFrame({"constant": np.ones(10)}),
                      pd.DataFrame({"x": range(10), "twice": np.arange(10) * 2}),
                      pd.DataFrame({"x": range(3)})]:
            with self.assertRaises(AssociationDataError):
                CovariateModel.from_frame(frame)


if __name__ == "__main__":
    unittest.main()
