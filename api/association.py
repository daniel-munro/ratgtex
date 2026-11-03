"""Small, CPU-only nominal additive tests; no permutations or multiple testing.

Match tensorQTL's centered covariate projection, mean-imputed hardcalls and
two-sided t test. Float64 arithmetic can differ slightly from its float32 GPU
implementation and from rounded release exports.
"""
from dataclasses import dataclass

import numpy as np
from scipy.stats import t as student_t


class AssociationDataError(ValueError):
    """Release inputs cannot support a faithful nominal test."""


@dataclass(frozen=True)
class CovariateModel:
    samples: tuple
    basis: np.ndarray
    dof: int

    @classmethod
    def from_frame(cls, frame):
        if not frame.index.is_unique or not frame.columns.is_unique:
            raise AssociationDataError("Duplicate covariate sample or covariate IDs.")
        try:
            values = frame.to_numpy(dtype=float)
        except (TypeError, ValueError) as error:
            raise AssociationDataError("Covariates must be numeric.") from error
        n, k = values.shape
        if not np.isfinite(values).all() or n <= k + 2:
            raise AssociationDataError("Invalid covariates or insufficient residual degrees of freedom.")
        centered = values - values.mean(axis=0)
        if k and np.linalg.matrix_rank(centered) != k:
            raise AssociationDataError("Centered covariates are rank deficient.")
        basis = np.linalg.qr(centered, mode="reduced")[0] if k else np.empty((n, 0))
        basis.setflags(write=False)
        return cls(tuple(frame.index), basis, n - k - 2)

    def test(self, phenotypes, genotypes):
        """Rows are phenotypes, columns are samples in exactly self.samples order.

        Missing genotypes are mean-imputed for inference only. Plotting retains
        null hardcalls. Invalid/constant vectors produce explicit untestable
        results rather than invented p-values or dropped-sample regressions.
        """
        y = np.atleast_2d(np.asarray(phenotypes, dtype=float))
        g = np.asarray(genotypes, dtype=float)
        if y.shape[1] != len(self.samples) or g.shape != (len(self.samples),):
            raise AssociationDataError("Association input dimensions do not match covariates.")
        if np.isinf(g).any() or np.any(np.isfinite(g) & ~np.isin(g, [0, 1, 2])):
            raise AssociationDataError("Expected biallelic diploid hardcall dosages.")
        empty = {"pValue": None, "slope": None, "slopeSE": None,
                 "associationSampleCount": len(g), "degreesOfFreedom": self.dof,
                 "imputedGenotypeCount": int(np.isnan(g).sum())}
        if np.isnan(g).all():
            return [{**empty, "associationStatus": "no_called_genotypes"} for _ in y]
        g = np.where(np.isnan(g), np.nanmean(g), g)
        g -= g.mean()
        gr = g - self.basis @ (self.basis.T @ g)
        gg = gr @ gr
        if gg <= np.finfo(float).eps * max(1, g @ g) * len(g):
            return [{**empty, "associationStatus": "no_residual_genotype_variation"} for _ in y]
        results = []
        for values in y:
            if not np.isfinite(values).all():
                results.append({**empty, "associationStatus": "nonfinite_phenotype"})
                continue
            centered = values - values.mean()
            residual = centered - self.basis @ (self.basis.T @ centered)
            yy = residual @ residual
            if yy <= np.finfo(float).eps * max(1, centered @ centered) * len(g):
                results.append({**empty, "associationStatus": "no_residual_phenotype_variation"})
                continue
            slope = float(gr @ residual / gg)
            errors = residual - slope * gr
            se = float(np.sqrt((errors @ errors) / (self.dof * gg)))
            pvalue = float(2 * student_t.sf(abs(slope / se), self.dof)) if se else 0.0
            results.append({**empty, "pValue": pvalue, "slope": slope, "slopeSE": se,
                            "associationStatus": "ok"})
        return results
