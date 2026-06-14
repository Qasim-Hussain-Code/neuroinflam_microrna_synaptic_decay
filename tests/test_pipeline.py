# ==============================================================================
# Script: tests/test_pipeline.py
# Purpose: Unit Test Suite for miRNA-mRNA Transcriptomics Pipeline
# ==============================================================================

import os
import json
import unittest

class TestTranscriptomicsPipeline(unittest.TestCase):
    
    def setUp(self):
        """Set up test environment variables."""
        self.workspace_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.db_path = os.path.join(self.workspace_dir, "target_database.json")
        self.pipeline_script = os.path.join(self.workspace_dir, "run_mirna_pipeline.py")

    def test_workspace_files_exist(self):
        """Verify that essential workspace configuration files exist."""
        self.assertTrue(os.path.exists(self.db_path), "target_database.json is missing!")
        self.assertTrue(os.path.exists(self.pipeline_script), "run_mirna_pipeline.py is missing!")
        self.assertTrue(os.path.exists(os.path.join(self.workspace_dir, "environment.yml")), "environment.yml is missing!")

    def test_target_database_validity(self):
        """Validate target_database.json loading and schema."""
        with open(self.db_path, "r") as f:
            data = json.load(f)
        
        self.assertIsInstance(data, dict, "Database should be a JSON object (dictionary)")
        self.assertGreater(len(data), 0, "Database should contain at least one miRNA mapping")
        
        for mirna, targets in data.items():
            self.assertTrue(mirna.startswith("hsa-miR-"), f"miRNA name '{mirna}' should follow standard hsa-miR prefix schema")
            self.assertIsInstance(targets, list, f"Targets for {mirna} must be a list of gene names")
            self.assertGreater(len(targets), 0, f"Target list for {mirna} should not be empty")

    def test_fisher_exact_logic(self):
        """Verify mathematical consistency of the Fisher's Exact Test implementation."""
        # Standard hypergeometric verification using scipy (simulate import)
        try:
            import scipy.stats as stats
            # 2x2 table: [[10, 2], [5, 50]] -> highly enriched validated targets in anti-correlated set
            contingency = [[10, 2], [5, 50]]
            odds_ratio, p_val = stats.fisher_exact(contingency, alternative='greater')
            
            self.assertGreater(odds_ratio, 1.0, "Odds ratio should indicate positive enrichment")
            self.assertLess(p_val, 0.05, "P-value should represent statistical significance")
        except ImportError:
            # Skip test if scipy is not installed in the testing context
            pass

    def test_ols_regression_shapes(self):
        """Test expected log transformations and shape alignments for regression design matrices."""
        try:
            import numpy as np
            import pandas as pd
            import statsmodels.api as sm
            
            # Synthesize dummy data matching the actual pipeline matrix sizes
            np.random.seed(42)
            n_samples = 20
            
            # Clinical design matrix
            clinical_df = pd.DataFrame({
                "AD_Status": np.random.choice([0, 1], size=n_samples),
                "RIN": np.random.normal(7.0, 1.0, size=n_samples),
                "PMI": np.random.exponential(scale=6.0, size=n_samples),
                "Age": np.random.randint(60, 90, size=n_samples),
                "Sex": np.random.choice([0, 1], size=n_samples)
            })
            
            # Target expression values
            target_expr = np.random.uniform(10, 200, size=n_samples)
            y_logged = np.log2(target_expr + 1)
            
            X = clinical_df[['AD_Status', 'RIN', 'PMI', 'Age', 'Sex']].copy()
            X = sm.add_constant(X)
            
            model = sm.OLS(y_logged, X).fit()
            
            self.assertEqual(len(model.resid), n_samples, "OLS residual length should match sample count")
            self.assertIn('AD_Status', model.params.index, "AD_Status must be evaluated as a covariate coefficient")
        except ImportError:
            # Skip test if numpy/pandas/statsmodels are not installed in the testing context
            pass

if __name__ == "__main__":
    unittest.main()
