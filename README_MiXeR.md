## MiXer Pipeline Execution Order

1. **Step 1: Univariate Fit 1 (`01_fit1_univariate.sh`)**
   - Runs single-trait architecture estimation for ADHD and iron traits.
   
2. **Step 2: Univariate Fit 2 (`02_fit2_univariate.sh`)**
   - Refines parameters from Fit 1.
   - **Output:** Generates JSON parameter files required for bivariate analysis (`fit2_*.json`).

3. **Step 3: Bivariate Test 2 (`03_bivariate_test2.sh`)**
   - Uses the JSON parameters generated in Step 2 (`--load-params-file`) to estimate genetic overlap and correlation between ADHD and iron traits.
