## MiXer Pipeline Execution Order

1. **Step 1: Univariate Fit 1**
   - Runs single-trait architecture estimation for ADHD and iron traits.
   
2. **Step 2: Univariate Fit 2**
   - Refines parameters from Fit 1.
   - **Output:** Generates JSON parameter files required for bivariate analysis (`fit2_*.json`).

3. **Step 3: Bivariate Test 2**
   - Uses the JSON parameters generated in Step 2 (`--load-params-file`) to estimate genetic overlap and correlation between ADHD and iron traits.
