# adm 0.5.0

## Bug fixes
-  `fit_abund_qrf` now works when a partition is used; it was calling the internal result wrapper with an outdated argument order. Without a partition it now also returns `predictors` and `metadata`, as the other `fit_abund_*` functions do.
-  `fit_abund_glm` and `fit_abund_gam` no longer fail when `hold_out_set` is used (`object 'train_set' not found`).
-  `adm_uncertainty` now accepts `...` (as documented), so the extra arguments needed for CNN models (`x`, `y`, `rasters`, `sample_size`, `custom_architecture`) can be passed.

## Testing
-  New tests for `fit_abund_*` (RAF, GLM, GAM, NET, SVM, GBM, XGB, QRF), `adm_uncertainty`, `model_selection`, `data_abund_pdp`, `p_abund_pdp`, `p_abund_bpdp`, `generate_dnn_architecture`, `generate_cnn_architecture`, `res_calculate`, `family_selector`, `croppin_hood`, `cnn_make_samples` and `get_partition_samples`.
-  The DNN/CNN tuning and prediction tests were shortened (smaller data, networks and epochs) and the heaviest ones are skipped on CRAN, reducing the test time from about 16 to 2 minutes.
-  The test-coverage workflow now uploads the `covr` report (`cobertura.xml`) to Codecov explicitly; previously the wrong file was uploaded and coverage was reported as 0.

## Documentation and CRAN preparation
-  Documented previously undocumented arguments (`hold_out_set`, `samples_list`, `weight_decay`, `optimizer`, `loss_function`, `framework`, `train_quantiles`, `eval_quantile`, `pred_quantile`, `sample_prop`, `na.rm`, and the arguments of `get_partition_samples`) and removed the non-existent `hold_out_evaluation` argument from the XGB help pages.
-  `grf` and `quantregForest` were added to `Suggests`.
-  The DNN sections of the "Modeling workflow" vignette are only evaluated when torch is installed; vignette titles now match their index entries.
-  Declared global variables to avoid R CMD check notes; updated the package title, CITATION, and README links.

# adm 0.0.2
-  `fit_abund_qrf` & `tune_abund_qrf` were implemented to construct Quantile Regression Random Forests models 
-  `p_abund_pdp` was improved to depict exactly training range values when projection data are used [#209](https://github.com/sjevelazco/adm/pull/209)
-  `p_abund_pdp` for GLM was fixed when polynomials are used [#209](https://github.com/sjevelazco/adm/pull/209)
-  `p_abund_bpdp` for XBT was fixed when polynomials are used [#209](https://github.com/sjevelazco/adm/pull/209)

