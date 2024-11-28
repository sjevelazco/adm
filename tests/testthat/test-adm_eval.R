# Creating some dummy data for the tests
set.seed(96) # set the random seed for reproducibility
obs <- rnorm(100) # Observed values
pred <- rnorm(100) # Predicted values

# Test1: Check if the function provides a tibble as output
test_that("Function output is a tibble", {
  expect_true(class(adm_eval(obs, pred))[1] == "tbl_df")
})

# Test 2: Check if the function provides a tibble with 6 variables as output
test_that("Function output has 6 variables", {
  expect_equal(ncol(adm_eval(obs, pred)), 6)
})

# Test3: Check if the output variables are numeric
test_that("Function output variables are numeric", {
  result_vars <- sapply(adm_eval(obs, pred), class)
  names(result_vars) <- NULL
  expect_identical(result_vars, rep("numeric", 6))
})

# Test 4: if function works with integer values as well
obs_int <- as.integer(obs)
pred_int <- as.integer(pred)
test_that("Function works with integer values", {
  expect_true(class(adm_eval(obs, pred))[1] == "tbl_df")
})

# Test 6: Check if function throws an error when obs and pred are non-numeric
obs_char <- as.character(obs)
pred_char <- as.character(pred)
test_that("Function throws an error when obs and pred are non-numeric", {
  expect_error(adm_eval(obs_char, pred_char))
})
