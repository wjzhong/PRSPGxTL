test_that("Test the parameter_tuning function", {
  data(ped)
  data(parameter_tuning_output)
  data(parameter_tuning_fixG_output)
  data(inner_layer_CV_output)
  data(inner_layer_CV_fixG_output)
  # when fixG = FALSE
  result_list_no_fixG = parameter_tuning(X = ped, s1_results = inner_layer_CV_output, fixG = FALSE)
  expect_equal(result_list_no_fixG, parameter_tuning_output)
  # when fixG = TRUE
  result_list_fixG = parameter_tuning(X = ped, s1_results = inner_layer_CV_fixG_output, fixG = TRUE)
  expect_equal(result_list_fixG, parameter_tuning_fixG_output)

})
