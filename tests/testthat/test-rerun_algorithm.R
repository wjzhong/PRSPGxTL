test_that("Test the rerun_algorithm function", {
  data(sum_stats)
  data(ped)
  data(G)
  data(bim)
  data(parameter_tuning_output)
  data(parameter_tuning_fixG_output)
  data(rerun_algorithm_output)
  data(rerun_algorithm_fixG_output)
  # when fixG = FALSE
  result_list_no_fixG = rerun_algorithm(sum_stats = sum_stats, ped = ped, G = G, bim = bim, initial = "PRS", num_snp = 2349, s2_results = parameter_tuning_output, fixG = FALSE)
  expect_equal(result_list_no_fixG, rerun_algorithm_output)
  # when fixG = TRUE
  result_list_fixG = rerun_algorithm(sum_stats = sum_stats, ped = ped, G = G, bim = bim, initial = "PRS", num_snp = 2350, s2_results = parameter_tuning_fixG_output, fixG = TRUE)
  expect_equal(result_list_fixG, rerun_algorithm_fixG_output)

})
