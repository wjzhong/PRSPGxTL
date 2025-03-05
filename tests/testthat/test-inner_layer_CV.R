test_that("Test the inner_layer_CV function", {
  data(sum_stats)
  data(ped)
  data(G)
  data(bim)
  data(inner_layer_CV_output)
  data(inner_layer_CV_fixG_output)
  # when fixG = FALSE
  result_list_no_fixG = inner_layer_CV(sum_stats = sum_stats, ped = ped, G = G, bim = bim, initial = "PRS", num_snp = 2349, covar = NULL, fixG = FALSE, verbose = TRUE)
  expect_equal(result_list_no_fixG, inner_layer_CV_output)
  # when fixG = TRUE
  result_list_fixG = inner_layer_CV(sum_stats = sum_stats, ped = ped, G = G, bim = bim, initial = "PRS", num_snp = 2350, covar = NULL, fixG = TRUE, verbose = TRUE)
  expect_equal(result_list_fixG, inner_layer_CV_fixG_output)

})
