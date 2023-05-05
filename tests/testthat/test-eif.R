test_that("EIF function produces correct output", {
  # create dummy data
  qn <- list(upshift = c(1, 2, 3, 4, 5), noshift = c(2, 3, 4, 5, 6))
  hn <- list(noshift = c(1, 1, 1, 1, 1))
  y <- c(5, 6, 7, 8, 9)
  fluc_mod_out <- list(
    qn_shift_star = c(2, 3, 4, 5, 6),
    qn_noshift_star = c(3, 4, 5, 6, 7)
  )

  # call eif function
  eif_out <- eif(y, qn, qn_unscaled = qn, hn, estimator = "tmle", fluc_mod_out = fluc_mod_out)

  # test that psi is equal to 4.0
  expect_equal(eif_out$psi, 4.0)

  # test that variance is equal to 0.5
  expect_equal(eif_out$var, 0.5)

  # test that CI is equal to c(3.8147, 5.1853)
  expect_equal(eif_out$CI, c(2.6, 5.4), tolerance = 0.1)

  # test that p-value is equal to 0.0795892
  expect_equal(eif_out$p_value, 0)
})
