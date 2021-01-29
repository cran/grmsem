#checks function grmsem.biher
context("check grmsem.biher cholesky")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.biher.R, quadvariate model, post, BIHER, estimate, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_BIHER_estimate_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_BIHER_estimate_cholesky.log")
  out.var.large <-
    grmsem.var(
      fit.large
    )
  out.biher.large <-
    grmsem.biher(
      ph.large,
      out.var.large
    )
  expect_equal(out.biher.large$BIHER[lower.tri(out.biher.large$BIHER, diag=TRUE)],
     c(0.348435105464312,
       0.473774083211919,
       0.387619745102565,
       0.827229948691183,
       0.539408694426588,
       0.680809964894649,
       0.805887225954664,
       0.536683278911405,
       0.793259100311895,
       0.658852500426340
     )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_BIHER_estimate_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.biher.R, quadvariate model, post, BIHER, se, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_BIHER_se_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_BIHER_se_cholesky.log")
  out.var.large <-
    grmsem.var(
      fit.large
    )
  out.biher.large <-
    grmsem.biher(
      ph.large,
      out.var.large
    )
  expect_equal(out.biher.large$BIHER.se[lower.tri(out.biher.large$BIHER.se, diag=TRUE)],
     c(0.0223563431285246,
       0.0328161247855043,
       0.0281493271200213,
       0.0845186244016231,
       0.0201728254788040,
       0.0271184181311658,
       0.0358886195462430,
       0.0204209874026193,
       0.0454838253142663,
       0.0207824907963596
     )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_BIHER_se_cholesky.log")
})
