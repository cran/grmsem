#checks function grmsem.stpar
context("check grmsem.stpar cholesky")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.stpar.R, quadvariate model, post, estimates, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_stpar_estimates_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_stpar_estimates_cholesky.log")
  out.st.large <-
    grmsem.stpar(
      fit.large
    )
  expect_equal(out.st.large$stand.model.out$stand.estimates,
     c(
       -0.59650222859900726,
       -0.42935496633778553,
       -0.43400309802078807,
       -0.26807173533010270,
       -0.64670535345279234,
       -0.43422589120592708,
       -0.43383040727545397,
       0.46670170979284126,
       -0.01672253961698979,
       0.70006232754093711,
       -0.80261141985173468,
       -0.30677292905967118,
       -0.46870085195182920,
       0.00717181200472886,
       0.55074210724505890,
       -0.05749115116651057,
       -0.01945506045804978,
       0.42695781912080505,
       -0.02430075801309616,
       -0.49853922168890219
     )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_stpar_estimates_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.stpar.R, quadvariate model, post, se, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_stpar_se_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_stpar_se_cholesky.log")
  out.st.large <-
    grmsem.stpar(
      fit.large
    )
  expect_equal(out.st.large$stand.model.out$stand.se,
     c(0.01623373289496773,
       0.02045616010410049,
       0.01954225127039980,
       0.02488761869025469,
       0.01173931254606051,
       0.01476161982141781,
       0.01898325881012201,
       0.01321678040019163,
       0.02469467697215491,
       0.01271022118383219,
       0.01206493903615692,
       0.01294071728838618,
       0.01208167418646985,
       0.01169742339426823,
       0.00961757282672022,
       0.01001635298645605,
       0.01196498706965341,
       0.00833052760770650,
       0.01239105245801407,
       0.01009380549231467
       )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_stpar_se_cholesky.log")
})

