#checks function grmsem.var in detail for VA, VE, RG and RE, Cholesky model
context("check grmsem.var cholesky")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.var.R, quadvariate model, post, varA estimates, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_varA_estimates_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_varA_estimates_cholesky.log")
  out.var.large <-
    grmsem.var(
      fit.large
    )
  expect_equal(out.var.large$VA[lower.tri(out.var.large$VA, diag=TRUE)],
     c(0.348435105464312,
       0.239790183207579,
       0.243363850244518,
       0.148268404937456,
       0.539408694426588,
       0.419874650532940,
       0.350760725810872,
       0.536683278911405,
       0.264288496143899,
       0.658852500426340
     )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_varA_estimates_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.var.R, quadvariate model, post, varA se, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_varA_se_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_varA_se_cholesky.log")
  out.var.large <-
    grmsem.var(
      fit.large
    )
  expect_equal(out.var.large$VA.se[lower.tri(out.var.large$VA.se, diag=TRUE)],
     c(0.0223563431285246,
       0.0166091494940618,
       0.0176733221572812,
       0.0151486798167347,
       0.0201728254788040,
       0.0167246910635207,
       0.0156204464284418,
       0.0204209874026193,
       0.0151537521428407,
       0.0207824907963596
     )
  )
  # end redirect
  sink()
  file.create("test_quadvariate_post_varA_se_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.var.R, quadvariate model, post, varE estimates, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_varE_estimates_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_varE_estimates_cholesky.log")
  out.var.large <-
    grmsem.var(
      fit.large
    )
  expect_equal(out.var.large$VE[lower.tri(out.var.large$VE, diag=TRUE)],
     c(0.63082432667768484,
       0.23052880866394745,
       0.35363256017525463,
       -0.00533727709684121,
       0.35576623985080269,
       0.10077361275028011,
       -0.01144931707420637,
       0.36572811826646967,
       -0.01123158690246034,
       0.21910455860835426
     )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_varE_estimates_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.var.R, quadvariate model, post, varE se, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_varE_se_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_varE_se_cholesky.log")
  out.var.large <-
    grmsem.var(
      fit.large
    )
  expect_equal(out.var.large$VE.se[lower.tri(out.var.large$VE.se, diag=TRUE)],
     c(0.01987809158560160,
       0.01221034662027264,
       0.01367955503481543,
       0.00928786746380551,
       0.01205485995085201,
       0.00908511057966976,
       0.00709329268292308,
       0.01241352510175409,
       0.00722055256084113,
       0.00835835565666557
     )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_varE_se_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.var.R, quadvariate model, post, RG estimates, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_RG_estimates_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_RG_estimates_cholesky.log")
  out.var.large <-
    grmsem.var(
      fit.large
    )
  expect_equal(out.var.large$RG[lower.tri(out.var.large$RG, diag=TRUE)],
               c(1.0000000000000000,
                 0.5531099563075157,
                 0.5627766672047456,
                 0.3094521924476351,
                 1.0000000000000000,
                 0.7803719525793050,
                 0.5883799833946660,
                 1.0000000000000002,
                 0.4444522695821256,
                 0.9999999999999999
               )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_RG_estimates_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.var.R, quadvariate model, post, RG se, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_RG_se_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_RG_se_cholesky.log")
  out.var.large <-
    grmsem.var(
      fit.large
    )
  expect_equal(out.var.large$RG.se[lower.tri(out.var.large$RG.se, diag=TRUE)],
               c(0.000000000000000e+00,
                 2.488356157530102e-02,
                 2.233206934683672e-02,
                 3.056074657092669e-02,
                 0.000000000000000e+00,
                 1.557784315695838e-02,
                 1.976462146594987e-02,
                 3.008488916437748e-18,
                 2.206251169793274e-02,
                 4.557190749931257e-18
                 )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_RG_se_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.var.R, quadvariate model, post, RE estimates, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_RE_estimates_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_RE_estimates_cholesky.log")
  out.var.large <-
    grmsem.var(
      fit.large
    )
  expect_equal(out.var.large$RE[lower.tri(out.var.large$RE, diag=TRUE)],
     c(1.0000000000000002,
       0.4866185078283961,
       0.7362390574671048,
       -0.0143562108318341,
       1.0000000000000002,
       0.2793735495990584,
       -0.0410082757556142,
       1.0000000000000000,
       -0.0396767635446130,
       1.0000000000000000
     )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_RE_estimates_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.var.R, quadvariate model, post, RE se, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_RE_se_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_RE_se_cholesky.log")
  out.var.large <-
    grmsem.var(
      fit.large
    )
  expect_equal(out.var.large$RE.se[lower.tri(out.var.large$RE.se, diag=TRUE)],
   c(1.11145194953163e-17,
     1.77078910165045e-02,
     1.09956785037082e-02,
     2.49918960270910e-02,
     0.00000000000000e+00,
     2.16459138874405e-02,
     2.54997036635452e-02,
     5.39777721839092e-18,
     2.55531091683535e-02,
     7.97145087103191e-18
    )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_RE_se_cholesky.log")
})
