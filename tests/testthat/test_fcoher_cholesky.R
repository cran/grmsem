#checks function grmsem.fcoher
context("check grmsem.fcoher cholesky")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fcoher.R, quadvariate model, post, Vi estimates, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_Vi_estimates_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_Vi_estimates_cholesky.log")
  out.fcoher.large <-
    grmsem.fcoher(
      fit.large
    )
  expect_equal(out.fcoher.large$fcoher.model.out$Vi,
     c(3.48435105464312e-01,
       1.65021638351001e-01,
       1.69977027793781e-01,
       6.30921498951806e-02,
       3.74387056075587e-01,
       1.70151586195351e-01,
       1.65239264090511e-01,
       1.96554664922273e-01,
       2.45514836675710e-04,
       4.30275571603973e-01,
       6.30824326677685e-01,
       8.42445818535660e-02,
       1.98242176668627e-01,
       4.51576225008531e-05,
       2.71521657997237e-01,
       2.98267944443768e-03,
       3.32306200251720e-04,
       1.64503262153405e-01,
       5.18457207737159e-04,
       2.18208637577865e-01
     )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_Vi_estimates_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fcoher.R, quadvariate model, post, Vi se, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_Vi_se_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_Vi_se_cholesky.log")
  out.fcoher.large <-
    grmsem.fcoher(
      fit.large
    )
  expect_equal(out.fcoher.large$fcoher.model.out$Vi.se,
     c(0.022356343128524598,
       0.018143375933269982,
       0.017985361259552038,
       0.012737851260687780,
       0.016008453197127074,
       0.012501712413665340,
       0.016461613560407592,
       0.012013785230538128,
       0.000773995728760167,
       0.019402101014056160,
       0.019878091585601596,
       0.007571274679006054,
       0.010823625553937823,
       0.000157211273211377,
       0.009214386386701651,
       0.001094511165146149,
       0.000436332757909888,
       0.006060648150202995,
       0.000564331332294193,
       0.008415001538399556
     )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_Vi_se_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fcoher.R, quadvariate model, post, FCOHER, estimates, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_FCOHER_estimates_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_FCOHER_estimates_cholesky.log")
  out.fcoher.large <-
    grmsem.fcoher(
      fit.large
    )
  expect_equal(out.fcoher.large$fcoher.model.out$FCOHER,
               c(1.000000000000000000,
                 0.305930623766501930,
                 0.316717577150081031,
                 0.095760659410648138,
                 0.694069376233498070,
                 0.317042831184310370,
                 0.250798568698737956,
                 0.366239591665608599, 
                 0.000372640062103184,
                 0.653068131828510645,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA
               )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_FCOHER_estimates_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fcoher.R, quadvariate model, post, FCOHER, se, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_FCOHER_se_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_FCOHER_se_cholesky.log")
  out.fcoher.large <-
    grmsem.fcoher(
      fit.large
    )
  expect_equal(out.fcoher.large$fcoher.model.out$FCOHER.se,
               c(0.00000000000000000,
                 0.02752669131138025,
                 0.02513593511759606,
                 0.01891418005841962,
                 0.02752669131138025,
                 0.02445917762717525,
                 0.02318344319982057,
                 0.02178638599939254,
                 0.00117446524384851,
                 0.02353685897972738,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA
               )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_FCOHER_se_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fcoher.R, quadvariate model, post, FCOENV, estimates, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_FCOENV_estimates_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_FCOENV_estimates_cholesky.log")
  out.fcoher.large <-
    grmsem.fcoher(
      fit.large
    )
  expect_equal(out.fcoher.large$fcoher.model.out$FCOENV,
               c(NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 1.000000000000000000,
                 0.236797572161134756,
                 0.542047949740050772,
                 0.000206100789448072,
                 0.763202427838865161,
                 0.008155455638946806,
                 0.001516655802884102,
                 0.449796594621002299,
                 0.002366254773657597,
                 0.995910988634010175
               )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_FCOENV_estimates_cholesky.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fcoher.R, quadvariate model, post, FCOENV, se, cholesky
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("checked_quadvariate_post_FCOENV_se_cholesky", {
  # redirect function output to file
  sink(file = "test_quadvariate_post_FCOENV_se_cholesky.log")
  out.fcoher.large <-
    grmsem.fcoher(
      fit.large
    )
  expect_equal(out.fcoher.large$fcoher.model.out$FCOENV.se,
               c(NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 0.000000000000000000,
                 0.017233975006478573,
                 0.016190895955562890,
                 0.000717577856904393,
                 0.017233975006478573,
                 0.002999598633011726,
                 0.001994418395796650,
                 0.016229481928463105,
                 0.002578211634564575,
                 0.003450609958129195
               )
  )
  # end redirect
  sink()
  file.remove("test_quadvariate_post_FCOENV_se_cholesky.log")
})
