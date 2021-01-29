#checks function grmsem.fit
context("check grmsem.fit cholesky estimates")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fit, trivariate model Cholesky, estimates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("small_trivariate_fit_cholesky_estimates", {
  # redirect function output to file
  sink(file = "test_small_trivariate_fit_cholesky_estimates.log")
  out <-
    grmsem.fit(
      ph.small,
      G.small,
      LogL = TRUE,
      estSE = FALSE,
      model = "Cholesky",
      compl.ph = TRUE
    )
  expect_equal(out$model.fit$estimates,
               c(0.5618621379116694,
                 0.5172726442442133,
                 0.6022142551121289,
                 0.7791053766476204,
                 0.3313129429567049,
                 0.4201453782260941,
                 0.8266180695166832,
                 0.2344298585212282,
                 0.4124524873920009,
                 0.3862238003986560,
                 0.3410508919960418,
                 0.2511360676160088
                 ),tolerance = 1E-4)
  # end redirect
  sink()
  file.remove("test_small_trivariate_fit_cholesky_estimates.log")
})
#checks function grmsem.fit
context("check grmsem.fit cholesky estimates")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fit, trivariate model Cholesky, se
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("small_trivariate_fit_cholesky_se", {
  # redirect function output to file
  sink(file = "test_small_trivariate_fit_cholesky_se.log")
  out <-
    grmsem.fit(
      ph.small,
      G.small,
      LogL = TRUE,
      estSE = TRUE,
      model = "Cholesky",
      compl.ph = TRUE
    )
  expect_equal(out$model.out$se,
               c(0.18730759858464494,
                 0.25149989646178950,
                 0.20040024800980732,
                 0.13696073278746754,
                 0.14755530556363713,
                 0.07000366546433208,
                 0.11621362028808112,
                 0.12464040728331546,
                 0.12534750576063711,
                 0.12634060929289900,
                 0.10625213311061614,
                 0.05125420979001577
                 ),tolerance = 1E-4)
  # end redirect
  sink()
  file.remove("test_small_trivariate_fit_cholesky_se.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fit, trivariate model Cholesky, estimates, allow for missingness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("small_trivariate_fit_cholesky_estimates_miss", {
  # redirect function output to file
  sink(file = "test_small_trivariate_fit_cholesky_estimates_miss.log")
  out <-
    grmsem.fit(
      ph.small,
      G.small,
      LogL = TRUE,
      estSE = FALSE,
      model = "Cholesky",
      compl.ph =FALSE
    )
  expect_equal(out$model.fit$estimates,
               c(0.5618621379116694,
                 0.5172726442442133,
                 0.6022142551121289,
                 0.7791053766476204,
                 0.3313129429567049,
                 0.4201453782260941,
                 0.8266180695166832,
                 0.2344298585212282,
                 0.4124524873920009,
                 0.3862238003986560,
                 0.3410508919960418,
                 0.2511360676160088
                 ),tolerance = 1E-4)
  # end redirect
  sink()
  file.remove("test_small_trivariate_fit_cholesky_estimates_miss.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fit, trivariate model Cholesky, se, allow for missingness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("small_trivariate_fit_cholesky_se_miss", {
  # redirect function output to file
  sink(file = "test_small_trivariate_fit_cholesky_se_miss.log")
  out <-
    grmsem.fit(
      ph.small,
      G.small,
      LogL = TRUE,
      estSE = TRUE,
      model = "Cholesky",
      compl.ph =FALSE
    )
  expect_equal(out$model.out$se,
               c(0.18730759858464494,
                 0.25149989646178950,
                 0.20040024800980732,
                 0.13696073278746754,
                 0.14755530556363713,
                 0.07000366546433208,
                 0.11621362028808112,
                 0.12464040728331546,
                 0.12534750576063711,
                 0.12634060929289900,
                 0.10625213311061614,
                 0.05125420979001577
                 ),tolerance = 1E-4)
  # end redirect
  sink()
  file.remove("test_small_trivariate_fit_cholesky_se_miss.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fit, trivariate model IP, estimates, allow for missingness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("small_trivariate_fit_IP_estimates_miss", {
  # redirect function output to file
  sink(file = "test_small_trivariate_fit_IP_estimates_miss.log")
  out <-
    grmsem.fit(
      ph.small,
      G.small,
      LogL = TRUE,
      estSE = FALSE,
      model = "IP",
      compl.ph =FALSE
    )
  expect_equal(out$model.fit$estimates,
               c(3.979949090981518e-01,
                 7.684376223992839e-01,
                 8.281862465883754e-01,
                 3.570360579626421e-01,
                 6.130402546360608e-01,
                 -7.194085747449750e-06,
                 6.152942244128992e-01,
                 3.413718305279933e-01,
                 5.773314519588433e-01,
                 5.774061651980960e-01,
                 2.103480438112928e-01,
                 -3.062940401510094e-06
                 ),tolerance = 1E-4)
  # end redirect
  sink()
  file.remove("test_small_trivariate_fit_IP_estimates_miss.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fit, trivariate model IPC, estimates, allow for missingness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("small_trivariate_fit_IPC_estimates_miss", {
  # redirect function output to file
  sink(file = "test_small_trivariate_fit_IPC_estimates_miss.log")
  out <-
    grmsem.fit(
      ph.small,
      G.small,
      LogL = TRUE,
      estSE = FALSE,
      model = "IPC",
      compl.ph =FALSE
    )
  expect_equal(out$model.fit$estimates,
               c(0.419696414830982967,
                 0.709395512363203640,
                 0.805618660396214992,
                 0.369915317514803998,
                 0.610645681990165778,
                 0.000276807199098239,
                 0.828285210532323179,
                 0.230876171746976139,
                 0.412568565221357542,
                 0.387506689465861254,
                 0.341753246738070393,
                 0.250531022348870147
                 ),tolerance = 1E-4)
  # end redirect
  sink()
  file.remove("test_small_trivariate_fit_IPC_estimates_miss.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.fit, trivariate model DS, estimates, allow for missingness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("small_trivariate_fit_DS_estimates_miss", {
  # redirect function output to file
  sink(file = "test_small_trivariate_fit_DS_estimates_miss.log")
  out <-
    grmsem.fit(
      ph.small,
      G.small,
      LogL = TRUE,
      estSE = FALSE,
      model = "DS",
      compl.ph =FALSE
    )
  expect_equal(out$model.fit$estimates,
               c(0.3161060549591701,
                 0.2912541558379716,
                 0.3389444285765710,
                 0.8754990115296768,
                 0.5705238767868135,
                 0.6497818696466298,
                 0.6828572814905260,
                 0.1932361501829878,
                 0.3403783654034593,
                 0.2034536049882292,
                 0.2277078829351579,
                 0.3487737422396363
                 ),tolerance = 1E-4)
  # end redirect
  sink()
  file.remove("test_small_trivariate_fit_DS_estimates_miss.log")
})
