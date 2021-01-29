#checks function grmsem.var VA miss, across Cholesky, IP and IPC model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.var, trivariate model Cholesky, VA, allow for missingness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("small_trivariate_var_cholesky_VA_miss", {
  # redirect function output to file
  sink(file = "test_small_trivariate_var_cholesky_VA_miss.log")
  out <-
    grmsem.fit(
      ph.small,
      G.small,
      LogL = TRUE,
      estSE = TRUE,
      model = "Cholesky",
      compl.ph =FALSE
    )
  var.out <- grmsem.var(out)
  expect_equal(var.out$VA[lower.tri(var.out$VA, diag = TRUE)],
               c(0.3156890620186718,
                 0.2906359137782761,
                 0.3383613888581842,
                 0.8745761764046308,
                 0.5696366553539253,
                 0.6489524140756368
                 ),tolerance = 1E-4)
  # end redirect
  sink()
  file.remove("test_small_trivariate_var_cholesky_VA_miss.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.var, trivariate model IP, VA, allow for missingness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("small_trivariate_var_IP_VA_miss", {
  # redirect function output to file
  sink(file = "test_small_trivariate_var_IP_VA_miss.log")
  out <-
    grmsem.fit(
      ph.small,
      G.small,
      LogL = TRUE,
      estSE = TRUE,
      model = "IP",
      compl.ph =FALSE
    )
  var.out <- grmsem.var(out)
  expect_equal(var.out$VA[lower.tri(var.out$VA, diag = TRUE)],
               c(0.2858746943535492,
                 0.3058342616744029,
                 0.3296139099272800,
                 0.9663147333229107,
                 0.6364094702321582,
                 0.6858924590898962
                 ),tolerance = 1E-4)
  # end redirect
  sink()
  file.remove("test_small_trivariate_var_IP_VA_miss.log")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: grmsem.var, trivariate model IPC, VA, allow for missingness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("small_trivariate_var_IPC_VA_miss", {
  # redirect function output to file
  sink(file = "test_small_trivariate_var_IPC_VA_miss.log")
  out <-
    grmsem.fit(
      ph.small,
      G.small,
      LogL = TRUE,
      estSE = TRUE,
      model = "IPC",
      compl.ph =FALSE
    )
  var.out <- grmsem.var(out)
  expect_equal(var.out$VA[lower.tri(var.out$VA, diag = TRUE)],
               c(0.3129824227540588,
                 0.2977307532360248,
                 0.3381152634892307,
                 0.8761301418942868,
                 0.5715022623611307,
                 0.6490215026008175
                 ),tolerance = 1E-4)
  # end redirect
  sink()
  file.remove("test_small_trivariate_var_IPC_VA_miss.log")
})
