#checks input functions
context("check input functions")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test: compare input data imported in different ways
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("test_input", {
  G.small.fromgz <- grm.input("small.gcta.grm.gz")
  G.small.frombin <- grm.bin.input("small.gcta.grm.bin")
  expect_equal(G.small.fromgz, G.small)
  expect_equal(G.small.fromgz, G.small.frombin, tolerance = 1E-7)
})

