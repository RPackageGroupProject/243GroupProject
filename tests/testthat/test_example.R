
test_that("there are no NA's in initialization", {
  expect_equal(sum(is.na(initialization(30,30))),0)
  #expect_equal(sum(is.na(initialization(30,30))),1) # this obviously fails (but I did it to check)
})
