context("Screening")

f.test.pval = function(pval) {
  out = pval >= 0 && pval <= 1
  return(out)
}

data('hfdata')
N = 10
rets = hfdata[,1:10]

test_that("Alpha", {
  tmp = alphaScreening(rets, control = list(nCore = 1))
  expect_equal(length(tmp$n), N)
  expect_equal(length(tmp$alpha), N)
  expect_equal(length(tmp$npeer), N)
  expect_equal(length(tmp$lambda), N)
  expect_equal(dim(tmp$dalpha), c(N,N))
  expect_equal(dim(tmp$tstat), c(N,N))
  expect_equal(dim(tmp$pval), c(N,N))

  tmp = alphaScreening(rets, control = list(nCore = 1, hac = TRUE))
  expect_equal(length(tmp$n), N)
  expect_equal(length(tmp$alpha), N)
  expect_equal(length(tmp$npeer), N)
  expect_equal(length(tmp$lambda), N)
  expect_equal(dim(tmp$dalpha), c(N,N))
  expect_equal(dim(tmp$tstat), c(N,N))
  expect_equal(dim(tmp$pval), c(N,N))
})

test_that("Sharpe", {
  tmp = sharpeScreening(rets, control = list(nCore = 1))
  expect_equal(length(tmp$n), N)
  expect_equal(length(tmp$sharpe), N)
  expect_equal(length(tmp$npeer), N)
  expect_equal(length(tmp$lambda), N)
  expect_equal(dim(tmp$dsharpe), c(N,N))
  expect_equal(dim(tmp$tstat), c(N,N))
  expect_equal(dim(tmp$pval), c(N,N))

  tmp = sharpeScreening(rets, control = list(nCore = 1, hac = TRUE))
  expect_equal(length(tmp$n), N)
  expect_equal(length(tmp$sharpe), N)
  expect_equal(length(tmp$npeer), N)
  expect_equal(length(tmp$lambda), N)
  expect_equal(dim(tmp$dsharpe), c(N,N))
  expect_equal(dim(tmp$tstat), c(N,N))
  expect_equal(dim(tmp$pval), c(N,N))
})

test_that("Modified Sharpe", {
  tmp = msharpeScreening(rets, level = 0.9, control = list(nCore = 1))
  expect_equal(length(tmp$n), N)
  expect_equal(length(tmp$msharpe), N)
  expect_equal(length(tmp$npeer), N)
  expect_equal(length(tmp$lambda), N)
  expect_equal(dim(tmp$dmsharpe), c(N,N))
  expect_equal(dim(tmp$tstat), c(N,N))
  expect_equal(dim(tmp$pval), c(N,N))

  tmp = msharpeScreening(rets, level = 0.95, control = list(nCore = 1, hac = TRUE))
  expect_equal(length(tmp$n), N)
  expect_equal(length(tmp$msharpe), N)
  expect_equal(length(tmp$npeer), N)
  expect_equal(length(tmp$lambda), N)
  expect_equal(dim(tmp$dmsharpe), c(N,N))
  expect_equal(dim(tmp$tstat), c(N,N))
  expect_equal(dim(tmp$pval), c(N,N))
})

test_that("Sharpe screenings store raw ratio differences for product tests", {
  pair <- rets[, 1:2]
  sr.diff <- unname(sharpe(pair)[1] - sharpe(pair)[2])
  msr.diff <- unname(msharpe(pair, level = 0.9)[1] -
                       msharpe(pair, level = 0.9)[2])

  sr <- sharpeScreening(pair, control = list(nCore = 1, ttype = 2))
  msr <- msharpeScreening(pair, level = 0.9,
                          control = list(nCore = 1, ttype = 2))
  expect_equal(unname(sr$dsharpe[1, 2]), sr.diff)
  expect_equal(unname(msr$dmsharpe[1, 2]), msr.diff)

  sr.xy <- sharpeScreening(pair[, 1], Y = pair[, 2],
                           control = list(nCore = 1, ttype = 2))
  msr.xy <- msharpeScreening(pair[, 1], Y = pair[, 2], level = 0.9,
                             control = list(nCore = 1, ttype = 2))
  expect_equal(unname(sr.xy$dsharpe[1, 1]), sr.diff)
  expect_equal(unname(msr.xy$dmsharpe[1, 1]), msr.diff)
})
