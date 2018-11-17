context("coverage")

test_that(".bestNcbCoverageCombination", {
    m1 <- sparseMatrix(i=c(1, 2, 2, 3, 5, 5, 6, 9, 8, 9, 8),
                       j=c(1, 5, 3, 4, 4, 5, 4, 4, 5, 5, 1),
                       x=c(1, 1, 3, 1, 1, 2, 2, 2, 3, 1, 2))
    expect_equal(.bestNcbCoverageCombination(m1),
                 cbind(index=c(5:4, 1, 3), fragments=c(5:4, 1, 1),
                       bonds=c(4, 2, 1, 0)))
    expect_equal(.bestNcbCoverageCombination(m1,
                                                        maximise="fragments"),
                 cbind(index=c(5:4, 1, 3), fragments=c(5:4, 1, 1),
                       bonds=c(4, 2, 1, 0)))
    expect_equal(.bestNcbCoverageCombination(m1, maximise="bonds"),
                 cbind(index=c(4, 1, 3), fragments=c(4, 2, 2),
                       bonds=c(4, 2, 1)))
    expect_equal(.bestNcbCoverageCombination(m1, intensity=1:5),
                 cbind(index=c(5:3, 1), fragments=c(5:4, 1, 1),
                       bonds=c(4, 2, 0, 1)))
    expect_equal(.bestNcbCoverageCombination(m1, minN=3),
                 cbind(index=c(5:4), fragments=c(5, 4), bonds=c(4, 2)))
    expect_equal(.bestNcbCoverageCombination(m1, n=3),
                 cbind(index=c(5:4, 1), fragments=c(5:4, 1), bonds=c(4, 2, 1)))
    expect_equal(.bestNcbCoverageCombination(m1, minN=5),
                 cbind(index=5, fragments=5, bonds=4))
})

test_that(".highestNcbBondCoverage", {
    m1 <- sparseMatrix(i=c(1:10, 2), j=c(1, 4, 3, 4, 4, 5, 4, 5, 5, 5, 1),
                       x=rep(1, 11))
    expect_error(.highestNcbBondCoverage(as.matrix(1:10)))
    expect_error(.highestNcbBondCoverage(m))
    expect_equal(.highestNcbBondCoverage(m1, intensity=1:5),
                 c(index=5, bonds=4))
    expect_equal(.highestNcbBondCoverage(m1, intensity=5:1),
                 c(index=4, bonds=4))
    expect_equal(.highestNcbBondCoverage(t(m1)),
                 c(index=2, bonds=2))
})

test_that(".highestNcbFragmentCoverage", {
    m1 <- sparseMatrix(i=c(1:10, 2), j=c(1, 5, 3, 4, 4, 5, 4, 5, 5, 5, 1),
                       x=c(1, 1, 3, 1, 1, 2, 2, 2, 3, 1, 1))
    expect_error(.highestNcbFragmentCoverage(as.matrix(1:10)))
    expect_error(.highestNcbFragmentCoverage(m))
    expect_equal(.highestNcbFragmentCoverage(m1, intensity=1:5),
                 c(index=5, fragments=6))
    expect_equal(.highestNcbFragmentCoverage(m1[,1:3]),
                 c(index=1, fragments=2))
    expect_equal(.highestNcbFragmentCoverage(m1[,1:3], intensity=1:3),
                 c(index=3, fragments=2))
    expect_equal(.highestNcbFragmentCoverage(t(m1)),
                 c(index=2, fragments=2))
    expect_equal(.highestNcbFragmentCoverage(t(m1)),
                 c(index=2, fragments=2))
})

test_that(".removeNcbCombinations", {
    m1 <- sparseMatrix(i=c(1, 2, 2, 3, 5, 5, 6, 9, 8, 9, 8),
                       j=c(1, 5, 3, 4, 4, 5, 4, 4, 5, 5, 1),
                       x=c(1, 1, 3, 1, 1, 2, 2, 2, 3, 1, 2))
    r1 <- sparseMatrix(i=c(1, 2, 3, 5, 6, 9),
                       j=c(1, 3, 4, 4, 4, 4),
                       x=c(1, 2, 1, 1, 2, 2), dims=c(9, 5))
    r2 <- sparseMatrix(i=c(1, 2, 2, 5, 8, 9, 8),
                       j=c(1, 5, 3, 5, 5, 5, 1),
                       x=c(1, 1, 3, 2, 3, 1, 2))
    expect_error(.removeNcbCombinations(matrix(1:10, ncol=2), 2))
    expect_error(.removeNcbCombinations(m, 2))
    expect_equal(.removeNcbCombinations(m1, 5), r1)
    expect_equal(.removeNcbCombinations(m1, 4), r2)
})
