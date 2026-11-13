reference_group_average <- function(probability, group) {
  onehot <- unname(stats::model.matrix(~group-1))
  numerator <- Matrix::crossprod(onehot, probability %*% onehot)
  count <- probability
  count@x[] <- 1
  denominator <- Matrix::crossprod(onehot, count %*% onehot)
  result <- as.matrix(numerator / denominator)
  result[is.nan(result)] <- 0
  unname(result)
}

reference_aggregation_fixture <- function(complex_ligand = FALSE, two_cells = FALSE) {
  L <- if (two_cells) c(.2,.1) else c(.2,.1,.6,.9,.3,.4,.7)
  R <- if (two_cells) c(.5,1) else c(.5,1,.3,.2,.7,.8,.4)
  expression <- Matrix::Matrix(rbind(L=L,R=R), sparse=TRUE,
    dimnames=list(c("L","R"),paste0("c",seq_along(L))))
  coordinates <- cbind(x=seq_along(L),y=0)
  rownames(coordinates) <- colnames(expression)
  group <- if (two_cells) factor(c("A","B")) else factor(c("A","A","B","B","C","C","C"))
  complex <- if (complex_ligand) data.frame(subunit_1="L",row.names="LC") else
    data.frame(subunit_1=character(),row.names=character())
  cofactor <- data.frame(cofactor1=character(),row.names=character())
  lr <- data.frame(ligand=if(complex_ligand) "LC" else "L",receptor="R",
    interaction_name="test",annotation="Secreted Signaling")
  list(expression=expression,coordinates=coordinates,group=group,lr=lr,
       complex=complex,cofactor=cofactor)
}

run_reference_fixture <- function(fixture, ...) {
  do.call(spatialess,c(fixture,list(tol=0,interaction.range=10,
    scale.distance=1,contact.range=1,min.percent=0,min.cells.sr=1L,
    use.AGAN=FALSE),list(...)))
}

fixture_cell_probability <- function(fixture) {
  x <- as.matrix(fixture$expression)/max(fixture$expression)
  ligand <- x["L",]
  if (fixture$lr$ligand == "LC") ligand <- exp(log(ligand))
  receptor <- x["R",]
  product <- outer(ligand,receptor)
  distance <- as.matrix(stats::dist(fixture$coordinates))
  distance[distance==0] <- min(distance[distance>0])
  Matrix::Matrix(product/(.5+product)*(1/distance),sparse=TRUE)
}

test_that("native reference aggregation uses the literal Matrix reduction order", {
  set.seed(1318)
  cells <- 127L
  group <- factor(sample(rep(LETTERS[1:5],length.out=cells)))
  probability <- Matrix::sparseMatrix(i=sample.int(cells,1400,TRUE),
    j=sample.int(cells,1400,TRUE),x=10^stats::runif(1400,-20,-2),
    dims=c(cells,cells))
  got <- spatialess_reference_aggregate_cpp(probability,as.integer(group),nlevels(group))
  expected <- reference_group_average(probability,group)
  expect_identical(unname(got$probability),expected)
  expect_identical(unname(got$sum/got$count),expected)
  expect_equal(sum(got$count),length(probability@x),tolerance=0)
})

test_that("default observed and permutation arithmetic matches the Matrix reference", {
  fixture <- reference_aggregation_fixture()
  seed <- 9L
  nboot <- 20L
  got <- run_reference_fixture(fixture,nboot=nboot,seed.use=seed)
  probability <- fixture_cell_probability(fixture)
  observed <- reference_group_average(probability,fixture$group)
  at <- cbind(got$records$sender_group,got$records$receiver_group)
  expect_identical(as.double(got$records$probability),as.double(observed[at]))
  expect_identical(got$parameters$aggregation,"reference")
  expect_identical(got$diagnostics$aggregation,"reference")

  # Reuse exactly the primary API's global RNG setup, rather than a new null.
  set.seed(seed)
  if (requireNamespace("future.apply",quietly=TRUE)) {
    invisible(future.apply::future_lapply(1L,identity,future.seed=TRUE))
  }
  permutations <- replicate(nboot,sample.int(ncol(fixture$expression),ncol(fixture$expression)))
  reject <- integer(nrow(got$records))
  for (b in seq_len(nboot)) {
    null <- reference_group_average(probability,fixture$group[permutations[,b]])
    reject <- reject + (null[at] > observed[at])
  }
  expect_identical(as.double(got$records$pvalue),as.double(reject/nboot))
})

test_that("reference memory guards never silently select native stream sums", {
  fixture <- reference_aggregation_fixture()
  expect_error(run_reference_fixture(fixture,max_reference_entries=1),
               "Reference aggregation exceeds")
  for (limit in list(0,NA_real_,Inf,"large")) {
    expect_error(run_reference_fixture(fixture,max_reference_entries=limit),
                 "max_reference_entries")
  }
  streamed <- run_reference_fixture(fixture,aggregation="stream",max_reference_entries=1)
  expect_identical(streamed$parameters$aggregation,"stream")
  expect_equal(streamed$diagnostics$reference_cell_group_entries,0,tolerance=0)
  probability <- fixture_cell_probability(fixture)
  expect_error(spatialess_reference_aggregate_cpp(probability,
    as.integer(fixture$group),nlevels(fixture$group),1),"max_reference_entries")
})

test_that("the primary stream preserves direct genes versus named one-subunit complexes", {
  for (complex_ligand in c(FALSE,TRUE)) {
    fixture <- reference_aggregation_fixture(complex_ligand,two_cells=TRUE)
    got <- run_reference_fixture(fixture,nboot=0)
    expected <- reference_group_average(fixture_cell_probability(fixture),fixture$group)
    at <- cbind(got$records$sender_group,got$records$receiver_group)
    expect_identical(as.double(got$records$probability),as.double(expected[at]))
  }
})
