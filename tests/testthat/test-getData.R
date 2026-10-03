# The 30 datasets this package has provided since its initial release; this
# list guards the Zenodo manifest against accidental row loss.
expected_datasets <- c(
    "E.MTAB.386_eset", "GSE12418_eset", "GSE12470_eset", "GSE13876_eset",
    "GSE14764_eset", "GSE17260_eset", "GSE18520_eset",
    "GSE19829.GPL570_eset", "GSE19829.GPL8300_eset", "GSE20565_eset",
    "GSE2109_eset", "GSE26193_eset", "GSE26712_eset", "GSE30009_eset",
    "GSE30161_eset", "GSE32062.GPL6480_eset", "GSE32063_eset",
    "GSE44104_eset", "GSE49997_eset", "GSE51088_eset", "GSE6008_eset",
    "GSE6822_eset", "GSE8842_eset", "GSE9891_eset", "PMID15897565_eset",
    "PMID17290060_eset", "PMID19318476_eset", "TCGA.mirna.8x15kv2_eset",
    "TCGA_eset", "TCGA.RNASeqV2_eset"
)

manifest <- read.csv(system.file("extdata", "zenodo-manifest.csv",
                                 package = "curatedOvarianData"))

test_that("manifest is complete and well-formed", {
    expect_setequal(manifest$dataset, expected_datasets)
    expect_false(anyDuplicated(manifest$dataset) > 0)
    expect_identical(manifest$filename, paste0(manifest$dataset, ".rda"))
    expect_true(all(grepl("^[0-9a-f]{32}$", manifest$md5)))
    expect_true(all(grepl(
        "^https://zenodo\\.org/records/[0-9]+/files/.+\\?download=1$",
        manifest$url)))
    expect_true(all(manifest$size_bytes > 0))
})

test_that("manifest, fixtures, and data/ stubs are in lock-step", {
    fixtures <- sub("\\.rda$", "", list.files(
        system.file("extdata", "testdata", package = "curatedOvarianData"),
        pattern = "\\.rda$"))
    expect_setequal(fixtures, expected_datasets)
    stubs <- data(package = "curatedOvarianData")$results[, "Item"]
    expect_setequal(stubs, expected_datasets)
})

test_that("no-argument call lists the datasets", {
    expect_setequal(curatedOvarianData(), expected_datasets)
})

test_that("getter returns ExpressionSets from offline fixtures", {
    eset <- curatedOvarianData("GSE30161_eset", test = TRUE)
    expect_s4_class(eset, "ExpressionSet")
    expect_gt(ncol(eset), 0)
    esets <- curatedOvarianData(c("GSE8842_eset", "GSE30161_eset"),
                                test = TRUE)
    expect_type(esets, "list")
    expect_named(esets, c("GSE8842_eset", "GSE30161_eset"))
    expect_s4_class(esets[[1]], "ExpressionSet")
})

test_that("unknown dataset names give an informative error", {
    expect_error(curatedOvarianData("NOT_A_DATASET"), "Unknown dataset")
    expect_error(curatedOvarianData("NOT_A_DATASET", test = TRUE),
                 "Unknown dataset")
})

test_that("data() stubs create a binding", {
    e <- new.env()
    data("GSE30161_eset", package = "curatedOvarianData", envir = e)
    expect_true(exists("GSE30161_eset", envir = e))
})

test_that("full download works (opt-in; set RUN_FULL_DOWNLOAD_TESTS=1)", {
    skip_on_bioc()
    skip_if(!nzchar(Sys.getenv("RUN_FULL_DOWNLOAD_TESTS")))
    skip_if_offline("zenodo.org")
    # smallest dataset in the collection
    eset <- curatedOvarianData("GSE30009_eset")
    expect_s4_class(eset, "ExpressionSet")
    # second call must hit the cache (no download message)
    expect_silent(suppressMessages(
        eset2 <- curatedOvarianData("GSE30009_eset")))
    expect_identical(dim(eset), dim(eset2))
})
