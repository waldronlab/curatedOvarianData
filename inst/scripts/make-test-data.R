## Build the small offline test fixtures in inst/extdata/testdata/: one per
## dataset, subset to the first <=200 features and <=100 samples with full
## phenoData retained, so that examples, tests, vignettes, and R CMD check
## run without network access (see the test= argument of the getter and
## .stubLoad in R/getData.R).
##
## Run from the package root with the full data/*.rda files present
## (git tag: pre-zenodo-refactor):
##   Rscript inst/scripts/make-test-data.R
##
## NOTE: this script is intentionally duplicated in curatedCRCData.

suppressPackageStartupMessages(library(Biobase))

dir.create("inst/extdata/testdata", recursive = TRUE, showWarnings = FALSE)
rdas <- sort(list.files("data", pattern = "\\.rda$"))
stopifnot(length(rdas) > 0L)

for (f in rdas) {
    nm <- sub("\\.rda$", "", f)
    env <- new.env(parent = emptyenv())
    obj <- load(file.path("data", f), envir = env)
    stopifnot(identical(obj, nm))
    eset <- get(nm, envir = env)
    ## keep the first ~200 features plus genes used in the vignettes
    ## (survival forest plots, probe-expansion examples)
    fn <- featureNames(eset)
    keep <- unique(c(intersect(c("CXCL12", "CXCR4"), fn),
                     head(grep("///", fn, value = TRUE), 5L),
                     head(fn, 200L)))
    fx <- eset[keep, seq_len(min(100L, ncol(eset)))]
    assign(nm, fx)
    save(list = nm, file = file.path("inst/extdata/testdata", f),
         compress = "xz")
    message("fixture ", nm, ": ", nrow(fx), " x ", ncol(fx))
}
message("total fixture size: ",
        round(sum(file.size(list.files("inst/extdata/testdata",
                                       full.names = TRUE))) / 2^10), " KB")
