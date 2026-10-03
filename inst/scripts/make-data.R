## Prepare the per-dataset .rda files uploaded to Zenodo, along with the
## download manifest, the data/ stub files, and data/datalist.
##
## The datasets themselves were produced by the curatedOvarianData curation
## pipeline (see the package vignette and inst/extdata/template_ov.csv);
## this script re-serializes the .rda files that shipped inside the package
## through version 1.51.x (git tag: pre-zenodo-refactor) with xz compression
## for hosting on Zenodo.
##
## Run from the package root with the full data/*.rda files present:
##   Rscript inst/scripts/make-data.R            # before upload (placeholder URLs)
##   Rscript inst/scripts/make-data.R <RECORD>   # after publishing, to finalize URLs
## Already-compressed files in the output directory are left untouched, so
## re-running with the record ID only rewrites the manifest. Do NOT
## regenerate the .rda files after uploading (checksums must match Zenodo).
##
## NOTE: this script is intentionally duplicated in curatedCRCData.

pkg <- unname(read.dcf("DESCRIPTION")[, "Package"])
record <- commandArgs(trailingOnly = TRUE)
record <- if (length(record)) record[1L] else "RECORD_ID"

outdir <- file.path("..", paste0(pkg, "-zenodo-upload"))
stubdir <- file.path(outdir, "stubs")
dir.create(stubdir, recursive = TRUE, showWarnings = FALSE)

rdas <- sort(list.files("data", pattern = "\\.rda$"))
stopifnot(length(rdas) > 0L)

rows <- lapply(rdas, function(f) {
    nm <- sub("\\.rda$", "", f)
    out <- file.path(outdir, f)
    if (!file.exists(out)) {
        env <- new.env(parent = emptyenv())
        obj <- load(file.path("data", f), envir = env)
        stopifnot(identical(obj, nm))
        save(list = nm, envir = env, file = out, compress = "xz")
        message("compressed ", f)
    }
    writeLines(sprintf(
        'delayedAssign("%s",\n    %s:::.stubLoad("%s"),\n    assign.env = environment())',
        nm, pkg, nm),
        file.path(stubdir, paste0(nm, ".R")))
    data.frame(dataset = nm, filename = f,
        url = sprintf("https://zenodo.org/records/%s/files/%s?download=1",
                      record, f),
        md5 = unname(tools::md5sum(out)), size_bytes = file.size(out))
})
manifest <- do.call(rbind, rows)

write.csv(manifest, file.path(outdir, "zenodo-manifest.csv"),
          row.names = FALSE, quote = FALSE)
writeLines(manifest$dataset, file.path(stubdir, "datalist"))
writeLines(sprintf("%s  %12d  %s", manifest$md5, manifest$size_bytes,
                   manifest$filename),
           file.path(outdir, "upload-manifest.txt"))

message(nrow(manifest), " datasets; total upload size ",
        round(sum(manifest$size_bytes) / 2^20), " MB; written to ", outdir)
if (identical(record, "RECORD_ID"))
    message("Upload the .rda files to a new Zenodo record, publish it, ",
            "verify the md5s shown by Zenodo against upload-manifest.txt, ",
            "then re-run this script with the record ID and copy ",
            "zenodo-manifest.csv to inst/extdata/.")
