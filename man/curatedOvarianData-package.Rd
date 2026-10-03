\name{curatedOvarianData-package}
\alias{curatedOvarianData-package}
\docType{package}
\title{Clinically Annotated Data for the Ovarian Cancer Transcriptome}
\description{The curatedOvarianData package provides manually curated
  clinical data, uniformly processed expression data, and convenience
  functions for gene expression analysis in patients with ovarian cancer.}
\details{
Datasets are hosted on Zenodo and downloaded on first use with
\code{\link{curatedOvarianData}()}, which caches them locally with
\pkg{BiocFileCache}. Legacy \code{data()} access still works but is
deprecated and downloads through the same mechanism. See the package
vignette for rule-based creation of a list of filtered datasets.

In this package, redundant probe sets are collapsed so that each gene
is represented by the probe set with maximum mean. Historical
alternative versions (NormalizerVcuratedOvarianData, where each gene is
the mean of probesets after removing "noisy" probesets, and
FULLVcuratedOvarianData, with no collapsing but a map in
featureData(eset)) are no longer distributed.

In the "Available sample meta-data" sections of each dataset, please refer to the following key.

For "sample_type": tumor / metastatic / adjacentnormal / healthy /
cellline: "healthy" should be only from individuals without cancer,
"adjacentnormal" from individuals with cancer, "metastatic" for
non-primary tumors.

For "histological_type": ser=serous / endo=endometrioid / clearcell / mucinous,
undifferentiated / other / mix.  Other includes sarcomatoid,
adenocarcinoma, dysgerminoma.

For "primarysite" and for "arrayedsite": ov|ft|other.  ov=ovary;ft=fallopian tube

For "summarygrade": low = 1, 2, LMP.  High= 3,4,23.

For "summarystage": early = 1,2, 12. late=3,4,23,34.

For "tumorstage": FIGO Stage (I-IV, but coded here as 1-4 to ensure
correct ordering in factors).  If multiple stages given (eg 34), use the
highest.

For "substage": substage (abcd).  For cases like ab, bc, use highest given.

For "grade": Grade (1-3): If multiple given, ie 12, 23, use highest
given.  Most ovarian cancer studies use FIGO grading, with a couple
exceptions in this package (Yoshihara and Tothill).

For "pltx": (y/n): patient treated with platin.

For "tax": (y/n): patient treated with taxol.

For "neo": (y/n): patient treated with neoadjuvant treatment.

For "primary_therapy_outcome_success":
completeresponse|partialresponse|progressivedisease|stabledisease:
response to any kind of therapy (including radiation only).

For "days_to_tumor_recurrence": time to recurrence or last follow-up in
days

For "recurrence_status": recurrence censoring variable (recurrence /
norecurrence)

For "days_to_death": time to death or last follow-up in days

For "vital_status": Overall survival censoring variable (living /
deceased)

For "os_binary": dichotomized overall survival variable as defined by
study authors (short / long).

For "relapse_binary": dichotomized relapse variable as defined by study
authors (short / long)

For "site_of_tumor_first_recurrence": (metastasis / locoregional / none
/ locoregional_plus_metastatic). none for no recurrence, na for unknown

For "primary_therapy_outcome_success": (completeresponse /
partialresponse / progressivedisease / stabledisease) Response to any kind
of therapy (including radiation only).

For "debulking": amount of residual disease (optimal = <1cm,
suboptimal=>1cm).

For "percent_normal_cells":
Estimated percentage of normal cells.  An integer 0-100, or can be >70,
<70, etc.

For "percent_stromal_cells":
Estimated percentage of stromal cells.  An integer 0-100, or can be >70,
<70, etc.

For "percent_tumor_cells":
Estimated percentage of tumor cells.  An integer 0-100, or can be >70,
<70, etc.

For "batch": batch variable when available.  Hybridization date when
Affymetrix CEL files are available.

For "uncurated_author_metadata":
Original uncurated data, with each field separated by ///.


}
\author{
Benjamin F. Ganzfried, Steve Skates, Markus Riester,
Victoria Wang, Thomas Risch, Benjamin Haibe-Kains, Curtis
Huttenhower, Svitlana Tyekucheva, Jie Ding, Ina Jazic,
Michael Birrer, Giovanni Parmigiani, Levi Waldron

Department of Biostatistics and Computational Biology, Dana-Farber Cancer Institute, Harvard School of Public Health

Maintainer: Levi Waldron <lwaldron.research@gmail.com>
}

\examples{
##List all datasets:
curatedOvarianData()
##Load one dataset (test=TRUE uses a small offline subset):
eset <- curatedOvarianData("GSE30161_eset", test = TRUE)
##
##See the actual template used for syntax checking of clinical metadata:
template.file <- system.file("extdata/template_ov.csv", package = "curatedOvarianData")
template <- read.csv(template.file, as.is=TRUE)
head(template)
}
