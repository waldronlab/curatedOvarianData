\name{curatedOvarianData-package}
\alias{curatedOvarianData}
\docType{package}
\title{Clinically Annotated Data for the Ovarian Cancer Transcriptome}
\description{The curatedOvarianData package provides manually curated
  clinical data, uniformly processed expression data, and convenience
  functions for gene expression analysis in patients with ovarian cancer.}
\details{
\tabular{ll}{
Package: \tab curatedOvarianData\cr
Type: \tab Package\cr
Version: \tab \Sexpr{packageDescription("curatedOvarianData")$Version}\cr
Date: \tab \Sexpr{packageDescription("curatedOvarianData")$Date}\cr
License: \tab Artistic-2.0\cr
Depends: \tab R (>= 2.10.0), affy\cr
}

Please see http://bcb.dfci.harvard.edu/ovariancancer for alterative
versions of this package, differing in how redundant probe sets are
dealt with.  In the curatedOvarianData version, each gene is represented
by the gene with maximum mean.  In NormalizerVcuratedOvarianData, each
gene is represented by the mean of the probesets after removing "noisy"
probesets (see the Normalizer function of the Sleipnir library for
computational biology), and in FULLVcuratedOvarianData, no collapsing of
probe sets is done, but a map is provided to allow the user to do so
by their method of choice through featureData(eset).

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

Maintainer: Levi Waldron <levi@jimmy.harvard.edu>
}

\examples{
##List all datasets:
data(package="curatedOvarianData")
##
##See the actual template used for syntax checking of clinical metadata:
template.file <- system.file("extdata/template_ov.csv", package = "curatedOvarianData")
template <- read.csv(template.file, as.is=TRUE)
head(template)
}
