\name{curatedOvarianData-package}
\alias{curatedOvarianData}
\docType{package}
\title{Ovarian Cancer Gene Expression Analysis}
\description{The curatedOvarianData package provides relevant functions and data for gene expression analysis in patients with ovarian cancer.}
\details{
\tabular{ll}{
Package: \tab curatedOvarianData\cr
Type: \tab Package\cr
Version: \tab \Sexpr{packageDescription("curatedOvarianData")$Version}\cr
Date: \tab \Sexpr{packageDescription("curatedOvarianData")$Date}\cr
License: \tab Artistic-2.0\cr
Depends: \tab R (>= 2.10.0), affy\cr
}

In the "Available sample meta-data" sections of each dataset, please refer to the following key.

For "subtype": ser=serous;endo=endometrioid;endo_clearcell=mixture of ser+endo.  Other includes sarcomatoid, endometroid, papillary.

For "primarysite" and for "arrayedsite": ov=ovary;ft=fallopian tube.

For "summarygrade": low = 1, 2, LMP.  High= 3,4,23.

For "summarystage": early = 1,2, 12. late=3,4,23,34.

For "tumorstage": FIGO Stage (I-IV, but coded here as 1-4 to ensure
correct ordering in factors).  If multiple stages given (eg 34), use the
highest.

For "T": TNM Stage (1-3).  If multiple stages given (eg 34), use the highest.

For "substage": substage (abcd).  For cases like ab, bc, use highest given.

For "G": Grade (1-4): If multiple given, ie 12, 23, use highest given.

For "N": N (0/1): degree of spread to regional lymph nodes.

For "M": M (0/1): presence of metastasis.

For "pltx": patient treated with platin.

For "tax": patient treated with taxol.

For "neo": patient treated with neoadjuvant treatment.

For "primary_therapy_outcome_success": response to any kind of therapy (including radiation only).

For "chemo_response": platinum resistance: refractory=3mo or less, resistant=6mo or less, sensitive=12mo or higher.

For "inferred_chemo_response": inferred platinum resistance: refractory=death in 6mo or less, sensitive=survival for 12mo or more. 

For "debulking": amount of residual disease (optimal = <1mm, suboptimal=>1mm).
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
}
