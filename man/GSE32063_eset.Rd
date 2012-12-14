\name{GSE32063_eset}
\alias{GSE32063_eset}
\docType{data}
\title{Japanese dataset B: Immune-activation as a therapeutic direction for patients with high-risk ovarian cancer based on gene expression signature}
\description{Two hundred sixty patients who were diagnosed as advanced-stage high-grade serous ovarian cancer were analyzed in this study.
}
\usage{data(GSE32063_eset)}
\format{

Authors:  Yoshihara K, Tsunoda T, Shigemizu D, Fujiwara H, Hatae M, Fujiwara H, Masuzaki H, Katabuchi H, Kawakami Y, Okamoto A, Nogawa T, Matsumura N, Udagawa Y, Saito T, Itamochi H, Takano M, Miyagi E, Sudo T, Ushijima K, Iwase H, Seki H, Terao Y, Enomoto T, Mikami M, Akazawa K, Tsuda H, Moriya T, Tajima A, Inoue I, Tanaka K

Lab: 	The Japanese Serous Ovarian Cancer Study G

Contact: Kosuke Yoshihara <yoshikou@med.niigata-u.ac.jp>

Title: High-risk ovarian cancer based on 126-gene expression signature is uniquely characterized by downregulation of antigen presentation pathway.

Url: http://pubget.com/search?q=22241791

PubMedID: 22241791
}
\details{
Platform used: Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)

Processed data here are those provided by the original study authors.

assayData: 41108 features, 40 samples.

GEO_platform_accession: GPL6480


---------------------------------------------------------------------------------

Available sample meta-data: 

---------------------------------------------------------------------------------

sample_type
tumor 
   40 

subtype
ser 
 40 

summarygrade
high  low 
  17   23 

summarystage
late 
  40 

T
 3  4 
31  9 

substage
 b  c 
 3 28 

G
 2  3 
23 17 

pltx
 y 
40 

tax
 y 
40 

days_to_death   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    210     705    1155    1346    1792    3330 

vital_status
deceased   living 
      22       18 

debulking
   optimal suboptimal 
        19         21 

}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE32063
}
\examples{
data(GSE32063_eset)
if(require(affy)){
  summary(GSE32063_eset$vital_status)
  time <- GSE32063_eset$days_to_death / 365
  cens <- ifelse(GSE32063_eset$vital_status=="deceased",1,0)
  library(survival)
  fit <- survfit(Surv(time,cens)~-1)
  fit          #for summary of survival
  summary(fit)
  plot(fit,xlab="Time (years)",ylab="Survivor function")
  inverse.fit <- survfit(Surv(time,1-cens) ~ -1)
  inverse.fit  #for summary of follow-up time
}
}
\keyword{datasets}
