### Changes in version 1.5.2
* Fixed bugs in `createEsets.R`

### Changes in version 1.3.2
* Updated all gene maps

### Changes in version 1.3.1
* Added missing grade to GSE26913
* Improved order of gene filtering in `createEsetList.R`

### Changes in version 1.3.0
* Added datasets: GSE26193, GSE49997, GSE44104, GSE51088
* Curated debulking status for Mok (GSE18520)
* Added sample_type classifications
* Refined platform_accession and metadata column names
* Upgraded to Ensembl release 76 for gene maps

### Changes in version 1.2.2
* Added "Non-unique gene symbols" section to vignette
* Updated `createEsetList.R` and patient selection config

### Changes in version 1.2.1
* Fixed curation of TCGA sample type "11"

### Changes in version 1.2.0
* Added RNA-seq TCGA data
* Updated gene maps

### Changes in version 1.1.1
* Added options `remove.retracted` and `remove.subsets` to `createEsets.R`

### Changes in version 1.1.0
* Added p-values to figure captions
* Renamed package variants in `inst/extdata/*`
* Updated metadata for GSE18520
* Added GSE8842 (Agilent array)
* Added DATABASE publication citation

### Changes in version 1.0.1
* Curated grade for GSE17260

### Changes in version 1.0.0
* Automated per-dataset documentation
* Expanded vignette
* Added missing grade/summarygrade fields
* Enriched `ExpressionSets` with metadata
* Fixed TCGA experiment name typo
* Fixed `featureData` in FULLV version

### Changes in version 0.99.21
* Improved merging logic for technical replicate metadata

### Changes in version 0.99.20
* Added stromal/tumor/normal cell percentages
* Added `uncurated_author_metadata`
* Removed obsolete metadata fields
* Renamed metadata fields
* Defined TCGA optimal debulking as <10mm

### Changes in version 0.99.10
* Added TCGA-mirna-8x15kv2 and GSE30161

### Changes in version 0.99.2
* Added tumorstage and summarystage to GSE26712

### Changes in version 0.99.1
* Fixed expO short title
* Updated staging system in template
* Curated stage/grade fields in multiple datasets
* Added tumorstage to GSE2109
* Removed curation files for uncurated datasets
* Fixed stage data in GSE6008
* Updated BLAST reference genome to GRCh37.68
* Updated BioMart gene maps
* Addressed missing publication metadata
* Made `test_counts.R` work
* Handled metadata disagreement in technical replicates
* Added batch variable to GSE17290060

### Changes in version 0.5.0
* Disallowed "unknown" in templates; use NA instead
* Added syntax-checking unit test

### Changes in version 0.4.2
* Actually included ages for GSE26712

### Changes in version 0.4.1
* Added ages to GSE26712
* Used `getPlatformsBiomaRt` for GSE6822

### Changes in version 0.4.0
* Added GSE30009 (Taqman qPCR dataset)

### Changes in version 0.3.0
* Added GSE32062.GPL6480 and GSE32063
* GSE13876: removed redundant recurrence status
* Resolved missing `ExperimentData` for GSE2109 and GSE6822
* Added subtype information to multiple datasets
* Converted dataset check test to R CMD check unit test

### Changes in version 0.2.4
* TCGA: replaced missing batch variable
* GSE9891: corrected survival data for three patients

### Changes in version 0.2.3
* Renamed batch variables for consistency
* Consistently label missing debulking status as NA
* Replaced missing probesets in FULL version
* Added unit test for sample counts (not yet functional)

### Changes in version 0.2.2
* Replace GSE17260, which was missing in 0.2.1

### Changes in version 0.2.1
* In the TCGA dataset, about 40 censored patients had `days_to_death` and `days_to_recurrence` labelled as NA instead of follow-up time.
* TCGA released updated clinical information for six patients
* Minor changes to GSE18520 based on a review, including updates to debulking status, stage, and grade.

### Changes in version 0.2.0
* Added E-MTAB-386 (FFPE samples with survival data)
* Put MaxMean probe mappings in `eset@featureData`
* Put official Bioconductor platform names in `eset@annotation`
