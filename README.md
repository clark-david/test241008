ICDPIC -- International Classification of Diseases Programs for Injury Categorization was originally developed in 2008 using Stata and using ICD Version 9 Clinical Modification (ICD-9-CM) diagnosis codes, providing an easy way to convert these codes to standard injury severity scores and categories.  After the introduction of ICD-10-CM to US hospitals in 2015, an update to accommodate this change was developed using R, and was eventually published on CRAN as package icdpicr.

Package icdpicr2 is being developed as a successor to icdpicr.  The "Implementation" vignette in this package gives more detail about the history of injury severity scoring, ICDPIC, and icdpicr; it discusses issues with the current version of icdpic and the aims for icdpicr2.

Package icdpicr2 consists primarily of a single function, cat_trauma2.  This function first reads in user data in a specified format.  It then calculates AIS, ISS, NISS, mortality predictions from the TQP and NIS models, and injury mechanisms.  It returns the original data file with these fields added.  Further details about the options available in cat_trauma2 are provided in the help file for this function.

