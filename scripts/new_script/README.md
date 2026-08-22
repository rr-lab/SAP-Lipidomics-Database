# Legacy lipid preprocessing and spatial-adjustment scripts

This directory contains the data-processing workflow through SpATS only.
It excludes GWAS, class composition, ratio, network, and figure-generation scripts.

Order is not fully executable as a single pipeline without checking paths and inputs. The original legacy files remain the authoritative record.

1. Raw QC inspection
2. Legacy SERRF correction
3. Post-SERRF wrangling
4. Legacy SpATS correction
5. Raw annotation audit
6. Corrected-name reaggregation
7. Duplicate intensity-profile audit
8. Assignment of nonfocused classes
9. Removal of non-lipid annotations
10. SERRF input preparation
11. SERRF feature-traceability audit
12. Pre-SERRF filtering/input build
13. Post-SERRF annotation and reaggregation
14. Post-SERRF SpATS fitting
