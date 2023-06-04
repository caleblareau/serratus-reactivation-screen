# serratus-reactivation-screen

Last updated: June 3, 2023

## Step 1 - pull libraries mapping to viruses using ViralZone and Serratus

The `serratus_data_setup` folder contains the workflow, including all relevant substeps to systematically query Serratus through the API based on viruses annotated through ViralZone.

## Step 2 - annotate libraries with relevant meta data. 

The `analysis_code` folder contains necessary `R` scripts for post-processing the outputs of the `serratus_data_setup` to reproduce relevant supplementary tables
and figures associated with our computational screen. These are annotated with the `GEOquery::getGEO` function. 

## Other notes

The `metadata` folder contains thousands of `.soft` files used to pull the submitter-associated metadata. 
These were downloaded programmatically using the `GEOquery::getGEO()` function in the `analysis_code/02_hhv6_followup.R` folder. 

The `output` folder represents finalized, processed data, including plots used in the `Lareau et al.` HHV-6 manuscript.

## Mapping to supplemental tables

Here, we explicitly map the files in `output` used to create the extended data tables found in the manuscript. 

```
Extended Data Table 1 = output/Serratus_hits_all_viruses.tsv
Extended Data Table 2 = output/Serratus_libraries_all_hits_reactivation.tsv
Extended Data Table 4 = output/tcell_hits_per_virus.tsv
Extended Data Table 7 = output/iPSC_hits_per_virus.tsv

```

### Other useful tables
```
All iPSC SRAs = output/all_iPSC_SRAs.tsv
All Tcell SRAs = output/all_tcell_SRAs.tsv
```
<br>
