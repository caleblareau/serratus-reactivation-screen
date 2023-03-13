# serratus-reactivation-screen


# Workflow for building table of all human viral hits using Serratus

(assumes you are in the `data` folder of this repository)


## 1) Set up a set of viruses to query
Here, we use the list provided by ViralZone (`https://viralzone.expasy.org/resources/Table_human_viruses.txt?`), accessed on August 1, 2022. At the time of processing, there were 129 viruses pulled from this resource. 
This file is downloaded locally to the repository at `data/Table_human_viruses.txt`

## 2) Setup Serratus API calls

```
Rscript step2_create_urls.R
```

## 3) Execute Serratus API calls

Now that we've built all possible nuccore IDs pertaining to our viruses of interest, we just run a bunch of serial `wget` commands. 

```
sh step3_run_wget_cmds.sh
```

## 4) Remove entries with 0 reads

Since most of the `NC` IDs aren't the version in Serratus, they will yield 0 lines from the API call. We'll quickly clean these up:

Code taken/modified from this [StackOverflow Post](https://unix.stackexchange.com/questions/327371/fast-way-to-delete-files-with-less-than-x-lines)
```
find . -type f -exec awk -v x=2 'NR==x{exit 1}' {} \; -exec rm -f {} \;
```

To keep things tidy, we'll move the resulting files to 

```
mkdir nc_pulls
mv NC*.txt nc_pulls
```

## 5) Aggregate summary table

Here, we loop over all the `NC` files per virus and count them. Further, we annotate whether the virus is a reactivation candidate. 

```
Rscript step5_aggregate_count.R
```

This generates `output/Serratus_hits_all_viruses.tsv`, which is used for the `Supplementary Table 1` in the revised manuscript. 