# SC.SkagerrakMTB

This repository contains data and code associated with the publication:

**"Navigating Past Oceans: Comparing Metabarcoding and Metagenomics of Marine Ancient Sediment Environmental DNA"**  
*Holman LE, Zampirolo G, Gyllencreutz R, Scourse J, Frøslev T, Carøe C, Gopalakrishnan S, Pedersen MW, Bohmann K (2025)*  
Published in *Molecular Ecology Resources*  
📄 [Read the paper](https://doi.org/10.1111/1755-0998.14086)

---

## Repository Structure

File/Folder | Description
--- | ---
`0.1.Bioinformatics.R` | R-based processing and formatting of metabarcoding eDNA data for downstream analysis
`0.1.TaxonomicAssignment.sh` | Shell script for assigning taxonomy to metabarcoding reads using reference databases
`0.2.FilePrep..R` | Formats and merges metabarcoding read tables for analysis
`0.3_Giulia_Damage_profiles_MTG.R` | Generates DNA damage profiles for metagenomic data
`0.4.1_Validation_Outlier.sh` | Shell script to detect and inspect outlier samples
`0.4.2_SupplementaryMTG.R` | R script to generate supplementary analyses for MTG datasets
`1.Analysis.R` | Main analysis script: diversity estimates, comparative analyses, and figure generation
`filtered_metazoan_all.csv` | Filtered metazoan taxa table (metagenomics)
`filtered_genus_all.csv` | Filtered genus-level taxa table (metagenomics)
`MTB.metadata.csv` | Metadata file for metabarcoding samples
`MTB.reads.csv` | Read count table for metabarcoding samples
`MTG.metadataAge.csv` | Metadata and age model data for metagenomics samples
`MTG.reads.csv` | Read count table (metagenomics, CSV version)
`MTG.reads.xlsx` | Read count table (metagenomics, Excel version)
`rawData/` | Raw or minimally processed input files (e.g. unfiltered outputs)
`cleanedData/` | Cleaned and filtered data used in analyses
`figures/` | Output figures included in the manuscript and supplementary materials
`taxonomy/` | Taxonomy files for metabarcoding
`README.md` | This file

---

