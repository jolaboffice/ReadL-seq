# ReadL-seq

Single-molecule DNA lesion mapping via nanopore signal-drop detection



## Overview

ReadL-seq is a nanopore-based framework for detecting DNA damage at single-read resolution.  
DNA lesions are converted into biotin-labeled tracts through polymerase-mediated repair synthesis and are detected via characteristic signal drops and alignment features in nanopore sequencing data.

Signal deviations and alignment features are integrated into a cumulative metric (BiotinScore), enabling robust identification of DNA damage sites from individual reads without amplification or ensemble averaging.



## Installation

Required tools:

- Remora v2.1.3 (Oxford Nanopore Technologies)
- Dorado
- minimap2
- samtools
- Python (numpy, pandas)



## Remora modification

This pipeline requires minor modifications to the Remora source code.

After installing Remora:

1. Locate `io.py` in the Remora installation directory and replace the corresponding function with the version provided in:

scripts/remora_patch/io.py

2. Locate `data_chunks.py` in the Remora installation directory and insert the function provided in:

scripts/remora_patch/data_chunks.py

Remora is not included in this repository and must be installed separately.


## Usage

Run the example pipeline:

```bash
git clone https://github.com/jolaboffice/ReadL-seq.git
cd ReadL-seq
bash run/run_example.sh
```

This script performs:

1. Basecalling (Dorado)  
2. Alignment (minimap2)  
3. BiotinScore-based analysis using `biotin_ssb_pipeline.py`  

## Workflow structure

Two execution modes are provided:

### Stepwise workflow

- extract_perbase_signal.py  
- mean_sigma.py  
- compute_delta_signal.py  
- biotin_scoring.py  

These scripts perform signal extraction, baseline generation, Δsignal calculation, and BiotinScore computation as separate steps.

The stepwise workflow primarily reflects the experimental baseline approach, but can also be used with an ONT-predicted baseline if provided in the same format as the output of `mean_sigma.py`.



### Integrated workflow

- biotin_ssb_pipeline.py  

This script performs the full BiotinScore pipeline in a single step using the ONT-predicted baseline derived from the Remora k-mer level table.



## Method summary

Signal deviation is defined as:

Δsignal_i = μ_i − signal_i  

where μ_i is the baseline signal at position i.

BiotinScore is calculated as the sum of Δsignal values within a window spanning −3 to +100 bases relative to the 3′ aligned end of each read.

Biotin-labeled reads exhibit sustained positive Δsignal values across this region due to consecutive biotin incorporation and extended 3′ soft-clipped regions.


## Output

The pipeline generates:

- Per-position CSV files containing coverage and biotin-positive hit counts  
- Per-read XLSX file containing BiotinScore values and classification results  


## Example data

The repository includes example λ DNA data (POD5 and FASTA) and corresponding output files generated from the example pipeline.



## Notes

- The file `resources/levels.txt` contains the ONT k-mer signal level table obtained from Remora.  
- This table provides expected normalized signal values for each k-mer and is used for ONT-predicted baseline estimation.  
- This enables control sequencing–free analysis.  


## License

MIT License


## Third-party dependency

This repository uses and modifies functions from Remora (Oxford Nanopore Technologies),  
which is distributed under the ONT Public Licence.

