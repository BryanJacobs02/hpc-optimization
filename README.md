# RNA-Seq Pipeline Optimization in Murine Colon Crypts
*Final project for [MCB 516A: Bioinformatics and Functional Genomic Analysis](https://infosci.arizona.edu/course/mcb-516a-bioinformatics-and-functional-genomic-analysis) — University of Arizona*
  
  
## Project Intro/Objective
This project explores how changes to high-performance computing (HPC) parameters during RNA-seq preprocessing affect downstream analysis of gene expression. Using a multiomics dataset of FACS-sorted mouse colonic epithelial cells, the study specifically investigates how varying fastp parameters including sliding window size, PHRED score threshold, and maximum read length influence alignment quality and transcriptomic interpretation. The biological context centers on the transition from stem cells to fully differentiated enterocytes in the colon, a process known to be regulated post-transcriptionally.

By aligning reads with STAR and comparing expression profiles across parameter settings, we evaluated the robustness and sensitivity of differential expression results. Despite modifications in preprocessing, overall trends remained consistent, suggesting that the dataset and analysis pipeline are relatively resilient to reasonable shifts in HPC parameterization. This work emphasizes the importance of optimizing preprocessing settings to balance data quality, computational efficiency, and biological insight in transcriptomic workflows.  

## Authors
**Bryan Jacobs**  
**Katharine Johnson**  
**Alexander Borowiec**  
  
  
## Languages/Packages:
* R
   * DESeq2
   * Enhanced Volcano
   * NMF
   * dplyr
* Bash  
  

## Tools, Techniques and Softwares
* University of Arizona High Performance Computing - Ocelote Cluster
* SLURM
* FastQC
* fastp
* STAR
* featureCounts
* MultiQC
* Principal Component Analysis
* Heatmaps
* Volcano Plots
* g:Profiler
* GSEA
* Cytoscape Enrichment Maps
  
  
## Repository Structure
- **`data/`**: Contains `.csv` files of converted gene lists necessary for functional analysis. Please contact me for access to featureCounts files derived from HPC.
- **`code/`**: Contains `.r` script files for differential expression and functional analysis as well as `.slurm` script files for high performance computing workflows.
- **`deliverables/`**: Contains `.pdf` files of final report and slide show.
- **`README.md`**
  
  
## How To Run
#### For Simple Viewing
1. Download and open desired `.pdf` files from the `deliverables/` folder.

#### To Run R Scripts
1. Clone the repository.
2. Contact me at bryanj82202@gmail.com for access to featureCounts files (too large for GitHub).
3. Download the desired `.r` file and open with a compatible software.
4. Install the packages listed at the top of the script.
5. Run the code as usual.
