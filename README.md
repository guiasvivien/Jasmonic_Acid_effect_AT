# Jasmonic Acid Effect on *Arabidopsis thaliana* (JA_AT Study)

## Overview
This project investigates the prolonged effects of Jasmonic Acid (JA) treatment on *Arabidopsis thaliana* using publicly available transcriptomic data from **Hickman et al. (2017)**. The primary objective is to assess whether **timepoint 16 (T16)** of JA treatment resembles **timepoint 0 (T0)**, potentially indicating the resolution of the stress response.

## Data Source
The dataset originates from the study by **Hickman et al.**, which analyzed JA-mediated transcriptional responses in *Arabidopsis thaliana*. The full article can be accessed here:  
🔗 **[PubMed - Hickman et al. (2017)](https://pubmed.ncbi.nlm.nih.gov/28827376/)**

## Project Goals
- Extract and analyze key transcriptomic data from **eight selected timepoints**.
- Evaluate whether the expression profile at **T16** returns to its baseline (**T0**).
- Implement a computational workflow to streamline data processing and visualization.

## Implementation
- A custom **Python/R script** was developed to extract and analyze relevant gene expression data.
- The script **filters and processes** transcriptomic datasets, focusing on stress-response genes.
- Statistical and visualization tools were applied to compare gene expression across timepoints.

## How to Use
1. Clone the repository:
   ```bash
   git clone https://github.com/guiasvivien/Jasmonic_Acid_AT.git
   cd Jasmonic_Acid_AT
