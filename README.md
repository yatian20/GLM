# Efficient Analysis of Latent Spaces in Heterogeneous Networks

## Project Description
This repository contains the implementation code for the paper "Efficient Analysis of Latent Spaces in Heterogeneous Networks".

## Repository Structure

### Core Algorithms (`algorithms/`)
- `latent_vectors_est.R`: Implements Algorithms A.1-A.2 and 1-2 proposed in the paper, and the algorithm for Section L in the Supplementary Material, and the compared method of MASE under COSIE model.
- `latent_dimensions_est.R`: Implements Algorithm B.1 proposed in the paper.

### Simulation Studies (`simulation/`)
- `experiments.R`: Data generation and analysis for Section 5 (under three distributions and Cases (A)-(C), fixing T = 5 as an example).
- `results/`: Output files by `experiments.R` across T = 5/10/20/40/80.
- `plot.R`: Visualization code for simulation results, producing Figures 1-3 and S1-S3 in the paper.
- `figures/`: Output files by `plot.R`.

### Additional Simulation Studies (`simulation2/`)
- `experiments.R`: Data generation and analysis for Sections B.3, G.2, G.3 and L.2 in the Supplementary Material. 
- `results/`: Output files by `experiments.R`.
- `plot.R`: Visualization code for additional simulation results, producing Figures S4-S8, S13-14, and Tables S1-S2 in the Supplementary Material.
- `figures/`: Output files by `plot.R`.

### Real Data Analysis (`real_data/`)
- `raw_data/`: Lazega lawyers data and documentation (original source: https://www.stats.ox.ac.uk/~snijders/siena/Lazega_lawyers_data.htm). 
- `analyze.R`: Data processing and analysis for Section 6. 
- `analyze_SecL.R`: Data processing and analysis for Supplementary Material Section L. 
- `results/`: Output files by `analyze.R` and `analyze_SecL.R`.
- `plot.R`: Visualization code for data analysis, producing Figures 4-6 in the paper.
- `plot2.R`: Supplementary visualization code for data analysis, producing Figures S10-S12 and Tables S3-S4 in the Supplementary Material.
- `plot_SecL.R`: Supplementary visualization code for real data analysis link prediction, producing Figure S15 in the Supplementary Material.
- `figures/`: Output files by `plot.R`, `plot2.R`, `plot_SecL.R`.
