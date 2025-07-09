# Efficient Analysis of Latent Spaces in Heterogeneous Networks

## Project Description
This repository contains the implementation code for the paper "Efficient Analysis of Latent Spaces in Heterogeneous Networks" by Yuang Tian, Jiajin Sun, and Yinqiu He.

## Repository Structure

### Core Algorithms (`algorithms/`)
- `latent_vectors_est.R`: Implements Algorithms A.1-A.2 and 1-2 proposed in the paper, and a compared method MASE under COSIE model.
- `latent_dimensions_est.R`: Implements Algorithm B.1 proposed in the paper.

### Simulation Studies (`simulation/`)
- `experiments.R`: Data generation and analysis for Section 5 (under three distributions and Cases (A)-(C), fixing T = 5 as an example).
- `results/`: Output files by `experiments.R` across T = 5/10/20/40/80.
- `plot.R`: Visualization code for simulation results, producing Figures 1-3 and S1-S3 in the paper.
- `figures/`: Output files by `plot.R`.

### Additional Simulation Studies (`simulation2/`)
- `experiments.R`: Data generation and analysis for Sections B.3, G.2 and G.3 in the Supplementary Material. 
- `results/`: Output files by `experiments.R`.
- `plot.R`: Visualization code for additional simulation results, producing Figures S4-S8 and Tables S1-S2 in the Supplementary Material.
- `figures/`: Output files by `plot.R`.

### Real Data Analysis (`real_data/`)
- `raw_data/`: Lazega lawyers data (original source: https://www.stats.ox.ac.uk/~snijders/siena/Lazega_lawyers_data.htm). 
- `analyze.R`: Data processing and analysis for Section 6. 
- `results/`: Output files by `analyze.R`.
- `plot.R`: Visualization code for data analysis, producing Figures 4-6 in the paper.
- `plot2.R`: Supplementary visualization code for data analysis, producing Figures S10-S12 and Tables S3-S4 in the Supplementary Material.
- `figures/`: Output files by `plot.R` and `plot2.R`.
