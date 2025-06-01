# Efficient Analysis of Latent Spaces in Heterogeneous Networks

## Project Description
This repository contains the implementation code for the paper "Efficient Analysis of Latent Spaces in Heterogeneous Networks" by Yuang Tian, Jiajin Sun, and Yinqiu He.

## Repository Structure

### Core Algorithms (`algorithms/`)
- `latent_vectors_est.R`: Implements Algorithms A.1-A.2 and 1-2 in the paper, and MASE under COSIE model.
- `latent_dimensions_est.R`: Implements Algorithm B.1 in the paper.

### Simulation Studies (`simulation/`)
- `experiments.R`: Data generation and analysis under three distributions and Cases (A)-(C) (only present T = 5 as an example).
- `results/`: Output files by `experiments.R`.
- `plot.R`: Visualization code for simulation results, producing Figures 1-3 and S1-S3 in the paper.
- `figures/`: Output files by `plot.R`.

### Real Data Analysis (`real_data/`)
- `raw_data/`: Lazega lawyers data (original source: https://www.stats.ox.ac.uk/~snijders/siena/Lazega_lawyers_data.htm). 
- `analyze.R`: Data processing and analysis. 
- `results/`: Output files by `analyze.R`.
- `plot.R`: Visualization code for data analysis, producing Figures 4-6 in the paper.
- `figures/`: Output files by `plot.R`.
