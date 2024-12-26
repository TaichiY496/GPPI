# GPPI: Gaussian Process Phase Interpolation

## Outline
This repository provides the codes of the paper below:
```
Gaussian Process Phase Interpolation for estimating the asymptotic phase of a limit cycle oscillator from time series data
Taichi Yamamoto, Hiroya Nakao, and Ryota Kobayashi  
Chaos, Solitons & Fractals, Volume 191, February 2025, 115913
```
[View paper](https://doi.org/10.1016/j.chaos.2024.115913)

## Contents
The codes implement 2 methods for estimating asymptotic phase:
- GPPI method by T.Yamamoto (proposed)
- DPR method by N.Namura (for comparison)

## Requirement
The environment in original experiment is as below.
- MATLAB 2023b.
- Statistics and Machine Learning Toolbox

## Usage
`test.mlx` reproduces the experiment of Stuart-Landau oscillator in our paper. (see Sec. 4.2.)

## How to cite our work
If you use GPPI part of this code, please cite our paper:
```
@article{YAMAMOTO2025115913,
  title = {Gaussian Process Phase Interpolation for estimating the asymptotic phase of a limit cycle oscillator from time series data},
  journal = {Chaos, Solitons & Fractals},
  volume = {191},
  pages = {115913},
  year = {2025},
  issn = {0960-0779},
  doi = {https://doi.org/10.1016/j.chaos.2024.115913},
  url = {https://www.sciencedirect.com/science/article/pii/S0960077924014656},
  author = {Taichi Yamamoto and Hiroya Nakao and Ryota Kobayashi},
  keywords = {Synchronization, Limit cycle oscillators, Phase reduction, Machine learning, Gaussian process regression}
}
```

If you use DPR part of this code, please cite the paper by N.Namura. et al.:
```
@article{PhysRevE.106.014204,
  title = {Estimating asymptotic phase and amplitude functions of limit-cycle oscillators from time series data},
  author = {Namura, Norihisa and Takata, Shohei and Yamaguchi, Katsunori and Kobayashi, Ryota and Nakao, Hiroya},
  journal = {Phys. Rev. E},
  volume = {106},
  issue = {1},
  pages = {014204},
  numpages = {12},
  year = {2022},
  month = {Jul},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevE.106.014204},
  url = {https://link.aps.org/doi/10.1103/PhysRevE.106.014204}
}
```

## The program was developed by
- Taichi Yamamoto (GPPI part)
- Norihisa Namura (DPR part).

We thank Norihisa Namura for sharing us the code of DPR method.
