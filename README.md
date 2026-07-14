# HHO-Nonlocal-Kirchhoff

MATLAB implementation of the Hybrid High-Order (HHO) method for nonlinear nonlocal problems of Kirchhoff type.

## Overview

This repository contains the MATLAB implementation accompanying the paper

> **A Hybrid High-Order Finite Element Method for a Nonlocal Nonlinear Problem of Kirchhoff Type**

by **Gouranga Mallik**.

The code implements the Hybrid High-Order (HHO) method for nonlinear nonlocal Kirchhoff problems and reproduces the numerical results presented in the paper.

## Original HHO Implementation

This work is based on the publicly available MATLAB implementation

**HHO-Lapl-OM**

developed by **Jérôme Droniou, Daniel Anderson**, and collaborators.

Original repository: https://github.com/jdroniou/HHO-Lapl-OM

The original implementation solves the Poisson model problem using the Hybrid High-Order method and is distributed under the GNU Lesser General Public License (LGPL v3).

The present repository extends this framework to nonlinear nonlocal problems of Kirchhoff type. While the original repository implements the HHO method for the Poisson model problem, this repository develops the HHO discretization, nonlinear solution algorithms, and numerical experiments for the nonlocal Kirchhoff problem.

## Main Contributions

Building upon the original HHO implementation, this repository introduces

- implementation of the nonlinear nonlocal Kirchhoff model;
- nonlinear solution algorithms based on Newton's method and the Sherman–Morrison–Woodbury formula;
- modifications to the local and global assembly procedures;
- numerical experiments presented in the accompanying paper.

## Running the Code

Run

```matlab
runHHO_NonLocal_Kirchhoff_Test.m
```
to reproduce the numerical experiments presented in the paper.

## Paper

The preprint accompanying this repository is available on arXiv:

**A Hybrid High-Order Finite Element Method for a Nonlocal Nonlinear Problem of Kirchhoff Type**

Arxiv version:
https://arxiv.org/pdf/2510.15574

Applied Numerical Mathematics - Journal Publication:
https://doi.org/10.1016/j.apnum.2026.07.005

## Citation

If you use this software or any part of its source code, or build upon the methods described in the accompanying paper, please cite:

```bibtex
@article{Mallik2025,
  author  = {Gouranga Mallik},
  title   = {A Hybrid High-Order Finite Element Method for a Nonlocal Nonlinear Problem of Kirchhoff Type},
  journal = {Applied Numerical Mathematics},
  year    = {2026},
  url     = {https://doi.org/10.1016/j.apnum.2026.07.005}
}
```

Once the paper is published, this citation will be updated with the journal information.

## License

This project is distributed under the **GNU Lesser General Public License v3.0 (LGPL-3.0)**.

Please refer to the LICENSE file for details.


