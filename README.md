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

## Citation

If you use this code in your research, please cite

```bibtex
% The BibTeX entry will be added after publication.
```

## License

This project is distributed under the **GNU Lesser General Public License v3.0 (LGPL-3.0)**.

Please refer to the LICENSE file for details.


