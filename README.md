# NUFFT2D

MATLAB research code for fast direct inversion of two-dimensional type-II nonuniform discrete Fourier transforms (NUDFTs).

## Mathematical problem

For two-dimensional type-II NUDFT, this project uses

```math
f_j
=
\sum_{k^{[x]}=0}^{n^{[x]}-1}
\sum_{k^{[y]}=0}^{n^{[y]}-1}
e^{-2\pi \mathrm{i}
\left(k^{[x]}x_j+k^{[y]}y_j\right)}
c_{k^{[x]},k^{[y]}},
\qquad
0 \leq j \leq M-1.
```

where the *sample points* $\{(x_{j}, y_{j})\}$  are distributed arbitrarily in $[0,1)^{2}$ and the *frequencies* $\{(k^{[x]}, k^{[y]})\}$ are distributed on a Cartesian grid of contiguous integers.
The *coefficients* $\{c_{k^{[x]}, k^{[y]}}\}$ and the target values $\{f_{j}\}$ are assumed to be in $\mathbb{C}$.
The number of sample points and the number of frequencies are $M$ an $N = n^{[x]} n^{[y]}$ respectively. 

The inversion routines recover the coefficient vector or image from the nonuniform samples.
Rectangular systems ($M \geq N$) are supported and are treated as least-squares problems.

## Installation

Clone the repository:

```bash
git clone https://github.com/JingyuLiuMath/NUFFT2D.git
cd NUFFT2D
```

### Build FFTW and FINUFFT

Detailed dependency instructions are in [`extern/readme.md`](extern/readme.md).

### Set the MATLAB path

Start MATLAB in the repository root and run:

```matlab
nufft2d_startup();
```

## Running the tests

From the repository root:

```matlab
nufft2d_startup();

% 2D random and polar sampling examples.
run("test/typeII_2d/test_typeII_2d_rand.m");
run("test/typeII_2d/test_typeII_2d_polar.m");
```

## Reproducing the experiments

Experiment families are stored in:

- `experiments/typeII_2d_rand/` — 2D random-point experiments.
- `experiments/typeII_2d_polar/` — 2D polar-grid experiments.

## Developers

- [Yingzhou Li](https://yingzhouli.com/), School of Mathematical Sciences, Fudan University.
- [Jingyu Liu](https://jingyuliumath.github.io/), School of Mathematical Sciences, Fudan University.

## Related papers

- Yingzhou Li and Jingyu Liu, A Superfast Direct Solver for 2D Type-II Inverse Nonuniform Discrete Fourier Transform Based on Hierarchically Semiseparable Matrix. [[arXiv]](https://arxiv.org/abs/2607.00928)
