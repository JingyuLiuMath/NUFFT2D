%% Basic Settings.
clear;
close all;
warning off;

p = 5;

tol = 1e-8;
maxit = 100;
fprintf("tol: %.4e\n", tol);
fprintf("maxit: %d\n", maxit);

%% Generate points x and omega.
nx = 2^p;
ny = 2^p;
N = nx * ny;
M = 2 * N;
fprintf("M: %d, N: %d\n", M, N);

% Random
x = rand(M, 2);

%% Solution: RHS (approximately) in range(A).
fprintf("RHS (approximately) in range(A)\n");
u_ex = randn(N, 1) + randn(N, 1) * 1i;
f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);

tic;
[u_solve, flag, relres, iter, resvec] = INUDFT2_2D_CG(x, nx, ny, f_nufft, tol, maxit);
t_solve = toc;
fprintf("Solve time: %.4e\n", t_solve);

r = f_nufft - MY_NUFFT2_2D(u_solve, x, nx, ny);
rel_res = norm(r) / norm(f_nufft);
fprintf("Rel res: %e\n", rel_res);
