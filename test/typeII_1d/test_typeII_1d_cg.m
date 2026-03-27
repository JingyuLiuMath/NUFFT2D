%% Basic Settings.
clear;
close all;
warning off;

p = 10;

tol = 1e-8;
maxit = 100;
fprintf("tol: %.4e\n", tol);
fprintf("maxit: %d\n", maxit);

%% Generate points x and omega.
N = 2^p;
M = 2 * N;
fprintf("M: %d, N: %d\n", M, N);

% Perturbation.
% perturb_size = 0.2;
% x = (0 : (M - 1))' / M;
% x = x + (rand(M, 1) * 2 - 1) * perturb_size / M;
% x = sort(x);

% Random
x = sort(rand(M, 1));

%% Solution: RHS (approximately) in range(A).
fprintf("RHS (approximately) in range(A)\n");
u_ex = randn(N, 1) + randn(N, 1) * 1i;
f_nufft = MY_NUFFT2(u_ex, x, N);

tic;
[u_solve, flag, relres, iter, resvec] = INUDFT2_CG(x, N, f_nufft, tol, maxit);
t_solve = toc;
fprintf("Solve time: %.4e\n", t_solve);

r = f_nufft - MY_NUFFT2(u_solve, x, N);
rel_res = norm(r) / norm(f_nufft);
fprintf("Rel res: %e\n", rel_res);

e = u_ex - u_solve;
rel_acc = norm(e) / norm(u_ex);
fprintf("Rel acc: %.4e\n", rel_acc);
