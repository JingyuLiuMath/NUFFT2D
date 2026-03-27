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
n = 2^p;
nx = n;
ny = n;
N = nx * ny;

% Random
x = PolarGrid(n);
M = size(x, 1);

fprintf("M: %d, N: %d\n", M, N);

figure();
plot(x(:, 1), x(:, 2), "Marker", ".", "LineStyle", "none");
axis equal;
xlim([0, 1]);
ylim([0, 1]);

%% NUDFT2_Matrix.
A = NUDFT2_2D_Matrix(x, nx, ny);
kappa_A = cond(A);
fprintf("cond(A): %.4e\n", kappa_A);

%% Solution: RHS (approximately) in range(A).
fprintf("RHS (approximately) in range(A)\n");
u_ex = randn(N, 1) + randn(N, 1) * 1i;
f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);

tic;
[u_solve, flag, relres, iter, resvec] = INUDFT2_2D_CG(x, nx, ny, f_nufft, tol, maxit);
t_solve = toc;
fprintf("Solve time: %.4e\n", t_solve);
fprintf("Iter number: %d\n", iter);

r = f_nufft - MY_NUFFT2_2D(u_solve, x, nx, ny);
rel_res = norm(r) / norm(f_nufft);
fprintf("Rel res: %e\n", rel_res);

e = u_ex - u_solve;
rel_acc = norm(e) / norm(u_ex);
fprintf("Rel acc: %.4e\n", rel_acc);
