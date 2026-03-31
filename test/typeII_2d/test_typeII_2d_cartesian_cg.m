%% Basic Settings.
clear;
close all;
warning off;

p = 5;

tol_cg = 1e-8;
maxit = 100;

fprintf("tol_cg: %.1e\n", tol_cg);
fprintf("maxit: %d\n", maxit);

%% Generate points x and omega.
n = 2^p;
nx = n;
ny = n;
N = nx * ny;

x = CartesianGrid(n);
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
fprintf("cond(A): %.1e\n", kappa_A);

%% MRI Reconstruction.
P = phantom('Modified Shepp-Logan', n);
u_ex = reshape(P, N, []);
f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);

tic;
[u_solve, flag, relres, iter, resvec] = INUDFT2_2D_CG(x, nx, ny, f_nufft, tol_cg, maxit);
t_solve = toc;
fprintf("Solve time: %.1e\n", t_solve);
fprintf("Iter number: %d\n", iter);
u_solve = real(u_solve);

r = f_nufft - MY_NUFFT2_2D(u_solve, x, nx, ny);
rel_res = norm(r) / norm(f_nufft);
fprintf("Rel res: %e\n", rel_res);

e = u_ex - u_solve;
rel_acc = norm(e) / norm(u_ex);
fprintf("Rel acc: %.1e\n", rel_acc);

P_reconstruct = reshape(u_solve, n, n);

figure;
subplot(1, 2, 1);
imagesc(P);
colormap gray;
axis image;
axis off;
title('Original Shepp-Logan Phantom', 'FontSize', 12);
colorbar;
subplot(1, 2, 2);
imagesc(P_reconstruct);
colormap gray;
axis image;
axis off;
title('Reconstructed Phantom', 'FontSize', 12);
colorbar;
sgtitle('Reconstruction Comparison', ...
    'FontSize', 14, 'FontWeight', 'bold');

