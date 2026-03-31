%% Basic Settings.
clear;
close all;
warning off;

p = 5;

tol_hss = 1e-5;

fprintf("tol_hss: %.1e\n", tol_hss);

%% Generate points x and omega.
n = 2^p;
nx = n;
ny = n;
N = nx * ny;
M = 2 * N;

x = rand(M, 2);

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

%% NUDFT2.
min_points = p * n;
fprintf("min_points: %d\n", min_points);

tic;
A = NUDFT2_2D(nx, ny);
A.Construct_ID_Proxy(x, min_points, tol_hss);
t_construct = toc;
fprintf("Construct time: %.1e\n", t_construct);

r = A.Rank();
fprintf("HSS rank: %d\n", r);

mem_exact = M * N;
mem = A.Storage();
fprintf("Memory (GB): %.1e\n", double_to_gb(mem));
mem_ratio = mem / mem_exact;
fprintf("Mem ratio: %.1e\n", mem_ratio);

%% Apply.
u_ex = randn(N, 1) + randn(N, 1) * 1i;

tic;
f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);
t_nufft = toc;
fprintf("NUFFT time: %.1e\n", t_nufft);

tic;
f = A.Apply(u_ex);
t_apply = toc;
fprintf("Apply time: %.1e\n", t_apply);

df = f - f_nufft;
rel_err = norm(df) / norm(f_nufft);
fprintf("Rel err: %.1e\n", rel_err);

%% URV Factorization.
tic;
A.URV_Factor();
t_factor = toc;
fprintf("Factor time: %.1e\n", t_factor);

%% MRI Reconstruction.
P = phantom('Modified Shepp-Logan', n);
u_ex = reshape(P, N, []);
f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);

tic;
u_solve = A.URV_Solve(f_nufft);
t_solve = toc;
fprintf("Solve time: %.1e\n", t_solve);
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
