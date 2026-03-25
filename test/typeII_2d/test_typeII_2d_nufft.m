%% Basic Settings.
clear;
close all;

p = 5;

%% Generate points x and omega.
nx = 2^p;
ny = 2^p;
N = nx * ny;
M = 2 * N;
fprintf("M: %d, N: %d\n", M, N);

% Random
xy = rand(M, 2);

%% Apply.
u_ex = randn(N, 1) + randn(N, 1) * 1i;
A = NUDFT2_2D_Matrix(xy, nx, ny);

tic;
f_ex = A * u_ex;
t_ex = toc;

tic;
f = MY_NUFFT2_2D(u_ex, xy, nx, ny);
t_nufft = toc;

df = f - f_ex;
rel_apply_err_norm = norm(df) / norm(f_ex);
fprintf("Relative application error: %.4e\n", rel_apply_err_norm);
fprintf("Exact comp time: %.4e\n", t_ex);
fprintf("NUFFT comp time: %.4e\n", t_nufft);
