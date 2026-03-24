%% Basic Settings.
clear;
close all;

p = 10;

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
x = rand(M, 1);

%% Apply.
u_ex = randn(N, 1) + randn(N, 1) * 1i;
A = NUDFT2_Matrix(x, N);

tic;
f_ex = A * u_ex;
t_ex = toc;

tic;
f = MY_NUFFT2(u_ex, x, N);
t_nufft = toc;

df = f - f_ex;
rel_apply_err_norm = norm(df) / norm(f_ex);
fprintf("Relative application error: %.4e\n", rel_apply_err_norm);
fprintf("Exact comp time: %.4e\n", t_ex);
fprintf("NUFFT comp time: %.4e\n", t_nufft);
