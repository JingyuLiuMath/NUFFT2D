%% Basic Settings.
clear;
close all;

p = 10;

%% Generate points x and omega.
N = 2^p;
M = 2 * N;
fprintf("M: %d, N: %d\n", M, N);

% Random
omega = rand(N, 1) * N - 0.5;

%% Apply.
u_ex = randn(N, 1) + randn(N, 1) * 1i;
A = NUDFT1_Matrix(omega, M);

tic;
f_ex = A * u_ex;
t_ex = toc;

tic;
f = MY_NUFFT1(u_ex, omega, M, "matlab");
t_nufft = toc;

df = f - f_ex;
rel_apply_err_norm = norm(df) / norm(f_ex);
fprintf("Relative application error: %.4e\n", rel_apply_err_norm);
fprintf("Exact comp time: %.4e\n", t_ex);
fprintf("NUFFT comp time: %.4e\n", t_nufft);
