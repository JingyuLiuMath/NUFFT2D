%% Basic Settings.
clear;
close all;

p = 5;

%% Generate points x and omega.
mx = 2^p;
my = 2^(p + 1);
M = mx * my;
n = 2^p;
N = M / 2;
fprintf("M: %d, N: %d\n", M, M);

% Random.
omega = rand(N, 2) * n - 0.5;

%% Apply.
u_ex = randn(N, 1) + randn(N, 1) * 1i;
A = NUDFT1_2D_Matrix(omega, mx, my);

tic;
f_ex = A * u_ex;
t_ex = toc;

tic;
f = MY_NUFFT1_2D(u_ex, omega, mx, my);
t_nufft = toc;

df = f - f_ex;
rel_apply_err_norm = norm(df) / norm(f_ex);
fprintf("Relative application error: %.4e\n", rel_apply_err_norm);
fprintf("Exact comp time: %.4e\n", t_ex);
fprintf("NUFFT comp time: %.4e\n", t_nufft);
