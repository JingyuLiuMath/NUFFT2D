%% Basic Settings.
clear;
close all;
warning off;

p = 10;

min_points = 128;
tol_hss = 1e-6;

tol_cg = 1e-12;
maxit_cg = 50;

%% Generate points x and omega.
N = 2^p;
M = 2 * N;
fprintf("M: %d, N: %d\n", M, N);

% Random
x = sort(rand(M, 1));

%% NUDFT2.
c_ex = randn(N, 1) + 1i * randn(N, 1);
f_ex = MY_NUFFT2(c_ex, x, N);

run_INUDFT2(...
    x, N, ...
    min_points, tol_hss, ...
    tol_cg, maxit_cg, ...
    c_ex, f_ex);
