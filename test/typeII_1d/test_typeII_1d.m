%% Basic Settings.
clear;
close all;
warning off;

p = 10;

min_points = 128;
tol = 1e-6;

tol_cg = 1e-12;
maxit_cg = 50;

fprintf("min_points: %d\n", min_points);
fprintf("tol: %.1e\n", tol);

%% Generate points x and omega.
N = 2^p;
M = 2 * N;
fprintf("M: %d, N: %d\n", M, N);

% Random
x = sort(rand(M, 1));

%% NUDFT2.
run_INUDFT2(x, N, min_points, tol, tol_cg, maxit_cg);
