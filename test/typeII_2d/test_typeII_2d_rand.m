%% Basic Settings.
clear;
close all;
warning off;

p = 5;

min_points = 256;
tol_hss = 1e-5;

tol_cg = 1e-12;
maxit_cg = 20;

%% Generate points x and omega.
n = 2^p;
nx = n;
ny = n;
N = nx * ny;
M = 2 * N;

xy = rand(M, 2);

fprintf("M: %d, N: %d\n", M, N);

figure();
plot(xy(:, 1), xy(:, 2), "Marker", ".", "LineStyle", "none");
axis equal;
xlim([0, 1]);
ylim([0, 1]);

%% NUDFT2_Matrix.
% A = NUDFT2_2D_Matrix(x, nx, ny);
% kappa_A = cond(A);
% fprintf("cond(A): %.1e\n", kappa_A);

%% NUDFT2.
P = phantom('Modified Shepp-Logan', n);
c_ex = reshape(P, N, []);
f_ex = MY_NUFFT2_2D(c_ex, xy, nx, ny);

run_INUDFT2_2D(...
    xy, n, ...
    min_points, tol_hss, ...
    tol_cg, maxit_cg, ...
    c_ex, f_ex);
