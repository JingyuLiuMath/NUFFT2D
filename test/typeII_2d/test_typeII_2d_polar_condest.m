%% Basic Settings.
clear;
close all;
warning off;

p = 5;

tol_hss = 1e-4;
tol_pcg = 1e-12;
maxit_pcg = 100;
maxit_condest = 8;
restarts = 2;

fprintf("tol_hss: %.1e\n", tol_hss);
fprintf("tol_pcg: %.1e\n", tol_pcg);
fprintf("maxit_pcg: %d\n", maxit_pcg);
fprintf("maxit_condest: %d\n", maxit_condest);
fprintf("restarts: %d\n", restarts);

%% Generate points.
n = 2^p;
nx = n;
ny = n;
N = nx * ny;

x = PolarGrid(n);
M = size(x, 1);

fprintf("M: %d, N: %d\n", M, N);

%% Exact matrix and exact condition number target.
A_exact = NUDFT2_2D_Matrix(x, nx, ny);
ATA_exact = A_exact' * A_exact;
kappa_exact = sqrt(cond(ATA_exact, 1));
kappa_exact_2norm = cond(A_exact);
kappa_condest_1norm = sqrt(condest(ATA_exact));
fprintf("Exact sqrt(cond_1(A^H A)): %.6e\n", kappa_exact);
fprintf("Exact cond_2(A): %.6e\n", kappa_exact_2norm);
fprintf("MATLAB condest sqrt(condest(A^H A)): %.6e\n", kappa_condest_1norm);

%% HSS object.
min_points = 256;
fprintf("min_points: %d\n", min_points);

tic;
A = NUDFT2_2D(nx, ny);
A.Construct_ID_Proxy(x, min_points, tol_hss);
t_construct = toc;
fprintf("Construct time: %.1e\n", t_construct);

tic;
A.URV_Factor();
t_factor = toc;
fprintf("Factor time: %.1e\n", t_factor);

%% Condition number estimate.
tic;
[kappa_est, info] = Condest_NUDFT2_2D(A, x, nx, ny, maxit_pcg, tol_pcg, maxit_condest, restarts);
t_condest = toc;
fprintf("Condest time: %.1e\n", t_condest);

fprintf("Estimated sqrt(cond_1(A^H A)): %.6e\n", kappa_est);
fprintf("Exact sqrt(cond_1(A^H A)): %.6e\n", kappa_exact);
fprintf("Exact cond_2(A): %.6e\n", kappa_exact_2norm);
fprintf("MATLAB condest sqrt(condest(A^H A)): %.6e\n", kappa_condest_1norm);
fprintf("Relative error: %.6e\n", abs(kappa_est - kappa_exact) / kappa_exact);
