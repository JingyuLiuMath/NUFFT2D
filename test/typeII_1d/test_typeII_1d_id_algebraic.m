%% Basic Settings.
clear;
close all;
warning off;

p = 10;

min_points = 128;
rank_or_tol = 1e-6;

fprintf("min_points: %d\n", min_points);
if rank_or_tol >= 1
    fprintf("rank: %d\n", rank_or_tol);
else
    fprintf("tol: %.4e\n", rank_or_tol);
end

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
x = sort(rand(M, 1));

%% NUDFT2.
tic;
A = NUDFT2(N);
A.Construct_ID_Algebraic(x, min_points, rank_or_tol);
t_construct = toc;
fprintf("Construct time: %.4e\n", t_construct);

r = A.Rank();
fprintf("HSS rank: %d\n", r);

mem_exact = M * N;
mem = A.Storage();
mem_ratio = mem / mem_exact;
fprintf("Mem ratio: %.4e\n", mem_ratio);

%% Apply.
u_ex = randn(N, 1) + randn(N, 1) * 1i;

tic;
f_nufft = MY_NUFFT2(u_ex, x, N);
t_nufft = toc;
fprintf("NUFFT time: %.4e\n", t_nufft);

tic;
f = A.Apply(u_ex);
t_apply = toc;
fprintf("Apply time: %.4e\n", t_apply);

df = f - f_nufft;
rel_err = norm(df) / norm(f_nufft);
fprintf("Rel err: %.4e\n", rel_err);

%% URV Factorization.
tic;
A.URV_Factor();
t_factor = toc;
fprintf("Factor time: %.4e\n", t_factor);

%% Solution: RHS (approximately) in range(A).
fprintf("RHS (approximately) in range(A)\n");
u_ex = randn(N, 1) + randn(N, 1) * 1i;
f_nufft = MY_NUFFT2(u_ex, x, N);

tic;
u_solve = A.URV_Solve(f_nufft);
t_solve = toc;
fprintf("Solve time: %.4e\n", t_solve);

r = f_nufft - MY_NUFFT2(u_solve, x, N);
rel_res = norm(r) / norm(f_nufft);
fprintf("Rel res: %e\n", rel_res);
