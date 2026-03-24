%% Basic Settings.
clear;
close all;
warning off;

p = 5;

min_points = 64;
rank_or_tol = 400;
fprintf("min_points: %d\n", min_points);
if rank_or_tol >= 1
    fprintf("rank: %d\n", rank_or_tol);
else
    fprintf("tol: %.4e\n", rank_or_tol);
end

%% Generate points x and omega.
nx = 2^p;
ny = 2^p;
N = nx * ny;
M = 8 * N;
fprintf("M: %d, N: %d\n", M, N);

% Perturbation.
% perturb_size = 0.2;
% x = (0 : (M - 1))' / M;
% x = x + (rand(M, 1) * 2 - 1) * perturb_size / M;
% x = sort(x);

% Random
x = rand(M, 2);

%% NUDFT2.
tic;
A = NUDFT2_2D(nx, ny);
A.Construct_BlackBox(x, min_points, rank_or_tol);
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
f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);
t_nufft = toc;
fprintf("NUFFT time: %.4e\n", t_nufft);

tic;
f = A.Apply(u_ex);
t_apply = toc;
fprintf("Apply time: %.4e\n", t_apply);

df = f - f_nufft;
rel_err = norm(df) / norm(f_nufft);
fprintf("Rel err: %.4e\n", rel_err);
