function result = run_INUDFT2_2D(xy, n, ...
    min_points, tol_hss, ...
    tol_cg, maxit_cg)

M = size(xy, 1);
nx = n;
ny = n;
N = nx * ny;
fprintf("Basic info.\n");
fprintf("  M: %d\n", M);
fprintf("  n: %d, N: %d\n", n, N);
fprintf("  M/N: %.1e\n", M / N);

fprintf("  min_points: %d\n", min_points);
fprintf("  tol_hss: %.1e\n", tol_hss);

fprintf("  tol_cg: %.1e\n", tol_cg);
fprintf("  maxit_cg: %d\n", maxit_cg);

result = struct();
result.M = M;
result.n = n;
result.N = N;
result.min_points = min_points;
result.tol_hss = tol_hss;
result.tol_cg = tol_cg;
result.maxit_cg = maxit_cg;

P = phantom('Modified Shepp-Logan', n);
c_ex = reshape(P, N, []);
result.c_ex = c_ex;
f_ex = MY_NUFFT2_2D(c_ex, xy, nx, ny);

% NUDFT2.
fprintf("NUDFT2.\n")
tic;
A = NUDFT2_2D(nx, ny);
A.Construct_ID_Proxy(xy, min_points, tol_hss);
result.t_construct = toc;
fprintf("  t_construct: %e\n", result.t_construct);

result.hss_rank = A.Rank();
fprintf("  HSS rank: %d\n", result.hss_rank);

mem_exact = M * N;
result.mem = A.Storage();
mem_ratio = result.mem / mem_exact;
fprintf("  Mem ratio: %.1e\n", mem_ratio);

f = A.Apply(c_ex);
result.rel_err = norm(f - f_ex) / norm(f_ex);
fprintf("  rel_err: %.1e\n", result.rel_err);

tic;
A.URV_Factor();
result.t_factor = toc;
fprintf("  t_factor: %e\n", result.t_factor);

% Direct solution.
fprintf("Direct solution.\n");
tic;
c_direct = A.URV_Solve(f_ex);
result.t_direct = toc;
c_direct = real(c_direct);
result.c_direct = c_direct;

fprintf("  t_direct: %.1e\n", result.t_direct);
result.rel_res_direct = norm(f_ex - MY_NUFFT2_2D(c_direct, xy, nx, ny)) / norm(f_ex);
fprintf("  rel_res_direct: %.1e\n", result.rel_res_direct);
result.rel_err_direct = norm(c_ex - c_direct) / norm(c_ex);
fprintf("  rel_err_direct: %.1e\n", result.rel_err_direct);

% Iterative solution.
fprintf("Iterative solution.\n");

% fprintf("  Without precond.\n")
% tic;
% [c_cg, result.flag_cg, ~, result.iter_cg, ~]  = INUDFT2_2D_CG(xy, nx, ny, f_ex, tol_cg, maxit_cg);
% result.t_cg = toc;
% c_cg = real(c_cg);
% result.c_cg = c_cg;
% 
% fprintf("    t_cg: %.1e\n", result.t_cg);
% fprintf("    iter_cg: %d\n", result.iter_cg);
% result.rel_res_cg = norm(f_ex - MY_NUFFT2_2D(c_cg, xy, nx, ny)) / norm(f_ex);
% fprintf("    rel_res_cg: %.1e\n", result.rel_res_cg);
% result.rel_err_cg = norm(c_ex - c_cg) / norm(c_ex);
% fprintf("    rel_err_cg: %.1e\n", result.rel_err_cg);

fprintf("  With precond.\n");
tic;
[c_pcg, result.flag_pcg, ~, result.iter_pcg] = INUDFT2_2D_PCG(A, xy, nx, ny, f_ex, tol_cg, maxit_cg);
result.t_pcg = toc;
c_pcg = real(c_pcg);
result.c_pcg = c_pcg;

fprintf("    t_pcg: %.1e\n", result.t_pcg);
fprintf("    iter_pcg: %d\n", result.iter_pcg);
result.rel_res_pcg = norm(f_ex - MY_NUFFT2_2D(c_pcg, xy, nx, ny)) / norm(f_ex);
fprintf("    rel_res_pcg: %.1e\n", result.rel_res_pcg);
result.rel_err_pcg = norm(c_ex - c_pcg) / norm(c_ex);
fprintf("    rel_err_pcg: %.1e\n", result.rel_err_pcg);


end