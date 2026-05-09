function result = run_INUDFT2(x, N, ...
    min_points, tol, ...
    tol_cg, maxit_cg, ...
    c_ex, f_ex)

M = size(x, 1);

fprintf("Basic info.\n");
fprintf("  N: %d\n", N);
fprintf("  M: %d\n", M);

fprintf("  min_points: %d\n", min_points);
fprintf("  tol: %.1e\n", tol);

fprintf("  tol_cg: %.1e\n", tol_cg);
fprintf("  maxit_cg: %d\n", maxit_cg);

result = struct();
result.N = N;
result.M = M;
result.min_points = min_points;
result.tol = tol;
result.tol_cg = tol_cg;
result.maxit_cg = maxit_cg;

% NUDFT2.
fprintf("NUDFT2.\n")
tic;
A = NUDFT2(N);
A.Construct_fADI(x, min_points, tol);
result.t_construct = toc;
fprintf("  t_construct: %e\n", result.t_construct);

r = A.Rank();
fprintf("  HSS rank: %d\n", r);

mem_exact = M * N;
mem = A.Storage();
mem_ratio = mem / mem_exact;
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

fprintf("  t_direct: %.1e\n", result.t_direct);
result.rel_res_direct = norm(f_ex - MY_NUFFT2(c_direct, x, N)) / norm(f_ex);
fprintf("  rel_res_direct: %.1e\n", result.rel_res_direct);
result.rel_err_direct = norm(c_ex - c_direct) / norm(c_ex);
fprintf("  rel_err_direct: %.1e\n", result.rel_err_direct);

% Iterative solution.
fprintf("Iterative solution.\n");

% fprintf("  Without precond.\n")
% tic;
% [c_cg, result.flag_cg, ~, result.iter_cg, ~]  = INUDFT2_CG(x, N, f_ex, tol_cg, maxit_cg);
% result.t_cg = toc;
% 
% fprintf("    t_cg: %.1e\n", result.t_cg);
% fprintf("    iter_cg: %d\n", result.iter_cg);
% result.rel_res_cg = norm(f_ex - MY_NUFFT2(c_cg, x, N)) / norm(f_ex);
% fprintf("    rel_res_cg: %.1e\n", result.rel_res_cg);
% result.rel_err_cg = norm(c_ex - c_cg) / norm(c_ex);
% fprintf("    rel_err_cg: %.1e\n", result.rel_err_cg);

fprintf("  With precond.\n");
tic;
[c_pcg, result.flag_pcg, ~, result.iter_pcg] = INUDFT2_PCG(A, x, N, f_ex, tol_cg, maxit_cg);
result.t_pcg = toc;

fprintf("    t_pcg: %.1e\n", result.t_pcg);
fprintf("    iter_pcg: %d\n", result.iter_pcg);
result.rel_res_pcg = norm(f_ex - MY_NUFFT2(c_pcg, x, N)) / norm(f_ex);
fprintf("    rel_res_pcg: %.1e\n", result.rel_res_pcg);
result.rel_err_pcg = norm(c_ex - c_pcg) / norm(c_ex);
fprintf("    rel_err_pcg: %.1e\n", result.rel_err_pcg);

end