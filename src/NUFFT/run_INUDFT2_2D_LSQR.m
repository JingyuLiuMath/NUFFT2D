function result = run_INUDFT2_2D_LSQR(xy, n, ...
    tol_cg, maxit_cg, ...
    c_ex, f_ex)

M = size(xy, 1);
nx = n;
ny = n;
N = nx * ny;
fprintf("Basic info.\n");
fprintf("  M: %d\n", M);
fprintf("  n: %d, N: %d\n", n, N);
fprintf("  M/N: %.1e\n", M / N);

fprintf("  tol_cg: %.1e\n", tol_cg);
fprintf("  maxit_cg: %d\n", maxit_cg);

result = struct();
result.M = M;
result.n = n;
result.N = N;
result.tol_cg = tol_cg;
result.maxit_cg = maxit_cg;

% Iterative solution.
fprintf("Iterative solution.\n");

fprintf("  Without precond.\n")
tic;
[c_cg, result.flag_cg, ~, result.iter_cg, ~]  = INUDFT2_2D_CG(xy, nx, ny, f_ex, tol_cg, maxit_cg);
result.t_cg = toc;
c_cg = real(c_cg);
result.c_cg = c_cg;

fprintf("    t_cg: %.1e\n", result.t_cg);
fprintf("    iter_cg: %d\n", result.iter_cg);
result.rel_res_cg = norm(f_ex - MY_NUFFT2_2D(c_cg, xy, nx, ny)) / norm(f_ex);
fprintf("    rel_res_cg: %.1e\n", result.rel_res_cg);
result.rel_err_cg = norm(c_ex - c_cg) / norm(c_ex);
fprintf("    rel_err_cg: %.1e\n", result.rel_err_cg);

end