function result = run_INUDFT2_2D_Condest(xy, n, ...
    min_points, tol_hss, ...
    tol_cg, maxit_cg, ...
    maxit_condest, restarts)

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

fprintf("  maxit_condest: %.1e\n", maxit_condest);
fprintf("  restarts: %d\n", restarts);

result = struct();
result.M = M;
result.n = n;
result.N = N;
result.min_points = min_points;
result.tol_hss = tol_hss;
result.maxit_condest = maxit_condest;
result.restarts = restarts;

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

tic;
A.URV_Factor();
result.t_factor = toc;
fprintf("  t_factor: %e\n", result.t_factor);

% Condition number estimate.
fprintf("Condition number estimate.\n")

tic;
[result.kappa1_est, ~] = Condest_NUDFT2_2D(A, xy, nx, ny, maxit_cg, tol_cg, maxit_condest, restarts);
t_cond_est = toc;
fprintf("  kappa1_est: %.1e\n", result.kappa1_est);
fprintf("  Estimate time: %.1e\n", t_cond_est);

if n <= 2^7
    tic;
    A_exact = NUDFT2_2D_Matrix(xy, nx, ny);
    ATA_exact = A_exact' * A_exact;
    result.kappa1_ex = sqrt(cond(ATA_exact, 1));
    result.kappa2_ex = cond(A_exact);
    t_cond_exact = toc;
    fprintf("  kappa1_ex: %.1e\n", result.kappa1_ex);
    fprintf("  kappa2_ex: %.1e\n", result.kappa2_ex);
    fprintf("  Exact time: %.1e\n", t_cond_exact);
else
    fprintf("Skip explicit matrix cond.\n");
end

end