clear;
close all;
warning off;

p_list = 5 : 8;
% p_list = 3 : 5;
num_n = length(p_list);

tol_hss = 1e-4;
tol_cg = 1e-12;
maxit = 500;
% maxit = 100;

fprintf("tol_hss: %e\n", tol_hss);
fprintf("tol_cg: %.4e\n", tol_cg);
fprintf("maxit: %d\n", maxit);

for it = 1 : num_n
    p = p_list(it);
    n = 2^p;
    N = n^2;

    x = PolarGrid(n);
    M = size(x, 1);

    fprintf("\n\n");
    fprintf("M: %d, N: %d\n", M, N);

    nx = n;
    ny = n;

    min_points = n * p;
    fprintf("min_points: %d\n", min_points);

    % NUDFT2.
    tic;
    A = NUDFT2_2D(nx, ny);
    A.Construct_ID_Algebraic(x, min_points, tol_hss);
    t_construct = toc;
    fprintf("Construct time: %.4e\n", t_construct);

    hss_rank = A.Rank();
    fprintf("HSS rank: %d\n", hss_rank);

    mem_exact = M * N;
    mem = A.Storage();
    mem_ratio = mem / mem_exact;
    fprintf("Mem ratio: %.4e\n", mem_ratio);

    % URV Factorization.
    tic;
    A.URV_Factor();
    t_factor = toc;
    fprintf("Factor time: %.4e\n", t_factor);

    % MRI Reconstruction.
    P = phantom('Modified Shepp-Logan', n);
    u_ex = reshape(P, N, []);
    f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);

    fprintf("Direct Solver\n");
    u_solve = A.URV_Solve(f_nufft);
    u_solve = real(u_solve);
    r = f_nufft - MY_NUFFT2_2D(u_solve, x, nx, ny);
    rel_res = norm(r) / norm(f_nufft);
    fprintf("Rel res: %e\n", rel_res);

    e = u_ex - u_solve;
    rel_acc = norm(e) / norm(u_ex);
    fprintf("Rel acc: %.4e\n", rel_acc);

    fprintf("Iterative Solver\n");
    tic;
    [u_solve, flag, relres, iter, resvec] = INUDFT2_2D_PCG(A, x, nx, ny, f_nufft, tol_cg, maxit);
    t_iter = toc;
    fprintf("Iterative solve time: %.4e\n", t_iter);
    fprintf("Iter number: %d\n", iter);
    u_solve = real(u_solve);

    r = f_nufft - MY_NUFFT2_2D(u_solve, x, nx, ny);
    rel_res = norm(r) / norm(f_nufft);
    fprintf("Rel res: %e\n", rel_res);

    e = u_ex - u_solve;
    rel_acc = norm(e) / norm(u_ex);
    fprintf("Rel acc: %.4e\n", rel_acc);

    P_reconstruct = reshape(u_solve, n, n);

    t_pre = t_construct + t_factor;
    t_total = t_pre + t_iter;
    fprintf("Total time: %e\n", t_total);

    save("./data/typeII_2d_results_" + string(p) + ".mat", ...
        "t_construct", ...
        "hss_rank", ...
        "mem", "mem_exact", "mem_ratio", ...
        ...
        "t_factor", ...
        "t_iter", ...
        "iter", ...
        "rel_res", ...
        "rel_acc", ...
        "P", "P_reconstruct", ...
        ...
        "t_pre", ...
        "t_total");
end