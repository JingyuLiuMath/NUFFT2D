clear;
close all;
warning off;

p_list = (5 : 9)';
num_n = length(p_list);

tol_hss = 1e-2;
% tol_hss = 1e-4;
tol_cg = 1e-12;
maxit = 500;

fprintf("tol_hss: %.1e\n", tol_hss);
fprintf("tol_cg: %.1e\n", tol_cg);
fprintf("maxit: %d\n", maxit);

for it = 1 : num_n
    p = p_list(it);
    n = 2^p;
    N = n^2;
    M = 2 * N;

    alpha = 6;
    x = PolarGrid(n, 1, alpha);
    M = size(x, 1);

    fprintf("\n\n\n\n");
    fprintf("alpha: %d\n", alpha);
    fprintf("M: %d, N: %d\n", M, N);

    nx = n;
    ny = n;

    min_points = n * p;
    fprintf("min_points: %d\n", min_points);

    % NUDFT2.
    tic;
    A = NUDFT2_2D(nx, ny);
    A.Construct_ID_Proxy(x, min_points, tol_hss);
    t_construct = toc;
    fprintf("Construct time: %.1e\n", t_construct);

    hss_rank = A.Rank();
    fprintf("HSS rank: %d\n", hss_rank);

    mem_exact = M * N;
    mem = A.Storage();
    fprintf("Memory (GB): %.1e\n", double_to_gb(mem));
    mem_ratio = mem / mem_exact;
    fprintf("Mem ratio: %.1e\n", mem_ratio);

    % URV Factorization.
    tic;
    A.URV_Factor();
    t_factor = toc;
    fprintf("Factor time: %.1e\n", t_factor);

    % MRI Reconstruction.
    P = phantom('Modified Shepp-Logan', n);
    u_ex = reshape(P, N, []);
    f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);

    fprintf("\nDirect Solver\n");
    tic;
    u_direct = A.URV_Solve(f_nufft);
    t_direct = toc;
    fprintf("Time: %.1e\n", t_direct);
    u_direct = real(u_direct);
    
    r_direct = f_nufft - MY_NUFFT2_2D(u_direct, x, nx, ny);
    rel_res_direct = norm(r_direct) / norm(f_nufft);
    fprintf("Rel res: %.1e\n", rel_res_direct);

    e_direct = u_ex - u_direct;
    rel_acc_direct = norm(e_direct) / norm(u_ex);
    fprintf("Rel acc: %.1e\n", rel_acc_direct);

    P_reconstruct_direct = reshape(u_direct, n, n);

    fprintf("\nIterative Solver\n");
    tic;
    [u_iter, flag, relres, iter, resvec] = INUDFT2_2D_PCG(A, x, nx, ny, f_nufft, tol_cg, maxit);
    t_iter = toc;
    fprintf("Time: %.1e\n", t_iter);
    fprintf("Iter number: %d\n", iter);
    u_iter = real(u_iter);
    
    r_iter = f_nufft - MY_NUFFT2_2D(u_iter, x, nx, ny);
    rel_res_iter = norm(r_iter) / norm(f_nufft);
    fprintf("Rel res: %.1e\n", rel_res_iter);
    
    e_iter = u_ex - u_iter;
    rel_acc_iter = norm(e_iter) / norm(u_ex);
    fprintf("Rel acc: %.1e\n", rel_acc_iter);
    
    P_reconstruct_iter = reshape(u_iter, n, n);

    t_pre = t_construct + t_factor;
    t_total = t_pre + t_iter;
    fprintf("Total time: %.1e\n", t_total);

    save("./data/typeII_2d_results_" + string(p) + "_tol_" + string(tol_hss) + ".mat", ...
        "M", "N", ...
        ...
        "t_construct", ...
        "hss_rank", ...
        "mem", "mem_exact", "mem_ratio", ...
        ...
        "t_factor", ...
        ...
        "P", ...
        ...
        "t_direct", "rel_res_direct", "rel_acc_direct", "P_reconstruct_direct", ...
        ...
        "t_iter", "iter", "rel_res_iter", "rel_acc_iter", "P_reconstruct_iter", ...
        ...
        "t_pre", "t_total");
end