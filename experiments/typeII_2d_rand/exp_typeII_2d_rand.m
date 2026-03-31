clear;
close all;
warning off;

p_list = 5 : 9;
% p_list = 3 : 5;
num_n = length(p_list);

tol_hss = 1e-4;

fprintf("tol_hss: %.1e\n", tol_hss);

for it = 1 : num_n
    p = p_list(it);
    n = 2^p;
    N = n^2;
    M = 3 * N;

    x = rand(M, 2);
    
    fprintf("\n\n");
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

    % Apply.
    u_ex = randn(N, 1) + randn(N, 1) * 1i;

    tic;
    f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);
    t_nufft = toc;
    fprintf("NUFFT time: %.1e\n", t_nufft);

    tic;
    f = A.Apply(u_ex);
    t_apply = toc;
    fprintf("Apply time: %.1e\n", t_apply);

    df = f - f_nufft;
    rel_err = norm(df) / norm(f_nufft);
    fprintf("Rel err: %.1e\n", rel_err);

    % URV Factorization.
    tic;
    A.URV_Factor();
    t_factor = toc;
    fprintf("Factor time: %.1e\n", t_factor);

    % MRI Reconstruction.
    P = phantom('Modified Shepp-Logan', n);
    u_ex = reshape(P, N, []);
    f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);

    tic;
    u_solve = A.URV_Solve(f_nufft);
    t_solve = toc;
    fprintf("Solve time: %.1e\n", t_solve);
    u_solve = real(u_solve);

    r = f_nufft - MY_NUFFT2_2D(u_solve, x, nx, ny);
    rel_res = norm(r) / norm(f_nufft);
    fprintf("Rel res: %.1e\n", rel_res);

    e = u_ex - u_solve;
    rel_acc = norm(e) / norm(u_ex);
    fprintf("Rel acc: %.1e\n", rel_acc);

    P_reconstruct = reshape(u_solve, n, n);

    t_total = t_construct + t_factor + t_solve;
    fprintf("Total time: %.1e\n", t_total);

    save("./data/typeII_2d_results_" + string(p) + ".mat", ...
        "M", "N", ...
        ...
        "t_construct", ...
        "hss_rank", ...
        "mem", "mem_exact", "mem_ratio", ...
        ...
        "t_apply", "t_nufft", ...
        "rel_err", ...
        ...
        "t_factor", ...
        "t_solve", ...
        "rel_res", ...
        "rel_acc", ...
        "P", "P_reconstruct", ...
        ...
        "t_total");
end