% Experiment for condition-number estimation of typeII 2D polar setting.
%
% Rule:
%   p <= 7      : explicit matrix condition number + Condest_NUDFT2_2D.
%   8 <= p <= 9 : Condest_NUDFT2_2D only.

clear;
p_list = (5 : 9)';
num_n = length(p_list);

close all;
warning off;

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

for it = 1 : num_n
    p = p_list(it);
    n = 2^p;
    nx = n;
    ny = n;
    N = n^2;

    alpha = 5;
    x = PolarGrid(n, 1, alpha);
    M = size(x, 1);

    fprintf("\n\n\n\n");
    fprintf("p: %d, alpha: %d\n", p, alpha);
    fprintf("M: %d, N: %d\n", M, N);

    min_points = n * p;
    fprintf("min_points: %d\n", min_points);

    % Build HSS-based NUDFT object.
    tic;
    A = NUDFT2_2D(nx, ny);
    A.Construct_ID_Proxy(x, min_points, tol_hss);
    t_construct = toc;
    fprintf("Construct time: %.1e\n", t_construct);

    hss_rank = A.Rank();
    fprintf("HSS rank: %d\n", hss_rank);

    mem_exact = M * N;
    mem = A.Storage();
    mem_ratio = mem / mem_exact;
    fprintf("Memory (GB): %.1e\n", double_to_gb(mem));
    fprintf("Mem ratio: %.1e\n", mem_ratio);

    tic;
    A.URV_Factor();
    t_factor = toc;
    fprintf("Factor time: %.1e\n", t_factor);

    % Condest_NUDFT2_2D always computed for p in [5, 9].
    tic;
    [cond_number_est, info_condest] = Condest_NUDFT2_2D(A, x, nx, ny, ...
        maxit_pcg, tol_pcg, maxit_condest, restarts);
    t_condest = toc;
    fprintf("Estimated cond number: %.1e\n", cond_number_est);
    fprintf("Condest time: %.1e\n", t_condest);

    % Explicit condition number when p <= 7.
    cond_number_exact = NaN;
    cond_number_exact_2norm = NaN;
    t_exact = NaN;
    if p <= 7
        tic;
        A_exact = NUDFT2_2D_Matrix(x, nx, ny);
        ATA_exact = A_exact' * A_exact;
        cond_number_exact = sqrt(cond(ATA_exact, 1));
        cond_number_exact_2norm = cond(A_exact);
        t_exact = toc;

        fprintf("Exact sqrt(cond_1(A^H A)): %.1e\n", cond_number_exact);
        fprintf("Exact cond_2(A): %.1e\n", cond_number_exact_2norm);
        fprintf("Exact computation time: %.1e\n", t_exact);
    else
        fprintf("Skip explicit matrix cond for p=%d (rule: only condest for 8<=p<=9).\n", p);
    end

    save("./data/typeII_2d_polar_plsqr_bad_condest_" + string(p) + ...
        "_tol_" + string(tol_hss) + ".mat", ...
        "p", "n", "alpha", "M", "N", "nx", "ny", "x", ...
        "tol_hss", "tol_pcg", "maxit_pcg", "maxit_condest", "restarts", ...
        "min_points", ...
        "t_construct", "hss_rank", "mem", "mem_exact", "mem_ratio", ...
        "t_factor", ...
        "cond_number_est", "info_condest", "t_condest", ...
        "cond_number_exact", "cond_number_exact_2norm", "t_exact");
end
