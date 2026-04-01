clear;
close all;
warning off;

p_list = (5 : 9)';
num_n = length(p_list);

tol_cg = 1e-12;
maxit = 500;

fprintf("tol_cg: %.1e\n", tol_cg);
fprintf("maxit: %d\n", maxit);

for it = 1 : num_n
    p = p_list(it);
    n = 2^p;
    N = n^2;
    M = 2 * N;

    x = rand(M, 2);

    fprintf("\n\n\n\n");
    fprintf("M: %d, N: %d\n", M, N);
    
    nx = n;
    ny = n;

    min_points = n * p;
    fprintf("min_points: %d\n", min_points);

    % MRI Reconstruction.
    P = phantom('Modified Shepp-Logan', n);
    u_ex = reshape(P, N, []);
    f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);

    fprintf("\nIterative Solver\n");
    tic;
    [u_solve, flag, relres, iter, resvec] = INUDFT2_2D_CG(x, nx, ny, f_nufft, tol_cg, maxit);
    t_iter = toc;
    fprintf("Time: %.1e\n", t_iter);
    fprintf("Iter number: %d\n", iter);
    u_solve = real(u_solve);

    r = f_nufft - MY_NUFFT2_2D(u_solve, x, nx, ny);
    rel_res = norm(r) / norm(f_nufft);
    fprintf("Rel res: %.1e\n", rel_res);

    e = u_ex - u_solve;
    rel_acc = norm(e) / norm(u_ex);
    fprintf("Rel acc: %.1e\n", rel_acc);

    P_reconstruct = reshape(u_solve, n, n);

    save("./data/typeII_2d_results_" + string(p) + "_alpha_" + string(alpha) + ".mat", ...
        "M", "N", ...
        ...
        "P", ...
        ...
        "t_iter", "iter", "rel_res", "rel_acc", "P_reconstruct");
end