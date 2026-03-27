clear;
close all;
warning off;

p_list = 5 : 8;
% p_list = 3 : 5;
num_n = length(p_list);

tol = 1e-12;
maxit = 500;
% maxit = 100;
fprintf("tol: %.4e\n", tol);
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

    % MRI Reconstruction.
    P = phantom('Modified Shepp-Logan', n);
    u_ex = reshape(P, N, []);
    f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);

    tic;
    [u_solve, flag, relres, iter, resvec] = INUDFT2_2D_CG(x, nx, ny, f_nufft, tol, maxit);
    t_iter = toc;
    fprintf("Solve time: %.4e\n", t_iter);
    fprintf("Iter number: %d\n", iter);
    u_solve = real(u_solve);

    r = f_nufft - MY_NUFFT2_2D(u_solve, x, nx, ny);
    rel_res = norm(r) / norm(f_nufft);
    fprintf("Rel res: %e\n", rel_res);

    e = u_ex - u_solve;
    rel_acc = norm(e) / norm(u_ex);
    fprintf("Rel acc: %.4e\n", rel_acc);

    P_reconstruct = reshape(u_solve, n, n);

    save("./data/typeII_2d_results_" + string(p) + ".mat", ...        ...
        "t_iter", ...
        "iter", ...
        "rel_res", ...
        "rel_acc", ...
        "P", "P_reconstruct");
end