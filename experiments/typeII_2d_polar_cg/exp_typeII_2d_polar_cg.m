clear;
close all;
warning off;

% p_list = 5 : 9;
p_list = 3 : 5;
num_n = length(p_list);

tol = 1e-8;
% maxit = 500;
maxit = 100;
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

    fprintf("RHS (approximately) in range(A)\n");
    u_ex = randn(N, 1) + randn(N, 1) * 1i;
    f_nufft = MY_NUFFT2_2D(u_ex, x, nx, ny);
    
    tic;
    [u_solve, flag, relres, iter, resvec] = INUDFT2_2D_CG(x, nx, ny, f_nufft, tol, maxit);
    t_iter = toc;
    fprintf("Solve time: %.4e\n", t_iter);
    fprintf("Iter number: %d\n", iter);
    
    r = f_nufft - MY_NUFFT2_2D(u_solve, x, nx, ny);
    rel_res = norm(r) / norm(f_nufft);
    fprintf("Rel res: %e\n", rel_res);

    save("./data/typeII_2d_results_" + string(p) + ".mat", ...
        "t_iter", ...
        "rel_res", ...
        "iter");
end