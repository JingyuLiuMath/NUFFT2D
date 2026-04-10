clear;
close all;
warning off;

p_list = (5 : 9)';
num_n = length(p_list);

beta_list = [0.5, 0.6, 0.7]';
num_beta = length(beta_list);

tol_cg = 1e-12;
maxit = 500;

fprintf("tol_cg: %.1e\n", tol_cg);
fprintf("maxit: %d\n", maxit);

for it_beta = 1 : num_beta
    for it_p = 1 : num_n
        beta = beta_list(it_beta);
        p = p_list(it_p);

        n = 2^p;
        nx = n;
        ny = n;
        N = nx * ny;
        x = PolarGrid(n, 1, beta * p);
        M = size(x, 1);

        fprintf("\n\n\n\n");
        fprintf("p: %d\n", p);
        fprintf("beta: %.1e\n", beta);
        fprintf("M: %d, N: %d\n", M, N);
        fprintf("M/N: %.1e\n", M / N);

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

        save("./data/typeII_2d_points_" + string(p) + "_" + string(beta) + ".mat", "x");
        
        save("./data/typeII_2d_results_" + string(p) + "_" + string(beta) + ".mat", ...
            "M", "N", ...
            ...
            "P", ...
            ...
            "t_iter", "iter", "rel_res", "rel_acc", "P_reconstruct", ...
            ...
            "flag", "resvec");
    end
end