clear;
close all;
warning off;

p_list = (5 : 9)';
num_n = length(p_list);

% beta_list = [0.6, 0.7]';
beta_list = [0.7]';
num_beta = length(beta_list);

tol_hss_list = [1e-4]';
num_tol_hss = length(tol_hss_list);

tol_cg = 1e-12;
maxit = 500;
maxit_condest = 8;
restarts = 2;

fprintf("tol_cg: %.1e\n", tol_cg);
fprintf("maxit: %d\n", maxit);
fprintf("maxit_condest: %d\n", maxit_condest);
fprintf("restarts: %d\n", restarts);

warmup = 1;
for it_tol_hss = 1 : num_tol_hss
    for it_beta = 1 : num_beta
        for it_p = 1 : num_n
            tol_hss = tol_hss_list(it_tol_hss);
            beta = beta_list(it_beta);
            p = p_list(it_p);

            n = 2^p;
            nx = n;
            ny = n;
            N = nx * ny;
            
            load("../typeII_2d_polar_lsqr/data/typeII_2d_points_" + string(p) + "_" + string(beta) + ".mat");
            M = size(x, 1);

            fprintf("\n\n\n\n");
            fprintf("p: %d\n", p);
            fprintf("beta: %.1e\n", beta);
            fprintf("tol_hss: %.1e\n", tol_hss);
            fprintf("M: %d, N: %d\n", M, N);
            fprintf("M/N: %.1e\n", M / N);

            min_points = n * p;
            fprintf("min_points: %d\n", min_points);

            % NUDFT2.
            if warmup == 1
                A = NUDFT2_2D(nx, ny);
                A.Construct_ID_Proxy(x, min_points, tol_hss);
                A.URV_Factor();
                warmup = 0;
            end

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

            cond_number_est = NaN;
            cond_number_exact = NaN;
            cond_number_exact_2norm = NaN;
            t_cond_est = NaN;
            t_cond_exact = NaN;
            if tol_hss == 1e-4
                tic;
                [cond_number_est, info_condest] = Condest_NUDFT2_2D(A, x, nx, ny, maxit, tol_cg, maxit_condest, restarts);
                t_cond_est = toc;
                fprintf("Estimated cond number: %.1e\n", cond_number_est);
                fprintf("Condest time: %.1e\n", t_cond_est);

                if p <= 7
                    tic;
                    A_exact = NUDFT2_2D_Matrix(x, nx, ny);
                    ATA_exact = A_exact' * A_exact;
                    cond_number_exact = sqrt(cond(ATA_exact, 1));
                    cond_number_exact_2norm = cond(A_exact);
                    t_cond_exact = toc;
                    fprintf("Exact sqrt(cond_1(A^H A)): %.1e\n", cond_number_exact);
                    fprintf("Exact cond_2(A): %.1e\n", cond_number_exact_2norm);
                    fprintf("Exact computation time: %.1e\n", t_cond_exact);
                else
                    fprintf("Skip explicit matrix cond for p=%d (rule: only condest for 8<=p<=9).\n", p);
                end
                clear A_exact ATA_exact;
            end

            save("./data/typeII_2d_results_" + string(p) + "_" + string(beta) + "_tol_" + string(tol_hss) + ".mat", ...
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
                "t_iter", "iter", "rel_res_iter", "rel_acc_iter", "flag", "resvec", "P_reconstruct_iter", ...
                ...
                "t_pre", "t_total", ...
                ...
                "cond_number_est", "cond_number_exact", "cond_number_exact_2norm", ...
                ...
                "t_cond_est", "t_cond_exact");
            clear A;
            clear P u_ex f_nufft;
            clear u_direct r_direct e_direct P_reconstruct_direct;
            clear u_iter r_iter e_iter P_reconstruct_iter;
        end
    end
end