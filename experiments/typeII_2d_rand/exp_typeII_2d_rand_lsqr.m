exp_typeII_2d_rand_settings;

for it_alpha = 1 : num_alpha
    for it_p = 1 : num_n
        alpha = alpha_list(it_alpha);
        p = p_list(it_p);
        n = 2^p;

        load("./data/typeII_2d_points" ...
            + "_" + string(p) ...
            + "_" + string(alpha) ...
            + ".mat");

        fprintf("\n\n\n\n");
        fprintf("p: %d\n", p);

        result = run_INUDFT2_2D_LSQR(x, n, tol_cg, maxit_cg);

        save("./data/typeII_2d_results_lsqr" ...
            + "_" + string(p) ...
            + "_" + string(alpha) ...
            + ".mat", "result");
    end
end