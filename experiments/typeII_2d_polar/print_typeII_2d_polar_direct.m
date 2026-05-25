exp_typeII_2d_polar_settings;

originalPath = path;
addpath("../")

for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    beta_string = beta_string_list(it_beta);
    beta_display = beta_display_list(it_beta);
    beta_factor_display = beta_factor_display_list(it_beta);
    caption_name = "Construction time, " + ...
        "factorization time, solve time, " + ...
        "relative residual " + ...
        "and relative error " + ...
        "of the direct solver " + ...
        "on a polar grid. " + ...
        "\(M \approx " + beta_factor_display + "N \log_{4}(N)\).";
    label_name = "typeII_2d_polar_direct_" + beta_string;

    result_list = [];
    for it_n = 1 : num_n
        p = p_list(it_n);

        result_list_tmp = [];
        for it_tol = 1 : num_tol_hss
            tol_hss = tol_hss_list(it_tol);
            load("./data/typeII_2d_results" ...
                + "_" + string(p) ...
                + "_" + string(beta) ...
                + "_" + string(tol_hss) ...
                + ".mat");
            result_list_tmp = [result_list_tmp, result];
        end
        result_list = [result_list; result_list_tmp];
    end

    print_table_direct(result_list, ...
        caption_name, label_name, ...
        tol_hss_display_list);
end

path(originalPath);
