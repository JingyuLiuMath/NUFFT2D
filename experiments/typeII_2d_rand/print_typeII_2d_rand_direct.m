exp_typeII_2d_rand_settings;

originalPath = path;
addpath("../");

for it_alpha = 1 : num_alpha
    alpha = alpha_list(it_alpha);
    alpha_string = alpha_string_list(it_alpha);
    alpha_display = alpha_display_list(it_alpha);

    caption_name = "Construction time, " + ...
        "factorization time, solve time, " + ...
        "relative residual " + ...
        "and relative error " + ...
        "of the direct solver " + ...
        "on a random grid. " + ...
        "\(M = " + alpha_display + "N\).";
    label_name = "typeII_2d_rand_direct_" + alpha_string;

    result_list = [];
    for it_n = 1 : num_n
        p = p_list(it_n);

        result_list_tmp = [];
        for it_hss_tol = 1 : num_tol_hss
            tol_hss = tol_hss_list(it_hss_tol);
            load("./data/typeII_2d_results" ...
                + "_" + string(p) ...
                + "_" + string(alpha) ...
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
