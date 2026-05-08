exp_typeII_2d_polar_settings;

originalPath = path;
addpath("../");

for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    beta_string = beta_string_list(it_beta);
    beta_display = beta_display_list(it_beta);

    result_list = [];
    for it_n = 1 : num_n
        p = p_list(it_n);

        tol_hss = tol_hss_list(1);
        load("./data/typeII_2d_results_condest" ...
            + "_" + string(p) ...
            + "_" + string(beta) ...
            + "_" + string(tol_hss) ...
            + ".mat");
        result_list = [result_list; result];
    end

    my_figure_prefix = figure_prefix + "_" + beta_string;
    plot_figure_condest(result_list, ...
        my_figure_prefix);
end

path(originalPath);
