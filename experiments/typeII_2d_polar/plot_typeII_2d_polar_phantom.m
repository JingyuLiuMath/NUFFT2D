exp_typeII_2d_polar_settings;

originalPath = path;
addpath("../");

p = p_list(num_n);
n = 2^p;
fprintf("n: %d\n", n);

selected_row = round(n * 832 / 1024);
fprintf("selected_row: %d\n", selected_row);

for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    beta_string = beta_string_list(it_beta);

    load("./data/typeII_2d_points" ...
        + "_" + string(p) ...
        + "_" + string(beta) ...
        + ".mat");
    load("./data/typeII_2d_results_lsqr" ...
        + "_" + string(p) ...
        + "_" + string(beta) ...
        + ".mat");
    result_lsqr = result;

    figure_name = figure_prefix ...
        + "_" + beta_string ...
        + "_phantom";
    plot_phantom(n, c_ex, selected_row, figure_name + "_exact");
    plot_phantom(n, result_lsqr.c_cg, selected_row, figure_name + "_cg");

    for it_hss_tol = 1 : num_tol_hss
        tol_hss = tol_hss_list(it_hss_tol);
        load("./data/typeII_2d_results_face_splitting" ...
            + "_" + string(p) ...
            + "_" + string(beta) ...
            + "_" + string(tol_hss) ...
            + ".mat");

        figure_name_tol = figure_name + "_" + string(tol_hss);
        plot_phantom(n, result.c_direct, selected_row, ...
            figure_name_tol + "_direct");
        plot_phantom(n, result.c_pcg, selected_row, ...
            figure_name_tol + "_pcg");
    end
end

path(originalPath);
