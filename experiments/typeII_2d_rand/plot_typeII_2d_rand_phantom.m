exp_typeII_2d_rand_settings;

originalPath = path;
addpath("../");

p = p_list(num_n);
n = 2^p;
fprintf("n: %d\n", n);

selected_row = ceil(n * 0.7);
fprintf("selected_row: %d\n", selected_row);


for it_alpha = 1 : num_alpha
    alpha = alpha_list(it_alpha);
    alpha_string = alpha_string_list(it_alpha);

    load("./data/typeII_2d_points" ...
        + "_" + string(p) ...
        + "_" + string(alpha) ...
        + ".mat");
    load("./data/typeII_2d_results_lsqr" ...
        + "_" + string(p) ...
        + "_" + string(alpha) ...
        + ".mat");
    result_lsqr = result;

    figure_name = figure_prefix ...
        + "_" + alpha_string ...
        + "_phantom";
    plot_phantom(n, c_ex, selected_row, figure_name + "_exact");
    plot_phantom(n, result_lsqr.c_cg, selected_row, figure_name + "_cg");

    for it_hss_tol = 1 : num_tol_hss
        tol_hss = tol_hss_list(it_hss_tol);
        load("./data/typeII_2d_results_face_splitting" ...
            + "_" + string(p) ...
            + "_" + string(alpha) ...
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
