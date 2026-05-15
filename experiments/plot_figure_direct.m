function plot_figure_direct(result_list, ...
    figure_prefix, ...
    tol_hss_display_list)

num_n = size(result_list, 1);
num_hss_tol = size(result_list, 2);

N_list = zeros(num_n, 1);
hss_rank_list = zeros(num_n, num_hss_tol);
t_construct_list = zeros(num_n, num_hss_tol);
t_factor_list = zeros(num_n, num_hss_tol);
t_solve_list = zeros(num_n, num_hss_tol);
rel_res_list = zeros(num_n, num_hss_tol);
rel_err_list = zeros(num_n, num_hss_tol);

for it_n = 1 : num_n
    for it_hss_tol = 1 : num_hss_tol
        curr_result = result_list(it_n, it_hss_tol);
        hss_rank_list(it_n, it_hss_tol) = curr_result.hss_rank;
        t_construct_list(it_n, it_hss_tol) = curr_result.t_construct;
        t_factor_list(it_n, it_hss_tol) = curr_result.t_factor;
        t_solve_list(it_n, it_hss_tol) = curr_result.t_direct;

        rel_res_list(it_n, it_hss_tol) = curr_result.rel_res_direct;
        rel_err_list(it_n, it_hss_tol) = curr_result.rel_err_direct;
    end
    N_list(it_n) = curr_result.N;
end

display_name_list = [];
for it_hss_tol = 1 : num_hss_tol
    tol_hss_display_string = tol_hss_display_list(it_hss_tol);
    display_name = "\varepsilon = " + tol_hss_display_string;
    display_name_list = [display_name_list; display_name];
end

scaling_type = "O(N^{0.5} \log (N))";
xlabel_name = "$N$";
ylabel_name = "$r_{\mathrm{h}}$";
title_name = "HSS rank";
my_name = "_hss_rank";
figure_name = figure_prefix + my_name;
plot_scaling(N_list, hss_rank_list, ...
    scaling_type, ...
    display_name_list, ...
    xlabel_name, ylabel_name, title_name, figure_name);

scaling_type = "O(N^{1.5} \log^{3} (N))";
xlabel_name = "$N$";
ylabel_name = "$t_{\mathrm{c}}$ (s)";
title_name = "Construction time";
my_name = "_t_construct";
figure_name = figure_prefix + my_name;
plot_scaling(N_list, t_construct_list, ...
    scaling_type, ...
    display_name_list, ...
    xlabel_name, ylabel_name, title_name, figure_name);

scaling_type = "O(N^{1.5} \log^{3} (N))";
xlabel_name = "$N$";
ylabel_name = "$t_{\mathrm{f}}$ (s)";
title_name = "Factorization time";
my_name = "_t_factor";
figure_name = figure_prefix + my_name;
plot_scaling(N_list, t_factor_list, ...
    scaling_type, ...
    display_name_list, ...
    xlabel_name, ylabel_name, title_name, figure_name);

scaling_type = "O(N \log^{3} (N))";
xlabel_name = "$N$";
ylabel_name = "$t_{\mathrm{s}}$ (s)";
title_name = "Solution time";
my_name = "_t_solve";
figure_name = figure_prefix + my_name;
plot_scaling(N_list, t_solve_list, ...
    scaling_type, ...
    display_name_list, ...
    xlabel_name, ylabel_name, title_name, figure_name);

xlabel_name = "$N$";
ylabel_name = "$r_{\mathrm{s}}$ (s)";
title_name = "Relative residual";
my_name = "r_s";
figure_name = figure_prefix + my_name;
plot_err(N_list, rel_res_list, ...
    display_name_list, ...
    xlabel_name, ylabel_name, title_name, figure_name);

xlabel_name = "$N$";
ylabel_name = "$e_{\mathrm{s}}$ (s)";
title_name = "Relative error";
my_name = "e_s";
figure_name = figure_prefix + my_name;
plot_err(N_list, rel_err_list, ...
    display_name_list, ...
    xlabel_name, ylabel_name, title_name, figure_name);

end

function plot_scaling(N_list, t_list, ...
    scaling_type, ...
    display_name_list, xlabel_name, ylabel_name, title_name, figure_name)

num_n = size(t_list, 1);
num_param = size(t_list, 2);

figure();
for it_param = 1 : num_param
    loglog(N_list, t_list(:, it_param), ...
        "LineWidth", 2, ...
        "Marker", "x", ...
        "Markersize", 12, ...
        "DisplayName", "$" + display_name_list(it_param) + "$");
    hold on;
end

switch title_name
    case "HSS rank"
        scaling_factor = 1.5;
    case "Construction time"
        scaling_factor = 2;
    case "Factorization time"
        scaling_factor = 4;
    case "Solution time"
        scaling_factor = 10;
end

switch scaling_type
    case "O(N^{0.5} \log (N))"
        ref_line = sqrt(N_list) .* log2(N_list);
        ref_line = ref_line / ref_line(1) * max(t_list(1, :)) * scaling_factor;
    case "O(N \log (N))"
        ref_line = N_list .* log2(N_list);
        ref_line = ref_line / ref_line(1) * max(t_list(1, :)) * scaling_factor;
    case "O(N \log^2 (N))"
        ref_line = N_list .* (log2(N_list).^2);
        ref_line = ref_line / ref_line(1) * max(t_list(1, :)) * scaling_factor;
    case "O(N \log^{3} (N))"
        ref_line = N_list .* log2(N_list).^3;
        ref_line = ref_line / ref_line(1) * max(t_list(1, :)) * scaling_factor;
    case "O(N^{1.5} \log^{3} (N))"
        ref_line = N_list.^(1.5) .* log2(N_list).^3;
        ref_line = ref_line / ref_line(1) * max(t_list(1, :)) * scaling_factor;
end
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "DisplayName", "$" + scaling_type + "$");
xlabel(xlabel_name, "Interpreter", "latex");
ylabel(ylabel_name, "Interpreter", "latex");
title(title_name, "Interpreter", "latex");
switch title_name
    case "HSS rank"
        yticks([10^2 10^3 10^4]);
    case "Construction time"
        yticks([10^0 10^2 10^4]);
    case "Factorization time"
        yticks([10^(-2) 10^0 10^2 10^4]);
    case "Solution time"
        yticks([10^(-2) 10^0 10^2]);
end

legend("Location", "southeast", "Interpreter", "latex");
set(gca, 'FontSize', 24);
saveas(gcf, figure_name + ".png", "png");
saveas(gcf, figure_name + ".eps", "epsc");

end

function plot_err(N_list, err_list, ...
    display_name_list, xlabel_name, ylabel_name, title_name, figure_name)

num_n = size(err_list, 1);
num_param = size(err_list, 2);

figure();
for it_param = 1 : num_param
    loglog(N_list, err_list(:, it_param), ...
        "LineWidth", 2, ...
        "Marker", "x", ...
        "Markersize", 12, ...
        "DisplayName", "$" + display_name_list(it_param) + "$");
    hold on;
end
xlabel(xlabel_name, "Interpreter", "latex");
ylabel(ylabel_name, "Interpreter", "latex");
title(title_name, "Interpreter", "latex");
legend("Location", "southeast", "Interpreter", "latex");
set(gca, 'FontSize', 24);
saveas(gcf, figure_name + ".png", "png");
saveas(gcf, figure_name + ".eps", "epsc");

end
