function plot_figure_condest(result_list, ...
    figure_prefix)

num_n = size(result_list, 1);

N_list = zeros(num_n, 1);
N_ex_list = [];
kappa1_est_list = zeros(num_n, 1);
kappa1_ex_list = [];
kappa2_ex_list = [];

for it_n = 1 : num_n
    curr_result = result_list(it_n);
    kappa1_est_list(it_n) = curr_result.kappa1_est;
    if curr_result.n <= 2^7
        kappa1_ex_list = [kappa1_ex_list; curr_result.kappa1_ex];
        kappa2_ex_list = [kappa2_ex_list; curr_result.kappa2_ex];
        N_ex_list = [N_ex_list; curr_result.N];
    end
    N_list(it_n) = curr_result.N;
end

figure()
loglog(N_list, kappa1_est_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Markersize", 12, ...
    "DisplayName", "Estimated $\kappa_{1}$");
hold on;
loglog(N_ex_list, kappa1_ex_list, ...
    "LineWidth", 2, ...
    "Marker", "+", ...
    "Markersize", 12, ...
    "DisplayName", "$\kappa_{1}$");
loglog(N_ex_list, kappa2_ex_list, ...
    "LineWidth", 2, ...
    "Marker", "square", ...
    "Markersize", 12, ...
    "DisplayName", "$\kappa_{2}$");

figure_name = figure_prefix + "_cond";
xlabel("$N$", "Interpreter", "latex");
ylabel("$\kappa$", "Interpreter", "latex");
title("Condition number", "Interpreter", "latex");
legend("Location", "southeast", "Interpreter", "latex");
set(gca, 'FontSize', 24);
saveas(gcf, figure_name + ".png", "png");
saveas(gcf, figure_name + ".eps", "epsc");


end
