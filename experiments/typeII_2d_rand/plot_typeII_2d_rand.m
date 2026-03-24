clear;
close all;

p_list = 3 : 5;
num_n = length(p_list);

n_list = 2.^p_list;
N_list = 2.^n_list;

hss_rank_list = zeros(num_n, 1);
mem_list = zeros(num_n, 1);
mem_ratio_list = zeros(num_n, 1);

t_construct_list = zeros(num_n, 1);
t_apply_list = zeros(num_n, 1);
t_nufft_list = zeros(num_n, 1);
t_factor_list = zeros(num_n, 1);
t_solve_list = zeros(num_n, 1);
t_total_list = zeros(num_n, 1);

rel_err_list = zeros(num_n, 1);
rel_res_list = zeros(num_n, 1);

for it = 1 : num_n
    p = p_list(it);

    load("./data/typeII_2d_results_" + string(p) + ".mat");

    hss_rank_list(it) = hss_rank;
    mem_list(it) = mem;
    mem_ratio_list(it) = mem_ratio;

    t_construct_list(it) = t_construct;
    t_apply_list(it) = t_apply;
    t_nufft_list(it) = t_nufft;
    t_factor_list(it) = t_factor;
    t_solve_list(it) = t_solve;
    t_total_list(it) = t_total;

    rel_err_list(it) = rel_err;
    rel_res_list(it) = rel_res;
end

figure();
loglog(N_list, hss_rank_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "HSS Rank");
hold on;
title("HSS rank");
xlabel("$N$", "Interpreter", "latex");
ylabel("$k$", "Interpreter", "latex");
legend("Location", "southeast");
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/hss_rank.png", "png");
saveas(gcf, "./figure/hss_rank.eps", "epsc");


figure();
loglog(N_list, mem_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "Memory");
hold on;
title("Memory");
xlabel("$N$", "Interpreter", "latex");
ylabel("Memory");
legend("Location", "southeast");
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/memory.png", "png");
saveas(gcf, "./figure/memory.eps", "epsc");

figure();
loglog(N_list, mem_ratio_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "Memory Saving");
hold on;
title("Memory Saving");
xlabel("$N$", "Interpreter", "latex");
ylabel("Memory Ratio");
legend("Location", "southeast");
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/memory_saving.png", "png");
saveas(gcf, "./figure/memory_saving.eps", "epsc");

figure();
loglog(N_list, t_construct_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "t_{c}");
hold on;
title("Construct Time");
xlabel("$N$", "Interpreter", "latex");
ylabel("time (s)");
legend("Location", "southeast");
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/t_construct.png", "png");
saveas(gcf, "./figure/t_construct.eps", "epsc");

figure();
loglog(N_list, t_apply_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "t_{a}");
hold on;
loglog(N_list, t_nufft_list, ...
    "LineWidth", 2, ...
    "Marker", "+", ...
    "Color", "b", ...
    "DisplayName", "t_{nufft}");
title("Apply Time");
xlabel("$N$", "Interpreter", "latex");
ylabel("time (s)");
legend("Location", "southeast");
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/t_apply.png", "png");
saveas(gcf, "./figure/t_apply.eps", "epsc");

figure();
loglog(N_list, t_factor_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "t_{f}");
hold on;
title("Factor Time");
xlabel("$N$", "Interpreter", "latex");
ylabel("time (s)");
legend("Location", "southeast");
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/t_factor.png", "png");
saveas(gcf, "./figure/t_factor.eps", "epsc");

figure();
loglog(N_list, t_solve_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "t_{s}");
hold on;
title("Solve Time");
xlabel("$N$", "Interpreter", "latex");
ylabel("time (s)");
legend("Location", "southeast");
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/t_solve.png", "png");
saveas(gcf, "./figure/t_solve.eps", "epsc");

figure();
loglog(N_list, rel_err_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "e_{a}");
hold on;
title("Relatve Error");
xlabel("$N$", "Interpreter", "latex");
ylabel("Error");
legend("Location", "southeast");
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/rel_err.png", "png");
saveas(gcf, "./figure/rel_err.eps", "epsc");

figure();
loglog(N_list, rel_res_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "r_{s}");
hold on;
title("Relatve Residual");
xlabel("$N$", "Interpreter", "latex");
ylabel("Residual");
legend("Location", "southeast");
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/rel_res.png", "png");
saveas(gcf, "./figure/rel_res.eps", "epsc");
