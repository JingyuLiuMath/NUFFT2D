clear;
close all;

p_list = (4 : 9)';
num_n = length(p_list);

n_list = 2.^p_list;
N_list = n_list.^2;

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
ref_line = sqrt(N_list) .* log2(N_list);
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "b", ...
    "DisplayName", "$O(\sqrt{N} \log N)$");
title("HSS Rank");
xlabel("$N$", "Interpreter", "latex");
ylabel("$k$", "Interpreter", "latex");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/hss_rank.png", "png");
saveas(gcf, "./figure/hss_rank.eps", "epsc");

figure();
loglog(N_list, hss_rank_list ./ N_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "HSS Rank / N");
hold on;
title("HSS Rank / N");
xlabel("$N$", "Interpreter", "latex");
ylabel("$k$", "Interpreter", "latex");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/hss_rank_ratio.png", "png");
saveas(gcf, "./figure/hss_rank_ratio.eps", "epsc");

figure();
loglog(N_list, mem_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "Memory");
hold on;
ref_line = N_list .* sqrt(N_list) .* log2(N_list);
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "b", ...
    "DisplayName", "$O(N^{3/ 2} \log N)$");
title("Memory");
xlabel("$N$", "Interpreter", "latex");
ylabel("Memory");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
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
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/memory_saving.png", "png");
saveas(gcf, "./figure/memory_saving.eps", "epsc");

figure();
loglog(N_list, t_construct_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "$t_{c}$");
hold on;
ref_line = N_list.^(5 / 2) .* log2(N_list);
ref_line = ref_line / ref_line(1) * t_construct_list(1) / 2;
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "b", ...
    "DisplayName", "$O(N^{5 / 2} \log N)$");
title("Construct Time");
xlabel("$N$", "Interpreter", "latex");
ylabel("time (s)");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/t_construct.png", "png");
saveas(gcf, "./figure/t_construct.eps", "epsc");

figure();
loglog(N_list, t_apply_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "$t_{a}$");
hold on;
loglog(N_list, t_nufft_list, ...
    "LineWidth", 2, ...
    "Marker", "+", ...
    "Color", "m", ...
    "DisplayName", "$t_{nufft}$");
ref_line = N_list.^(3 / 2) .* log2(N_list);
ref_line = ref_line / ref_line(1) * t_apply_list(1) / 2;
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "b", ...
    "DisplayName", "$O(N^{3/2} \log N)$");
title("Apply Time");
xlabel("$N$", "Interpreter", "latex");
ylabel("time (s)");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/t_apply.png", "png");
saveas(gcf, "./figure/t_apply.eps", "epsc");

figure();
loglog(N_list, t_factor_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "$t_{f}$");
hold on;
ref_line = (N_list.^2) .* (log2(N_list).^2);
ref_line = ref_line / ref_line(1) * t_factor_list(1) / 2;
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "b", ...
    "DisplayName", "$O(N^{2} \log^2 N)$");
title("Factor Time");
xlabel("$N$", "Interpreter", "latex");
ylabel("time (s)");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/t_factor.png", "png");
saveas(gcf, "./figure/t_factor.eps", "epsc");

figure();
loglog(N_list, t_solve_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "$t_{s}$");
hold on;
ref_line = N_list.^(3 / 2) .* log2(N_list);
ref_line = ref_line / ref_line(1) * t_solve_list(1) / 2;
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "b", ...
    "DisplayName", "$O(N^{3/2} \log N)$");
title("Solve Time");
xlabel("$N$", "Interpreter", "latex");
ylabel("time (s)");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/t_solve.png", "png");
saveas(gcf, "./figure/t_solve.eps", "epsc");

figure();
loglog(N_list, rel_err_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "$e_{a}$");
hold on;
title("Relatve Error");
xlabel("$N$", "Interpreter", "latex");
ylabel("Error");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/rel_err.png", "png");
saveas(gcf, "./figure/rel_err.eps", "epsc");

figure();
loglog(N_list, rel_res_list, ...
    "LineWidth", 2, ...
    "Marker", "x", ...
    "Color", "r", ...
    "DisplayName", "$r_{s}$");
hold on;
title("Relatve Residual");
xlabel("$N$", "Interpreter", "latex");
ylabel("Residual");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/rel_res.png", "png");
saveas(gcf, "./figure/rel_res.eps", "epsc");
