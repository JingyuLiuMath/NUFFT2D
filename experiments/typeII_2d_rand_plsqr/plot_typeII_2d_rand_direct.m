clear;
close all;

p_list = (5 : 9)';
% p_list = (5 : 8)';
num_n = length(p_list);

tol_hss_list = [1e-2, 1e-4];
% tol_hss_list = [1e-2];
num_tol = length(tol_hss_list);

n_list = 2.^p_list;
N_list = n_list.^2;

M_list = zeros(num_n, num_tol);

hss_rank_list = zeros(num_n, num_tol);
mem_list = zeros(num_n, num_tol);

t_construct_list = zeros(num_n, num_tol);
t_factor_list = zeros(num_n, num_tol);
t_solve_list = zeros(num_n, num_tol);

rel_res_list = zeros(num_n, num_tol);
rel_acc_list = zeros(num_n, num_tol);

for it_tol = 1 : num_tol
    tol_hss = tol_hss_list(it_tol);
    fprintf("tol_hss: %.1e\n", tol_hss);
    for it_n = 1 : num_n
        p = p_list(it_n);

        load("./data/typeII_2d_results_" + string(p) + "_tol_" + string(tol_hss) + ".mat");

        M_list(it_n, it_tol) = M;

        hss_rank_list(it_n, it_tol) = hss_rank;
        mem_list(it_n, it_tol) = double_to_gb(mem);

        t_construct_list(it_n, it_tol) = t_construct;
        t_factor_list(it_n, it_tol) = t_factor;
        t_solve_list(it_n, it_tol) = t_direct;

        rel_res_list(it_n, it_tol) = rel_res_direct;
        rel_acc_list(it_n, it_tol) = rel_acc_direct;
    end
end

M_list

figure();
for it_tol = 1 : num_tol
    tol_hss = tol_hss_list(it_tol);
    if tol_hss == 1e-2
        marker = "x";
        color = "#0072BD";
        display_name = "tol = 1e-2";
    elseif tol_hss == 1e-4
        marker = "+";
        color = "#D95319";
        display_name = "tol = 1e-4";
    end
    loglog(N_list, hss_rank_list(:, it_tol), ...
        "LineWidth", 2, ...
        "Marker", marker, ...
        "Color", color, ...
        "DisplayName", display_name);
    hold on;
end
ref_line = sqrt(N_list) .* log2(N_list);
ref_line = ref_line / ref_line(1) * hss_rank_list(1, 2) * 1.2;
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "#77AC30", ...
    "DisplayName", "$O(\sqrt{N} \log N)$");
title("HSS Rank");
xlabel("$N$", "Interpreter", "latex");
ylabel("$k$", "Interpreter", "latex");
yticks([10^2, 10^3, 10^4]);
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/rand_hss_rank.png", "png");
saveas(gcf, "./figure/rand_hss_rank.eps", "epsc");

figure();
for it_tol = 1 : num_tol
    tol_hss = tol_hss_list(it_tol);
    if tol_hss == 1e-2
        marker = "x";
        color = "#0072BD";
        display_name = "tol = 1e-2";
    elseif tol_hss == 1e-4
        marker = "+";
        color = "#D95319";
        display_name = "tol = 1e-4";
    end
    loglog(N_list, mem_list(:, it_tol), ...
        "LineWidth", 2, ...
        "Marker", marker, ...
        "Color", color, ...
        "DisplayName", display_name);
    hold on;
end
ref_line = M_list(:, 1) .* sqrt(N_list) .* log2(N_list);
ref_line = ref_line / ref_line(1) * mem_list(1, 2) * 1.2;
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "#77AC30", ...
    "DisplayName", "$O(M \sqrt{N} \log N)$");
title("Memory");
xlabel("$N$", "Interpreter", "latex");
ylabel("Memory (GB)");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/rand_memory.png", "png");
saveas(gcf, "./figure/rand_memory.eps", "epsc");

figure();
for it_tol = 1 : num_tol
    tol_hss = tol_hss_list(it_tol);
    if tol_hss == 1e-2
        marker = "x";
        color = "#0072BD";
        display_name = "tol = 1e-2";
    elseif tol_hss == 1e-4
        marker = "+";
        color = "#D95319";
        display_name = "tol = 1e-4";
    end
    loglog(N_list, t_construct_list(:, it_tol), ...
        "LineWidth", 2, ...
        "Marker", marker, ...
        "Color", color, ...
        "DisplayName", display_name);
    hold on;
end
ref_line = M_list(:, 1) .* N_list .* log2(N_list).^2;
ref_line = ref_line / ref_line(1) * t_construct_list(1, 2) * 1.2;
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "#77AC30", ...
    "DisplayName", "$O(M N \log^{2} N)$");
title("Construct Time");
xlabel("$N$", "Interpreter", "latex");
ylabel("time (s)");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/rand_t_construct.png", "png");
saveas(gcf, "./figure/rand_t_construct.eps", "epsc");

figure();
for it_tol = 1 : num_tol
    tol_hss = tol_hss_list(it_tol);
    if tol_hss == 1e-2
        marker = "x";
        color = "#0072BD";
        display_name = "tol = 1e-2";
    elseif tol_hss == 1e-4
        marker = "+";
        color = "#D95319";
        display_name = "tol = 1e-4";
    end
    loglog(N_list, t_factor_list(:, it_tol), ...
        "LineWidth", 2, ...
        "Marker", marker, ...
        "Color", color, ...
        "DisplayName", display_name);
    hold on;
end
ref_line = M_list(:, 1) .* N_list .* log2(N_list).^2;
ref_line = ref_line / ref_line(1) * t_factor_list(1, 2) * 1.2;
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "#77AC30", ...
    "DisplayName", "$O(M N \log^{2} N)$");
title("Factor Time");
xlabel("$N$", "Interpreter", "latex");
ylabel("time (s)");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/rand_t_factor.png", "png");
saveas(gcf, "./figure/rand_t_factor.eps", "epsc");

figure();
for it_tol = 1 : num_tol
    tol_hss = tol_hss_list(it_tol);
    if tol_hss == 1e-2
        marker = "x";
        color = "#0072BD";
        display_name = "tol = 1e-2";
    elseif tol_hss == 1e-4
        marker = "+";
        color = "#D95319";
        display_name = "tol = 1e-4";
    end
    loglog(N_list, t_solve_list(:, it_tol), ...
        "LineWidth", 2, ...
        "Marker", marker, ...
        "Color", color, ...
        "DisplayName", display_name);
    hold on;
end
ref_line = M_list(:, 1) .* sqrt(N_list) .* log2(N_list);
ref_line = ref_line / ref_line(1) * t_solve_list(1, 2) * 1.2;
loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "#77AC30", ...
    "DisplayName", "$O(M \sqrt{N} \log N)$");
title("Solve Time");
xlabel("$N$", "Interpreter", "latex");
ylabel("time (s)");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/rand_t_solve.png", "png");
saveas(gcf, "./figure/rand_t_solve.eps", "epsc");

figure();
for it_tol = 1 : num_tol
    tol_hss = tol_hss_list(it_tol);
    if tol_hss == 1e-2
        marker = "x";
        color = "#0072BD";
        display_name = "tol = 1e-2";
    elseif tol_hss == 1e-4
        marker = "+";
        color = "#D95319";
        display_name = "tol = 1e-4";
    end
    semilogx(N_list, rel_res_list(:, it_tol), ...
        "LineWidth", 2, ...
        "Marker", marker, ...
        "Color", color, ...
        "DisplayName", display_name);
    hold on;
end
title("Relatve Residual");
xlabel("$N$", "Interpreter", "latex");
ylabel("Residual");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/rand_rel_res.png", "png");
saveas(gcf, "./figure/rand_rel_res.eps", "epsc");

figure();
for it_tol = 1 : num_tol
    tol_hss = tol_hss_list(it_tol);
    if tol_hss == 1e-2
        marker = "x";
        color = "#0072BD";
        display_name = "tol = 1e-2";
    elseif tol_hss == 1e-4
        marker = "+";
        color = "#D95319";
        display_name = "tol = 1e-4";
    end
    semilogx(N_list, rel_acc_list(:, it_tol), ...
        "LineWidth", 2, ...
        "Marker", marker, ...
        "Color", color, ...
        "DisplayName", display_name);
    hold on;
end
title("Relatve Accuracy");
xlabel("$N$", "Interpreter", "latex");
ylabel("Accuracy");
lgd = legend("Location", "southeast");
lgd.Interpreter = 'latex';
set(gca, 'FontSize', 18);
saveas(gcf, "./figure/rand_rel_acc.png", "png");
saveas(gcf, "./figure/rand_rel_acc.eps", "epsc");
