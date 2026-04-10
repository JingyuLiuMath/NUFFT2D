clear;
close all;

p_list = (5 : 9)';
% p_list = (5 : 8)';
num_n = length(p_list);

tol_hss_list = [1e-2, 1e-4];
% tol_hss_list = [1e-2];
num_tol = length(tol_hss_list);

beta_list = [0.6, 0.7, 0.8]';
num_beta = length(beta_list);

n_list = 2.^p_list;
N_list = n_list.^2;

M_list = zeros(num_n, num_beta);

hss_rank_list = zeros(num_n, num_beta);
mem_list = zeros(num_n, num_beta);

t_construct_list = zeros(num_n, num_beta);
t_factor_list = zeros(num_n, num_beta);
t_solve_list = zeros(num_n, num_beta);

rel_res_list = zeros(num_n, num_beta);
rel_acc_list = zeros(num_n, num_beta);

colors = ["#0072BD", "#D95319", "#77AC30", "#EDB120", "#7E2F8E", "#4DBEEE"];
markers = ["o", "+", "x", "s", "d", "^"];

tol_hss = 1e-4;  % Use 1e-4 for resvec plots

for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    
    for it_n = 1 : num_n
        p = p_list(it_n);

        load("./data/typeII_2d_results_" + string(p) + "_" + string(beta) + "_tol_" + string(tol_hss) + ".mat");

        M_list(it_n, it_beta) = M;

        hss_rank_list(it_n, it_beta) = hss_rank;
        mem_list(it_n, it_beta) = double_to_gb(mem);

        t_construct_list(it_n, it_beta) = t_construct;
        t_factor_list(it_n, it_beta) = t_factor;
        t_solve_list(it_n, it_beta) = t_direct;

        rel_res_list(it_n, it_beta) = rel_res_direct;
        rel_acc_list(it_n, it_beta) = rel_acc_direct;
    end
end

M_list

%% Plot: Residual History for each beta (different p)
for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    
    figure();
    hold on;
    
    for it_p = 1 : num_n
        p = p_list(it_p);
        n = n_list(it_p);
        N = N_list(it_p);
        
        load("./data/typeII_2d_results_" + string(p) + "_" + string(beta) + "_tol_" + string(tol_hss) + ".mat");
        
        iterations = 1:length(resvec);
        semilogy(iterations, resvec, ...
            "LineWidth", 1.5, ...
            "Marker", markers(it_p), ...
            "MarkerSize", 6, ...
            "Color", colors(it_p), ...
            "DisplayName", sprintf("$p = %d$ ($N = %d^2$)", p, n));
    end
    
    % Calculate average beta
    beta_av = mean(M_list(:, it_beta) ./ N_list ./ p_list);
    
    title(sprintf("Residual History (beta = %.2f, $M \\approx %.2f N \\log_4 N$)", beta, beta_av), ...
        "Interpreter", "latex");
    xlabel("Iteration", "Interpreter", "latex");
    ylabel("Residual", "Interpreter", "latex");
    lgd = legend("Location", "northeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 16);
    grid on;
    hold off;
    
    % Save figure
    saveas(gcf, "./figure/resvec_beta_" + string(beta) + ".png", "png");
    saveas(gcf, "./figure/resvec_beta_" + string(beta) + ".eps", "epsc");
end

%% Summary plot: all beta on one figure (for a fixed p)
for it_p = 1 : num_n
    p = p_list(it_p);
    n = n_list(it_p);
    N = N_list(it_p);
    
    figure();
    hold on;
    
    for it_beta = 1 : num_beta
        beta = beta_list(it_beta);
        
        load("./data/typeII_2d_results_" + string(p) + "_" + string(beta) + "_tol_" + string(tol_hss) + ".mat");
        
        iterations = 1:length(resvec);
        beta_av = M_list(it_p, it_beta) / N / p;
        
        semilogy(iterations, resvec, ...
            "LineWidth", 1.5, ...
            "Marker", markers(it_beta), ...
            "MarkerSize", 6, ...
            "Color", colors(it_beta), ...
            "DisplayName", sprintf("$\\beta = %.2f$, $M \\approx %.2f N \\log_4 N$", beta, beta_av));
    end
    
    title(sprintf("Residual History ($p = %d$, $N = %d^2$)", p, n), ...
        "Interpreter", "latex");
    xlabel("Iteration", "Interpreter", "latex");
    ylabel("Residual", "Interpreter", "latex");
    lgd = legend("Location", "northeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 16);
    grid on;
    hold off;
    
    % Save figure
    saveas(gcf, "./figure/resvec_p_" + string(p) + ".png", "png");
    saveas(gcf, "./figure/resvec_p_" + string(p) + ".eps", "epsc");
end

%% Original plots (only for first beta to maintain compatibility)
beta = beta_list(1);

for it_n = 1 : num_n
    p = p_list(it_n);

    load("./data/typeII_2d_results_" + string(p) + "_" + string(beta) + "_tol_1e-2.mat");

    M_list(it_n, 1) = M;

    hss_rank_list(it_n, 1) = hss_rank;
    mem_list(it_n, 1) = double_to_gb(mem);

    t_construct_list(it_n, 1) = t_construct;
    t_factor_list(it_n, 1) = t_factor;
    t_solve_list(it_n, 1) = t_direct;

    rel_res_list(it_n, 1) = rel_res_direct;
    rel_acc_list(it_n, 1) = rel_acc_direct;
end

for it_n = 1 : num_n
    p = p_list(it_n);

    load("./data/typeII_2d_results_" + string(p) + "_" + string(beta) + "_tol_1e-4.mat");

    M_list(it_n, 2) = M;

    hss_rank_list(it_n, 2) = hss_rank;
    mem_list(it_n, 2) = double_to_gb(mem);

    t_construct_list(it_n, 2) = t_construct;
    t_factor_list(it_n, 2) = t_factor;
    t_solve_list(it_n, 2) = t_direct;

    rel_res_list(it_n, 2) = rel_res_direct;
    rel_acc_list(it_n, 2) = rel_acc_direct;
end

figure();
for it_tol = 1 : 2
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
saveas(gcf, "./figure/polar_bad_hss_rank.png", "png");
saveas(gcf, "./figure/polar_bad_hss_rank.eps", "epsc");

figure();
for it_tol = 1 : 2
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
saveas(gcf, "./figure/polar_bad_memory.png", "png");
saveas(gcf, "./figure/polar_bad_memory.eps", "epsc");

figure();
for it_tol = 1 : 2
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
saveas(gcf, "./figure/polar_bad_t_construct.png", "png");
saveas(gcf, "./figure/polar_bad_t_construct.eps", "epsc");

figure();
for it_tol = 1 : 2
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
saveas(gcf, "./figure/polar_bad_t_factor.png", "png");
saveas(gcf, "./figure/polar_bad_t_factor.eps", "epsc");

figure();
for it_tol = 1 : 2
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
saveas(gcf, "./figure/polar_bad_t_solve.png", "png");
saveas(gcf, "./figure/polar_bad_t_solve.eps", "epsc");

figure();
for it_tol = 1 : 2
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
saveas(gcf, "./figure/polar_bad_rel_res.png", "png");
saveas(gcf, "./figure/polar_bad_rel_res.eps", "epsc");

figure();
for it_tol = 1 : 2
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
saveas(gcf, "./figure/polar_bad_rel_acc.png", "png");
saveas(gcf, "./figure/polar_bad_rel_acc.eps", "epsc");

fprintf("Figures saved to ./figure/\n");
