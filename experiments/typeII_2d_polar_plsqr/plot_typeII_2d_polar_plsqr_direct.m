clear;
close all;

p_list = (5 : 9)';
num_n = length(p_list);

beta_list = [0.6, 0.7]';
num_beta = length(beta_list);

tol_hss_list = [1e-2, 1e-4];
num_tol = length(tol_hss_list);

n_list = 2.^p_list;
N_list = n_list.^2;

M_list = zeros(num_n, num_tol, num_beta);

hss_rank_list = zeros(num_n, num_tol, num_beta);
mem_list = zeros(num_n, num_tol, num_beta);

t_construct_list = zeros(num_n, num_tol, num_beta);
t_factor_list = zeros(num_n, num_tol, num_beta);
t_solve_list = zeros(num_n, num_tol, num_beta);

rel_res_list = zeros(num_n, num_tol, num_beta);
rel_acc_list = zeros(num_n, num_tol, num_beta);

cond_est_list = zeros(num_n, num_tol, num_beta);
cond_exact_list = zeros(num_n, num_tol, num_beta);
cond_exact_2norm_list = zeros(num_n, num_tol, num_beta);

%% Load all data
for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        fprintf("beta: %.1f, tol_hss: %.1e\n", beta, tol_hss);
        for it_n = 1 : num_n
            p = p_list(it_n);

            load("./data/typeII_2d_results_" + string(p) + "_" + string(beta) + "_" + string(tol_hss) + ".mat");

            M_list(it_n, it_tol, it_beta) = M;

            hss_rank_list(it_n, it_tol, it_beta) = hss_rank;
            mem_list(it_n, it_tol, it_beta) = double_to_gb(mem);

            t_construct_list(it_n, it_tol, it_beta) = t_construct;
            t_factor_list(it_n, it_tol, it_beta) = t_factor;
            t_solve_list(it_n, it_tol, it_beta) = t_direct;

            rel_res_list(it_n, it_tol, it_beta) = rel_res_direct;
            rel_acc_list(it_n, it_tol, it_beta) = rel_acc_direct;

            cond_est_list(it_n, it_tol, it_beta) = cond_number_est;
            cond_exact_list(it_n, it_tol, it_beta) = cond_number_exact;
            cond_exact_2norm_list(it_n, it_tol, it_beta) = cond_number_exact_2norm;
        end
    end
end


coeff = mean(M_list(:, 1, :) ./ N_list ./ p_list);


%% Plot for each beta
markers = ["x", "+"];
colors = ["#0072BD", "#D95319"];
ref_color = "#77AC30";

for it_beta = 1 : num_beta
    beta = beta_list(it_beta);

    % HSS Rank
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "$\varepsilon$ = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "$\varepsilon$ = 1e-4";
        end
        loglog(N_list, hss_rank_list(:, it_tol, it_beta), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    % Reference line: O(sqrt(N) * log(N))
    ref_line = sqrt(N_list) .* log2(N_list);
    ref_line = ref_line / ref_line(1) * hss_rank_list(1, 2, it_beta) * 1.1;
    loglog(N_list, ref_line, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Color", ref_color, ...
        "DisplayName", "$O(\sqrt{N} \\log N)$");
    if beta == 0.7
        title(sprintf("HSS Rank ($M \\approx 0.55 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    elseif beta == 0.6
        title(sprintf("HSS Rank ($M \\approx 0.47 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$r_{\\mathrm{h}}$", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/polar_hss_rank_beta_" + string(beta) + ".png", "png");
    saveas(gcf, "./figure/polar_hss_rank_beta_" + string(beta) + ".eps", "epsc");

    % Memory
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "$\varepsilon$ = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "$\varepsilon$ = 1e-4";
        end
        loglog(N_list, mem_list(:, it_tol, it_beta), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    % Reference line: O(M * sqrt(N) * log(N))
    % ref_line = M_list(:, 1, it_beta) .* sqrt(N_list) .* log2(N_list);
    ref_line = N_list .* log2(N_list).^3;
    ref_line = ref_line / ref_line(1) * mem_list(1, 2, it_beta) * 1.1;
    loglog(N_list, ref_line, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Color", ref_color, ...
        "DisplayName", "$O(N \log^{3} N)$");
    if beta == 0.7
        title(sprintf("Memory ($M \\approx 0.55 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    elseif beta == 0.6
        title(sprintf("Memory ($M \\approx 0.47 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$m_{\\mathrm{h}}$ (GB)", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/polar_memory_beta_" + string(beta) + ".png", "png");
    saveas(gcf, "./figure/polar_memory_beta_" + string(beta) + ".eps", "epsc");

    % Construct Time
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "$\varepsilon$ = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "$\varepsilon$ = 1e-4";
        end
        loglog(N_list, t_construct_list(:, it_tol, it_beta), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    % Reference line: O(M * N * log^2(N))
    % ref_line = M_list(:, 1, it_beta) .* N_list .* log2(N_list).^2;
    ref_line = N_list.^(1.5) .* log2(N_list).^3;
    ref_line = ref_line / ref_line(1) * t_construct_list(1, 2, it_beta) * 1.1;
    loglog(N_list, ref_line, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Color", ref_color, ...
        "DisplayName", "$O(N^{3 / 2} \log^{3} N)$");
    if beta == 0.7
        title(sprintf("Construct Time ($M \\approx 0.55 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    elseif beta == 0.6
        title(sprintf("Construct Time ($M \\approx 0.47 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$t_{\\mathrm{c}}$ (s)", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/polar_t_construct_beta_" + string(beta) + ".png", "png");
    saveas(gcf, "./figure/polar_t_construct_beta_" + string(beta) + ".eps", "epsc");

    % Factor Time
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "$\varepsilon$ = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "$\varepsilon$ = 1e-4";
        end
        loglog(N_list, t_factor_list(:, it_tol, it_beta), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    % Reference line: O(M * N * log^2(N))
    % ref_line = M_list(:, 1, it_beta) .* N_list .* log2(N_list).^2;
    ref_line = N_list.^(1.5) .* log2(N_list).^3;
    ref_line = ref_line / ref_line(1) * t_factor_list(1, 2, it_beta) * 1.1;
    loglog(N_list, ref_line, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Color", ref_color, ...
        "DisplayName", "$O(N^{3 / 2} \log^{3} N)$");
    if beta == 0.7
        title(sprintf("Factor Time ($M \\approx 0.55 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    elseif beta == 0.6
        title(sprintf("Factor Time ($M \\approx 0.47 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$t_{\\mathrm{f}}$ (s)", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/polar_t_factor_beta_" + string(beta) + ".png", "png");
    saveas(gcf, "./figure/polar_t_factor_beta_" + string(beta) + ".eps", "epsc");

    % Solve Time
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "$\varepsilon$ = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "$\varepsilon$ = 1e-4";
        end
        loglog(N_list, t_solve_list(:, it_tol, it_beta), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    % Reference line: O(M * sqrt(N) * log(N))
    % ref_line = M_list(:, 1, it_beta) .* sqrt(N_list) .* log2(N_list);
    ref_line = N_list .* log2(N_list).^3;
    ref_line = ref_line / ref_line(1) * t_solve_list(1, 2, it_beta) * 1.1;
    loglog(N_list, ref_line, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Color", ref_color, ...
        "DisplayName", "$O(N \log^{3} N)$");
    if beta == 0.7
        title(sprintf("Solve Time ($M \\approx 0.55 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    elseif beta == 0.6
        title(sprintf("Solve Time ($M \\approx 0.47 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$t_{\\mathrm{s}}$ (s)", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/polar_t_solve_beta_" + string(beta) + ".png", "png");
    saveas(gcf, "./figure/polar_t_solve_beta_" + string(beta) + ".eps", "epsc");

    % Relative Residual
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "$\varepsilon$ = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "$\varepsilon$ = 1e-4";
        end
        loglog(N_list, rel_res_list(:, it_tol, it_beta), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    if beta == 0.7
        title(sprintf("Relative Residual ($M \\approx 0.55 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    elseif beta == 0.6
        title(sprintf("Relative Residual ($M \\approx 0.47 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$r_{\\mathrm{s}}$", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/polar_rel_res_beta_" + string(beta) + ".png", "png");
    saveas(gcf, "./figure/polar_rel_res_beta_" + string(beta) + ".eps", "epsc");

    % Relative Accuracy
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "$\varepsilon$ = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "$\varepsilon$ = 1e-4";
        end
        loglog(N_list, rel_acc_list(:, it_tol, it_beta), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    if beta == 0.7
        title(sprintf("Relative Accuracy ($M \\approx 0.55 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    elseif beta == 0.6
        title(sprintf("Relative Accuracy ($M \\approx 0.47 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$e_{\\mathrm{s}}$", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/polar_rel_acc_beta_" + string(beta) + ".png", "png");
    saveas(gcf, "./figure/polar_rel_acc_beta_" + string(beta) + ".eps", "epsc");

    % Condition Number (only tol_hss == 1e-4)
    figure();
    it_tol = find(tol_hss_list == 1e-4, 1);

    valid_ind = ~isnan(cond_est_list(:, it_tol, it_beta));
    if any(valid_ind)
        loglog(N_list(valid_ind), cond_est_list(valid_ind, it_tol, it_beta), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", "#0072BD", ...
            "DisplayName", "Estimated cond (1-norm)");
        hold on;
    end

    valid_ind = ~isnan(cond_exact_2norm_list(:, it_tol, it_beta));
    if any(valid_ind)
        loglog(N_list(valid_ind), cond_exact_2norm_list(valid_ind, it_tol, it_beta), ...
            "LineWidth", 2, ...
            "Marker", "+", ...
            "LineStyle", "--", ...
            "Color", "#D95319", ...
            "DisplayName", "Exact cond (2-norm)");
        hold on;
    end

    valid_ind = ~isnan(cond_exact_list(:, it_tol, it_beta));
    if any(valid_ind)
        loglog(N_list(valid_ind), cond_exact_list(valid_ind, it_tol, it_beta), ...
            "LineWidth", 2, ...
            "Marker", "o", ...
            "LineStyle", ":", ...
            "Color", "#77AC30", ...
            "DisplayName", "Exact cond (1-norm)");
        hold on;
    end

    if beta == 0.7
        title(sprintf("Condition Number ($M \\approx 0.55 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    elseif beta == 0.6
        title(sprintf("Condition Number ($M \\approx 0.47 N \\log_{4} (N) $)", beta), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$\kappa$", "Interpreter", "latex");
    lgd = legend("Location", "northwest");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    hold off;
    saveas(gcf, "./figure/polar_cond_beta_" + string(beta) + ".png", "png");
    saveas(gcf, "./figure/polar_cond_beta_" + string(beta) + ".eps", "epsc");
end

fprintf("Figures saved to ./figure/\n");
