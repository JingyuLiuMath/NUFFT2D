clear;
close all;

p_list = (5 : 9)';
num_n = length(p_list);

alpha_list = [1.5, 2.0]';
num_alpha = length(alpha_list);

tol_hss_list = [1e-2, 1e-4]';
num_tol = length(tol_hss_list);

n_list = 2.^p_list;
N_list = n_list.^2;

M_list = zeros(num_n, num_tol, num_alpha);

hss_rank_list = zeros(num_n, num_tol, num_alpha);
mem_list = zeros(num_n, num_tol, num_alpha);

t_construct_list = zeros(num_n, num_tol, num_alpha);
t_factor_list = zeros(num_n, num_tol, num_alpha);
t_solve_list = zeros(num_n, num_tol, num_alpha);

rel_res_list = zeros(num_n, num_tol, num_alpha);
rel_acc_list = zeros(num_n, num_tol, num_alpha);

cond_est_list = zeros(num_n, num_tol, num_alpha);
cond_exact_list = zeros(num_n, num_tol, num_alpha);
cond_exact_2norm_list = zeros(num_n, num_tol, num_alpha);

%% Load all data
for it_alpha = 1 : num_alpha
    alpha = alpha_list(it_alpha);
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        fprintf("alpha: %.1f, tol_hss: %.1e\n", alpha, tol_hss);
        for it_n = 1 : num_n
            p = p_list(it_n);

            load("./data/typeII_2d_results_" + string(p) + "_" + string(alpha) + "_tol_" + string(tol_hss) + ".mat");

            M_list(it_n, it_tol, it_alpha) = M;

            hss_rank_list(it_n, it_tol, it_alpha) = hss_rank;
            mem_list(it_n, it_tol, it_alpha) = double_to_gb(mem);

            t_construct_list(it_n, it_tol, it_alpha) = t_construct;
            t_factor_list(it_n, it_tol, it_alpha) = t_factor;
            t_solve_list(it_n, it_tol, it_alpha) = t_direct;

            rel_res_list(it_n, it_tol, it_alpha) = rel_res_direct;
            rel_acc_list(it_n, it_tol, it_alpha) = rel_acc_direct;

            cond_est_list(it_n, it_tol, it_alpha) = cond_number_est;
            cond_exact_list(it_n, it_tol, it_alpha) = cond_number_exact;
            cond_exact_2norm_list(it_n, it_tol, it_alpha) = cond_number_exact_2norm;
        end
    end
end

%% Plot for each alpha
markers = ["x", "+"];
colors = ["#0072BD", "#D95319"];
ref_color = "#77AC30";

for it_alpha = 1 : num_alpha
    alpha = alpha_list(it_alpha);

    % HSS Rank
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "tol = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "tol = 1e-4";
        end
        loglog(N_list, hss_rank_list(:, it_tol, it_alpha), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    % Reference line: O(sqrt(N) * log(N))
    ref_line = sqrt(N_list) .* log2(N_list);
    ref_line = ref_line / ref_line(1) * hss_rank_list(1, 2, it_alpha) * 1.1;
    loglog(N_list, ref_line, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Color", ref_color, ...
        "DisplayName", "$O(\sqrt{N} \log N)$");
    if alpha == 2.0
        title(sprintf("HSS Rank ($M = 2 N$)"), "Interpreter", "latex");
    elseif alpha == 1.5
        title(sprintf("HSS Rank ($M = 1.5 N$)"), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$r_{\mathrm{h}}$", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/rand_hss_rank_alpha_" + string(alpha) + ".png", "png");
    saveas(gcf, "./figure/rand_hss_rank_alpha_" + string(alpha) + ".eps", "epsc");

    % Memory
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "tol = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "tol = 1e-4";
        end
        loglog(N_list, mem_list(:, it_tol, it_alpha), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    % Reference line: O(M * sqrt(N) * log(N))
    % ref_line = M_list(:, 1, it_alpha) .* sqrt(N_list) .* log2(N_list);
    ref_line = N_list .* log2(N_list).^3;
    ref_line = ref_line / ref_line(1) * mem_list(1, 2, it_alpha) * 1.1;
    loglog(N_list, ref_line, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Color", ref_color, ...
        "DisplayName", "$O(N \log^{3} N)$");
    if alpha == 2.0
        title(sprintf("Memory ($M = 2 N$)"), "Interpreter", "latex");
    elseif alpha == 1.5
        title(sprintf("Memory ($M = 1.5 N$)"), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$m_{\mathrm{h}}$ (GB)", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/rand_memory_alpha_" + string(alpha) + ".png", "png");
    saveas(gcf, "./figure/rand_memory_alpha_" + string(alpha) + ".eps", "epsc");

    % Construct Time
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "tol = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "tol = 1e-4";
        end
        loglog(N_list, t_construct_list(:, it_tol, it_alpha), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    % Reference line: O(M * N * log^2(N))
    % ref_line = M_list(:, 1, it_alpha) .* N_list .* log2(N_list).^2;
    ref_line = N_list.^(1.5) .* log2(N_list).^3;
    ref_line = ref_line / ref_line(1) * t_construct_list(1, 2, it_alpha) * 1.1;
    loglog(N_list, ref_line, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Color", ref_color, ...
        "DisplayName", "$O(N^{1.5} \log^{3} N)$");
    if alpha == 2.0
        title(sprintf("Construct Time ($M = 2 N$)"), "Interpreter", "latex");
    elseif alpha == 1.5
        title(sprintf("Construct Time ($M = 1.5 N$)"), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$t_{\mathrm{c}}$ (s)", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/rand_t_construct_alpha_" + string(alpha) + ".png", "png");
    saveas(gcf, "./figure/rand_t_construct_alpha_" + string(alpha) + ".eps", "epsc");

    % Factor Time
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "tol = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "tol = 1e-4";
        end
        loglog(N_list, t_factor_list(:, it_tol, it_alpha), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    % Reference line: O(M * N * log^2(N))
    % ref_line = M_list(:, 1, it_alpha) .* N_list .* log2(N_list).^2;
    ref_line = N_list.^(1.5) .* log2(N_list).^3;
    ref_line = ref_line / ref_line(1) * t_factor_list(1, 2, it_alpha) * 1.1;
    loglog(N_list, ref_line, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Color", ref_color, ...
        "DisplayName", "$O(N^{1.5} \log^{3} N)$");
    if alpha == 2.0
        title(sprintf("Factor Time ($M = 2 N$)"), "Interpreter", "latex");
    elseif alpha == 1.5
        title(sprintf("Factor Time ($M = 1.5 N$)"), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$t_{\mathrm{f}}$ (s)", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/rand_t_factor_alpha_" + string(alpha) + ".png", "png");
    saveas(gcf, "./figure/rand_t_factor_alpha_" + string(alpha) + ".eps", "epsc");

    % Solve Time
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "tol = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "tol = 1e-4";
        end
        loglog(N_list, t_solve_list(:, it_tol, it_alpha), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    % Reference line: O(M * sqrt(N) * log(N))
    % ref_line = M_list(:, 1, it_alpha) .* sqrt(N_list) .* log2(N_list);
    ref_line = N_list .* log2(N_list).^3;
    ref_line = ref_line / ref_line(1) * t_solve_list(1, 2, it_alpha) * 1.1;
    loglog(N_list, ref_line, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Color", ref_color, ...
        "DisplayName", "$O(N \log^{3} N)$");
    if alpha == 2.0
        title(sprintf("Solve Time ($M = 2 N$)"), "Interpreter", "latex");
    elseif alpha == 1.5
        title(sprintf("Solve Time ($M = 1.5 N$)"), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$t_{\mathrm{s}}$ (s)", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/rand_t_solve_alpha_" + string(alpha) + ".png", "png");
    saveas(gcf, "./figure/rand_t_solve_alpha_" + string(alpha) + ".eps", "epsc");

    % Relative Residual
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "tol = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "tol = 1e-4";
        end
        loglog(N_list, rel_res_list(:, it_tol, it_alpha), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    if alpha == 2.0
        title(sprintf("Relative Residual ($M = 2 N$)"), "Interpreter", "latex");
    elseif alpha == 1.5
        title(sprintf("Relative Residual ($M = 1.5 N$)"), "Interpreter", "latex");
    end
    title(sprintf("Relative Residual ($M \\approx %.1f N$)", alpha), "Interpreter", "latex");
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$r_{\mathrm{s}}$", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/rand_rel_res_alpha_" + string(alpha) + ".png", "png");
    saveas(gcf, "./figure/rand_rel_res_alpha_" + string(alpha) + ".eps", "epsc");

    % Relative Accuracy
    figure();
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        if tol_hss == 1e-2
            display_name = "tol = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "tol = 1e-4";
        end
        loglog(N_list, rel_acc_list(:, it_tol, it_alpha), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", colors(it_tol), ...
            "DisplayName", display_name);
        hold on;
    end
    if alpha == 2.0
        title(sprintf("Relative Accuracy ($M = 2 N$)"), "Interpreter", "latex");
    elseif alpha == 1.5
        title(sprintf("Relative Accuracy ($M = 1.5 N$)"), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$e_{\mathrm{s}}$", "Interpreter", "latex");
    lgd = legend("Location", "southeast");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    saveas(gcf, "./figure/rand_rel_acc_alpha_" + string(alpha) + ".png", "png");
    saveas(gcf, "./figure/rand_rel_acc_alpha_" + string(alpha) + ".eps", "epsc");

    % Condition Number (only tol_hss == 1e-4)
    figure();
    it_tol = find(tol_hss_list == 1e-4, 1);

    valid_ind = ~isnan(cond_est_list(:, it_tol, it_alpha));
    if any(valid_ind)
        loglog(N_list(valid_ind), cond_est_list(valid_ind, it_tol, it_alpha), ...
            "LineWidth", 2, ...
            "Marker", markers(it_tol), ...
            "Color", "#0072BD", ...
            "DisplayName", "Estimated cond (1-norm)");
        hold on;
    end

    valid_ind = ~isnan(cond_exact_2norm_list(:, it_tol, it_alpha));
    if any(valid_ind)
        loglog(N_list(valid_ind), cond_exact_2norm_list(valid_ind, it_tol, it_alpha), ...
            "LineWidth", 2, ...
            "Marker", "+", ...
            "LineStyle", "--", ...
            "Color", "#D95319", ...
            "DisplayName", "Exact cond (2-norm)");
        hold on;
    end

    valid_ind = ~isnan(cond_exact_list(:, it_tol, it_alpha));
    if any(valid_ind)
        loglog(N_list(valid_ind), cond_exact_list(valid_ind, it_tol, it_alpha), ...
            "LineWidth", 2, ...
            "Marker", "o", ...
            "LineStyle", ":", ...
            "Color", "#77AC30", ...
            "DisplayName", "Exact cond (1-norm)");
        hold on;
    end
    if alpha == 2.0
        title(sprintf("Condition Number ($M = 2 N$)"), "Interpreter", "latex");
    elseif alpha == 1.5
        title(sprintf("Condition Number ($M = 1.5 N$)"), "Interpreter", "latex");
    end
    xlabel("$N$", "Interpreter", "latex");
    ylabel("$\kappa$", "Interpreter", "latex");
    lgd = legend("Location", "northwest");
    lgd.Interpreter = 'latex';
    set(gca, 'FontSize', 18);
    % grid on;
    hold off;
    saveas(gcf, "./figure/rand_cond_alpha_" + string(alpha) + ".png", "png");
    saveas(gcf, "./figure/rand_cond_alpha_" + string(alpha) + ".eps", "epsc");
end

fprintf("Figures saved to ./figure/\n");
