clear;
close all;

p_list = (5 : 9)';
num_n = length(p_list);

alpha_list = [2.0, 1.5]';
num_alpha = length(alpha_list);

tol_hss_list = [1e-2, 1e-4];
num_tol = length(tol_hss_list);

n_list = 2.^p_list;
N_list = n_list.^2;

M_list = zeros(num_n, num_alpha);
t_construct_list = zeros(num_n, num_tol, num_alpha);
t_factor_list    = zeros(num_n, num_tol, num_alpha);
t_direct_list    = zeros(num_n, num_tol, num_alpha);
rel_res_direct_list = zeros(num_n, num_tol, num_alpha);
rel_acc_direct_list = zeros(num_n, num_tol, num_alpha);

for it_alpha = 1 : num_alpha
    alpha = alpha_list(it_alpha);
    
    for it_n = 1 : num_n
        p = p_list(it_n);
    
        for it_tol = 1 : num_tol
            tol_hss = tol_hss_list(it_tol);
            load("./data/typeII_2d_results_" + string(p) + "_" + string(alpha) + "_" + string(tol_hss) + ".mat");
            
            M_list(it_n, it_alpha) = M;
            t_construct_list(it_n, it_tol, it_alpha) = t_construct;
            t_factor_list(it_n, it_tol, it_alpha)    = t_factor;
            t_direct_list(it_n, it_tol, it_alpha)    = t_direct;
            rel_res_direct_list(it_n, it_tol, it_alpha) = rel_res_direct;
            rel_acc_direct_list(it_n, it_tol, it_alpha) = rel_acc_direct;
        end
    end
end

%% Print table for each alpha
k = num_tol;  % Number of rows for multirow

for it_alpha = 1 : num_alpha
    alpha = alpha_list(it_alpha);
    
    fprintf("\n\n")
    
    fprintf("\\begin{table}[tbhp]\n")
    fprintf("\\centering\n")
    fprintf("\\begin{tabular}{cc c ccccc}\n")
    fprintf("\\toprule\n")
    fprintf("\\(N\\) & \\(M\\) & \\(\\varepsilon\\) & \\(t_{\\mathrm{c}}\\) & \\(t_{\\mathrm{f}}\\) & \\(t_{\\mathrm{s}}\\) & \\(r_{\\mathrm{s}}\\) & \\(e_{\\mathrm{s}}\\) \\\\ \n")
    for it_n = 1 : num_n
        fprintf("\\midrule\n")
        n = n_list(it_n);
        tmp_str = string(n) + "^2";
        fprintf("\\multirow{%d}{*}{\\(%s\\)} ", k, tmp_str);
        tmp_str = string(M_list(it_n, it_alpha));
        fprintf("& \\multirow{%d}{*}{\\(%s\\)} ", k, tmp_str);
        for it_tol = 1 : num_tol
            tol_hss = tol_hss_list(it_tol);
            if tol_hss == 1e-2
                display_eps = "1e-2";
            elseif tol_hss == 1e-4
                display_eps = "1e-4";
            else
                display_eps = string(tol_hss);
            end
            if it_tol == 1
                fprintf(" & %s ", display_eps);
            else
                fprintf(" & & %s ", display_eps);
            end
            fprintf("& %.1e ", t_construct_list(it_n, it_tol, it_alpha));
            fprintf("& %.1e ", t_factor_list(it_n, it_tol, it_alpha));
            fprintf("& %.1e ", t_direct_list(it_n, it_tol, it_alpha));
            fprintf("& %.1e ", rel_res_direct_list(it_n, it_tol, it_alpha));
            fprintf("& %.1e ", rel_acc_direct_list(it_n, it_tol, it_alpha));
            fprintf("\\\\ \n");
        end
    end
    fprintf("\\bottomrule\n")
    fprintf("\\end{tabular}\n")
    if alpha == 1.5
        fprintf("\\caption{Construction time, factorization time, solve time, relative residual and relative error of the direct solver on a random grid. \\(\\alpha = 1.5 \\). \\(M = 1.5 N\\).}\n")
        fprintf("\\label{tab:typeII_2d_rand_direct_1.5}\n")
    elseif alpha == 2.0
        fprintf("\\caption{Construction time, factorization time, solve time, relative residual and relative error of the direct solver on a random grid. \\(\\alpha = 2\\). \\(M = 2 N\\).}\n")
        fprintf("\\label{tab:typeII_2d_rand_direct_2}\n")
    end
    fprintf("\\end{table}\n")
    fprintf("\n\n")
end
