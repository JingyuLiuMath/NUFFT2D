clear;
close all;

p_list = (5 : 9)';
num_n = length(p_list);

beta_list = [0.7, 0.6]';
num_beta = length(beta_list);

tol_hss_list = [1e-2, 1e-4];
num_tol = length(tol_hss_list);

n_list = 2.^p_list;
N_list = n_list.^2;

M_list = zeros(num_n, num_beta);
t_construct_list = zeros(num_n, num_tol, num_beta);
t_factor_list    = zeros(num_n, num_tol, num_beta);
t_direct_list    = zeros(num_n, num_tol, num_beta);
rel_res_direct_list = zeros(num_n, num_tol, num_beta);
rel_acc_direct_list = zeros(num_n, num_tol, num_beta);

for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    
    for it_n = 1 : num_n
        p = p_list(it_n);
        
        for it_tol = 1 : num_tol
            tol_hss = tol_hss_list(it_tol);
            load("./data/typeII_2d_results_" + string(p) + "_" + string(beta) + "_" + string(tol_hss) + ".mat");
            
            M_list(it_n, it_beta) = M;
            t_construct_list(it_n, it_tol, it_beta) = t_construct;
            t_factor_list(it_n, it_tol, it_beta)    = t_factor;
            t_direct_list(it_n, it_tol, it_beta)    = t_direct;
            rel_res_direct_list(it_n, it_tol, it_beta) = rel_res_direct;
            rel_acc_direct_list(it_n, it_tol, it_beta) = rel_acc_direct;
        end
    end
end

%% Print table for each beta
k = num_tol;  % Number of rows for multirow

for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    
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
        tmp_str = string(M_list(it_n, it_beta));
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
            fprintf("& %.1e ", t_construct_list(it_n, it_tol, it_beta));
            fprintf("& %.1e ", t_factor_list(it_n, it_tol, it_beta));
            fprintf("& %.1e ", t_direct_list(it_n, it_tol, it_beta));
            fprintf("& %.1e ", rel_res_direct_list(it_n, it_tol, it_beta));
            fprintf("& %.1e ", rel_acc_direct_list(it_n, it_tol, it_beta));
            fprintf("\\\\ \n");
        end
    end
    fprintf("\\bottomrule\n")
    fprintf("\\end{tabular}\n")
    if beta == 0.7
        fprintf("\\caption{Construction time, factorization time, solve time, relative residual and relative error of the direct solver on a polar grid. \\(\\beta = 0.7\\). \\(M \\approx 0.55 N \\log_{4} (N) \\).}\n")
        fprintf("\\label{tab:typeII_2d_polar_direct_0.7}\n")
    elseif beta == 0.6
        fprintf("\\caption{Construction time, factorization time, solve time, relative residual and relative error of the direct solver on a polar grid. \\(\\beta = 0.6\\). \\(M \\approx 0.47 N \\log_{4} (N) \\).}\n")
        fprintf("\\label{tab:typeII_2d_polar_direct_0.6}\n")
    end
    fprintf("\\end{table}\n")
    fprintf("\n\n")
end
