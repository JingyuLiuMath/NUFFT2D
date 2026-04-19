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

t_pre_list_lsqr = zeros(num_n, num_beta);
t_iter_list_lsqr = zeros(num_n, num_beta);
n_iter_list_lsqr = zeros(num_n, num_beta);
rel_res_list_lsqr = zeros(num_n, num_beta);
rel_acc_list_lsqr = zeros(num_n, num_beta);

t_pre_list_plsqr = zeros(num_n, num_tol, num_beta);
t_iter_list_plsqr = zeros(num_n, num_tol, num_beta);
n_iter_list_plsqr = zeros(num_n, num_tol, num_beta);
rel_res_list_plsqr = zeros(num_n, num_tol, num_beta);
rel_acc_list_plsqr = zeros(num_n, num_tol, num_beta);

beta_av_list = zeros(num_beta, 1);
for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    
    for it_n = 1 : num_n
        p = p_list(it_n);

        load("./typeII_2d_polar_lsqr/data/typeII_2d_results_" + string(p) + "_" + string(beta) + ".mat");

        M_list(it_n, it_beta) = M;

        t_iter_list_lsqr(it_n, it_beta) = t_iter;
        n_iter_list_lsqr(it_n, it_beta) = iter;
        rel_res_list_lsqr(it_n, it_beta) = rel_res;
        rel_acc_list_lsqr(it_n, it_beta) = rel_acc;

        for it_tol = 1 : num_tol
            tol_hss = tol_hss_list(it_tol);
            load("./typeII_2d_polar_plsqr/data/typeII_2d_results_" + string(p) + "_" + string(beta) + "_" + string(tol_hss) + ".mat");

            assert(M == M_list(it_n, it_beta), "M does not match for LSQR and PLSQR.");

            t_pre_list_plsqr(it_n, it_tol, it_beta) = t_pre;
            t_iter_list_plsqr(it_n, it_tol, it_beta) = t_iter;
            n_iter_list_plsqr(it_n, it_tol, it_beta) = iter;
            rel_res_list_plsqr(it_n, it_tol, it_beta) = rel_res_iter;
            rel_acc_list_plsqr(it_n, it_tol, it_beta) = rel_acc_iter;
        end
    end
end

for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    fprintf("beta = %.1e: ", beta);
    beta_av = mean(M_list(:, it_beta) ./ N_list ./ p_list);
    beta_av_list(it_beta) = beta_av;
    fprintf("beta_av = %.1e: ", beta_av);
    disp(M_list(:, it_beta) ./ N_list ./ p_list);
end

%% Print table for each beta
k = num_tol + 1;  % Number of rows for multirow

for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    beta_av = beta_av_list(it_beta);
    
    fprintf("\n\n")
    
    fprintf("\\begin{table}[tbhp]\n")
    fprintf("\\centering\n")
    fprintf("\\begin{tabular}{cc c ccccc}\n")
    fprintf("\\toprule\n")
    fprintf("\\(N\\) & \\(M\\) & Method & \\(t_{\\pre}\\) & \\(t_{\\iter}\\) & \\(n_{\\iter}\\) & \\(r_{\\solve}\\) & \\(e_{\\solve}\\) \\\\ \n")
    for it_n = 1 : num_n
        fprintf("\\midrule\n")
        n = n_list(it_n);
        tmp_str = string(n) + "^2";
        fprintf("\\multirow{%d}{*}{\\(%s\\)} ", k, tmp_str);
        tmp_str = string(M_list(it_n, it_beta));
        fprintf("& \\multirow{%d}{*}{\\(%s\\)} ", k, tmp_str);
        fprintf("& LSQR ");
        fprintf("& - ");
        fprintf("& %.1e ", t_iter_list_lsqr(it_n, it_beta));
        fprintf("& %d ", n_iter_list_lsqr(it_n, it_beta));
        fprintf("& %.1e ", rel_res_list_lsqr(it_n, it_beta));
        fprintf("& %.1e ", rel_acc_list_lsqr(it_n, it_beta));
        fprintf("\\\\ \n");
        for it_tol = 1 : num_tol
            tol_hss = tol_hss_list(it_tol);
            fprintf(" ");
            if tol_hss == 1e-2
                display_name = "\\(\\varepsilon\\) = 1e-2";
            elseif tol_hss == 1e-4
                display_name = "\\(\\varepsilon\\) = 1e-4";
            end
            fprintf("& & PLSQR (" + display_name + ") ");
            fprintf("& %.1e ", t_pre_list_plsqr(it_n, it_tol, it_beta));
            fprintf("& %.1e ", t_iter_list_plsqr(it_n, it_tol, it_beta));
            fprintf("& %d ", n_iter_list_plsqr(it_n, it_tol, it_beta));
            fprintf("& %.1e ", rel_res_list_plsqr(it_n, it_tol, it_beta));
            fprintf("& %.1e ", rel_acc_list_plsqr(it_n, it_tol, it_beta));
            fprintf("\\\\ \n");
        end
    end
    fprintf("\\bottomrule\n")
    fprintf("\\end{tabular}\n")
    if beta == 0.7
        fprintf("\\caption{Time cost, iteration number, relative residual and relative error of LSQR and PLSQR on a polar grid. \\(\\beta = 0.7\\). \\(M \\approx 0.55 N \\log_{4} (N)\\).}\n", beta_av)
        fprintf("\\label{tab:typeII_2d_polar_iterative_0.7}\n")
    elseif beta == 0.6
        fprintf("\\caption{Time cost, iteration number, relative residual and relative error of LSQR and PLSQR on a polar grid. \\(\\beta = 0.6\\). \\(M \\approx 0.47 N \\log_{4} (N)\\).}\n", beta_av)
        fprintf("\\label{tab:typeII_2d_polar_iterative_0.6}\n")
    end
    fprintf("\\end{table}\n")
    fprintf("\n\n")
end
