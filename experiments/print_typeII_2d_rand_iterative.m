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

t_pre_list_lsqr = zeros(num_n, num_alpha);
t_iter_list_lsqr = zeros(num_n, num_alpha);
n_iter_list_lsqr = zeros(num_n, num_alpha);
rel_res_list_lsqr = zeros(num_n, num_alpha);
rel_acc_list_lsqr = zeros(num_n, num_alpha);

t_pre_list_plsqr = zeros(num_n, num_tol, num_alpha);
t_iter_list_plsqr = zeros(num_n, num_tol, num_alpha);
n_iter_list_plsqr = zeros(num_n, num_tol, num_alpha);
rel_res_list_plsqr = zeros(num_n, num_tol, num_alpha);
rel_acc_list_plsqr = zeros(num_n, num_tol, num_alpha);

for it_alpha = 1 : num_alpha
    alpha = alpha_list(it_alpha);
    
    for it_n = 1 : num_n
        p = p_list(it_n);

        load("./typeII_2d_rand_lsqr/data/typeII_2d_results_" + string(p) + "_" + string(alpha) + ".mat");

        M_list(it_n, it_alpha) = M;

        t_iter_list_lsqr(it_n, it_alpha) = t_iter;
        n_iter_list_lsqr(it_n, it_alpha) = iter;
        rel_res_list_lsqr(it_n, it_alpha) = rel_res;
        rel_acc_list_lsqr(it_n, it_alpha) = rel_acc;

        for it_tol = 1 : num_tol
            tol_hss = tol_hss_list(it_tol);
            load("./typeII_2d_rand_plsqr/data/typeII_2d_results_" + string(p) + "_" + string(alpha) + "_" + string(tol_hss) + ".mat");

            assert(M == M_list(it_n, it_alpha), "M does not match for LSQR and PLSQR.");

            t_pre_list_plsqr(it_n, it_tol, it_alpha) = t_pre;
            t_iter_list_plsqr(it_n, it_tol, it_alpha) = t_iter;
            n_iter_list_plsqr(it_n, it_tol, it_alpha) = iter;
            rel_res_list_plsqr(it_n, it_tol, it_alpha) = rel_res_iter;
            rel_acc_list_plsqr(it_n, it_tol, it_alpha) = rel_acc_iter;
        end
    end
end

for it_alpha = 1 : num_alpha
    alpha = alpha_list(it_alpha);
    fprintf("alpha = %.1e\n", alpha);
    disp(M_list(:, it_alpha) ./ N_list);
end

%% Print table for each alpha
k = num_tol + 1;  % Number of rows for multirow

for it_alpha = 1 : num_alpha
    alpha = alpha_list(it_alpha);
    
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
        tmp_str = string(M_list(it_n, it_alpha));
        fprintf("& \\multirow{%d}{*}{\\(%s\\)} ", k, tmp_str);
        fprintf("& LSQR ");
        fprintf("& - ");
        fprintf("& %.1e ", t_iter_list_lsqr(it_n, it_alpha));
        fprintf("& %d ", n_iter_list_lsqr(it_n, it_alpha));
        fprintf("& %.1e ", rel_res_list_lsqr(it_n, it_alpha));
        fprintf("& %.1e ", rel_acc_list_lsqr(it_n, it_alpha));
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
            fprintf("& %.1e ", t_pre_list_plsqr(it_n, it_tol, it_alpha));
            fprintf("& %.1e ", t_iter_list_plsqr(it_n, it_tol, it_alpha));
            fprintf("& %d ", n_iter_list_plsqr(it_n, it_tol, it_alpha));
            fprintf("& %.1e ", rel_res_list_plsqr(it_n, it_tol, it_alpha));
            fprintf("& %.1e ", rel_acc_list_plsqr(it_n, it_tol, it_alpha));
            fprintf("\\\\ \n");
        end
    end
    fprintf("\\bottomrule\n")
    fprintf("\\end{tabular}\n")
    if alpha == 1.5
        fprintf("\\caption{Time cost, iteration number, relative residual and relative error of LSQR and PLSQR on a random grid. \\(\\alpha = 1.5 \\). \\(M = 1.5 N\\).}\n")
        fprintf("\\label{tab:typeII_2d_rand_iterative_1.5}\n")
    elseif alpha == 2.0
        fprintf("\\caption{Time cost, iteration number, relative residual and relative error of LSQR and PLSQR on a random grid. \\(\\alpha = 2\\). \\(M = 2 N\\).}\n")
        fprintf("\\label{tab:typeII_2d_rand_iterative_2}\n")
    end
    fprintf("\\end{table}\n")
    fprintf("\n\n")
end
