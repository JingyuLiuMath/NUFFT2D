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

M_list = zeros(num_n, 1);

t_pre_list_lsqr = zeros(num_n, 1);
t_iter_list_lsqr = zeros(num_n, 1);
n_iter_list_lsqr = zeros(num_n, 1);
rel_res_list_lsqr = zeros(num_n, 1);
rel_acc_list_lsqr = zeros(num_n, 1);

t_pre_list_plsqr = zeros(num_n, num_tol);
t_iter_list_plsqr = zeros(num_n, num_tol);
n_iter_list_plsqr = zeros(num_n, num_tol);
rel_res_list_plsqr = zeros(num_n, num_tol);
rel_acc_list_plsqr = zeros(num_n, num_tol);

for it_n = 1 : num_n
    p = p_list(it_n);

    load("./typeII_2d_rand_lsqr/data/typeII_2d_results_" + string(p) + ".mat");

    M_list(it_n) = M;

    t_iter_list_lsqr(it_n) = t_iter;
    n_iter_list_lsqr(it_n) = iter;
    rel_res_list_lsqr(it_n) = rel_res;
    rel_acc_list_lsqr(it_n) = rel_acc;

    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        load("./typeII_2d_rand_plsqr/data/typeII_2d_results_" + string(p) + "_tol_" + string(tol_hss) + ".mat");

        t_pre_list_plsqr(it_n, it_tol) = t_pre;
        t_iter_list_plsqr(it_n, it_tol) = t_iter;
        n_iter_list_plsqr(it_n, it_tol) = iter;
        rel_res_list_plsqr(it_n, it_tol) = rel_res_iter;
        rel_acc_list_plsqr(it_n, it_tol) = rel_acc_iter;
    end
end

[M_list ./ N_list]

fprintf("\n\n")

fprintf("\\begin{table}[tbhp]\n")
fprintf("\\centering\n")
fprintf("\\begin{tabular}{c c ccccc}\n")
fprintf("\\toprule\n")
fprintf("$N$ & Method & $t_{\\pre}$ & $t_{\\iter}$ & $n_{\\iter}$ & $r_{\\solve}$ & $e_{\\solve}$ \\\\ \n")
for it_n = 1 : num_n
    fprintf("\\midrule\n")
    n = n_list(it_n);
    tmp_str = string(n) + "^2";
    fprintf("\\multirow{3}{*}{$%s$} ", tmp_str);
    fprintf("& LSQR ");
    fprintf("& - ");
    fprintf("& %.1e ", t_iter_list_lsqr(it_n));
    fprintf("& %d ", n_iter_list_lsqr(it_n));
    fprintf("& %.1e ", rel_res_list_lsqr(it_n));
    fprintf("& %.1e ", rel_acc_list_lsqr(it_n));
    fprintf("\\\\ \n");
    for it_tol = 1 : num_tol
        tol_hss = tol_hss_list(it_tol);
        fprintf(" ");
        if tol_hss == 1e-2
            display_name = "tol = 1e-2";
        elseif tol_hss == 1e-4
            display_name = "tol = 1e-4";
        end
        fprintf("& PLSQR (" + display_name + ") ");
        fprintf("& %.1e ", t_pre_list_plsqr(it_n, it_tol));
        fprintf("& %.1e ", t_iter_list_plsqr(it_n, it_tol));
        fprintf("& %d ", n_iter_list_plsqr(it_n, it_tol));
        fprintf("& %.1e ", rel_res_list_plsqr(it_n, it_tol));
        fprintf("& %.1e ", rel_acc_list_plsqr(it_n, it_tol));
        fprintf("\\\\ \n");
    end
end
fprintf("\\bottomrule\n")
fprintf("\\end{tabular}\n")
fprintf("\\caption{Time cost, iteration number, relative residual and relative error of LSQR and PLSQR for 2D type-II NUDFT with a random grid. \\(M = 2 N\\).}\n")
fprintf("\\label{tab:typeII_2d_rand_iterative}\n")
fprintf("\\end{table}\n")
fprintf("\n\n")
