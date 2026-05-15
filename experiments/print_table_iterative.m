function print_table_iterative(...
    result_list_lsqr, ...
    result_list, ...
    caption_name, label_name, ...
    tol_hss_display_list)

num_n = size(result_list, 1);

fprintf("\n");

fprintf("\\begin{table}[tbhp]\n")
fprintf("\\centering\n")
fprintf("\\begin{tabular}{cc c ccccc}\n")
fprintf("\\toprule\n")
fprintf("\\(N\\) ");
fprintf("& \\(M\\) ");
fprintf("& method ");
fprintf("& \\(t_{\\pre}\\) ");
fprintf("& \\(t_{\\iter}\\) ");
fprintf("& \\(n_{\\iter}\\) ");
fprintf("& \\(r_{\\solve}\\) ");
fprintf("& \\(e_{\\solve}\\) ");
fprintf("\\\\ \n");

for it_n = 1 : num_n
    curr_result_lsqr = result_list_lsqr(it_n);
    curr_result = result_list(it_n, :);
    print_result_iterative(curr_result_lsqr, curr_result, tol_hss_display_list);
end
fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\caption{%s}\n", caption_name);
fprintf("\\label{tab:%s}\n", label_name);
fprintf("\\end{table}\n");
fprintf("\n");


end

function print_result_iterative(result_lsqr, result_list, tol_hss_display_list)

num_hss_tol = size(result_list, 2);
fprintf("\\midrule\n");
n = result_list(1).n;
tmp_str = string(n) + "^2";
fprintf("\\multirow{%d}{*}{\\(%s\\)} ", num_hss_tol + 1, tmp_str);
M = result_list(1).M;
fprintf("& \\multirow{%d}{*}{%d} ", num_hss_tol + 1, M);

fprintf("& lsqr ");
fprintf("& - ");
fprintf("& %.1e ", result_lsqr.t_cg);
fprintf("& %d ", result_lsqr.iter_cg);
fprintf("& %.1e ", result_lsqr.rel_res_cg);
fprintf("& %.1e ", result_lsqr.rel_err_cg);
fprintf("\\\\ \n");
for it_hss_tol = 1 : num_hss_tol
    tol_hss = result_list(it_hss_tol).tol_hss;
    tol_hss_display = tol_hss_display_list(it_hss_tol);
    fprintf(" & & plsqr (\\(\\varepsilon = %s\\)) ", tol_hss_display);
    t_pre = result_list(it_hss_tol).t_construct + result_list(it_hss_tol).t_factor;
    fprintf("& %.1e ", t_pre);
    fprintf("& %.1e ", result_list(it_hss_tol).t_pcg);
    fprintf("& %d ", result_list(it_hss_tol).iter_pcg);
    fprintf("& %.1e ", result_list(it_hss_tol).rel_res_pcg);
    fprintf("& %.1e ", result_list(it_hss_tol).rel_err_pcg);
    fprintf("\\\\ \n");
end

end