clear;
close all;

p_list = (4 : 9)';
num_n = length(p_list);

n_list = 2.^p_list;
N_list = n_list.^2;

hss_rank_list = zeros(num_n, 1);
mem_list = zeros(num_n, 1);
mem_ratio_list = zeros(num_n, 1);

t_construct_list = zeros(num_n, 1);
t_apply_list = zeros(num_n, 1);
t_nufft_list = zeros(num_n, 1);
t_factor_list = zeros(num_n, 1);
t_solve_list = zeros(num_n, 1);
t_total_list = zeros(num_n, 1);

rel_err_list = zeros(num_n, 1);
rel_res_list = zeros(num_n, 1);

for it = 1 : num_n
    p = p_list(it);

    load("./data/typeII_2d_results_" + string(p) + ".mat");

    hss_rank_list(it) = hss_rank;
    mem_list(it) = mem;
    mem_ratio_list(it) = mem_ratio;

    t_construct_list(it) = t_construct;
    t_apply_list(it) = t_apply;
    t_nufft_list(it) = t_nufft;
    t_factor_list(it) = t_factor;
    t_solve_list(it) = t_solve;
    t_total_list(it) = t_total;

    rel_err_list(it) = rel_err;
    rel_res_list(it) = rel_res;
end
