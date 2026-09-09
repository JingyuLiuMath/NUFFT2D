clear;
close all;
warning off;

if ispc
    p_list = (4 : 7)';
else
    p_list = (5 : 9)';
end
num_n = length(p_list);

alpha_list = [1.5];
alpha_string_list = ["1.5"];
alpha_display_list = ["1.5"];
num_alpha = length(alpha_list);

if ispc
    tol_hss_list = [1e-3];
    tol_hss_display_list = ["10^{-3}"];
else
    tol_hss_list = [1e-2; 1e-4];
    tol_hss_display_list = ["10^{-2}"; "10^{-4}"];
end
num_tol_hss = length(tol_hss_list);

if ispc
    n_leaf = 8;
    N_leaf = 64;
else
    n_leaf = 16;
    N_leaf = 256;
end

if ispc
    tol_cg = 1e-8;
    maxit_cg = 50;
else
    tol_cg = 1e-12;
    maxit_cg = 500;
end

figure_prefix = "./figure/rand";
