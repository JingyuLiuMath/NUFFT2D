clear;
close all;
warning off;

p_list = (5 : 9)';
num_n = length(p_list);

alpha_list = [2.0; 1.5];
alpha_string_list = ["2.0"; "1.5"];
alpha_display_list = ["2"; "1.5"];
num_alpha = length(alpha_list);

tol_hss_list = [1e-3; 1e-6];
tol_hss_display_list = ["10^{-3}"; "10^{-6}"];
num_tol_hss = length(tol_hss_list);

min_points = 256;

tol_cg = 1e-12;
maxit_cg = 500;

figure_prefix = "./figure/rand";
