clear;
close all;

p_list = [9]';
num_n = length(p_list);

beta_list = [0.6]';
num_beta = length(beta_list);

tol_hss_list = [1e-4];
num_tol = length(tol_hss_list);


for it_beta = 1 : num_beta
    beta = beta_list(it_beta);
    
    for it_n = 1 : num_n
        p = p_list(it_n);

        load("./typeII_2d_polar_lsqr/data/typeII_2d_results_" + string(p) + "_" + string(beta) + ".mat");

        for it_tol = 1 : num_tol
            tol_hss = tol_hss_list(it_tol);
            load("./typeII_2d_polar_plsqr/data/typeII_2d_results_" + string(p) + "_" + string(beta) + "_tol_" + string(tol_hss) + ".mat");

            norm_P = norm(P, "fro");
            fprintf("rel err of direct: %.1e\n", norm(P_reconstruct_direct - P, "fro") / norm_P);
            fprintf("rel err of lsqr: %.1e\n", norm(P_reconstruct - P, "fro") / norm_P);
            fprintf("rel err of plsqr: %.1e\n", norm(P_reconstruct_iter - P, "fro") / norm_P);

            figure;
            subplot(2, 2, 1);
            imagesc(P);
            colormap gray;
            axis image;
            axis off;
            title('Original', 'FontSize', 18);
            colorbar;
            subplot(2, 2, 2);
            imagesc(P_reconstruct_direct);
            colormap gray;
            axis image;
            axis off;
            title('Direct Solver', 'FontSize', 18);
            colorbar;
            subplot(2, 2, 3);
            imagesc(P_reconstruct);
            colormap gray;
            axis image;
            axis off;
            title('LSQR', 'FontSize', 18);
            colorbar;
            subplot(2, 2, 4);
            imagesc(P_reconstruct_iter);
            colormap gray;
            axis image;
            axis off;
            title('PLSQR', 'FontSize', 18);
            colorbar;

            saveas(gcf, "./figure/polar_image_" + string(p) + "_" + string(beta) + "_tol_" + string(tol_hss) + ".png", "png");
            saveas(gcf, "./figure/polar_image_" + string(p) + "_" + string(beta) + "_tol_" + string(tol_hss) + ".eps", "epsc");
        end
    end
end
