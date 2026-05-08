exp_typeII_2d_polar_settings;

for it_beta = 1 : num_beta
    for it_p = 1 : num_n
        beta = beta_list(it_beta);
        p = p_list(it_p);

        n = 2^p;
        nx = n;
        ny = n;
        N = nx * ny;
        x = PolarGrid(n, 1, beta * p);
        M = size(x, 1);

        save("./data/typeII_2d_points" ...
            + "_" + string(p) ...
            + "_" + string(beta) ...
            + ".mat", "x");
    end
end