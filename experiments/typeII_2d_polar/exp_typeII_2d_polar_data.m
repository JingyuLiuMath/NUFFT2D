exp_typeII_2d_polar_settings;

for it_beta = 1 : num_beta
    for it_p = 1 : num_n
        beta = beta_list(it_beta);
        p = p_list(it_p);

        n = 2^p;
        nx = n;
        ny = n;
        N = nx * ny;
        xy = PolarGrid(n, 1, beta * p);
        M = size(xy, 1);

        P = phantom('Modified Shepp-Logan', n);
        c_ex = reshape(P, N, []);
        f_ex = MY_NUFFT2_2D(c_ex, xy, nx, ny);

        save("./data/typeII_2d_points" ...
            + "_" + string(p) ...
            + "_" + string(beta) ...
            + ".mat", "xy", "c_ex", "f_ex", "n");
    end
end