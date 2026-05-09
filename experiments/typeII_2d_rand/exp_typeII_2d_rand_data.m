exp_typeII_2d_rand_settings;

for it_alpha = 1 : num_alpha
    for it_p = 1 : num_n
        alpha = alpha_list(it_alpha);
        p = p_list(it_p);

        n = 2^p;
        nx = n;
        ny = n;
        N = nx * ny;
        M = round(alpha * N);
        xy = rand(M, 2);

        P = phantom('Modified Shepp-Logan', n);
        c_ex = reshape(P, N, []);
        f_ex = MY_NUFFT2_2D(c_ex, xy, nx, ny);

        save("./data/typeII_2d_points" ...
            + "_" + string(p) ...
            + "_" + string(alpha) ...
            + ".mat", "xy", "c_ex", "f_ex", "n");
    end
end