function k_z_w = NUFFT2_Kernel(z, w, N)

arguments (Input)
    z (:, 1) double;
    w (:, 1) double;
    N (1, 1) double;
end

arguments (Output)
    k_z_w (:, :) double;
end

phi = z ./ w.';
tmp_ind = abs(phi(:) - 1) <= eps;
k_z_w = (phi.^N - 1) ./ (phi - 1) / N;
k_z_w(tmp_ind) = 1;

end