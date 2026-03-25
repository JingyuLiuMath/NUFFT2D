function k_z_w = NUFFT2_Kernel(z, w, N)
% NUFFT2_Kernel

arguments (Input)
    z (:, 1) double;
    w (:, 1) double;
    N (1, 1) double;
end

arguments (Output)
    k_z_w (:, :) double;
end

% xi = z * w';
% theta = angle(xi);
% k_z_w = 1 ./ (xi - 1);
% [row_ind, col_ind] = find(abs(1 - xi) < 1e-9);
% sub_theta = theta(row_ind, col_ind);
% k_z_w(row_ind, col_ind) = 1 ./ expm1(sub_theta * 1i);
% k_z_w = ((z.^N - 1) .* k_z_w) / N;

theta_z = angle(z);
xi = z * w';
theta = angle(xi);
k_z_w = expm1(theta_z * N * 1i) ./ expm1(theta * 1i) / N;

% xi = z * w';
% theta = angle(xi);
% k_z_w = expm1(theta * N * 1i) ./ expm1(theta * 1i) / N;
% k_z_w = k_z_w .* (w.^N).';

end