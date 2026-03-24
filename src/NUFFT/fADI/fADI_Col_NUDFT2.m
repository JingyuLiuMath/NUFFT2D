function [col_sk, V, k] = fADI_Col_NUDFT2(hat_c, N, k, alpha, beta)
% fADI_Col_NUDFT2 A * X - X * B = M * N'.

n = length(hat_c);
rho = size(N, 2);

Y = zeros(n, k * rho);
conj_hat_c = conj(hat_c);
conj_alpha = conj(alpha);
conj_beta = conj(beta);

y = N ./ (conj_hat_c - conj_alpha(1));
Y(:, 1) = y;
d = zeros(k * rho, 1);
d(1) = (beta(1) - alpha(1));
for j = 1 : (k - 1)
    y = (conj_hat_c - conj_beta(j)) .* Y(:, j);
    y = y ./ (conj_hat_c - conj_alpha(j + 1));
    Y(:, j + 1) = y;
    d(j + 1) = beta(j + 1) - alpha(j + 1);
end

[col_sk, V, k] = LowRank_ID(d .* Y', k);

end