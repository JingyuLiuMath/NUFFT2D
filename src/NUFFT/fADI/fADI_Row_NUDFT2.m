function [row_sk, U, k] = fADI_Row_NUDFT2(gamma, M, k, alpha, beta)
% fADI_Row_NUDFT2 A * X - X * B = M * N'.

m = length(gamma);
rho = size(M, 2);

W = zeros(m, k * rho);

w = M ./ (gamma - beta(1));
W(:, 1) = w;
d = zeros(k * rho, 1);
d(1) = (beta(1) - alpha(1));
for j = 1 : (k - 1)
    w = (gamma - alpha(j)) .* W(:, j);
    w = w ./ (gamma - beta(j + 1));
    W(:, j + 1) = w;
    d(j + 1) = beta(j + 1) - alpha(j + 1);
end

[row_sk, U, k, ~] = LowRank_Row_ID(W .* d.', k);

end