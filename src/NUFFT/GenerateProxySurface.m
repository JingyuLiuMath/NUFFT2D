function x_proxy = GenerateProxySurface(N, x_start, x_end, proxy_size)
% GenerateProxySurface Generate proxy surface points.

% Jingyu Liu, November 19, 2024.

arguments (Input)
    N (1, 1) double;
    x_start (1, 1) double;
    x_end (1, 1) double;
    proxy_size (1, 1) double;
end

arguments (Output)
    x_proxy (:, 1) double;
end

tmp_x = 1 / N / 4 + (0 : (proxy_size - 1))' / N;
tmp_x_1 = tmp_x + (2 * rand(proxy_size, 1) - 1) / N / 8;
tmp_x_2 = tmp_x + (2 * rand(proxy_size, 1) - 1) / N / 8;

x_start_proxy = x_start - tmp_x_1;
x_end_proxy = x_end + tmp_x_2;
x_proxy = [x_start_proxy; x_end_proxy];

end