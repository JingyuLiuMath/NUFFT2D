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
delta = phi - 1;
near = N * abs(delta) <= 0.1;

k_z_w = zeros(size(phi), 'like', phi);

% Original formula away from coincident points.
p = phi(~near);
k_z_w(~near) = (p.^N - 1) ./ delta(~near) / N;

% Local polynomial expansion near coincident points.
d = delta(near);
term = ones(size(d), 'like', phi);
s = term;

for j = 1 : min(12, N - 1)
    term = term .* (((N - j) / (j + 1)) .* d);
    s = s + term;
end

k_z_w(near) = s;

end
