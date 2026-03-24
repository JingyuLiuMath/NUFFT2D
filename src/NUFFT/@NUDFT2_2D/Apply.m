function f = Apply(A, u)
% Apply

% Jingyu Liu, December 8, 2024.

arguments (Input)
    A NUDFT2_2D;
    u (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

n1 = A.nx_;
n2 = A.ny_;
N = A.N_;

u = reshape(u, n1, n2, []);
u = fft2(u);
u = reshape(u, N, []);
u = u(A.omega_perm_, :);
f = A.AFinv_HSS_.Apply(u);
f = f(A.xy_inv_perm_, :);

end