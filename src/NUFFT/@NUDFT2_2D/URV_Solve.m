function u = URV_Solve(A, f)
% URV_Solve

arguments (Input)
    A NUDFT2_2D;
    f (:, :) double;
end

arguments (Output)
    u (:, :) double;
end

f = f(A.xy_perm_, :);
u = A.AFinv_HSS_.URV_Solve(f);
u = u(A.omega_inv_perm_, :);
u = reshape(u, A.nx_, A.ny_, []);
u = ifft2(u);
u = reshape(u, A.N_, []);

end