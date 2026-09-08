function f = Adjoint_Solve(A, u)

arguments (Input)
    A NUDFT2_2D;
    u (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

u = reshape(u, A.nx_, A.ny_, []);
u = fft2(u) / A.N_;
u = reshape(u, A.N_, []);
u = u(A.omega_perm_, :);
f = URV_Adjoint_Solve(A.Factor_, u);
f = f(A.xy_inv_perm_, :);

end
