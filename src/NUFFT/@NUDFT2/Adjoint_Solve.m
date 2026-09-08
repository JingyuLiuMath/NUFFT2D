function f = Adjoint_Solve(A, u)

arguments (Input)
    A NUDFT2;
    u (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

u = fft(u) / A.N_;
f = URV_Adjoint_Solve(A.Factor_, u);
f = f(A.x_inv_perm_, :);

end
