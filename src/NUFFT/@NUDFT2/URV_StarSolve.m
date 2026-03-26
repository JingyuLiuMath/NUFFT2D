function f = URV_StarSolve(A, u)

arguments (Input)
    A NUDFT2;
    u (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

u = fft(u) / A.N_;
f = A.AFinv_HSS_.URV_StarSolve(u);
f = f(A.x_inv_perm_, :);

end