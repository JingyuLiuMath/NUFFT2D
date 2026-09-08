function c = Apply_Adjoint(A, f)

arguments (Input)
    A NUDFT2;
    f (:, :) double;
end

arguments (Output)
    c (:, :) double;
end

f = f(A.x_perm_, :);
c = A.AFinv_HSS_.Apply_Adjoint(f);
c = ifft(c) * A.N_;

end
