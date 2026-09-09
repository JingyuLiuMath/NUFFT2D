function f = Apply(A, c)

arguments (Input)
    A NUDFT2;
    c (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

c = fft(c, [], 1);
f = A.AFinv_HSS_.Apply(c);
f = f(A.x_inv_perm_, :);

end