function f = Apply(A, c)
% Apply

% Jingyu Liu, December 8, 2024.

arguments (Input)
    A NUDFT2;
    c (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

c = fft(c);
f = A.AFinv_HSS_.Apply(c);
f = f(A.x_inv_perm_, :);

end