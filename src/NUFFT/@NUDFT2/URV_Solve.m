function u = URV_Solve(A, f)
% URV_Solve

% Jingyu Liu, December 8, 2024.

arguments (Input)
    A NUDFT2;
    f (:, :) double;
end

arguments (Output)
    u (:, :) double;
end

f = f(A.x_perm_, :);
c = A.AFinv_HSS_.URV_Solve(f);
u = ifft(c);

end