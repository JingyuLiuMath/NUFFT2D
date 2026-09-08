function f = URV_Adjoint_Solve(Fac, u)
% URV_Adjoint_Solve Apply the adjoint of the URV solution operator.
% For full column rank A, return the minimum-norm solution of A' * f = u.
% For a wide reduced root, apply the adjoint of URV_Solve's basic solution.

arguments (Input)
    Fac (1, 1) struct;
    u (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

% Upward.
data = URV_Adjoint_Solve_Upward(Fac, u, true);

% Root.
m = size(Fac.A_re_re, 1);
n = size(Fac.A_re_re, 2);
re_size = min(m, n);
f = zeros(m, size(u, 2));
f(1 : re_size, :) = Fac.A_re_re(1 : re_size, 1 : re_size)' ...
    \ data.u_sk(1 : re_size, :);
f = Fac.P * f;

% Downward.
f = URV_Adjoint_Solve_Downward(Fac, data, f, true);

end
