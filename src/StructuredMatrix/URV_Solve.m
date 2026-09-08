function u = URV_Solve(Fac, f)
% URV_Solve Solve A * u = f using the factor tree from A.URV_Factor().
% For full column rank A, return the least-squares solution.
% A wide reduced root uses a basic solution with zero free coordinates.

arguments (Input)
    Fac (1, 1) struct;
    f (:, :) double;
end

arguments (Output)
    u (:, :) double;
end

% Upward.
[data, f] = URV_Solve_Upward(Fac, f, true);

% Root.
m = size(Fac.A_re_re, 1);
n = size(Fac.A_re_re, 2);
re_size = min(m, n);
f = Fac.P' * f;
u = zeros(n, size(f, 2));
u(1 : re_size, :) = Fac.A_re_re(1 : re_size, 1 : re_size) ...
    \ f(1 : re_size, :);

% Downward.
u = URV_Solve_Downward(Fac, data, u, zeros(0, size(f, 2)), true);

end
