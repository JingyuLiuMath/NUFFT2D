function [u, flag, relres, iter, resvec] = INUDFT2_PCG(A, x, N, f, tol, maxit)

arguments (Input)
    A;
    x (:, 1) double;
    N (1, 1) double;
    f (:, :) double;
    tol (1, 1) double = 1e-6;
    maxit (1, 1) double = min(N, 20);
end

arguments (Output)
    u (:, :) double;
    flag (1, :) double;
    relres (1, :) double;
    iter (1, :) double;
    resvec (1, :) cell;
end

num_rhs = size(f, 2);

% pcg.
f = MY_NUFFT1(f, omega, N);
op_A = @(v) MY_NUFFT1(MY_NUFFT2(v, x, N), omega, N);
op_M = @(v) A.URV_Solve(A.URV_StarSolve(v));

u = zeros(N, num_rhs);
flag = zeros(1, num_rhs);
relres = zeros(1, num_rhs);
iter = zeros(1, num_rhs);
resvec = cell(1, num_rhs);
for it = 1 : num_rhs
    [u(:, it), flag(it), relres(it), iter(it), resvec{it}] = pcg(op_A, f(:, it), tol, maxit, op_M);
end

end