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

op_A = @(v, opt) afun(v, opt, xy, nx, ny);
op_M = @(v, opt) mfun(v, opt, A);

u = zeros(N, num_rhs);
flag = zeros(1, num_rhs);
relres = zeros(1, num_rhs);
iter = zeros(1, num_rhs);
resvec = cell(1, num_rhs);
for it = 1 : num_rhs
    [u(:, it), flag(it), relres(it), iter(it), resvec{it}] = pcg(op_A, f(:, it), tol, maxit, op_M);
end

end

function y = afun(v, opt, x, N)

if strcmp(opt,'notransp')
    y = MY_NUFFT2(v, x, N);
else
    omega = -N * x;
    y = MY_NUFFT1(v, omega, N);
end

end

function y = mfun(v, opt, A)

if strcmp(opt,'notransp')
    y = A.URV_Solve(v);
else
    y = A.URV_StarSolve(v);
end

end