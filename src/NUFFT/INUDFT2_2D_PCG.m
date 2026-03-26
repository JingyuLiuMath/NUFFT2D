function [u, flag, relres, iter, resvec] = INUDFT2_2D_PCG(A, xy, nx, ny, f, tol, maxit)

arguments (Input)
    A NUDFT2_2D;
    xy (:, 2) double;
    nx (1, 1) double;
    ny (1, 1) double;
    f (:, :) double;
    tol (1, 1) double = 1e-6;
    maxit (1, 1) double = min(nx, 20);
end

arguments (Output)
    u (:, :) double;
    flag (1, :) double;
    relres (1, :) double;
    iter (1, :) double;
    resvec (1, :) cell;
end

N = nx * ny;
num_rhs = size(f, 2);

% pcg.
omega = -[xy(:, 1) * nx, xy(:, 2) * ny];
f = MY_NUFFT1_2D(f, omega, nx, ny);
op_A = @(v) MY_NUFFT1_2D(MY_NUFFT2_2D(v, xy, nx, ny), omega, nx, ny);
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