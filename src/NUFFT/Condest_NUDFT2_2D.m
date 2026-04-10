function [kappa1_est, info] = Condest_NUDFT2_2D(A, xy, nx, ny, maxit_pcg, tol_pcg, maxit_condest, restarts)
% Condest_NUDFT2_2D Estimate condition number of 2D NUDFT operator.
%
% This routine follows operator definitions used in INUDFT2_2D_PCG.
% It estimates cond_1(A^H A) = ||A^H A||_1 * ||(A^H A)^{-1}||_1
% through Hager algorithm via Condest_Implict.
%
% Note:
%   The inverse action is approximated by PCG on (A^H A) x = b, with an
%   A-based preconditioner.

arguments (Input)
    A NUDFT2_2D;
    xy (:, 2) double;
    nx (1, 1) double;
    ny (1, 1) double;
    maxit_pcg (1, 1) double {mustBeInteger, mustBePositive} = 100;
    tol_pcg (1, 1) double {mustBePositive} = 1e-10;
    maxit_condest (1, 1) double {mustBeInteger, mustBePositive} = 8;
    restarts (1, 1) double {mustBeInteger, mustBePositive} = 2;
end

arguments (Output)
    kappa1_est (1, 1) double;
    info (1, 1) struct;
end

N = nx * ny;
M = size(xy, 1);

op_A = @(v, opt) nufft_normal_op(v, opt, xy, nx, ny);
op_A_inv = @(v, opt) nufft_normal_inv_op(v, opt, A, xy, nx, ny, maxit_pcg, tol_pcg);

[kappa1_est, info_cond] = Condest_Implict(op_A, op_A_inv, N, maxit_condest, restarts);

raw_kappa1_est = kappa1_est;
kappa1_est = sqrt(kappa1_est);

info = info_cond;
info.M = M;
info.N = N;
info.nx = nx;
info.ny = ny;
info.raw_kappa1_est = raw_kappa1_est;
info.target = "sqrt(cond_1(A^H*A))";

end

function y = nufft_normal_op(v, opt, xy, nx, ny)

y = MY_NUFFT2_2D(v, xy, nx, ny);
omega = -[xy(:, 1) * nx, xy(:, 2) * ny];
y = MY_NUFFT1_2D(y, omega, nx, ny);

end

function y = nufft_normal_inv_op(v, opt, A, xy, nx, ny, maxit_pcg, tol_pcg)

op_B = @(u) nufft_normal_op(u, 'notransp', xy, nx, ny);
op_M = @(u) A.URV_Solve(A.URV_StarSolve(u));

y = zeros(size(v));
for it = 1 : size(v, 2)
    y(:, it) = pcg(op_B, v(:, it), tol_pcg, maxit_pcg, op_M);
end

end