function [Q, B] = RangeFinderImplict(m, n, ...
    AX_fun, XA_fun, rank_or_tol, power_iter)
% RangeFinderImplict
% rank_or_tol >= 1:
% Find orthogonal Q (m-by-l) and B (l-by-n) of a m-by-n A such that Q
% approximates the leading l left singular vectors.
% l = k + p where k = rank_or_tol and p = oversampling_number.
% rank_or_tol < 1:
% Find orthonormal Q (m-by-l) and B (l-by-n) of a m-by-n matrix A such that
% norm(A - Q * B, "fro") / norm(A, "fro") < rank_or_tol.

% Only X -> AX and X -> XA is available.

arguments (Input)
    m (1, 1) double;
    n (1, 1) double;
    AX_fun function_handle;
    XA_fun function_handle;
    rank_or_tol (1, 1) double;
    power_iter (1, 1) double = 0;
end

arguments (Output)
    Q (:, :) double;
    B (:, :) double;
end

min_szA = min(m, n);

if min_szA == 0
    Q = zeros(m, 0);
    B = zeros(0, n);
    return;
end

if rank_or_tol >= 1
    oversampling_number = min(min_szA, 5);
    k = min(rank_or_tol, min_szA);
    Omega = randn(n, k + oversampling_number);
    [Q, ~] = qr(AX_fun(Omega), "econ");
    for it = 1 : power_iter
        [Q, ~] = qr(XA_fun(Q')', "econ");
        [Q, ~] = qr(AX_fun(Q), "econ");
    end
    % Q = Q(:, 1 : k);
    B = XA_fun(Q');
else
    Q = zeros(m, 0);
    B = zeros(0, n);

    % Keep the test samples independent of the construction samples.
    Omega_test = randn(n, 20);
    Y_test = AX_fun(Omega_test);
    norm_test = norm(Y_test, "fro");
    if norm_test == 0
        return;
    end

    block_size0 = min(max(floor(0.1 * min_szA), 5), min_szA);
    approx_err = 1;

    while size(Q, 2) < min_szA && approx_err > rank_or_tol
        block_size = min(block_size0, min_szA - size(Q, 2));
        Omega_i = randn(n, block_size);

        Y_i = AX_fun(Omega_i) - Q * (B * Omega_i);
        [Q_i, ~, ~, ~] = MyQRSketch(Y_i, block_size);

        for it = 1 : power_iter
            Z_i = XA_fun(Q_i')' - B' * (Q' * Q_i);
            [Z_i, ~, ~, ~] = MyQRSketch(Z_i, block_size);

            Y_i = AX_fun(Z_i) - Q * (B * Z_i);
            [Q_i, ~, ~, ~] = MyQRSketch(Y_i, block_size);
        end

        Y_i = Q_i - Q * (Q' * Q_i);
        Y_i = Y_i - Q * (Q' * Y_i);
        [Q_i, ~, ~, ~] = MyQRSketch(Y_i, block_size);

        if isempty(Q_i)
            break;
        end

        B_i = XA_fun(Q_i');
        Q = [Q, Q_i];
        B = [B; B_i];

        approx_err = norm(Y_test - Q * (B * Omega_test), "fro") ...
            / norm_test;
    end

    if approx_err > rank_or_tol
        error("RangeFinderImplict:NotConverged", ...
            "The sampled relative residual %g exceeds tolerance %g.", ...
            approx_err, rank_or_tol);
    end
end

end