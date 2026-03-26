function K = NUFFT2_Kernel_Real(x, y, N)

arguments (Input)
    x (:, 1) double;
    y (:, 1) double;
    N (1, 1) double;
end

arguments (Output)
    K (:, :) double;
end

phi = -2 * pi * 1i * (x - y.');
Num = expm1(-2 * pi * 1i * x * N) / N;
Dom = expm1(phi);
K = Num ./ Dom;

if max(isinf(K(:))) == 1
    [I, J] = find(abs(K) == inf);
    % subK = K(I, J);
    % abs_subK = abs(subK);
    % subD = Dom(I, J);
    % keyboard;
    for rit = 1 : length(I)
        x_i = x(I(rit));
        for cit = 1 : length(J)
            y_j = y(J(cit));
            num = expm1(-2 * pi * 1i * x_i * N) / N;
            dom = expm1(-2 * pi * 1i * (x_i - y_j));
            if abs(num) <= 1e-15 && abs(dom) <= 1e-15
                v = 1;
            else
                v = 0;
            end
            K(I(rit), J(cit)) = v;
        end
    end
end

end