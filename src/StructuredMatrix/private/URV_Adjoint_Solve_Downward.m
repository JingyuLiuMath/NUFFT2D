function f = URV_Adjoint_Solve_Downward(Fac, data, f_sk, is_root)
% Restore row coordinates by applying the adjoint right-hand-side transforms.

if is_root
    f = f_sk;
else
    f = zeros(size(Fac.P, 2), size(f_sk, 2));
    f(1 : Fac.re_size, :) = data.f_re;
    f((Fac.re_size + 1) : (Fac.re_size + Fac.row_sk_size), :) = f_sk;
    f = Fac.P * f;
    if size(Fac.Omega, 1) > 0
        f = Fac.Omega * f;
    end
end

if Fac.leaf == 0
    f_sk = f;
    f = zeros(Fac.row_size, size(f_sk, 2));
    row_offset = 0;
    sk_offset = 0;
    for i = 1 : Fac.num_children
        current_size = Fac.children{i}.row_size;
        current_sk_size = Fac.children{i}.row_sk_size;
        f((row_offset + 1) : (row_offset + current_size), :) ...
            = URV_Adjoint_Solve_Downward(Fac.children{i}, data.children{i}, ...
            f_sk((sk_offset + 1) : (sk_offset + current_sk_size), :), false);
        row_offset = row_offset + current_size;
        sk_offset = sk_offset + current_sk_size;
    end
end

end
