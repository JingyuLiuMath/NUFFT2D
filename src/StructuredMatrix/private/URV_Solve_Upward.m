function [data, f] = URV_Solve_Upward(Fac, f, is_root)
% Transform the right-hand side and pass skeleton rows to the parent.

data.children = cell(1, Fac.num_children);
if Fac.leaf == 0
    row_size = 0;
    for i = 1 : Fac.num_children
        row_size = row_size + Fac.children{i}.row_sk_size;
    end
    f_sk = zeros(row_size, size(f, 2));
    row_offset = 0;
    sk_offset = 0;
    for i = 1 : Fac.num_children
        current_size = Fac.children{i}.row_size;
        current_sk_size = Fac.children{i}.row_sk_size;
        [data.children{i}, f_sk((sk_offset + 1) : (sk_offset + current_sk_size), :)] ...
            = URV_Solve_Upward(Fac.children{i}, ...
            f((row_offset + 1) : (row_offset + current_size), :), false);
        row_offset = row_offset + current_size;
        sk_offset = sk_offset + current_sk_size;
    end
    f = f_sk;
end

if ~is_root
    if size(Fac.Omega, 1) > 0
        f = Fac.Omega' * f;
    end
    f = Fac.P' * f;
    data.f_re = f(1 : Fac.re_size, :);
    f = f((Fac.re_size + 1) : (Fac.re_size + Fac.row_sk_size), :);
end

end
