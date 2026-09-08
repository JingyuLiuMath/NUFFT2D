function u = URV_Solve_Downward(Fac, data, u_sk, fhvec, is_root)
% Back-substitute redundant coordinates and propagate sibling interactions.

if is_root
    u = u_sk;
else
    u = zeros(size(Fac.Q, 2), size(u_sk, 2));
    u(1 : Fac.re_size, :) = Fac.A_re_re \ ...
        (data.f_re - Fac.A_re_sk * u_sk - Fac.U_re * fhvec);
    u((Fac.re_size + 1) : (Fac.re_size + Fac.col_sk_size), :) = u_sk;
    u = Fac.Q * u;
end

if Fac.leaf == 0
    col_offset = 0;
    for i = 1 : Fac.num_children
        current_size = Fac.children{i}.col_sk_size;
        data.children{i}.u_sk = u((col_offset + 1) : (col_offset + current_size), :);
        data.children{i}.uhvec = Fac.children{i}.V_sk' * data.children{i}.u_sk;
        col_offset = col_offset + current_size;
    end

    u = zeros(Fac.col_size, size(u_sk, 2));
    col_offset = 0;
    for i = 1 : Fac.num_children
        if is_root
            fhvec_i = zeros(Fac.children{i}.row_rank, size(u_sk, 2));
        else
            fhvec_i = Fac.children{i}.Rmat * fhvec;
        end
        for j = 1 : Fac.num_children
            if i ~= j
                fhvec_i = fhvec_i + Fac.Bmat{i, j} * data.children{j}.uhvec;
            end
        end
        current_size = Fac.children{i}.col_size;
        u((col_offset + 1) : (col_offset + current_size), :) ...
            = URV_Solve_Downward(Fac.children{i}, data.children{i}, ...
            data.children{i}.u_sk, fhvec_i, false);
        col_offset = col_offset + current_size;
    end
end

end
