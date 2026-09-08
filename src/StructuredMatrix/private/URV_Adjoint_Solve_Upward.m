function data = URV_Adjoint_Solve_Upward(Fac, u, is_root)
% Apply the adjoint of back-substitution, including sibling interactions.
% fhvec accumulates adjoint contributions to incoming row-basis coefficients.

data.children = cell(1, Fac.num_children);
data.fhvec = zeros(Fac.row_rank, size(u, 2));
if Fac.leaf == 0
    col_size = 0;
    col_offset = 0;
    for i = 1 : Fac.num_children
        current_size = Fac.children{i}.col_size;
        data.children{i} = URV_Adjoint_Solve_Upward(Fac.children{i}, ...
            u((col_offset + 1) : (col_offset + current_size), :), false);
        if ~is_root
            data.fhvec = data.fhvec ...
                + Fac.children{i}.Rmat' * data.children{i}.fhvec;
        end
        col_offset = col_offset + current_size;
        col_size = col_size + Fac.children{i}.col_sk_size;
    end

    u = zeros(col_size, size(u, 2));
    col_offset = 0;
    for i = 1 : Fac.num_children
        uhvec = zeros(Fac.children{i}.col_rank, size(u, 2));
        for j = 1 : Fac.num_children
            if i ~= j
                uhvec = uhvec + Fac.Bmat{j, i}' * data.children{j}.fhvec;
            end
        end
        current_size = Fac.children{i}.col_sk_size;
        u((col_offset + 1) : (col_offset + current_size), :) ...
            = data.children{i}.u_sk + Fac.children{i}.V_sk * uhvec;
        col_offset = col_offset + current_size;
    end
end

if is_root
    data.u_sk = u;
else
    u = Fac.Q' * u;
    data.f_re = Fac.A_re_re' \ u(1 : Fac.re_size, :);
    data.u_sk = u((Fac.re_size + 1) : (Fac.re_size + Fac.col_sk_size), :) ...
        - Fac.A_re_sk' * data.f_re;
    data.fhvec = data.fhvec - Fac.U_re' * data.f_re;
end

end
