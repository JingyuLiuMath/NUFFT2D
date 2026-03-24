function BBC_ConstructGenerators(A, level, target_rank, tol)
% BBC_ConstructGenerators

% Jingyu Liu, December 5, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
    target_rank (1, 1) double;
    tol (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 0
        A.BBC_MergeAuxiliaryMatrix();
    end
    row_target_rank = min(target_rank, size(A.BBC_Y_, 1));
    col_target_rank = min(target_rank, size(A.BBC_Z_, 1));
    P = NullBasis(A.BBC_Omega_, row_target_rank);
    [A.Umat_, A.row_rank_] = ColBasis(A.BBC_Y_ * P, row_target_rank,  tol);
    Q = NullBasis(A.BBC_Psi_, col_target_rank);
    [A.Vmat_, A.col_rank_] = ColBasis(A.BBC_Z_ * Q, col_target_rank, tol);
    Y_OmegaInv = A.BBC_Y_ / A.BBC_Omega_;
    Z_PsiInv = A.BBC_Z_ / A.BBC_Psi_;
    A.BBC_A_ = Y_OmegaInv - A.Umat_ * (A.Umat_' * Y_OmegaInv) ...
        + A.Umat_ * (A.Umat_' ...
        * (Z_PsiInv - A.Vmat_ * (A.Vmat_' * Z_PsiInv))');
    if A.leaf_ == 1
        A.Amat_ = A.BBC_A_;
    else
        % Assign R and W.
        row_offset = 0;
        col_offset = 0;
        for i = 1 : A.num_children_
            current_row_size = A.children_{i}.row_rank_;
            current_col_size = A.children_{i}.col_rank_;
            A.Rmat_{i} = A.Umat_(...
                (row_offset + 1) : (row_offset + current_row_size), :);
            A.Wmat_{i} = A.Vmat_(...
                (col_offset + 1) : (col_offset + current_col_size), :);
            row_offset = row_offset + current_row_size;
            col_offset = col_offset + current_col_size;
        end

        % Assign B.
        row_offset = 0;
        for i = 1 : A.num_children_
            current_row_size = A.children_{i}.row_rank_;
            col_offset = 0;
            for j = 1 : A.num_children_
                current_col_size = A.children_{j}.col_rank_;
                A.Bmat_{i, j} = A.BBC_A_(...
                    (row_offset + 1) : (row_offset + current_row_size), ...
                    (col_offset + 1) : (col_offset + current_col_size));
                col_offset = col_offset + current_col_size;
            end
            row_offset = row_offset + current_row_size;
        end
    end
elseif A.leaf_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.BBC_ConstructGenerators(level, target_rank, tol);
    end
end

end