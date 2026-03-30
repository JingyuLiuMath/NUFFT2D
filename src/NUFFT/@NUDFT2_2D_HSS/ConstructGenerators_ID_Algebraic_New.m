function ConstructGenerators_ID_Algebraic_New(A, ...
    H, ...
    level, rank_or_tol)

arguments (Input)
    A NUDFT2_2D_HSS;
    H (:, :) double;
    level (1, 1) double;
    rank_or_tol (1, 1) double;
end

M = A.row_global_size_;
nx = A.nx_;
ny = A.ny_;
N = A.col_global_size_;

if A.level_ == level
    if A.leaf_ == 1
        A.row_ind_ = A.row_offset_ + (1 : A.row_size_)';
        A.col_ind_ = A.col_offset_ + (1 : A.col_size_)';
    elseif A.leaf_ == 0
        A.row_ind_ = A.row_sk_offset_ + (1 : A.row_sk_size_)';
        A.col_ind_ = A.col_sk_offset_ + (1 : A.col_sk_size_)';
        for i = 1 : A.num_children_
            c_I = A.children_{i}.row_sk_offset_ + (1 : A.children_{i}.row_sk_size_)';
            for j = [1 : (i - 1), (i + 1) : A.num_children_]
                c_J = A.children_{j}.col_sk_offset_ + (1 : A.children_{j}.col_sk_size_)';
                A.Bmat_{i, j} = H(c_I, c_J);
            end
        end
    end

    % Construct U.
    Jc = (1 : size(H, 2))';
    Jc(A.col_ind_) = [];
    A_I_Jc = H(A.row_ind_, Jc);
    [row_sk, U, A.row_rank_, row_re] = LowRank_Row_ID(A_I_Jc, rank_or_tol);

    % Construct V using proxy surface.
    Ic = (1 : size(H, 1))';
    Ic(A.row_ind_) = [];
    A_Ic_J = H(Ic, A.col_ind_);
    [col_sk, V, A.col_rank_, col_re] = LowRank_ID(A_Ic_J, rank_or_tol);

    if A.leaf_ == 1
        % Assign U and V.
        A.Umat_ = U;
        A.Vmat_ = V;

        % Assign full mat.
        A.Amat_ = H(A.row_ind_, A.col_ind_);
    else
        % Assign R and W.
        offset = 0;
        for i = 1 : A.num_children_
            current_size = size(A.children_{i}.row_ind_, 1);
            A.Rmat_{i} = U((offset + 1) : (offset + current_size), :);
            offset = offset + current_size;
        end

        offset = 0;
        for i = 1 : A.num_children_
            current_size = size(A.children_{i}.col_ind_, 1);
            A.Wmat_{i} = V((offset + 1) : (offset + current_size), :);
            offset = offset + current_size;
        end

        % Clear data.
        for i = 1 : A.num_children_
            A.children_{i}.row_ind_ = [];
            A.children_{i}.col_ind_ = [];
        end
    end

    % Update row and col.
    A.row_re_ = A.row_ind_(row_re);
    A.col_re_ = A.col_ind_(col_re);
    A.row_ind_ = A.row_ind_(row_sk);
    A.col_ind_ = A.col_ind_(col_sk);
else
    for i = 1 : A.num_children_
        A.children_{i}.ConstructGenerators_ID_Algebraic_New(...
            H, ...
            level, rank_or_tol);
    end
end

end