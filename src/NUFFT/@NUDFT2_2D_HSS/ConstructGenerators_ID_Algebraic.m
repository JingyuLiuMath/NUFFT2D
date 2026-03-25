function ConstructGenerators_ID_Algebraic(A, ...
    H, ...
    row_active, col_active, ...
    level, rank_or_tol)
% ConstructGenerators

% Jingyu Liu, November 24, 2024.

arguments (Input)
    A NUDFT2_2D_HSS;
    H (:, :) double;
    row_active (:, 1) double;
    col_active (:, 1) double;
    level (1, 1) double;
    rank_or_tol (1, 1) double;
end

M = A.row_global_size_;
nx = A.nx_;
ny = A.ny_;
N = A.col_global_size_;

if A.level_ == level
    if A.leaf_ == 0
        % Merge data from children and construct B.
        A.row_ind_ = [];
        A.col_ind_ = [];
        for i = 1 : A.num_children_
            A.row_ind_ = [A.row_ind_; A.children_{i}.row_ind_];
            A.col_ind_ = [A.col_ind_; A.children_{i}.col_ind_];
            for j = [1 : (i - 1), (i + 1) : A.num_children_]
                A.Bmat_{i, j} = H(A.children_{i}.row_ind_, A.children_{j}.col_ind_);
            end
        end
    end

    % Construct U.
    Jc = (1 : N)';
    Jc(A.self_col_ind_) = [];
    Jc = intersect(Jc, col_active);
    A_I_Jc = H(A.row_ind_, Jc);
    [row_sk, U, A.row_rank_] = LowRank_Row_ID(A_I_Jc, rank_or_tol);

    % Construct V.
    Ic = (1 : M)';
    Ic(A.self_row_ind_) = [];
    Ic = intersect(Ic, row_active);
    A_Ic_J = H(Ic, A.col_ind_);
    [col_sk, V, A.col_rank_] = LowRank_ID(A_Ic_J, rank_or_tol);

    % figure();
    % plot(xy_I(:, 1), xy_I(:, 2), "rx", "DisplayName", "row pts");
    % hold on;
    % plot(col_pos_J(:, 1) / nx, col_pos_J(:, 2) / ny, "go", "DisplayName", "col pts");
    % keyboard;

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
    A.row_ind_ = A.row_ind_(row_sk);
    A.col_ind_ = A.col_ind_(col_sk);
else
    for i = 1 : A.num_children_
        A.children_{i}.ConstructGenerators_ID_Algebraic(...
            H, ...
            row_active, col_active, ...
            level, rank_or_tol);
    end
end

end