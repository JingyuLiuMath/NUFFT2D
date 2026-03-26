function ConstructRootGenerators3(A, H)
% ConstructGenerators

% Jingyu Liu, November 20, 2024.

arguments (Input)
    A NUDFT2_2D_HSS;
    H (:, :) double;
end

M = A.row_global_size_;
nx = A.nx_;
ny = A.ny_;
N = A.col_global_size_;

if A.leaf_ == 1
    A.row_ind_ = A.row_offset_ + (1 : A.row_size_)';
    A.col_ind_ = A.col_offset_ + (1 : A.col_size_)';
    
    % Assign full mat.
    A.Amat_ = H(A.row_ind_, A.col_ind_);

    % Clear.
    A.row_ind_ = [];
    A.col_ind_ = [];
else
    % Construct B.
    for i = 1 : A.num_children_
        for j = [1 : (i - 1), (i + 1) : A.num_children_]
            A.Bmat_{i, j} = H(A.children_{i}.row_ind_, A.children_{j}.col_ind_);
        end
    end

    % Clear data.
    for i = 1 : A.num_children_
        A.children_{i}.row_ind_ = [];
        A.children_{i}.col_ind_ = [];
    end
end

end