function ConstructRootGenerators4(A, H)

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
        c_I = A.children_{i}.row_sk_offset_ + (1 : A.children_{i}.row_sk_size_)';
        for j = [1 : (i - 1), (i + 1) : A.num_children_]
            c_J = A.children_{j}.col_sk_offset_ + (1 : A.children_{j}.col_sk_size_)';
            A.Bmat_{i, j} = H(c_I, c_J);
        end
    end

    % Clear data.
    for i = 1 : A.num_children_
        A.children_{i}.row_ind_ = [];
        A.children_{i}.col_ind_ = [];
    end
end

end