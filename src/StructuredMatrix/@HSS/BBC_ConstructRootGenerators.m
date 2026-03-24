function BBC_ConstructRootGenerators(A)
% BBC_ConstructRootGenerators

% Jingyu Liu, December 5, 2024.

arguments (Input)
    A HSS;
end

if A.leaf_ == 1
    A.Amat_ = A.BBC_Y_ / A.BBC_Omega_;
else
    A.BBC_MergeAuxiliaryMatrix();
    A.BBC_A_ = A.BBC_Y_ / A.BBC_Omega_;

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

% Clear.
A.BBC_Omega_ = [];
A.BBC_Y_ = [];
A.BBC_Psi_ = [];
A.BBC_Z_ = [];
A.BBC_A_ = [];

end