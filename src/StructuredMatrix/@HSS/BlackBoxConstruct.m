function BlackBoxConstruct(A, op_A, op_Astar, ...
    row_leaf_size, col_leaf_size, target_rank, tol)
% BlackBoxConstruct

% Jingyu Liu, December 4, 2024.

arguments (Input)
    A HSS;
    op_A function_handle;
    op_Astar function_handle;
    row_leaf_size (1, 1) double;
    col_leaf_size (1, 1) double;
    target_rank (1, 1) double;
    tol (1, 1) double;
end

% Settings.
p = 5;
r = target_rank + p;
s_row = max(r + row_leaf_size, 3 * r);
s_col = max(r + col_leaf_size, 3 * r);

% Sampling.
Omega = randn(A.col_size_, s_col);
Psi = randn(A.row_size_, s_row);

Y = op_A(Omega);
Z = op_Astar(Psi);

A.BBC_FillAuxiliaryMatrix(Omega, Y, Psi, Z);

% Recursive construction.
for level = A.max_level_ : -1 : 1
    A.BBC_ConstructGenerators(level, r, tol);
end

A.BBC_ConstructRootGenerators();

% Eliminate B{i, i} matrices.
A.BBC_EliminateRootBMat();

for level = 1 : 1 : A.max_level_
    A.BBC_EliminateBMat(level);
end

end