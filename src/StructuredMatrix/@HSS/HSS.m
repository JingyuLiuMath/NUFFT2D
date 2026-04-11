classdef (Abstract) HSS < handle
    % HSS

    properties
        % *****************************************************************
        % PROPERTY: Tree.
        level_ (1, 1) double = 0;
        leaf_ (1, 1) double = 0;
        max_level_ (1, 1) double = 0;
        children_ (1, :) cell;
        num_children_ (1, 1) double = 0;
        % -----------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: Matrix.
        row_global_size_ (1, 1) double = 0;
        col_global_size_ (1, 1) double = 0;
        row_size_ (1, 1) double = 0;
        col_size_ (1, 1) double = 0;
        row_rank_ (1, 1) double = 0;  % size(U, 2), rank of the col space.
        col_rank_ (1, 1) double = 0;  % size(V, 2), rank of the row space.
        row_offset_ (1, 1) double = 0;
        col_offset_ (1, 1) double = 0;
        % -----------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: Generators.
        % Basis matrices.
        Umat_ (:, :) double;
        Vmat_ (:, :) double;

        % Transfer matrices.
        Rmat_ (1, :) cell;
        Wmat_ (1, :) cell;

        % Interaction matrices.
        Bmat_ (:, :) cell;

        % Full matrices.
        Amat_ (:, :) double;
        %------------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: Vectors.
        uvec_ (:, :) double;
        uhvec_ (:, :) double;
        fhvec_ (:, :) double;
        fvec_ (:, :) double;
        vec_col_size_ (:, :) double;
        %------------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: URV.
        % Matrices.
        URV_m_old_ (1, 1) double;
        URV_m_new_ (1, 1) double;
        URV_re_size_ (1, 1) double;
        URV_row_sk_size_ (1, 1) double;
        URV_col_sk_size_ (1, 1) double;
        URV_Omega_ (:, :) double;
        URV_V_sk_ (:, :);
        URV_Q_ (:, :) double;  % Col zeroing-out matrix.
        URV_P_ (:, :) double;  % Row zeroing-out matrix.
        URV_A_re_re_ (:, :) double;
        URV_A_re_sk_ (:, :) double;
        URV_A_sk_sk_ (:, :) double;
        URV_U_re_ (:, :) double;
        URV_U_sk_ (:, :) double;
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Vectors.
        URV_f_sk_ (:, :) double;
        URV_f_re_ (:, :) double;
        URV_f_reduced_row_ (:, :) double;
        URV_u_sk_ (:, :) double;
        URV_u_re_ (:, :) double;
        %------------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: Black-Box Construction.
        BBC_Omega_ (:, :) double;
        BBC_Y_ (:, :) double;
        BBC_Psi_ (:, :) double;
        BBC_Z_ (:, :) double;
        BBC_A_ (:, :) double;
        %------------------------------------------------------------------

    end

    methods
        % ****************************************************************
        % METHOD: Utilization
        U = Matrix_U(A);
        V = Matrix_V(A);
        mem = Storage(A);
        r = Rank(A);
        FillVector_Col(A, u);
        f = FetchVector_Row(A);
        FillVector_Row(A, f);
        u = FetchVector_Col(A);
        [max_row_leaf_size, max_col_leaf_size] = MaxLeafSize(A);
        %-----------------------------------------------------------------

        % ****************************************************************
        % METHOD: Apply.
        f = Apply(A, u)
        Apply_Upward(A, level);
        Apply_Root(A);
        Apply_Downward(A, level);

        u = StarApply(A, f);
        StarApply_Upward(A, level);
        StarApply_Root(A);
        StarApply_Downward(A, level);
        %-----------------------------------------------------------------

        % ****************************************************************
        % METHOD: URV factorization and the least squares solver.
        % Factorization.
        URV_Factor(A);
        URV_Eliminate(A, level);
        URV_Merge(A);
        URV_RootFactor(A);
        % Solution.
        u = URV_Solve(A, f);
        URV_Solve_Upward(A, level);
        URV_Solve_Merge(A);
        URV_Solve_Root(A);
        URV_Solve_Split(A);
        URV_Solve_Downward(A, level);
        %-----------------------------------------------------------------

        % ****************************************************************
        % METHOD: Block-Box Construction.
        BlackBoxConstruct(A, op_A, op_Astar, ...
            row_leaf_size, col_leaf_size, target_rank);
        BBC_FillAuxiliaryMatrix(A, Omega, Y, Psi, Z);
        BBC_ConstructGenerators(A, level, target_rank, tol);
        BBC_MergeAuxiliaryMatrix(A);
        BBC_ConstructRootGenerators(A);
        BBC_EliminateRootBMat(A);
        BBC_EliminateBMat(A, level);
        %-----------------------------------------------------------------

    end

end