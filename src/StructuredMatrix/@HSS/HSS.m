classdef (Abstract) HSS < handle
    % HSS HSS matrix representation.
    % Algorithm workspaces and factorizations are stored in separate structs.

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
        row_rank_ (1, 1) double = 0;  % Column basis rank.
        col_rank_ (1, 1) double = 0;  % Row basis rank.
        row_offset_ (1, 1) double = 0;
        col_offset_ (1, 1) double = 0;
        % -----------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: Generators.
        % Leaf basis matrices.
        Umat_ (:, :) double;
        Vmat_ (:, :) double;

        % Transfer matrices.
        Rmat_ (:, :) double;
        Wmat_ (:, :) double;

        % Sibling interaction matrices.
        Bmat_ (:, :) cell;

        % Dense diagonal block at a leaf.
        Amat_ (:, :) double;
        % -----------------------------------------------------------------
    end

    methods
        % *****************************************************************
        % METHOD: Utilization.
        r = Rank(A);
        mem = Storage(A);
        U = Matrix_U(A);
        V = Matrix_V(A);
        % -----------------------------------------------------------------

        % *****************************************************************
        % METHOD: Apply.
        f = Apply(A, u);
        u = Apply_Adjoint(A, f);
        % -----------------------------------------------------------------

        % *****************************************************************
        % METHOD: Recompression.
        Recompress(A, tol);
        % -----------------------------------------------------------------

        % *****************************************************************
        % METHOD: URV factorization.
        Fac = URV_Factor(A);
        % -----------------------------------------------------------------
    end

end
