function URV_Solve_Upward(A, level)
% URV_Solve_Upward

% Jingyu Liu, November 23, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 0
        A.URV_Solve_MergeVector();
    end
    
    if ~isempty(A.URV_Omega_)
        A.fvec_ = A.URV_Omega_' * A.fvec_;
        A.fvec_(1 : A.URV_m_new_, :) ...
            = A.URV_P_' * A.fvec_(1 : A.URV_m_new_, :);
        A.URV_f_re_ = A.fvec_(1 : A.URV_re_size_, :);
        A.URV_f_sk_ = A.fvec_(...
            (A.URV_re_size_ + 1) : (A.URV_re_size_ + A.URV_row_sk_size_), :);
        % Full QR.
        % A.URV_f_reduced_row_ = A.fvec_(...
        %     (A.URV_re_size_ + A.URV_row_sk_size_ + 1) : A.URV_m_old_, :);
    else
        A.fvec_ = A.URV_P_' * A.fvec_;
        A.URV_f_re_ = A.fvec_(1 : A.URV_re_size_, :);
        A.URV_f_sk_ = A.fvec_(...
            (A.URV_re_size_ + 1) : (A.URV_re_size_ + A.URV_row_sk_size_), :);
    end

    % Clear.
    A.fvec_ = [];
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.URV_Solve_Upward(level);
    end
end

end