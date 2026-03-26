function URV_StarSolve_Upward(A, level)
% URV_StarSolve_Upward

% Jingyu Liu, December 9, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
end

if A.level_ == level  
    if A.leaf_ == 0
        A.URV_StarSolve_MergeVector();
    end
    
    A.uvec_ = A.URV_Q_' * A.uvec_;
    A.URV_u_re_ = A.uvec_(1 : A.URV_re_size_, :);
    if isempty(A.URV_u_sk_)
        A.URV_u_sk_ = A.uvec_(...
            (A.URV_re_size_ + 1) : (A.URV_re_size_ + A.URV_col_sk_size_), :);
    else
        A.URV_u_sk_ = A.URV_u_sk_ + A.uvec_(...
            (A.URV_re_size_ + 1) : (A.URV_re_size_ + A.URV_col_sk_size_), :);
    end
    A.URV_f_re_ = A.URV_A_re_re_' \ A.URV_u_re_;
    
    % Clear.
    A.uvec_ = [];
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.URV_StarSolve_Upward(level);
    end
end

end