function URV_Solve_Downward(A, level)
% URV_Solve_Downward

% Jingyu Liu, November 23, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
end

if A.level_ == level
    % Update f_re.
    A.URV_f_re_ = A.URV_f_re_ ...
        - A.URV_A_re_sk_ * A.URV_u_sk_ ...
        - A.URV_U_re_ * A.fhvec_;
    u_re = A.URV_A_re_re_ \ A.URV_f_re_;
    A.uvec_ = A.URV_Q_ * [u_re; A.URV_u_sk_];

    % Clear.
    A.URV_f_re_ = [];
    A.URV_u_sk_ = [];

    if A.leaf_ == 0
        A.URV_Solve_SplitVector();
    end
elseif A.leaf_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.URV_Solve_Downward(level);
    end
end

end