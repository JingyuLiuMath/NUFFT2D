function URV_StarSolve_Downward(A, level)
% URV_StarSolve_Downward

% Jingyu Liu, December 9, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
end

if A.level_ == level
    % Update f_re.
    if ~isempty(A.URV_Omega_)
        A.fvec_ = A.URV_P_ * [A.URV_f_re_; A.URV_f_sk_];
        % Full QR.
        % A.fvec_ = A.URV_Omega_ * [
        %     A.fvec_;
        %     zeros(A.URV_m_old_ - A.URV_m_new_, A.vec_col_size_)];
        % Econ QR.
        A.fvec_ = A.URV_Omega_ * A.fvec_;
    else
        A.fvec_ = A.URV_P_ * [A.URV_f_re_; A.URV_f_sk_];
    end

    % Clear.
    A.URV_f_re_ = [];
    A.URV_f_sk_ = [];

    if A.leaf_ == 0
        A.URV_StarSolve_SplitVector();
    end
elseif A.leaf_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.URV_StarSolve_Downward(level);
    end
end

end