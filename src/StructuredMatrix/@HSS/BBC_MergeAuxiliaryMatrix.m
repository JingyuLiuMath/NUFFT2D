function BBC_MergeAuxiliaryMatrix(A)
% BBC_MergeAuxiliaryMatrix

% Jingyu Liu, December 5, 2024.

arguments (Input)
    A HSS;
end

A.BBC_Omega_ = [];
A.BBC_Y_ = [];
A.BBC_Psi_ = [];
A.BBC_Z_ = [];
for i = 1 : A.num_children_
    A.BBC_Omega_ = [A.BBC_Omega_;...
        A.children_{i}.Vmat_' * A.children_{i}.BBC_Omega_];
    A.BBC_Y_ = [A.BBC_Y_; ...
        A.children_{i}.Umat_' ...
        * (A.children_{i}.BBC_Y_ ...
        - A.children_{i}.BBC_A_ * A.children_{i}.BBC_Omega_)];
    A.BBC_Psi_ = [A.BBC_Psi_;...
        A.children_{i}.Umat_' * A.children_{i}.BBC_Psi_];
    A.BBC_Z_ = [A.BBC_Z_; ...
        A.children_{i}.Vmat_' ...
        * (A.children_{i}.BBC_Z_ ...
        - A.children_{i}.BBC_A_' * A.children_{i}.BBC_Psi_)];

    % Clear.
    A.children_{i}.BBC_Omega_ = [];
    A.children_{i}.BBC_Y_ = [];
    A.children_{i}.BBC_Psi_ = [];
    A.children_{i}.BBC_Z_ = [];
    A.children_{i}.BBC_A_ = [];
    if A.children_{i}.leaf_ == 0
        A.children_{i}.Umat_ = [];
        A.children_{i}.Vmat_ = [];
    end

end