function Construct_FaceSplitting(A, Ax, Ay, row_index)
% Construct generators on the existing 2D tree from two 1D HSS roots.

arguments (Input)
    A NUDFT2_2D_HSS;
    Ax NUDFT2_HSS;
    Ay NUDFT2_HSS;
    row_index (1, 1) struct;
end

% D(parent)|child = D(child)*E + U(child)*F.
Fx = FaceSplitting_Transfers(Ax);
Fy = FaceSplitting_Transfers(Ay);

% Generators. Recompression is performed separately.
ConstructGenerators_FaceSplitting(A, Ax, Ay, Fx, Fy, row_index);

end
