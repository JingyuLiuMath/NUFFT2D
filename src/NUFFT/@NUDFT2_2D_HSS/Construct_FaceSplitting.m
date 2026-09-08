function [p, q] = Construct_FaceSplitting(A, Ax, Ay, px, py)
% Construct_FaceSplitting Assemble raw 2D generators from two 1D HSS roots.
% px and py map source rows to the same original row ordering.
% The input is a fresh root. Recompression is performed separately.

arguments (Input)
    A NUDFT2_2D_HSS;
    Ax NUDFT2_HSS;
    Ay NUDFT2_HSS;
    px (:, 1) double;
    py (:, 1) double;
end

arguments (Output)
    p (:, 1) double;
    q (:, 1) double;
end

[p, q, index] = BuildTree_FaceSplitting(A, Ax, Ay, px, py);
index.p = p;
[~, Fx] = FaceSplitting_Transfers(Ax);
[~, Fy] = FaceSplitting_Transfers(Ay);
ConstructGenerators_FaceSplitting(A, Ax, Ay, Fx, Fy, index, true);

end
