function Construct_FaceSplitting(A, Ax, Ay)

arguments (Input)
    A NUDFT2_2D_HSS;
    Ax NUDFT2_HSS;
    Ay NUDFT2_HSS;
end

Fx = FaceSplitting_Transfers(Ax);
Fy = FaceSplitting_Transfers(Ay);
ConstructGenerators_FaceSplitting(A, Ax, Ay, Fx, Fy);

end
