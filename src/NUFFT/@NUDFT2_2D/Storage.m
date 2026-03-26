function mem = Storage(A)
% Storage

arguments (Input)
    A NUDFT2_2D;
end

arguments (Output)
    mem (1, 1) double;
end

mem = A.AFinv_HSS_.Storage();

end