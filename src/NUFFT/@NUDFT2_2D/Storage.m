function mem = Storage(A)
% Storage

% Jingyu Liu, November 19, 2024.

arguments (Input)
    A NUDFT2_2D;
end

arguments (Output)
    mem (1, 1) double;
end

mem = A.AFinv_HSS_.Storage();

end