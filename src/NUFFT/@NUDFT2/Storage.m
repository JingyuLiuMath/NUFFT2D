function mem = Storage(A)

arguments (Input)
    A NUDFT2;
end

arguments (Output)
    mem (1, 1) double;
end

mem = A.AFinv_HSS_.Storage();

end