function r = Rank(A)

arguments (Input)
    A NUDFT2;
end

arguments (Output)
    r(1, 1) double;
end

r = A.AFinv_HSS_.Rank();

end