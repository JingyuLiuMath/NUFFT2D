function Recompress(A, tol)
% Recompress

arguments (Input)
    A HSS;
    tol (1, 1) double;
end

if A.leaf_ == 1
    return;
end

Recompress_Upward(A, true);
Recompress_Downward(A, zeros(0), zeros(0), tol, true);

end
