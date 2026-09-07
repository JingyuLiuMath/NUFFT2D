function Compress(A, tol)
% HSS recompression: bottom-up QR followed by top-down SVD.
% J. Xia, SIMAX 33 (2012), Section 5, Algorithm 3.
% https://doi.org/10.1137/110827788
% A is an unfactored HSS root, and 0 <= tol < 1.
% tol is the relative singular-value threshold for each local compression.

arguments (Input)
    A HSS;
    tol (1, 1) double;
end

% Upward.
for level = A.max_level_ : -1 : 1
    A.Compress_Upward(level);
end

% Root.
A.Compress_Root(tol);

% Downward.
for level = 1 : 1 : A.max_level_
    A.Compress_Downward(level, tol);
end

end
