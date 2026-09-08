function Fac = URV_Factor(A)
% URV_Factor Compute a separate URV factor tree.
% Fac stores the tree, coupling generators and the factors needed to solve.
% At the root, P and A_re_re contain the QR factors of the reduced matrix.
% Other nodes also store Omega, Q, A_re_sk, U_re, V_sk and elimination sizes.
% The HSS matrix is preserved; no right-hand-side workspace is stored.

arguments (Input)
    A HSS;
end

arguments (Output)
    Fac (1, 1) struct;
end

[Fac, data] = URV_Eliminate(A, true);
[Fac.P, Fac.A_re_re] = qr(data.Amat);

end
