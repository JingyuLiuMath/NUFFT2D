function u = Solve(A, f)

arguments (Input)
    A NUDFT2;
    f (:, :) double;
end

arguments (Output)
    u (:, :) double;
end

f = f(A.x_perm_, :);
c = URV_Solve(A.Factor_, f);
u = ifft(c);

end
