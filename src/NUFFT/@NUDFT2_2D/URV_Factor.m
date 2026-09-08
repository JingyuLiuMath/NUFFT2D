function URV_Factor(A)

arguments (Input)
    A NUDFT2_2D;
end

A.Factor_ = A.AFinv_HSS_.URV_Factor();

end
