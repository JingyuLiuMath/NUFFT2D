function URV_Factor(A)

arguments (Input)
    A NUDFT2;
end

A.Factor_ = A.AFinv_HSS_.URV_Factor();

end
