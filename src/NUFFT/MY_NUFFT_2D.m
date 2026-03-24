function f = MY_NUFFT_2D(c, xy, omega)

f = exp(-2 * pi * 1i * xy * omega.') * c;

end