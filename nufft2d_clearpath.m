function nufft2d_clearpath()
% nufft2d_clearpath

file_path = mfilename('fullpath2d');
tmp = strfind(file_path, 'nufft');
file_path = file_path(1:(tmp(end)-1));
rmpath(genpath([file_path 'src']));
rmpath(genpath([file_path 'extern']));

end