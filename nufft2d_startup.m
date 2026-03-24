function nufft2d_startup()
% nufft2d_startup

file_path = mfilename('fullpath');
tmp = strfind(file_path, 'nufft2d');
file_path = file_path(1:(tmp(end)-1));
addpath(genpath([file_path 'src']));
addpath(genpath([file_path 'extern']));

end