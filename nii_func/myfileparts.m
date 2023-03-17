function [pathname, filename, ext] = myfileparts(filename)
% FUNCTION [pathname, filename, ext] = myfileparts(filename)
%
% Robust against .nii.gz file extension

[pathname, filename, ext2] = fileparts(filename);
[~, filename, ext1]        = fileparts(filename);
ext                        = [ext1 ext2];
