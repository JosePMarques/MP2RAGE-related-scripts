function Vol = spm_vol_gz(FileName)
% FUNCTION Vol = spm_vol_gz(FileName)
%
% A wrapper around spm_vol that unzips .nii.gz files in a temp-folder

[~, ~, Ext] = myfileparts(FileName);
switch Ext
    case '.nii.gz'
        FileName = char(gunzip(FileName, tempdir));
    case '.nii'
    otherwise
        error('Unknown file extenstion %s in %s', Ext, FileName)
end

Vol = spm_vol(FileName);
