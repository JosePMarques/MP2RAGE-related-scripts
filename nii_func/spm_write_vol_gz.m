function Vol = spm_write_vol_gz(Hdr, data)
% FUNCTION Vol = spm_write_vol_gz(Hdr, data)
%
% A wrapper around spm_write_vol that writes and zips .nii files

[fpath, ~, Ext] = myfileparts(Hdr.fname);
Hdr             = rmfield(Hdr, 'pinfo');
[~,~]           = mkdir(fpath);
switch Ext
    case '.nii.gz'
        Hdr.fname = spm_file(Hdr.fname, 'ext','');
        Vol       = spm_write_vol(Hdr, data);
        gzip(Hdr.fname)
        delete(Hdr.fname)
    case '.nii'
        Vol = spm_write_vol(Hdr, data);
    otherwise
        error('Unknown file extenstion %s in %s', Ext, FileName)
end
