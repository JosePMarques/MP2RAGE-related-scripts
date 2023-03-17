function bids_T1B1correct_job(BIDSroot, InvEff, B1Scaling, Realign, FWHM, Target, ...
    Correct, Fingerprint, B1correctM0, MP2RAGE, B1map, B1Ref, R1mapname)

% Performs the work for bids_T1B1correct()
%
% Marcel Zwiers, 17/03/2023

Debug = 0;

% Load the headers & data
B1Src     = spm_vol_gz(fullfile(B1map.folder, B1map.name));
B1Src_Vol = spm_read_vols(B1Src);
if Realign || FWHM ~= 0
    B1Ref_ = spm_vol_gz(fullfile(B1Ref.folder, B1Ref.name));
end

% Smooth the B1-map in order to avoid influence of salt & pepper border noise
if FWHM ~= 0
    B1Ref_Vol   = spm_read_vols(B1Ref_);
    PixDim      = spm_imatrix(B1Src.mat);
    Pre_Smooth  = double(B1Ref_Vol).^1 .* exp(1i*double(B1Src_Vol) / B1Scaling);
    Post_Smooth = smooth3D(Pre_Smooth, FWHM, abs(PixDim(7:9)));
    B1Src_Vol   = angle(Post_Smooth) * B1Scaling;
    B1Src.fname = spm_file(B1Src.fname, 'suffix', '_smooth');
    B1Src       = spm_write_vol(B1Src, B1Src_Vol);
end

% Realign & reslice the B1 reference image with the INV2 image
INV2Ref = [];
if Realign
    INV2Ref   = spm_vol_gz(MP2RAGE.filenameINV2);
    x         = spm_coreg(INV2Ref, B1Ref_);                 % Coregister B1Ref with INV2Ref
    R         = B1Src.mat \ spm_matrix(x) * INV2Ref.mat;    % R = Mapping from voxels in INV2Ref to voxels in B1Src
    B1img.img = NaN(INV2Ref.dim);
    for z = 1:INV2Ref.dim(3)                                % Reslice the B1Src volume at the coordinates of each coregistered transverse slice of INV2Ref
        B1img.img(:,:,z) = spm_slice_vol(B1Src, R * spm_matrix([0 0 z]), INV2Ref.dim(1:2), 1);
    end
else
    B1img.img = B1Src_Vol;
end

if Debug
    figure(452)
    subplot(121)
    Orthoview(B1img.img / B1Scaling, [], [-0.8 1.1])
    title('After Coregisteration and Smoothing')
    subplot(122)
    Orthoview(B1Src_Vol / B1Scaling, [], [-0.8 1.1])
    title('Before Coregisteration')
end

% Perform the unbiased B1-map estimation and the UNI image correction
B1img.img  = double(B1img.img) / B1Scaling;
MP2RAGEimg = load_untouch_nii(MP2RAGE.filenameUNI);
MP2RAGESrc = spm_vol_gz(MP2RAGE.filenameINV2);       % A bit redundant to unzip the image twice, but hey

% Clean-up the temporarily unzipped/smoothed images
for TempVol = [B1Src, B1Ref_, INV2Ref, MP2RAGESrc]
    if startsWith(TempVol.fname, tempdir) && isfile(TempVol.fname)
        delete(TempVol.fname)
        if endsWith(TempVol.fname, '_smooth.nii')
            delete(strrep(TempVol.fname, '_smooth.nii', '.nii'))
        end
    end
end

if Fingerprint == false

    [~, MP2RAGECorr] = T1B1correctpackageTFL(B1img, MP2RAGEimg, [], MP2RAGE, [], InvEff);

    % Compute the M0- and R1-map
    INV2img            = load_untouch_nii(MP2RAGE.filenameINV2);
    [~, M0map , R1map] = T1M0estimateMP2RAGE(MP2RAGECorr, INV2img, MP2RAGE, InvEff);

else

    INV1img = load_untouch_nii(MP2RAGE.filenameINV1);
    INV2img = load_untouch_nii(MP2RAGE.filenameINV2);
    [INV1img, INV2img] = Correct_INV1INV2_withMP2RAGEuni(INV1img, INV2img, MP2RAGEimg, 0);
    MP2RAGE.invEff  =  InvEff;
    [~, M0map.img, R1map.img] = MP2RAGE_dictionaryMatching(MP2RAGE, INV1img.img, INV2img.img, B1img.img, [0.002,0.005], 1, B1img.img ~= 0);

    if Correct    
         [MP2RAGE.Intensity, MP2RAGE.T1vector] = MP2RAGE_lookuptable(2, MP2RAGE.TR, MP2RAGE.TIs, MP2RAGE.FlipDegrees, MP2RAGE.NZslices, MP2RAGE.TRFLASH, 'normal', MP2RAGE.invEff);
         MP2RAGECorr     = MP2RAGEimg;
         MP2RAGECorr.img = reshape(interp1(MP2RAGE.T1vector, MP2RAGE.Intensity, 1./R1map.img(:)), size(R1map.img));
         MP2RAGECorr.img(isnan(MP2RAGECorr.img)) =- 0.5;
         MP2RAGECorr.img = round(4095*(MP2RAGECorr.img + 0.5));
    end

end

if B1correctM0 ~= 0
    M0map.img = M0map.img ./ flipdim(B1img.img, B1correctM0);
end

% data is only valid where B1 was mapped
R1map.img(B1img.img == 0) = 0;
M0map.img(B1img.img == 0) = 0;

% Save the R1-map image
R1Hdr       = MP2RAGESrc;
R1Hdr.fname = R1mapname;
spm_write_vol_gz(R1Hdr, R1map.img);

% Read & enrich the UNI json-file and write it as a R1-map json-file
[UNIpath, UNIname, UNIext]    = myfileparts(MP2RAGE.filenameUNI);
jsonR1map                     = jsondecode(fileread(fullfile(UNIpath, [UNIname '.json'])));
jsonR1map.BasedOn             = {MP2RAGE.filenameUNI, MP2RAGE.filenameINV1, MP2RAGE.filenameINV2, fullfile(B1map.folder, B1map.name), fullfile(B1Ref.folder, B1Ref.name)};
jsonR1map.SeriesDescription   = [jsonR1map.ProtocolName '_B1_bias_corrected'];
jsonR1map.InversionEfficiency = InvEff;
jsonR1map.NumberShots         = MP2RAGE.NZslices;
jsonR1map.EchoSpacing         = MP2RAGE.TRFLASH;
jsonR1map.InversionTime       = MP2RAGE.TIs;
jsonR1map.FlipAngle           = MP2RAGE.FlipDegrees;
fid = fopen(spm_file(spm_file(R1Hdr.fname,'ext',''), 'ext','.json'), 'w');
fprintf(fid, '%s', jsonencode(jsonR1map));
fclose(fid);

% Save the M0-map & json file
M0Hdr       = MP2RAGESrc;
M0Hdr.fname = strrep(R1mapname, '_R1map.nii', '_M0map.nii');
spm_write_vol_gz(M0Hdr, M0map.img);
fid = fopen(spm_file(spm_file(M0Hdr.fname,'ext',''), 'ext','.json'), 'w');
fprintf(fid, '%s', jsonencode(jsonR1map));
fclose(fid);

% Save the B1-map & json file
B1Hdr       = MP2RAGESrc;
B1Hdr.fname = strrep(R1mapname, '_R1map.nii', '_B1map.nii');
spm_write_vol_gz(B1Hdr, B1img.img);
fid = fopen(spm_file(spm_file(B1Hdr.fname,'ext',''), 'ext','.json'), 'w');
fprintf(fid, '%s', jsonencode(jsonR1map));
fclose(fid);

% Save the corrected UNI image & json file
if Correct
    CorrHdr       = MP2RAGESrc;
    CorrHdr.fname = strrep(R1mapname, '_R1map.nii', '_desc-B1corr_UNIT1.nii');
    spm_write_vol_gz(CorrHdr, MP2RAGECorr.img);
    fid = fopen(spm_file(spm_file(CorrHdr.fname,'ext',''), 'ext','.json'), 'w');
    fprintf(fid, '%s', jsonencode(jsonR1map));
    fclose(fid);
end

% Adapt the scans.tsv file
if ~startsWith(Target, 'derivatives')
    SubSes = split(UNIname, '_');
    if contains(UNIname, '_ses-')
        SubSes = SubSes(1:2);
    else
        SubSes = SubSes(1);
    end
    scansfile = fullfile(BIDSroot, SubSes{:}, sprintf('%sscans.tsv', sprintf('%s_', SubSes{:})));
    if isfile(scansfile)
        ScansTable  = readtable(scansfile, 'FileType','text', 'ReadRowNames',true, 'Delimiter','\t', 'PreserveVariableNames',true, 'DatetimeType','text');
        [~, Source] = fileparts(UNIpath);
        UNIscan     = [Source '/' UNIname UNIext];
        if any(contains(ScansTable.Properties.RowNames, UNIscan))
            UNIdata = ScansTable(UNIscan,:).Variables;
        else
            UNIdata = repmat({'n/a'}, 1, size(ScansTable.Variables, 2));
        end
        [~, R1mapfname, R1mapext] = myfileparts(R1mapname);
        R1mapscan                 = ['anat/' R1mapfname R1mapext];
        ScansTable(R1mapscan, :)  = UNIdata;
        fprintf('Updating %s:\n--> %s%s\n\n', scansfile, R1mapscan, sprintf('\t%s',UNIdata{:}))
        writetable(ScansTable, scansfile, 'FileType','text', 'WriteRowNames',true, 'Delimiter','\t')
    end
end
