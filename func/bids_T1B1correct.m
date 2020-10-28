function bids_T1B1correct(BIDSroot, NrShots, EchoSpacing, Expression, subjects, InvEff, B1Scaling, Realign, Correct)

% FUNCTION bids_T1B1correct(BIDSroot, [NrShots], [EchoSpacing], [Expression], [subjects], [InvEff], [B1Scaling], [Realign], [Correct])
%
% A BIDS-aware wrapper ('bidsapp') around 'T1B1correctpackageTFL' function that reads and writes BIDS compliant data.
% The MP2RAGE images are assumed to be stored with a suffix in the filename (e.g. as "sub-001_acq-MP2RAGE_inv1.nii.gz").
% NB: Fieldmaps intended for MP2RAGE are not accomodated for.
%
% 'T1B1correctpackageTFL' removes residual B1 bias from T1-maps estimated from the MP2RAGE data as suggested in:
%
%   Marques, J.P., Gruetter, R., 2013. New Developments and Applications of the MP2RAGE Sequence -
%   Focusing the Contrast and High Spatial Resolution R1 Mapping. PLoS ONE 8. doi:10.1371/journal.pone.0069294
% 
% INPUT
%   BIDSroot        - The root directory of the BIDS repository with all the subject sub-directories
%   NrShots         - The number of shots in the inner loop, i.e. SlicesPerSlab * [PartialFourierInSlice-0.5 0.5].
%                     The json file doesn't usually reliably contain this information. Default = "ReconMatrixPE"
%   EchoSpacing     - The echo spacing in secs that typically is not given in the json file. Default = 2 * TE
%   Expression      - A structure with 'uni', 'inv1', 'inv2', 'B1map' and 'B1mag' fields for selecting the
%                     corresponding MP2RAGE and B1-map images. A suffix (e.g. '_uni') must be included.
%                     Default = struct('uni',  ['extra_data' filesep '*_uni.nii*'], ...
%                                      'inv1', ['extra_data' filesep '*_inv1.nii*'], ...
%                                      'inv2', ['extra_data' filesep '*_inv2.nii*'], ...
%                                      'B1map',['extra_data' filesep '*_B1map.nii*'], ...
%                                      'B1mag',['extra_data' filesep '*_mod-B1map_*magnitude.nii*']);
%                     The 'B1mag' field is only needed if Realign==true (see below)
%   subjects        - Directory list of BIDS subjects that are processed. All are subjects processed
%                     if left empty (default), i.e. then subjects = dir(fullfile(bidsroot, 'sub-*'))
%   InvEff          - The inversion efficiency of the adiabatic inversion. Ideally it should be 1 but in the first
%                     implementation of the MP2RAGE it was measured to be ~0.96. Default = 0.96
%   B1Scaling       - Relative scaling factor of the B1-map, i.e. the nr for which the B1-map the nominal
%                     flip-angle is the actual flip-angle. Default = 900
%   Realign         - Uses the magnitude image of the B1-scan to realign and reslice the B1-map to the space of the
%                     MP2RAGE image if Realign==true (default). Otherwise it is assumed this is already the case
%                     and this processing step is skipped. NB: Realign requires the SPM12 software on your Matlab-path
%   Correct         - If Correct==true (default) a B1 bias corrected UNI image is saved in BIDSroot/derivatives/MP2RAGE
%
% EXAMPLES
%   >> bids_T1B1correct('/project/3015046.06/bids')
%   >> bids_T1B1correct('/project/3015046.06/bids', [round(224*3/8) round(224*4/8)], 7.5e-3)
%   >> bids_T1B1correct('/project/3015046.06/bids', round(224*4/8), [], ...
%         struct('uni',  'anat/*_uni.nii*', ...
%                'inv1', 'anat/*_inv1.nii*', ...
%                'inv2', 'anat/*_inv2.nii*', ...
%                'B1map','extra_data/*_b1.nii*', ...
%                'B1mag','extra_data/*_mod-B1map_*magnitude.nii*'));
%   >> bids_T1B1correct('/project/3015046.06/bids', [], [], ...
%         struct('uni',  'extra_data/*_acq-Prot1_*_UNI.nii.gz', ...
%                'inv1', 'extra_data/*_acq-Prot1_*_INV1.nii.gz', ...
%                'inv2', 'extra_data/*_acq-Prot1_*_INV2.nii.gz', ...
%                'B1map','extra_data/*FLIPANGLEMAP*.nii.gz', ...
%                'B1mag','extra_data/*FLIPANGLEMAP*_magn.nii.gz'), ...
%         dir('/project/3015046.06/bids/sub-00*'));
%
% See also: DemoForR1Correction, T1B1correctpackageTFL, T1B1correctpackage
%
% Marcel Zwiers, 8/10/2020


%% Parse the input arguments
if nargin<4 || isempty(Expression)
    Expression = struct('uni',  ['extra_data' filesep '*_uni.nii*'], ...
                        'inv1', ['extra_data' filesep '*_inv1.nii*'], ...
                        'inv2', ['extra_data' filesep '*_inv2.nii*'], ...
                        'B1map',['extra_data' filesep '*_B1map.nii*'], ...
                        'B1mag',['extra_data' filesep '*_mod-B1map_*magnitude.nii*']);
end
if nargin<5 || isempty(subjects)
    subjects = dir(fullfile(BIDSroot, 'sub-*'));
end
if nargin<6 || isempty(InvEff)
    InvEff = 0.96;
end
if nargin<7 || isempty(B1Scaling)
    B1Scaling = 900;
end
if nargin<8 || isempty(Realign)
    Realign = true;
end
if nargin<9 || isempty(Correct)
    Correct = true;
end
assert(contains(Expression.uni, '_'), ...
    'The output will not be BIDS-compliant because the uni-expression "%s" does not seem to contain a suffix (e.g. "_uni")', Expression.uni)


%% Get all the MP2RAGE and B1map images
MP2RAGE = {};
HG = figure('Name', 'bids_T1B1correct');
for subject = subjects'
    
    sessions = dir(fullfile(subject.folder, subject.name, 'ses-*'));
    if isempty(sessions)
        sessions(1).folder = fullfile(subject.folder, subject.name);
        sessions(1).name   = '.';
    end
    
    for session = sessions'
        
        n = numel(MP2RAGE) + 1;
        fprintf('%s (%i):\n', fullfile(session.folder, session.name), n)
        
        uni      = dir(fullfile(session.folder, session.name, Expression.uni));
        inv1     = dir(fullfile(session.folder, session.name, Expression.inv1));
        inv2     = dir(fullfile(session.folder, session.name, Expression.inv2));
        B1map{n} = dir(fullfile(session.folder, session.name, Expression.B1map));
        B1mag{n} = dir(fullfile(session.folder, session.name, Expression.B1mag));
        if isempty(uni) || isempty(B1map{n}) || (Realign && isempty(B1mag{n}))
            fprintf('Could not find UNI & B1-map images with search terms:\n%s\n%s\n%s\n\n', fullfile(subject.name, session.name, Expression.uni), fullfile(subject.name, session.name, Expression.B1map), fullfile(subject.name, session.name, Expression.B1mag))
            continue
        elseif numel(uni) > 1
            warning('Too many UNI-images found in:\n%s\n%s\n', uni.folder, sprintf('%s\n', uni.name))
            continue
        elseif numel(B1map{n}) > 1
            warning('Too many B1map-images found in:\n%s\n%s\n', B1map{n}.folder, sprintf('%s\n', B1map{n}.name))
            continue
        elseif Realign && numel(B1mag{n}) > 1
            warning('Too many B1map magnitude-images found in:\n%s\n%s\n', B1mag{n}.folder, sprintf('%s\n', B1mag{n}.name))
            continue
        end
        
        MP2RAGE{n}   = PopulateMP2RAGEStructure(uni, inv1, inv2, EchoSpacing, NrShots);

        suffix       = split(strtok(Expression.uni,'.'), '_');              % ASSUMPTION ALERT: The MP2RAGE image is stored with a (custom) suffix
        T1mapname{n} = fullfile(session.folder, session.name, 'anat', strrep(uni.name, ['_' suffix{end}], '_T1map'));   % Corrected B1-map
        
        fprintf('%s\n%s\n%s\n%s\n--> %s\n\n', uni.name, inv1.name, inv2.name, B1map{n}.name, T1mapname{n})
        
        % Check the properties of this MP2RAGE protocol... this happens to be a very B1 insensitive protocol
        plotMP2RAGEproperties(MP2RAGE{n}, HG)
    
    end
    
end


%% Process all the images
for n = 1:numel(MP2RAGE)
    
    fprintf('\n--> Processing: %s (%i/%i\n)', MP2RAGE{n}.filenameUNI, n, numel(MP2RAGE))
    
    % Realign and reslice the B1-map to the INV2 image
    if Realign
        B1Src     = spm_vol_gz(fullfile(B1map{n}.folder, B1map{n}.name));
        B1SrcMag  = spm_vol_gz(fullfile(B1mag{n}.folder, B1mag{n}.name));
        INV2Ref   = spm_vol_gz(MP2RAGE{n}.filenameINV2);
        x         = spm_coreg(INV2Ref, B1SrcMag);               % Coregister B1SrcMag with INV2Ref
        R         = B1Src.mat \ spm_matrix(x) * INV2Ref.mat;    % R = Mapping from voxels in INV2Ref to voxels in B1Src
        B1img.img = NaN(INV2Ref.dim);
        for z = 1:INV2Ref.dim(3)                                % Reslice the B1Src volume at the coordinates of each coregistered transverse slice of INV2Ref
            B1img.img(:,:,z) = spm_slice_vol(B1Src, R * spm_matrix([0 0 z]), INV2Ref.dim(1:2), 1);
        end
    else
        B1img = load_untouch_nii(fullfile(B1map{n}.folder, B1map{n}.name));
    end
    
    % Perform the unbiased B1-map estimation and the UNI image correction
    B1img.img            = double(B1img.img) / B1Scaling;
    MP2RAGEimg           = load_untouch_nii(MP2RAGE{n}.filenameUNI);
    [T1map, MP2RAGECorr] = T1B1correctpackageTFL(B1img, MP2RAGEimg, [], MP2RAGE{n}, [], InvEff);
    
    % Save the T1-map image
    save_untouch_nii(T1map, T1mapname{n})
    
    % Read & enrich the UNI json-file and write it as a T1map json-file
    [UNIpath, UNIname, UNIext]    = myfileparts(MP2RAGE{n}.filenameUNI);
    jsonT1map                     = jsondecode(fileread(fullfile(UNIpath, [UNIname '.json'])));
    jsonT1map.SeriesDescription   = [jsonT1map.ProtocolName '_B1_bias_corrected'];
    jsonT1map.InversionEfficiency = InvEff;
    jsonT1map.NumberOfShots       = MP2RAGE{n}.NZslices;
    jsonT1map.EchoSpacing         = MP2RAGE{n}.TRFLASH;
    jsonT1map.InversionTime       = MP2RAGE{n}.TIs;
    jsonT1map.FlipAngle           = MP2RAGE{n}.FlipDegrees;
    
    [T1path, T1name] = myfileparts(T1mapname{n});
    fid = fopen(fullfile(T1path, [T1name '.json']), 'w');
    fprintf(fid, '%s', jsonencode(jsonT1map));
    fclose(fid);

    % Save the corrected UNI image & json file
    if Correct
        MP2RAGECorrpath = strrep(UNIpath, fullfile(BIDSroot,filesep), fullfile(BIDSroot,'derivatives','MP2RAGE',filesep));
        MP2RAGECorrname = [UNIname 'B1corr'];
        [~,~]           = mkdir(MP2RAGECorrpath);
        save_untouch_nii(MP2RAGECorr, fullfile(MP2RAGECorrpath, [MP2RAGECorrname UNIext]))
        
        fid = fopen(fullfile(MP2RAGECorrpath, [MP2RAGECorrname '.json']), 'w');
        fprintf(fid, '%s', jsonencode(jsonT1map));
        fclose(fid);
    end
    
    % Adapt the scans.tsv file
    subses = split(UNIname, '_');
    if contains(UNIname, '_ses-')
        subses = subses(1:2);
    else
        subses = subses(1);
    end
    scansfile = fullfile(BIDSroot, subses{:}, sprintf('%sscans.tsv', sprintf('%s_', subses{:})));
    if isfile(scansfile)
        scanstable  = readtable(scansfile, 'FileType','text', 'ReadRowNames',true, 'Delimiter','\t', 'PreserveVariableNames',true);
        [~, source] = fileparts(UNIpath);
        UNIscan     = [source '/' UNIname UNIext];
        if any(contains(scanstable.Properties.RowNames, UNIscan))
            UNIdata = scanstable(UNIscan,:).Variables;
        else
            UNIdata = repmat({'n/a'}, 1, size(scanstable.Variables, 2));
        end
        [~, T1mapfname, T1mapext] = myfileparts(T1mapname{n});
        T1mapscan                 = ['anat/' T1mapfname T1mapext];
        scanstable(T1mapscan, :)  = UNIdata;
        fprintf('Updating %s:\n--> %s%s\n\n', scansfile, T1mapscan, sprintf('\t%s',UNIdata{:}))
        writetable(scanstable, scansfile, 'FileType','text', 'WriteRowNames',true, 'Delimiter','\t')
    end
    
end


function MP2RAGEstructure = PopulateMP2RAGEStructure(uni, inv1, inv2, EchoSpacing, NrShots)
%
%   uni         - The directory item of the UNI file
%   inv1        - The directory item of the INV1 file
%   inv2        - The directory item of the INV2 file
%   EchoSpacing - The echo spacing in secs that typically is not given on the json file. Default: twice the echo time
%   NrShots     - The number of shots in the inner loop, the json file doesn't usually accomodate this. Default: ReconMatrixPE

jsonINV1 = jsondecode(fileread(fullfile(inv1.folder, [strtok(inv1.name,'.') '.json'])));
jsonINV2 = jsondecode(fileread(fullfile(inv2.folder, [strtok(inv2.name,'.') '.json'])));

MP2RAGEstructure.B0           =  jsonINV1.MagneticFieldStrength;                 % In Tesla
MP2RAGEstructure.TR           =  jsonINV1.RepetitionTime;                        % MP2RAGE TR in seconds
MP2RAGEstructure.TIs          = [jsonINV1.InversionTime jsonINV2.InversionTime]; % Inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGEstructure.FlipDegrees  = [jsonINV1.FlipAngle     jsonINV2.FlipAngle];     % Flip angle of the two readouts in degrees
MP2RAGEstructure.filenameUNI  = fullfile(uni.folder, uni.name);                  % Standard MP2RAGE T1w image
MP2RAGEstructure.filenameINV1 = fullfile(inv1.folder, inv1.name);
MP2RAGEstructure.filenameINV2 = fullfile(inv2.folder, inv2.name);

if nargin<4 || isempty(EchoSpacing)
    MP2RAGEstructure.TRFLASH = jsonINV1.EchoTime * 2;       % TR of the GRE readout in seconds
else
    MP2RAGEstructure.TRFLASH = EchoSpacing;
end

if nargin<5 || isempty(NrShots)
    assert(isfield(jsonINV1,'ReconMatrixPE'), 'The json-file does not contain "NrSHots"-info (i.e. "ReconMatrixPE") belonging to:\n%s', MP2RAGEstructure.filenameINV1)
    MP2RAGEstructure.NZslices = jsonINV1.ReconMatrixPE;     % Slices Per Slab * [PartialFourierInSlice-0.5 0.5] OR Base Resolution * [PartialFourierInPE-0.5 0.5]/iPATpe + [RefLines/2 RefLines/2]*(1-1/iPATpe )
else
    MP2RAGEstructure.NZslices = NrShots;
end


function Vol = spm_vol_gz(FileName)
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


function [pathname, filename, ext] = myfileparts(filename)
%
% Robust against .nii.gz file extension

[pathname, filename, ext2] = fileparts(filename);
[~, filename, ext1]        = fileparts(filename);
ext                        = [ext1 ext2];
