function bids_T1B1correct(BIDSroot, NrShots, EchoSpacing, Expression, subjects, InvEff, B1Scaling, Realign, Target, Correct)

% FUNCTION bids_T1B1correct(BIDSroot, [NrShots], [EchoSpacing], [Expression], [subjects], [InvEff], [B1Scaling], [Realign], [Target], [Correct])
%
% A BIDS-aware wrapper ('bidsapp') around 'T1B1correctpackageTFL' function that reads and writes BIDS compliant data.
% The MP2RAGE images are assumed to be stored with a suffix in the filename (e.g. as "sub-001_acq-MP2RAGE_inv1.nii.gz"
% or as BIDS v1.5 images, e.g. as "sub-001_inv-1__MP2RAGE.nii.gz").
% NB: Fieldmaps intended for MP2RAGE are not accomodated for.
%
% 'T1B1correctpackageTFL' removes residual B1 bias from T1-maps estimated from the MP2RAGE data as suggested in:
%
%   Marques, J.P., Gruetter, R., 2013. New Developments and Applications of the MP2RAGE Sequence -
%   Focusing the Contrast and High Spatial Resolution R1 Mapping. PLoS ONE 8. doi:10.1371/journal.pone.0069294
%
%   Marques, J.P., Kober, T., Krueger, G., van der Zwaag, W., Van de Moortele, P., Gruetter, R., 2010.
%   MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field.
%   NeuroImage 49, doi.org/10.1016/j.neuroimage.2009.10.002
% 
% INPUT
%   BIDSroot        - The root directory of the BIDS repository with all the subject sub-directories
%   NrShots         - The number of shots in the inner loop, i.e. SlicesPerSlab * [PartialFourierInSlice-0.5 0.5].
%                     The json file doesn't usually / reliably contain this information, but it should be available
%                     from the scan protocol (adviced). Default: NrShots = nr of slices in the z-dimension.
%                     NrShots has previously also been referred to as `NZslices`.
%   EchoSpacing     - The RepetitionTimeExcitation value in secs that is not always given in the json file, but it
%                     should be available from the scan protocol (adviced). Default: EchoSpacing = 2*TE
%                     EchoSpacing has previously also been referred to as `TRFLASH`.
%   Expression      - A structure with 'uni', 'inv1', 'inv2', 'B1map' and 'B1Ref' search fields for selecting the
%                     corresponding MP2RAGE images in the sub-directory. A suffix needs to be included in the
%                     uni-expression (e.g. '_UNIT1'). The B1map/Ref fields must each select 1 image for correction
%                     of the MP2RAGE images. The 'B1Ref' field is only used if Realign==true (see further below)
%                     Default = struct('uni',   ['anat' filesep '*_UNIT1.nii*'], ...
%                                      'inv1',  ['anat' filesep '*_inv-1*_MP2RAGE.nii*'], ...
%                                      'inv2',  ['anat' filesep '*_inv-2*_MP2RAGE.nii*'], ...
%                                      'B1map', ['fmap' filesep '*_acq-famp*_TB1*.nii*'], ...
%                                      'B1Ref', ['fmap' filesep '*_acq-anat*_TB1*.nii*']);
%                     NB: Paths before '::' are prepended to the sub-directories, such that
%                     'derivatives/SIEMENS::anat/*_UNIT1.nii*' will perform a search for UNIT1 images in
%                     e.g. 'bidsroot/derivatives/SIEMENS/sub-01/ses-01/anat/' and 'anat/*_inv-1*_MP2RAGE.nii*'
%                     will search for MP2RAGE images in e.g. 'bidsroot/sub-01/ses-01/anat/'
%   subjects        - Directory list of BIDS subjects that are processed. All subjects are processed
%                     if left empty (default), i.e. then subjects = dir(fullfile(bidsroot, 'sub-*'))
%   InvEff          - The inversion efficiency of the adiabatic inversion. Ideally it should be 1 but in the first
%                     implementation of the MP2RAGE it was measured to be ~0.96. Default = 0.96
%   B1Scaling       - Relative scaling factor of the B1-map, i.e. the nr for which the B1-map the nominal
%                     flip-angle is the actual flip-angle. Default = 900
%   Realign         - Uses the magnitude image of the B1-scan to realign and reslice the B1-map to the space of the
%                     MP2RAGE image if Realign==true (default). Otherwise it is assumed this is already the case
%                     and this processing step is skipped. NB: Realign requires the SPM12 software on your Matlab-path
%   Target          - The target sub-directory in which the corrected B1-map is saved, e.g. 'anat'
%                     (default = 'derivatives')
%   Correct         - If Correct==true (default) a B1 bias corrected UNI image is saved in BIDSroot/derivatives/MP2RAGE
%
% EXAMPLES
%   >> bids_T1B1correct('/project/3015046.06/bids')
%   >> bids_T1B1correct('/project/3015046.06/bids', [round(224*3/8) round(224*4/8)], 7.5e-3)    % NB: NrSHots and EchoSpacing are typically available from the scan protocol and are adviced to be used
%   >> bids_T1B1correct('/project/3015046.06/bids', round(224*4/8), [], ...
%         struct('uni',  'derivatives/SIEMENS::anat/*_uni.nii*', ...
%                'inv1', 'anat/*_inv1.nii*', ...
%                'inv2', 'anat/*_inv2.nii*', ...
%                'B1map','extra_data/*_b1.nii*'), ...
%         [], [], [], false);
%   >> bids_T1B1correct('/project/3015046.06/bids', [], [], ...
%         struct('uni',  'extra_data/*_acq-Prot1_*_UNI.nii.gz', ...
%                'inv1', 'extra_data/*_acq-Prot1_*_INV1.nii.gz', ...
%                'inv2', 'extra_data/*_acq-Prot1_*_INV2.nii.gz', ...
%                'B1map','fmap/*_TB1map.nii.gz', ...
%                'B1Ref','fmap/*_TB1TFL.nii.gz'), ...
%         dir('/project/3015046.06/bids/sub-00*'));
%
% See also: DemoForR1Correction, T1B1correctpackageTFL, T1B1correctpackage
%
% Marcel Zwiers, 17/03/2021


%% Parse the input arguments
if nargin<2
    NrShots = [];
end
if nargin<3
    EchoSpacing = [];
end
if nargin<4 || isempty(Expression)
    Expression = struct('uni',   ['anat' filesep '*_UNIT1.nii*'], ...
                        'inv1',  ['anat' filesep '*_inv-1*_MP2RAGE.nii*'], ...
                        'inv2',  ['anat' filesep '*_inv-2*_MP2RAGE.nii*'], ...
                        'B1map', ['fmap' filesep '*_acq-famp*_TB1*.nii*'], ...
                        'B1Ref', ['fmap' filesep '*_acq-anat*_TB1*.nii*']);
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
if nargin<9 || isempty(Target)
    Target = 'derivatives';
end
if nargin<10 || isempty(Correct)
    Correct = true;
end

assert(contains(Expression.uni, '_'), 'The output will not be BIDS-compliant because the uni-expression "%s" does not seem to contain a suffix (e.g. "_UNIT1")', Expression.uni)
suffix = split(strtok(Expression.uni,'.'), '_');              % ASSUMPTION ALERT: The MP2RAGE image is stored with a (custom) suffix
suffix = suffix{end};


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
        
        fprintf('Indexing (%i): %s\n', numel(MP2RAGE) + 1, fullfile(session.folder, session.name))
        
        uni    = getfiles(session, Expression.uni);
        inv1   = getfiles(session, Expression.inv1);
        inv2   = getfiles(session, Expression.inv2);
        B1map_ = getfiles(session, Expression.B1map);
        B1Ref_ = getfiles(session, Expression.B1Ref);
        if isempty(uni) || isempty(B1map_)
            fprintf('Could not find UNI & B1-map images with search terms "%s" and "%s"\n', Expression.uni, Expression.B1map)
            continue
        elseif Realign && isempty(B1Ref_)
            fprintf('Could not find B1-reference image for realignment with search term "%s\n', Expression.B1Ref)
            continue
        elseif ~isequal(numel(uni), numel(inv1), numel(inv2))
            warning('Unequal number of UNI (%i), INV1 (%i) and INV2 (%i) images found with search terms "%s", "%s" and "%s"\n', numel(uni), numel(inv1), numel(inv2), expression.uni, Expression.inv1, Expression.inv2)
            continue
        elseif numel(B1map_) ~= 1
            warning('Ambiguous (%i instead of 1) B1map-images found using "%s"\n', numel(B1map_), Expression.B1map)
            disp(char(B1map_.name))
            continue
        elseif Realign && numel(B1Ref_) ~= 1
            warning('Ambiguous (%i instead of 1) B1Ref-images found using "%s"\n', numel(B1Ref_), Expression.B1Ref)
            disp(char(B1Ref_.name))
            continue
        end
        
        for n = 1:numel(uni)
            
            index        = numel(MP2RAGE) + 1;
            [MP2RAGE{index}, EchoSpacing, NrShots] = PopulateMP2RAGEStructure(uni(n), inv1(n), inv2(n), EchoSpacing, NrShots);
            B1map{index} = B1map_;         % NB: The same B1-map is used for all MP2RAGE images
            B1Ref{index} = B1Ref_;         % NB: The same B1-ref is used for all MP2RAGE images            
            if strcmp(Target, 'derivatives')
                T1mapname{index} = fullfile(BIDSroot, 'derivatives', 'MP2RAGE', subject.name, session.name, 'anat', strrep(uni(n).name, ['_' suffix], '_T1map'));   % Corrected T1-map
            else
                T1mapname{index} = fullfile(session.folder, session.name, Target, strrep(uni(n).name, ['_' suffix], '_T1map'));   % Corrected T1-map
            end

            fprintf('%s\n%s\n%s\n%s\n--> %s\n\n', uni(n).name, inv1(n).name, inv2(n).name, B1map{index}.name, T1mapname{index})

            % Check the properties of this MP2RAGE protocol... this happens to be a very B1 insensitive protocol
            plotMP2RAGEproperties(MP2RAGE{index}, HG)

        end
    end
    
end


%% Process all the images
for n = 1:numel(MP2RAGE)
    
    fprintf('\n--> Processing (%i/%i): %s\n', n, numel(MP2RAGE), MP2RAGE{n}.filenameUNI)
    
    % Realign and reslice the B1-map to the INV2 image
    if Realign
        B1Src     = spm_vol_gz(fullfile(B1map{n}.folder, B1map{n}.name));
        B1Ref_    = spm_vol_gz(fullfile(B1Ref{n}.folder, B1Ref{n}.name));
        INV2Ref   = spm_vol_gz(MP2RAGE{n}.filenameINV2);
        x         = spm_coreg(INV2Ref, B1Ref_);                 % Coregister B1Ref with INV2Ref
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
    [T1path, T1name] = myfileparts(T1mapname{n});
    [~,~]            = mkdir(T1path);
    save_untouch_nii(T1map, T1mapname{n})
    
    % Read & enrich the UNI json-file and write it as a T1map json-file
    [UNIpath, UNIname, UNIext]    = myfileparts(MP2RAGE{n}.filenameUNI);
    jsonT1map                     = jsondecode(fileread(fullfile(UNIpath, [UNIname '.json'])));
    jsonT1map.BasedOn             = {MP2RAGE{n}.filenameUNI, MP2RAGE{n}.filenameINV1, MP2RAGE{n}.filenameINV2, fullfile(B1map{n}.folder, B1map{n}.name), fullfile(B1Ref{n}.folder, B1Ref{n}.name)};
    jsonT1map.SeriesDescription   = [jsonT1map.ProtocolName '_B1_bias_corrected'];
    jsonT1map.InversionEfficiency = InvEff;
    jsonT1map.NumberShots         = MP2RAGE{n}.NZslices;
    jsonT1map.EchoSpacing         = MP2RAGE{n}.TRFLASH;
    jsonT1map.InversionTime       = MP2RAGE{n}.TIs;
    jsonT1map.FlipAngle           = MP2RAGE{n}.FlipDegrees;
    
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
        scanstable  = readtable(scansfile, 'FileType','text', 'ReadRowNames',true, 'Delimiter','\t', 'PreserveVariableNames',true, 'DatetimeType','text');
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


function [MP2RAGEstructure, EchoSpacing, NrShots] = PopulateMP2RAGEStructure(uni, inv1, inv2, EchoSpacing, NrShots)
%
%   uni         - The directory item of the UNI file
%   inv1        - The directory item of the INV1 file
%   inv2        - The directory item of the INV2 file
%   EchoSpacing - The RepetitionTimeExcitation value in secs that typically is not given on the json file. Default: twice the echo time
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
    if isfield(jsonINV1, 'RepetitionTimeExcitation')
        EchoSpacing = jsonINV1.RepetitionTimeExcitation;                         % TR of the GRE readout in seconds
    else
        EchoSpacing = jsonINV1.EchoTime * 2;                                     % 2 X EchoTime can be used as a surrogate
    end
    disp(['Extracted EchoSpacing: ' num2str(EchoSpacing)])
end
MP2RAGEstructure.TRFLASH = EchoSpacing;

if nargin<5 || isempty(NrShots)
    UNIimg  = load_untouch_nii(MP2RAGEstructure.filenameUNI);                    % A bit overkill to load the whole image but loading just the header is only available as an internal function :-(
    NrShots = size(UNIimg.img, 3);
    disp(['Extracted NrShots: ' num2str(NrShots)])
end
MP2RAGEstructure.NZslices = NrShots;


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


function data = getfiles(session, expression)
%
% Uses dir(expression) to return the matching files in the session

expression = split(expression, '::');
if numel(expression)>1
    bidsroot       = split(session.folder, 'sub-');
    session.folder = fullfile(bidsroot{1}, expression{1}, ['sub-' bidsroot{2}]);
end

data = dir(fullfile(session.folder, session.name, expression{end}));
