function bids_T1B1correct(bidsroot, NrShots, EchoSpacing, expression, invEFF, B1Scaling, realign)

% FUNCTION bids_T1B1correct(bidsroot, [NrShots], [EchoSpacing], [expression], [invEFF], [B1Scaling], [realign])
%
% A BIDS-aware wrapper ('bidsapp') around 'T1B1correctpackageTFL' that reads and writes BIDS compliant data.
% The MP2RAGE images are assumed to be stored with a suffix in the filename (e.g. as "sub-001_acq-MP2RAGE_inv1.nii.gz"). 
%
% The 'T1B1correctpackageTFL' removes residual B1 bias from T1-maps estimated from the MP2RAGE data as suggested in:
%
%   Marques, J.P., Gruetter, R., 2013. New Developments and Applications of the MP2RAGE Sequence -
%   Focusing the Contrast and High Spatial Resolution R1 Mapping. PLoS ONE 8. doi:10.1371/journal.pone.0069294
% 
% In this script it is assumed that the B1 maps have been coregistered to the space of the MP2RAGE
% image and that they have the B1 has in the process been interpolated to the same resolution.
%
% INPUT
%   bidsroot        - The root of the BIDS directory with all the subject directories
%   NrShots         - The number of shots in the inner loop, the json file doesn't usually contain this ("ReconMatrixPE")
%   EchoSpacing     - The echo spacing in secs that typically is not given in the json file. If it is not
%                     given we take it to be twice the echo time
%   expression      - A structure with 'uni', 'inv1', 'inv2' and 'B1map' fields for selecting the
%                     corresponding MP2RAGE and B1-map images. A suffix (e.g. '_uni') must be included.
%                     Default = struct('uni',  ['extra_data' filesep '*_uni.nii*'], ...
%                                      'inv1', ['extra_data' filesep '*_inv1.nii*'], ...
%                                      'inv2', ['extra_data' filesep '*_inv2.nii*'], ...
%                                      'B1map',['extra_data' filesep '*_B1map.nii*'])
%   invEFF          - Default = 0.96
%   B1Scaling       - Relative scaling factor of the B1-map, i.e. the nr for which the B1-map the nominal
%                     flip-angle is the actual flip-angle. Default = 900
%   realign         - Realigns and reslices the B1-map to the space of the UNI image if realign==True (default).
%                     Otherwise it is assumed this is already the case and this processing step is skipped.
%                     NB: The realign option requires SPM12 software package on your Matlab-path
%
% EXAMPLES
%   >> bids_T1B1correct('/project/3015046.06/bids')
%
% See also: DemoForR1Correction, T1B1correctpackageTFL, T1B1correctpackage
%
% Marcel Zwiers, 25/9/2020


%% Parse the input arguments
if nargin<4 || isempty(expression)
    expression = struct('uni',  ['extra_data' filesep '*_uni.nii*'], ...
                        'inv1', ['extra_data' filesep '*_inv1.nii*'], ...
                        'inv2', ['extra_data' filesep '*_inv2.nii*'], ...
                        'B1map',['extra_data' filesep '*_B1map.nii*']);
end
assert(contains(expression.uni, '_'), ...
    'The output will not be bids-compliant because the uni-expression "%s" does not seem to contain a suffix (e.g. "_uni")', expression.uni)
if nargin<5 || isempty(invEFF)
    invEFF = 0.96;
end
if nargin<6 || isempty(B1Scaling)
    B1Scaling = 900;
end


%% Get all the MP2RAGE images
MP2RAGE  = [];
subjects = dir(fullfile(bidsroot, 'sub-*'));
for subject = subjects'
    
    sessions = dir(fullfile(subject.folder, subject.name, 'ses-*'));
    if isempty(sessions)
        sessions(1).folder = fullfile(subject.folder, subject.name);
        sessions(1).name   = '.';
    end
    
    for session = sessions'
        
        n = numel(MP2RAGE) + 1;
        fprintf('%s (%i):\n', fullfile(session.folder, session.name), n)
        
        uni      = dir(fullfile(session.folder, session.name, expression.uni));
        inv1     = dir(fullfile(session.folder, session.name, expression.inv1));
        inv2     = dir(fullfile(session.folder, session.name, expression.inv2));
        B1map{n} = dir(fullfile(session.folder, session.name, expression.B1map));
        if isempty(uni) || isempty(B1map{n})
            fprintf('Could not find UNI & B1-map images with search terms:\n%s\n%s\n\n', fullfile(subject.name, session.name, expression.uni),  fullfile(subject.name, session.name, expression.B1map))
            continue
        elseif numel(uni) > 1
            warning('Too many UNI-images found in:\n%s\n%s\n', uni.folder, sprintf('%s\n', uni.name))
            continue
        elseif numel(B1map{n}) > 1
            warning('Too many B1map-images found in:\n%s\n%s\n', B1map{n}.folder, sprintf('%s\n', B1map{n}.name))
            continue
        end
        
        MP2RAGE(n).filenameUNI  = fullfile(uni.folder, uni.name);           % Standard MP2RAGE T1w image
        MP2RAGE(n).filenameINV1 = fullfile(inv1.folder, inv1.name);
        MP2RAGE(n).filenameINV2 = fullfile(inv2.folder, inv2.name);
        MP2RAGE(n)              = PopulateMP2RAGEStructure(MP2RAGE(n), EchoSpacing, NrShots);
        
        suffix       = split(strtok(expression.uni,'.'), '_');              % ASSUMPTION ALERT: The MP2RAGE image is stored with a (custom) suffix
        T1mapname{n} = fullfile(session.folder, session.name, 'anat', strrep(uni.name, ['_' suffix{end}], '_T1map'));   % Corrected B1-map
        
        fprintf('%s\n%s\n%s\n%s\n--> %s\n\n', uni.name, inv1.name, inv2.name, B1map{n}.name, T1mapname{n})
        
        % Check the properties of this MP2RAGE protocol... this happens to be a very B1 insensitive protocol
        plotMP2RAGEproperties(MP2RAGE(n))
    
    end
    
end


%% Process all the images
for n = 1:numel(MP2RAGE)
    
    fprintf('\n--> Processing: %s (%i/%i\n)', MP2RAGE(n).filenameUNI, n, numel(MP2RAGE))
    
    % Get the UNI and T1map fileparts (robust for .nii.gz file extensions)
    [UNIbase, UNIext]     = strtok(MP2RAGE(n).filenameUNI, '.');
    [UNIpath, UNIfname]   = fileparts(UNIbase);
    [T1mapbase, T1mapext] = strtok(T1mapname{n}, '.');
    [~, T1mapfname]       = fileparts(T1mapbase);
    
    % Realign and reslice the B1-map to the UNI image
    if realign
        B1Src     = spm_vol_gz(fullfile(B1map{n}.folder, B1map{n}.name));
        UNIRef    = spm_vol_gz(MP2RAGE(n).filenameUNI);
        x         = spm_coreg(UNIRef, B1Src);
        R         = B1Src.mat \ spm_matrix(x) * UNIRef.mat;     % R = Mapping from voxels in UNIRef to voxels in B1Src
        B1img.img = NaN(Ref.dim);
        for z = 1:Ref.dim(3)
            B1img.img(:,:,z) = spm_slice_vol_gz(B1Src, R * spm_matrix([0 0 z]), UNIRef.dim(1:2), 1);
        end
    else
        B1img = load_untouch_nii(fullfile(B1map{n}.folder, B1map{n}.name));
    end
    
    % Perform the correction
    B1img.img  = double(B1img.img) / B1Scaling;
    MP2RAGEimg = load_untouch_nii(MP2RAGE(n).filenameUNI);
    T1map      = T1B1correctpackageTFL(B1img, MP2RAGEimg, [], MP2RAGE(n), [], invEFF);
    
    % Save the T1-map image and copy-over & enrich the UNI json-file
    save_untouch_nii(T1map, T1mapname{n})
    copyfile(fullfile(UNIbase,'.json'), [strtok(T1mapname{n},'.') '.json'])
    
    % Adapt the IntendedFor fieldmap values (TODO)
    
    % Adapt the scans.tsv file
    subses = split(UNIfname, '_');
    if contains(UNIfname, '_ses-')
        subses = subses(1:2);
    else
        subses = subses(1);
    end
    scansfile = fullfile(bidsroot, subses{:}, sprintf('%sscans.tsv', sprintf('%s_', subses{:})));
    if exist(scansfile, 'file')
        scanstable  = readtable(scansfile, 'FileType','text', 'ReadRowNames',true, 'Delimiter','\t', 'PreserveVariableNames',true);
        [~, source] = fileparts(UNIpath);
        UNIscan     = [source '/' UNIfname UNIext];
        if any(contains(scanstable.Properties.RowNames, UNIscan))
            UNIdata = scanstable(UNIscan,:).Variables;
        else
            UNIdata = repmat({'n/a'}, 1, size(scanstable.Variables, 2));
        end
        T1mapscan                = ['anat/' T1mapfname T1mapext];
        scanstable(T1mapscan, :) = UNIdata;
        fprintf('Updating %s:\n--> %s%s\n\n', scansfile, T1mapscan, sprintf('\t%s',UNIdata{:}))
        writetable(scanstable, scansfile, 'FileType','text', 'WriteRowNames',true, 'Delimiter','\t')
    end
    
end


function MP2RAGEstructure = PopulateMP2RAGEStructure(MP2RAGEstructure, EchoSpacing, NrShots)

% MP2RAGEstructure = PopulateMP2RAGEStructure(MP2RAGEstructure, EchoSpacing, NrShots)
%
% INPUT
%   MP2RAGEstructure.filenameUNI  - the UNI file
%   MP2RAGEstructure.filenameINV1 - the INV1 file
%   MP2RAGEstructure.filenameINV2 - the INV2 file
%   EchoSpacing  - the echo spacing in secs that typically is not given on the
%                  json file. if it is not given we take it to be twice the echo time
%   NrShots      - refers to the number of shots in the inner loop, the json
%                  file doesn't usually accomodate this

jsonINV1 = jsondecode(fileread([strtok(MP2RAGEstructure.filenameINV1,'.') '.json']));
jsonINV2 = jsondecode(fileread([strtok(MP2RAGEstructure.filenameINV2,'.') '.json']));

MP2RAGEstructure.B0          = jsonINV1.MagneticFieldStrength;                  % in Tesla
MP2RAGEstructure.TR          = jsonINV1.RepetitionTime;                         % MP2RAGE TR in seconds
MP2RAGEstructure.TIs         = [jsonINV1.InversionTime jsonINV2.InversionTime]; % inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGEstructure.FlipDegrees = [jsonINV1.FlipAngle jsonINV2.FlipAngle];         % Flip angle of the two readouts in degrees

if nargin<2 || isempty(EchoSpacing)
    MP2RAGEstructure.TRFLASH = jsonINV1.EchoTime * 2;       % TR of the GRE readout
else
    MP2RAGEstructure.TRFLASH = EchoSpacing;
end

if nargin<3 || isempty(NrShots)
    assert(isfield(jsonINV1,'ReconMatrixPE'), 'The json-file does not contain "NrSHots"-info (i.e. "ReconMatrixPE") beloning to:\n%s', MP2RAGEstructure.filenameINV1)
    MP2RAGEstructure.NZslices = jsonINV1.ReconMatrixPE;     % Slices Per Slab * [PartialFourierInSlice-0.5 0.5] OR Base Resolution * [PartialFourierInPE-0.5 0.5]/iPATpe + [RefLines/2 RefLines/2]*(1-1/iPATpe )
else
    MP2RAGEstructure.NZslices = NrShots;
end
