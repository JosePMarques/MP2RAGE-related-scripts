function bids_T1B1correct(BIDSroot, NrShots, EchoSpacing, Expression, subjects, InvEff, B1Scaling, Realign, ...
    FWHM, Target, Correct, Fingerprint, B1correctM0, Cluster)

% FUNCTION bids_T1B1correct(BIDSroot, [NrShots], [EchoSpacing], [Expression], [subjects], [InvEff], [B1Scaling], ...
%    [Realign], [FWHM], [Target], [Correct], [Fingerprint], [B1correctM0], [Cluster])
%
% A BIDS-aware wrapper ('bidsapp') around 'T1B1correctpackageTFL' function that reads and writes BIDS compliant data.
% The MP2RAGE images are assumed to be stored with a suffix in the filename (e.g. as "sub-001_acq-MP2RAGE_inv1.nii.gz"
% or as BIDS v1.5 images, e.g. as "sub-001_inv-1__MP2RAGE.nii.gz").
% NB: Fieldmaps intended for MP2RAGE are not accomoda.gzted for.
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
%                                      'B1map', ['derivatives/SIEMENS::fmap' filesep '*_acq-famp*_TB1*.nii*'], ...
%                                      'B1Ref', ['derivatives/SIEMENS::fmap' filesep '*_acq-anat*_TB1*.nii*']);
%                     NB: Paths before '::' are prepended to the sub-directories, such that
%                     'derivatives/SIEMENS::anat/*_UNIT1.nii*' will perform a search for UNIT1 images in
%                     e.g. 'bidsroot/derivatives/SIEMENS/sub-01/ses-01/anat/' and 'anat/*_inv-1*_MP2RAGE.nii*'
%                     will search for MP2RAGE images in e.g. 'bidsroot/sub-01/ses-01/anat/'
%   subjects        - A wildcard expression to select the BIDS subjects that are processed. All subjects are
%                     processed if left empty (default), i.e. then subjects = 'sub-*'
%   InvEff          - The inversion efficiency of the adiabatic inversion. Ideally it should be 1 but in the first
%                     implementation of the MP2RAGE it was measured to be ~0.96. Default = 0.96
%   B1Scaling       - Relative scaling factor of the B1-map, i.e. the nr for which the B1-map the nominal
%                     flip-angle is the actual flip-angle. Default = 900
%   Realign         - Uses the magnitude image of the B1-scan to realign and reslice the B1-map to the space of the
%                     MP2RAGE image if Realign==true (default). Otherwise it is assumed this is already the case
%                     and this processing step is skipped. NB: Realign requires the SPM12 software on your Matlab-path
%   FWHM            - FWHM used to gaussian smooth the B1 map and ensure that skull regions are populated with valid
%                     B1 values. Default FWHM is 12. Pass FWHM = 0 if no smoothing is desired
%   Target          - The target sub-directory in which the corrected B1-map is saved, e.g. 'anat'
%                     (default = 'derivatives/MP2RAGE_scripts')
%   Correct         - If Correct==true (default) a B1 bias corrected UNI image is saved in the BIDS derivatives folder
%   Fingerprint     - If Fingerprint==true (not default) The T1 mapping and correction is performed using a
%                     fingerprinting approach
%   B1correctM0     - if 0 no correction is applied (default), 
%                     if 1,2,3 applies a flip on the left right direction of the transmit field as an approximation
%                     of the receive field. the number 1 2 3 refers to the left-right dimension only use this option
%                     if anatomical data is well centred in the field of view (as provided by auto-align) and if the
%                     Bias Field correction was enabled on the scanner
%   Cluster         - Uses Fieldtrip's qsubfeval function to submit the job(s) to the cluster if TRUE
%
% EXAMPLES
%   >> bids_T1B1correct('/project/3015046.06/bids')
%   >> bids_T1B1correct('/project/3015046.06/bids', [round(224*3/8) round(224*4/8)], 7.5e-3)    % NB: NrSHots and EchoSpacing are typically available from the scan protocol and are adviced to be used
%   >> bids_T1B1correct('/project/3015046.06/bids', round(224*4/8), [], ...
%         struct('uni',  'derivatives/SIEMENS::anat/*_uni.nii*', ...
%                'inv1', 'anat/*_inv1.nii*', ...
%                'inv2', 'anat/*_inv2.nii*', ...
%                'B1map','extra_data/*_b1.nii*'));
%   >> bids_T1B1correct('/project/3015046.06/bids', [], [], ...
%         struct('uni',  'extra_data/*_acq-Prot1_*_UNI.nii.gz', ...
%                'inv1', 'extra_data/*_acq-Prot1_*_INV1.nii.gz', ...
%                'inv2', 'extra_data/*_acq-Prot1_*_INV2.nii.gz', ...
%                'B1map','fmap/*_TB1map.nii.gz', ...
%                'B1Ref','fmap/*_TB1TFL.nii.gz'), ...
%         'sub-00*', [], [], false);
%
% See also: DemoForR1Correction, T1B1correctpackageTFL, T1B1correctpackage
%
% Marcel Zwiers, 10/02/2023

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
                        'B1map', ['derivatives/SIEMENS::fmap' filesep '*_acq-famp*_TB1*.nii*'], ...
                        'B1Ref', ['derivatives/SIEMENS::fmap' filesep '*_acq-anat*_TB1*.nii*']);
end
if nargin<5 || isempty(subjects)
    subjects = 'sub-*';
end
if ~contains(subjects, '*')
    subjects = [subjects '*'];
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
if nargin<9 || isempty(FWHM)
    FWHM = 12;
end
if nargin<10 || isempty(Target)
    Target = 'derivatives/MP2RAGE_scripts';
end
if nargin<11 || isempty(Correct)
    Correct = true;
end
if nargin<12 || isempty(Fingerprint)
    Fingerprint = false;
end
if nargin<13 || isempty(B1correctM0)
    B1correctM0 = 0;
end
if nargin<14 || isempty(Cluster)
    Cluster = false;
end

assert(contains(Expression.uni, '_'), 'The output will not be BIDS-compliant because the uni-expression "%s" does not seem to contain a suffix (e.g. "_UNIT1")', Expression.uni)
suffix = split(strtok(Expression.uni,'.'), '_');        % ASSUMPTION: The MP2RAGE image is stored with a unique suffix
suffix = suffix{end};


%% Get all the MP2RAGE and B1map images
MP2RAGE = {};
HG = figure('Name', 'bids_T1B1correct');
for subject = dir(fullfile(BIDSroot, subjects))'

    assert(startsWith(subject.name, 'sub-'), ['Your "' subjects '" subjects input argument did not select BIDS subject directories. See the help for usage'])
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
            fprintf('Could not find UNI & B1-map images with search terms "%s" and "%s"\n\n', Expression.uni, Expression.B1map)
            continue
        elseif Realign && isempty(B1Ref_)
            fprintf('Could not find B1-reference image for realignment with search term "%s\n\n', Expression.B1Ref)
            continue
        elseif ~isequal(numel(uni), numel(inv1), numel(inv2))
            warning('Unequal number of UNI (%i), INV1 (%i) and INV2 (%i) images found with search terms "%s", "%s" and "%s"\n\n', numel(uni), numel(inv1), numel(inv2), Expression.uni, Expression.inv1, Expression.inv2)
            continue
        elseif numel(B1map_) ~= 1
            warning('Ambiguous (%i instead of 1) B1map-images found using "%s"\n\n', numel(B1map_), Expression.B1map)
            disp(char(B1map_.name))
            continue
        elseif Realign && numel(B1Ref_) ~= 1
            warning('Ambiguous (%i instead of 1) B1Ref-images found using "%s"\n\n', numel(B1Ref_), Expression.B1Ref)
            disp(char(B1Ref_.name))
            continue
        end

        for n = 1:numel(uni)

            index        = numel(MP2RAGE) + 1;
            [MP2RAGE{index}, EchoSpacing, NrShots] = PopulateMP2RAGEStructure(uni(n), inv1(n), inv2(n), EchoSpacing, NrShots);
            B1map{index} = B1map_;         % NB: The same B1-map is used for all MP2RAGE images
            B1Ref{index} = B1Ref_;         % NB: The same B1-ref is used for all MP2RAGE images
            if startsWith(Target, 'derivatives')
                R1mapname{index} = fullfile(BIDSroot, Target, subject.name, session.name, 'anat', strrep(uni(n).name, ['_' suffix], '_R1map'));   % Corrected R1-map
            else
                R1mapname{index} = fullfile(session.folder, session.name, Target, strrep(uni(n).name, ['_' suffix], '_R1map'));   % Corrected R1-map
            end

            % fprintf('%s\n%s\n%s\n%s\n%s\n--> %s\n\n', uni(n).name, inv1(n).name, inv2(n).name, B1map{index}.name, B1Ref{index}.name, R1mapname{index})

            % Check the properties of this MP2RAGE protocol... this happens to be a very B1 insensitive protocol
            plotMP2RAGEproperties(MP2RAGE{index}, HG)

        end
    end

end


%% Process all the images
for n = 1:numel(MP2RAGE)
    fprintf('\n--> Processing (%i/%i): %s\n', n, numel(MP2RAGE), R1mapname{n})
    if Cluster
        qsubfeval('bids_T1B1correct_job', BIDSroot, InvEff, B1Scaling, Realign, FWHM, Target, Correct, ...
            Fingerprint, B1correctM0, MP2RAGE{n}, B1map{n}, B1Ref{n}, R1mapname{n}, 'timreq',60*60, 'memreq',4*1024^3);
    else
        bids_T1B1correct_job(BIDSroot, InvEff, B1Scaling, Realign, FWHM, Target, Correct, ...
            Fingerprint, B1correctM0, MP2RAGE{n}, B1map{n}, B1Ref{n}, R1mapname{n})
    end
end
if Cluster
    disp('Now wait for the jobs to finish... (see `man qstat`)')
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
    Hdr     = spm_vol_gz(MP2RAGEstructure.filenameUNI);
    NrShots = Hdr.dim(3);
    if startsWith(Hdr.fname, tempdir)
        delete(Hdr.fname)
    end
    disp(['Extracted NrShots: ' num2str(NrShots)])
end
MP2RAGEstructure.NZslices = NrShots;


function data = getfiles(session, expression)
%
% Uses dir(expression) to return the matching files in the session

expression = split(expression, '::');
if numel(expression)>1
    bidsroot       = split(session.folder, 'sub-');
    session.folder = fullfile(bidsroot{1}, expression{1}, ['sub-' bidsroot{2}]);
end

data = dir(fullfile(session.folder, session.name, expression{end}));
