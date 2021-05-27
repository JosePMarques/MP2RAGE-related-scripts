function bids_RobustCombination(bidsroot, regularization, expression, subjects, target)

% FUNCTION bids_RobustCombination(bidsroot, [regularization], [expression], [subjects], [target])
%
% A BIDS-aware wrapper ('bidsapp') around 'RobustCombination' that reads and writes BIDS compliant
% data. The MP2RAGE images are assumed to be stored with a suffix in the filename (e.g. as
% "sub-001_acq-MP2RAGE_inv1.nii.gz") or stored according to the BIDS v1.5 standard (e.g.
% "sub-001_inv-1_MP2RAGE.nii.gz"). NB: Fieldmaps intended for MP2RAGE are not accomodated for.
%
% INPUT
%   bidsroot        - The root directory of the BIDS repository with all the subject sub-directories
%   regularization  - A noise level regularization term, (default = find level interactively)
%   expression      - A structure with 'uni', 'inv1' and 'inv2' search fields for selecting the
%                     corresponding MP2RAGE images in the sub-directory. A suffix needs to be included
%                     in the uni-expression (e.g. '_UNIT1').
%                     Default = struct('uni',  ['anat' filesep '*_UNIT1.nii*'], ...
%                                      'inv1', ['anat' filesep '*_inv-1*_MP2RAGE.nii*'], ...
%                                      'inv2', ['anat' filesep '*_inv-2*_MP2RAGE.nii*']);
%                     NB: Paths before '::' are prepended to the sub-directories, such that
%                     'derivatives/SIEMENS::anat/*_UNIT1.nii*' will perform a search for UNIT1 images in
%                     e.g. 'bidsroot/derivatives/SIEMENS/sub-01/ses-01/anat/' and 'anat/*_inv-1*_MP2RAGE.nii*'
%                     will search for MP2RAGE images in e.g. 'bidsroot/sub-01/ses-01/anat/'
%   subjects        - Directory list of BIDS subjects that are processed. All subjects are processed if
%                     left empty (default), i.e. then subjects = dir(fullfile(bidsroot, 'sub-*'))
%   target          - The target sub-directory in which the combined image is saved, e.g. 'anat'
%                     (default = 'derivatives')
%
% EXAMPLES
%   >> bids_RobustCombination('/project/3015046.06/bids')
%   >> bids_RobustCombination('/project/3015046.06/bids', 15)
%   >> bids_RobustCombination('/project/3015046.06/bids', [], ...
%         struct('uni','anat/*_uni.nii*', 'inv1','anat/*_inv1.nii*', 'inv2','anat/*_inv2.nii*'))
%   >> bids_RobustCombination('/project/3015046.06/bids', [], ...
%         struct('uni', 'derivatives/SIEMENS::anat/*UNIimages*_MP2RAGE.nii.gz', ...
%                'inv1','anat/*INV1*_MP2RAGE.nii.gz', ...
%                'inv2','anat/*INV2*_MP2RAGE.nii.gz'))
%   >> bids_RobustCombination('/project/3015046.06/bids', [], ...
%         struct('uni', 'extra_data/*_acq-Prot1_*_UNI.nii*', ...
%                'inv1','extra_data/*_acq-Prot1_*_INV1.nii*', ...
%                'inv2','extra_data/*_acq-Prot1_*_INV2.nii*'))
%   >> bids_RobustCombination('/project/3015046.06/bids', [], [], dir('/project/3015046.06/bids/sub-00*'), 'derivatives')
%
% See also: DemoRemoveBackgroundNoise, RobustCombination
%
% Marcel Zwiers, 17/03/2021


%% Parse the input arguments
if nargin<2
    regularization = [];
end
if nargin<3 || isempty(expression)
    expression = struct('uni',  ['anat' filesep '*_UNIT1.nii*'], ...
                        'inv1', ['anat' filesep '*_inv-1*_MP2RAGE.nii*'], ...
                        'inv2', ['anat' filesep '*_inv-2*_MP2RAGE.nii*']);
end
if nargin<4 || isempty(subjects)
    subjects = dir(fullfile(bidsroot, 'sub-*'));
end
if nargin<5 || isempty(target)
    target = 'derivatives';
end

assert(contains(expression.uni, '_'), 'The output will not be bids-compliant because the uni-expression "%s" does not seem to contain a suffix (e.g. "_UNIT1")', expression.uni)
suffix = split(expression.uni, '_');                                   % ASSUMPTION ALERT: The MP2RAGE image is stored with a (custom) suffix
suffix = suffix{end};


%% Get all the MP2RAGE images
MP2RAGE = [];
for subject = subjects'
    
    sessions = dir(fullfile(subject.folder, subject.name, 'ses-*'));
    if isempty(sessions)
        sessions(1).folder = fullfile(subject.folder, subject.name);
        sessions(1).name   = '.';
    end
    
    for session = sessions'
        
        fprintf('%s:\n', fullfile(session.folder, session.name))
        
        uni  = getdata(session, expression.uni);
        inv1 = getdata(session, expression.inv1);
        inv2 = getdata(session, expression.inv2);
        if isempty(uni) || isempty(inv1) || isempty(inv2)
            fprintf('Could not find UNI, INV1 & INV2 images with search terms:\n%s\n%s\n%s\n\n', fullfile(subject.name, session.name, expression.uni), fullfile(subject.name, session.name, expression.inv1),  fullfile(subject.name, session.name, expression.inv2))
            continue
        elseif ~isequal(numel(uni), numel(inv1), numel(inv2))
            warning('Unequal number of UNI (%i), INV1 (%i) & INV2 (%i) images found in:\n%s\n', numel(uni), numel(inv1), numel(inv2), fullfile(subject.name, session.name, expression.uni))
            continue
        end

        for n = 1:numel(uni)
            index                       = numel(MP2RAGE) + 1;
            MP2RAGE(index).filenameUNI  = fullfile(uni(n).folder, uni(n).name);                    % Standard MP2RAGE T1w image
            MP2RAGE(index).filenameINV1 = fullfile(inv1(n).folder, inv1(n).name);                  % Inversion Time 1 MP2RAGE T1w image
            MP2RAGE(index).filenameINV2 = fullfile(inv2(n).folder, inv2(n).name);                  % Inversion Time 2 MP2RAGE T1w image
            T1name                      = strrep(uni(n).name, ['_' strtok(suffix,'.')], '_T1w');   % Guess the suffix from the search expression
            if strcmp(target, 'derivatives')
                MP2RAGE(index).filenameOUT = fullfile(bidsroot, 'derivatives', 'MP2RAGE', subject.name, session.name, 'anat', T1name);  % T1w image without background noise;
            else
                MP2RAGE(index).filenameOUT = fullfile(session.folder, session.name, target, T1name);                                    % T1w image without background noise;
            end

            fprintf('%s\n%s\n%s\n--> %s\n\n', uni(n).name, inv1(n).name, inv2(n).name, MP2RAGE(index).filenameOUT)
        end

    end
    
end


%% Get a good regularization value from the first MP2RAGE image
HG = figure('Name', 'bids_RobustCombination');
if isempty(regularization) && ~isempty(MP2RAGE)
    [~, regularization] = RobustCombination(rmfield(MP2RAGE(1),'filenameOUT'), [], HG);
end


%% Process all the MP2RAGE images
for n = 1:numel(MP2RAGE)
    
    % Save a combined image
    [outpath, outname, outext] = fileparts(MP2RAGE(n).filenameOUT);
    if ~isfolder(outpath)
        mkdir(outpath)
    end
    RobustCombination(MP2RAGE(n), regularization, HG);

    % Adapt the scans.tsv file
    if ~strcmp(target, 'derivatives')
    
        subses = split(outname,'_');
        if contains(outname, '_ses-')
            subses = subses(1:2);
        else
            subses = subses(1);
        end
        
        scansfile = fullfile(bidsroot, subses{:}, sprintf('%sscans.tsv', sprintf('%s_', subses{:})));
        if isfile(scansfile)
            scanstable                 = readtable(scansfile, 'FileType','text', 'ReadRowNames',true, 'Delimiter','\t', 'PreserveVariableNames',true, 'DatetimeType','text');
            [UNIpath, UNIname, UNIext] = fileparts(MP2RAGE(n).filenameUNI);
            [~, source]                = fileparts(UNIpath);            
            UNIscan                    = [source '/' UNIname UNIext];
            if any(contains(scanstable.Properties.RowNames, UNIscan))
                UNIdata = scanstable(UNIscan,:).Variables;
            else
                UNIdata = repmat({'n/a'}, 1, size(scanstable.Variables, 2));
            end
            T1scan = [target '/' outname outext];
            scanstable(T1scan, :) = UNIdata;
            fprintf('Updating %s:\n--> %s%s\n\n', scansfile, T1scan, sprintf('\t%s',UNIdata{:}))
            writetable(scanstable, scansfile, 'FileType','text', 'WriteRowNames',true, 'Delimiter','\t')
        end
    
    end
    
end


function data = getdata(session, expression)

expression = split(expression, '::');
if numel(expression)>1
    bidsroot       = split(session.folder, 'sub-');
    session.folder = fullfile(bidsroot{1}, expression{1}, ['sub-' bidsroot{2}]);
end

data = dir(fullfile(session.folder, session.name, expression{end}));
