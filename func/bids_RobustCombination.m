function bids_RobustCombination(bidsroot, regularization, expression, subjects, target)

% FUNCTION bids_RobustCombination(bidsroot, [regularization], [expression], [subjects], [target])
%
% A BIDS-aware wrapper ('bidsapp') around 'RobustCombination' that reads and writes BIDS compliant
% data. The MP2RAGE images are assumed to be stored with a suffix in the filename (e.g. as
% "sub-001_acq-MP2RAGE_inv1.nii.gz"). NB: Fieldmaps intended for MP2RAGE are not accomodated for.
%
% INPUT
%   bidsroot        - The root directory of the BIDS repository with all the subject sub-directories
%   regularization  - A noise level regularization term, (default = find level interactively)
%   expression      - A structure with 'uni', 'inv1' and 'inv2' fields for selecting the
%                     corresponding MP2RAGE images. The suffix (e.g. '_uni') needs to be included.
%                     Default = struct('uni', ['extra_data' filesep '*_uni.nii*'], ...
%                                      'inv1',['extra_data' filesep '*_inv1.nii*'], ...
%                                      'inv2',['extra_data' filesep '*_inv2.nii*'])
%   subjects        - Directory list of BIDS subjects that are processed. All are subjects processed
%                     if left empty (default), i.e. then subjects = dir(fullfile(bidsroot, 'sub-*'))
%   target          - The target sub-directory in which the combined image is saved (default = 'anat')
%                     If 'derivatives' then the output is saved in the bids 'derivatives' root-folder
%
% EXAMPLES
%   >> bids_RobustCombination('/project/3015046.06/bids')
%   >> bids_RobustCombination('/project/3015046.06/bids', 15)
%   >> bids_RobustCombination('/project/3015046.06/bids', [], ...
%         struct('uni','anat/*_uni.nii*', 'inv1','anat/*_inv1.nii*', 'inv2','anat/*_inv2.nii*'))
%   >> bids_RobustCombination('/project/3015046.06/bids', [], ...
%         struct('uni', 'anat/*UNIimages*_MP2RAGE.nii.gz', ...
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
% Marcel Zwiers, 24/9/2020


%% Parse the input arguments
if nargin<2
    regularization = [];
end
if nargin<3 || isempty(expression)
    expression = struct('uni', ['extra_data' filesep '*_uni.nii*'], ...
                        'inv1',['extra_data' filesep '*_inv1.nii*'], ...
                        'inv2',['extra_data' filesep '*_inv2.nii*']);
end
if nargin<4 || isempty(subjects)
    subjects = dir(fullfile(bidsroot, 'sub-*'));
end
if nargin<5 || isempty(target)
    target = 'anat';
end
assert(contains(expression.uni, '_'), ...
    'The output will not be bids-compliant because the uni-expression "%s" does not seem to contain a suffix (e.g. "_inv1")', expression.uni)


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
        
        uni  = dir(fullfile(session.folder, session.name, expression.uni));
        inv1 = dir(fullfile(session.folder, session.name, expression.inv1));
        inv2 = dir(fullfile(session.folder, session.name, expression.inv2));
        if isempty(uni) || isempty(inv1) || isempty(inv2)
            fprintf('Could not find UNI, INV1 & INV2 images with search terms:\n%s\n%s\n%s\n\n', fullfile(subject.name, session.name, expression.uni), fullfile(subject.name, session.name, expression.inv1),  fullfile(subject.name, session.name, expression.inv2))
            continue
        elseif numel(uni)>1 || numel(inv1)>1 || numel(inv2)>1
            warning('Too many UNI, INV1 & INV2 images found in:\n%s\n%s\n', fullfile(subject.name, session.name, expression.uni), sprintf('%s\n', uni.name, inv1.name, inv2.name))
            continue
        end
        
        index                       = numel(MP2RAGE) + 1;
        MP2RAGE(index).filenameUNI  = fullfile(uni.folder, uni.name);                               % Standard MP2RAGE T1w image
        MP2RAGE(index).filenameINV1 = fullfile(inv1.folder, inv1.name);                             % Inversion Time 1 MP2RAGE T1w image
        MP2RAGE(index).filenameINV2 = fullfile(inv2.folder, inv2.name);                             % Inversion Time 2 MP2RAGE T1w image
        suffix                      = split(expression.uni, '_');                                   % ASSUMPTION ALERT: The MP2RAGE image is stored with a (custom) suffix
        T1name                      = strrep(uni.name, ['_' strtok(suffix{end},'.')], '_T1w');      % Guess the suffix from the search expression
        if strcmp(target, 'derivatives')
            MP2RAGE(index).filenameOUT = fullfile(bidsroot, 'derivatives', 'MP2RAGE', subject.name, session.name, T1name);      % T1w image without background noise;
        else
            MP2RAGE(index).filenameOUT = fullfile(session.folder, session.name, target, T1name);                                % T1w image without background noise;
        end
        
        fprintf('%s\n%s\n%s\n--> %s\n\n', uni.name, inv1.name, inv2.name, MP2RAGE(index).filenameOUT)        

    end
    
end


%% Get a good regularization value from the first MP2RAGE image
HG = figure('Name', 'bids_RobustCombination');
if isempty(regularization)
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
            scanstable                 = readtable(scansfile, 'FileType','text', 'ReadRowNames',true, 'Delimiter','\t', 'PreserveVariableNames',true);
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
