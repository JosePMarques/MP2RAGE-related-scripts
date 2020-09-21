function bids_RobustCombination(bidsroot, regularization, expression, source, target)

% FUNCTION bids_RobustCombination(bidsroot, regularization)
%
% INPUT
%   bidsroot        - The root of the BIDS directory with all the subject directories
%   regularization  - A noise level regularization term, see also: DemoRemoveBackgroundNoise.m
%
% A BIDS-aware (a 'bidsapp') wrapper around 'RobustCombination' that reads and writes BIDS
% compliant data
%
% Marcel Zwiers, 18/9/2020

if nargin<2
    regularization = [];
end
if nargin<3 || isempty(expression)
    expression = struct('uni', ['extra_data' filesep '*_uni.nii*'], ...
                        'inv1',['extra_data' filesep '*_inv1.nii*'], ...
                        'inv2',['extra_data' filesep '*_inv2.nii*']);
end
if nargin<4 || isempty(source)
    source = 'extra_data';
end
if nargin<4 || isempty(target)
    target = 'anat';
end

% Get all the MP2RAGE images
MP2RAGE  = [];
subjects = dir(fullfile(bidsroot, 'sub-*'));
for subject = subjects'
    sessions = dir(fullfile(bidsroot, subject.name, 'ses-*'));
    if isempty(sessions)
        sessions(1).name = '.';
    end
    for session = sessions'
        for expr = expression(:)'
            fprintf('%s:\n', fullfile(bidsroot, subject.name, session.name))
            uni  = dir(fullfile(bidsroot, subject.name, session.name, expr.uni));
            inv1 = dir(fullfile(bidsroot, subject.name, session.name, expr.inv1));
            inv2 = dir(fullfile(bidsroot, subject.name, session.name, expr.inv2));
            if isempty(uni) || isempty(inv1) || isempty(inv2)
                fprintf('Could not find UNI, INV1 & INV2 images with search terms:\n%s\n%s\n%s\n\n', fullfile(subject.name, session.name, expr.uni), fullfile(subject.name, session.name, expr.inv1),  fullfile(subject.name, session.name, expr.inv2))
                continue
            elseif numel(uni)>1 || numel(inv1)>1 || numel(inv2)>1
                warning('Too many UNI, INV1 & INV2 images found in:\n%s\n%s\n', fullfile(subject.name, session.name, expr.uni), sprintf('%s\n', uni.name, inv1.name, inv2.name))
                continue
            end
            index                       = numel(MP2RAGE) + 1;
            MP2RAGE(index).filenameUNI  = fullfile(uni.folder, uni.name);                               % standard MP2RAGE T1w image;
            MP2RAGE(index).filenameINV1 = fullfile(inv1.folder, inv1.name);                             % Inversion Time 1 MP2RAGE T1w image;
            MP2RAGE(index).filenameINV2 = fullfile(inv2.folder, inv2.name);                             % Inversion Time 2 MP2RAGE T1w image;
            suffix                      = split(expr.uni, '_');                                         % ASSUMPTION ALERT: The MP2RAGE image is stored with a (custom) suffix
            T1name                      = strrep(uni.name, ['_' strtok(suffix{end}, '.')], '_T1w');     % Guess the suffix from the search expression
            if strcmp(target, 'derivatives')
                MP2RAGE(index).filenameOUT = fullfile(bidsroot, 'derivatives', 'MP2RAGE', subject.name, session.name, T1name);      % T1w image without background noise;
            else
                MP2RAGE(index).filenameOUT = fullfile(bidsroot, subject.name, session.name, target, T1name);                        % T1w image without background noise;
            end
            fprintf('%s\n%s\n%s\n--> %s\n\n', uni.name, inv1.name, inv2.name, MP2RAGE(index).filenameOUT)
        end
    end
end

% Get a good regularization value from the first MP2RAGE image
if isempty(regularization)
    [~, regularization] = RobustCombination(rmfield(MP2RAGE(1),'filenameOUT'), regularization, true);
end

% Process all the MP2RAGE images
for n = 1:numel(MP2RAGE)
    if ~exist(fileparts(MP2RAGE(n).filenameOUT), 'dir')
        mkdir(fileparts(MP2RAGE(n).filenameOUT))
    end
    RobustCombination(MP2RAGE(n), regularization, false);
end
