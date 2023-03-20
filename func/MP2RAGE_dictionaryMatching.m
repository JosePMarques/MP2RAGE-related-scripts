function [T1, PD, R1] = MP2RAGE_dictionaryMatching(MP2RAGE, INV1,INV2,B1map,varargin)
% function that computes T1 and PD maps using a dictionary matching
% approach.
% input variables are:
% MP2RAGE is a structure that should have the following fields
%       MP2RAGE.TR
%       MP2RAGE.TIs,
%       MP2RAGE.NZslices,
%       MP2RAGE.TRFLASH,
%       MP2RAGE.invEff
% INV1  is the Inv1 image, it should be a real value (positive and negative
%       values) its intensities should have been corrected using "" function
% INV2  is the Inv1 image, it should be a real value (positive and negative
%       values) its intensities should have been corrected using "" function
% B1map  is the normalized B1 map, it should have alredy been coregistered
%       to the MP2RAGE spaceimages, it can be left empty
% varargin{1}(1) represents the step size to of the R1 vector
% varargin{1}(2) represents the step size to of the B1 vector ()
% varargin{2} is the input that defines if the T1 maps will be discretized (0) or
% interpolated (1) the interpolation requires is done by default using the
% three highest correlation points of the dictionary matching;
% varargin{3} is a mask to be used (only at massk==1 will the match be performed)

% in the future B1vector or R1 vector could be provided as an alternatice

ProbabilityEstimate = 1;
interpPt = 3; % numebr of points to use on interpolation 

dims = size(INV1);
mask = ones(dims);

if isempty(B1map)
    B1map = ones(dims);
end
if ~isfield(MP2RAGE,'invEff')
    MP2RAGE.invEff = 0.96;
end
if nargin == 4
    deltaR1 = 0.01;
    deltaB1 = 0.01;
else
    if nargin >=5
        if length(varargin{1})==1
            deltaR1 = varargin{1};
            deltaB1 = 0.01;
        else
            deltaR1 = varargin{1}(1);
            deltaB1 = varargin{1}(2);
        end
    end
    if nargin >= 6
        if ~isempty(varargin{2})
            ProbabilityEstimate = varargin{2};
        end
        if nargin >= 7
            mask = varargin{3};
        end
    end
end

R1vector = 0.2:deltaR1:5;
B1vector = min(B1map(mask==1)):deltaB1:max(B1map(mask==1));

% make the MP2RAGE data into a table
x = zeros(prod(dims),2);
x(:,1) = INV1(:);
x(:,2) = INV2(:);

% initialize variables
c = zeros(prod(dims),1);
idx = zeros(prod(dims),1);
T1 = zeros(prod(dims),1);
R1 = zeros(prod(dims),1);
PD = zeros(prod(dims),1);
% IntensityBeforeComb = zeros(length(R1vector),1);

%% make the matching per B1 value

count = 0;
fprintf('\nFinger printing B1-values: %f -> %f\n', B1vector([1 end]))
for B1 = B1vector

    count = count + 1;
    
    %% create dictionary for the specific B1 value
    j = 0;
    for R1val = R1vector
        j = j + 1;
        Signal(j,1:2) = 1*MPRAGEfunc(2, MP2RAGE.TR, MP2RAGE.TIs, MP2RAGE.NZslices, MP2RAGE.TRFLASH, B1*MP2RAGE.FlipDegrees, 'normal', 1./R1val, MP2RAGE.invEff);
    end

    dictionary = Signal(:,:) ./ vecnorm(Signal(:,:),2,2);
    %% dictionary Matching

    ind_B1 = find(and(and(B1map>=B1, B1map<B1+deltaB1), mask==1));
    if rem(count, 10) == 0
        fprintf('%f %i\n', B1, numel(ind_B1))
    else
        fprintf('%f %i, ', B1, numel(ind_B1))
    end
    % https://bitbucket.org/asslaender/nyu_mrf_recon/src/master/example/MRF_recon_example.m
    % Dictionary matching

    if ProbabilityEstimate == 0

        for q = ind_B1'
            [c(q,1), idx(q,1)] = max(x(q,:) * (dictionary'), [], 2);
        end
        PD(ind_B1) = c(ind_B1) ./ vecnorm(Signal(idx(ind_B1),:),2,2);
        % PD = c ./ D.normalization(idx).';
        % T1(ind_B1) = 1./R1vector(idx(ind_B1));
        R1(ind_B1) = R1vector(idx(ind_B1));

    else
        % initialize interpolation points
        weight = zeros(interpPt+1,1);
        % idxall = zeros(interpPt+1,1);     % Is defined again in the for loop below..

        for q = ind_B1'

            m = ones(1,length(dictionary));
            idxall = [0 0 0];
            [c((q),1), idx((q),1)] = max(x((q),:) * (dictionary'), [], 2);
            weight(1) = [c((q),1)];
            idxall(1) = idx((q),1);
            % look for following best dictionary matches
            m(idxall(1)) = 0;
            for n = 2:(interpPt+1)
                [weight(n), idxall(n)] = max((x((q),:)) * (dictionary') .* m, [], 2);
                m(idxall(n)) = 0;
            end
            % perform weighted average
            weightNorm(1:interpPt) = ((weight(1:interpPt)-weight(n)) ./ (weight(1)-weight(n))).^2;
            R1((q)) = sum(weightNorm.*R1vector(idxall(1:interpPt))) ./ sum(weightNorm);

        end % loop on pixels
        PD(ind_B1) = c(ind_B1) ./ vecnorm(Signal(idx(ind_B1),:),2,2);

    end     % loop on B1

end

PD = reshape(PD, dims);
R1 = reshape(R1, dims);
T1 = 1 ./ R1;
T1(mask==0) = 0;
