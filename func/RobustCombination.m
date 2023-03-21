function [MP2RAGEimgRobustPhaseSensitive, multiplyingFactor] = RobustCombination(MP2RAGE, regularization, HG)

% FUNCTION [MP2RAGEimgRobustPhaseSensitive, multiplyingFactor] = RobustCombination(MP2RAGE, [regularization], [HG])
%
% This script creates MP2RAGE T1w images without the strong background noise in air regions.
%
% MP2RAGE is a structure that should have the following fields:
% MP2RAGE.filenameUNI
% MP2RAGE.filenameINV1
% MP2RAGE.filenameINV2
% MP2RAGE.filenameOUT - it does not have to exist, only if you want to save the output file.
%
% If you have already done your 'homework' with your datasets using the same protocol you can
% then just use this function shows one possible implementation of the methods suggested in:
%
% O'Brien, et al, 2014.
% Robust T1-Weighted Structural Brain Imaging and Morphometry at 7T Using MP2RAGE
% PLOS ONE 9, e99676. doi:10.1371/journal.pone.0099676
% http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099676
%
% Although in the original paper the method only worked on raw multichannel data, here that
% constraint has been overcome and the correction can be implemented if both SOS images of the
% two inversion times exist and a MP2RAGE T1w image that has been calculated directly from the
% multichannel data as initially proposed in Marques et al, Neuroimage, 2009


%% Parse the input arguments

if nargin<3
    HG = gcf;
    set(HG, 'Name','RobustCombination')
end
if nargin<2 || isempty(regularization)
    multiplyingFactor = 1;
    HG                = gcf;
else
    multiplyingFactor = regularization;
end
FinalChoice = 'n';


%% Define relevant functions

MP2RAGErobustfunc = @(INV1, INV2, beta)(conj(INV1).*INV2-beta) ./ (INV1.^2 + INV2.^2 + 2*beta);
rootsquares_pos   = @(a, b, c)(-b + sqrt(b.^2 - 4*a.*c)) ./ (2*a);
rootsquares_neg   = @(a, b, c)(-b - sqrt(b.^2 - 4*a.*c)) ./ (2*a);


%% Load Data

disp(['Loading images from: ' fileparts(MP2RAGE.filenameUNI)])
UNIhdr  = spm_vol_gz(MP2RAGE.filenameUNI);
UNIimg  = double(spm_read_vols(UNIhdr));
INV1hdr = spm_vol_gz(MP2RAGE.filenameINV1);
INV1img = double(spm_read_vols(INV1hdr));
INV2hdr = spm_vol_gz(MP2RAGE.filenameINV2);
INV2img = double(spm_read_vols(INV2hdr));

if min(UNIimg(:))>=0 && max(UNIimg(:))>=0.51
    % converts MP2RAGE to -0.5 to 0.5 scale - assumes that it is getting only positive values
    UNIimg = (UNIimg - max(UNIimg(:))/2) ./ max(UNIimg(:));
    integerformat = 1;
else
    integerformat = 0;
end


%% Compute correct INV1 dataset

% Gives the correct polarity to INV1
INV1img = sign(UNIimg) .* INV1img;

% Because the MP2RAGE INV1 and INV2 is a summ of squares data, while the
% MP2RAGEimg is a phase sensitive coil combination.. some more maths has to
% be performed to get a better INV1 estimate which here is done by assuming
% both INV2 is closer to a real phase sensitive combination

INV1pos = rootsquares_pos(-UNIimg, INV2img, -INV2img.^2 .* UNIimg);
INV1neg = rootsquares_neg(-UNIimg, INV2img, -INV2img.^2 .* UNIimg);

INV1final = INV1img;
INV1final(abs(INV1img-INV1pos) >  abs(INV1img-INV1neg)) = INV1neg(abs(INV1img-INV1pos) >  abs(INV1img-INV1neg));
INV1final(abs(INV1img-INV1pos) <= abs(INV1img-INV1neg)) = INV1pos(abs(INV1img-INV1pos) <= abs(INV1img-INV1neg));


%% Visualize the data

pos = round(3/5 * size(INV1final));
if isgraphics(HG)
    figureJ(200)
    subplot(411)
    Orthoview(INV1pos, pos, [-200 200])
    title('positive root')
    
    subplot(412)
    Orthoview(INV1neg, pos, [-200 200])
    title('negative root')
    
    subplot(413)
    Orthoview(INV1img, pos, [-200 200])
    title('Phase Corrected Sum of Squares  root')
    
    subplot(414)
    Orthoview(INV1final, pos, [-200 200])
    title('INV1 final')
end


%% Calculate lambda
% Usually the multiplicative factor shouldn't be greater than 10, but that is not the ase when the image is
% bias field corrected, in which case the noise estimated at the edge of the image might not be a good measure

while ~strcmpi(FinalChoice, 'y')
    
    noiselevel = multiplyingFactor*mean(mean(mean(INV2img(1:end, end-10:end, end-10:end))));
    
    % MP2RAGEimgRobustScanner = MP2RAGErobustfunc(INV1img, INV2img, noiselevel.^2);
    MP2RAGEimgRobustPhaseSensitive = MP2RAGErobustfunc(INV1final, INV2img, noiselevel.^2);
    
    if isgraphics(HG)
        
        % Robust Image view
        range = [-0.5 0.40];
        
        subplot(211)
        Orthoview(UNIimg, pos, range)
        title('MP2RAGE UNI-Image')
        
        % subplot(312)
        % Orthoview(MP2RAGEimgRobustScanner, pos, range)
        % title('MP2RAGE Robust Scanner')
        
        subplot(212)
        Orthoview(MP2RAGEimgRobustPhaseSensitive, pos, range)
        title('MP2RAGE Robust')
        ylabel(['Noise level = ' num2str(multiplyingFactor)])
        
        if isempty(regularization)
            FinalChoice = input('Is it a satisfactory noise level?? (y/n) [n]: ', 's');
            if strcmpi(FinalChoice,'y')
                fprintf('Final regularization noise level = %g\n\n', multiplyingFactor)
            else
                multiplyingFactor = input(['New regularization noise level (current = ' num2str(multiplyingFactor) '): ']);
                HG = gcf;                                   % Make sure we have a figure to plot the noise level
            end
        else
            FinalChoice = 'y';
        end
        
    else
        
        FinalChoice = 'y';
        
    end
    
end


%% Save data if filenameOUT is given

if isfield(MP2RAGE, 'filenameOUT') && ~isempty(MP2RAGE.filenameOUT)

    % Save a nifti-file
    disp(['Saving: ' MP2RAGE.filenameOUT])
    if integerformat==0
        UNIimg = MP2RAGEimgRobustPhaseSensitive;
    else
        UNIimg = round(4095 * (MP2RAGEimgRobustPhaseSensitive + 0.5));
    end
    UNIhdr.fname = MP2RAGE.filenameOUT;
    spm_write_vol_gz(UNIhdr, UNIimg);
    
    % Read & enrich the UNI json-file and write it as a T1w json-file
    [UNIpath, UNIname] = myfileparts(MP2RAGE.filenameUNI);
    jsonUNIfile        = fullfile(UNIpath, [UNIname '.json']);
    if isfile(jsonUNIfile)
        jsonT1                     = jsondecode(fileread(jsonUNIfile));
        jsonT1.BasedOn             = {MP2RAGE.filenameUNI, MP2RAGE.filenameINV1, MP2RAGE.filenameINV2};
        jsonT1.SeriesDescription   = [jsonT1.ProtocolName '_MP2RAGE_denoised_background'];
        jsonT1.NoiseRegularization = multiplyingFactor;
        
        [T1path, T1fname] = myfileparts(MP2RAGE.filenameOUT);
        fid = fopen(fullfile(T1path, [T1fname '.json']), 'w');
        fprintf(fid, '%s', jsonencode(jsonT1));
        fclose(fid);
    end
    
end

% Clean-up the temporarily unzipped/smoothed images
for TempVol = [INV1hdr, INV2hdr, UNIhdr]
    if startsWith(TempVol.fname, tempdir) && isfile(TempVol.fname)
        delete(TempVol.fname)
    end
end
