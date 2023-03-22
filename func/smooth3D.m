function S=smooth3D(X,FWHM,vox)

% will apply a 3D gaussian smoothing filter to a 3D image, call with
% S=smooth(X,FWFHM,voxdims)
% FWHM is in mm - (FWHM = 2.35 * stdev)

%zero fill to required just in case it is not even elements
dimsorig=size(X);
dim = round(dimsorig(1:3)./2)*2;

if dimsorig(1:3) ~= dim

    X = zeropad_odd_dimension(X,'pre',dim);

end



dim=size(X);
%create gaussian smoothing matrix

x1=((-dim(1)/2:dim(1)/2-1))*vox(1);
x2=((-dim(2)/2:dim(2)/2-1))*vox(2);
x3=((-dim(3)/2:dim(3)/2-1))*vox(3);


[X1,X2,X3]=ndgrid(x1,x2,x3);

dev=(FWHM/2.35);

GAUSS=exp(-(X1.^2+X2.^2+X3.^2)/(2*dev^2));
GAUSS=GAUSS/(sum(abs(GAUSS(:))));


S=zeros(dim);

S=(ifft3s(bsxfun(@times,fft3s(X),fft3s(GAUSS))));

if dimsorig(1:3) ~= dim

    S = zeropad_odd_dimension(S,'post',dimsorig);

end




function out = fft3s(data)
% function out = fft3s(data)
%
% Do 3D FFT

out = fftshift(fftshift(fftshift(fft(fft(fft(fftshift(fftshift(fftshift(data,1),2),3),[],1),[],2),[],3),1),2),3);
end;

function out = ifft3s(data)
% function out = fft3s(data)
%
% Do 3D iFFT

out = ifftshift(ifftshift(ifftshift(ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(data,1),2),3),[],1),[],2),[],3),1),2),3);
end;

function output = zeropad_odd_dimension(input,mode,matrixSize_o)

matrixSize = size(input);

% determine if a dimension needs to be zeropadded
padsize     = zeros(size(matrixSize));
for kd = 1:length(matrixSize)
    if mod(matrixSize(kd),2) == 1 && kd < 4
        padsize(kd) = 1;
    end
end

switch mode
    case 'pre'
        % zero padding if the dimension of the matrix is an odd number
        output = padarray(input, padsize, 0,'post');

    case 'post'
        % remove zero padding
        output = input(1:matrixSize_o(1),1:matrixSize_o(2),1:matrixSize_o(3),:,:);

end

end
end

