function [INV1final, INV2final] = Correct_INV1INV2_withMP2RAGEuni(INV1img,INV2img,MP2RAGEimg, CorrectBoth )
% usage:
% function [INV1final, INV2final] = Correct_INV1INV2_withMP2RAGEuni(INV1img,INV2img,MP2RAGEimg )
% all inputs and outputs are structures where the image informtion is
% stored on the .img field


% just ensuring they are doubles for the time being
INV1img.img = double(INV1img.img);
INV2img.img = double(INV2img.img);
MP2RAGEimg.img = double(MP2RAGEimg.img);

if min(MP2RAGEimg.img(:))>=0 && max(MP2RAGEimg.img(:))>=0.51
    % converts MP2RAGE to -0.5clear to 0.5 scale - assumes that it is getting only positive values
    %     MP2RAGEimg.img = (MP2RAGEimg.img - max(MP2RAGEimg.img(:))/2) ./ max(MP2RAGEimg.img(:));
    MP2RAGEimg.img = (MP2RAGEimg.img - 4095/2) ./ max(4095);
    integerformat = 1;
else
    integerformat = 0;
end

if isempty(CorrectBoth)
    CorrectBoth = 0;
end;
%% Compute correct INV1 & INV2 datasets by using the phase sensitivity information available on thee UNI image

MP2RAGErobustfunc = @(INV1, INV2, beta)(conj(INV1).*INV2-beta) ./ (INV1.^2 + INV2.^2 + 2*beta);
rootsquares_pos   = @(a, b, c)(-b + sqrt(b.^2 - 4*a.*c)) ./ (2*a);
rootsquares_neg   = @(a, b, c)(-b - sqrt(b.^2 - 4*a.*c)) ./ (2*a);


% Gives the correct polarity to INV1
INV1img.img = sign(MP2RAGEimg.img) .* INV1img.img;

% Because the MP2RAGE INV1 and INV2 is a summ of squares data, while the
% MP2RAGEimg is a phase sensitive coil combination.. some more maths has to
% be performed to get a better INV2 and INV1 estimate, which here is done by assuming
% both the bigger value at any pixel will be the closest approximations to a phase sensitive combination

INV1pos = rootsquares_pos(-MP2RAGEimg.img, INV2img.img, -INV2img.img.^2 .* MP2RAGEimg.img);
INV1neg = rootsquares_neg(-MP2RAGEimg.img, INV2img.img, -INV2img.img.^2 .* MP2RAGEimg.img);

INV2pos = rootsquares_pos(-MP2RAGEimg.img, INV1img.img, -INV1img.img.^2 .* MP2RAGEimg.img);
INV2neg = rootsquares_neg(-MP2RAGEimg.img, INV1img.img, -INV1img.img.^2 .* MP2RAGEimg.img);

if CorrectBoth==1
    %making the correction when INV2img>abs(INV1img)
    INV1final = INV1img;
    INV1final.img(abs(INV1img.img-INV1pos) >  abs(INV1img.img-INV1neg)) = INV1neg(abs(INV1img.img-INV1pos) >  abs(INV1img.img-INV1neg));
    INV1final.img(abs(INV1img.img-INV1pos) <= abs(INV1img.img-INV1neg)) = INV1pos(abs(INV1img.img-INV1pos) <= abs(INV1img.img-INV1neg));
    %     INV1final.img(abs(INV1img.img) > abs(INV2img.img) ) = INV1img.img(abs(INV1img.img) > abs(INV2img.img) );



    %making the correction when abs(INV1img)>INV2img
    INV2final = INV2img;
    INV2final.img(abs(INV2img.img-INV2pos) >  abs(INV2img.img-INV2neg)) = INV2neg(abs(INV2img.img-INV2pos) >  abs(INV2img.img-INV2neg));
    INV2final.img(abs(INV2img.img-INV2pos) <= abs(INV2img.img-INV2neg)) = INV2pos(abs(INV2img.img-INV2pos) <= abs(INV2img.img-INV2neg));
    %     INV2final.img(abs(INV2img.img) >= abs(INV1img.img) ) = INV2img.img(abs(INV2img.img) >= abs(INV1img.img) );

    % REINFORCE THE POLARITY
    INV2final.img = abs(INV2final.img);
    INV1final.img= sign(MP2RAGEimg.img) .* abs(INV1final.img);


elseif CorrectBoth== 0

    INV2final = INV2img;

    INV1final = INV1img;
    INV1final.img(abs(INV1img.img-INV1pos) >  abs(INV1img.img-INV1neg)) = INV1neg(abs(INV1img.img-INV1pos) >  abs(INV1img.img-INV1neg));
    INV1final.img(abs(INV1img.img-INV1pos) <= abs(INV1img.img-INV1neg)) = INV1pos(abs(INV1img.img-INV1pos) <= abs(INV1img.img-INV1neg));
    %     INV1final.img(abs(INV1img.img) > abs(0.6 * INV2img.img) ) = INV1img.img(abs(INV1img.img) > abs(0.6*INV2img.img) );

else
    % currently the  only condition is that the difference is smaller, but no
    % requirement is made on the polarity being respected:
    % equal to the positive if closer and polarity respected
    % equal to negative if closer and polarity is respected
    INV2final = INV2img;
    INV1final = INV1img;
    INV1final.img(and(abs(INV1img.img-INV1pos) >  abs(INV1img.img-INV1neg) , INV1neg.*INV1img.img > 0 )) = INV1neg(and(abs(INV1img.img-INV1pos) >  abs(INV1img.img-INV1neg), INV1neg.*INV1img.img > 0 ));
    INV1final.img(and(abs(INV1img.img-INV1pos) <= abs(INV1img.img-INV1neg), INV1neg.*INV1img.img > 0 )) = INV1pos(and(abs(INV1img.img-INV1pos) <= abs(INV1img.img-INV1neg), INV1neg.*INV1img.img > 0 ));
    %     INV1final.img(abs(INV1img.img) > abs(0.6 * INV2img.img) ) = INV1img.img(abs(INV1img.img) > abs(0.6*INV2img.img) );

end
% only applies the correction if one of the images has less than 60 % of
% the SNR of the higher SNR image
INV1final.img(abs(INV1img.img) > abs(0.6 * INV2img.img) ) = INV1img.img(abs(INV1img.img) > abs(0.6 * INV2img.img) );
INV2final.img(abs(INV2img.img) >= abs(0.6 * INV1img.img) ) = INV2img.img(abs(INV2img.img) >= abs(0.6 * INV1img.img) );
