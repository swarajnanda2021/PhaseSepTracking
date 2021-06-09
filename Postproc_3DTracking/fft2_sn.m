function [Ftransform] = fft2_sn(Data,index)

if index==1
    Dy = Data' - mean(Data'); % Average over length for each timestep
else
    Dy = (Data - mean(Data))'; %Average over time for each lengthstep
end



NFFT_timey  = 2 ^ (nextpow2( size(Dy,1) )); % Calculate fft data output size as a power of 2 (for fft efficiency)
NFFT_spacey = 2 ^ (nextpow2( size(Dy,2) ));
ry        = size(Dy,1); % ry and cy are window dimensions
cy        = size(Dy,2);
wy        = window2(ry,cy,@hann); %@rectwin,%@hamming, @hann % A windowing function
pow_windy = (sum(sum (sqrt(wy.^2))))/(ry*cy);
% complete 2-D FFT of Dy
Ftransform_jb = fft2(Dy.*wy,NFFT_timey,NFFT_spacey);
% fft2 = Y2Dy_jb;
Ftransform_jb = Ftransform_jb(:,:)/(0.5*pow_windy*NFFT_timey*NFFT_spacey*mean(mean(Data)));
Ftransform = fftshift(Ftransform_jb); %shift origin of frequency to center



end

