function [realx]=Amplitude_To_TimeSequence(dt,amplitude_spectrum,fmax,number,length_wavelet)
% dt;%时间采样/s
% signal;%输入时域信号
% fmax;%最大显示频率/Hz
% amplitude_spectrum;%输入信号signal的振幅谱

N=2^ceil(log2(number));
Amplitude=zeros(N/2,1);
fmax_number=ceil(fmax/(1/(N*dt)));
Amplitude(1:fmax_number)=amplitude_spectrum*N/2;
Amplitude=[Amplitude;flip(Amplitude)];
xifft=ifft(Amplitude);
realx=real(xifft);
n=length(realx);
realx=[realx(n/2:end);realx(1:n/2)];
nn=(n-length_wavelet+1)/2;
realx=realx(nn+2:n-nn+2)';
% realx=realx./max(realx);

end
