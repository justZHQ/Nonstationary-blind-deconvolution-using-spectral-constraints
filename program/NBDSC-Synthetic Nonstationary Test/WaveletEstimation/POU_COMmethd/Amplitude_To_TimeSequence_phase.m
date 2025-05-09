function [signal_rot]=Amplitude_To_TimeSequence_phase(dt,amplitude_spectrum,fmax,number,length_wavelet,phase)
% dt;%时间采样/s
% signal;%输入时域信号
% fmax;%最大显示频率/Hz
% amplitude_spectrum;%输入信号signal的振幅谱

N=2^ceil(log2(number));
Amplitude=zeros(N/2,1);
fmax_number=ceil(fmax/(1/(N*dt)));
Amplitude(1:fmax_number)=amplitude_spectrum*N/2;
% Amplitude(1:fmax_number)=amplitude_spectrum*N/2.*exp(phase*i);
Amplitude=[Amplitude;flip(Amplitude)];
xifft=ifft(Amplitude);
realx=real(xifft);
n=length(realx);
realx=[realx(n/2:end);realx(1:n/2)];
nn=(n-length_wavelet+1)/2;
realx=realx(nn+2:n-nn+2)';

% realx=realx./max(realx);
hilbert_realx=hilbert(realx);

% phase=0;
% signal_rot=signal*cos(phase)+hilbert_signal*sin(phase);
signal_rot=real(hilbert_realx)*cos(phase)+imag(hilbert_realx)*sin(phase);

% figure;
% plot(real(hilbert_realx),'b');
% hold on;
% plot(signal_rot,'r');


end
