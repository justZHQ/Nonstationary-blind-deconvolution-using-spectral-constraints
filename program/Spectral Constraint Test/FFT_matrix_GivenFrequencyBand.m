function [f,F,F_cut]=FFT_matrix_GivenFrequencyBand(dt,Wavelength,FFT_length,FrequencyBand)
% dt;%时间采样/s
% signal;%输入时域信号
% fmax;%最大显示频率/Hz
% f;%输出频率序列/Hz
% amplitude_spectrum;%输入信号signal的振幅谱
% F=fft(signal,2^ceil(log2(1000)));


N=FFT_length;
F=dftmtx(N);
df=1/(N*dt);
f=(0:N-1)*df;
fmax_number=ceil(FrequencyBand./df);
F=F(fmax_number(1):fmax_number(2),:);
f=f(fmax_number(1):fmax_number(2));
F_cut=F(:,1:Wavelength);
