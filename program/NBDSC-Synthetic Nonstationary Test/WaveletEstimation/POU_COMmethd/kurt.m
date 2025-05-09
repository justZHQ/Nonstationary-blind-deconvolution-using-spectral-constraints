function [signal_rot,kurt]=kurt(signal,phase)
% dt;%时间采样/s
% L;%子波长度；%wavelength=2*L+1
% rick_mainf;%雷克子波的主频/HZ
% fa;%相位
% fmax;%最大显示频率/Hz
% t;%输出时间序列/s
% ricker;%输出时间域雷克子波
% f;%输出频率序列
% amplitude_spectrum;%输出雷克子波的振幅谱

% t=-L*dt:dt:L*dt;
% s1=(1-2*(pi*rick_mainf*t).^2).*exp(-(pi*rick_mainf*t).^2);
hilbert_signal=hilbert(signal);
% signal_rot=signal*cos(phase)+hilbert_signal*sin(phase);
signal_rot=real(hilbert_signal)*cos(phase)+imag(hilbert_signal)*sin(phase);
% signal_rot=cos(signal*phase)+sin(hilbert_signal*phase);
kurt=length(signal_rot)*sum(signal_rot.^4)/(sum(signal_rot.^2))^2-3;
% [f,amplitude_spectrum]=Amplitude_spectrum_my(dt,ricker,fmax);
