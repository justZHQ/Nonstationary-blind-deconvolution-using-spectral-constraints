function [signal_rot,phase_rot]=maxkurt(signal,phase)
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
n=length(phase);
kurtvec=zeros(n,1);
for i=1:n
% hilbert_signal=hilbert(signal);
% signal_rot=signal*cos(fa)+hilbert_signal*sin(phase(i));
% kurt=length(signal)*sum(signal.^4)/(sum(signal.^2))^2-3;
[~,kurtvec(i)]=kurt(signal,phase(i));
end
[~,mi]=max(kurtvec);
phase_rot=phase(mi);
[signal_rot,~]=kurt(signal,phase_rot);