function [t,klauder,f,amplitude_spectrum]=Klauder_my(dt,L,F_L,F_H,T,fa,fmax)
% dt;%时间采样/s
% L;%子波长度；%wavelength=2*L+1
% F_L,F_H;%klauder子波的低截频和高截频/HZ，满足F_L<F_H
% T;%周期/s
% fa;%相位
% fmax;%最大显示频率/Hz
% t;%输出时间序列/s
% klauder;%输出时间域klauder子波
% f;%输出频率序列
% amplitude_spectrum;%输出klauder子波的振幅谱
t=-L*dt:dt:L*dt;
t((length(t)+1)/2)=0.000000001;
k=(F_H-F_L)/T;
f0=(F_H+F_L)/2;
s1=sin(pi*k.*t.*(T-t))./(pi*k.*t).*cos(2*pi*f0.*t);
s1=s1/max(s1);
s2=hilbert(s1);
klauder=real(s2)*cos(fa)+imag(s2)*sin(fa);
[f,amplitude_spectrum]=Amplitude_spectrum_my(dt,klauder,fmax);
