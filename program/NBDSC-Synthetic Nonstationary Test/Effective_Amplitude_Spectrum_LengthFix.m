function [BegNum,Bandf,df,ff,log_ampx,ampx]=Effective_Amplitude_Spectrum_LengthFix(dt,DATA,Begf,Endf,number)
% dt;%时间采样/s
% signal;%输入时域信号
% fmax;%最大显示频率/Hz
% f;%输出频率序列/Hz
% amplitude_spectrum;%输入信号signal的振幅谱
% F=fft(signal,2^ceil(log2(number)));
% N=length(F);
% amplitude_spectrum=abs(F);%abs(Cn)=abs(X(k))=sqrt(realX(k)^2+imagX(k)^2)
% f=(0:N-1)*(1/(N*dt));%频率采样间隔等于基波频率f0=1/T;T=dt*N=N/fs
% 
% 
% fmax_number=ceil(fmax/(1/(N*dt)));
% f=f(1:fmax_number);
% amplitude_spectrum=amplitude_spectrum(1:fmax_number)*2/N;%频率与真实振幅之间的关系——振幅谱；
% F=F(1:fmax_number);
% end

[~,dm]=size(DATA);
F=abs(fft(DATA(:,1),2^ceil(log2(number))));
N=length(F);
f=(0:N-1)*(1/(N*dt));%频率采样间隔等于基波频率f0=1/T;T=dt*N=N/fs
df=f(2)-f(1);
% frequency band parameter slection
BegNum=floor(Begf/df)+1; 
EndNum=floor(Endf/df)+1;
Bandf=((BegNum-1):(EndNum-1))*df;
% Bandf=Bandf';
log_ampx=[];
ampx=[];
for i=1:dm
AbsSrcf=abs(fft(DATA(:,i),2^ceil(log2(number))));
AbsSrcf=AbsSrcf;%鍔熺巼璋?
BandAbsSrcf=AbsSrcf(BegNum:EndNum);
log_ampx=[log_ampx,log(BandAbsSrcf)];
ampx=[ampx,BandAbsSrcf];
end

ff=f(BegNum:EndNum);