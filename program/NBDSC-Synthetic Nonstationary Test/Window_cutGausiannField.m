function [signal_afterwin,number,Gaussian,nstart,nend]=Window_cutGausiannField(GausiannVariance,Winsize,step,signal)

% Winsize=1000;
% average=150;
% startwin=100;
% signal=randn(1000,1);
% step=100;
% GausiannVariance=20;
% 
% 
% Winsize=length(signal);
T=-Winsize/2+1:Winsize/2;
% T=T';

Gaussian=exp(-T'.^2/(2*GausiannVariance^2));
% figure;
% plot(Gaussian);
% qwavelet=[];
% qwaveletwin=[];
signal_afterwin=[];
% Gaussian_win=[];
number=floor((length(signal(:))-Winsize)/step)+1;
nstart=[];
nend=[];
for i=1:number
%     averagei=startwin+(i-1)*step;
%     Ti=T-averagei;
%     Gaussian=exp(-Ti.^2/(2*GausiannVariance^2));
%     signal_afterwin=[signal_afterwin,signal.*Gaussian];
%     Gaussian_win=[Gaussian_win,Gaussian];
starti=(i-1)*step+1;
endi=(i-1)*step+Winsize;
mid=round((starti+endi)/2);
% qwavelet=[qwavelet,qmat(starti:endi,mid)];
% qwaveletwin=[qwaveletwin,Gaussian.*qmat(starti:endi,mid)];
signal_afterwin=[signal_afterwin,Gaussian.*signal(starti:endi)];
nstart=[nstart;starti];
nend=[nend;endi];

end