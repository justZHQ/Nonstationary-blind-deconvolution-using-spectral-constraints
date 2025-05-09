clear;clc;close all;
addpath('SpectralConstraints');
addpath(genpath(fullfile(pwd, 'WaveletEstimation')));
addpath(genpath(fullfile(pwd, 'AlternateIteration')));
%--------------------------generate A and R-------------------%
n=2000;m=1;dt=0.001;
wave_L=51;
fmax=200;fa=0;
t=0:dt:(n-1)*dt;

%% give wave and generate A
[t_wavelet,wave,f_wavelet,amplitude_spectrum_wavelet]=Ricker_my(dt,wave_L,40,fa,fmax);
length_wavelet=length(wave);

%% generate R 
load R.mat
gcf1=figure;
set(gcf1,'position',[800 600 800 200]);
plot(t,R,'k','linewidth',1.2);
ylabel('Amplitude');
xlabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
set(gca,'position',[0.12 0.26 0.86 0.6]);
% xlim([0 n*dt]);
annotation('textbox',[.04 .76 .1 .2], ...
    'String','a)','EdgeColor','none','FontSize',14,'FontWeight','bold');
%% generate Y 
Qeff=ones(size(R))*150;
S=zeros(size(R));
qmat=qmatrix(Qeff,t,wave,t_wavelet);%build the qmatrix
S=qmat*R;
SNR=10;SNoise=awgn(S,SNR,'measured');
SNRDB=snr(S,SNoise-S)
load SNoise.mat

gcf2=figure;
set(gcf2,'position',[800 600 800 200]);
plot(t,S,'k','linewidth',1.2);
ylabel('Amplitude');
xlabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
set(gca,'position',[0.12 0.26 0.86 0.6]);
xlim([0 n*dt]);
% text(-0.15,0.85,'b)','FontSize',14,'FontWeight','bold');
annotation('textbox',[.04 .76 .1 .2], ...
    'String','b)','EdgeColor','none','FontSize',14,'FontWeight','bold');
%% POU Processing
GausiannVariance=10000000000;
step=250;L=step;
Winsize=250;

sigma = 0.01;
numBasis=80;
molecular=[10,20,30,40,50,60,70];

[molecular_vecter,molecular_Signal]=POU_my(numBasis,molecular,sigma,dt,S);
gcf3=figure;
set(gcf3,'position',[800 600 800 200]);
plot(t,molecular_vecter');
ylabel('Amplitude');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
set(gca,'position',[0.12 0.26 0.86 0.6]);
ylim([0 1.5]);
txt1='${\psi _1}(t)$';txt2='${\psi _2}(t)$';txt3='${\psi _3}(t)$';txt4='${\psi _4}(t)$';txt5='${\psi _5}(t)$';txt6='${\psi _6}(t)$';txt7='${\psi _7}(t)$';txt8='${\psi _8}(t)$';
text(0.07,1.2,txt1,'FontSize',14,'interpreter','latex');text(0.32,1.2,txt2,'FontSize',14,'interpreter','latex');text(0.57,1.2,txt3,'FontSize',14,'interpreter','latex');
text(0.82,1.2,txt4,'FontSize',14,'interpreter','latex');text(1.07,1.2,txt5,'FontSize',14,'interpreter','latex');text(1.32,1.2,txt6,'FontSize',14,'interpreter','latex');
text(1.57,1.2,txt7,'FontSize',14,'interpreter','latex');text(1.82,1.2,txt8,'FontSize',14,'interpreter','latex');


[qwavelet,qwaveletwin,signal_afterwin,N,Gausiann,nstart,nend]=Window_cutGausiann(GausiannVariance,Winsize,step,S,qmat);
[~,~,signalNoise_afterwin,~,~,~,~]=Window_cutGausiann(GausiannVariance,Winsize,step,SNoise,qmat);
[~,~,R_afterwin,~,~,~,~]=Window_cutGausiann(GausiannVariance,Winsize,step,R,qmat);

Rwindows=R_afterwin;
Wwindows=qwaveletwin;
Swindows=signal_afterwin;
SwindowsNoise=signalNoise_afterwin;

%% Amplitude spectrum of each POU seismic trace
number=1024;
Begfmin=0;
Endfmax=fmax;
[~,~,~,~,W_log_ampx,W_ampx]=Effective_Amplitude_Spectrum_LengthFix(dt,Wwindows,Begfmin,Endfmax,number);
[BegNum,Bandf,df,ff,log_ampx,ampx]=Effective_Amplitude_Spectrum_LengthFix(dt,Swindows,Begfmin,Endfmax,number);
[~,~,~,~,~,ampxNoise]=Effective_Amplitude_Spectrum_LengthFix(dt,SwindowsNoise,Begfmin,Endfmax,number);
fnumber=length(ff);

%% Amplitude spectrum estimation by AWPSI-COM method
Begf_high=5;
Endf_low=50;
Begf_high=10;
Endf_low=60;
epsilon=0.01;
[Begf,Endf]=Set_wi(df,epsilon,W_ampx);
% Begf=zeros(N,1);
Begf=0*ones(N,1);Begf=Begf*df;Endf=Endf*df;
Iterations=10;

p=0.26*ones(N,1);
% p=0.2*ones(N,1);
GausiannVariance_AWPSI=4;

lamda1=0;
lamda2=0;
lamda1AWPSI=6000;
lamda2AWPSI=6000;

[log_ampx_Correction,ampx_Correction_COMAWPSI]=Fast_COM(lamda1AWPSI,lamda2AWPSI,ampx,df,p,Iterations,Begf,Endf,Begf_high,Endf_low,GausiannVariance_AWPSI);
[~,ampxNoise_Correction_COMAWPSI]=Fast_COM(lamda1AWPSI,lamda2AWPSI,ampxNoise,df,p,Iterations,Begf,Endf,Begf_high,Endf_low,GausiannVariance_AWPSI);
[~,ampx_Correction_COM]=Fast_COM(lamda1,lamda2,ampx,df,p,Iterations,Begf,Endf,Begf_high,Endf_low,GausiannVariance_AWPSI);
[~,ampxNoise_Correction_COM]=Fast_COM(lamda1,lamda2,ampxNoise,df,p,Iterations,Begf,Endf,Begf_high,Endf_low,GausiannVariance_AWPSI);
ampx_Correction_COMAWPSI=ampx_Correction_COMAWPSI.^2;
ampxNoise_Correction_COMAWPSI=ampxNoise_Correction_COMAWPSI.^2;
ampx_Correction_COM=ampx_Correction_COM.^2;
ampxNoise_Correction_COM=ampxNoise_Correction_COM.^2;
W_ampx=W_ampx.^2;
ampx=ampx.^2;

W_ampxPlot=zeros(fnumber,2*N+1);
ampxPlot=zeros(fnumber,2*N+1);
ampxNoisePlot=zeros(fnumber,2*N+1);
ampx_Correction_COMAWPSIPlot=zeros(fnumber,2*N+1);
ampxNoise_Correction_COMAWPSIPlot=zeros(fnumber,2*N+1);
ampx_Correction_COMPlot=zeros(fnumber,2*N+1);
ampxNoise_Correction_COMPlot=zeros(fnumber,2*N+1);

W_ampxPlot(:,2:2:2*N+1)=W_ampx./max(W_ampx);
ampxPlot(:,2:2:2*N+1)=ampx./max(ampx);
ampxNoisePlot(:,2:2:2*N+1)=ampxNoise./max(ampxNoise);
ampx_Correction_COMAWPSIPlot(:,2:2:2*N+1)=ampx_Correction_COMAWPSI./max(ampx_Correction_COMAWPSI);
ampxNoise_Correction_COMAWPSIPlot(:,2:2:2*N+1)=ampxNoise_Correction_COMAWPSI./max(ampxNoise_Correction_COMAWPSI);
ampx_Correction_COMPlot(:,2:2:2*N+1)=ampx_Correction_COM./max(ampx_Correction_COM);
ampxNoise_Correction_COMPlot(:,2:2:2*N+1)=ampxNoise_Correction_COM./max(ampxNoise_Correction_COM);

zwin=0:1:2*N;z=zwin;
%%
gcf4=figure;
set(gcf4,'position',[800 400 1000 400]);
subplot(8,5,1);
wigb_my(-ampxPlot,1,zwin,ff);title({'Seismic trace';'power spectrum'});
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
xticklabels({'0','100'});
yticks(zwin);
yticklabels({'','no.1','','no.2','','no.3','','no.4','','no.5','','no.6','','no.7','','no.8',''});
set(gca,'position',[0.07 0.12 0.14 0.77]);
annotation('textbox',[.025 .78 .1 .2], ...
    'String','a)','EdgeColor','none','FontSize',14,'FontWeight','bold');
xlim([-30 160]);
for i=1:N-1
    BegNuma=floor(Begf(i)/df)+1; 
    BegNumb=floor(Begf(i+1)/df)+1; 
    EndNuma=floor(Endf(i)/df)+1;
    EndNumb=floor(Endf(i+1)/df)+1;
    
    Begf_highNuma=floor(Begf_high/df)+1;
    Begf_highNumb=floor(Begf_high/df)+1;
    Endf_lowNuma=floor(Endf_low/df)+1;
    Endf_lowNumb=floor(Endf_low/df)+1;
    
    hold on;
% 获取a和b点的时间值
time_a = z(2*i);
time_b = z(2*i+2);
% 获取a和b点的频率值
Begfreq_a = ff(BegNuma);
Begfreq_b = ff(BegNumb);
Endfreq_a = ff(EndNuma);
Endfreq_b = ff(EndNumb);

Begf_highfreq_a = ff(Begf_highNuma);
Begf_highfreq_b = ff(Begf_highNumb);
Endf_lowfreq_a = ff(Endf_lowNuma);
Endf_lowfreq_b = ff(Endf_lowNumb);
% 绘制红色线连接a和b点
plot([Begfreq_a, Begfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);
plot([Endfreq_a, Endfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);

plot([Endf_lowfreq_a, Endf_lowfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);
plot([Begf_highfreq_a, Begf_highfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);

txt1=strcat('$$f_',num2str(i),'^1$$');
txt2=strcat('$$f_',num2str(i),'^2$$');
text(Begfreq_a-22,time_a +.3,txt1,'Interpreter', 'latex','FontSize',12,'Color','red');
text(Endfreq_a+1,time_a +.3,txt2,'Interpreter', 'latex','FontSize',12,'Color','red');
if i==N-1
txt1=strcat('$$f_{',num2str(i+1),'}^1$$');
txt2=strcat('$$f_{',num2str(i+1),'}^2$$');

text(Begfreq_b-22,time_b +.3,txt1,'Interpreter', 'latex','FontSize',12,'Color','red');
text(Endfreq_b+1,time_b +.3,txt2,'Interpreter', 'latex','FontSize',12,'Color','red');

text(Begf_highfreq_b+1,time_b +.8,'$$f^1$$','Interpreter', 'latex','FontSize',12,'Color','green');
text(Endf_lowfreq_b+1,time_b +.8,'$$f^2$$','Interpreter', 'latex','FontSize',12,'Color','green');
end
end
subplot(8,5,2);
wigb_my(-W_ampxPlot,1,zwin,ff);title({'Ture wavelet';'power spectrum'});
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
xticklabels({'0','100'});
yticks(zwin);
yticklabels({''});
set(gca,'position',[0.22 0.12 0.14 0.77]);
xlim([-30 160]);
for i=1:N-1
[maxa,peak_a]=max(W_ampxPlot(:,2*i));
[maxab,peak_b]=max(W_ampxPlot(:,2*i+2));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
time_aa =  2*i;
time_bb = 2*i+2;
plot([peakfreq_a, peakfreq_b],[time_aa-W_ampxPlot(peak_a,time_aa)./max(W_ampxPlot(:,time_aa))-2, time_bb-W_ampxPlot(peak_b,time_bb)./max(W_ampxPlot(:,time_bb))-2], '-*y', 'LineWidth', 2);
end
subplot(8,5,3);
wigb_my(-ampx_Correction_COMPlot,1,zwin,ff);title('Result of COM');
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
xticklabels({'0','100'});
yticks(zwin);
yticklabels({''});
set(gca,'position',[0.37 0.12 0.14 0.77]);
xlim([-30 160]);
for i=1:N-1
[maxa,peak_a]=max(ampx_Correction_COMPlot(:,2*i));
[maxab,peak_b]=max(ampx_Correction_COMPlot(:,2*i+2));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
time_aa =  2*i;
time_bb = 2*i+2;
plot([peakfreq_a, peakfreq_b],[time_aa-ampx_Correction_COMPlot(peak_a,time_aa)./max(ampx_Correction_COMPlot(:,time_aa))-2, time_bb-ampx_Correction_COMPlot(peak_b,time_bb)./max(ampx_Correction_COMPlot(:,time_bb))-2], '-*y', 'LineWidth', 2);
end
subplot(8,5,4);
wigb_my(-ampx_Correction_COMAWPSIPlot,1,zwin,ff);title('Result of AWPSI-COM');
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
xticklabels({'0','100'});
yticks(zwin);
yticklabels({''});
set(gca,'position',[0.52 0.12 0.14 0.77]);
xlim([-30 160]);
for i=1:N-1
[maxa,peak_a]=max(ampx_Correction_COMAWPSIPlot(:,2*i));
[maxab,peak_b]=max(ampx_Correction_COMAWPSIPlot(:,2*i+2));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
time_aa =  2*i;
time_bb = 2*i+2;
plot([peakfreq_a, peakfreq_b],[time_aa-ampx_Correction_COMAWPSIPlot(peak_a,time_aa)./max(ampx_Correction_COMAWPSIPlot(:,time_aa))-2, time_bb-ampx_Correction_COMAWPSIPlot(peak_b,time_bb)./max(ampx_Correction_COMAWPSIPlot(:,time_bb))-2], '-*y', 'LineWidth', 2);
end
for i=1:N
    subplot(N,5,i*5);hold on;box on;
    plot(ff,W_ampxPlot(:,2*i)./max(W_ampxPlot(:,2*i)),'r','linewidth',1.2);
    plot(ff,ampx_Correction_COMAWPSIPlot(:,2*i)./max(ampx_Correction_COMAWPSIPlot(:,2*i)),'b','linewidth',1.2);
    plot(ff,ampx_Correction_COMPlot(:,2*i)./max(ampx_Correction_COMPlot(:,2*i)),'-.k','linewidth',1.2);
%     ylabel('no.1');
    set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5);
    set(gca,'TickLength',[0.008 0.001]);
    yticklabels({''});
    if i==1
        title('Power spectrum');
        legend('True wavelet', 'AWPSI-COM','COM', 'Location', 'none', 'Position', [0.72,0.95,0.1,0.07],'Orientation','horizontal');
        legend('boxoff');
    end
     if i~=N
    xticklabels({''});
%     xlim([15 80]);
     else
         xlabel('Frequency/Hz');
%          xlim([15 80]);
%          xticks([20 70]);
%           xlim([15 80]);
     end
    set(gca,'position',[0.7 0.92-i*0.1 0.2 0.08]);
    text(-35,0.7,['no.',num2str(i)],'FontSize',12);
    text(100,0.75,['COM:',num2str(round(corr(W_ampxPlot(:,2*i)./max(W_ampxPlot(:,2*i)),ampx_Correction_COMPlot(:,2*i)./max(ampx_Correction_COMPlot(:,2*i))),3))],'FontSize',8);
    text(100,0.35,['AWPSI-COM:',num2str(round(corr(W_ampxPlot(:,2*i)./max(W_ampxPlot(:,2*i)),ampx_Correction_COMAWPSIPlot(:,2*i)./max(ampx_Correction_COMAWPSIPlot(:,2*i))),3))],'FontSize',8);
end

%%
gcf5=figure;
set(gcf5,'position',[800 400 1000 400]);
subplot(8,5,1);
wigb_my(-ampxNoisePlot,1,zwin,ff);title({'Seismic trace';'power spectrum'});
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
xticklabels({'0','100'});
yticks(zwin);
yticklabels({'','no.1','','no.2','','no.3','','no.4','','no.5','','no.6','','no.7','','no.8',''});
set(gca,'position',[0.07 0.12 0.14 0.77]);
annotation('textbox',[.025 .78 .1 .2], ...
    'String','b)','EdgeColor','none','FontSize',14,'FontWeight','bold');
xlim([-30 160]);
for i=1:N-1
    BegNuma=floor(Begf(i)/df)+1; 
    BegNumb=floor(Begf(i+1)/df)+1; 
    EndNuma=floor(Endf(i)/df)+1;
    EndNumb=floor(Endf(i+1)/df)+1;
    
    Begf_highNuma=floor(Begf_high/df)+1;
    Begf_highNumb=floor(Begf_high/df)+1;
    Endf_lowNuma=floor(Endf_low/df)+1;
    Endf_lowNumb=floor(Endf_low/df)+1;
    
    hold on;
% 获取a和b点的时间值
time_a = z(2*i);
time_b = z(2*i+2);
% 获取a和b点的频率值
Begfreq_a = ff(BegNuma);
Begfreq_b = ff(BegNumb);
Endfreq_a = ff(EndNuma);
Endfreq_b = ff(EndNumb);

Begf_highfreq_a = ff(Begf_highNuma);
Begf_highfreq_b = ff(Begf_highNumb);
Endf_lowfreq_a = ff(Endf_lowNuma);
Endf_lowfreq_b = ff(Endf_lowNumb);
% 绘制红色线连接a和b点
plot([Begfreq_a, Begfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);
plot([Endfreq_a, Endfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);

plot([Endf_lowfreq_a, Endf_lowfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);
plot([Begf_highfreq_a, Begf_highfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);

txt1=strcat('$$f_',num2str(i),'^1$$');
txt2=strcat('$$f_',num2str(i),'^2$$');
text(Begfreq_a-22,time_a +.3,txt1,'Interpreter', 'latex','FontSize',12,'Color','red');
text(Endfreq_a+1,time_a +.3,txt2,'Interpreter', 'latex','FontSize',12,'Color','red');
if i==N-1
txt1=strcat('$$f_{',num2str(i+1),'}^1$$');
txt2=strcat('$$f_{',num2str(i+1),'}^2$$');

text(Begfreq_b-22,time_b +.3,txt1,'Interpreter', 'latex','FontSize',12,'Color','red');
text(Endfreq_b+1,time_b +.3,txt2,'Interpreter', 'latex','FontSize',12,'Color','red');

text(Begf_highfreq_b+1,time_b +.8,'$$f^1$$','Interpreter', 'latex','FontSize',12,'Color','green');
text(Endf_lowfreq_b+1,time_b +.8,'$$f^2$$','Interpreter', 'latex','FontSize',12,'Color','green');
end
end
subplot(8,5,2);
wigb_my(-W_ampxPlot,1,zwin,ff);title({'Ture wavelet';'power spectrum'});
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
xticklabels({'0','100'});
yticks(zwin);
yticklabels({''});
set(gca,'position',[0.22 0.12 0.14 0.77]);
xlim([-30 160]);
for i=1:N-1
[maxa,peak_a]=max(W_ampxPlot(:,2*i));
[maxab,peak_b]=max(W_ampxPlot(:,2*i+2));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
time_aa =  2*i;
time_bb = 2*i+2;
plot([peakfreq_a, peakfreq_b],[time_aa-W_ampxPlot(peak_a,time_aa)./max(W_ampxPlot(:,time_aa))-2, time_bb-W_ampxPlot(peak_b,time_bb)./max(W_ampxPlot(:,time_bb))-2], '-*y', 'LineWidth', 2);
end

subplot(8,5,3);
wigb_my(-ampxNoise_Correction_COMPlot,1,zwin,ff);title('Result of COM');
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
xticklabels({'0','100'});
yticks(zwin);
yticklabels({''});
set(gca,'position',[0.37 0.12 0.14 0.77]);
xlim([-30 160]);
for i=1:N-1
[maxa,peak_a]=max(ampxNoise_Correction_COMPlot(:,2*i));
[maxab,peak_b]=max(ampxNoise_Correction_COMPlot(:,2*i+2));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
time_aa =  2*i;
time_bb = 2*i+2;
plot([peakfreq_a, peakfreq_b],[time_aa-ampxNoise_Correction_COMPlot(peak_a,time_aa)./max(ampxNoise_Correction_COMPlot(:,time_aa))-2, time_bb-ampxNoise_Correction_COMPlot(peak_b,time_bb)./max(ampxNoise_Correction_COMPlot(:,time_bb))-2], '-*y', 'LineWidth', 2);
end

subplot(8,5,4);
wigb_my(-ampxNoise_Correction_COMAWPSIPlot,1,zwin,ff);title('Result of AWPSI-COM');
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
xticklabels({'0','100'});
yticks(zwin);
yticklabels({''});
set(gca,'position',[0.52 0.12 0.14 0.77]);
xlim([-30 160]);
for i=1:N-1
[maxa,peak_a]=max(ampxNoise_Correction_COMAWPSIPlot(:,2*i));
[maxab,peak_b]=max(ampxNoise_Correction_COMAWPSIPlot(:,2*i+2));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
time_aa =  2*i;
time_bb = 2*i+2;
plot([peakfreq_a, peakfreq_b],[time_aa-ampxNoise_Correction_COMAWPSIPlot(peak_a,time_aa)./max(ampxNoise_Correction_COMAWPSIPlot(:,time_aa))-2, time_bb-ampxNoise_Correction_COMAWPSIPlot(peak_b,time_bb)./max(ampxNoise_Correction_COMAWPSIPlot(:,time_bb))-2], '-*y', 'LineWidth', 2);
end
for i=1:N
    subplot(N,5,i*5);hold on;box on;
    plot(ff,W_ampxPlot(:,2*i)./max(W_ampxPlot(:,2*i)),'r','linewidth',1.2);
    plot(ff,ampxNoise_Correction_COMAWPSIPlot(:,2*i)./max(ampxNoise_Correction_COMAWPSIPlot(:,2*i)),'b','linewidth',1.2);
    plot(ff,ampxNoise_Correction_COMPlot(:,2*i)./max(ampxNoise_Correction_COMPlot(:,2*i)),'-.k','linewidth',1.2);
%     ylabel('no.1');
    set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5);
    set(gca,'TickLength',[0.008 0.001]);
    yticklabels({''});
    if i==1
        title('Power spectrum');
        legend('True wavelet', 'AWPSI-COM','COM', 'Location', 'none', 'Position', [0.72,0.95,0.1,0.07],'Orientation','horizontal');
        legend('boxoff');
    end
     if i~=N
    xticklabels({''});
%     xlim([15 80]);
     else
         xlabel('Frequency/Hz');
%          xlim([15 80]);
%          xticks([20 70]);
%           xlim([15 80]);
     end
    set(gca,'position',[0.7 0.92-i*0.1 0.2 0.08]);
%     text(120,0.7,['no.',num2str(i)],'FontSize',12);
    text(-35,0.7,['no.',num2str(i)],'FontSize',12);
    text(100,0.75,['COM:',num2str(round(corr(W_ampxPlot(:,2*i)./max(W_ampxPlot(:,2*i)),ampx_Correction_COMPlot(:,2*i)./max(ampx_Correction_COMPlot(:,2*i))),3))],'FontSize',8);
    text(100,0.35,['AWPSI-COM:',num2str(round(corr(W_ampxPlot(:,2*i)./max(W_ampxPlot(:,2*i)),ampx_Correction_COMAWPSIPlot(:,2*i)./max(ampx_Correction_COMAWPSIPlot(:,2*i))),3))],'FontSize',8);
end

phase_rot_fit=zeros(N,1);
%% Time domain wavelet obtained by kurtosis maximizing constant phase rotation based on AWPSI
FFT_Length=number;
WaveLength_Give=length_wavelet;
EstimationWave=[];
EstimationWaveNoise=[];
% figure;
for i=1:N
    [EstimationWavei]=Amplitude_To_TimeSequence_phase(dt,sqrt(ampx_Correction_COMAWPSI(:,i)),fmax,FFT_Length,WaveLength_Give,phase_rot_fit(i));
    [EstimationWaveiNoise]=Amplitude_To_TimeSequence_phase(dt,sqrt(ampxNoise_Correction_COMAWPSI(:,i)),fmax,FFT_Length,WaveLength_Give,phase_rot_fit(i));
    EstimationWave=[EstimationWave,EstimationWavei'];
    EstimationWaveNoise=[EstimationWaveNoise,EstimationWaveiNoise'];
end

% wave_L=31;length_wavelet=wave_L*2+1;
qwaveletShift=zeros(length_wavelet,N);
EstimationWaveShift=zeros(length_wavelet,N);
EstimationWaveNoiseShift=zeros(length_wavelet,N);

for i=1:N
    maxfqwavelet=find(qwavelet(:,i)==max(qwavelet(:,i)),1);
    maxfEstimationWave=find(EstimationWave(:,i)==max(EstimationWave(:,i)),1); 
    maxfEstimationWaveNoise=find(EstimationWaveNoise(:,i)==max(EstimationWaveNoise(:,i)),1); 
    Difference=0;
    qwaveletShift(:,i)=qwavelet(maxfqwavelet-wave_L-Difference:maxfqwavelet+wave_L-Difference,i);
    EstimationWaveShift(:,i)=EstimationWave(maxfEstimationWave-wave_L:maxfEstimationWave+wave_L,i);
    EstimationWaveNoiseShift(:,i)=EstimationWaveNoise(maxfEstimationWaveNoise-wave_L:maxfEstimationWaveNoise+wave_L,i);
end

RwindowsPlot=zeros(Winsize,2*N+1);
SwindowsPlot=zeros(Winsize,2*N+1);
SwindowsNoisePlot=zeros(Winsize,2*N+1);
WwindowsPlot=zeros(length_wavelet,2*N+1);
RwindowsPlot(:,2:2:2*N+1)=Rwindows;
SwindowsPlot(:,2:2:2*N+1)=Swindows;
SwindowsNoisePlot(:,2:2:2*N+1)=SwindowsNoise;
WwindowsPlot(:,2:2:2*N+1)=qwaveletShift;
twindows=0:dt:(Winsize-1)*dt;
scal=0.6;
zwin=0:1:2*N;
gcf6=figure;
set(gcf6,'position',[800 400 800 500]);
subplot(1,3,1);
wigb_Nofill(-WwindowsPlot,scal,zwin,t_wavelet);title('True wavelet');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticks(zwin);
yticklabels({'','no.1','','no.2','','no.3','','no.4','','no.5','','no.6','','no.7','','no.8',''});
set(gca,'position',[0.1 0.1 0.2 0.82]);
subplot(1,3,2);
wigb_Nofill(-SwindowsPlot,scal,zwin,twindows);title('Seismic trace');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticklabels({''});
set(gca,'position',[0.32 0.1 0.3 0.82]);
subplot(1,3,3);
wigb_Nofill(-SwindowsNoisePlot,scal,zwin,twindows);title('Noisy seismic trace');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticklabels({''});
set(gca,'position',[0.64 0.1 0.3 0.82]);
%%
Swindows=Swindows./max(Swindows,[],1);
EstimationWaveShift=EstimationWaveShift./max(EstimationWaveShift,[],1);
SwindowsNoise=SwindowsNoise./max(SwindowsNoise,[],1);
EstimationWaveNoiseShift=EstimationWaveNoiseShift./max(EstimationWaveNoiseShift,[],1);
[Power]=PowerGet(Rwindows);
%%
Wavelength=length(wave);
FFT_length=length(Swindows(:,1));
FrequencyBand=[20,120];
FrequencyBand=[5,120];
[f,F,F_cut]=FFT_matrix_GivenFrequencyBand(dt,Wavelength,FFT_length,FrequencyBand);
F_wave=F_cut*EstimationWaveShift;
p=abs(F_wave);

FrequencyBandNoise=[20,150];
FrequencyBandNoise=[5,200];
[fNoise,~,F_cutNoise]=FFT_matrix_GivenFrequencyBand(dt,Wavelength,FFT_length,FrequencyBandNoise);
F_waveNoise=F_cutNoise*EstimationWaveNoiseShift;
pNoise=abs(F_waveNoise);

% kSC=norm(abs(F_cut*wave').^2-p.^2,2)/norm(wave,2);
L0=1.05;eta=2;
iter_alternate=30;iter_fista=100;opt=1;
lambda=0.8;% 1 norm of r 0.8
lambda0=0.0015;% SpectralConstraints of wavelet
lambdaNoise=1;% 1 norm of r 0.8
lambda0Noise=0.003;% SpectralConstraints of wavelet
lambda1=0.2;% 1 norm of wavelet noiseadd 0.1
lambda2SC=0;% 2 norm of wavelet noiseadd 0.1
kSC=norm(abs(F_cut*EstimationWaveShift).^2,2)/norm(EstimationWaveShift,2)
lambda2=lambda0*kSC*100;% 2 norm of wavelet noiseadd 0.1
kSCNoise=norm(abs(F_cutNoise*EstimationWaveNoiseShift).^2,2)/norm(EstimationWaveNoiseShift,2)
lambda2Noise=lambda0Noise*kSCNoise*100;% 2 norm of wavelet noiseadd 0.1

[R_hat,InverseWave,Convergence1]=Alternate_iteration_SpectralConstraints_AWPSI(F_cut,p,0,lambda1,lambda2,Swindows,EstimationWaveShift,L0,eta,lambda,iter_alternate,iter_fista,opt);
[R_hatPc,InverseWavePc,Convergence2]=Alternate_iteration_SpectralConstraints_AWPSI(F_cut,p,lambda0,lambda1,lambda2SC,Swindows,EstimationWaveShift,L0,eta,lambda,iter_alternate,iter_fista,opt);
[R_hatNoise,InverseWaveNoise,~]=Alternate_iteration_SpectralConstraints_AWPSI(F_cutNoise,pNoise,0,lambda1,lambda2Noise,SwindowsNoise,EstimationWaveNoiseShift,L0,eta,lambda,iter_alternate,iter_fista,opt);
[R_hatPcNoise,InverseWavePcNoise,~]=Alternate_iteration_SpectralConstraints_AWPSI(F_cutNoise,pNoise,lambda0Noise,lambda1,lambda2SC,SwindowsNoise,EstimationWaveNoiseShift,L0,eta,lambdaNoise,iter_alternate,iter_fista,opt);

R_hat=R_hat.*Power;
R_hatPc=R_hatPc.*Power;
R_hatNoise=R_hatNoise.*Power;
R_hatPcNoise=R_hatPcNoise.*Power;
%%
TrueWavelet=qwaveletShift./max(qwaveletShift,[],1);
InitialWavelet=EstimationWaveShift./max(EstimationWaveShift,[],1);
InitialWaveletNoise=EstimationWaveNoiseShift./max(EstimationWaveNoiseShift,[],1);

TrueWaveletPlot=zeros(length_wavelet,2*N+1);
RwindowsPlot=zeros(Winsize,2*N+1);
TrueWaveletPlot(:,2:2:2*N+1)=TrueWavelet;
RwindowsPlot(:,2:2:2*N+1)=Rwindows;

InitialWaveletPlot=zeros(length_wavelet,2*N+1);
InverseWavePlot=zeros(length_wavelet,2*N+1);
InverseWavePcPlot=zeros(length_wavelet,2*N+1);
R_hatPlot=zeros(Winsize,2*N+1);
R_hatPcPlot=zeros(Winsize,2*N+1);
InitialWaveletPlot(:,2:2:2*N+1)=InitialWavelet;
InverseWavePlot(:,2:2:2*N+1)=InverseWave;
InverseWavePcPlot(:,2:2:2*N+1)=InverseWavePc;
R_hatPlot(:,2:2:2*N+1)=R_hat;
R_hatPcPlot(:,2:2:2*N+1)=R_hatPc;

InitialWaveletNoisePlot=zeros(length_wavelet,2*N+1);
InverseWaveNoisePlot=zeros(length_wavelet,2*N+1);
InverseWavePcNoisePlot=zeros(length_wavelet,2*N+1);
R_hatNoisePlot=zeros(Winsize,2*N+1);
R_hatPcNoisePlot=zeros(Winsize,2*N+1);
InitialWaveletNoisePlot(:,2:2:2*N+1)=InitialWaveletNoise;
InverseWaveNoisePlot(:,2:2:2*N+1)=InverseWaveNoise;
InverseWavePcNoisePlot(:,2:2:2*N+1)=InverseWavePcNoise;
R_hatNoisePlot(:,2:2:2*N+1)=R_hatNoise;
R_hatPcNoisePlot(:,2:2:2*N+1)=R_hatPcNoise;

ampInverse=abs(F_cut*InverseWave).^2;
ampInversePc=abs(F_cut*InverseWavePc).^2;
ampInverseNoise=abs(F_cutNoise*InverseWaveNoise).^2;
ampInversePcNoise=abs(F_cutNoise*InverseWavePcNoise).^2;
ampTrue=abs(F_cut*TrueWavelet).^2;

gcf7=figure;
set(gcf7,'position',[800 400 800 500]);
subplot(1,2,1);
wigb_Nofill(-TrueWaveletPlot,scal,zwin,t_wavelet);title('True wavelet');
% title('TureWave');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticks(zwin);
yticklabels({'','no.1','','no.2','','no.3','','no.4','','no.5','','no.6','','no.7','','no.8',''});
set(gca,'position',[0.1 0.1 0.2 0.82]);
subplot(1,2,2);
wigb_Nofill(-RwindowsPlot,scal,zwin,twindows);title('True reflectivity');
% title('Rwindows');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
% yticks(zwin);
yticklabels({''});
set(gca,'position',[0.32 0.1 0.3 0.82]);

%%
gcf8=figure;
set(gcf8,'position',[800 400 1000 500]);
subplot(8,4,1);
wigb_Nofill(-RwindowsPlot,scal,zwin,twindows);title('True reflectivity');
% title('TureWave');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticks(zwin);
yticklabels({'','no.1','','no.2','','no.3','','no.4','','no.5','','no.6','','no.7','','no.8',''});
set(gca,'position',[0.05 0.1 0.2 0.82]);
subplot(8,4,2);
wigb_Nofill(-R_hatPlot,scal,zwin,twindows);title('Estimated reflectivity');
% title('Rwindows');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
% yticks(zwin);
yticklabels({''});
set(gca,'position',[0.265 0.1 0.2 0.82]);
subplot(8,4,3);
wigb_Nofill(-InverseWavePlot,scal,zwin,t_wavelet);title('Estimated wavelet');
% title('TureWave');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticks(zwin);
yticklabels({''});
% yticklabels({'','no.1','','no.2','','no.3','','no.4','','no.5','','no.6','','no.7','','no.8',''});
set(gca,'position',[0.48 0.1 0.2 0.82]);
for i=1:N
annotation('textbox',[.62 .732-0.082*i .1 .2], ...
    'String',num2str(round(corr(TrueWavelet(:,i),InverseWave(:,i)),3)),'EdgeColor','none','FontSize',12,'FontWeight','bold','Color','red');
end
% subplot(1,2,2);
for i=1:N
    subplot(N,4,i*4);hold on;box on;
%     plot(f,p(:,i).^2./max(p(:,i).^2),'r','linewidth',1.2);
    plot(f,ampTrue(:,i)./max(ampTrue(:,i)),'r','linewidth',1.2);
%     plot(f,ampInverse(:,i)./max(ampInverse(:,i)),'k','linewidth',1.2);
    plot(f,ampInverse(:,i)./max(ampInverse(:,i)),'b','linewidth',1.2);
%     ylabel('no.1');
    set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5);
    set(gca,'TickLength',[0.008 0.001]);
    yticklabels({''});
    if i==1
        title('Power spectrum');
        legend('True wavelet', 'Estimated wavelet', 'Location', 'none', 'Position', [0.72,0.95,0.1,0.07],'Orientation','horizontal');
        legend('boxoff');
    end
     if i~=N
    xticklabels({''});
    xlim([15 100]);
     else
         xlabel('Frequency/Hz');
%          xlim([15 80]);
         xticks([20 70]);
          xlim([15 100]);
     end
    set(gca,'position',[0.7 0.938-i*0.105 0.2 0.09]);
    text(65,0.7,['no.',num2str(i)],'FontSize',12);
end
annotation('textbox',[.02 .78 .1 .2], ...
    'String','a)','EdgeColor','none','FontSize',14,'FontWeight','bold');

gcf9=figure;
set(gcf9,'position',[800 400 1000 500]);
subplot(8,4,1);
wigb_Nofill(-RwindowsPlot,scal,zwin,twindows);title('True reflectivity');
% title('TureWave');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticks(zwin);
yticklabels({'','no.1','','no.2','','no.3','','no.4','','no.5','','no.6','','no.7','','no.8',''});
set(gca,'position',[0.05 0.1 0.2 0.82]);
subplot(8,4,2);
wigb_Nofill(-R_hatPcPlot,scal,zwin,twindows);title('Estimated reflectivity');
% title('Rwindows');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
% yticks(zwin);
yticklabels({''});
set(gca,'position',[0.265 0.1 0.2 0.82]);
subplot(8,4,3);
wigb_Nofill(-InverseWavePcPlot,scal,zwin,t_wavelet);title('Estimated wavelet');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
% yticks(zwin);
yticklabels({''});
set(gca,'position',[0.48 0.1 0.2 0.82]);
for i=1:N
annotation('textbox',[.62 .732-0.082*i .1 .2], ...
    'String',num2str(round(corr(TrueWavelet(:,i),InverseWavePc(:,i)),3)),'EdgeColor','none','FontSize',12,'FontWeight','bold','Color','red');
end
for i=1:N
    subplot(N,4,i*4);hold on;box on;
%     plot(f,p(:,i).^2./max(p(:,i).^2),'r','linewidth',1.2);
    plot(f,ampTrue(:,i)./max(ampTrue(:,i)),'r','linewidth',1.2);
%     plot(f,ampInverse(:,i)./max(ampInverse(:,i)),'k','linewidth',1.2);
    plot(f,ampInversePc(:,i)./max(ampInversePc(:,i)),'b','linewidth',1.2);
%     ylabel('no.1');
    set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5);
    set(gca,'TickLength',[0.008 0.001]);
    yticklabels({''});
    if i==1
        title('Power spectrum');
        legend('True wavelet', 'Estimated wavelet', 'Location', 'none', 'Position', [0.72,0.95,0.1,0.07],'Orientation','horizontal');
        legend('boxoff');
    end
     if i~=N
    xticklabels({''});
    xlim([15 100]);
     else
         xlabel('Frequency/Hz');
%          xlim([15 80]);
         xticks([20 70]);
          xlim([15 100]);
     end
    set(gca,'position',[0.7 0.938-i*0.105 0.2 0.09]);
    text(65,0.7,['no.',num2str(i)],'FontSize',12);
end
annotation('textbox',[.02 .78 .1 .2], ...
    'String','b)','EdgeColor','none','FontSize',14,'FontWeight','bold');
%%
twin=.2*1;tinc=.01*3;%defines the Gaussian windows
tinc=.1*2;
gdb=60;%truncation factor in fgabor
p=1; %Leave this at 1
% normflag=0;
tsmo=0.1;%temporal smoother in seconds
fsmo=10;%frequency smoother in Hz
ihyp=0;%flag for hyperbolic smoothing. 2 gets AWPSI, 1 gets Hyperbolic, 0 gets boxcar
order=10;%order of the Burg spectrum if iburg is 1
phase=0;%gabor decon stab factor and phase flag
ipow=1;%1 means the output trace will be balanced in power to the input, 0 means not balancing.
transforms=1;
taperpct=60;
stabg1=0.00001;
stabg2=0.0001;
stabg3=0.001;
[r21,tvs_op1]=gabordeconb(S,t,twin,tinc,tsmo,fsmo,ihyp,order,stabg1,phase,p,gdb,taperpct);
[r22,tvs_op2]=gabordeconb(S,t,twin,tinc,tsmo,fsmo,ihyp,order,stabg2,phase,p,gdb,taperpct);
[r23,tvs_op3]=gabordeconb(S,t,twin,tinc,tsmo,fsmo,ihyp,order,stabg3,phase,p,gdb,taperpct);
%%
R_result=[];
R_result=[S,R_result,reshape(Rwindows,[],1),reshape(R_hat,[],1),reshape(R_hatPc,[],1),r21,r22,r23];
R_resultPlot=zeros(n,13);
% R_resultPlot(:,1:2:13)=R_result;
R_resultPlot(:,1:2:13)=R_result./[norm(R_result(:,1),2),norm(R_result(:,2),2),norm(R_result(:,3),2),norm(R_result(:,4),2),norm(R_result(:,5),2),norm(R_result(:,6),2),norm(R_result(:,7),2)];
yno=1:1:13;
gcf10=figure;
set(gcf10,'position',[800 300 800 300]);
wigb_Nofill(-R_resultPlot,0.6,yno,t);
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticks(yno(1:end));
yticklabels({'Synthetic','','True','','NBDM','','NBDSC','','GD-stab=0.00001','','GD-stab=0.0001','','GD-stab=0.001'});
set(gca,'position',[0.12 0.2 0.86 0.76]);
annotation('textbox',[.07 .81 .1 .2], ...
    'String','a)','EdgeColor','none','FontSize',14,'FontWeight','bold'); 
%%
gcf11=figure;
set(gcf11,'position',[800 400 1000 500]);
subplot(8,4,1);
wigb_Nofill(-RwindowsPlot,scal,zwin,twindows);title('True reflectivity');
% title('TureWave');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticks(zwin);
yticklabels({'','no.1','','no.2','','no.3','','no.4','','no.5','','no.6','','no.7','','no.8',''});
set(gca,'position',[0.05 0.1 0.2 0.82]);
subplot(8,4,2);
wigb_Nofill(-R_hatNoisePlot,scal,zwin,twindows);title('Estimated reflectivity');
% title('Rwindows');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
% yticks(zwin);
yticklabels({''});
set(gca,'position',[0.265 0.1 0.2 0.82]);
subplot(8,4,3);
wigb_Nofill(-InverseWaveNoisePlot,scal,zwin,t_wavelet);title('Estimated wavelet');
% title('TureWave');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticks(zwin);
yticklabels({''});
% yticklabels({'','no.1','','no.2','','no.3','','no.4','','no.5','','no.6','','no.7','','no.8',''});
set(gca,'position',[0.48 0.1 0.2 0.82]);
for i=1:N
annotation('textbox',[.62 .732-0.082*i .1 .2], ...
    'String',num2str(round(corr(TrueWavelet(:,i),InverseWaveNoise(:,i)),3)),'EdgeColor','none','FontSize',12,'FontWeight','bold','Color','red');
end
for i=1:N
    subplot(N,4,i*4);hold on;box on;
%     plot(fNoise,pNoise(:,i).^2./max(pNoise(:,i).^2),'r','linewidth',1.2);
    plot(f,ampTrue(:,i)./max(ampTrue(:,i)),'r','linewidth',1.2);
%     plot(f,ampInverse(:,i)./max(ampInverse(:,i)),'k','linewidth',1.2);
    plot(fNoise,ampInverseNoise(:,i)./max(ampInverseNoise(:,i)),'b','linewidth',1.2);
%     ylabel('no.1');
    set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5);
    set(gca,'TickLength',[0.008 0.001]);
    yticklabels({''});
    if i==1
        title('Power spectrum');
        legend('True wavelet', 'Estimated wavelet', 'Location', 'none', 'Position', [0.72,0.95,0.1,0.07],'Orientation','horizontal');
        legend('boxoff');
    end
     if i~=N
    xticklabels({''});
%     xlim([15 80]);
     else
         xlabel('Frequency/Hz');
%          xlim([15 80]);
%          xticks([20 70]);
%           xlim([15 80]);
     end
%     set(gca,'position',[0.64 0.938-i*0.105 0.3 0.09]);
    set(gca,'position',[0.7 0.938-i*0.105 0.2 0.09]);
    text(120,0.7,['no.',num2str(i)],'FontSize',12);
end
annotation('textbox',[.02 .78 .1 .2], ...
    'String','a)','EdgeColor','none','FontSize',14,'FontWeight','bold');

gcf12=figure;
set(gcf12,'position',[800 400 1000 500]);
subplot(8,4,1);
wigb_Nofill(-RwindowsPlot,scal,zwin,twindows);title('True reflectivity');
% title('TureWave');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticks(zwin);
yticklabels({'','no.1','','no.2','','no.3','','no.4','','no.5','','no.6','','no.7','','no.8',''});
set(gca,'position',[0.05 0.1 0.2 0.82]);
subplot(8,4,2);
wigb_Nofill(-R_hatPcNoisePlot,scal,zwin,twindows);title('Estimated reflectivity');
% title('Rwindows');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
% yticks(zwin);
yticklabels({''});
set(gca,'position',[0.265 0.1 0.2 0.82]);
subplot(8,4,3);
wigb_Nofill(-InverseWavePcNoisePlot,scal,zwin,t_wavelet);title('Estimated wavelet');
% title('TureWave');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticks(zwin);
yticklabels({''});
% yticklabels({'','no.1','','no.2','','no.3','','no.4','','no.5','','no.6','','no.7','','no.8',''});
set(gca,'position',[0.48 0.1 0.2 0.82]);
for i=1:N
annotation('textbox',[.62 .732-0.082*i .1 .2], ...
    'String',num2str(round(corr(TrueWavelet(:,i),InverseWavePcNoise(:,i)),3)),'EdgeColor','none','FontSize',12,'FontWeight','bold','Color','red');
end
for i=1:N
    subplot(N,4,i*4);hold on;box on;
%     plot(fNoise,pNoise(:,i).^2./max(pNoise(:,i).^2),'r','linewidth',1.2);
    plot(f,ampTrue(:,i)./max(ampTrue(:,i)),'r','linewidth',1.2);
%     plot(f,ampInverse(:,i)./max(ampInverse(:,i)),'k','linewidth',1.2);
    plot(fNoise,ampInversePcNoise(:,i)./max(ampInversePcNoise(:,i)),'b','linewidth',1.2);
%     ylabel('no.1');
    set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5);
    set(gca,'TickLength',[0.008 0.001]);
    yticklabels({''});
    if i==1
        title('Power spectrum');
        legend('True wavelet', 'Estimated wavelet', 'Location', 'none', 'Position', [0.72,0.95,0.1,0.07],'Orientation','horizontal');
        legend('boxoff');
    end
     if i~=N
    xticklabels({''});
%     xlim([15 80]);
     else
         xlabel('Frequency/Hz');
%          xlim([15 80]);
%          xticks([20 70]);
%           xlim([15 80]);
     end
    set(gca,'position',[0.7 0.938-i*0.105 0.2 0.09]);
    text(120,0.7,['no.',num2str(i)],'FontSize',12);
end
annotation('textbox',[.02 .78 .1 .2], ...
    'String','b)','EdgeColor','none','FontSize',14,'FontWeight','bold');
%%
stabgNoise1=0.1;
stabgNoise2=1;
stabgNoise3=10;
[r2Noise1,tvs_opNoise1]=gabordeconb(SNoise,t,twin,tinc,tsmo,fsmo,ihyp,order,stabgNoise1,phase,p,gdb,taperpct);
[r2Noise2,tvs_opNoise2]=gabordeconb(SNoise,t,twin,tinc,tsmo,fsmo,ihyp,order,stabgNoise2,phase,p,gdb,taperpct);
[r2Noise3,tvs_opNoise3]=gabordeconb(SNoise,t,twin,tinc,tsmo,fsmo,ihyp,order,stabgNoise3,phase,p,gdb,taperpct);
%%
R_result=[];
R_result=[SNoise,R_result,reshape(Rwindows,[],1),reshape(R_hatNoise,[],1),reshape(R_hatPcNoise,[],1),r2Noise1,r2Noise2,r2Noise3];
R_resultPlot=zeros(n,13);
R_resultPlot(:,1:2:13)=R_result./[norm(R_result(:,1),2),norm(R_result(:,2),2),norm(R_result(:,3),2),norm(R_result(:,4),2),norm(R_result(:,5),2),norm(R_result(:,6),2),norm(R_result(:,7),2)];
yno=1:1:13;
gcf13=figure;hold on;
set(gcf13,'position',[800 300 800 300]);
wigb_Nofill(-R_resultPlot,0.6,yno,t);
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0.008 0.001]);
yticks(yno(1:end));
yticklabels({'Synthetic','','True','','NBDM','','NBDSC','','GD-stab=0.1','','GD-stab=1','','GD-stab=10'});
set(gca,'position',[0.12 0.2 0.86 0.76]);
annotation('textbox',[.07 .81 .1 .2], ...
    'String','b)','EdgeColor','none','FontSize',14,'FontWeight','bold'); 