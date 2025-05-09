clear all;
close all;
clc;
delt=0;Q=300;
dt=0.001;%鏃堕棿閲囨牱闂撮殧/s
L=35;
fmax=200;%鏈?澶ф樉绀洪鐜?/Hz
fa=0;%鐩镐綅
seismic_length=800;
ASSW_Ormsby=0;
ASSW_Klauder=0;
for i_plot=1:500
if ASSW_Ormsby<0.008&&ASSW_Klauder<0.0185&&i_plot>1
    break;
else
close all;
[t_wavelet,wavelet,f_wavelet,amplitude_spectrum_wavelet]=Ricker_my(dt,L,40,fa,fmax);
waveletname='Ricker';
r_Gaussian=normrnd(0,.3,1,seismic_length);
t_seismic=0:dt:(seismic_length-1)*dt;
r=r_Gaussian.^3/max(abs(r_Gaussian))/1;
Rname='姝ｆ?佸垎甯?';
r1=r;
save('..\reflection coefficients and wavelet\r1.mat','r1');
%%
seismic_original=conv(r,wavelet);
seismic=seismic_original(L+1:end-L);
Conv_Length=length(seismic_original);
[f_seismic,F_seismic,amplitude_spectrum_seismic,N_length]=Amplitude_spectrum_my_LengthFix(dt,seismic_original,fmax,Conv_Length);
[f_wavelet,F_wavelet,amplitude_spectrum_wavelet,~]=Amplitude_spectrum_my_LengthFix(dt,wavelet,fmax,Conv_Length);
mu=0;
amplitude_spectrum_wavelet=amplitude_spectrum_wavelet';
%%
% amplitude_spectrum_wavelet=amplitude_spectrum_wavelet.*exp(-f_wavelet'/800);
amplitude_spectrum_wavelet=amplitude_spectrum_wavelet.*exp(-f_wavelet'*pi*delt/Q);
[wavelet]=Amplitude_To_TimeSequence(dt,amplitude_spectrum_wavelet,fmax,Conv_Length,length(wavelet));

seismic_original=conv(r,wavelet);
seismic=seismic_original(L+1:end-L);
Conv_Length=length(seismic_original);
[f_seismic,F_seismic,amplitude_spectrum_seismic,N_length]=Amplitude_spectrum_my_LengthFix(dt,seismic_original,fmax,Conv_Length);
[f_wavelet,F_wavelet,amplitude_spectrum_wavelet,~]=Amplitude_spectrum_my_LengthFix(dt,wavelet,fmax,Conv_Length);

mu=0;
amplitude_spectrum_wavelet=amplitude_spectrum_wavelet';
amplitude_spectrum_seismic=amplitude_spectrum_seismic';
m=inv(amplitude_spectrum_wavelet'*amplitude_spectrum_wavelet+mu)*amplitude_spectrum_wavelet'*amplitude_spectrum_seismic;
amplitude_spectrum_seismic=amplitude_spectrum_seismic/m;
%% COM parameters and estimated ASSW    
mu=0; %正则化参数 
Iterations=50; %迭代次数
p=0.4;
Contraction_start=5;
Contraction_end=120;
Contraction=zeros(size(amplitude_spectrum_seismic));
Contraction(Contraction_start:Contraction_end)=Contraction_mapping(f_seismic(Contraction_start:Contraction_end),amplitude_spectrum_seismic(Contraction_start:Contraction_end),p,Iterations,mu);
%%
Signal_Wavelet=amplitude_spectrum_wavelet;
Signal_Seismic=amplitude_spectrum_seismic;
%% POU decomposition
dt_POU=2;
GaussianLength=100;
GausiannVariance=40;
GivenSignalLength=300;
Scale=40;
TimeOffset=150;
EffectiveAtoms_Start=2;
EffectiveAtoms_End=5;
number=EffectiveAtoms_End-EffectiveAtoms_Start+1;
[NormalizedGaborAtomicLibrary,TimeShiftedSignal_Wavelet,GaborTimeShiftedSignalTruncation_Wavelet]=Gabor(dt_POU,GaussianLength,GausiannVariance,GivenSignalLength,Scale,TimeOffset,Signal_Wavelet);
[GaborIntransformation_Wavelet]=InGabor(EffectiveAtoms_Start,EffectiveAtoms_End,GaborTimeShiftedSignalTruncation_Wavelet);
[~,TimeShiftedSignal_Seismic,GaborTimeShiftedSignalTruncation_Seismic]=Gabor(dt_POU,GaussianLength,GausiannVariance,GivenSignalLength,Scale,TimeOffset,Signal_Seismic);
[GaborIntransformation_Seismic]=InGabor(EffectiveAtoms_Start,EffectiveAtoms_End,GaborTimeShiftedSignalTruncation_Seismic);
Wavelet=[];
Seismic=[];
for i=EffectiveAtoms_Start:EffectiveAtoms_End
    Wavelet=[Wavelet,GaborTimeShiftedSignalTruncation_Wavelet(:,i)];
    Seismic=[Seismic,GaborTimeShiftedSignalTruncation_Seismic(:,i)];
end

%% POU-COM parameters and estimated ASSW   
df=f_seismic(2)-f_seismic(1);
Iterations=20;
mu=0;
VecterCut_start=[5;10;20;45];
VecterCut_end=[50;75;100;120];
Vecterp=[0.3;0.3;0.3;0.3];
flip=[1,1,1,1];
[Vectercp,Vecteralpha,Vecterbeta,k]=Parameter_determination_ALL(EffectiveAtoms_Start,EffectiveAtoms_End,df,Seismic,Vecterp,VecterCut_start,VecterCut_end,mu,flip);
data=Seismic;
for i=1:Iterations
    [InversionResult,AfterFitting]=Gaborfitting_ALL(data,df,Vecterp,VecterCut_start,VecterCut_end,Vectercp,Vecteralpha,Vecterbeta,k,EffectiveAtoms_Start,EffectiveAtoms_End,flip);
    data=AfterFitting;
end
RestoringSignal=InversionResult;
%% zero-phase seismic wavelet (ZSW)
[wavelet]=Amplitude_To_TimeSequence(dt,Signal_Wavelet,fmax,Conv_Length,length(wavelet));
[realx_pou]=Amplitude_To_TimeSequence(dt,RestoringSignal,fmax,Conv_Length,length(wavelet));
[realx_Contraction]=Amplitude_To_TimeSequence(dt,Contraction,fmax,Conv_Length,length(wavelet));
%% misfit function
ASSW_COM=Contraction./norm(Contraction,2)-Signal_Wavelet./norm(Signal_Wavelet,2);
ASSW_POU=RestoringSignal./norm(RestoringSignal,2)-Signal_Wavelet./norm(Signal_Wavelet,2);
ASSW_COM_time=realx_Contraction./norm(realx_Contraction,2)-wavelet./norm(wavelet,2);
ASSW_POU_time=realx_pou./norm(realx_pou,2)-wavelet./norm(wavelet,2);
%% 
F1=figure;
set(F1,'position',[20 10 750 980]);
subplot(7,4,[1,2]);
hold on;box on;
plot(t_seismic,r,'b','linewidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(a)'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.855 0.39 0.135]);
subplot(7,4,[3,4]);
hold on;box on;
plot(t_seismic,seismic,'b','LineWidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(b)'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.58 0.855 0.39 0.135]);
subplot(7,4,[5,6,9,10,13,14,17,18]);
hold on;box on;
plot(f_wavelet,amplitude_spectrum_wavelet,'r','LineWidth',1.5);
plot(f_seismic,amplitude_spectrum_seismic,'b','LineWidth',1);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(c)'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.515 0.39 0.26]);

for i=1:number
    xx=0:0:0;
    subplot(7,4,4*(i-1)+7);
    hold on;box on;
    Synthetic_seismic_record=plot(Seismic(:,i),'b','LineWidth',1);
    Seismic_wavelet=plot(Wavelet(:,i),'r','LineWidth',1.5);
    plot(AfterFitting(:,i)*k(i),'--k','LineWidth',1.5);
    ylabel('Magnitude');
    if i==number
    xlabel('Frequency/Hz');
    else
    set(gca,'xtick',xx);
    end
    set(gca,'FontName','Arial','FontSize',9.5,'linewidth',1.5,'FontWeight','bold');
    set(gca,'TickLength',[0 0.0001]); 
    set(gca,'position',[0.575 0.725-(i-1)*0.072 0.17 0.05]);
    subplot(7,4,4*(i-1)+8);
    hold on;box on;
    plot(Wavelet(:,i),'r','LineWidth',2);plot(AfterFitting(:,i)*k(i),'--k','LineWidth',2);
    if i==number
    xlabel('Frequency/Hz');
    else
    set(gca,'xtick',xx);
    end    
    set(gca,'FontName','Arial','FontSize',9.5,'linewidth',1.5,'FontWeight','bold');
    set(gca,'TickLength',[0 0.0001]); 
    set(gca,'position',[0.805 0.725-(i-1)*0.072 0.17 0.05]);
end

subplot(7,4,[21,22]);
hold on;box on;
plot(f_wavelet,Signal_Wavelet,'r','LineWidth',1.5);
plot(f_seismic,Contraction,'--g','LineWidth',1.5);
plot(f_seismic,RestoringSignal,'--k','LineWidth',1.5);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(e)'});
legend('Ture','COM','POU-COM','location','NorthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.298 0.39 0.135]);
subplot(7,4,[23,24]);
hold on;box on;
plot(t_wavelet,wavelet,'r','LineWidth',1.5);
plot(t_wavelet,realx_Contraction,'--g','LineWidth',1.5);
plot(t_wavelet,realx_pou,'--k','LineWidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(f)'});
legend('Ture','COM','POU-COM','location','NorthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.585 0.298 0.39 0.135]);
subplot(7,4,[25,26]);
hold on;box on;
plot(f_seismic,ASSW_COM,'b','linewidth',1.5);
plot(f_seismic,ASSW_POU,'--r','linewidth',1.5);
plot(f_seismic,zeros(size(ASSW_COM)),'k','linewidth',1.5);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(g)'});
legend('COM','POU-COM','location','SouthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.08 0.39 0.135]);
subplot(7,4,[27,28]);
hold on;box on;
plot(t_wavelet,ASSW_COM_time,'b','linewidth',1.5);
plot(t_wavelet,ASSW_POU_time,'--r','linewidth',1.5);
plot(t_wavelet,zeros(size(ASSW_POU_time)),'k','linewidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(h)'});
legend('COM','POU-COM','location','SouthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.585 0.08 0.39 0.135]);
annotation(F1,'textbox',...
    [0.749928571428572 0.44183673527168 0.0657142841815949 0.0295918361568938],...
    'String','(d)',...
    'LineStyle','none',...
    'Interpreter','none',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off');
Ricker1=wavelet;

%Ricker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ormsby

L=35;
[t_wavelet,wavelet,f_wavelet,amplitude_spectrum_wavelet]=Ormsby_my(dt,L,30,60,110,140,fa,fmax);
waveletname='Ormsby';
%%
seismic_original=conv(r,wavelet);
seismic=seismic_original(L+1:end-L);
Conv_Length=length(seismic_original);
[f_seismic,F_seismic,amplitude_spectrum_seismic,N_length]=Amplitude_spectrum_my_LengthFix(dt,seismic_original,fmax,Conv_Length);
[f_wavelet,F_wavelet,amplitude_spectrum_wavelet,~]=Amplitude_spectrum_my_LengthFix(dt,wavelet,fmax,Conv_Length);
mu=0;
amplitude_spectrum_wavelet=amplitude_spectrum_wavelet';
%%
% amplitude_spectrum_wavelet=amplitude_spectrum_wavelet.*exp(-f_wavelet'/800);
amplitude_spectrum_wavelet=amplitude_spectrum_wavelet.*exp(-f_wavelet'*pi*delt/Q);
[wavelet]=Amplitude_To_TimeSequence(dt,amplitude_spectrum_wavelet,fmax,Conv_Length,length(wavelet));

seismic_original=conv(r,wavelet);
seismic=seismic_original(L+1:end-L);
Conv_Length=length(seismic_original);
[f_seismic,F_seismic,amplitude_spectrum_seismic,N_length]=Amplitude_spectrum_my_LengthFix(dt,seismic_original,fmax,Conv_Length);
[f_wavelet,F_wavelet,amplitude_spectrum_wavelet,~]=Amplitude_spectrum_my_LengthFix(dt,wavelet,fmax,Conv_Length);

mu=0;
amplitude_spectrum_wavelet=amplitude_spectrum_wavelet';
amplitude_spectrum_seismic=amplitude_spectrum_seismic';
m=inv(amplitude_spectrum_wavelet'*amplitude_spectrum_wavelet+mu)*amplitude_spectrum_wavelet'*amplitude_spectrum_seismic;
amplitude_spectrum_seismic=amplitude_spectrum_seismic/m;
%% COM parameters and estimated ASSW    
mu=0; %正则化参数 
Iterations=5; %迭代次数
p=0.4;
Contraction_start=30;
Contraction_end=155;
Contraction=zeros(size(amplitude_spectrum_seismic));
Contraction(Contraction_start:Contraction_end)=Contraction_mapping(f_seismic(Contraction_start:Contraction_end),amplitude_spectrum_seismic(Contraction_start:Contraction_end),p,Iterations,mu);
%%
Signal_Wavelet=amplitude_spectrum_wavelet;
Signal_Seismic=amplitude_spectrum_seismic;
%% POU decomposition
dt_POU=2;
GaussianLength=100;
GausiannVariance=40;
GivenSignalLength=300;
Scale=40;
TimeOffset=150;
EffectiveAtoms_Start=3;
EffectiveAtoms_End=6;
number=EffectiveAtoms_End-EffectiveAtoms_Start+1;
[NormalizedGaborAtomicLibrary,TimeShiftedSignal_Wavelet,GaborTimeShiftedSignalTruncation_Wavelet]=Gabor(dt_POU,GaussianLength,GausiannVariance,GivenSignalLength,Scale,TimeOffset,Signal_Wavelet);
[GaborIntransformation_Wavelet]=InGabor(EffectiveAtoms_Start,EffectiveAtoms_End,GaborTimeShiftedSignalTruncation_Wavelet);
[~,TimeShiftedSignal_Seismic,GaborTimeShiftedSignalTruncation_Seismic]=Gabor(dt_POU,GaussianLength,GausiannVariance,GivenSignalLength,Scale,TimeOffset,Signal_Seismic);
[GaborIntransformation_Seismic]=InGabor(EffectiveAtoms_Start,EffectiveAtoms_End,GaborTimeShiftedSignalTruncation_Seismic);
Wavelet=[];
Seismic=[];
for i=EffectiveAtoms_Start:EffectiveAtoms_End
    Wavelet=[Wavelet,GaborTimeShiftedSignalTruncation_Wavelet(:,i)];
    Seismic=[Seismic,GaborTimeShiftedSignalTruncation_Seismic(:,i)];
end
%% POU-COM parameters and estimated ASSW   
df=f_seismic(2)-f_seismic(1);
Iterations=20;
mu=0;
VecterCut_start=[30;40;60;90];
VecterCut_end=[90;125;145;152];
Vecterp=[0.4;0.3;0.8;0.5];
flip=[1,1,-1,-1];
[Vectercp,Vecteralpha,Vecterbeta,k]=Parameter_determination_ALL(EffectiveAtoms_Start,EffectiveAtoms_End,df,Seismic,Vecterp,VecterCut_start,VecterCut_end,mu,flip);
data=Seismic;
for i=1:Iterations
    [InversionResult,AfterFitting]=Gaborfitting_ALL(data,df,Vecterp,VecterCut_start,VecterCut_end,Vectercp,Vecteralpha,Vecterbeta,k,EffectiveAtoms_Start,EffectiveAtoms_End,flip);
    data=AfterFitting;
end
RestoringSignal=InversionResult;
%% zero-phase seismic wavelet (ZSW)
[wavelet]=Amplitude_To_TimeSequence(dt,Signal_Wavelet,fmax,Conv_Length,length(wavelet));
[realx_pou]=Amplitude_To_TimeSequence(dt,RestoringSignal,fmax,Conv_Length,length(wavelet));
[realx_Contraction]=Amplitude_To_TimeSequence(dt,Contraction,fmax,Conv_Length,length(wavelet));
%% misfit function
ASSW_COM=Contraction./norm(Contraction,2)-Signal_Wavelet./norm(Signal_Wavelet,2);
ASSW_POU=RestoringSignal./norm(RestoringSignal,2)-Signal_Wavelet./norm(Signal_Wavelet,2);
ASSW_COM_time=realx_Contraction./norm(realx_Contraction,2)-wavelet./norm(wavelet,2);
ASSW_POU_time=realx_pou./norm(realx_pou,2)-wavelet./norm(wavelet,2);
i_plot
ASSW_Ormsby=max(abs(ASSW_POU))
%% 
F1=figure;
set(F1,'position',[20 10 750 980]);
subplot(7,4,[1,2]);
hold on;box on;
plot(t_seismic,r,'b','linewidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(a)'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.855 0.39 0.135]);
subplot(7,4,[3,4]);
hold on;box on;
plot(t_seismic,seismic,'b','LineWidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(b)'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.58 0.855 0.39 0.135]);
subplot(7,4,[5,6,9,10,13,14,17,18]);
hold on;box on;
plot(f_wavelet,amplitude_spectrum_wavelet,'r','LineWidth',1.5);
plot(f_seismic,amplitude_spectrum_seismic,'b','LineWidth',1);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(c)'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.515 0.39 0.26]);

for i=1:number
    xx=0:0:0;
    subplot(7,4,4*(i-1)+7);
    hold on;box on;
    Synthetic_seismic_record=plot(Seismic(:,i),'b','LineWidth',1);
    Seismic_wavelet=plot(Wavelet(:,i),'r','LineWidth',1.5);
    plot(AfterFitting(:,i)*k(i),'--k','LineWidth',1.5);
    ylabel('Magnitude');
    if i==number
    xlabel('Frequency/Hz');
    else
    set(gca,'xtick',xx);
    end
    set(gca,'FontName','Arial','FontSize',9.5,'linewidth',1.5,'FontWeight','bold');
    set(gca,'TickLength',[0 0.0001]); 
    set(gca,'position',[0.575 0.725-(i-1)*0.072 0.17 0.05]);
    subplot(7,4,4*(i-1)+8);
    hold on;box on;
    plot(Wavelet(:,i),'r','LineWidth',2);plot(AfterFitting(:,i)*k(i),'--k','LineWidth',2);
    if i==number
    xlabel('Frequency/Hz');
    else
    set(gca,'xtick',xx);
    end    
    set(gca,'FontName','Arial','FontSize',9.5,'linewidth',1.5,'FontWeight','bold');
    set(gca,'TickLength',[0 0.0001]); 
    set(gca,'position',[0.805 0.725-(i-1)*0.072 0.17 0.05]);
end

subplot(7,4,[21,22]);
hold on;box on;
plot(f_wavelet,Signal_Wavelet,'r','LineWidth',1.5);
plot(f_seismic,Contraction,'--g','LineWidth',1.5);
plot(f_seismic,RestoringSignal,'--k','LineWidth',1.5);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(e)'});
legend('Ture','COM','POU-COM','location','NorthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.298 0.39 0.135]);
subplot(7,4,[23,24]);
hold on;box on;
plot(t_wavelet,wavelet,'r','LineWidth',1.5);
plot(t_wavelet,realx_Contraction,'--g','LineWidth',1.5);
plot(t_wavelet,realx_pou,'--k','LineWidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(f)'});
legend('Ture','COM','POU-COM','location','NorthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.585 0.298 0.39 0.135]);
subplot(7,4,[25,26]);
hold on;box on;
plot(f_seismic,ASSW_COM,'b','linewidth',1.5);
plot(f_seismic,ASSW_POU,'--r','linewidth',1.5);
plot(f_seismic,zeros(size(ASSW_COM)),'k','linewidth',1.5);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(g)'});
legend('COM','POU-COM','location','SouthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.08 0.39 0.135]);
subplot(7,4,[27,28]);
hold on;box on;
plot(t_wavelet,ASSW_COM_time,'b','linewidth',1.5);
plot(t_wavelet,ASSW_POU_time,'--r','linewidth',1.5);
plot(t_wavelet,zeros(size(ASSW_POU_time)),'k','linewidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(h)'});
legend('COM','POU-COM','location','SouthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.585 0.08 0.39 0.135]);
annotation(F1,'textbox',...
    [0.749928571428572 0.44183673527168 0.0657142841815949 0.0295918361568938],...
    'String','(d)',...
    'LineStyle','none',...
    'Interpreter','none',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off');
Ormsby1=wavelet;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POU-COM schematic diagram

F=figure;
set(F,'position',[20 50 750 700]);
subplot(2,2,1);
hold on;box on;
plot(f_wavelet,amplitude_spectrum_wavelet,'r','LineWidth',1.5);
plot(f_seismic,amplitude_spectrum_seismic,'b','LineWidth',1);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(a)'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.09 0.62 0.4 0.364]);
subplot(2,2,2);
hold on;box on;
Synthetic_seismic_record=plot(amplitude_spectrum_seismic,'b','LineWidth',1);
for i=1:length(NormalizedGaborAtomicLibrary(1,:))
    if i>=EffectiveAtoms_Start && i<=EffectiveAtoms_End
        Effective_Gabor_atoms=plot(f_wavelet,0.04*NormalizedGaborAtomicLibrary(TimeOffset+1:TimeOffset+length(f_wavelet),i),'g','LineWidth',1);
    else
        Gabor_atoms=plot(f_wavelet,0.04*NormalizedGaborAtomicLibrary(TimeOffset+1:TimeOffset+length(f_wavelet),i),'k','LineWidth',1);
    end
end
Seismic_wavelet=plot(f_wavelet,amplitude_spectrum_wavelet,'r','LineWidth',1.5);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(b)'});
xlim([0,200]);
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001])
set(gca,'position',[0.583 0.62 0.4 0.364]);
for i=1:number
    xx=0:0:0;
    subplot(8,4,4*(i-1)+17);
    hold on;box on;
    Synthetic_seismic_record=plot(Seismic(:,i),'b','LineWidth',1);
    Seismic_wavelet=plot(Wavelet(:,i),'r','LineWidth',1.5);
    plot(AfterFitting(:,i)*k(i),'--k','LineWidth',1.5);
    ylabel('Magnitude');
    if i==number
    xlabel('Frequency/Hz');
    else
    set(gca,'xtick',xx);
    end
    set(gca,'FontName','Arial','FontSize',9.5,'linewidth',1.5,'FontWeight','bold');
    set(gca,'TickLength',[0 0.0001]); 
    set(gca,'position',[0.085 0.412-(i-1)*0.1008 0.17 0.07]);
    subplot(8,4,4*(i-1)+18);
    hold on;box on;
    plot(Wavelet(:,i),'r','LineWidth',2);plot(AfterFitting(:,i)*k(i),'--k','LineWidth',2);
    if i==number
    xlabel('Frequency/Hz');
    else
    set(gca,'xtick',xx);
    end
    set(gca,'FontName','Arial','FontSize',9.5,'linewidth',1.5,'FontWeight','bold');
    set(gca,'TickLength',[0 0.0001]); 
    set(gca,'position',[0.315 0.412-(i-1)*0.1008 0.17 0.07]);
end
subplot(2,2,4);
hold on;box on;
plot(f_wavelet,Signal_Wavelet,'r','LineWidth',1.5);
plot(f_seismic,RestoringSignal,'--k','LineWidth',1.5);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(d)'});
ylim([0,0.02]);
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.583 0.11 0.4 0.3861]);

annotation(F,'textbox',...
    [0.258678571428571 0.0208928571428571 0.04775 0.02625],'String',{'(c)'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off');
%%

%Ormsby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Klauder

L=40;
[t_wavelet,wavelet,f_wavelet,amplitude_spectrum_wavelet]=Klauder_my(dt,L,30,140,6,fa,fmax);
waveletname='Klauder';
%%
seismic_original=conv(r,wavelet);
seismic=seismic_original(L+1:end-L);
Conv_Length=length(seismic_original);
[f_seismic,F_seismic,amplitude_spectrum_seismic,N_length]=Amplitude_spectrum_my_LengthFix(dt,seismic_original,fmax,Conv_Length);
[f_wavelet,F_wavelet,amplitude_spectrum_wavelet,~]=Amplitude_spectrum_my_LengthFix(dt,wavelet,fmax,Conv_Length);
mu=0;
amplitude_spectrum_wavelet=amplitude_spectrum_wavelet';
%%
% amplitude_spectrum_wavelet=amplitude_spectrum_wavelet.*exp(-f_wavelet'/800);
amplitude_spectrum_wavelet=amplitude_spectrum_wavelet.*exp(-f_wavelet'*pi*delt/Q);
[wavelet]=Amplitude_To_TimeSequence(dt,amplitude_spectrum_wavelet,fmax,Conv_Length,length(wavelet));

seismic_original=conv(r,wavelet);
seismic=seismic_original(L+1:end-L);
Conv_Length=length(seismic_original);
[f_seismic,F_seismic,amplitude_spectrum_seismic,N_length]=Amplitude_spectrum_my_LengthFix(dt,seismic_original,fmax,Conv_Length);
[f_wavelet,F_wavelet,amplitude_spectrum_wavelet,~]=Amplitude_spectrum_my_LengthFix(dt,wavelet,fmax,Conv_Length);

mu=0;
amplitude_spectrum_wavelet=amplitude_spectrum_wavelet';
amplitude_spectrum_seismic=amplitude_spectrum_seismic';
m=inv(amplitude_spectrum_wavelet'*amplitude_spectrum_wavelet+mu)*amplitude_spectrum_wavelet'*amplitude_spectrum_seismic;
amplitude_spectrum_seismic=amplitude_spectrum_seismic/m;
%% COM parameters and estimated ASSW    
mu=0; %正则化参数 
Iterations=5; %迭代次数
p=0.4;
Contraction_start=20;
Contraction_end=160;
Contraction=zeros(size(amplitude_spectrum_seismic));
Contraction(Contraction_start:Contraction_end)=Contraction_mapping(f_seismic(Contraction_start:Contraction_end),amplitude_spectrum_seismic(Contraction_start:Contraction_end),p,Iterations,mu);
%%
Signal_Wavelet=amplitude_spectrum_wavelet;
Signal_Seismic=amplitude_spectrum_seismic;
%% POU decomposition
dt_POU=2;
GaussianLength=100;
GausiannVariance=40;
GivenSignalLength=300;
Scale=40;
TimeOffset=150;
EffectiveAtoms_Start=3;
EffectiveAtoms_End=6;
number=EffectiveAtoms_End-EffectiveAtoms_Start+1;
[NormalizedGaborAtomicLibrary,TimeShiftedSignal_Wavelet,GaborTimeShiftedSignalTruncation_Wavelet]=Gabor(dt_POU,GaussianLength,GausiannVariance,GivenSignalLength,Scale,TimeOffset,Signal_Wavelet);
[GaborIntransformation_Wavelet]=InGabor(EffectiveAtoms_Start,EffectiveAtoms_End,GaborTimeShiftedSignalTruncation_Wavelet);
[~,TimeShiftedSignal_Seismic,GaborTimeShiftedSignalTruncation_Seismic]=Gabor(dt_POU,GaussianLength,GausiannVariance,GivenSignalLength,Scale,TimeOffset,Signal_Seismic);
[GaborIntransformation_Seismic]=InGabor(EffectiveAtoms_Start,EffectiveAtoms_End,GaborTimeShiftedSignalTruncation_Seismic);
Wavelet=[];
Seismic=[];
for i=EffectiveAtoms_Start:EffectiveAtoms_End
    Wavelet=[Wavelet,GaborTimeShiftedSignalTruncation_Wavelet(:,i)];
    Seismic=[Seismic,GaborTimeShiftedSignalTruncation_Seismic(:,i)];
end
%% POU-COM parameters and estimated ASSW   
df=f_seismic(2)-f_seismic(1);
Iterations=20;
mu=0;
VecterCut_start=[25;35;65;100];
VecterCut_end=[80;120;150;160];
Vecterp=[0.4;0.2;0.05;0.6];
flip=[1,1,-1,-1];
[Vectercp,Vecteralpha,Vecterbeta,k]=Parameter_determination_ALL(EffectiveAtoms_Start,EffectiveAtoms_End,df,Seismic,Vecterp,VecterCut_start,VecterCut_end,mu,flip);
data=Seismic;
for i=1:Iterations
    [InversionResult,AfterFitting]=Gaborfitting_ALL(data,df,Vecterp,VecterCut_start,VecterCut_end,Vectercp,Vecteralpha,Vecterbeta,k,EffectiveAtoms_Start,EffectiveAtoms_End,flip);
    data=AfterFitting;
end
RestoringSignal=InversionResult;
%% zero-phase seismic wavelet (ZSW)
[wavelet]=Amplitude_To_TimeSequence(dt,Signal_Wavelet,fmax,Conv_Length,length(wavelet));
[realx_pou]=Amplitude_To_TimeSequence(dt,RestoringSignal,fmax,Conv_Length,length(wavelet));
[realx_Contraction]=Amplitude_To_TimeSequence(dt,Contraction,fmax,Conv_Length,length(wavelet));
%% misfit function
ASSW_COM=Contraction./norm(Contraction,2)-Signal_Wavelet./norm(Signal_Wavelet,2);
ASSW_POU=RestoringSignal./norm(RestoringSignal,2)-Signal_Wavelet./norm(Signal_Wavelet,2);
ASSW_COM_time=realx_Contraction./norm(realx_Contraction,2)-wavelet./norm(wavelet,2);
ASSW_POU_time=realx_pou./norm(realx_pou,2)-wavelet./norm(wavelet,2);
ASSW_Klauder=max(abs(ASSW_POU))
%% 
F1=figure;
set(F1,'position',[20 10 750 980]);
subplot(7,4,[1,2]);
hold on;box on;
plot(t_seismic,r,'b','linewidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(a)'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.855 0.39 0.135]);
subplot(7,4,[3,4]);
hold on;box on;
plot(t_seismic,seismic,'b','LineWidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(b)'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.58 0.855 0.39 0.135]);
subplot(7,4,[5,6,9,10,13,14,17,18]);
hold on;box on;
plot(f_wavelet,amplitude_spectrum_wavelet,'r','LineWidth',1.5);
plot(f_seismic,amplitude_spectrum_seismic,'b','LineWidth',1);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(c)'});
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.515 0.39 0.26]);

for i=1:number
    xx=0:0:0;
    subplot(7,4,4*(i-1)+7);
    hold on;box on;
    Synthetic_seismic_record=plot(Seismic(:,i),'b','LineWidth',1);
    Seismic_wavelet=plot(Wavelet(:,i),'r','LineWidth',1.5);
    plot(AfterFitting(:,i)*k(i),'--k','LineWidth',1.5);
    ylabel('Magnitude');
    if i==number
    xlabel('Frequency/Hz');
    else
    set(gca,'xtick',xx);
    end
    set(gca,'FontName','Arial','FontSize',9.5,'linewidth',1.5,'FontWeight','bold');
    set(gca,'TickLength',[0 0.0001]); 
    set(gca,'position',[0.575 0.725-(i-1)*0.072 0.17 0.05]);
    subplot(7,4,4*(i-1)+8);
    hold on;box on;
    plot(Wavelet(:,i),'r','LineWidth',2);plot(AfterFitting(:,i)*k(i),'--k','LineWidth',2);
    if i==number
    xlabel('Frequency/Hz');
    else
    set(gca,'xtick',xx);
    end    
    set(gca,'FontName','Arial','FontSize',9.5,'linewidth',1.5,'FontWeight','bold');
    set(gca,'TickLength',[0 0.0001]); 
    set(gca,'position',[0.805 0.725-(i-1)*0.072 0.17 0.05]);
end

subplot(7,4,[21,22]);
hold on;box on;
plot(f_wavelet,Signal_Wavelet,'r','LineWidth',1.5);
plot(f_seismic,Contraction,'--g','LineWidth',1.5);
plot(f_seismic,RestoringSignal,'--k','LineWidth',1.5);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(e)'});
legend('Ture','COM','POU-COM','location','NorthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.298 0.39 0.135]);
subplot(7,4,[23,24]);
hold on;box on;
plot(t_wavelet,wavelet,'r','LineWidth',1.5);
plot(t_wavelet,realx_Contraction,'--g','LineWidth',1.5);
plot(t_wavelet,realx_pou,'--k','LineWidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(f)'});
legend('Ture','COM','POU-COM','location','NorthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.585 0.298 0.39 0.135]);
subplot(7,4,[25,26]);
hold on;box on;
plot(f_seismic,ASSW_COM,'b','linewidth',1.5);
plot(f_seismic,ASSW_POU,'--r','linewidth',1.5);
plot(f_seismic,zeros(size(ASSW_COM)),'k','linewidth',1.5);
ylabel('Magnitude');
xlabel({'Frequency/Hz';'(g)'});
legend('COM','POU-COM','location','SouthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.08 0.39 0.135]);
subplot(7,4,[27,28]);
hold on;box on;
plot(t_wavelet,ASSW_COM_time,'b','linewidth',1.5);
plot(t_wavelet,ASSW_POU_time,'--r','linewidth',1.5);
plot(t_wavelet,zeros(size(ASSW_POU_time)),'k','linewidth',1.5);
ylabel('Amplitude');
xlabel({'Time/s';'(h)'});
legend('COM','POU-COM','location','SouthEast','FontName','Arial');
set(gca,'FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.585 0.08 0.39 0.135]);
annotation(F1,'textbox',...
    [0.749928571428572 0.44183673527168 0.0657142841815949 0.0295918361568938],...
    'String','(d)',...
    'LineStyle','none',...
    'Interpreter','none',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off');
Klauder1=wavelet;
end
end
%%
save('..\reflection coefficients and wavelet\Ricker1.mat','Ricker1');
save('..\reflection coefficients and wavelet\Ormsby1.mat','Ormsby1');
save('..\reflection coefficients and wavelet\Klauder1.mat','Klauder1');