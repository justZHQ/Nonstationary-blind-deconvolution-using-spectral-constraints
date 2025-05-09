clear;clc;close all;

n=500;m=60;dt=0.001;
wave_L=40;
fmax=200;fa=1*pi/4;

[t_wavelet,wave_object,f_wavelet,amplitude_spectrum_wavelet1]=Ormsby_my(dt,wave_L,10,40,90,120,0,fmax);
[~,wave0,~,amplitude_spectrum_wavelet2]=Ricker_my(dt,wave_L,50,0,fmax);
[~,w_give,~,~]=Ricker_my(dt,wave_L,50,fa,fmax);w_give=w_give';
[~,wave_object_fa,~,~]=Ormsby_my(dt,wave_L,10,40,90,120,fa,fmax);wave_object_fa=wave_object_fa';
t_wavelet=t_wavelet*1000;
wave_object=wave_object';
wave0=wave0';
wave_n=length(wave_object);


Wavelength=length(wave_object);
FFT_length=512;
FrequencyBand=[5,130];
[f,F,F_cut]=FFT_matrix_GivenFrequencyBand(dt,Wavelength,FFT_length,FrequencyBand);
p=abs(F_cut*wave_object);




iter=250;
lambda1=1;

lambda2=0;
e=0.00001;
% e=0.0001;
k=[3,50];

[loss,x]=grad_new(F_cut,p,wave0,w_give,lambda1,lambda2,e,iter,k);
Fwinverse=F_cut*x;
gcf1=figure;
set(gcf1,'position',[800 600 400 250]);
plot(loss, 'k','LineWidth', 2);
ylabel('Loss');
xlabel('Number of iterations');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
annotation('textbox',[.002 .79 .1 .2], ...
    'String','a)','EdgeColor','none','FontSize',14,'FontWeight','bold');

gcf2=figure;
set(gcf2,'position',[800 600 400 250]);
hold on;box on;
plot(f,p.^2,'r','LineWidth', 2);
plot(f,abs(Fwinverse(:,1)).^2,'k','LineWidth', 2);
plot(f,abs(Fwinverse(:,2)).^2,'--g','LineWidth', 2);
plot(f,abs(Fwinverse(:,3)).^2,'-.b','LineWidth', 2);
ylabel('Magnitude');
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
legend('$\bf{p}$','Iter=0','Iter=10','Iter=200','location','NorthEast','FontName','Arial','interpreter','latex');
annotation('textbox',[.002 .79 .1 .2], ...
    'String','b)','EdgeColor','none','FontSize',14,'FontWeight','bold');

gcf3=figure;
set(gcf3,'position',[800 600 400 250]);
hold on;box on;
% plot(t_wavelet,wave_object_fa,'r','LineWidth', 2);
% plot(t_wavelet,w_give,'r','LineWidth', 2);
plot(t_wavelet,x(:,1),'k','LineWidth', 2);
plot(t_wavelet,x(:,2),'--g','LineWidth', 2);
plot(t_wavelet,x(:,3),'-.b','LineWidth', 2);
xlim([-25 25]);
% plot(t_wavelet,x(:,4),'-.g','LineWidth', 2);
% plot(t_wavelet,w_give,'k','LineWidth', 2);
ylabel('Amplitude');
xlabel('Time/ms');
legend('Iter=0','Iter=10','Iter=200','location','NorthEast','FontName','Arial','interpreter','latex');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
annotation('textbox',[.03 .79 .1 .2], ...
    'String','c)','EdgeColor','none','FontSize',14,'FontWeight','bold');
%%
lambda2=1200;
e=0.000005;
k=[10,200];
% lambda2=100;
[loss,x]=grad_new(F_cut,p,wave0,w_give,lambda1,lambda2,e,iter,k);
Fwinverse=F_cut*x;
gcf4=figure;
set(gcf4,'position',[800 600 400 250]);
plot(loss, 'k','LineWidth', 2);
ylabel('Loss');
xlabel('Number of iterations');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
annotation('textbox',[.002 .79 .1 .2], ...
    'String','d)','EdgeColor','none','FontSize',14,'FontWeight','bold');

gcf5=figure;
set(gcf5,'position',[800 600 400 250]);
hold on;box on;
plot(f,p.^2,'r','LineWidth', 2);
plot(f,abs(Fwinverse(:,1)).^2,'k','LineWidth', 2);
plot(f,abs(Fwinverse(:,2)).^2,'--g','LineWidth', 2);
plot(f,abs(Fwinverse(:,3)).^2,'-.b','LineWidth', 2);
ylabel('Magnitude');
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
legend('$\bf{p}$','Iter=0','Iter=10','Iter=200','location','NorthEast','FontName','Arial','interpreter','latex');
annotation('textbox',[.002 .79 .1 .2], ...
    'String','e)','EdgeColor','none','FontSize',14,'FontWeight','bold');

gcf6=figure;
set(gcf6,'position',[800 600 400 250]);
hold on;box on;
% plot(t_wavelet,wave_object_fa,'r','LineWidth', 2);
plot(t_wavelet,w_give,'m','LineWidth', 2);
plot(t_wavelet,x(:,1),'k','LineWidth', 2);
plot(t_wavelet,x(:,2),'--g','LineWidth', 2);
plot(t_wavelet,x(:,3),'-.b','LineWidth', 2);
xlim([-25 25]);
ylim([-0.67 1]);
% plot(t_wavelet,x(:,4),'-.g','LineWidth', 2);
% plot(t_wavelet,w_give,'k','LineWidth', 2);
ylabel('Amplitude');
xlabel('Time/ms');
legend('$\bf{\tilde w}$','Iter=0','Iter=10','Iter=200','location','NorthEast','FontName','Arial','interpreter','latex');
set(gca,'FontName','Arial','FontSize',12,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
annotation('textbox',[.03 .79 .1 .2], ...
    'String','f)','EdgeColor','none','FontSize',14,'FontWeight','bold');
