function [AMP_Wavelet,ASSW_POUCOM,EstimationWave_POUCOM]=WaveletEstimation(wavelet,dt,Y,fmax,FFT_Length,WaveLength_Give)

[phase_rot]=ConstPhaseRotation(Y);
%%
[~,trcNum]=size(Y);
WaveLength_Give=WaveLength_Give*2+1;

Fa0 = 0;
for k = 1 : trcNum
    a=Y(:,k);
    [f_seismic,~,Fak,~]=Amplitude_spectrum_my_LengthFix(dt,a,fmax,FFT_Length);
     Fa0 = Fa0 + Fak;
end
Fa=Fa0;

[~,~,amplitude_spectrum_wavelet,~]=Amplitude_spectrum_my_LengthFix(dt,wavelet,fmax,FFT_Length);
amplitude_spectrum_seismic=Fa;

mu=0;
amplitude_spectrum_wavelet=amplitude_spectrum_wavelet';
m=inv(amplitude_spectrum_wavelet'*amplitude_spectrum_wavelet+mu)*amplitude_spectrum_wavelet'*amplitude_spectrum_seismic;
amplitude_spectrum_seismic=amplitude_spectrum_seismic/m;

%%
AMP_Wavelet=amplitude_spectrum_wavelet;
AMP_Seismic=amplitude_spectrum_seismic;
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
[NormalizedGaborAtomicLibrary,TimeShiftedSignal_Wavelet,GaborTimeShiftedSignalTruncation_Wavelet]=Gabor(dt_POU,GaussianLength,GausiannVariance,GivenSignalLength,Scale,TimeOffset,AMP_Wavelet);
[GaborIntransformation_Wavelet]=InGabor(EffectiveAtoms_Start,EffectiveAtoms_End,GaborTimeShiftedSignalTruncation_Wavelet);
[~,TimeShiftedSignal_Seismic,GaborTimeShiftedSignalTruncation_Seismic]=Gabor(dt_POU,GaussianLength,GausiannVariance,GivenSignalLength,Scale,TimeOffset,AMP_Seismic);
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
ASSW_POUCOM=RestoringSignal;
phase_rot=-0*pi/6;
[EstimationWave_POUCOM]=Amplitude_To_TimeSequence_phase(dt,RestoringSignal,fmax,FFT_Length,WaveLength_Give,-phase_rot);
