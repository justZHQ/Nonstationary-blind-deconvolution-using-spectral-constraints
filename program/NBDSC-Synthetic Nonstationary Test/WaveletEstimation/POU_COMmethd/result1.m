function [GaborIntransformationAfterFitting,TimeShiftedSignalTruncation,InversionResult,m]=result(TimeShiftedSignal,EffectiveAtoms_Start,EffectiveAtoms_End,GaborTimeShiftedSignalTruncation_AfterFitting,Time,mu)
%%
GaussianLength=(length(GaborTimeShiftedSignalTruncation_AfterFitting(:,1))-1)/2;
G=zeros(Time(EffectiveAtoms_End)-Time(EffectiveAtoms_Start)+2*GaussianLength+1+1,EffectiveAtoms_End-EffectiveAtoms_Start+1);
GaborIntransformationAfterFitting=zeros(size(G(:,1)));
for j=EffectiveAtoms_Start:EffectiveAtoms_End
    t=Time(j)-Time(EffectiveAtoms_Start)+1;
    G(t:t+2*GaussianLength+1-1,j-EffectiveAtoms_Start+1)=GaborTimeShiftedSignalTruncation_AfterFitting(:,j-EffectiveAtoms_Start+1);
    GaborIntransformationAfterFitting=GaborIntransformationAfterFitting+G(:,j-EffectiveAtoms_Start+1);
end
TimeShiftedSignalTruncation=TimeShiftedSignal(Time(EffectiveAtoms_Start):Time(EffectiveAtoms_End)+2*GaussianLength+1);
% m=inv(G'*G+mu*eye(EffectiveAtoms_End-EffectiveAtoms_Start+1))*G'*TimeShiftedSignalTruncation;

n=EffectiveAtoms_End-EffectiveAtoms_Start+1;
m=zeros(1,n);
GG=zeros(size(G(:,1)));
for i=1:n
    GG=GG+G(:,i);
end
for i=1:n
    m(i)=inv(GG'*GG+mu)*GG'*TimeShiftedSignalTruncation;
end    
    
    
InversionResult=zeros(size(G(:,1)));
for i=1:EffectiveAtoms_End-EffectiveAtoms_Start+1
    InversionResult=G(:,i)*m(i)+InversionResult;
end