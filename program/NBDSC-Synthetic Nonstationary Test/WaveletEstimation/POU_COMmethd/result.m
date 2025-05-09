function [InversionResult,m]=result(EffectiveAtoms_Start,EffectiveAtoms_End,GaborTimeShiftedSignalTruncation_AfterFitting,mu,Seismic)
%%
% GaussianLength=(length(GaborTimeShiftedSignalTruncation_AfterFitting(:,1))-1)/2;
% G=zeros(Time(EffectiveAtoms_End)-Time(EffectiveAtoms_Start)+2*GaussianLength+1+1,EffectiveAtoms_End-EffectiveAtoms_Start+1);
% GaborIntransformationAfterFitting=zeros(size(G(:,1)));%%%
% 
% S=zeros(size(G));
% for j=EffectiveAtoms_Start:EffectiveAtoms_End
%     t=Time(j)-Time(EffectiveAtoms_Start)+1;
%     G(t:t+2*GaussianLength+1-1,j-EffectiveAtoms_Start+1)=GaborTimeShiftedSignalTruncation_AfterFitting(:,j-EffectiveAtoms_Start+1);
%     GaborIntransformationAfterFitting=GaborIntransformationAfterFitting+G(:,j-EffectiveAtoms_Start+1);
%     S(t:t+2*GaussianLength+1-1,j-EffectiveAtoms_Start+1)=Seismic(:,j-EffectiveAtoms_Start+1);%%%
% end
% TimeShiftedSignalTruncation=TimeShiftedSignal(Time(EffectiveAtoms_Start):Time(EffectiveAtoms_End)+2*GaussianLength+1);
% % m=inv(G'*G+mu*eye(EffectiveAtoms_End-EffectiveAtoms_Start+1))*G'*TimeShiftedSignalTruncation;
% 
n=EffectiveAtoms_End-EffectiveAtoms_Start+1;
% m=zeros(1,n);
% SS=zeros(size(G(:,1)));
G=GaborTimeShiftedSignalTruncation_AfterFitting;
S=Seismic;
SS=zeros(size(G(:,1)));
for i=1:n
    m(i)=inv(G(:,i)'*G(:,i)+mu)*G(:,i)'*S(:,i);
    SS=SS+S(:,i);
%     
% ss=0;
% pp=0;
% for j=1:length(G(:,1))
%     ss=ss+S(j,i)^2;
%     pp=pp+G(j,i)^2;
% end
% m(i)=ss/pp;
    
end




% figure;
% hold on;
% plot(TimeShiftedSignalTruncation);
% plot(SS);

% m=inv(G'*G+mu*eye(EffectiveAtoms_End-EffectiveAtoms_Start+1))*G'*SS;
InversionResult=zeros(size(G(:,1)));
for i=1:EffectiveAtoms_End-EffectiveAtoms_Start+1
    InversionResult=G(:,i)*m(i)+InversionResult;
end