function [GaborIntransformation]=InGabor(EffectiveAtoms_Start,EffectiveAtoms_End,GaborTimeShiftedSignalTruncation)

% GaborIntransformation=zeros(size(GaborTimeShiftedSignal(:,1)));
% for j=EffectiveAtoms_Start:EffectiveAtoms_End
%     GaborIntransformation=GaborIntransformation+GaborTimeShiftedSignal(:,j);
% end
% GaborIntransformation=GaborIntransformation(Time(EffectiveAtoms_Start):Time(EffectiveAtoms_End)+GaussianLength);
% TimeShiftedSignalTruncation=TimeShiftedSignal(Time(EffectiveAtoms_Start):Time(EffectiveAtoms_End)+GaussianLength);

%%
G=zeros(length(GaborTimeShiftedSignalTruncation(:,1)),EffectiveAtoms_End-EffectiveAtoms_Start+1);
GaborIntransformation=zeros(size(G(:,1)));
for j=EffectiveAtoms_Start:EffectiveAtoms_End
%     t=Time(j)-Time(EffectiveAtoms_Start)+1;
%     G(t:t+2*GaussianLength+1-1,j)=GaborTimeShiftedSignalTruncation(:,j);
    GaborIntransformation=GaborIntransformation+GaborTimeShiftedSignalTruncation(:,j);    
end
% TimeShiftedSignalTruncation=TimeShiftedSignal(Time(EffectiveAtoms_Start):Time(EffectiveAtoms_End)+2*GaussianLength+1);
