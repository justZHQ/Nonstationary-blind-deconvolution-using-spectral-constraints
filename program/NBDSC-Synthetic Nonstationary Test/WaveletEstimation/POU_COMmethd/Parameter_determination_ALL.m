function [Vecter_cp,Vecter_alpha,Vecter_beta,k]=Parameter_determination_ALL(EffectiveAtoms_Start,EffectiveAtoms_End,df,Vecter_data,Vecter_p,VecterCut_start,VecterCut_end,mu,flip)

VecterLength=length(Vecter_data(1,:));
Vecter_cp=zeros(VecterLength,1);
Vecter_alpha=zeros(VecterLength,1);
Vecter_beta=zeros(VecterLength,1);


AfterFitting=zeros(size(Vecter_data));
for i=1:VecterLength
    Signal=Vecter_data(:,i);
    if flip(i)==-1
        Signal=Signal(VecterCut_end(i):-1:VecterCut_start(i));
        [Vecter_cp(i),Vecter_alpha(i),Vecter_beta(i)]=Parameter_determination_EachP(df,Signal,Vecter_p(i),mu);
        [~,AfterFitting(:,i)]=Gaborfitting_EachP_Flip(Vecter_data(:,i),df,Vecter_p(i),VecterCut_start(i),VecterCut_end(i),Vecter_cp(i),Vecter_alpha(i),Vecter_beta(i));
    elseif flip(i)==1
        Signal=Signal(VecterCut_start(i):VecterCut_end(i));    
        [Vecter_cp(i),Vecter_alpha(i),Vecter_beta(i)]=Parameter_determination_EachP(df,Signal,Vecter_p(i),mu);
        [~,AfterFitting(:,i)]=Gaborfitting_EachP(Vecter_data(:,i),df,Vecter_p(i),VecterCut_start(i),VecterCut_end(i),Vecter_cp(i),Vecter_alpha(i),Vecter_beta(i));
    end   
end

[~,k]=result(EffectiveAtoms_Start,EffectiveAtoms_End,AfterFitting,0,Vecter_data);
