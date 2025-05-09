function [InputSignalCut,AfterFitting]=Gaborfitting_EachP_Flip(InputSignal,df,p,Cut_start,Cut_end,cp,alpha,beta)
InputSignalCut=InputSignal(Cut_start:Cut_end);
n=length(InputSignalCut);
S_new=abs(InputSignalCut).^p+0;
S0=S_new./sum(S_new)*df;
F=zeros(n,1);
F(1)=S0(1);
for i=2:n
    F(i)=F(i-1)+S0(i)*df;
end
P=(cp*F.^alpha).*(1-F).^beta;
InputSignal=zeros(size(InputSignal));
InputSignal(Cut_end:-1:Cut_start)=P;
AfterFitting=InputSignal;