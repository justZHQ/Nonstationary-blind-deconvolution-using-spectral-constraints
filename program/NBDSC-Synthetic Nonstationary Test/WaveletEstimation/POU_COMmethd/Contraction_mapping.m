function [Sk_output]=Contraction_mapping(f_seismic,amplitude_spectrum_seismic,p,Iterations,mu)
n=length(f_seismic); df=f_seismic(2)-f_seismic(1);
S_new=abs(amplitude_spectrum_seismic).^p;
S0=S_new./sum(S_new)*df;
F=zeros(n,1);
F(1)=S0(1);
for i=2:n
    F(i)=F(i-1)+S0(i)*df;
end
% S=log(S0');
% A1=ones(n,1);A2=log(F);A3=log(1-F);
% A=[A1,A2,A3];
% m=inv(A'*A+mu*eye(3))*A'*S;
% cp=exp(m(1));alpha=m(2);beta=m(3);

S0=abs(amplitude_spectrum_seismic);
S=log(S0);
A1=ones(n,1);A2=log(F);A3=log(1-F);
A=[A1,A2,A3];
m=inv(A'*A+mu*eye(3))*A'*S;
cp=exp(m(1));alpha=m(2);beta=m(3);



Sk=S0;
for k=1:Iterations
    Sk=Sk.^p;
    Sk=Sk./sum(Sk)*df;
    F=zeros(n,1); 
    F(1)=Sk(1);
    for i=2:n
        F(i)=F(i-1)+Sk(i)*df;
    end
    P=(cp*F.^alpha).*(1-F).^beta;
    Sk=P; 
    Sk=Sk;
end
Sk_output=Sk;

% 
% InputSignalCut=InputSignal(Cut_start:Cut_end);
% n=length(InputSignalCut);
% S_new=abs(InputSignalCut).^p;
% S0=S_new./sum(S_new)*df;
% F=zeros(n,1);
% F(1)=S0(1);
% for i=2:n
%     F(i)=F(i-1)+S0(i)*df;
% end
% P=(cp*F.^alpha).*(1-F).^beta;
% InputSignal=zeros(size(InputSignal));
% InputSignal(Cut_start:Cut_end)=P;
% AfterFitting=InputSignal;