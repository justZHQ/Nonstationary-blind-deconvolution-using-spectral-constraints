function [log_ampx_Correction,ampx_Correction]=Fast_COM(lamda1,lamda2,ampx,df,p,Iterations,Begf,Endf,Begf_high,Endf_low,GausiannVariance)

% n=length(f_seismic); df=f_seismic(2)-f_seismic(1);
[~,n]=size(ampx);
A_all=[];
y=[];
A_allcut=[];
log_ampx=2*log(ampx);

BegNum_high=floor(Begf_high/df)+1; 
EndNum_low=floor(Endf_low/df)+1; 
ma=EndNum_low-BegNum_high+1;
for j=1:n
    BegNum=floor(Begf(j)/df)+1; 
    EndNum=floor(Endf(j)/df)+1;
    S_new=abs(ampx(BegNum:EndNum,j)).^p(j);
    S0=S_new./sum(S_new)*df;
    m=EndNum-BegNum+1;
    F=zeros(m,1);
    F(1)=S0(1);
    for i=2:m
        F(i)=F(i-1)+S0(i)*df;
    end
    A1=ones(m,1);A2=log(F);A3=log(1-F);
    A0=[A1,A2,A3];
    A_allj=zeros(m,3*n);
    A_allj(:,(j-1)*3+1:j*3)=A0;
    A_all=[A_all;A_allj];
    y=[y;log_ampx(BegNum:EndNum,j)];
    A_allcutj=zeros(EndNum_low-BegNum_high+1,3*n);
    A_allcutj(:,(j-1)*3+1:j*3)=A0(BegNum_high:EndNum_low,:);
    A_allcut=[A_allcut;A_allcutj]; 
end

%%
L1=zeros(ma-1,ma);
for i=1:ma-1
    L1(i,i)=-1;
    L1(i,i+1)=1;
end

L2=zeros(ma-2,ma);
for i=1:ma-2
    L2(i,i)=1;
    L2(i,i+1)=-2;
    L2(i,i+2)=1;
end

% Para=[3 2 1.5 1 1 1 1 1 1.5 2 3]*lamda2;
% Para=[0.2, 0.2, 0.4, 0.6, 0.8, 1, 0.8, 0.6, 0.4, 0.2, 0.2]*lamda2;
% Lam2=[];
% for i=1:n
%     Lam=[Para(i); Para(i); Para(i)];
%     Lam2=[Lam;Lam2];
% end
% lamda2=diag(Lam2);

% GausiannVariance=5;
averge=5;
T=-floor(n/2)+averge:floor(n/2)+averge;
Para=exp(-T.^2/(2*GausiannVariance^2))*lamda2;
% figure;
% plot(Para);
Lam2=[];
for i=1:n
%     if k~=0
%         Lam=[Para(i); Para(i)];
%     else
        Lam=[Para(i); Para(i); Para(i)];
%     end
    Lam2=[Lam;Lam2];
end
lamda2=diag(Lam2);

%%
Gall=zeros(3*n,3*n);
K=1;
N=n-K;

for i=1:N
    Pi=zeros(ma,ma*n);
    for j=1:ma
%         if i==1
        Pi(j,(i-1)*ma+j)=-1;
%         end
%         Pi(j,j)=-1;
        Pi(j,(i-1+K)*ma+j)=1;
    end
    Gii=L1*Pi*A_allcut;
        Gi=L2*Pi*A_allcut;
%     Gall=Gall+(Gii'*Gii);
        Gall=Gall+lamda1*(Gi'*Gi)+lamda2*(Gii'*Gii);
end

y=A_all'*y;
invG=inv(Gall+1*(A_all'*A_all)+0.0*eye(3*n));
M=invG*y;
M=reshape(M,3,n);

Sk=zeros(size(ampx));
for k=1:Iterations
    for i=1:n
        BegNum=floor(Begf(i)/df)+1; 
        EndNum=floor(Endf(i)/df)+1;
        m=EndNum-BegNum+1;
        cp=exp(M(1,i));alpha=M(2,i);beta=M(3,i);
        Sk(BegNum:EndNum,i)=ampx(BegNum:EndNum,i).^p(i);
        Sk(BegNum:EndNum,i)=Sk(BegNum:EndNum,i)./sum(Sk(BegNum:EndNum,i))*df;
        F=zeros(m,n);
        F(1,i)=Sk(BegNum,i);
        for j=2:m
            F(j,i)=F(j-1,i)+Sk(BegNum+j-1,i)*df;
        end
        P=(cp*F(:,i).^alpha).*(1-F(:,i)).^beta;
        Sk(BegNum:EndNum,i)=P;
%         Sk(:,i)=Sk(:,i);
    end
end
ampx_Correction=Sk;
log_ampx_Correction=log(ampx_Correction);

