function [R_hat,InverseWave,Convergence]=Alternate_iteration(Y,e_wave,L0,eta,lambda,iter_alternate,iter_fista,opt)
%%
[n,m]=size(Y);
length_wavelet=length(e_wave);
wave_L=(length_wavelet-1)/2;

R0 = 1e-3*ones(size(Y));
%----------------------estimate A and R--------------------%
% param.lambda=0.1;
% param.lambda2=1.5;
param.lambda=0.001;
param.lambda2=10;
param.mode=2;
% param.ols = 'ture';
R_hat = zeros(size(Y));
Convergence=zeros(1,iter_alternate);

wavelet=e_wave./max(e_wave,[],1);
tic
for iter = 1:iter_alternate
%% A
    [R_hat,~]= fista_backtracking_lasso(wavelet,Y,R0,L0,eta,lambda,iter_fista,opt);
%     [R_hat,~] = fista_backtracking_lasso_LateraConstraints(nt,dt,flow,fhigh,PRE_F,PRE_B,wavelet,Y,R0,L0,eta,lambda1,lambda2,iter_fista,opt,eps);
    disp('Lasso I done.');
%% R
    RR=[];
    for j=1:m
        Rmid=zeros(length_wavelet+n-1,length_wavelet);
        for i=1:length(Rmid(1,:))
            Rmid(i:n+i-1,i)=R_hat(:,j);
        end
        Rmid=Rmid(wave_L+1:end-wave_L,:);
        RR=[RR;Rmid];
    end    
    InverseWave = mexLasso(reshape(Y,[],1),RR,param);   
    InverseWave=full(InverseWave);
    disp('Lasso I done.');
%%   
    wavelet_hat=InverseWave/max(InverseWave);
    R_hat = R_hat*max(InverseWave);   
%% 
    Y_now=conv2(wavelet_hat,1,R_hat,'same');
    Convergence(iter)=0.5*sum(sum((Y-Y_now).^2))   
    if norm(wavelet_hat-wavelet)/norm(wavelet) < 1e-8 && norm(R0-R_hat,'fro')/norm(R0,'fro') < 1e-8
        break;
    else
        wavelet=wavelet_hat;
        R0 = R_hat;
    end
    iter
end
toc
