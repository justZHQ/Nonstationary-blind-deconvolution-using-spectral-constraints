function [R_hat,InverseWave,Convergence]=Alternate_iteration_SpectralConstraints_AWPSI(F_cut,b,lambda0,lambda1,lambda2,Y,e_wave,L0,eta,lambda,iter_alternate,iter_fista,opt)
%%
%%
[~,m]=size(Y);
% length_wavelet=length(e_wave);
% wave_L=(length_wavelet-1)/2;

R0 = 1e-3*ones(size(Y));
%----------------------estimate A and R--------------------%
R_hat = zeros(size(Y));
Convergence=zeros(1,iter_alternate);

wavelet=e_wave./max(e_wave,[],1);
% tic
for iter = 1:iter_alternate
%% A
    [R_hat,~]= fista_backtracking_lasso(wavelet,Y,R0,L0,eta,lambda,iter_fista,opt);
%     disp('Lasso I done.');
%%  
    [InverseWave,~] = SpectralConstraints_FISTA_con(R_hat,F_cut,Y,b,wavelet,L0,eta,lambda0,lambda1,lambda2,iter,opt,eps);
%     disp('Lasso I done.');
%%   
    wavelet_hat=InverseWave./max(InverseWave,[],1);
%     R_hat = R_hat./max(InverseWave,[],1);   
%% 
Y_now=zeros(size(Y));
for i=1:m
    Y_now(:,i)=conv(R_hat(:,i),wavelet_hat(:,i),'same');
end
%     Y_now=conv2(wavelet_hat,2,R_hat,'same');
    Convergence(iter)=0.5*sum(sum((Y-Y_now).^2));   
    if norm(wavelet_hat-wavelet)/norm(wavelet) < 1e-4 && norm(R0-R_hat,'fro')/norm(R0,'fro') < 1e-4
        break;
    else
        wavelet=wavelet_hat;
        R0 = R_hat;
    end
%     iter
end
% toc
