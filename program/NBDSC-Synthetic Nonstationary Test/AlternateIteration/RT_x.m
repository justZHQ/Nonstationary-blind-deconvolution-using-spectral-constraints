function [y] = RT_x(R,x,length_wavelet)
% input:
% x          is the reflection coefficient to be inverted (n*m)
% wavelet    is the seismic wavelet constituting the Toeplitz wavelet matrix A

% output:
% y          is the convolution result of wavelet and reflection
%            coefficient, that is, the super gather of seismic records (nm*1)
wave_L=(length_wavelet-1)/2;
[n,m]=size(R);
y=[];
for k=1:m
    RR=[];
    %     for j=1:m
    Rmid=zeros(length_wavelet+n-1,length_wavelet);
    for i=1:length(Rmid(1,:))
        Rmid(i:n+i-1,i)=R(:,k);
    end
    RR=Rmid(wave_L+1:end-wave_L,:);
    y=[y,RR'*x(:,k)];
    %         RR=[RR;Rmid];
    %     end
%     aa=1;
end


% r=r(wave_L+1:end-wave_L,:);
% r_flip=fliplr(r);
% [n,m]=size(r);
% y=zeros(n,m);
% for i=1:m
%     y(:,i)=conv(r_flip(:,i),x(:,i),'same');
% end

% y=conv2(wavelet_flip,1,x,'same');