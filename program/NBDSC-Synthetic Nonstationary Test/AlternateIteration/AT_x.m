function [y] = AT_x(wavelet,x)
% input:
% x          is the reflection coefficient to be inverted (n*m)
% wavelet    is the seismic wavelet constituting the Toeplitz wavelet matrix A

% output:
% y          is the convolution result of wavelet and reflection
%            coefficient, that is, the super gather of seismic records (nm*1)

wavelet_flip=fliplr(wavelet);
[n,m]=size(x);
y=zeros(n,m);
for i=1:m
    y(:,i)=conv(x(:,i),wavelet_flip(:,i),'same');
end

% y=conv2(wavelet_flip,1,x,'same');