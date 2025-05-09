function [y] = A_x(wavelet,x)
% input:
% x          is the reflection coefficient to be inverted (n*m)
% wavelet    is the seismic wavelet constituting the Toeplitz wavelet matrix A

% output:
% y          is the convolution result of wavelet and reflection
%            coefficient, that is, the super gather of seismic records (nm*1)

[n,m]=size(x);
y=zeros(n,m);
for i=1:m
    y(:,i)=conv(x(:,i),wavelet(:,i),'same');
end


% y=conv2(wavelet,1,x,'same');
