function [gradall]=gradFft(lambda0,b,FR,FI,x)
[n,m]=size(x);
f_number=length(b(:,1));
% [~,m]=size(b);

% i=1;
% j=2;
% xj=x(:,j);
% FRi=FR(i,:);
% FIi=FI(i,:);
% bi=b(i,j);
% a=(xj'*FRi'*FRi*xj+xj'*FIi'*FIi*xj-bi^2);
% ab=(FRi'*FRi+FIi'*FIi)*xj;
% c=a*ab;

gradi=@(x,FRi,FIi,bi) (x'*FRi'*FRi*x+x'*FIi'*FIi*x-bi^2)*(FRi'*FRi+FIi'*FIi)*x;
gradall=zeros(n,m);
if lambda0~=0
    for j=1:m
        grad=zeros(n,1);
        for i=1:f_number
%             grad=grad+gradi(x(:,j),FR(:,i),FI(:,i),b(i,j));
            grad=grad+gradi(x(:,j),FR(i,:),FI(i,:),b(i,j));
        end
        gradall(:,j)=grad;
    end
end