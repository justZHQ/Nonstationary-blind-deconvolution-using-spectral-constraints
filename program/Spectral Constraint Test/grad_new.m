function [loss,x]=grad_new(F,b,x0,x_object,lambda1,lambda2,e,iter,k)
%% f = 1/2*lambda1*|| |Fx|^2-|b|^2 ||_2^2+1/2*lambda2*|| |x|^2 ||_2^2
m=length(b);
n=length(x0);
FR=real(F);
FI=imag(F);

gradi=@(x,FRi,FIi,bi) (x'*FRi'*FRi*x+x'*FIi'*FIi*x-bi^2)*(FRi'*FRi+FIi'*FIi)*x;

xk=x0;
x=[xk];

nk=1;

loss=[];
for j=1:iter
    grad=zeros(n,1);
    for i=1:m
        grad=grad+gradi(xk,FR(i,:),FI(i,:),b(i));
    end
    grad=lambda1*grad+lambda2*(xk-x_object);
    xk=xk-e*grad;
    if j==k(nk)
    x=[x,xk];
    if nk<length(k)
        nk=nk+1;
    end
    end
    loss=[loss norm(abs(F*xk)-b)]
end
