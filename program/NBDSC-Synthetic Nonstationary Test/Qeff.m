function [q_eff]=Qeff(q,t)
 
n=length(q);
dt=t(2)-t(1);
q_eff=q(1);
for i=1:n-1
    eff=sum(1./q(1:i)*dt)/t(i+1);
    q_eff=[q_eff;1/eff];    
end