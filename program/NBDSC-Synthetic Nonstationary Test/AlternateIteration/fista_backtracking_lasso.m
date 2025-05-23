function [x_s,convergence] = fista_backtracking_lasso(wavelet,b,x0,L0,eta,lambda,iter,opt,eps)
% Solves the following problem via [F]ISTA:
% minimize 1/2*|| Ax - b ||_2^2 + \lambda || x ||_1

% INPUT
%=======================================
% A           n*n
% b           n*m 
% x0......... initial point
% L0 ........ initial choice of stepsize
% eta ....... the constant in which the stepsize is multiplied
% lambda .... coefficient of l1-norm
% iter ...... iteration number
% opt ....... 0 for ISTA or 1 for FISTA
% eps ....... stop criterion
% OUTPUT
%=======================================
% x_s ....... sequences {xk} generated by [F]ISTA

%% 
Ax=@(x) A_x(wavelet,x);

%% 
ATx=@(x) AT_x(wavelet,x);

%% f = 1/2*|| Ax - b ||_2^2
f = @(x) 0.5 * norm(Ax(x)-b,'fro')^2;

%% g = \lambda || x ||_1
g = @(x) lambda * sum(sum(abs(x)));

%% the gradient of f
grad = @(x) ATx((Ax(x)-b));

%% computer F
F = @(x) 0.5*(norm(Ax(x)-b,'fro'))^2 + lambda*sum(sum(abs(x)));

%% shrinkage operator
S = @(tau, g) max(0, g - tau) + min(0, g + tau);

%% projection
P = @(L, y) S(lambda/L, y - (1/L)*grad(y));

%% computer Q
Q = @(L, x, y) f(y) + sum(sum((x-y).*grad(y))) + 0.5*L*norm(x-y,'fro') + g(x);



% x_s = [];
convergence=[];
x_old = x0;
y_old = x0;
L_new = L0;
t_old = 1;
%% MAIN LOOP
for ii = 1:iter
    % find i_k
    j = 1;
    while true
        L_bar = eta^j * L_new;
        if F(P(L_bar, y_old)) <= Q(L_bar, P(L_bar, y_old), y_old)
            L_new = L_new * eta^j;
            break
        else
            j = j + 1;
        end
    end
    x_new = P(L_new, y_old);
    t_new = 0.5 * (1+sqrt(1+4*t_old^2));
    del = opt * (t_old-1)/(t_new);
    y_new = x_new + del*(x_new-x_old);
    % record x_s
%     x_s = [x_s, x_new];
    % check stop criteria
    % e = norm(x_new-x_old,1)/numel(x_new);
    % if e < eps
        % break
    % end
    % update
    x_old = x_new;
    t_old = t_new;
    y_old = y_new;
    convergence=[convergence;F(x_new)];
end
x_s=x_new;
