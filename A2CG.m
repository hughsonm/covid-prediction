function [ xmin, nfuncevals ] = A2CG(f_handle,x0,method,grad_tol,kmax)
%A2CG Solves a supplied function using conjugate-gradient technique

% Methods are:
% Fletcher-Reeves
% Polak-Ribiere
% Hestenes-Steifel

% Pick x0

nfuncevals = 0;
N = length(x0);
gnew = f_handle(x0,1);
nfuncevals = nfuncevals + 1;
gold = gnew;
sold = gold;
%snew = sold;
xi = x0;
kk = 1;
while((norm(gnew)>grad_tol) && (kk < kmax))
    %clc;disp(norm(gnew));
    kk = kk+1;
    for ii = 1:N
        if(ii==1)
            snew = -gnew;
        else
            switch(method)
                case 0
                    % Fletcher-Reeves
                    gamma = (gnew'*gnew)/(gold'*gold);
                case 1
                    % Hestenes-Stiefel
                    gamma = ((gnew-gold)'*gnew)/((gnew-gold)'*snew);
                case 2
                    % Polak-Ribiere
                    gamma = ((gnew-gold)'*gnew)/(gold'*gold);
                otherwise
                    % Default to Fletcher-Reeves
                    gamma = (gnew'*gnew)/(gold'*gold);
            end
            snew = -gnew + gamma*sold;
        end
        % Minimize in that direction
        [xi,~,fev] = A1P4SwannGoldSearch(...
            f_handle,...
            snew/norm(snew),...
            xi,...
            1e-4,...
            0.000001,...
            2);
        nfuncevals = nfuncevals + fev;
        gold = gnew;
        gnew = f_handle(xi,1);
        nfuncevals = nfuncevals + 1;
        sold = snew;
        if(norm(gnew)<grad_tol)
            break;
        end
    end
end
xmin = xi;
end

