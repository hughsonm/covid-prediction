function [ minx, mingrad, nfuncevals ] = A1P4SwannGoldSearch(...
    f_handle,...
    dd,...
    x0,...
    tolerance,...
    delta,...
    magnification)
%A1P4SwannGoldSearch Finds the minimum of supplied function along direction
%dd
nfuncevals = 0;
%figure(24);
%subplot(2,1,1);
% Error checking
if(length(dd) ~= length(x0)), error('x0 and dd must be same length');end
if(norm(dd) == 0), error('dd must contain a non-zero direction'); end

% Evaluate function at initial point and plot.
f0 = f_handle(x0,0); nfuncevals = nfuncevals + 1;
%hold off;stem(0,f0); hold on;

% Evaluate function at delta along d and plot.
f1 = f_handle(x0 + delta*dd,0);nfuncevals = nfuncevals + 1;
%stem(delta,f1);

ready_for_gold = 0;
search_dir = 1;

if(f1>f0)
    % Try the other direction, this one was no good
    f2 = f_handle(x0-delta*dd,0);nfuncevals = nfuncevals + 1;
    %stem(-delta,f2);
    if(f2>f0)
        % This direction is also no good.
        % Begin golden section search between f1 and f2
        ready_for_gold = 1;
        xl = x0 - delta*dd;
        xh = x0+delta*dd;
    else
        % Okay so we should search in the negative direction
        search_dir = -1;
        dd = -dd;
        xim1 = x0;
        fim1 = f0;
        xi   = x0 + delta*dd;
        fi   = f2;
    end
else
    % We got the direction right, first try.
    xim1 = x0;
    fim1 = f0;
    xi   = x0+ delta*dd;
    fi   = f1;
end

% Perform Swann steps in the direction dd, plotting all function
% evaluations.
cur_step_size = delta;
lambda = delta;
while(~ready_for_gold)
    min_was_found = 0;
    cur_step_size = cur_step_size * magnification;
    cur_step = dd*cur_step_size;
    lambda = lambda + cur_step_size*search_dir;
    xip1 = xi + cur_step;
    fip1 = f_handle(xip1,0);nfuncevals = nfuncevals + 1;
    %stem(lambda,fip1);
    if(fip1 > fi)
        min_was_found = 1;
        % Set up golden section search
        xl = xim1;
        xh = xip1;        
    else
        xim1 = xi;
        xi = xip1;
        fim1 = fi;
        fi = fip1;
    end        
    ready_for_gold = min_was_found;
end

% Prepare for golden section search.
%subplot(2,1,2);
dx = (xh-xl);
PHI = (sqrt(5)-1)/2;

xb = xl + PHI*dx;
xa = xl + (1-PHI)*dx;

interval_size = norm(dx);
interval_bottom = 0;

fa = f_handle(xa,0);nfuncevals = nfuncevals + 1;
fb = f_handle(xb,0);nfuncevals = nfuncevals + 1;
%hold off;stem(interval_size*(1-PHI),fa);hold on;
%stem(interval_size*PHI,fb);

%stem(0,0);
%stem(interval_size,0);

while(interval_size>tolerance)
    % Pick the new interval, perform new function evaluation, and plot.
    if(fa<fb)        
        xh = xb;
        dx = xh-xl;
        xb = xa;
        xa = xl + (1-PHI)*dx;        
        fb = fa;
        fa = f_handle(xa,0);nfuncevals = nfuncevals + 1;   
        interval_size = interval_size*PHI;
        %stem(interval_bottom + interval_size*(1-PHI),fa);
    else
        xl = xa;
        dx = xh-xl;
        xa = xb;
        xb = xl + PHI*dx;
        fa = fb;
        fb = f_handle(xb,0);nfuncevals = nfuncevals + 1;
        interval_bottom = interval_bottom + (1-PHI)*interval_size;
        interval_size = interval_size*PHI;
        %stem(interval_bottom + interval_size*(PHI),fb);
    end
end

% Return x-star
minx = (xl + xh)/2;
mingrad = f_handle(minx,1);

end

