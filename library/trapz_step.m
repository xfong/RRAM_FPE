function [ outMat, errMsg ] = trapz_step( dF, dt, t, y0, ypred, atol, rtol )
%TRAPZ_STEP Take a step in trapz

tn = t+dt;
NIter = 1;
y1p = ypred;
dy0 = dt.*feval(dF, t, y0);
slope1 = feval(dF, tn, y1p);
y1 = y0 + 0.5.*(dy0 + dt.*slope1);
%% Begin iterations
while 1
    if (NIter > 100)
        break;
    end
    if isrow(y1p)
        ytol = max(max([abs(y1p); abs(y1)]).*rtol, atol);
    else
        ytol = max(max([abs(y1p.'); abs(y1.')]).*rtol, atol);
    end
    yerr = abs(y1p - y1);
    if min(yerr <= ytol.') > 0 % Converged
        outMat = y1;
        errMsg = [0; NIter; max(yerr)];
        return
    end
    % If not converged...
    y1p = y1;
    dy1 = feval(dF, tn, y1p).*dt;
    y1 = y0 + 0.5.*(dy0 + dy1);
    NIter = NIter + 1; % Increment iteration count
end
outMat = y1;
errMsg = [-1; NIter; max(yerr)];
end