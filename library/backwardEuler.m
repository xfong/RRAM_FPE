function [ outMat, errMsg ] = backwardEuler( dF, dt, t0, y0, ypred, atol, rtol )
%BACKWARDEULER Calculate backward Euler step for solving ODE
%   Assumes a predicted point is given

tn = dt+t0;
y1p = ypred;
y1 = y0 + feval(dF, tn, y1p).*dt;
Niter = 1;
while 1
    if Niter > 100
        break;
    end
    if isrow(y1p)
        ytol = max([abs(y1p); abs(y1)]).*rtol + atol;
    else
        ytol = max([abs(y1p.'); abs(y1.')]).*rtol + atol;
    end
    yerr = abs(y1p - y1);
    chk = (yerr <= ytol.');
    if min(chk) > 0 % Converged
        outMat = y1;
        errMsg = [0; Niter; max(yerr)];
        return
    end
    % If not converged...
    y1p = y1;
    y1 = y0 + feval(dF, tn, y1p).*dt;
    Niter = Niter + 1; % Increment iteration count
end
% If never converge
outMat = y1;
errMsg = [-1; Niter; max(yerr)];
end