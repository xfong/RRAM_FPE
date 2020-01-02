function [ tpoints, y_out, errEst, errMsg ] = abm23_stepper( dF, tspan, y0, dtinit, atol, rtol )
%STEPPER Perform numerical time solution using AB2 predictor and
%trapezoidal corrector
%   Detailed explanation goes here

% Determine how we step through time points based on size of tspan array
if (length(tspan) == 1) % Given only end time point
    tstart = 0;
    tend = tspan;
    if (dtinit > 0)
        dt0 = dtinit;
    else
        dt0 = tend / 100;
    end
    FixDt = 0;
elseif (length(tspan) == 2) % Given start and end time points
    tstart = tspan(1);
    tend = tspan(2);
    if tstart == 0
        if (dtinit > 0)
            dt0 = dtinit;
        else
            dt0 = tend / 100;
        end
    else
        if (dtinit > 0)
            dt0 = dtinit;
        else
            dt0 = tstart / 100;
        end
    end
    FixDt = 0;
elseif (length(tspan) == 3) % Given start and end time points and interval between time points
    tstart = tspan(1);
    tend = tspan(2);
    dt0 = tspan(3);
    FixDt = 1;
end

tfinish = tend;
if tstart > 0
    tfinish = tstart;
elseif tstart < 0
    error('Time interval must start at 0!');
end
tarray = 0;
tcurr = 0;
tnext = dt0;
if isrow(y0)
    y0 = y0.';
end
ytmp = zeros(length(y0),2);
ytmp(:,1) = y0;
yearr = 0;
stepCount = 1;
delt_curr = dt0;
delt_prev_1 = dt0;
y_prev_1 = y0;
slope_prev_1 = y0;
delt_prev_2 = dt0;
y_prev_2 = y0;
convFlag = 0;
while tcurr < tfinish % Loop until time hits tfinish
    if convFlag ~= 0 % When converged at a particular time step
        convFlag = 0; % Reset flag
        stepCount = stepCount + 1; % Increment step count
         % Shift previously computed points one step back
        delt_prev_2 = delt_prev_1;
        delt_prev_1 = delt_curr;
        delt_curr = delt_next;
        % Update stamps
        tarray(stepCount) = tnext;
        tcurr = tnext;
        tnext = tnext + delt_curr;
        ytmp(:,stepCount) = ycorr;
        yearr(stepCount) = max(yerr);
        slope_prev_1 = slope_curr;
        if tcurr == tfinish
            break;
        end
    end
    if tnext > tfinish
        delt_curr = delt_curr - (tnext - tfinish);
        tnext = tfinish;
    end
    slope_curr = feval(dF, tcurr, ytmp(:,stepCount));
    if stepCount < 2 % When not enough steps to execute AB2 method...
        ypred = ytmp(:,stepCount) + delt_curr .* slope_curr;
        [ycorr, errRet] = backwardEuler(dF, delt_curr, tcurr, ytmp(:,stepCount), ypred, atol, rtol);
        if errRet(1) < 0 % Failed to cnverge to solution within step
            if FixDt ~= 0 % If not allowed to vary time step...
                error('Failed to converge with fixed dt!');
            else
                delt_curr = delt_curr ./ 8;
                continue;
            end
        end
        if isrow(ypred) && isrow(ycorr)
            maxy = max([abs(ypred);abs(ycorr)]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred);
            ychk = yerr ./ ytol;
        elseif iscolumn(ypred) && iscolumn(ycorr)
            maxy = max([abs(ypred) abs(ycorr)].');
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred).';
            ychk = yerr ./ ytol;
        elseif isrow(ypred) && iscolumn(ycorr)
            maxy = max([abs(ypred);abs(ycorr.')]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr.' - ypred);
            ychk = yerr ./ ytol;
        elseif iscolumn(ypred) && isrow(ycorr)
            maxy = max([abs(ypred.');abs(ycorr)]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred.');
            ychk = yerr ./ ytol;
        else
            error('ypred or ycorr are not vectors!');
        end
        facinv = max(ychk);
        if facinv == 0
            fac = 2;
        elseif facinv < 0.25
            fac = 0.96.*exp(log(facinv) ./ -2); % Be a little conservative when increasing time step
        else
            fac = exp(log(facinv) ./ -2);
        end
        if ( (facinv > 1) && (FixDt == 0) ) % When we should reduce time step, capture as convergence failed
            convFlag = 0;
            delt_curr = fac .* delt_curr;
            tnext = tcurr + delt_curr;
        else
            convFlag = 1;
            if FixDt ~= 0
                fac = 1;
            end
            delt_next = fac * delt_curr;
        end
    else % Enough steps to execute AB2 method
        ypred = ytmp(:,stepCount) + (0.5.*delt_curr./delt_prev_1).*(((delt_curr + 2 .* delt_prev_1).*slope_curr) - delt_curr .* slope_prev_1);
        [ycorr, errRet] = trapz_step(dF, delt_curr, tcurr, ytmp(:,stepCount), ypred, atol, rtol);
        if errRet(1) < 0 % Failed to converge to solution within step
            if FixDt ~= 0 % If not allowed to vary time step...
                error('Failed to converge with fixed dt!');
            else
                delt_curr = delt_curr ./ 8; % Immediately cut delta
                continue;
            end
        end
        if isrow(ypred) && isrow(ycorr)
            maxy = max([abs(ypred);abs(ycorr)]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred);
            ychk = yerr ./ ytol;
        elseif iscolumn(ypred) && iscolumn(ycorr)
            maxy = max([abs(ypred) abs(ycorr)].');
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred).';
            ychk = yerr ./ ytol;
        elseif isrow(ypred) && iscolumn(ycorr)
            maxy = max([abs(ypred);abs(ycorr.')]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr.' - ypred);
            ychk = yerr ./ ytol;
        elseif iscolumn(ypred) && isrow(ycorr)
            maxy = max([abs(ypred.');abs(ycorr)]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred.');
            ychk = yerr ./ ytol;
        else
            error('ypred or ycorr are not vectors!');
        end
        facinv = max(ychk);
        if facinv == 0
            fac = 2;
        elseif facinv < 0.125
            fac = 0.96.*exp(log(facinv) ./ -3); % Be a little conservative when increasing time step
        else
            fac = exp(log(facinv) ./ -3);
        end
        if ( (facinv > 1) && (FixDt == 0) ) % When we should reduce time step, capture as convergence failed
            convFlag = 0;
            delt_curr = fac .* delt_curr ./ 2;
            tnext = tcurr + delt_curr;
        else
            convFlag = 1;
            if FixDt ~= 0
                fac = 1;
            end
            delt_next = fac .* delt_curr;
        end
    end
end

if tstart == 0
    tpoints = tarray;
    y_out = ytmp;
    errEst = yearr;
    errMsg = 0;
    return;
end
tfinish = tend;
if stepCount >= 2;
    startInd = stepCount - 1;
    ytmp = ytmp(:, startInd:stepCount);
    tarray = tarray(startInd:stepCount);
    yearr = yearr(startInd:end);
    stepCount = 2;
end
convFlag = 0;
while tcurr < tfinish % Loop until time hits tfinish
    if convFlag ~= 0 % When converged at a particular time step
        convFlag = 0; % Reset flag
        stepCount = stepCount + 1; % Increment step count
         % Shift previously computed points one step back
        delt_prev_2 = delt_prev_1;
        delt_prev_1 = delt_curr;
        delt_curr = delt_next;
        % Update stamps
        tarray(stepCount) = tnext;
        tcurr = tnext;
        tnext = tnext + delt_curr;
        ytmp(:,stepCount) = ycorr;
        yearr(stepCount) = max(yerr);
        slope_prev_1 = slope_curr;
        if tcurr == tfinish
            break;
        end
    end
    if tnext > tfinish
        delt_curr = delt_curr - (tnext - tfinish);
        tnext = tfinish;
    end
    slope_curr = feval(dF, tcurr, ytmp(:,stepCount));
    if stepCount < 2 % When not enough steps to execute AB2 method...
        ypred = ytmp(:,stepCount) + delt_curr .* slope_curr;
        [ycorr, errRet] = backwardEuler(dF, delt_curr, tcurr, ytmp(:,stepCount), ypred, atol, rtol);
        if errRet(1) < 0 % Failed to cnverge to solution within step
            if FixDt ~= 0 % If not allowed to vary time step...
                error('Failed to converge with fixed dt!');
            else
                delt_curr = delt_curr ./ 8;
                continue;
            end
        end
        if isrow(ypred) && isrow(ycorr)
            maxy = max([abs(ypred);abs(ycorr)]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred);
            ychk = yerr ./ ytol;
        elseif iscolumn(ypred) && iscolumn(ycorr)
            maxy = max([abs(ypred) abs(ycorr)].');
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred).';
            ychk = yerr ./ ytol;
        elseif isrow(ypred) && iscolumn(ycorr)
            maxy = max([abs(ypred);abs(ycorr.')]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr.' - ypred);
            ychk = yerr ./ ytol;
        elseif iscolumn(ypred) && isrow(ycorr)
            maxy = max([abs(ypred.');abs(ycorr)]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred.');
            ychk = yerr ./ ytol;
        else
            error('ypred or ycorr are not vectors!');
        end
        facinv = max(ychk);
        if facinv == 0
            fac = 2;
        elseif facinv < 0.25
            fac = 0.96.*exp(log(facinv) ./ -2); % Be a little conservative when increasing time step
        else
            fac = exp(log(facinv) ./ -2);
        end
        if ( (facinv > 1) && (FixDt == 0) ) % When we should reduce time step, capture as convergence failed
            convFlag = 0;
            delt_curr = fac .* delt_curr;
            tnext = tcurr + delt_curr;
        else
            convFlag = 1;
            if FixDt ~= 0
                fac = 1;
            end
            delt_next = fac * delt_curr;
        end
    else % Enough steps to execute AB2 method
        ypred = ytmp(:,stepCount) + (0.5.*delt_curr./delt_prev_1).*(((delt_curr + 2 .* delt_prev_1).*slope_curr) - delt_curr .* slope_prev_1);
        [ycorr, errRet] = trapz_step(dF, delt_curr, tcurr, ytmp(:,stepCount), ypred, atol, rtol);
        if errRet(1) < 0 % Failed to converge to solution within step
            if FixDt ~= 0 % If not allowed to vary time step...
                error('Failed to converge with fixed dt!');
            else
                delt_curr = delt_curr ./ 8; % Immediately cut delta
                continue;
            end
        end
        if isrow(ypred) && isrow(ycorr)
            maxy = max([abs(ypred);abs(ycorr)]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred);
            ychk = yerr ./ ytol;
        elseif iscolumn(ypred) && iscolumn(ycorr)
            maxy = max([abs(ypred) abs(ycorr)].');
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred).';
            ychk = yerr ./ ytol;
        elseif isrow(ypred) && iscolumn(ycorr)
            maxy = max([abs(ypred);abs(ycorr.')]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr.' - ypred);
            ychk = yerr ./ ytol;
        elseif iscolumn(ypred) && isrow(ycorr)
            maxy = max([abs(ypred.');abs(ycorr)]);
            ytol = max(maxy.*rtol, atol);
            yerr = abs(ycorr - ypred.');
            ychk = yerr ./ ytol;
        else
            error('ypred or ycorr are not vectors!');
        end
        facinv = max(ychk);
        if facinv == 0
            fac = 2;
        elseif facinv < 0.125
            fac = 0.96.*exp(log(facinv) ./ -3); % Be a little conservative when increasing time step
        else
            fac = exp(log(facinv) ./ -3);
        end
        if ( (facinv > 1) && (FixDt == 0) ) % When we should reduce time step, capture as convergence failed
            convFlag = 0;
            delt_curr = fac .* delt_curr ./ 2;
            tnext = tcurr + delt_curr;
        else
            convFlag = 1;
            if FixDt ~= 0
                fac = 1;
            end
            delt_next = fac .* delt_curr;
        end
    end
end
tpoints = tarray(2:end);
y_out = ytmp(:,2:end);
errEst = yearr(2:end);
errMsg = 0;
% % if FixDt % If interval between time points is fixed
% %     if tstart > 0 % If we only want data after a perticular time point
% %         tpoints = tstart:dt0:tend; % For output
% %         tintern = 0:dt0:tstart; % Calculate internally till the time point to start from
% %     elseif tstart == 0 % If start time point is zero, we do not discard any data
% %         tintern = 0:dt0:tend;
% %     else % Error check the time point input
% %         error('tstart cannot be negative!');
% %     end
% %     y_intern = zeros(length(y0),length(tintern)); % Temporary buffer
% %     % Loading initial condition...
% %     if isrow(y0) % Make sure the array is a column vector
% %         y0 = y0.';
% %     elseif ~iscol(y0)
% %         error('Initial conditions must be a row or column vector!');
% %     end
% %     y_intern(:,1) = y0; % Results for corrector method
% %     y_intern_alt = y_intern; % Results for predictor method
% %     slope_n_min1 = zeros(size(y0(:,1))); % Temporary buffer
% %     slope_n = feval(@slope_func,0,y0(:,1)); % Temporary buffer
% %     for i0 = 2:length(tintern) % First round stepping through time points
% %         % Moving computed slopes according to time points
% %         slope_n_min2 = slope_n_min1;
% %         slope_n_min1 = slope_n;
% %         if i0 <= 3 % Insufficient time points have been generated for AB2 predictor
% %             % Use backward Euler predictor first
% %             [y_intern(:,i0), errMsg] = trapz_step(@slope_func, dt, tintern(i0), y_intern(:,i0-1), atol, rtol);
% %             slope_n = feval(@slope_func, tintern(i0), y_intern(:,i0));
% %             if errMsg(2) >= Niter
% %                 Niter = errMsg(2);
% %                 delmax = errMsg(3);
% %                 occTime = i0;
% %             end
% %             if errMsg(1) < 0
% %                 printf('Non-convergence detected at step %d\n', i0-1);
% %             end
% %             y_intern_alt(:,i0)=y_intern(:,i0);
% %         else
% %             y_intern_alt(:,i0) = y_intern_alt(:,i0-1)+(0.5.*dt0).*(3.*slope_n_min1 - slope_n_min2);
% %             [y_intern(:,i0), errMsg] = trapz_step(@slope_func, dt, tintern(i0), y_intern_alt(:,i0-1), atol, rtol);
% %             if errMsg(2) >= Niter
% %                 Niter = errMsg(2);
% %                 delmax = errMsg(3);
% %                 occTime = i0;
% %             end
% %             if errMsg(1) < 0
% %                 printf('Non-convergence detected at step %d\n', i0-1);
% %             end
% %             slope_n = feval(@slope_func, tintern(i0), y_intern(:,i0));
% %         end
% %     end
% %     if tstart == 0
% %         tpoints = tintern;
% %         y_out = y_intern;
% %         errEst = abs(y_intern_alt - y_intern);
% %         errMsg = 0;
% %         return
% %     end
% %     ytmp = y_intern(:,end-1:end);
% %     y_intern = zeros(length(ytmp),length(tpoints)+1);
% %     y_intern(:,1:2) = ytmp;
% %     y_intern_alt = y_intern;
% %     for i0 = 2:length(tpoints)
% %         slope_n_min2 = slope_n_min1;
% %         slope_n_min1 = slope_n;
% %         y_intern_alt(:,i0) = y_intern_alt(:,i0-1)+(0.5.*dt0).*(3.*slope_n_min1 - slope_n_min2);
% %         [y_intern(:,i0), errMsg] = trapz_step(@slope_func, dt, tpoints(i0), y_intern_alt(:,i0-1), atol, rtol);
% %         if errMsg(2) >= Niter
% %             Niter = errMsg(2);
% %             delmax = errMsg(3);
% %             occTime = i0;
% %         end
% %         if errMsg(1) < 0
% %             printf('Non-convergence detected at step %d\n', i0-1);
% %         end
% %         slope_n = feval(@slope_func, tintern(i0), y_intern(:,i0));
% %     end
% %     y_out = y_intern(:,2:end);
% %     errEst = abs(y_intern_alt(:,2:end)-y_intern(:,2:end));
% %     errMsg = 0;
% % else % When we do not have fixed time steps
% %     if tstart > 0 % Find out when we need to start recording results
% %         t_term = tstart;
% %     else
% %         t_term = tend;
% %     end
% %     tsave = zeros(1,1);
% %     ysave = zeros(size(y0));
% %     ysave(:,1) = y0;
% %     tstep = dt0;
% %     step_iter = 1;
% %     tcurr = tstep;
% %     while t_curr < t_term
% %         if step_iter < 3 % When we have insufficient points to obtain predictor using AB2 method
% %             % Half time step, and do one forward Euler steps
% %             
% %         else
% %         end
% %         Terminate if we hit time point
% %         if tcurr >= t_term
% %             break;
% %         end
% %         % Setup for next step at end of step
% %         step_iter = step_iter + 1;
% %         ysave(:,end+1) = yconv;
% %         tnext = tcurr + tstep;
% %         tsave(end+1) = tnext;
% %         tcurr = tnext;
% %     end
% % end

end