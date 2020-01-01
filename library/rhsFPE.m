function [dp_dt] = rhsFPE(timeStamp, pdfData, xdata, driftfcn, diffusionfcn, timeScale)
%RHSFPE This function calculates the RHS of the Fokker-Planck equation
%during solution using ode solver in MATLAB
%
%   driftfcn: Function handle that returns the vector theta for every
%   spatial point defined by xdata, at the time point given by timeStamp
%
%   diffusionfcn: Function handle that returns the scalar diffusion
%   parameter for every spatial point defined by xdata, at the time point
%   given by timeStamp
%
%   pdfData: the pdf at every spatial point defined by xdata, at the time
%   point timeStamp. With the exception of the initial value, this should
%   be calculated by the ode solver in MATLAB
driftArray = generateDriftTerm(driftfcn, pdfData, xdata, timeStamp);        % Calculate the drift flux
diffArray = generateDiffusionTerm(diffusionfcn, pdfData, xdata, timeStamp); % Calculate the diffusion flux
flowArray = driftArray + diffArray;                                         % Calculate total flux
tmpArray = zeros(size(xdata));                                              % Initialize temporary array for dp_dt
indEnd = length(tmpArray);
for ind=1:indEnd
    if (flowArray(ind) > 0)
        tmpArray(ind) = tmpArray(ind) - flowArray(ind);
        if (ind < indEnd)
            tmpArray(ind+1) = tmpArray(ind+1) + flowArray(ind);
        end
    else
        if (flowArray(ind) < 0)
            tmpArray(ind) = tmpArray(ind) + flowArray(ind);
            if (ind > 1)
                tmpArray(ind-1) = tmpArray(ind-1) - flowArray(ind);
            end
        end
    end
end
dp_dt = timeScale .* tmpArray;                                             % Calculate the RHS of the Fokker-Planck equation
end