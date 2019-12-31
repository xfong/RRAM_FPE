function [dp_dt] = rhsFPE(timeStamp, pdfData, xdata, driftfcn, diffusionfcn)
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
driftArray = generateDriftTerm(driftfcn, pdfData, xdata, timeStamp);        % Calculate the drift current flow
diffArray = generateDiffusionTerm(diffusionfcn, pdfData, xdata, timeStamp); % Calculate the diffusion current flow
flowArray = driftArray + diffArray;                                         % Calculate total current flow
flowSpline = createSplineFromData(flowArray, xdata);                        % Generate spline to fit total current flow
ppfit = fnder(flowSpline);                                                  % Calculate the divergence of the total current flow
dp_dt = -1.0 .* ppval(ppfit, xdata);                                        % Calculate the RHS of the Fokker-Planck equation
end