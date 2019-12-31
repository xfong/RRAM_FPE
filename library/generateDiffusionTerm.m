function [diffusionTerm] = generateDiffusionTerm(diffusionfcn,pdfData,xdata, timeStamp)
%GENERATEDIFFUSIONTERM This function determines the diffusion term that
%goes into the calculation of the flow of the PDF
%
%   diffusionfcn: Function handle that returns the scalar diffusion
%   parameter for every spatial point defined by xdata, at the time point
%   given by timeStamp
%
%   pdfData: the pdf at every spatial point defined by xdata, at the time
%   point timeStamp. With the exception of the initial value, this should
%   be calculated by the ode solver in MATLAB
diffusionArray = feval(diffusionfcn, xdata, timeStamp);    % Determine diffusion parameter for every point in space (defined in xdata) at a particular time point
diffSpline = createSplineFromData(diffusionArray, xdata);  % Generate spline for diffusion parameter
pdfSpline = createSplineFromData(pdfData, xdata);          % Generate spline for pdf at current time point
diffusionprime = fnder(diffSpline);                        % Generate function for spatial derivative of diffusion parameter
pdfprime = fnder(pdfSpline);                               % Generate function for spatial derivative of pdf parameter
% Calculate the diffusion term using chain rule for -grad(D.p);
ATmp = ppval(diffusionprime,xdata);
ATerm = ATmp.*pdfData;
BTmp = ppval(pdfprime,xdata);
BTerm = BTmp.*diffusionArray;
diffusionTerm = - ATerm - BTerm;
end