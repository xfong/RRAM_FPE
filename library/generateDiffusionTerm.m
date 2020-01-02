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
dx = xdata(2) - xdata(1);
ATmp = gradient(pdfData,dx);
diffusionTerm = - ATmp .* diffusionArray;
end