function [driftTerm] = generateDriftTerm(driftfcn,pdfData,xdata,timeStamp)
%GENERATEDRIFTTERM This function determines the drift term that goes into
%the calculation of the flow of the PDF
%
%   driftfcn: Function handle that returns the vector theta for every
%   spatial point defined by xdata, at the time point given by timeStamp
%
%   pdfData: the pdf at every spatial point defined by xdata, at the time
%   point timeStamp. With the exception of the initial value, this should
%   be calculated by the ode solver in MATLAB
driftArray = feval(driftfcn, xdata, timeStamp);  % Determine drift parameter at every point in space (defined by xdata) for current timeStamp
driftTerm = driftArray .* pdfData;               % Calculate the drift term
end