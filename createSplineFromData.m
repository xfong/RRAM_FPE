function [pp] = createSplineFromData(ydata,xdata)
%createSplineFromData Function takes as input ydata and xdata, and creates
%a spline that fits the data points
%   The function error checks the input arrays before calling spline to
%   generate the fitting spline

if (isrow(ydata))
    if (isrow(xdata))
        if (length(ydata) ~= length(xdata))
            error('Number of x- and y- data points do not match!');
        else
            pp = spline(xdata, ydata);
        end
    else
        if (length(ydata) ~= length(xdata))
            error('Number of x- and y- data points do not match!');
        else
            pp = spline(xdata, ydata.');
        end        
    end
else
    if (iscolumn(ydata))
        if (iscolumn(xdata))
            if (length(ydata) ~= length(xdata))
                error('Number of x- and y- data points do not match!');
            else
                pp = spline(xdata, ydata);
            end
        else
            if (length(ydata) ~= length(xdata))
                error('Number of x- and y- data points do not match!');
            else
                pp = spline(xdata, ydata.');
            end
        end
    else
        error('Each of ydata and xdata must be either a row or a column vector!');
    end
end

end