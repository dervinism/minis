function [proceed, startLim, endLim] = thrMinis(succeed, threshold, startLim, endLim, dt)

if succeed
    startLim = threshold;
    if threshold < endLim
        proceed = false;
    else
        proceed = true;
    end
else
    if threshold > startLim + dt
        endLim = threshold;
        proceed = false;
    elseif threshold > startLim
        endLim = startLim;
        proceed = false;
    else
        proceed = true;
    end
end