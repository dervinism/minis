function excludedTimes = calcExcludedTimes(sweepDuration, numSweeps, startPulse, endPulse, startGlitch, endGlitch, dt)
% excludedTimes is a helper subfunction of the optimiseMinis function. It
% calculates a time vector with values that should be excluded from further
% analysis.
%
% Calling syntax:
%   excludedTimes = calcExcludedTimes(sweepDuration, numSweeps,...
%   startPulse, endPulse, startGlitch, endGlitch, dt)
%
% The input variables are:
%   sweepDuration - a scalar value of recording sweep duration in
%       milliseconds (ms);
%   numSweeps - a total number of recording sweeps;
%   startPulse - a vector containing starting times of different pulses, ms;
%   endPulse - a vector containing ending times of different pulses, ms;
%   startGlitch - a vector containing starting times of various glitches, ms;
%   endGlitch - a vector containing ending times of various glitches, ms;
%   dt - a sampling interval, ms.
%



if sweepDuration > 0
    % Mark the beginning/end of each sweep for exclusion:
    if numSweeps >= 1
        numSweeps = round(numSweeps);
        sweepExcludedT = [(sweepDuration:sweepDuration:numSweeps*sweepDuration) - sweepDuration numSweeps*sweepDuration];
    else
        sweepExcludedT = [];
    end
    
    % Exclude pulses:
    if ~isempty(startPulse) && ~isempty(endPulse) && max(endPulse) > 0 && length(startPulse) == length(endPulse)
        pulseExcludedT = zeros(1,ceil(max(endPulse) + sweepDuration*(numSweeps+1)));
        for iSweep = 1:numSweeps
            for iPulse = 1:length(startPulse)
                if iSweep == 1 && iPulse == 1
                    pulse = (startPulse(iPulse) + dt: dt :endPulse(iPulse));
                    pulseExcludedT(1:length(pulse)) = pulse;
                else
                    lastEntry = find(pulseExcludedT,1,'last');
                    pulse = (startPulse(iPulse) + sweepDuration*(iSweep-1) + dt: dt :endPulse(iPulse) + sweepDuration*(iSweep-1));
                    pulseExcludedT(lastEntry + 1: lastEntry+length(pulse)) = pulse;
                end
            end
            lastEntry = find(pulseExcludedT,1,'last');
            pulseExcludedT(lastEntry+1:end) = [];
        end
    else
        pulseExcludedT = [];
    end
    
    % Exclude glitches:
    if ~isempty(startGlitch) && ~ isempty(endGlitch) && max(endGlitch) > 0 && length(startGlitch) == length(endGlitch)
        glitchExcludedT = zeros(1,ceil(max(endGlitch)));
        for iGlitch = 1:length(startGlitch)
            if iGlitch == 1
                glitch = (startGlitch(iGlitch) + dt: dt :endGlitch(iGlitch));
                glitchExcludedT(1:length(glitch)) = glitch;
            else
                lastEntry = find(glitchExcludedT,1,'last');
                glitch = (startGlitch(iGlitch) + dt: dt :endGlitch(iGlitch));
                glitchExcludedT(lastEntry + 1: lastEntry + length(glitch)) = glitch;
            end
        end
        lastEntry = find(glitchExcludedT,1,'last');
        glitchExcludedT(lastEntry+1:end) = [];
    else
        glitchExcludedT = [];
    end
    
    % Concatenate the excluded times:
    initExcludedT = [sweepExcludedT pulseExcludedT glitchExcludedT];
    uniqueExcludedT = unique(initExcludedT);
    excludedTimes = sort(uniqueExcludedT);
else
    excludedTimes = [];
end
end