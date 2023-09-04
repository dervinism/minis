function [croppedData, excludedTimes] = cropData(data, sweepDuration, numSweeps, startPulse, endPulse, startGlitch, endGlitch, dt)
% cropData is a helper subfunction of the optimiseMinis function. It
% removes pulse periods and glitches from the data and calculates a time
% vector with values that should be excluded from further analysis of the
% cropped data.
%
% The input variables are:
%   data - data vector (either mV or nA);
%   sweepDuration - a scalar value of recording sweep duration in miliseconds (ms);
%   numSweeps - a total number of recording sweeps;
%   startPulse - a vector containing starting times of different pulses, ms;
%   endPulse - a vector containing ending times of different pulses, ms;
%   startGlitch - a vector containing starting times of various glitches, ms;
%   endGlitch - a vector containing ending times of various glitches, ms;
%   dt - a sampling interval, ms.
%
% The output variables are:
%   croppedData - a vector with cropped data (either mV or nA);
%   excludedT - a vector of excluded times in the croppedData, ms.
%

if sweepDuration > 0
    if numSweeps >= 1
        numSweeps = round(numSweeps);
        sweepExcludedT = (sweepDuration:sweepDuration:numSweeps*sweepDuration);
    else
        sweepExcludedT = [];
    end
    
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
    
    initExcludedT = [pulseExcludedT glitchExcludedT];
    excludedTimes = unique(initExcludedT);
    
    
    % Crop data:
    if ~isempty(pulseExcludedT) || ~isempty(glitchExcludedT)
        cropStartIndex = round(startGlitch/dt);
        croppedLength = 0;
        sweepEndIndices = round(sweepExcludedT/dt);
        croppedSweepEndIndices = zeros(1,length(sweepEndIndices));
        excludedIndices = round(excludedTimes/dt);
        for iSweep = 1:numSweeps
            sweepExcludedIndices = excludedIndices(excludedIndices <= sweepEndIndices(iSweep));
            croppedLength = croppedLength + length(sweepExcludedIndices);
            croppedSweepEndIndices(iSweep) = sweepEndIndices(iSweep) - length(sweepExcludedIndices);
            if iSweep == 1
                cropStartIndex(cropStartIndex <= sweepEndIndices(iSweep)) = cropStartIndex(cropStartIndex <= sweepEndIndices(iSweep)) - croppedLength;
            else
                cropStartIndex(cropStartIndex > sweepEndIndices(iSweep-1) & cropStartIndex <= sweepEndIndices(iSweep)) =...
                    cropStartIndex(cropStartIndex > sweepEndIndices(iSweep-1) & cropStartIndex <= sweepEndIndices(iSweep)) - croppedLength;
            end
        end
        croppedData = data;
        croppedData(excludedIndices) = [];
        initExcludedT = [croppedSweepEndIndices*dt cropStartIndex*dt];
        uniqueExcludedT = unique(initExcludedT);
        excludedTimes = sort(uniqueExcludedT);
    else
        croppedData = data;
        excludedTimes = sweepExcludedT;
    end
else
    excludedTimes = [];
end
end