function [initialised, excludedTimes] = initExclTimesTarget(handles)

excludedTimes = [];

[initialised, startPulse, endPulse] = initPulseTarget(handles);
if ~initialised
    return
end

[initialised, startGlitch, endGlitch] = initGlitchTarget(handles);
if ~initialised
    return
end

excludedTimes = struct('startPulse', startPulse, 'endPulse', endPulse, 'startGlitch', startGlitch, 'endGlitch', endGlitch);