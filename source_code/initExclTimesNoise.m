function [initialised, excludedTimes] = initExclTimesNoise(handles)

excludedTimes = [];

[initialised, startPulse, endPulse] = initPulseNoise(handles);
if ~initialised
    return
end

[initialised, startGlitch, endGlitch] = initGlitchNoise(handles);
if ~initialised
    return
end

excludedTimes = struct('startPulse', startPulse, 'endPulse', endPulse, 'startGlitch', startGlitch, 'endGlitch', endGlitch);