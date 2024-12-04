function [initialised, startPulse, endPulse] = initPulseNoise(handles)

initialised = false;

errmsgNaN = 'Error: At least one of the noise pulse vectors is NaN';
errmsgVectorLengths = 'Error: Noise pulse vectors are of unequal lengths';
errmsgStartEnd = 'Error: At least one of the noise pulse end vector elements is smaller than its corresponding start pulse element';

try
  startPulseNoiseStr = get(handles.startPulseNoiseEdit, 'string');
catch
  startPulseNoiseStr = handles.startPulseNoiseEdit;
end
try
  endPulseNoiseStr = get(handles.endPulseNoiseEdit, 'string');
catch
  endPulseNoiseStr = handles.endPulseNoiseEdit;
end
if ~strcmpi(startPulseNoiseStr, '...') && ~strcmpi(endPulseNoiseStr, '...')
    try
        [startPulse, endPulse] = testExcludedTimes(startPulseNoiseStr, endPulseNoiseStr, errmsgNaN, errmsgVectorLengths, errmsgStartEnd);
    catch %#ok<CTCH>
        return
    end
else
    startPulse = [];
    endPulse = [];
end

initialised = true;