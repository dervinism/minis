function [initialised, startGlitch, endGlitch] = initGlitchNoise(handles)

initialised = false;

errmsgNaN = 'Error: At least one of the noise glitch vectors is NaN';
errmsgVectorLengths = 'Error: Noise glitch vectors are of unequal lengths';
errmsgStartEnd = 'Error: At least one of the noise glitch end vector elements is smaller than its corresponding start glitch element';

try
  startGlitchNoiseStr = get(handles.startGlitchNoiseEdit, 'string');
catch
  startGlitchNoiseStr = handles.startGlitchNoiseEdit;
end
try
  endGlitchNoiseStr = get(handles.endGlitchNoiseEdit, 'string');
catch
  endGlitchNoiseStr = handles.endGlitchNoiseEdit;
end
if ~strcmpi(startGlitchNoiseStr, '...') && ~strcmpi(endGlitchNoiseStr, '...')
    try
        [startGlitch, endGlitch] = testExcludedTimes(startGlitchNoiseStr, endGlitchNoiseStr, errmsgNaN, errmsgVectorLengths, errmsgStartEnd);
    catch %#ok<CTCH>
        return
    end
else
    startGlitch = [];
    endGlitch = [];
end

initialised = true;