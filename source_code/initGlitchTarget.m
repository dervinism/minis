function [initialised, startGlitch, endGlitch] = initGlitchTarget(handles)

initialised = false;

errmsgNaN = 'Error: At least one of the target glitch vectors is NaN';
errmsgVectorLengths = 'Error: Target glitch vectors are of unequal lengths';
errmsgStartEnd = 'Error: At least one of the target glitch end vector elements is smaller than its corresponding start glitch element';

if isfield(handles, 'figure1')
  startGlitchTargetStr = get(handles.startGlitchTargetEdit, 'string');
  endGlitchTargetStr = get(handles.endGlitchTargetEdit, 'string');
else
  startGlitchTargetStr = handles.startGlitchTargetEdit;
  endGlitchTargetStr = handles.endGlitchTargetEdit;
end
if ~strcmpi(startGlitchTargetStr, '...') && ~strcmpi(endGlitchTargetStr, '...')
    try
        [startGlitch, endGlitch] = testExcludedTimes(startGlitchTargetStr, endGlitchTargetStr, errmsgNaN, errmsgVectorLengths, errmsgStartEnd);
    catch %#ok<CTCH>
        return
    end
else
    startGlitch = [];
    endGlitch = [];
end

initialised = true;