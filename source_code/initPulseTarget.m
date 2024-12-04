function [initialised, startPulse, endPulse] = initPulseTarget(handles)

initialised = false;

errmsgNaN = 'Error: At least one of the target pulse vectors is NaN';
errmsgVectorLengths = 'Error: Target pulse vectors are of unequal lengths';
errmsgStartEnd = 'Error: At least one of the target pulse end vector elements is smaller than its corresponding start pulse element';

if isfield(handles, 'figure1')
  startPulseTargetStr = get(handles.startPulseTargetEdit, 'string');
  endPulseTargetStr = get(handles.endPulseTargetEdit, 'string');
else
  startPulseTargetStr = handles.startPulseTargetEdit;
  endPulseTargetStr = handles.endPulseTargetEdit;
end
if ~strcmpi(startPulseTargetStr, '...') && ~strcmpi(endPulseTargetStr, '...')
    try
        [startPulse, endPulse] = testExcludedTimes(startPulseTargetStr, endPulseTargetStr, errmsgNaN, errmsgVectorLengths, errmsgStartEnd);
    catch %#ok<CTCH>
        return
    end
else
    startPulse = [];
    endPulse = [];
end

initialised = true;