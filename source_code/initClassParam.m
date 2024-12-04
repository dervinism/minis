function classificationParameters = initClassParam(handles, Ampupbound)

% Amplitude arrays:
ampStepSize = .01;
amplitudeArray = round(100*(0: ampStepSize :Ampupbound))/100;
amplitudeArrayExt = [amplitudeArray amplitudeArray(end) + ampStepSize];

% Rise time arrays:
if isfield(handles, 'figure1')
  value = get(handles.RTbinSizeEdit, 'value');
  RTstepSize = get(handles.RTbinSizeEdit,'string');
  RTstepSize = str2double(strtrim(RTstepSize(value)));
else
  RTstepSize = str2double(handles.RTbinSizeEdit);
end
riseTimeArray = round(100*(0: RTstepSize :11))/100;
riseTimeArrayExt = [riseTimeArray riseTimeArray(end) + RTstepSize];

classificationParameters = struct('amplitudeArray', amplitudeArray, 'amplitudeArrayExt', amplitudeArrayExt, 'riseTimeArray', riseTimeArray,...
    'riseTimeArrayExt', riseTimeArrayExt);