function [initialised, detectionParameters] = initDetectParam(handles)

initialised = false;
detectionParameters = [];

if isfield(handles, 'figure1')
  maxTimeToPeak = str2double(get(handles.maxTimeToPeakEdit, 'string'));
  if isnan(maxTimeToPeak)
    msgbox('Error: Maximum time to peak is NaN', 'Error', 'Error');
    return
  end
  
  baselineDuration = str2double(get(handles.baselineDurationEdit,'string'));
  if isnan(baselineDuration)
    msgbox('Error: Baseline duration is NaN', 'Error', 'Error');
    return
  end
  
  peakIntegrationPeriod = str2double(get(handles.peakIntegrationPeriodEdit,'string'));
  if isnan(peakIntegrationPeriod)
    msgbox('Error: Peak integration period is NaN', 'Error', 'Error');
    return
  end
  
  Amplobound = str2double(get(handles.AmploboundEdit,'string'));
  if isnan(Amplobound)
    msgbox('Error: Amplitude lower bound is NaN', 'Error', 'Error');
    return
  end
  
  Ampupbound = str2double(get(handles.AmpupboundEdit,'string'));
  if isnan(Ampupbound)
    msgbox('Error: Amplitude upper bound is NaN', 'Error', 'Error');
    return
  end
  
  smoothWindow = str2double(get(handles.smoothWindowEdit,'string'));
  if isnan(smoothWindow)
    msgbox('Error: Gaussian smoothing window is NaN', 'Error', 'Error');
    return
  end
  
  pulseDuration = get(handles.pulseDurationEdit,'string');
  if ~strcmpi(pulseDuration, '...')
    pulseDuration = str2double(pulseDuration);
    if isnan(pulseDuration)
      msgbox('Error: Brief (second) pulse duration is NaN', 'Error', 'Error');
      return
    end
  end
  
  value = get(handles.RTintEdit,'value');
  RTinterval = get(handles.RTintEdit,'string');
  RTinterval = strtrim(RTinterval(value));
  downGoing = get(handles.downGoingCheckbox,'Value');
  voltageClamp = get(handles.voltageClampCheckbox,'Value');
  
else
  maxTimeToPeak = str2double(handles.maxTimeToPeakEdit);
  if isnan(maxTimeToPeak)
    msgbox('Error: Maximum time to peak is NaN', 'Error', 'Error');
    return
  end
  
  baselineDuration = str2double(handles.baselineDurationEdit);
  if isnan(baselineDuration)
    msgbox('Error: Baseline duration is NaN', 'Error', 'Error');
    return
  end
  
  peakIntegrationPeriod = str2double(handles.peakIntegrationPeriodEdit);
  if isnan(peakIntegrationPeriod)
    msgbox('Error: Peak integration period is NaN', 'Error', 'Error');
    return
  end
  
  Amplobound = str2double(handles.AmploboundEdit);
  if isnan(Amplobound)
    msgbox('Error: Amplitude lower bound is NaN', 'Error', 'Error');
    return
  end
  
  Ampupbound = str2double(handles.AmpupboundEdit);
  if isnan(Ampupbound)
    msgbox('Error: Amplitude upper bound is NaN', 'Error', 'Error');
    return
  end
  
  smoothWindow = str2double(handles.smoothWindowEdit);
  if isnan(smoothWindow)
    msgbox('Error: Gaussian smoothing window is NaN', 'Error', 'Error');
    return
  end
  
  pulseDuration = handles.pulseDurationEdit;
  if ~strcmpi(pulseDuration, '...')
    pulseDuration = str2double(pulseDuration);
    if isnan(pulseDuration)
      msgbox('Error: Brief (second) pulse duration is NaN', 'Error', 'Error');
      return
    end
  end
  
  RTinterval = handles.RTintEdit;
  if RTinterval == 1
    RTinterval = '10-90%';
  elseif RTinterval == 2
    RTinterval = '20-80%';
  end
  downGoing = handles.downGoingCheckbox;
  voltageClamp = handles.voltageClampCheckbox;
end

initialised = true;
detectionParameters = struct('SWstart', maxTimeToPeak, 'BLduration', baselineDuration, 'refractoryPeriod', peakIntegrationPeriod,...
  'Amplobound', Amplobound, 'Ampupbound', Ampupbound, 'smoothWindow', smoothWindow, 'smoothWindowLite', 8, 'pulseDuration', pulseDuration, 'RTinterval', RTinterval, 'downGoing',downGoing, 'voltageClamp',voltageClamp);