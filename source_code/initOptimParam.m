function [initialised, optimisationParameters] = initOptimParam(handles)

initialised = false;
optimisationParameters = [];

if isfield(handles, 'figure1')
  distTypeStr = get(handles.distTypeEdit, 'string');
  distTypeValue = get(handles.distTypeEdit, 'value');
  distType = strtrim(distTypeStr(distTypeValue));
  
  distBaselineStr = get(handles.distBaselineEdit, 'string');
  distBaselineValue = get(handles.distBaselineEdit, 'value');
  distBaseline = strtrim(distBaselineStr(distBaselineValue));
  
  options = handles.options;
  
  if handles.options.SDupbound
    SDupbound = str2double(get(handles.SDupbound, 'string'));
    if isnan(SDupbound)
      msgbox('Error: Standard deviation upper bound is NaN', 'Error', 'Error');
      return
    end
  else
    SDupbound = [];
  end
  
  if handles.options.SDlobound
    SDlobound = str2double(get(handles.SDlobound, 'string'));
    if isnan(SDlobound)
      msgbox('Error: Standard deviation lower bound is NaN', 'Error', 'Error');
      return
    end
  else
    SDlobound = [];
  end
  
  maxSAD = str2double(get(handles.maxErrEdit, 'string'));
  if isnan(maxSAD)
    msgbox('Error: Maximum combined SAD is NaN','Error','Error');
    return
  end
  
  maxAmpSAD = str2double(get(handles.maxAmpErrEdit, 'string'));
  if isnan(maxAmpSAD)
    msgbox('Error: Maximum amplitude SAD is NaN','Error','Error');
    return
  end
  
  maxRTSAD = str2double(get(handles.maxRTErrEdit, 'string'));
  if isnan(maxRTSAD)
    msgbox('Error: Maximum rise time SAD is NaN','Error','Error');
    return
  end
  
  maxMAD = str2double(get(handles.maxDevEdit, 'string'));
  if isnan(maxMAD)
    msgbox('Error: Maximum combined MAD is NaN','Error','Error');
    return
  end
  
  maxAmpMAD = str2double(get(handles.maxDevAmpEdit, 'string'));
  if isnan(maxAmpMAD)
    msgbox('Error: Maximum amplitude MAD is NaN','Error','Error');
    return
  end
  
  maxRTMAD = str2double(get(handles.maxDevRTEdit, 'string'));
  if isnan(maxRTMAD)
    msgbox('Error: Maximum rise time MAD is NaN','Error','Error');
    return
  end
  
  maxAmpBottomSAD = str2double(get(handles.maxAmpBottomErrEdit, 'string'));
  if isnan(maxAmpBottomSAD)
    msgbox('Error: Maximum amplitude top 50% SAD is NaN','Error','Error');
    return
  end
  
  maxAmpBottomMAD = str2double(get(handles.maxAmpBottomDevEdit, 'string'));
  if isnan(maxAmpBottomMAD)
    msgbox('Error: Maximum amplitude top 50% MAD is NaN','Error','Error');
    return
  end
  
  maxAmpMidSAD = str2double(get(handles.maxAmpMidErrEdit, 'string'));
  if isnan(maxAmpMidSAD)
    msgbox('Error: Maximum amplitude top 10% SAD is NaN','Error','Error');
    return
  end
  
  maxAmpMidMAD = str2double(get(handles.maxAmpMidDevEdit, 'string'));
  if isnan(maxAmpMidMAD)
    msgbox('Error: Maximum amplitude top 10% MAD is NaN','Error','Error');
    return
  end
  
  maxAmpTopSAD = str2double(get(handles.maxAmpTopErrEdit, 'string'));
  if isnan(maxAmpTopSAD)
    msgbox('Error: Maximum amplitude top 2% SAD is NaN','Error','Error');
    return
  end
  
  maxAmpTopMAD = str2double(get(handles.maxAmpTopDevEdit, 'string'));
  if isnan(maxAmpTopMAD)
    msgbox('Error: Maximum amplitude top 2% MAD is NaN','Error','Error');
    return
  end
  
else
  distType = strtrim(handles.distTypeEdit);
  
  distBaseline = strtrim(handles.distBaselineEdit);
  
  options = handles.options;
  
  if handles.options.SDupbound
    SDupbound = str2double(handles.SDupbound);
    if isnan(SDupbound)
      msgbox('Error: Standard deviation upper bound is NaN', 'Error', 'Error');
      return
    end
  else
    SDupbound = [];
  end
  
  if handles.options.SDlobound
    SDlobound = str2double(handles.SDlobound);
    if isnan(SDlobound)
      msgbox('Error: Standard deviation lower bound is NaN', 'Error', 'Error');
      return
    end
  else
    SDlobound = [];
  end
  
  maxSAD = str2double(handles.maxErrEdit);
  if isnan(maxSAD)
    msgbox('Error: Maximum combined SAD is NaN','Error','Error');
    return
  end
  
  maxAmpSAD = str2double(handles.maxAmpErrEdit);
  if isnan(maxAmpSAD)
    msgbox('Error: Maximum amplitude SAD is NaN','Error','Error');
    return
  end
  
  maxRTSAD = str2double(handles.maxRTErrEdit);
  if isnan(maxRTSAD)
    msgbox('Error: Maximum rise time SAD is NaN','Error','Error');
    return
  end
  
  maxMAD = str2double(handles.maxDevEdit);
  if isnan(maxMAD)
    msgbox('Error: Maximum combined MAD is NaN','Error','Error');
    return
  end
  
  maxAmpMAD = str2double(handles.maxDevAmpEdit);
  if isnan(maxAmpMAD)
    msgbox('Error: Maximum amplitude MAD is NaN','Error','Error');
    return
  end
  
  maxRTMAD = str2double(handles.maxDevRTEdit);
  if isnan(maxRTMAD)
    msgbox('Error: Maximum rise time MAD is NaN','Error','Error');
    return
  end
  
  maxAmpBottomSAD = str2double(handles.maxAmpBottomErrEdit);
  if isnan(maxAmpBottomSAD)
    msgbox('Error: Maximum amplitude top 50% SAD is NaN','Error','Error');
    return
  end
  
  maxAmpBottomMAD = str2double(handles.maxAmpBottomDevEdit);
  if isnan(maxAmpBottomMAD)
    msgbox('Error: Maximum amplitude top 50% MAD is NaN','Error','Error');
    return
  end
  
  maxAmpMidSAD = str2double(handles.maxAmpMidErrEdit);
  if isnan(maxAmpMidSAD)
    msgbox('Error: Maximum amplitude top 10% SAD is NaN','Error','Error');
    return
  end
  
  maxAmpMidMAD = str2double(handles.maxAmpMidDevEdit);
  if isnan(maxAmpMidMAD)
    msgbox('Error: Maximum amplitude top 10% MAD is NaN','Error','Error');
    return
  end
  
  maxAmpTopSAD = str2double(handles.maxAmpTopErrEdit);
  if isnan(maxAmpTopSAD)
    msgbox('Error: Maximum amplitude top 2% SAD is NaN','Error','Error');
    return
  end
  
  maxAmpTopMAD = str2double(handles.maxAmpTopDevEdit);
  if isnan(maxAmpTopMAD)
    msgbox('Error: Maximum amplitude top 2% MAD is NaN','Error','Error');
    return
  end
end

optimisationParameters = struct('distType', distType, 'distBaseline', distBaseline, 'SDupbound', SDupbound, 'SDlobound', SDlobound, 'SAD', maxSAD,...
  'AmpSAD', maxAmpSAD, 'RTSAD', maxRTSAD, 'MAD', maxMAD, 'AmpMAD', maxAmpMAD, 'RTMAD', maxRTMAD, 'AmpBottomSAD', maxAmpBottomSAD, 'AmpBottomMAD',...
  maxAmpBottomMAD, 'AmpMidSAD', maxAmpMidSAD, 'AmpMidMAD', maxAmpMidMAD, 'AmpTopSAD', maxAmpTopSAD, 'AmpTopMAD', maxAmpTopMAD, 'options', options);
initialised = true;