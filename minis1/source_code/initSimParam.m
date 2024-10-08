function [initialised, simulationParameters] = initSimParam(handles, pulseDuration, riseTimeArray, excludedTimes, filename)

initialised = false;
simulationParameters = [];

if isfield(handles, 'figure1')
  loSimAmp = str2double(get(handles.loSimAmpEdit,'String'));
  if isnan(loSimAmp)
    msgbox('Error: Amplitude lower bound is NaN', 'Error', 'Error');
    return
  end
  
  L = str2double(get(handles.LEdit,'string'));
  if isnan(L)
    msgbox('Error: Initial L (electrotonic length) is NaN','Error','Error');
    return
  end
  
  tau_m = get(handles.tau_mEdit,'string');
  if ~strcmpi(tau_m, '...')
    tau_m = str2double(tau_m);
    if isnan(tau_m)
      msgbox('Error: tau_m (passive membrane time constant; lower limit) is NaN','Error','Error');
      return
    end
  end
  
  if handles.options.tauRange
    tau_PSPm = get(handles.tau_PSPmEdit,'string');
    if ~strcmpi(tau_PSPm, '...')
      tau_PSPm = str2double(tau_PSPm);
      if isnan(tau_PSPm)
        msgbox('Error: tau_m (passive membrane time constant; upper limit) is NaN','Error','Error');
        return
      end
    else
      tau_PSPm = [];
    end
  else
    tau_PSPm = [];
  end
  
else
  loSimAmp = str2double(handles.loSimAmpEdit);
  if isnan(loSimAmp)
    msgbox('Error: Amplitude lower bound is NaN', 'Error', 'Error');
    return
  end
  
  L = str2double(handles.LEdit);
  if isnan(L)
    msgbox('Error: Initial L (electrotonic length) is NaN','Error','Error');
    return
  end
  
  tau_m = handles.tau_mEdit;
  if ~strcmpi(tau_m, '...')
    tau_m = str2double(tau_m);
    if isnan(tau_m)
      msgbox('Error: tau_m (passive membrane time constant; lower limit) is NaN','Error','Error');
      return
    end
  end
  
  if handles.options.tauRange
    tau_PSPm = handles.tau_PSPmEdit;
    if ~strcmpi(tau_PSPm, '...')
      tau_PSPm = str2double(tau_PSPm);
      if isnan(tau_PSPm)
        msgbox('Error: tau_m (passive membrane time constant; upper limit) is NaN','Error','Error');
        return
      end
    else
      tau_PSPm = [];
    end
  else
    tau_PSPm = [];
  end
end

C_m = 1;
R_i = 200;
d = 1E-4;
Q = 2000E-12;
nSeries = 100;

% Estimate tau_m lower bound:
if ~isempty(filename)
  if isfield(handles.options, 'estimateTauLo') && handles.options.estimateTauLo
    button = 'Yes';
  else
    if ~strcmpi(tau_m, '...')
      button = questdlg('Estimate tau_m based on impulses?','Estimate Parameter','Yes','No','No');
    else
      button = 'Yes';
    end
  end
  if strcmpi(button, 'Yes')
    waveform = initWaveform(handles, pulseDuration, riseTimeArray, excludedTimes, filename);
    tau_m = waveform.tau_m;
  end
end

% Estimate tau_m upper bound:
if ~isempty(filename)
  if handles.options.tauRange
    if isfield(handles.options, 'estimateTauHi') && handles.options.estimateTauHi
      button = 'Yes';
    else
      if ~isempty(tau_PSPm)
        button = questdlg('Do you want to estimate the upper limit on tau_m using postsynaptic potentials?','Estimate Parameter','Yes','No','Yes');
      else
        button = 'Yes';
      end
    end
    if strcmpi(button, 'Yes')
      tau_PSPm = [];
    end
  end
end

simulationParameters = struct('tau_m', tau_m, 'tau_PSPm', tau_PSPm, 'C_m', C_m, 'R_i', R_i, 'nSeries', nSeries, 'd', d, 'Q', Q, 'L', L,...
  'loSimAmp', loSimAmp);
initialised = true;