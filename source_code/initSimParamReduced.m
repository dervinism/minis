function [initialised, simulationParameters] = initSimParamReduced(handles)

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
  tau_m = str2double(tau_m);
  if isnan(tau_m)
    msgbox('Error: tau_m (passive membrane time constant; lower limit) is NaN','Error','Error');
    return
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
  tau_m = str2double(tau_m);
  if isnan(tau_m)
    msgbox('Error: tau_m (passive membrane time constant; lower limit) is NaN','Error','Error');
    return
  end
end

tau_PSPm = tau_m;

C_m = 1;
R_i = 200;
d = 1E-4;
Q = 2000E-12;
nSeries = 100;

simulationParameters = struct('tau_m', tau_m, 'tau_PSPm', tau_PSPm, 'C_m', C_m, 'R_i', R_i, 'nSeries', nSeries, 'd', d, 'Q', Q, 'L', L,...
  'loSimAmp', loSimAmp);
initialised = true;