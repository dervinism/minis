function waveform = initWaveform(handles, pulseDuration, riseTimeArrayExt, excludedTimes, filename)

waveform.estimate = 1;
waveform.riseTimeArrayExt = riseTimeArrayExt;

% Estimate tau_m:
if length(excludedTimes.startPulse) == 2 || length(excludedTimes.startPulse) == 1
    if isnan(pulseDuration)
        options.Resize = 'on';
        options.WindowStyle = 'normal';
        pulseDuration = inputdlg('Enter exact pulse duration (ms):','Current Pulse',1,{'0.5'},options);
        if ~isempty(pulseDuration)
            pulseDuration = cell2mat(pulseDuration);
            if isnan(pulseDuration)
                msgbox('Error: Zero pulse duration', 'Error', 'Error');
                return
            end
        else
            msgbox('Error: Zero pulse duration', 'Error', 'Error');
            return
        end
    end
    if length(excludedTimes.startPulse) == 2
        tau_m = doublePeel(pulseDuration, filename, false, (excludedTimes.startPulse(2)-.1)*1000, (excludedTimes.endPulse(2))*1000);
    elseif length(excludedTimes.startPulse) == 1
        tau_m = doublePeel(pulseDuration, filename, false, (excludedTimes.startPulse(1)-.1)*1000, (excludedTimes.endPulse(1))*1000);
    end
    if ~isreal(tau_m)
        if isfield(handles, 'figure1')
            tau_m = str2double(get(handles.tau_mEdit,'string'));
        else
            tau_m = str2double(handles.tau_mEdit);
        end
        if isnan(tau_m)
            tau_m = 10;
        end
    end
else
    if isfield(handles, 'figure1')
        tau_m = str2double(get(handles.tau_mEdit,'string'));
    else
        tau_m = str2double(handles.tau_mEdit);
    end
    if isnan(tau_m)
        tau_m = 10;
    end
end
waveform.tau_m = tau_m;