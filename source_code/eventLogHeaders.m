function headers = eventLogHeaders(RTinterval, includeSD, includePseudoNoiseSD, includeWaveform)
% EVENTLOGHEADERS returns the TSV column headers for a detected-event
% table, matching the column layout produced by detectMinis /
% detectMinisStandalone. The 22 basic detection columns are always
% included; the optional trailing column groups are selected explicitly
% (rather than inferred from a column count, which is ambiguous: a
% 25-column table means "22 basic + 3 SD" in some callers and "22 basic +
% 3 waveform" in others):
%   INCLUDESD           - append 'total SD', '15ms average SD', '30ms
%                          average SD' (3 columns).
%   INCLUDEPSEUDONOISESD - append 'pseudo noise SD' (1 column). Only
%                          meaningful together with INCLUDESD.
%   INCLUDEWAVEFORM      - append 'mean amplitude', 'median amplitude',
%                          'tau' (3 columns) -- the average-waveform
%                          statistics computed over all detected events.

if nargin < 3 || isempty(includePseudoNoiseSD)
    includePseudoNoiseSD = false;
end
if nargin < 4 || isempty(includeWaveform)
    includeWaveform = false;
end

if strcmpi(RTinterval, '10-90%')
    rt1 = '10%'; rt2 = '90%'; rtRange = '10-90% RT';
elseif strcmpi(RTinterval, '20-80%')
    rt1 = '20%'; rt2 = '80%'; rtRange = '20-80% RT';
else
    error('eventLogHeaders:RTinterval', 'Unrecognised RTinterval: %s', RTinterval);
end

headers = {'Peak potential', 'Peak Time', 'Peak Index', 'Amplitude', 'Baseline (BL)', 'BL start time', 'BL end time',...
    'BL start index', 'BL end index', 'Rise time (RT) length', 'RT', rtRange, [rt1 ' RT time mark'], '50% RT time mark',...
    [rt2 ' RT time mark'], [rt1 ' RT index'], '50% RT index', [rt2 ' RT index'], [rt1 ' RT potential'], '50% RT potential',...
    [rt2 ' RT potential'], 'decay time'};

if includeSD
    headers = [headers, {'total SD', '15ms average SD', '30ms average SD'}];
    if includePseudoNoiseSD
        headers = [headers, {'pseudo noise SD'}];
    end
end
if includeWaveform
    headers = [headers, {'mean amplitude', 'median amplitude', 'tau'}];
end
end
