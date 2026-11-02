function ld = saveTargetEventLogReduced(minisArray, incidenceRates, filename, RTinterval)
% INCIDENCERATES is a per-event vector (see incidenceRate.m), appended as
% the TSV's final column.

% Save the detected event data (TSV):
ld = fileparts(filename);
cd(ld);
headers = [eventLogHeaders(RTinterval, true, true, size(minisArray,2) == 29), {'Incidence rate (events/s)'}];
writeEventsTSV(filename, headers, [minisArray, incidenceRates(:)]);
