function ld = saveEventLog(minisTarget, minisNoise, targetIncidenceRates, noiseIncidenceRates, curationFigs, resultsPanel, RTinterval, graphicsFormats, ld)
% SAVEEVENTLOG saves the Detect-and-compare task's output: the target and
% noise event data (TSV), the two interactive event-summary/curation
% figures, and the combined results panel (amplitude/rise-time histograms,
% 2D histograms, frequency spectra, and the target average
% waveform/rise-time distribution, whichever were computed).
%
% CURATIONFIGS is a 2-element array: [targetCurationFig, noiseCurationFig].
% TARGETINCIDENCERATES/NOISEINCIDENCERATES are per-event vectors (see
% incidenceRate.m), appended as each TSV's final column.

% Save the event data files (TSV):
button = questdlg('Save the event data (TSV)?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')

    % Target events:
    [eventFilename, eventPathname, filterIndex] = uiputfile({'*.tsv','TSV files (*.tsv)'},'Save Target Event Log as', ld);
    if filterIndex
        ld = eventPathname;
        targetHeaders = [eventLogHeaders(RTinterval, true, true, size(minisTarget,2) == 29), {'Incidence rate (events/s)'}];
        writeEventsTSV(fullfile(eventPathname,eventFilename), targetHeaders, [minisTarget, targetIncidenceRates(:)]);
    end

    % Noise events (never includes waveform statistics; the noise pass is
    % run with waveform estimation disabled -- see compareMinis.m):
    [eventFilename, eventPathname, filterIndex] = uiputfile({'*.tsv','TSV files (*.tsv)'},'Save Noise Event Log as', ld);
    if filterIndex
        ld = eventPathname;
        noiseHeaders = [eventLogHeaders(RTinterval, true, false, false), {'Incidence rate (events/s)'}];
        writeEventsTSV(fullfile(eventPathname,eventFilename), noiseHeaders, [minisNoise, noiseIncidenceRates(:)]);
    end
end


% Save detected event (curation) figure files:
button = questdlg('Save the graphs showing detected events?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    % Target events:
    [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats, 'Save Target Event Graph as', ld);
    if filterIndex
        ld = figurePathname;
        figFullName = fullfile(figurePathname, figureFilename);
        if ishghandle(curationFigs(1))
            saveas(curationFigs(1), figFullName);
        else
            disp('Warning: target detected-event figure handle is no longer valid; skipping figure save.');
        end
    end
    % Noise events:
    [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats,'Save Noise Event Graph as', ld);
    if filterIndex
        ld = figurePathname;
        figFullName = fullfile(figurePathname, figureFilename);
        if ishghandle(curationFigs(2))
            saveas(curationFigs(2), figFullName);
        else
            disp('Warning: noise detected-event figure handle is no longer valid; skipping figure save.');
        end
    end
end


% Save the combined results panel:
if ~isempty(resultsPanel) && ishghandle(resultsPanel)
    button = questdlg('Save the results figure (histograms, spectra, waveform)?','Save File','Yes','No','Yes');
    if strcmpi(button, 'Yes')
        [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats, 'Save Results Figure as', ld);
        if filterIndex
            ld = figurePathname;
            figFullName = fullfile(figurePathname, figureFilename);
            saveas(resultsPanel, figFullName);
        end
    end
end
