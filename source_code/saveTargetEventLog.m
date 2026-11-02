function ld = saveTargetEventLog(minisArray, incidenceRates, curationFig, resultsPanel, RTinterval, graphicsFormats, ld)
% SAVETARGETEVENTLOG saves the Detect task's output: the detected-event
% data (TSV), the interactive event-summary/curation figure, and the
% combined results panel (histograms, 2D histogram, frequency spectrum,
% and average waveform/rise-time distribution, whichever were computed).
%
% INCIDENCERATES is a per-event vector (same row count as minisArray) --
% see incidenceRate.m -- appended as the TSV's final column.

% Save the detected event data (TSV):
button = questdlg('Save the event data (TSV)?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    [eventFilename, eventPathname, filterIndex] = uiputfile({'*.tsv','TSV files (*.tsv)'},'Save Event Log as', ld);
    if filterIndex
        ld = eventPathname;
        headers = [eventLogHeaders(RTinterval, true, true, size(minisArray,2) == 29), {'Incidence rate (events/s)'}];
        writeEventsTSV(fullfile(eventPathname,eventFilename), headers, [minisArray, incidenceRates(:)]);
    end
end

% Save the detected-event summary graph (the interactive curation figure):
button = questdlg('Save the graph showing detected events?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats,'Save Event Graph as', ld);
    if filterIndex
        ld = figurePathname;
        figFullName = fullfile(figurePathname, figureFilename);
        if ishghandle(curationFig)
            saveas(curationFig, figFullName);
        else
            disp('Warning: detected event figure handle is no longer valid; skipping figure save.');
        end
    end
end

% Save the combined results panel:
if ~isempty(resultsPanel) && ishghandle(resultsPanel)
    button = questdlg('Save the results figure (histograms, spectrum, waveform)?','Save File','Yes','No','Yes');
    if strcmpi(button, 'Yes')
        [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats,'Save Results Figure as', ld);
        if filterIndex
            ld = figurePathname;
            figFullName = fullfile(figurePathname, figureFilename);
            saveas(resultsPanel, figFullName);
        end
    end
end
