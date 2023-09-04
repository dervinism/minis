function ld = saveTargetEventLog(minisArray, F, waveform, filtering, RTinterval, graphicsFormats, ld)

% Save the detected event text file:
button = questdlg('Save the event text file?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    [eventFilename, eventPathname, filterIndex] = uiputfile({'*.txt','Text files (*.txt)'},'Save Event Log as', ld);
    if filterIndex
        ld = eventPathname;
        fid = fopen(fullfile(eventPathname,eventFilename),'wt+');
        if strcmpi(RTinterval, '10-90%')
            if size(minisArray,2) == 29
                fprintf(fid, '%16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\n',...
                    'Peak potential', 'Peak Time', 'Peak Index', 'Amplitude', 'Baseline (BL)', 'BL start time', 'BL end time', 'BL start index',...
                    'BL end index', 'Rise time (RT) length', 'RT', '10-90% RT', '10% RT time mark', '50% RT time mark', '90% RT time mark',...
                    '10% RT index', '50% RT index', '90% RT index', '10% RT potential', '50% RT potential', '90% RT potential', 'decay time',...
                    'total SD', '15ms average SD', '30ms average SD', 'pseudo noise SD', 'mean top 10%', 'median top 10%', 'tau');
                fprintf(fid,'%16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\n',minisArray');
            else
                fprintf(fid, '%16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\n',...
                    'Peak potential', 'Peak Time', 'Peak Index', 'Amplitude', 'Baseline (BL)', 'BL start time', 'BL end time', 'BL start index',...
                    'BL end index', 'Rise time (RT) length', 'RT', '10-90% RT', '10% RT time mark', '50% RT time mark', '90% RT time mark',...
                    '10% RT index', '50% RT index', '90% RT index', '10% RT potential', '50% RT potential', '90% RT potential', 'decay time',...
                    'total SD', '15ms average SD', '30ms average SD', 'pseudo noise SD');
                fprintf(fid,'%16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\n',minisArray');
            end
        elseif strcmpi(RTinterval, '20-80%')
            if size(minisArray,2) == 29
                fprintf(fid, '%16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\n',...
                    'Peak potential', 'Peak Time', 'Peak Index', 'Amplitude', 'Baseline (BL)', 'BL start time', 'BL end time', 'BL start index',...
                    'BL end index', 'Rise time (RT) length', 'RT', '20-80% RT', '20% RT time mark', '50% RT time mark', '80% RT time mark',...
                    '20% RT index', '50% RT index', '80% RT index', '20% RT potential', '50% RT potential', '80% RT potential', 'decay time',...
                    'total SD', '15ms average SD', '30ms average SD', 'pseudo noise SD', 'mean top 10%', 'median top 10%', 'tau');
                fprintf(fid,'%16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\n',minisArray');
            else
                fprintf(fid, '%16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\n',...
                    'Peak potential', 'Peak Time', 'Peak Index', 'Amplitude', 'Baseline (BL)', 'BL start time', 'BL end time', 'BL start index',...
                    'BL end index', 'Rise time (RT) length', 'RT', '20-80% RT', '20% RT time mark', '50% RT time mark', '80% RT time mark',...
                    '20% RT index', '50% RT index', '80% RT index', '20% RT potential', '50% RT potential', '80% RT potential', 'decay time',...
                    'total SD', '15ms average SD', '30ms average SD', 'pseudo noise SD');
                fprintf(fid,'%16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\n',minisArray');
            end
        end
        fclose(fid);
    end
end

% Save the detected event figure file:
button = questdlg('Save the graph showing detected events?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats,'Save Event Graph as', ld);
    if filterIndex
        ld = figurePathname;
        figFullName = fullfile(figurePathname, figureFilename);
        saveas(F(3), figFullName);
    end
end

% Save event distribution histograms:
button = questdlg('Save the event distribution histograms?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats,'Save Event Histograms as', ld);
    if filterIndex
        ld = figurePathname;
        figFullName = fullfile(figurePathname, figureFilename);
        saveas(F(1), figFullName);
    end
    [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats,'Save Event 2D Histogram as', ld);
    if filterIndex
        ld = figurePathname;
        figFullName = fullfile(figurePathname, figureFilename);
        saveas(F(2), figFullName);
    end
end

% Save the frequency spectrum graph:
if strcmpi(filtering.state,'on')
    button = questdlg('Save the frequency spectrum graph?','Save File','Yes','No','Yes');
    if strcmpi(button, 'Yes')
        [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats,'Save Frequency Spectrum as', ld);
        if filterIndex
            ld = figurePathname;
            figFullName = fullfile(figurePathname, figureFilename);
            saveas(F(4), figFullName);
        end
    end
end

% Save the the average waveform and related rise time distribution graph:
if ~isempty(waveform) && isfield(waveform,'F') && ~isempty(waveform.F)
    if strcmpi(RTinterval, '10-90%')
        qst = 'Save the average waveform graph and a related 10-90% rise time distribution graph?';
        msg = 'Save 10-90% Rise Time Distribution Graph as';
    elseif strcmpi(RTinterval, '20-80%')
        qst = 'Save the average waveform graph and a related 20-80% rise time distribution graph?';
        msg = 'Save 20-80% Rise Time Distribution Graph as';
    end
    button = questdlg(qst,'Save File','Yes','No','Yes');
    if strcmpi(button, 'Yes')
        [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats, 'Save Average Waveform as', ld);
        if filterIndex
            ld = figurePathname;
            figFullName = fullfile(figurePathname, figureFilename);
            saveas(waveform.F(1), figFullName);
        end
        [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats, msg, ld);
        if filterIndex
            ld = figurePathname;
            figFullName = fullfile(figurePathname, figureFilename);
            saveas(waveform.F(2), figFullName);
        end
    end
end