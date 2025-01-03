function ld = saveTargetEventLogReduced(minisArray, filename, RTinterval)

% Save the detected event text file:
ld = fileparts(filename);
cd(ld);
fid = fopen(filename,'wt+');
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