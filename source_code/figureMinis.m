function path = figureMinis(path)

global ratio initFile %#ok<*NUSED>



%% Load the data:
[dataFilename, dataPathname, filterIndex] = uigetfile({'*.mat','MAT files (*.mat)'}, 'Provide data file for figure retrieval', path);
if filterIndex
    dataFilename = fullfile(dataPathname, dataFilename);
    path = dataPathname;
    load(dataFilename); %#ok<LOAD>
    if ~exist('tHalfF', 'var')
        uiwait(msgbox('The supplied data file is corrupt. Please provide another file.','modal'));
        return
    end
    
    
    
    %% Plot the data:
    % Plot the average and median amplitude value of the top 10% of detected events:
    f2 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    axes('Position', [.04 .07 .92 .82]);
    l1 = line(tSweep,medianAmp10', 'LineStyle',':', 'LineWidth',1.5, 'Color','r');
    ax1 = gca;
    xlabel('Recording time (s)');
    ylabel('mV or nA');
    
    ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
    l2 = line(fSweep,averageAmp10', 'LineStyle',':', 'LineWidth',1.5, 'Color','g', 'Parent',ax2);
    xlabel('Files');
    axisvals1 = get(ax1,'ylim');
    axisvals2 = get(ax2,'ylim');
    yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
    set(ax1,'xlim',[0 tEnd]);
    set(ax2,'xlim',[initFile initFile+length(files)]);
    set(ax1,'ylim',yLimit);
    set(ax2,'ylim',yLimit);
    
    xTimeMarks = [APinfuse APblock minisInfuse];
    yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
    nTimeMarks = length(xTimeMarks);
    yTimeMarks1(nTimeMarks+1 : end) = [];
    yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
    yTimeMarks2(nTimeMarks+1 : end) = [];
    line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    l3 = zeros(1,nTimeMarks);
    for iMark = 1:nTimeMarks
        l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
    end
    
    set(f2, 'NumberTitle', 'off');
    set(f2, 'Name', 'Examine Amplitudes');
    title('Amplitudes of the top 10% of detected events');
    if nTimeMarks
        legend([l1 l2 l3(1)],'Median','Mean','Drug timing', 'Location','NorthEast');
    else
        legend([l1 l2],'Median','Mean', 'Location','NorthEast');
    end
    
    % Create the zoom object for the figure f2:
    figureZoomHandle = zoom;
    set(figureZoomHandle,'ActionPostCallback',@minisZoom);
    
    % Create the pan object for the figure f2:
    figurePanHandle = pan(f2);
    set(figurePanHandle,'ActionPostCallback',@minisPan);
    
    
    % Plot tau_m:
    f3 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    axes('Position', [.04 .07 .92 .82]);
    l1 = line(tHalfF,tau_m, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    ax1 = gca;
    xlabel('Recording time (s)');
    ylabel('Time constant (ms)');
    
    ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
    xlabel('Files');
    axisvals1 = get(ax1,'ylim');
    axisvals2 = get(ax2,'ylim');
    yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
    set(ax1,'xlim',[0 tEnd]);
    set(ax2,'xlim',[initFile initFile+length(files)]);
    set(ax1,'ylim',yLimit);
    set(ax2,'ylim',yLimit);
    
    xTimeMarks = [APinfuse APblock minisInfuse];
    yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
    nTimeMarks = length(xTimeMarks);
    yTimeMarks1(nTimeMarks+1 : end) = [];
    yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
    yTimeMarks2(nTimeMarks+1 : end) = [];
    line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    l2 = zeros(1,nTimeMarks);
    for iMark = 1:nTimeMarks
        l2(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
    end
    
    set(f3, 'NumberTitle', 'off');
    set(f3, 'Name', 'Examine PSP-based Effective Membrane Time Constant');
    title('Effective passive membrane time constant of the top 10% of detected events');
    if nTimeMarks
        legend([l1 l2(1)],'Effective membrane time constant','Drug timing', 'Location','NorthEast');
    else
        legend(l1,'Effective membrane time constant', 'Location','NorthEast');
    end
    
    % Create the zoom object for the figure f3:
    figureZoomHandle = zoom;
    set(figureZoomHandle,'ActionPostCallback',@minisZoom);
    
    % Create the pan object for the figure f3:
    figurePanHandle = pan(f3);
    set(figurePanHandle,'ActionPostCallback',@minisPan);
    
    
    
    % Plot SD (100ms window):
    f4 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    axes('Position', [.04 .07 .92 .82]);
    l1 = line(0.001*tSTD100, STD100, 'Color','b');
    ax1 = gca;
    xlabel('Recording time (s)');
    ylabel('mV or nA');
    
    ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
    l2 = line(swpSTD100,STDsmooth100, 'Color','g', 'Parent',ax2);
    xlabel('Files');
    axisvals1 = get(ax1,'ylim');
    axisvals2 = get(ax2,'ylim');
    yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
    set(ax1,'xlim',[0 tEnd]);
    set(ax2,'xlim',[initFile initFile+length(files)]);
    set(ax1,'ylim',yLimit);
    set(ax2,'ylim',yLimit);
    
    xTimeMarks = [APinfuse APblock minisInfuse];
    yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
    nTimeMarks = length(xTimeMarks);
    yTimeMarks1(nTimeMarks+1 : end) = [];
    yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
    yTimeMarks2(nTimeMarks+1 : end) = [];
    line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    l3 = zeros(1,nTimeMarks);
    for iMark = 1:nTimeMarks
        l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
    end
    
    set(f4, 'NumberTitle', 'off');
    set(f4, 'Name', 'Examine Standard Deviation (100 ms Window)');
    title('Standard deviation of the electrophysiological recording data (100 ms window)');
    if nTimeMarks
        legend([l1 l2 l3(1)],'Raw data','Smoothed data','Drug timing', 'Location','NorthEast');
    else
        legend([l1 l2],'Raw data','Smoothed data', 'Location','NorthEast');
    end
    
    % Create the zoom object for the figure f4:
    figureZoomHandle = zoom;
    set(figureZoomHandle,'ActionPostCallback',@minisZoom);
    
    % Create the pan object for the figure f4:
    figurePanHandle = pan(f4);
    set(figurePanHandle,'ActionPostCallback',@minisPan);
    
    
    
    % Plot SD (15ms window):
    f5 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    axes('Position', [.04 .07 .92 .82]);
    l1 = line(0.001*tSTD15, STD15, 'Color','b');
    ax1 = gca;
    xlabel('Recording time (s)');
    ylabel('mV or nA');
    
    ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
    l2 = line(swpSTD15,STDsmooth15, 'Color','g', 'Parent',ax2);
    xlabel('Files');
    axisvals1 = get(ax1,'ylim');
    axisvals2 = get(ax2,'ylim');
    yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
    set(ax1,'xlim',[0 tEnd]);
    set(ax2,'xlim',[initFile initFile+length(files)]);
    set(ax1,'ylim',yLimit);
    set(ax2,'ylim',yLimit);
    
    xTimeMarks = [APinfuse APblock minisInfuse];
    yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
    nTimeMarks = length(xTimeMarks);
    yTimeMarks1(nTimeMarks+1 : end) = [];
    yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
    yTimeMarks2(nTimeMarks+1 : end) = [];
    line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    l3 = zeros(1,nTimeMarks);
    for iMark = 1:nTimeMarks
        l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
    end
    
    set(f5, 'NumberTitle', 'off');
    set(f5, 'Name', 'Examine Standard Deviation (15 ms Window)');
    title('Standard deviation of the electrophysiological recording data (15 ms window)');
    if nTimeMarks
        legend([l1 l2 l3(1)],'Raw data','Smoothed data','Drug timing', 'Location','NorthEast');
    else
        legend([l1 l2],'Raw data','Smoothed data', 'Location','NorthEast');
    end
    
    % Create the zoom object for the figure f5:
    figureZoomHandle = zoom;
    set(figureZoomHandle,'ActionPostCallback',@minisZoom);
    
    % Create the pan object for the figure f5:
    figurePanHandle = pan(f5);
    set(figurePanHandle,'ActionPostCallback',@minisPan);
    
    
    
    % Plot the baseline:
    f6 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    axes('Position', [.04 .07 .92 .82]);
    l1 = line(0.001*tSTD100, mean100, 'Color','b');
    ax1 = gca;
    xlabel('Recording time (s)');
    ylabel('mV or nA');
    
    ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
    l2 = line(swpSTD100,meanSmooth100, 'Color','g', 'Parent',ax2);
    xlabel('Files');
    axisvals1 = get(ax1,'ylim');
    axisvals2 = get(ax2,'ylim');
    yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
    set(ax1,'xlim',[0 tEnd]);
    set(ax2,'xlim',[initFile initFile+length(files)]);
    set(ax1,'ylim',yLimit);
    set(ax2,'ylim',yLimit);
    
    xTimeMarks = [APinfuse APblock minisInfuse];
    yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
    nTimeMarks = length(xTimeMarks);
    yTimeMarks1(nTimeMarks+1 : end) = [];
    yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
    yTimeMarks2(nTimeMarks+1 : end) = [];
    line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    l3 = zeros(1,nTimeMarks);
    for iMark = 1:nTimeMarks
        l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
    end
    
    set(f6, 'NumberTitle', 'off');
    set(f6, 'Name', 'Examine Baseline Fluctuations (100 ms Window)');
    title('Baseline of the electrophysiological recording data (100 ms window)');
    if nTimeMarks
        legend([l1 l2 l3(1)],'Raw data','Smoothed data','Drug timing', 'Location','NorthEast');
    else
        legend([l1 l2],'Raw data','Smoothed data', 'Location','NorthEast');
    end
    
    % Create the zoom object for the figure f6:
    figureZoomHandle = zoom;
    set(figureZoomHandle,'ActionPostCallback',@minisZoom);
    
    % Create the pan object for the figure f6:
    figurePanHandle = pan(f6);
    set(figurePanHandle,'ActionPostCallback',@minisPan);
    
    
    
    % Plot pulse-based tau_m:
    f7 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    axes('Position', [.04 .07 .92 .82]);
    l1 = line(tHalfF,tauPulseEff, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    ax1 = gca;
    xlabel('Recording time (s)');
    ylabel('Time constant (ms)');
    
    ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
    l2 = line(fHalfF,tauPulse, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g', 'Parent',ax2);
    xlabel('Files');
    axisvals1 = get(ax1,'ylim');
    axisvals2 = get(ax2,'ylim');
    yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
    set(ax1,'xlim',[0 tEnd]);
    set(ax2,'xlim',[initFile initFile+length(files)]);
    set(ax1,'ylim',yLimit);
    set(ax2,'ylim',yLimit);
    
    xTimeMarks = [APinfuse APblock minisInfuse];
    yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
    nTimeMarks = length(xTimeMarks);
    yTimeMarks1(nTimeMarks+1 : end) = [];
    yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
    yTimeMarks2(nTimeMarks+1 : end) = [];
    line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    l3 = zeros(1,nTimeMarks);
    for iMark = 1:nTimeMarks
        l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
    end
    
    set(f7, 'NumberTitle', 'off');
    set(f7, 'Name', 'Examine Impulse-based Passive Membrane Time Constant');
    title('Impulse-based passive membrane time constant');
    if nTimeMarks
        legend([l2 l1 l3(1)],'Membrane time constant','Effective time constant','Drug timing', 'Location','NorthEast');
    else
        legend([l2 l1],'Membrane time constant','Effective time constant', 'Location','NorthEast');
    end
    
    % Create the zoom object for the figure f7:
    figureZoomHandle = zoom;
    set(figureZoomHandle,'ActionPostCallback',@minisZoom);
    
    % Create the pan object for the figure f7:
    figurePanHandle = pan(f7);
    set(figurePanHandle,'ActionPostCallback',@minisPan);
    
    
    
    % Plot capacitance:
    f8 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    axes('Position', [.04 .07 .92 .82]);
    l1 = line(tHalfF,capacitanceEff, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    ax1 = gca;
    xlabel('Recording time (s)');
    ylabel('Capacitance (pF)');
    
    ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
    l2 = line(fHalfF,capacitance, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g', 'Parent',ax2);
    xlabel('Files');
    axisvals1 = get(ax1,'ylim');
    axisvals2 = get(ax2,'ylim');
    yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
    set(ax1,'xlim',[0 tEnd]);
    set(ax2,'xlim',[initFile initFile+length(files)]);
    set(ax1,'ylim',yLimit);
    set(ax2,'ylim',yLimit);
    
    xTimeMarks = [APinfuse APblock minisInfuse];
    yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
    nTimeMarks = length(xTimeMarks);
    yTimeMarks1(nTimeMarks+1 : end) = [];
    yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
    yTimeMarks2(nTimeMarks+1 : end) = [];
    line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    l3 = zeros(1,nTimeMarks);
    for iMark = 1:nTimeMarks
        l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
    end
    
    set(f8, 'NumberTitle', 'off');
    set(f8, 'Name', 'Examine Total Capacitance');
    title('Total capacitance');
    if nTimeMarks
        legend([l2 l1 l3(1)],'Capacitance','Effective capacitance','Drug timing', 'Location','NorthEast');
    else
        legend([l2 l1],'Capacitance','Effective capacitance', 'Location','NorthEast');
    end
    
    % Create the zoom object for the figure f8:
    figureZoomHandle = zoom;
    set(figureZoomHandle,'ActionPostCallback',@minisZoom);
    
    % Create the pan object for the figure f8:
    figurePanHandle = pan(f8);
    set(figurePanHandle,'ActionPostCallback',@minisPan);
    
    
    
    % Plot neuron's input resistance:
    f9 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    axes('Position', [.04 .07 .92 .82]);
    l1 = line(tHalfF,Rseries, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    ax1 = gca;
    xlabel('Recording time (s)');
    ylabel('Resistance (Megaohms)');
    
    ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
    xlabel('Files');
    axisvals1 = get(ax1,'ylim');
    axisvals2 = get(ax2,'ylim');
    yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
    set(ax1,'xlim',[0 tEnd]);
    set(ax2,'xlim',[initFile initFile+length(files)]);
    set(ax1,'ylim',yLimit);
    set(ax2,'ylim',yLimit);
    
    xTimeMarks = [APinfuse APblock minisInfuse];
    yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
    nTimeMarks = length(xTimeMarks);
    yTimeMarks1(nTimeMarks+1 : end) = [];
    yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
    yTimeMarks2(nTimeMarks+1 : end) = [];
    line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
    l2 = zeros(1,nTimeMarks);
    for iMark = 1:nTimeMarks
        l2(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
    end
    
    set(f9, 'NumberTitle', 'off');
    set(f9, 'Name', 'Examine neuron''s input resistance');
    title('Neuron''s input resistance (assuming electrode''s resistance is balanced)');
    if nTimeMarks
        legend([l1 l2(1)],'R_N','Drug timing', 'Location','NorthEast');
    else
        legend(l1,'R_N', 'Location','NorthEast');
    end
    
    % Create the zoom object for the figure f9:
    figureZoomHandle = zoom;
    set(figureZoomHandle,'ActionPostCallback',@minisZoom);
    
    % Create the pan object for the figure f9:
    figurePanHandle = pan(f9);
    set(figurePanHandle,'ActionPostCallback',@minisPan);
    
    
    
    fclose all;
    disp('Task completed');
end