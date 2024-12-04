function manualFitMinisGUI(varargin)
% MANUALFITMINISGUI MATLAB code for manualFitMinisGUI.fig
%      MANUALFITMINISGUI, by itself, creates a new MANUALFITMINISGUI or raises the existing
%      singleton*.
%
%      H = MANUALFITMINISGUI returns the handle to a new MANUALFITMINIS or the handle to
%      the existing singleton*.
%
%      MANUALFITMINISGUI('CALLBACK',hObject,~,handles,...) calls the local
%      function named CALLBACK in MANUALFITMINISGUI.M with the given input arguments.
%
%      MANUALFITMINISGUI('Property','Value',...) creates a new MANUALFITMINISGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before manualFitMinisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to manualFitMinisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help manualFitMinisGUI

% Last Modified by GUIDE v2.5 27-Dec-2011 19:35:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @manualFitMinisGUI_OpeningFcn, ...
    'gui_OutputFcn',  @manualFitMinisGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:}); %#ok<NASGU>
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before manualFitMinisGUI is made visible.
function manualFitMinisGUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manualFitMinisGUI (see VARARGIN)

% Initialise variables:
handles.wd = varargin{1};
handles.V = varargin{2};
handles.minisSimulated = varargin{3};
handles.target2D = varargin{4};
handles.simulated2D = varargin{5};
handles.difference2D = handles.target2D - handles.simulated2D;
handles.difference1D = sum(handles.difference2D,1);
handles.lowestCost1D = sum(abs(handles.difference1D));
handles.lowestCost2D = sum(sum(abs(handles.difference2D)));
handles.minis2D = varargin{6};
handles.minis1D = sum(handles.minis2D,1);
handles.lengthRatio = varargin{7};
handles.noiseProperties = varargin{8};
handles.detectionParameters = varargin{9};
handles.simulationParameters = varargin{10};
classificationParameters = varargin{11};
handles.amplitudeArray = classificationParameters.amplitudeArray;
handles.amplitudeArrayExt = classificationParameters.amplitudeArrayExt;
handles.ampLim = (handles.amplitudeArray(2) - handles.amplitudeArray(1))/2;
handles.riseTimeArray = classificationParameters.riseTimeArray;
handles.riseTimeArrayExt = classificationParameters.riseTimeArrayExt;
handles.graphicsFormats = varargin{12};
handles.parallelCores = varargin{13};
handles.F = varargin{14};
handles.options = varargin{15};
handles.shapes = varargin{16};

% Initialise the fitting window:
axes(handles.axes1); %#ok<*MAXES>
plot(gca, handles.amplitudeArray, sum(handles.target2D,1), 'b.-');
xlim(gca, [handles.amplitudeArray(1)-handles.ampLim handles.amplitudeArray(end)+handles.ampLim]);
hold on
plot(gca, handles.amplitudeArray, sum(handles.simulated2D,1),'r.-');
plot(gca, handles.amplitudeArray, sum(handles.target2D-handles.simulated2D,1),'g.-');
plot(gca, handles.amplitudeArray, sum(handles.minis2D,1),'k.-');
title('Amplitude distribution fit');
xlabel('Amplitude(mV)');
ylabel('Number of Events');
legend('Target','Simulated','Target-simulated','Location','NorthEast');

% Update handles structure
guidata(hObject, handles);

iterate(hObject, handles);

% UIWAIT makes manualFitMinisGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function manualFitMinisGUI_OutputFcn(~, ~, ~)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FileMenu_Callback(~, ~, ~) %#ok<*DEFNU>
% hObject    handle to FileMenu (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function PrintMenuItem_Callback(~, ~, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(~, ~, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg('Close the manual distribution fitting window?','Close Window','Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, ~, handles)
% hObject    handle to cancelButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

handleIndices = find(handles.bins(:,2));
existingHandles = ishandle(handles.bins(handleIndices,3));
iexistingHandles = handleIndices(existingHandles);
delete([handles.bins(iexistingHandles,3)' handles.bins(iexistingHandles,4)']);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in clearButton.
function clearButton_Callback(hObject, ~, handles)
% hObject    handle to clearButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

if handles.bins(handles.i,2)
    delete([handles.bins(handles.i,3) handles.bins(handles.i,4)]);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on button press in enterButton.
function enterButton_Callback(hObject, ~, handles)
% hObject    handle to enterButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

set(handles.enterButton, 'UserData', 1);

% Adjust the minis distribution:
adjustmentIndices = find(handles.bins(:,2));
handles.minis1D(adjustmentIndices) = round(handles.lengthRatio*handles.bins(adjustmentIndices,2));
guidata(hObject, handles);

% Simulate minis:
[~, handles.V, ~, ~, handles.shapes] = simulateMinis('Zero', handles.simulationParameters, 1, 'Arbitrary', handles.minis1D, handles.noiseProperties,...
    'basal', handles.parallelCores, handles.amplitudeArray);
handles.minis2D = hist2d(handles.shapes(:,2), handles.shapes(:,3), handles.amplitudeArrayExt, handles.riseTimeArrayExt);
handles.minis2D = handles.minis2D(1:end-1,1:end-1);
handles.minis1D = sum(handles.minis2D,1);

% Detect simulated events:
waveform.estimate = 2;
waveform.tau_m = simulationParameters.tau_m;
if ishandle(handles.F(1)) && handles.F(1)
    close(handles.F(1));
end
filt.state = 'off';
[handles.minisSimulated, ~, ~, simWave, handles.F(1)] = detectMinis(handles.V, handles.noiseProperties.excludedTimes, handles.detectionParameters, filt,...
    waveform, handles.parallelCores, handles.options);
tau = ones(size(handles.minisSimulated,1),1)*simWave.parameters.tau;
handles.minisSimulated = [handles.minisSimulated tau];
handles.simulated2D = handles.lengthRatio*hist2d(handles.minisSimulated(:,4), handles.minisSimulated(:,12), handles.amplitudeArrayExt, handles.riseTimeArrayExt);
handles.simulated2D = handles.simulated2D(1:end-1,1:end-1);
handles.difference2D = handles.target2D - handles.simulated2D;
handles.difference1D = sum(handles.difference2D,1);

% Evaluate the fitness of the new distribution and save best examples:
handles.cost1D = sum(abs(handles.difference1D));
handles.cost2D = sum(sum(abs(handles.difference2D)));
if handles.cost1D > handles.lowestCost1D
    handles.lowestCost1D = handles.cost1D;
    handles.lowestCost2D = handles.cost2D;
    filename = num2str(handles.cost1D);
    filename(filename == '.') = '_';
    dataFilename = strcat(filename,'.mat');
    difference1D = handles.difference1D;
    difference2D = handles.difference2D;
    cost1D = handles.cost1D; 
    cost2D = handles.cost2D; 
    amplitudeArray = handles.amplitudeArray; 
    riseTimeArray = handles.riseTimeArray; 
    minis1D = handles.minis1D; 
    minis2D = handles.minis2D; 
    shapes = handles.shapes; 
    simulated2D = handles.simulated2D; 
    target2D = handles.target2D; 
    if ~isdeployed
        cd(cell2mat(handles.wd));
    end
    dataDir = 'data';
    if ~exist(dataDir,'dir')
        mkdir(dataDir);
    end
    cd(dataDir);
    save(dataFilename,'difference1D','difference2D','cost1D','cost2D','amplitudeArray','riseTimeArray','minis1D','minis2D','shapes','simulated2D','target2D');
    
    filename = strcat(filename,'.abf');
    if handles.noiseProperties.nchans_to_save == 1
        writeABF(handles.V,filename,1000/handles.noiseProperties.dt,{'mV'});
    elseif handles.noiseProperties.nchans_to_save == 2
        writeABF([handles.V; handles.noiseProperties.current],filename,1000/handles.noiseProperties.dt,{'mV';'pA'});
    else
        disp('error in noiseProperties.nchans_to_save');
    end
    cd ..
end

% Update the figures:
[handles.F(2),handles.F(3),handles.F(4)] = plotMinis(handles.amplitudeArray, handles.riseTimeArray, handles.target2D, handles.simulated2D, handles.minis2D,...
    handles.shapes,handles.F(2),handles.F(3),handles.F(4));

% Update the fitting window:
axes(handles.axes1);
cla;
plot(gca, handles.amplitudeArray, sum(handles.target2D,1), 'b.-');
xlim(gca, [handles.amplitudeArray(1)-handles.ampLim handles.amplitudeArray(end)+handles.ampLim]);
hold on
plot(gca, handles.amplitudeArray, sum(handles.simulated2D,1),'r.-');
plot(gca, handles.amplitudeArray, handles.difference1D,'g.-');
plot(gca, handles.amplitudeArray, handles.minis1D,'k.-');
title('Amplitude distribution fit');
xlabel('Amplitude(mV)');
ylabel('Number of Events');
legend('Target','Simulated','Target-simulated','Location','NorthEast');

% Update handles structure
guidata(hObject, handles);

iterate(hObject, handles);


% --- Executes on button press in saveButton.
function saveButton_Callback(~, ~, handles)
% hObject    handle to saveButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

if ~isdeployed
    cd(cell2mat(handles.wd));
end
dataDir = 'data';
if ~exist(dataDir,'dir')
    mkdir(dataDir);
end
cd(dataDir);

% Save the noise + simulated minis recording trace in abf format:
[filename, ~, filterIndex] = uiputfile({'*.abf','Axon ABF files (*.abf)'},'Save Noise + Simulated Minis Data as', 'simulation recording trace.abf');
if filterIndex
    if handles.noiseProperties.nchans_to_save == 1
        writeABF(handles.V,filename,1000/handles.noiseProperties.dt,{'mV'});
    elseif handles.noiseProperties.nchans_to_save == 2
        writeABF([handles.V; handles.noiseProperties.current],filename,1000/handles.noiseProperties.dt,{'mV';'pA'});
    else
        disp('error in handles.noiseProperties.nchans_to_save');
    end
end

% Save the event log file:
[filename, ~, filterIndex] = uiputfile({'*.txt','Text files (*.txt)'},'Save Detected Event Log as', 'detected events.txt');
if filterIndex
    fid = fopen(filename,'wt+');
    fprintf(fid, '%16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\n',...
        'Peak potential', 'Peak Time', 'Peak Index', 'Amplitude', 'Baseline (BL)', 'BL start time', 'BL end time', 'BL start index', 'BL end index',...
        'Rise time (RT) length', 'RT', '10-90% RT', '10% RT time mark', '50% RT time mark', '90% RT time mark',...
        '10% RT index', '50% RT index', '90% RT index', '10% RT potential', '50% RT potential', '90% RT potential', 'tau');
    fprintf(fid,'%16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\n',handles.minisSimulated');
    fclose(fid);
end

% Save the detected event figure file:
if handles.F(1)
    [filename, pathname, filterIndex] = uiputfile(handles.graphicsFormats,'Save detected event graph as', 'detected events');
    if filterIndex
        figfullname = fullfile(pathname, filename);
        saveas(handles.F(1), figfullname);
    end
end

% Save event distribution histograms:
[filename, pathname, filterIndex] = uiputfile(handles.graphicsFormats,'Save Event Histograms as', 'event histograms');
if filterIndex
    figFullName = fullfile(pathname, filename);
    saveas(handles.F(2), figFullName);
end

[filename, pathname, filterIndex] = uiputfile(handles.graphicsFormats,'Save Event 2D Histograms as', 'event 2D histograms');
if filterIndex
    figFullName = fullfile(pathname, filename);
    saveas(handles.F(3), figFullName);
end

[filename, pathname, filterIndex] = uiputfile(handles.graphicsFormats,'Save Simulated Minis Histograms as', 'minis histograms');
if filterIndex
    figFullName = fullfile(pathname, filename);
    saveas(handles.F(4), figFullName);
end

% Save the .MAT data file:
[filename, ~, filterIndex] = uiputfile({'*.mat','Matlab MAT files (*.mat)'},'Save Workspace Variables as', 'workspace variables.mat');
if filterIndex
    difference1D = handles.difference1D; 
    difference2D = handles.difference2D; 
    cost1D = handles.lowestCost1D; 
    cost2D = handles.lowestCost2D; 
    amplitudeArray = handles.amplitudeArray; 
    riseTimeArray = handles.riseTimeArray; 
    minis1D = handles.minis1D; 
    minis2D = handles.minis2D; 
    shapes = handles.shapes; 
    simulated2D = handles.simulated2D; 
    target2D = handles.target2D; 
    save(filename,'difference1D','difference2D','cost1D','cost2D','amplitudeArray','riseTimeArray','minis1D','minis2D','shapes','simulated2D','target2D');
end

cd ..


% --- Executes on button press in parameterButton.
function parameterButton_Callback(hObject, ~, handles)
% hObject    handle to parameterButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

try
    simulationParameters = changeParametersGUI(handles.simulationParameters);
    set(handles.parameterButton, 'UserData', simulationParameters);
catch %#ok<CTCH>
end
guidata(hObject, handles);


function iterate(hObject, handles)
% ITERATE is a helper subfunction of manualFitMinisGUI. It gathers the user
% graphical input and is iterated until the enter button is pressed.
%

handles.minis1D = sum(handles.minis2D,1);
handles.bins = zeros(length(handles.amplitudeArray),4);
handles.bins(:,1) = round(handles.amplitudeArray*100)/100;
set(handles.enterButton, 'UserData', 0);                                    % enterButton controls the while loop execution
set(handles.parameterButton, 'UserData', []);
while get(handles.enterButton, 'UserData') == 0
    try
        [x,y] = ginput(1);
    catch %#ok<CTCH>
        return
    end
    if ~isempty(get(handles.parameterButton, 'UserData'))
        handles.simulationParameters = get(handles.parameterButton, 'UserData');
        guidata(hObject, handles);
    end
    if get(handles.enterButton, 'UserData') == 0
        if ~isempty(x) && ~isempty(y)
            xpos = round(100*x)/100;                                        % dummy operation for matching "handles.bins(:,1) = round(handles.amplitudeArray*100)/100;"
            ypos = round(y);
            if xpos >= handles.amplitudeArray(1) && xpos <= handles.amplitudeArray(end) && ypos >= 0
                p1 = plot(xpos, ypos, 'k.', 'markersize', 1);
                p2 = plot(xpos, ypos, 'ko');
                handles.i = find(handles.bins(:,1) == xpos);
                handles.bins(handles.i,2) = ypos;
                if handles.bins(handles.i,3) && ishandle(handles.bins(handles.i,3))
                    delete([handles.bins(handles.i,3) handles.bins(handles.i,4)]);
                    handles.bins(handles.i,3:4) = [p1 p2];
                else
                    handles.bins(handles.i,3) = p1;
                    handles.bins(handles.i,4) = p2;
                end
            end
        end
        % Delete cleared or cancelled values based on existent handles:
        handleIndices = find(handles.bins(:,2));
        existingHandles = ishandle(handles.bins(handleIndices,3));
        if sum(existingHandles) < length(handleIndices)
            nonexistentHandle = find(~existingHandles);
            iNonexistentHandle = handleIndices(nonexistentHandle); %#ok<FNDSB>
            handles.bins(iNonexistentHandle,2:4) = 0;
        end
        guidata(hObject, handles);
    end
end