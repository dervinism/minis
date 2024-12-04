function varargout = preprocessMinis(varargin)
% PREPROCESSMINIS MATLAB code for preprocessMinis.fig
%      PREPROCESSMINIS, by itself, creates a new PREPROCESSMINIS or raises the existing
%      singleton*.
%
%      H = PREPROCESSMINIS returns the handle to a new PREPROCESSMINIS or the handle to
%      the existing singleton*.
%
%      PREPROCESSMINIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROCESSMINIS.M with the given input arguments.
%
%      PREPROCESSMINIS('Property','Value',...) creates a new PREPROCESSMINIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before preprocessMinis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to preprocessMinis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help preprocessMinis

% Last Modified by GUIDE v2.5 04-Apr-2022 11:37:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preprocessMinis_OpeningFcn, ...
                   'gui_OutputFcn',  @preprocessMinis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before preprocessMinis is made visible.
function preprocessMinis_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to preprocessMinis (see VARARGIN)

startPulse = varargin{1};
endPulse = varargin{2};
if ~isempty(startPulse) && ~isempty(endPulse)
    startPulse = regexprep(num2str(startPulse), '\s(\s*)\s', ',');
    set(handles.startPulseEdit,'String',startPulse);
    endPulse = regexprep(num2str(endPulse), '\s(\s*)\s', ',');
    set(handles.endPulseEdit,'String',endPulse);
end
handles.detectionParameters = varargin{3};
handles.classificationParameters = varargin{4};
handles.ld = varargin{5};
handles.graphicsFormats = varargin{6};
handles.parallelCores = varargin{7};
handles.minisHandles = varargin{8};

% Choose default command line output for preprocessMinis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = preprocessMinis_OutputFcn(hObject, ~, handles)  %#ok<STOUT>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Update handles structure
guidata(hObject, handles);



function loadFileInput_Callback(hObject, ~, handles) %#ok<*DEFNU>
% hObject    handle to loadFileInput (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

filename = get(hObject,'String');
set(handles.loadFileInput,'String',filename);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function loadFileInput_CreateFcn(hObject, ~, ~)
% hObject    handle to loadFileInput (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadFileButton.
function loadFileButton_Callback(hObject, ~, handles)
% hObject    handle to loadFileButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% Load an abf recording file:
dataDirectory = uigetdir(handles.ld, 'Choose Data Directory');
if dataDirectory ~= 0
    set(handles.loadFileInput,'String',dataDirectory);
    handles.ld = dataDirectory;
else
    return
end

% Update handles structure
guidata(hObject, handles);


function startPulseEdit_Callback(hObject, ~, handles)
% hObject    handle to startPulseEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

startPulseStr = get(hObject,'String');
if isempty(startPulseStr)
    startPulseStr = '...';
    set(handles.startPulseEdit,'String',startPulseStr);
    guidata(hObject, handles);
    return
end
startPulseStr = strrep(startPulseStr, ' ', '');
startPulseStrSplit = regexp(startPulseStr, ',*', 'split');
startPulse = strarray2numarray(startPulseStrSplit');
if sum(isnan(startPulse))
    msgbox('Error: At least one of the pulse vectors is NaN', 'Error', 'Error');
else
    set(handles.startPulseEdit,'String',startPulseStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function startPulseEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to startPulseEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endPulseEdit_Callback(hObject, ~, handles)
% hObject    handle to endPulseEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

endPulseStr = get(hObject,'String');
if isempty(endPulseStr)
    endPulseStr = '...';
    set(handles.endPulseEdit,'String',endPulseStr);
    guidata(hObject, handles);
    return
end
endPulseStr = strrep(endPulseStr, ' ', '');
endPulseStrSplit = regexp(endPulseStr, ',*', 'split');
endPulse = strarray2numarray(endPulseStrSplit');
if sum(isnan(endPulse))
    msgbox('Error: At least one of the pulse vectors is NaN', 'Error', 'Error');
else
    set(handles.endPulseEdit,'String',endPulseStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function endPulseEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to endPulseEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pulseAmplitudeEdit_Callback(hObject, ~, handles)
% hObject    handle to pulseAmplitudeEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

pulseAmplitudeStr = get(hObject,'String');
if isempty(pulseAmplitudeStr)
    pulseAmplitudeStr = '...';
    set(handles.pulseAmplitudeEdit,'String',pulseAmplitudeStr);
    guidata(hObject, handles);
    return
end

pulseAmplitude = str2double(pulseAmplitudeStr);
if isnan(pulseAmplitude)
    msgbox('Error: Brief (second) pulse amplitude is NaN', 'Error', 'Error');
else
    set(handles.pulseAmplitudeEdit,'String',pulseAmplitudeStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function pulseAmplitudeEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to pulseAmplitudeEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function APinfuseEdit_Callback(hObject, ~, handles)
% hObject    handle to APinfuseEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

APinfuseStr = get(hObject,'String');
if isempty(APinfuseStr)
    APinfuseStr = '...';
    set(handles.APinfuseEdit,'String',APinfuseStr);
    guidata(hObject, handles);
    return
end

APinfuse = str2double(APinfuseStr);
if isnan(APinfuse)
    msgbox('Error: AP blocker infusion time is NaN', 'Error', 'Error');
else
    set(handles.APinfuseEdit,'String',APinfuseStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function APinfuseEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to APinfuseEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function APblockEdit_Callback(hObject, ~, handles)
% hObject    handle to APblockEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

APblockStr = get(hObject,'String');
if isempty(APblockStr)
    APblockStr = '...';
    set(handles.APblockEdit,'String',APblockStr);
    guidata(hObject, handles);
    return
end

APblock = str2double(APblockStr);
if isnan(APblock)
    msgbox('Error: AP blocking time is NaN', 'Error', 'Error');
else
    set(handles.APblockEdit,'String',APblockStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function APblockEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to APblockEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minisInfuseEdit_Callback(hObject, ~, handles)
% hObject    handle to minisInfuseEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

minisInfuseStr = get(hObject,'String');
if isempty(minisInfuseStr)
    minisInfuseStr = '...';
    set(handles.minisInfuseEdit,'String',minisInfuseStr);
    guidata(hObject, handles);
    return
end

minisInfuse = str2double(minisInfuseStr);
if isnan(minisInfuse)
    msgbox('Error: Minis blocker infusion time is NaN', 'Error', 'Error');
else
    set(handles.minisInfuseEdit,'String',minisInfuseStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function minisInfuseEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to minisInfuseEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function initFileEdit_Callback(hObject, ~, handles)
% hObject    handle to initFileEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

initFileStr = get(hObject,'String');
if isempty(initFileStr)
    initFileStr = '...';
    set(handles.initFileEdit,'String',initFileStr);
    guidata(hObject, handles);
    return
end

initFile = str2double(initFileStr);
if isnan(initFile)
    msgbox('Error: Starting file number is NaN', 'Error', 'Error');
else
    set(handles.initFileEdit,'String',initFileStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function initFileEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to initFileEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in acceptButton.
function acceptButton_Callback(hObject, ~, handles)
% hObject    handle to acceptButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

try
    data = guidata(handles.preprocessMinis);
    errmsgNaN = 'Error: At least one of the pulse vectors is NaN';
    errmsgVectorLengths = 'Error: Pulse vectors are of unequal lengths';
    errmsgStartEnd = 'Error: At least one of the pulse end vector elements is smaller than its corresponding start pulse element';

    dataDirectory = get(data.loadFileInput,'String');
    if strcmpi(dataDirectory, '>>> <<<') || isempty(dataDirectory)
        errmsgNoPath = 'Error: no data path provided';
        msgbox(errmsgNoPath,'Error','Error');
        return
    end
    handles.ld = dataDirectory;

    % Pulse vector inputs:
    startPulseStr = get(handles.startPulseEdit,'String');
    endPulseStr = get(handles.endPulseEdit,'String');
    if ~strcmpi(startPulseStr, '...') && ~strcmpi(endPulseStr, '...')
        try
            [startPulse, endPulse] = testExcludedTimes(startPulseStr, endPulseStr, errmsgNaN, errmsgVectorLengths, errmsgStartEnd);
        catch %#ok<CTCH>
            return
        end
    else
        startPulse = [];
        endPulse = [];
    end
    handles.excludedTimes = struct('startPulse', startPulse, 'endPulse', endPulse);

    pulseAmplitudeStr = get(handles.pulseAmplitudeEdit,'String');
    if strcmpi(pulseAmplitudeStr, '...')
        pulseAmplitude = [];
    else
        pulseAmplitude = str2double(pulseAmplitudeStr);
        if isnan(pulseAmplitude)
            msgbox('Error: Brief (second) pulse amplitude is NaN','Error','Error');
            return
        end
    end

    % Drug infusion times:
    APinfuseStr = get(handles.APinfuseEdit,'String');
    if strcmpi(APinfuseStr, '...')
        APinfuse = [];
    else
        APinfuse = str2double(APinfuseStr);
        if isnan(APinfuse)
            msgbox('Error: AP blocker infusion time is NaN','Error','Error');
            return
        end
    end

    APblockStr = get(handles.APblockEdit,'String');
    if strcmpi(APblockStr, '...')
        APblock = [];
    else
        APblock = str2double(APblockStr);
        if isnan(APblock)
            msgbox('Error: AP blocking time is NaN','Error','Error');
            return
        end
    end

    minisInfuseStr = get(handles.minisInfuseEdit,'String');
    if strcmpi(minisInfuseStr, '...')
        minisInfuse = [];
    else
        minisInfuse = str2double(minisInfuseStr);
        if isnan(minisInfuse)
            msgbox('Error: Minis blocker infusion time is NaN','Error','Error');
            return
        end
    end

    global initFile
    initFileStr = get(handles.initFileEdit,'String');
    if strcmpi(initFileStr, '...')
        initFile = 0;
    else
        initFile = str2double(initFileStr);
        if isnan(initFile)
            msgbox('Error: Starting file number is NaN','Error','Error');
            return
        end
    end

    if strcmpi(handles.detectionParameters.pulseDuration, '...')
        options.Resize = 'on';
        options.WindowStyle = 'normal';
        pulseDuration = inputdlg('Enter exact pulse duration (ms):','Current Pulse',1,{'0.5'},options);
        if ~isempty(pulseDuration)
            pulseDuration = str2double(cell2mat(pulseDuration));
            if isnan(pulseDuration)
                msgbox('Error: Zero pulse duration', 'Error', 'Error');
                return
            end
        else
            msgbox('Error: Zero pulse duration', 'Error', 'Error');
            return
        end
        handles.detectionParameters.pulseDuration = pulseDuration;
    end

    % Preprocess the data:
    [data, data2, figures] = dataMinis(dataDirectory, handles.excludedTimes, handles.detectionParameters, handles.classificationParameters,...
        pulseAmplitude, APinfuse, APblock, minisInfuse, initFile, handles.parallelCores);

    button = questdlg('Save the pre-processing data text file?','Save File','Yes','No','Yes');
    if strcmpi(button, 'Yes')
        [eventFilename, eventPathname, filterIndex] = uiputfile({'*.txt','Text files (*.txt)'},'Save Pre-processing Log as', dataDirectory);
        if filterIndex
            handles.ld = eventPathname;
            fid = fopen(fullfile(eventPathname,eventFilename),'wt+');
            fprintf(fid, '%16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\n',...
                'File', 'Top 20% mean', 'Top 20% median', 'Top 10% mean', 'Top 10% median', 'Top 5% mean', 'Top 5% median', 'Top 2% mean', 'Top 2% median',...
                'Top 1% mean', 'Top 1% median', '100ms SD', '100ms smoothed SD', '15ms-mean SD', '15ms-mean smoothed SD', '15ms-med SD',...
                '15ms-med smoothed SD', '100ms BL', '100ms smoothed BL', '15ms BL', '15ms smoothed BL');
            fprintf(fid,'%16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\n', data');
            fprintf(fid,'\r\n');

            fprintf(fid, '%16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\n',...
                'File', 'Top 10% tau_m', 'Impulse tau_m', 'Effective tau_m', 'Capacitance', 'Effective capacitance', 'Pseudo series R');
            fprintf(fid,'%16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\n', data2');
            fprintf(fid,'\r\n');
            fclose(fid);
        end
    end

    % Save the figures:
    button = questdlg('Save the summary figures?','Save Figures','Yes','No','Yes');
    if strcmpi(button, 'Yes')
        [figureFilename, figurePathname, filterIndex] = uiputfile(handles.graphicsFormats, 'Save Top 10% Amplitude Graph as', dataDirectory);
        if filterIndex
            dataDirectory = figurePathname;
            handles.ld = dataDirectory;
            figFullName = fullfile(dataDirectory, figureFilename);
            saveas(figures(1), figFullName);
        end
        [figureFilename, figurePathname, filterIndex] = uiputfile(handles.graphicsFormats, 'Save PSP Membrane Time Constant Graph as', dataDirectory);
        if filterIndex
            dataDirectory = figurePathname;
            handles.ld = dataDirectory;
            figFullName = fullfile(dataDirectory, figureFilename);
            saveas(figures(2), figFullName);
        end
        [figureFilename, figurePathname, filterIndex] = uiputfile(handles.graphicsFormats, 'Save Standard Deviation (100ms Window) Graph as', dataDirectory);
        if filterIndex
            dataDirectory = figurePathname;
            handles.ld = dataDirectory;
            figFullName = fullfile(dataDirectory, figureFilename);
            saveas(figures(3), figFullName);
        end
        [figureFilename, figurePathname, filterIndex] = uiputfile(handles.graphicsFormats, 'Save Standard Deviation (15ms Window) Graph as', dataDirectory);
        if filterIndex
            dataDirectory = figurePathname;
            handles.ld = dataDirectory;
            figFullName = fullfile(dataDirectory, figureFilename);
            saveas(figures(4), figFullName);
        end
        [figureFilename, figurePathname, filterIndex] = uiputfile(handles.graphicsFormats, 'Save Baseline Graph as', dataDirectory);
        if filterIndex
            dataDirectory = figurePathname;
            handles.ld = dataDirectory;
            figFullName = fullfile(dataDirectory, figureFilename);
            saveas(figures(5), figFullName);
        end
        [figureFilename, figurePathname, filterIndex] = uiputfile(handles.graphicsFormats, 'Save Impulse Membrane Time Constant Graph as', dataDirectory);
        if filterIndex
            dataDirectory = figurePathname;
            handles.ld = dataDirectory;
            figFullName = fullfile(dataDirectory, figureFilename);
            saveas(figures(6), figFullName);
        end
        [figureFilename, figurePathname, filterIndex] = uiputfile(handles.graphicsFormats, 'Save Capacitance Graph as', dataDirectory);
        if filterIndex
            dataDirectory = figurePathname;
            handles.ld = dataDirectory;
            figFullName = fullfile(dataDirectory, figureFilename);
            saveas(figures(7), figFullName);
        end
        [figureFilename, figurePathname, filterIndex] = uiputfile(handles.graphicsFormats, 'Save Series Resistance Graph as', dataDirectory);
        if filterIndex
            dataDirectory = figurePathname;
            handles.ld = dataDirectory;
            figFullName = fullfile(dataDirectory, figureFilename);
            saveas(figures(8), figFullName);
        end
    end
catch me
    if strcmp(me.message, 'Interrupted by user')
        handles.minisHandles.state.String = '  State: Interrupted by user';
    else
        handles.minisHandles.state.String = '  State: Execution failed. See minisLogFile inside the working directory.';
        guidata(hObject, handles);
        disp('Program execution failed');
        rethrow(me);
    end
end

disp('Task completed');
handles.minisHandles.state.String = '  State: Task completed';

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in retrieveButton.
function retrieveButton_Callback(hObject, ~, handles)
% hObject    handle to retrieveButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

handles.ld = figureMinis(handles.ld);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function preprocessMinis_DeleteFcn(hObject, ~, handles)
% hObject    handle to detectionPanel (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

disp('minis is ready');
handles.minisHandles.state.String = '  State: Ready';
guidata(hObject, handles);
