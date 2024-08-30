function varargout = errorMinis(varargin)
% ERRORMINIS MATLAB code for errorMinis.fig
%      ERRORMINIS, by itself, creates a new ERRORMINIS or raises the existing
%      singleton*.
%
%      H = ERRORMINIS returns the handle to a new ERRORMINIS or the handle to
%      the existing singleton*.
%
%      ERRORMINIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ERRORMINIS.M with the given input arguments.
%
%      ERRORMINIS('Property','Value',...) creates a new ERRORMINIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before errorMinis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to errorMinis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help errorMinis

% Last Modified by GUIDE v2.5 04-Apr-2022 11:57:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @errorMinis_OpeningFcn, ...
                   'gui_OutputFcn',  @errorMinis_OutputFcn, ...
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


% --- Executes just before errorMinis is made visible.
function errorMinis_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to errorMinis (see VARARGIN)

excludedTimes = varargin{1};
if ~isempty(excludedTimes)
    startPulse = regexprep(num2str(excludedTimes.startPulse), '\s(\s*)\s', ',');
    set(handles.startPulseEdit,'String',startPulse);
    endPulse = regexprep(num2str(excludedTimes.endPulse), '\s(\s*)\s', ',');
    set(handles.endPulseEdit,'String',endPulse);
end
handles.detectionParameters = varargin{2};
handles.classificationParameters = varargin{3};
handles.ld = varargin{4};
handles.wd = varargin{5};
handles.parallelCores = varargin{6};
handles.minisHandles = varargin{7};

% Choose default command line output for errorMinis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = errorMinis_OutputFcn(hObject, ~, handles)  %#ok<STOUT>
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



% --- Executes on button press in acceptButton.
function acceptButton_Callback(hObject, ~, handles)
% hObject    handle to acceptButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

try
    data = guidata(handles.errorMinis);
    errmsgNaN = 'Error: At least one of the pulse vectors is NaN';
    errmsgVectorLengths = 'Error: Pulse vectors are of unequal lengths';
    errmsgStartEnd = 'Error: At least one of the pulse end vector elements is smaller than its corresponding start pulse element';

    dataDirectory = get(data.loadFileInput,'String');
    if strcmpi(dataDirectory, '>>> <<<') || isempty(dataDirectory)
        errmsgNoPath = 'Error: no data path provided';
        msgbox(errmsgNoPath,'Error','Error');
        return
    end
    if ~isdir(dataDirectory) %#ok<*ISDIR>
        errmsgNoPath = 'Error: the supplied path is not a directory';
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
    excludedTimes = struct('startPulse', startPulse, 'endPulse', endPulse);

    % Estimate the error margin for optimisation (initial stage):
    detectionParameters = handles.detectionParameters;
    classificationParameters = handles.classificationParameters;
    [Amps, AmpsNeg, RTs, RTsNeg, twoDs, twoDsNeg, fileNames, fileSweeps, SD, fileAmps, fileRTs, fileAmpsNeg, fileRTsNeg, sLength, sweepcount, handles.ld] = estErrBounds(...
        excludedTimes, detectionParameters, classificationParameters, handles.ld, handles.parallelCores);
    if isempty(Amps)
        return
    end

    % Save the initial data for later retrieval:
    button = questdlg('Save the intermediate error bound data for later retrieval?','Save File','Yes','No','Yes');
    if strcmpi(button, 'Yes')
        [eventFilename, eventPathname, filterIndex] = uiputfile({'*.mat', 'MAT files (*.mat)'},'Save data as', handles.ld);
        if filterIndex
            handles.ld = eventPathname;
            eventFilename = fullfile(eventPathname, eventFilename);
            save(eventFilename, 'Amps', 'AmpsNeg', 'RTs', 'RTsNeg', 'twoDs', 'twoDsNeg', 'fileNames', 'fileSweeps', 'SD', 'fileAmps', 'fileRTs',...
                'fileAmpsNeg', 'fileRTsNeg', 'sLength', 'detectionParameters', 'classificationParameters', 'sweepcount');
        end
    end

    % Resume error margin estimation:
    expand = '15-score'; % 3-score, 6-score, 9-score, or 15-score 50th centile
    [AmpsMean, RTsMean, TwoDsMean, AmpsMedian, RTsMedian, TwoDsMedian, AmpsMin, RTsMin, TwoDsMin, AmpsMax, RTsMax, TwoDsMax, AmpsPrct, RTsPrct, TwoDsPrct,...
        boundAmps, boundRTs, boundTwoDs, boundAmpsLinear, boundRTsLinear, boundTwoDsLinear, AmpsMeanBottom, AmpsMedianBottom, AmpsMinBottom, AmpsMaxBottom,...
        AmpsPrctBottom, boundAmpsBottom, boundAmpsLinearBottom, AmpsMeanMid, AmpsMedianMid, AmpsMinMid, AmpsMaxMid, AmpsPrctMid, boundAmpsMid, boundAmpsLinearMid,...
        AmpsMeanTop, AmpsMedianTop, AmpsMinTop, AmpsMaxTop, AmpsPrctTop, boundAmpsTop, boundAmpsLinearTop, sLength, slFactor, fileNames, fileSweeps, flFactor, pow, powWgh,...
        estCount, UFest, AmpsLSE, RTsLSE, TwoDsLSE, AmpsLinearLSE, RTsLinearLSE, TwoDsLinearLSE, AmpsLSEBottom, AmpsLinearLSEBottom, AmpsLSEMid, AmpsLinearLSEMid,...
        AmpsLSETop, AmpsLinearLSETop, SD, optimData, dataF, cntF, F, G, H, A, B, C, FBottom, ABottom, FMid, AMid, FTop, ATop] = estErrBounds2(Amps, AmpsNeg, RTs, RTsNeg,...
        twoDs, twoDsNeg, fileNames, fileSweeps, SD, fileAmps, fileRTs, fileAmpsNeg, fileRTsNeg, sLength, detectionParameters, classificationParameters,...
        sweepcount, expand);

    handles.ld = saveErrorBounds(AmpsMean, RTsMean, TwoDsMean, AmpsMedian, RTsMedian, TwoDsMedian, AmpsMin, RTsMin, TwoDsMin, AmpsMax, RTsMax, TwoDsMax,...
        AmpsPrct, RTsPrct, TwoDsPrct, boundAmps, boundRTs, boundTwoDs, boundAmpsLinear, boundRTsLinear, boundTwoDsLinear, AmpsMeanBottom, AmpsMedianBottom,...
        AmpsMinBottom, AmpsMaxBottom, AmpsPrctBottom, boundAmpsBottom, boundAmpsLinearBottom, AmpsMeanMid, AmpsMedianMid, AmpsMinMid, AmpsMaxMid, AmpsPrctMid,...
        boundAmpsMid, boundAmpsLinearMid, AmpsMeanTop, AmpsMedianTop, AmpsMinTop, AmpsMaxTop, AmpsPrctTop, boundAmpsTop, boundAmpsLinearTop, sLength,...
        slFactor, fileNames, fileSweeps, flFactor, pow, powWgh, estCount, UFest, AmpsLSE, RTsLSE, TwoDsLSE, AmpsLinearLSE, RTsLinearLSE, TwoDsLinearLSE,...
        AmpsLSEBottom, AmpsLinearLSEBottom, AmpsLSEMid, AmpsLinearLSEMid, AmpsLSETop, AmpsLinearLSETop, SD, optimData, dataF, cntF, F, G, H, A, B, C, FBottom,...
        ABottom, FMid, AMid, FTop, ATop, detectionParameters, classificationParameters, handles.ld, handles.wd, expand);
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

fclose all;
disp('Task completed');
handles.minisHandles.state.String = '  State: Task completed';

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in stopButton.
function stopButton_Callback(hObject, ~, handles)
% hObject    handle to stopButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

try
    % Retrieve data for estimating the error margin for optimisation:
    [dataFilename, dataPathname, filterIndex] = uigetfile({'*.mat','MAT files (*.mat)'}, 'Provide data file for resuming error bound estimation', handles.ld);
    if filterIndex
        dataFilename = fullfile(dataPathname, dataFilename);
        handles.ld = dataPathname;
        set(handles.loadFileInput,'String',dataFilename);
        load(dataFilename); %#ok<*LOAD>
        if ~exist('Amps', 'var') || ~exist('AmpsNeg', 'var') || ~exist('RTs', 'var') || ~exist('RTsNeg', 'var') || ~exist('twoDs', 'var')...
                || ~exist('twoDsNeg', 'var') || ~exist('fileNames', 'var') || ~exist('fileSweeps', 'var') || ~exist('SD', 'var') || ~exist('fileAmps', 'var')...
                || ~exist('fileRTs', 'var') || ~exist('fileAmpsNeg', 'var') || ~exist('fileRTsNeg', 'var') || ~exist('sLength', 'var')...
                || ~exist('detectionParameters', 'var') || ~exist('classificationParameters', 'var') || ~exist('sweepcount', 'var')
            msgbox('The supplied data file is corrupt','Error','Error');
            return
        end

        % Resume error margin estimation:
        expand = '15-score'; % 3-score, 6-score, or 9-score 50th centile
        [AmpsMean, RTsMean, TwoDsMean, AmpsMedian, RTsMedian, TwoDsMedian, AmpsMin, RTsMin, TwoDsMin, AmpsMax, RTsMax, TwoDsMax, AmpsPrct, RTsPrct, TwoDsPrct,...
            boundAmps, boundRTs, boundTwoDs, boundAmpsLinear, boundRTsLinear, boundTwoDsLinear, AmpsMeanBottom, AmpsMedianBottom, AmpsMinBottom, AmpsMaxBottom,...
            AmpsPrctBottom, boundAmpsBottom, boundAmpsLinearBottom, AmpsMeanMid, AmpsMedianMid, AmpsMinMid, AmpsMaxMid, AmpsPrctMid, boundAmpsMid, boundAmpsLinearMid,...
            AmpsMeanTop, AmpsMedianTop, AmpsMinTop, AmpsMaxTop, AmpsPrctTop, boundAmpsTop, boundAmpsLinearTop, sLength, slFactor, fileNames, fileSweeps, flFactor, pow, powWgh,...
            estCount, UFest, AmpsLSE, RTsLSE, TwoDsLSE, AmpsLinearLSE, RTsLinearLSE, TwoDsLinearLSE, AmpsLSEBottom, AmpsLinearLSEBottom, AmpsLSEMid, AmpsLinearLSEMid,...
            AmpsLSETop, AmpsLinearLSETop, SD, optimData, dataF, cntF, F, G, H, A, B, C, FBottom, ABottom, FMid, AMid, FTop, ATop] = estErrBounds2(Amps, AmpsNeg, RTs, RTsNeg,...
            twoDs, twoDsNeg, fileNames, fileSweeps, SD, fileAmps, fileRTs, fileAmpsNeg, fileRTsNeg, sLength, detectionParameters, classificationParameters,...
            sweepcount, expand); %#ok<NODEF>

        handles.ld = saveErrorBounds(AmpsMean, RTsMean, TwoDsMean, AmpsMedian, RTsMedian, TwoDsMedian, AmpsMin, RTsMin, TwoDsMin, AmpsMax, RTsMax, TwoDsMax,...
            AmpsPrct, RTsPrct, TwoDsPrct, boundAmps, boundRTs, boundTwoDs, boundAmpsLinear, boundRTsLinear, boundTwoDsLinear, AmpsMeanBottom, AmpsMedianBottom,...
            AmpsMinBottom, AmpsMaxBottom, AmpsPrctBottom, boundAmpsBottom, boundAmpsLinearBottom, AmpsMeanMid, AmpsMedianMid, AmpsMinMid, AmpsMaxMid, AmpsPrctMid,...
            boundAmpsMid, boundAmpsLinearMid, AmpsMeanTop, AmpsMedianTop, AmpsMinTop, AmpsMaxTop, AmpsPrctTop, boundAmpsTop, boundAmpsLinearTop, sLength,...
            slFactor, fileNames, fileSweeps, flFactor, pow, powWgh, estCount, UFest, AmpsLSE, RTsLSE, TwoDsLSE, AmpsLinearLSE, RTsLinearLSE, TwoDsLinearLSE,...
            AmpsLSEBottom, AmpsLinearLSEBottom, AmpsLSEMid, AmpsLinearLSEMid, AmpsLSETop, AmpsLinearLSETop, SD, optimData, dataF, cntF, F, G, H, A, B, C, FBottom,...
            ABottom, FMid, AMid, FTop, ATop, detectionParameters, classificationParameters, handles.ld, handles.wd, expand);
    else
        return
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

fclose all;
disp('Task completed');
handles.minisHandles.state.String = '  State: Task completed';

% Update handles structure
guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function errorMinis_DeleteFcn(hObject, ~, handles)
% hObject    handle to errorMinis (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

disp('minis is ready');
handles.minisHandles.state.String = '  State: Ready';
guidata(hObject, handles);