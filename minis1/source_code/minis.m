function varargout = minis(varargin)
%MINIS M-file for minis.fig
%      MINIS, by itself, creates a new MINIS or raises the existing
%      singleton*.
%
%      H = MINIS returns the handle to a new MINIS or the handle to
%      the existing singleton*.
%
%      MINIS('Property','Value',...) creates a new MINIS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to minis_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MINIS('CALLBACK') and MINIS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MINIS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help minis

% Last Modified by GUIDE v2.5 03-Apr-2022 23:49:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @minis_OpeningFcn, ...
    'gui_OutputFcn',  @minis_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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


% --- Executes just before minis is made visible.
function minis_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Declare initialisation
disp('minis is initialising');

% Suppress javacomponent warning:
warning('off', 'MATLAB:ui:javacomponent:FunctionToBeRemoved');

% Simulation Panel:
backColour = javaObjectEDT(java.awt.Color(.79, .82, .87));
pposSimulationPanel = getpixelposition(handles.simulationPanel);

% loAmp text:
pposLoAmpText = [.0975*pposSimulationPanel(3) pposSimulationPanel(2)+.48*pposSimulationPanel(4) .20*pposSimulationPanel(3) .25*pposSimulationPanel(4)];
labelStr = '<html><font style="font-size:1.3em;" face = "MS Sans Serif">Amplitude lower bound (mV):</font></html>';
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
setBackground(jLabel, backColour);
[handles.pposLoAmpText, pposLoAmpTextContainer] = javacomponent(jLabel, pposLoAmpText, gcf); %#ok<*JAVCM>
set(pposLoAmpTextContainer,'units','normalized');

% L text:
pposL = [.0975*pposSimulationPanel(3) pposSimulationPanel(2)+.15*pposSimulationPanel(4) .20*pposSimulationPanel(3) .25*pposSimulationPanel(4)];
labelStr = '<html><font style="font-size:1.3em;" face = "MS Sans Serif">Initial L (electrotonic length, &#955):</font></html>';
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
setBackground(jLabel, backColour);
[handles.L, LContainer] = javacomponent(jLabel, pposL, gcf);
set(LContainer,'units','normalized');

% tau_m text:
ppostau_m = [.505*pposSimulationPanel(3) pposSimulationPanel(2)+.44*pposSimulationPanel(4) .3*pposSimulationPanel(3) .25*pposSimulationPanel(4)];
labelStr = '<html><font style="font-size:1.3em;" face = "MS Sans Serif">&#964<sub>m</sub> (passive membrane time constant, ms; lower limit):</font></html>';
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
setBackground(jLabel, backColour);
[handles.tau_m, tau_mContainer] = javacomponent(jLabel, ppostau_m, gcf);
set(tau_mContainer,'units','normalized');

% tau_PSPm text:
ppostau_PSPm = [.505*pposSimulationPanel(3) pposSimulationPanel(2)+.12*pposSimulationPanel(4) .3*pposSimulationPanel(3) .25*pposSimulationPanel(4)];
labelStr = '<html><font style="font-size:1.3em;" face = "MS Sans Serif">&#964<sub>m</sub> (passive membrane time constant, ms; upper limit):</font></html>';
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
setBackground(jLabel, backColour);
[handles.tau_PSPm, tau_PSPmContainer] = javacomponent(jLabel, ppostau_PSPm, gcf);
set(tau_PSPmContainer,'units','normalized');


% Register source directory:
if (~isdeployed)
    currentFullFilename = mfilename('fullpath');
    currentPath = fileparts(currentFullFilename);
    addpath(genpath(currentPath));
else
    if ispc
        [~, currentPath] = system('cd');
    elseif isunix
        [~, currentPath] = system('pwd');
    end
end
handles.wd = currentPath;
handles.ld = currentPath;

% Initiate logging
if (~isdeployed)
    diary([currentPath filesep 'minisLogFile']);
    diary on
    disp(['Logging is on, see ' currentPath filesep 'minisLogFile']);
end

% Initialise default population option variables:
value = get(handles.distTypeEdit,'Value');
distType = get(handles.distTypeEdit,'String');
distType = distType(value);
if strcmpi(distType,'Log Normal') || strcmpi(distType,'Bimodal Log Normal') || strcmpi(distType,'Trimodal Log Normal')
    %          mu1_1       A_1  alpha_1 scale_1     mu2_1    B_1 beta_1    rho_1    mu1_2    A_2   alpha_2 scale_2  mu2_2      B_2   beta_2  rho_2     mu1_3    A_3  alpha_3  scale_3    mu2_3      B_3   beta_3    rho_3
    loBounds = [-0.1,      .01,     .01,      0,       -1,   .01,   .01,      -1,    -0.1,   .01,      .01, -10000,    -1,     .01,     .01,    -1,     -0.1,   .01,     .01,  -10000,      -1,     .01,     .01,      -1];
    upBounds = [   1,        1,     1.5,  10000,       10,    20,   1.5,       1,       1,     1,      1.5,  10000,    10,      20,     1.5,     1,        1,     1,     1.5,   10000,      10,      20,     1.5,       1];
elseif ~strcmpi(distType,'Skew-normal') && ~strcmpi(distType,'Bimodal skew-normal') && ~strcmpi(distType,'Trimodal skew-normal')
    %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  mu1_2  sigma1_2  scale_2  mu2_2  sigma2_2  rho_2  mu1_3  sigma1_3  scale_3  mu2_3  sigma2_3  rho_3    mu1_4  sigma1_4  scale_4  mu2_4  sigma2_4    rho_4
    loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,   -0.1,      .01,  -10000,    -1,      .01,      -1];
    upBounds = [   1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,   	  1,        1,   10000,    10,       10,       1];
elseif strcmpi(distType,'Skew-normal')
    %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  skew1_1  skew2_1
    loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,     -10,     -10, -10000,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,     -10,     -10,      -10,   -10,      -10,     -10];
    upBounds = [   1,        1,   10000,    10,       10,     1,      10,      10,  10000,    10,       10,     1,     2,        1,   10000,    10,       10,     1,      10,      10,       10,    10,       10,      10];
elseif ~strcmpi(distType,'Trimodal skew-normal')
    %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  mu1_2  sigma1_2  scale_2  mu2_2  sigma2_2  rho_2  skew1_1  skew2_1  skew1_2  skew2_2
    loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,   .01,      .01,  -10000,    -1,      .01,    -1,     -10,     -10,     -10,     -10,   .01,    -1,     -10,     -10,      -10,   -10,      -10,     -10];
    upBounds = [   1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,      10,      10,      10,      10,    10,     1,      10,      10,       10,    10,       10,      10];
else
    %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  mu1_2  sigma1_2  scale_2  mu2_2  sigma2_2  rho_2  mu1_3  sigma1_3  scale_3  mu2_3  sigma2_3  rho_3  skew1_1  skew2_1  skew1_2  skew2_2  skew1_3  skew2_3
    loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,     -10,     -10,     -10,     -10,     -10,     -10];
    upBounds = [   1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,      10,      10,      10,      10,      10,      10];
end
bounds = [loBounds; upBounds];

delete(gcp('nocreate'));
%pool = parpool;
%parallelCores = pool.NumWorkers;
parallelCores = 1;

handles.options = struct('bounds', bounds, 'nGenerations', 50, 'parallelCores', parallelCores, 'fullParallel', 0, 'tauRange', 0,...
    'figureDisplay', 1, 'cluster', 0, 'clusterProfile', 'local', 'cliff', 0, 'SDlobound', 0, 'SDupbound', 0);
  
% faf
if fafTest
    % Load settings file:
    settings = load('fafSettings.mat');
    settings = settings.settings;
    
    % File panel:
    set(handles.loadTargetFileInput, 'string', settings.loadTargetFileInput);
    set(handles.loadNoiseFileInput, 'string', settings.loadNoiseFileInput);

    % Detection Parameters panel:
    set(handles.maxTimeToPeakEdit, 'string', settings.maxTimeToPeakEdit);
    set(handles.baselineDurationEdit, 'string', settings.baselineDurationEdit);
    set(handles.peakIntegrationPeriodEdit, 'string', settings.peakIntegrationPeriodEdit);
    set(handles.AmploboundEdit, 'string', settings.AmploboundEdit);
    set(handles.AmpupboundEdit, 'string', settings.AmpupboundEdit);
    set(handles.smoothWindowEdit, 'string', settings.smoothWindowEdit);
    set(handles.RTintEdit, 'value', settings.RTintEdit);
    set(handles.RTbinSizeEdit, 'value', settings.RTbinSizeEdit);
    set(handles.startPulseTargetEdit, 'string', settings.startPulseTargetEdit);
    set(handles.endPulseTargetEdit, 'string', settings.endPulseTargetEdit);
    set(handles.startGlitchTargetEdit, 'string', settings.startGlitchTargetEdit);
    set(handles.endGlitchTargetEdit, 'string', settings.endGlitchTargetEdit);
    set(handles.startPulseNoiseEdit, 'string', settings.startPulseNoiseEdit);
    set(handles.endPulseNoiseEdit,'string',settings.endPulseNoiseEdit);
    set(handles.startGlitchNoiseEdit, 'string', settings.startGlitchNoiseEdit);
    set(handles.endGlitchNoiseEdit, 'string', settings.endGlitchNoiseEdit);
    set(handles.pulseDurationEdit, 'string', settings.pulseDurationEdit);
    set(handles.downGoingCheckbox, 'value', settings.downGoingCheckbox);
    set(handles.voltageClampCheckbox, 'value', settings.voltageClampCheckbox);

    % Optimisation Parameters panel:
    set(handles.distTypeEdit, 'value', settings.distTypeEdit);
    set(handles.distBaselineEdit, 'value', settings.distBaselineEdit);
    set(handles.SDlobound, 'string', settings.SDlobound);
    set(handles.SDupbound, 'string', settings.SDupbound);
    set(handles.maxErrEdit, 'string', settings.maxErrEdit);
    set(handles.maxAmpErrEdit, 'string', settings.maxAmpErrEdit);
    set(handles.maxRTErrEdit, 'string', settings.maxRTErrEdit);
    set(handles.maxDevEdit, 'string', settings.maxDevEdit);
    set(handles.maxDevAmpEdit, 'string', settings.maxDevAmpEdit);
    set(handles.maxDevRTEdit, 'string', settings.maxDevRTEdit);
    set(handles.maxAmpBottomErrEdit, 'string', settings.maxAmpBottomErrEdit);
    set(handles.maxAmpBottomDevEdit, 'String', settings.maxAmpBottomDevEdit);
    set(handles.maxAmpMidErrEdit, 'string', settings.maxAmpMidErrEdit);
    set(handles.maxAmpMidDevEdit, 'String', settings.maxAmpMidDevEdit);
    set(handles.maxAmpTopErrEdit, 'string', settings.maxAmpTopErrEdit);
    set(handles.maxAmpTopDevEdit,'String', settings.maxAmpTopDevEdit);

    % Simulation Parameters panel:
    set(handles.loSimAmpEdit, 'string', settings.loSimAmpEdit);
    set(handles.LEdit, 'string', settings.LEdit);
    set(handles.tau_mEdit, 'string', settings.tau_mEdit);
    set(handles.tau_PSPmEdit, 'string', settings.tau_PSPmEdit);

    handles.options = settings.options;

    parallelCores = str2double(handles.options.parallelCores);
    pool = gcp('nocreate');
    if isempty(pool)
        openWorkers = 0;
    else
        openWorkers = pool.NumWorkers;
    end
    scheduler = parcluster('local');
    if parallelCores > scheduler.NumWorkers
      parallelCores = scheduler.NumWorkers;
    end
    if openWorkers ~= parallelCores
      delete(gcp('nocreate'));
      if parallelCores > 1
        parpool(parallelCores);
        handles.options.parallelCores = parallelCores;
      end
    end
end

% Choose default command line output for minis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = minis_OutputFcn(hObject, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% handles    structure with handles and user data (see GUIDATA)

% Declare the end of the initialisation process
handles.state.String = '  State: Ready';
guidata(hObject, handles);
disp('minis is ready');

% Get default command line output from handles structure
varargout{1} = handles.output;










%% Files panel:
function loadTargetFileInput_Callback(hObject, ~, handles) %#ok<*DEFNU>
% hObject    handle to loadTargetFileInput (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

targetFilename = get(hObject,'String');
set(handles.loadTargetFileInput,'String',targetFilename);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function loadTargetFileInput_CreateFcn(hObject, ~, ~)
% hObject    handle to loadTargetFileInput (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadTargetFileButton.
function loadTargetFileButton_Callback(hObject, ~, handles)
% hObject    handle to loadTargetFileButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% Load a target minis+noise abf file:
[targetFilename, targetPathname, button] = uigetfile({'*.abf', 'Axon ABF files (*.abf)'}, 'Select a Target File', handles.ld);
if button
    targetFilename = fullfile(targetPathname, targetFilename);
    set(handles.loadTargetFileInput,'String',targetFilename);
    handles.ld = targetPathname;
else
    return
end

% Update handles structure
guidata(hObject, handles);


function loadNoiseFileInput_Callback(hObject, ~, handles)
% hObject    handle to loadNoiseFileInput (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

noiseFilename = get(hObject,'String');
set(handles.loadNoiseFileInput,'String',noiseFilename);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function loadNoiseFileInput_CreateFcn(hObject, ~, ~)
% hObject    handle to loadNoiseFileInput (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadNoiseFileButton.
function loadNoiseFileButton_Callback(hObject, ~, handles)
% hObject    handle to loadNoiseFileButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% Load a noise abf file:
[noiseFilename, noisePathname, button] = uigetfile({'*.abf', 'Axon ABF files (*.abf)'}, 'Noise file (TTX + gabazine + CPP + NBQX)', handles.ld);
if button
    noiseFilename = fullfile(noisePathname, noiseFilename);
    set(handles.loadNoiseFileInput,'String',noiseFilename);
    handles.ld = noisePathname;
else
    return
end

% Update handles structure
guidata(hObject, handles);










%% Detection Parameters panel:
function maxTimeToPeakEdit_Callback(hObject, ~, handles)
% hObject    handle to maxTimeToPeakEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

maxTimeToPeakEdit = str2double(get(hObject,'String'));
if isnan(maxTimeToPeakEdit)
    msgbox('Error: Maximum time to peak is NaN', 'Error', 'Error');
else
    set(handles.maxTimeToPeakEdit,'String',maxTimeToPeakEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxTimeToPeakEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxTimeToPeakEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function baselineDurationEdit_Callback(hObject, ~, handles)
% hObject    handle to baselineDurationEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

baselineDurationEdit = str2double(get(hObject,'String'));
if isnan(baselineDurationEdit)
    msgbox('Error: Baseline duration is NaN', 'Error', 'Error');
else
    set(handles.baselineDurationEdit,'String',baselineDurationEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function baselineDurationEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to baselineDurationEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function peakIntegrationPeriodEdit_Callback(hObject, ~, handles)
% hObject    handle to peakIntegrationPeriodEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

peakIntegrationPeriodEdit = str2double(get(hObject,'String'));
if isnan(peakIntegrationPeriodEdit)
    msgbox('Error: Peak integration period is NaN', 'Error', 'Error');
else
    set(handles.peakIntegrationPeriodEdit,'String',peakIntegrationPeriodEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function peakIntegrationPeriodEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to peakIntegrationPeriodEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function AmploboundEdit_Callback(hObject, ~, handles)
% hObject    handle to AmploboundEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

AmploboundEdit = str2double(get(hObject,'String'));
if isnan(AmploboundEdit)
    msgbox('Error: Amplitude lower bound is NaN', 'Error', 'Error');
else
    set(handles.AmploboundEdit,'String',AmploboundEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function AmploboundEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to AmploboundEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function AmpupboundEdit_Callback(hObject, ~, handles)
% hObject    handle to AmpupboundEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

AmpupboundEdit = str2double(get(hObject,'String'));
if isnan(AmpupboundEdit)
    msgbox('Error: Amplitude upper bound is NaN', 'Error', 'Error');
else
    set(handles.AmpupboundEdit,'String',AmpupboundEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function AmpupboundEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to AmpupboundEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function smoothWindowEdit_Callback(hObject, ~, handles)
% hObject    handle to smoothWindowEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

smoothWindowEdit = str2double(get(hObject,'String'));
if isnan(smoothWindowEdit)
    msgbox('Error: Gaussian smoothing window is NaN', 'Error', 'Error');
else
    set(handles.smoothWindowEdit,'String',smoothWindowEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function smoothWindowEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to smoothWindowEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in RTintEdit.
function RTintEdit_Callback(hObject, ~, handles)
% hObject    handle to RTintEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

value = get(hObject,'value');
set(handles.RTintEdit,'value',value);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function RTintEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to RTintEdit (see GCBO)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in RTbinSizeEdit.
function RTbinSizeEdit_Callback(hObject, ~, handles)
% hObject    handle to RTbinSizeEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

value = get(hObject,'value');
set(handles.RTbinSizeEdit,'value',value);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function RTbinSizeEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to RTbinSizeEdit (see GCBO)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function startPulseTargetEdit_Callback(hObject, ~, handles)
% hObject    handle to startPulseTargetEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

startPulseTargetStr = get(hObject,'String');
if isempty(startPulseTargetStr)
    startPulseTargetStr = '...';
    set(handles.startPulseTargetEdit,'String',startPulseTargetStr);
    guidata(hObject, handles);
    return
end
startPulseTargetStr = strrep(startPulseTargetStr, ' ', '');
startPulseTargetStrSplit = regexp(startPulseTargetStr, ',*', 'split');
startPulseTarget = strarray2numarray(startPulseTargetStrSplit');
if sum(isnan(startPulseTarget))
    errmsg = 'Error: At least one of the target pulse vectors is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.startPulseTargetEdit,'String',startPulseTargetStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function startPulseTargetEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to startPulseTargetEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function endPulseTargetEdit_Callback(hObject, ~, handles)
% hObject    handle to endPulseTargetEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

endPulseTargetStr = get(hObject,'String');
if isempty(endPulseTargetStr)
    endPulseTargetStr = '...';
    set(handles.endPulseTargetEdit,'String',endPulseTargetStr);
    guidata(hObject, handles);
    return
end
endPulseTargetStr = strrep(endPulseTargetStr, ' ', '');
endPulseTargetStrSplit = regexp(endPulseTargetStr, ',*', 'split');
endPulseTarget = strarray2numarray(endPulseTargetStrSplit');
if sum(isnan(endPulseTarget))
    errmsg = 'Error: At least one of the target pulse vectors is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.endPulseTargetEdit,'String',endPulseTargetStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function endPulseTargetEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to endPulseTargetEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function startGlitchTargetEdit_Callback(hObject, ~, handles)
% hObject    handle to startGlitchTargetEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

startGlitchTargetStr = get(hObject,'String');
if isempty(startGlitchTargetStr)
    startGlitchTargetStr = '...';
    set(handles.startGlitchTargetEdit,'String',startGlitchTargetStr);
    guidata(hObject, handles);
    return
end
startGlitchTargetStr = strrep(startGlitchTargetStr, ' ', '');
startGlitchTargetStrSplit = regexp(startGlitchTargetStr, ',*', 'split');
startGlitchTarget = strarray2numarray(startGlitchTargetStrSplit');
if sum(isnan(startGlitchTarget))
    errmsg = 'Error: At least one of the target glitch vectors is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.startGlitchTargetEdit,'String',startGlitchTargetStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function startGlitchTargetEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to startGlitchTargetEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function endGlitchTargetEdit_Callback(hObject, ~, handles)
% hObject    handle to endGlitchTargetEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

endGlitchTargetStr = get(hObject,'String');
if isempty(endGlitchTargetStr)
    endGlitchTargetStr = '...';
    set(handles.endGlitchTargetEdit,'String',endGlitchTargetStr);
    guidata(hObject, handles);
    return
end
endGlitchTargetStr = strrep(endGlitchTargetStr, ' ', '');
endGlitchTargetStrSplit = regexp(endGlitchTargetStr, ',*', 'split');
endGlitchTarget = strarray2numarray(endGlitchTargetStrSplit');
if sum(isnan(endGlitchTarget))
    errmsg = 'Error: At least one of the target glitch vectors is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.endGlitchTargetEdit,'String',endGlitchTargetStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function endGlitchTargetEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to endGlitchTargetEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function startPulseNoiseEdit_Callback(hObject, ~, handles)
% hObject    handle to startPulseNoiseEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

startPulseNoiseStr = get(hObject,'String');
if isempty(startPulseNoiseStr)
    startPulseNoiseStr = '...';
    set(handles.startPulseNoiseEdit,'String',startPulseNoiseStr);
    guidata(hObject, handles);
    return
end
startPulseNoiseStr = strrep(startPulseNoiseStr, ' ', '');
startPulseNoiseStrSplit = regexp(startPulseNoiseStr, ',*', 'split');
startPulseNoise = strarray2numarray(startPulseNoiseStrSplit');
if sum(isnan(startPulseNoise))
    errmsg = 'Error: At least one of the noise pulse vectors is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.startPulseNoiseEdit,'String',startPulseNoiseStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function startPulseNoiseEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to startPulseNoiseEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function endPulseNoiseEdit_Callback(hObject, ~, handles)
% hObject    handle to endPulseNoiseEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

endPulseNoiseStr = get(hObject,'String');
if isempty(endPulseNoiseStr)
    endPulseNoiseStr = '...';
    set(handles.endPulseNoiseEdit,'String',endPulseNoiseStr);
    guidata(hObject, handles);
    return
end
endPulseNoiseStr = strrep(endPulseNoiseStr, ' ', '');
endPulseNoiseStrSplit = regexp(endPulseNoiseStr, ',*', 'split');
endPulseNoise = strarray2numarray(endPulseNoiseStrSplit');
if sum(isnan(endPulseNoise))
    errmsg = 'Error: At least one of the noise pulse vectors is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.endPulseNoiseEdit,'String',endPulseNoiseStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function endPulseNoiseEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to endPulseNoiseEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function startGlitchNoiseEdit_Callback(hObject, ~, handles)
% hObject    handle to startGlitchNoiseEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

startGlitchNoiseStr = get(hObject,'String');
if isempty(startGlitchNoiseStr)
    startGlitchNoiseStr = '...';
    set(handles.startGlitchNoiseEdit,'String',startGlitchNoiseStr);
    guidata(hObject, handles);
    return
end
startGlitchNoiseStr = strrep(startGlitchNoiseStr, ' ', '');
startGlitchNoiseStrSplit = regexp(startGlitchNoiseStr, ',*', 'split');
startGlitchNoise = strarray2numarray(startGlitchNoiseStrSplit');
if sum(isnan(startGlitchNoise))
    errmsg = 'Error: At least one of the noise glitch vectors is NaN';
    msgbox(errmsg,'Error','Error');
    return
else
    set(handles.startGlitchNoiseEdit,'String',startGlitchNoiseStr);
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function startGlitchNoiseEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to startGlitchNoiseEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function endGlitchNoiseEdit_Callback(hObject, ~, handles)
% hObject    handle to endGlitchNoiseEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

endGlitchNoiseStr = get(hObject,'String');
if isempty(endGlitchNoiseStr)
    endGlitchNoiseStr = '...';
    set(handles.endGlitchNoiseEdit,'String',endGlitchNoiseStr);
    guidata(hObject, handles);
    return
end
endGlitchNoiseStr = strrep(endGlitchNoiseStr, ' ', '');
endGlitchNoiseStrSplit = regexp(endGlitchNoiseStr, ',*', 'split');
endGlitchNoise = strarray2numarray(endGlitchNoiseStrSplit');
if sum(isnan(endGlitchNoise))
    errmsg = 'Error: At least one of the noise glitch vectors is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.endGlitchNoiseEdit,'String',endGlitchNoiseStr);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function endGlitchNoiseEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to endGlitchNoiseEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pulseDurationEdit_Callback(hObject, ~, handles)
% hObject    handle to pulseDurationEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

pulseDurationEdit = get(hObject, 'String');
if isempty(pulseDurationEdit) || strcmpi(pulseDurationEdit, '...')
    pulseDurationEdit = '...';
    set(handles.pulseDurationEdit, 'String', pulseDurationEdit);
    guidata(hObject, handles);
    return
end

pulseDurationEdit = str2double(get(hObject,'String'));
if isnan(pulseDurationEdit)
    msgbox('Error: Brief (second) pulse duration is NaN', 'Error', 'Error');
else
    set(handles.pulseDurationEdit,'String',pulseDurationEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function pulseDurationEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to pulseDurationEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in downGoingCheckbox.
function downGoingCheckbox_Callback(hObject, ~, handles)
% hObject    handle to downGoingCheckbox (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in voltageClampCheckbox.
function voltageClampCheckbox_Callback(hObject, ~, handles)
% hObject    handle to voltageClampCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update handles structure
guidata(hObject, handles);












%% Optimisation Parameters panel:
% --- Executes on selection change in distTypeEdit.
function distTypeEdit_Callback(hObject, ~, handles)
% hObject    handle to distTypeEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

value = get(handles.distTypeEdit,'Value');
set(handles.distTypeEdit,'Value',value);
distType = get(handles.distTypeEdit,'String');
distType = distType(value);
if strcmpi(distType,'Log Normal') || strcmpi(distType,'Bimodal Log Normal') || strcmpi(distType,'Trimodal Log Normal')
    %          mu1_1       A_1  alpha_1 scale_1     mu2_1    B_1 beta_1    rho_1    mu1_2    A_2   alpha_2 scale_2  mu2_2      B_2   beta_2  rho_2     mu1_3    A_3  alpha_3  scale_3    mu2_3      B_3   beta_3    rho_3
    loBounds = [-0.1,      .01,     .01,      0,       -1,   .01,   .01,      -1,    -0.1,   .01,      .01, -10000,    -1,     .01,     .01,    -1,     -0.1,   .01,     .01,  -10000,      -1,     .01,     .01,      -1];
    upBounds = [   1,        1,     1.5,  10000,       10,    20,   1.5,       1,       1,     1,      1.5,  10000,    10,      20,     1.5,     1,        1,     1,     1.5,   10000,      10,      20,     1.5,       1];
elseif ~strcmpi(distType,'Skew-normal') && ~strcmpi(distType,'Bimodal skew-normal') && ~strcmpi(distType,'Trimodal skew-normal')
    %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  mu1_2  sigma1_2  scale_2  mu2_2  sigma2_2  rho_2  mu1_3  sigma1_3  scale_3  mu2_3  sigma2_3  rho_3    mu1_4  sigma1_4  scale_4  mu2_4  sigma2_4    rho_4
    loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,   -0.1,      .01,  -10000,    -1,      .01,      -1];
    upBounds = [   1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,   	  1,        1,   10000,    10,       10,       1];
elseif strcmpi(distType,'Skew-normal')
    %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  skew1_1  skew2_1
    loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,     -10,     -10, -10000,    -1,      .01,    -1,   .01,      .01,  -10000,    -1,      .01,    -1,     -10,     -10,      -10,   -10,      -10,     -10];
    upBounds = [   1,        1,   10000,    10,       10,     1,      10,      10,  10000,    10,       10,     1,     2,        1,   10000,    10,       10,     1,      10,      10,       10,    10,       10,      10];
elseif ~strcmpi(distType,'Trimodal skew-normal')
    %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  mu1_2  sigma1_2  scale_2  mu2_2  sigma2_2  rho_2  skew1_1  skew2_1  skew1_2  skew2_2
    loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,     -10,     -10,     -10,     -10,   .01,    -1,     -10,     -10,      -10,   -10,      -10,     -10];
    upBounds = [   1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,      10,      10,      10,      10,    10,     1,      10,      10,       10,    10,       10,      10];
else
    %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  mu1_2  sigma1_2  scale_2  mu2_2  sigma2_2  rho_2  mu1_3  sigma1_3  scale_3  mu2_3  sigma2_3  rho_3  skew1_1  skew2_1  skew1_2  skew2_2  skew1_3  skew2_3
    loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,     -10,     -10,     -10,     -10,     -10,     -10];
    upBounds = [   1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,      10,      10,      10,      10,      10,      10];
end
handles.options.bounds = [loBounds; upBounds];

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function distTypeEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to distTypeEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in distBaselineEdit.
function distBaselineEdit_Callback(hObject, ~, handles)
% hObject    handle to distBaselineEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

value = get(hObject,'value');
set(handles.distBaselineEdit,'value',value);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function distBaselineEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to distBaselineEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SDupbound_Callback(hObject, ~, handles)
% hObject    handle to SDupbound (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

SDupbound = get(hObject, 'string');
if (isempty(SDupbound) || strcmpi(SDupbound, '...')) && ~handles.options.SDupbound
    SDupbound = '...';
    set(handles.SDupbound, 'String', SDupbound);
    guidata(hObject, handles);
    return
end

SDupbound = str2double(SDupbound);
if isnan(SDupbound)
    msgbox('Error: Standard deviation upper bound is NaN', 'Error', 'Error');
else
    set(handles.SDupbound, 'String', SDupbound);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function SDupbound_CreateFcn(hObject, ~, ~)
% hObject    handle to SDupbound (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SDlobound_Callback(hObject, ~, handles)
% hObject    handle to SDlobound (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

SDlobound = get(hObject,'String');
if (isempty(SDlobound) || strcmpi(SDlobound, '...')) && ~handles.options.SDlobound
    SDlobound = '...';
    set(handles.SDlobound, 'String' ,SDlobound);
    guidata(hObject, handles);
    return
end

SDlobound = str2double(SDlobound);
if isnan(SDlobound)
    msgbox('Error: Standard deviation lower bound is NaN', 'Error', 'Error');
else
    set(handles.SDlobound, 'String', SDlobound);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function SDlobound_CreateFcn(hObject, ~, ~)
% hObject    handle to SDlobound (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxErrEdit_Callback(hObject, ~, handles)
% hObject    handle to maxErrEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

maxErrEdit = str2double(get(hObject,'String'));
if isnan(maxErrEdit)
    msgbox('Error: Maximum combined SAD is NaN','Error','Error');
else
    set(handles.maxErrEdit,'String', maxErrEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxErrEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxErrEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxAmpErrEdit_Callback(hObject, ~, handles)
% hObject    handle to maxAmpErrEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

maxAmpErrEdit = str2double(get(hObject,'String'));
if isnan(maxAmpErrEdit)
    msgbox('Error: Maximum amplitude SAD is NaN','Error','Error');
else
    set(handles.maxAmpErrEdit,'String', maxAmpErrEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxAmpErrEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxAmpErrEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxRTErrEdit_Callback(hObject, ~, handles)
% hObject    handle to maxRTErrEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

maxRTErrEdit = str2double(get(hObject,'String'));
if isnan(maxRTErrEdit)
    msgbox('Error: Maximum rise time SAD is NaN','Error','Error');
else
    set(handles.maxRTErrEdit,'String', maxRTErrEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxRTErrEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxRTErrEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxDevEdit_Callback(hObject, ~, handles)
% hObject    handle to maxDevEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

maxDevEdit = str2double(get(hObject,'String'));
if isnan(maxDevEdit)
    msgbox('Error: Maximum combined MAD is NaN','Error','Error');
else
    set(handles.maxDevEdit,'String', maxDevEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxDevEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxDevEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxDevAmpEdit_Callback(hObject, ~, handles)
% hObject    handle to maxDevAmpEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

maxDevAmpEdit = str2double(get(hObject,'String'));
if isnan(maxDevAmpEdit)
    msgbox('Error: Maximum amplitude MAD is NaN','Error','Error');
else
    set(handles.maxDevAmpEdit,'String', maxDevAmpEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxDevAmpEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxDevAmpEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxDevRTEdit_Callback(hObject, ~, handles)
% hObject    handle to maxDevRTEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

maxDevRTEdit = str2double(get(hObject,'String'));
if isnan(maxDevRTEdit)
    msgbox('Error: Maximum rise time MAD is NaN','Error','Error');
else
    set(handles.maxDevRTEdit,'String', maxDevRTEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxDevRTEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxDevRTEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxAmpBottomErrEdit_Callback(hObject, ~, handles)
% hObject    handle to maxAmpBottomErrEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

maxAmpBottomErrEdit = str2double(get(hObject,'String'));
if isnan(maxAmpBottomErrEdit)
    msgbox('Error: Maximum amplitude top 50% SAD is NaN','Error','Error');
else
    set(handles.maxAmpBottomErrEdit,'String', maxAmpBottomErrEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxAmpBottomErrEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxAmpBottomErrEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxAmpBottomDevEdit_Callback(hObject, ~, handles)
% hObject    handle to maxAmpBottomDevEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

maxAmpBottomDevEdit = str2double(get(hObject,'String'));
if isnan(maxAmpBottomDevEdit)
    msgbox('Error: Maximum amplitude top 50% MAD is NaN','Error','Error');
else
    set(handles.maxAmpBottomDevEdit,'String', maxAmpBottomDevEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxAmpBottomDevEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxAmpBottomDevEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxAmpMidErrEdit_Callback(hObject, ~, handles)
% hObject    handle to maxAmpMidErrEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

maxAmpMidErrEdit = str2double(get(hObject,'String'));
if isnan(maxAmpMidErrEdit)
    msgbox('Error: Maximum amplitude top 10% SAD is NaN','Error','Error');
else
    set(handles.maxAmpMidErrEdit,'String', maxAmpMidErrEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxAmpMidErrEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxAmpMidErrEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxAmpMidDevEdit_Callback(hObject, ~, handles)
% hObject    handle to maxAmpMidDevEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

maxAmpMidDevEdit = str2double(get(hObject,'String'));
if isnan(maxAmpMidDevEdit)
    msgbox('Error: Maximum amplitude top 10% MAD is NaN','Error','Error');
else
    set(handles.maxAmpMidDevEdit,'String', maxAmpMidDevEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxAmpMidDevEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxAmpMidDevEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxAmpTopErrEdit_Callback(hObject, ~, handles)
% hObject    handle to maxAmpTopErrEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

maxAmpTopErrEdit = str2double(get(hObject,'String'));
if isnan(maxAmpTopErrEdit)
    msgbox('Error: Maximum amplitude top 2% SAD is NaN','Error','Error');
else
    set(handles.maxAmpTopErrEdit,'String', maxAmpTopErrEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxAmpTopErrEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxAmpTopErrEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxAmpTopDevEdit_Callback(hObject, ~, handles)
% hObject    handle to maxAmpTopDevEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

maxAmpTopDevEdit = str2double(get(hObject,'String'));
if isnan(maxAmpTopDevEdit)
    msgbox('Error: Maximum amplitude top 2% MAD is NaN','Error','Error');
else
    set(handles.maxAmpTopDevEdit,'String', maxAmpTopDevEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function maxAmpTopDevEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to maxAmpTopDevEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end










%% Simulation Parameters panel:
function loSimAmpEdit_Callback(hObject, ~, handles)
% hObject    handle to loSimAmpEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loSimAmpEdit = str2double(get(hObject,'String'));
if isnan(loSimAmpEdit)
    msgbox('Error: Amplitude lower bound is NaN', 'Error', 'Error');
else
    set(handles.loSimAmpEdit,'String',loSimAmpEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loSimAmpEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to loSimAmpEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function LEdit_Callback(hObject, ~, handles)
% hObject    handle to LEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

LEdit = str2double(get(hObject,'String'));
if isnan(LEdit)
    msgbox('Error: Initial L (electrotonic length) is NaN','Error','Error');
else
    set(handles.LEdit,'String',LEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function LEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to LEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function tau_mEdit_Callback(hObject, ~, handles)
% hObject    handle to nSeriesEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

tau_mEdit = get(hObject, 'String');
if isempty(tau_mEdit) || strcmpi(tau_mEdit, '...')
    tau_mEdit = '...';
    set(handles.tau_mEdit, 'String', tau_mEdit);
    guidata(hObject, handles);
    return
end

tau_mEdit = str2double(tau_mEdit);
if isnan(tau_mEdit)
    msgbox('Error: tau_m (passive membrane time constant; lower limit) is NaN','Error','Error');
else
    set(handles.tau_mEdit, 'String', tau_mEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function tau_mEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to nSeriesEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function tau_PSPmEdit_Callback(hObject, ~, handles)
% hObject    handle to tau_PSPmEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

tau_PSPmEdit = get(hObject, 'String');
if isempty(tau_PSPmEdit) || strcmpi(tau_PSPmEdit, '...')
    tau_PSPmEdit = '...';
    set(handles.tau_PSPmEdit, 'String', tau_PSPmEdit);
    guidata(hObject, handles);
    return
end

tau_PSPmEdit = str2double(tau_PSPmEdit);
if isnan(tau_PSPmEdit)
    msgbox('Error: tau_m (passive membrane time constant; upper limit) is NaN','Error','Error');
else
    set(handles.tau_PSPmEdit, 'String', tau_PSPmEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function tau_PSPmEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to tau_PSPmEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end










%% Settings menu:
% --------------------------------------------------------------------
function settingsMenu_Callback(~, ~, ~)


% --------------------------------------------------------------------
function loadSettings_Callback(hObject, ~, handles)

[settingsFilename, settingsPathname, filterIndex] = uigetfile({'*.mat', 'Matlab MAT file (*.mat)'}, 'Load settings', handles.ld);
if filterIndex
    try
        settFullName = fullfile(settingsPathname, settingsFilename);
        settings = load(settFullName);
        settings = settings.settings;
        handles.ld = settingsPathname;
    catch %#ok<CTCH>
        msgbox('Some parameters were not loaded properly because your settings file might be corrupt.','Error','Error');
        fclose all;
        return
    end
else
    return
end

try
    % File panel:
    set(handles.loadTargetFileInput, 'string', settings.loadTargetFileInput);
    set(handles.loadNoiseFileInput, 'string', settings.loadNoiseFileInput);
    
    % Detection Parameters panel:
    set(handles.maxTimeToPeakEdit, 'string', settings.maxTimeToPeakEdit);
    set(handles.baselineDurationEdit, 'string', settings.baselineDurationEdit);
    set(handles.peakIntegrationPeriodEdit, 'string', settings.peakIntegrationPeriodEdit);
    set(handles.AmploboundEdit, 'string', settings.AmploboundEdit);
    set(handles.AmpupboundEdit, 'string', settings.AmpupboundEdit);
    set(handles.smoothWindowEdit, 'string', settings.smoothWindowEdit);
    set(handles.RTintEdit, 'value', settings.RTintEdit);
    set(handles.RTbinSizeEdit, 'value', settings.RTbinSizeEdit);
    set(handles.startPulseTargetEdit, 'string', settings.startPulseTargetEdit);
    set(handles.endPulseTargetEdit, 'string', settings.endPulseTargetEdit);
    set(handles.startGlitchTargetEdit, 'string', settings.startGlitchTargetEdit);
    set(handles.endGlitchTargetEdit, 'string', settings.endGlitchTargetEdit);
    set(handles.startPulseNoiseEdit, 'string', settings.startPulseNoiseEdit);
    set(handles.endPulseNoiseEdit,'string',settings.endPulseNoiseEdit);
    set(handles.startGlitchNoiseEdit, 'string', settings.startGlitchNoiseEdit);
    set(handles.endGlitchNoiseEdit, 'string', settings.endGlitchNoiseEdit);
    set(handles.pulseDurationEdit, 'string', settings.pulseDurationEdit);
    set(handles.downGoingCheckbox, 'value', settings.downGoingCheckbox);
    set(handles.voltageClampCheckbox, 'value', settings.voltageClampCheckbox);
    
    % Optimisation Parameters panel:
    set(handles.distTypeEdit, 'value', settings.distTypeEdit);
    set(handles.distBaselineEdit, 'value', settings.distBaselineEdit);
    set(handles.SDlobound, 'string', settings.SDlobound);
    set(handles.SDupbound, 'string', settings.SDupbound);
    set(handles.maxErrEdit, 'string', settings.maxErrEdit);
    set(handles.maxAmpErrEdit, 'string', settings.maxAmpErrEdit);
    set(handles.maxRTErrEdit, 'string', settings.maxRTErrEdit);
    set(handles.maxDevEdit, 'string', settings.maxDevEdit);
    set(handles.maxDevAmpEdit, 'string', settings.maxDevAmpEdit);
    set(handles.maxDevRTEdit, 'string', settings.maxDevRTEdit);
    set(handles.maxAmpBottomErrEdit, 'string', settings.maxAmpBottomErrEdit);
    set(handles.maxAmpBottomDevEdit, 'String', settings.maxAmpBottomDevEdit);
    set(handles.maxAmpMidErrEdit, 'string', settings.maxAmpMidErrEdit);
    set(handles.maxAmpMidDevEdit, 'String', settings.maxAmpMidDevEdit);
    set(handles.maxAmpTopErrEdit, 'string', settings.maxAmpTopErrEdit);
    set(handles.maxAmpTopDevEdit,'String', settings.maxAmpTopDevEdit);
    
    % Simulation Parameters panel:
    set(handles.loSimAmpEdit, 'string', settings.loSimAmpEdit);
    set(handles.LEdit, 'string', settings.LEdit);
    set(handles.tau_mEdit, 'string', settings.tau_mEdit);
    set(handles.tau_PSPmEdit, 'string', settings.tau_PSPmEdit);
    
    handles.options = settings.options;
    
    parallelCores = str2double(handles.options.parallelCores);
    pool = gcp('nocreate');
    if isempty(pool)
        openWorkers = 0;
    else
        openWorkers = pool.NumWorkers;
    end
    scheduler = parcluster('local');
    if parallelCores > scheduler.NumWorkers
        parallelCores = scheduler.NumWorkers;
    end
    if openWorkers ~= parallelCores
        delete(gcp('nocreate'));
        if parallelCores > 1
            parpool(parallelCores);
            handles.options.parallelCores = parallelCores;
        end
    end
catch %#ok<CTCH>
    msgbox('Some parameters were not loaded properly because your settings file might be corrupt.','Error','Error');
end
fclose all;

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function saveSettings_Callback(hObject, ~, handles)

% File panel:
settings.loadTargetFileInput = get(handles.loadTargetFileInput, 'string');
settings.loadNoiseFileInput = get(handles.loadNoiseFileInput, 'string');

% Detection Parameters panel:
settings.maxTimeToPeakEdit = get(handles.maxTimeToPeakEdit, 'string');
settings.baselineDurationEdit = get(handles.baselineDurationEdit, 'string');
settings.peakIntegrationPeriodEdit = get(handles.peakIntegrationPeriodEdit, 'string');
settings.AmploboundEdit = get(handles.AmploboundEdit, 'string');
settings.AmpupboundEdit = get(handles.AmpupboundEdit, 'string');
settings.smoothWindowEdit = get(handles.smoothWindowEdit, 'string');
settings.RTintEdit = get(handles.RTintEdit, 'value');
settings.RTbinSizeEdit = get(handles.RTbinSizeEdit, 'value');
settings.startPulseTargetEdit = get(handles.startPulseTargetEdit, 'string');
settings.endPulseTargetEdit = get(handles.endPulseTargetEdit, 'string');
settings.startGlitchTargetEdit = get(handles.startGlitchTargetEdit, 'string');
settings.endGlitchTargetEdit = get(handles.endGlitchTargetEdit, 'string');
settings.startPulseNoiseEdit = get(handles.startPulseNoiseEdit, 'string');
settings.endPulseNoiseEdit = get(handles.endPulseNoiseEdit, 'string');
settings.startGlitchNoiseEdit = get(handles.startGlitchNoiseEdit, 'string');
settings.endGlitchNoiseEdit = get(handles.endGlitchNoiseEdit, 'string');
settings.pulseDurationEdit = get(handles.pulseDurationEdit, 'string');
settings.downGoingCheckbox = get(handles.downGoingCheckbox, 'value');
settings.voltageClampCheckbox = get(handles.voltageClampCheckbox, 'value');

% Optimisation Parameters panel:
settings.distTypeEdit = get(handles.distTypeEdit, 'value');
settings.distBaselineEdit = get(handles.distBaselineEdit, 'value');
settings.SDlobound = get(handles.SDlobound, 'string');
settings.SDupbound = get(handles.SDupbound, 'string');
settings.maxErrEdit = get(handles.maxErrEdit, 'string');
settings.maxAmpErrEdit = get(handles.maxAmpErrEdit, 'string');
settings.maxRTErrEdit = get(handles.maxRTErrEdit, 'string');
settings.maxDevEdit = get(handles.maxDevEdit, 'string');
settings.maxDevAmpEdit = get(handles.maxDevAmpEdit, 'string');
settings.maxDevRTEdit = get(handles.maxDevRTEdit, 'string');
settings.maxAmpBottomErrEdit = get(handles.maxAmpBottomErrEdit, 'string');
settings.maxAmpBottomDevEdit = get(handles.maxAmpBottomDevEdit, 'String');
settings.maxAmpMidErrEdit = get(handles.maxAmpMidErrEdit, 'string');
settings.maxAmpMidDevEdit = get(handles.maxAmpMidDevEdit, 'String');
settings.maxAmpTopErrEdit = get(handles.maxAmpTopErrEdit, 'string');
settings.maxAmpTopDevEdit = get(handles.maxAmpTopDevEdit, 'String');

% Simulation Parameters panel:
settings.loSimAmpEdit = get(handles.loSimAmpEdit, 'string');
settings.LEdit = get(handles.LEdit, 'string');
settings.tau_mEdit = get(handles.tau_mEdit, 'string');
settings.tau_PSPmEdit = get(handles.tau_PSPmEdit, 'string');

settings.options = handles.options;

pathname = handles.ld;
[settingsFilename, settingsPathname, filterIndex] = uiputfile({'*.mat', 'Matlab MAT file (*.mat)'}, 'Save settings as', pathname);
if filterIndex
    save(fullfile(settingsPathname,settingsFilename), 'settings');
    handles.ld = settingsPathname;
end
fclose all;

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function resetSettings_Callback(hObject, ~, handles)
% handles    structure with handles and user data (see GUIDATA)

if fafTest
    % Load settings file:
    settings = load('fafSettings.mat');
    settings = settings.settings;
    
    % File panel:
    set(handles.loadTargetFileInput, 'string', settings.loadTargetFileInput);
    set(handles.loadNoiseFileInput, 'string', settings.loadNoiseFileInput);

    % Detection Parameters panel:
    set(handles.maxTimeToPeakEdit, 'string', settings.maxTimeToPeakEdit);
    set(handles.baselineDurationEdit, 'string', settings.baselineDurationEdit);
    set(handles.peakIntegrationPeriodEdit, 'string', settings.peakIntegrationPeriodEdit);
    set(handles.AmploboundEdit, 'string', settings.AmploboundEdit);
    set(handles.AmpupboundEdit, 'string', settings.AmpupboundEdit);
    set(handles.smoothWindowEdit, 'string', settings.smoothWindowEdit);
    set(handles.RTintEdit, 'value', settings.RTintEdit);
    set(handles.RTbinSizeEdit, 'value', settings.RTbinSizeEdit);
    set(handles.startPulseTargetEdit, 'string', settings.startPulseTargetEdit);
    set(handles.endPulseTargetEdit, 'string', settings.endPulseTargetEdit);
    set(handles.startGlitchTargetEdit, 'string', settings.startGlitchTargetEdit);
    set(handles.endGlitchTargetEdit, 'string', settings.endGlitchTargetEdit);
    set(handles.startPulseNoiseEdit, 'string', settings.startPulseNoiseEdit);
    set(handles.endPulseNoiseEdit,'string',settings.endPulseNoiseEdit);
    set(handles.startGlitchNoiseEdit, 'string', settings.startGlitchNoiseEdit);
    set(handles.endGlitchNoiseEdit, 'string', settings.endGlitchNoiseEdit);
    set(handles.pulseDurationEdit, 'string', settings.pulseDurationEdit);
    set(handles.downGoingCheckbox, 'value', settings.downGoingCheckbox);
    set(handles.voltageClampCheckbox, 'value', settings.voltageClampCheckbox);

    % Optimisation Parameters panel:
    set(handles.distTypeEdit, 'value', settings.distTypeEdit);
    set(handles.distBaselineEdit, 'value', settings.distBaselineEdit);
    set(handles.SDlobound, 'string', settings.SDlobound);
    set(handles.SDupbound, 'string', settings.SDupbound);
    set(handles.maxErrEdit, 'string', settings.maxErrEdit);
    set(handles.maxAmpErrEdit, 'string', settings.maxAmpErrEdit);
    set(handles.maxRTErrEdit, 'string', settings.maxRTErrEdit);
    set(handles.maxDevEdit, 'string', settings.maxDevEdit);
    set(handles.maxDevAmpEdit, 'string', settings.maxDevAmpEdit);
    set(handles.maxDevRTEdit, 'string', settings.maxDevRTEdit);
    set(handles.maxAmpBottomErrEdit, 'string', settings.maxAmpBottomErrEdit);
    set(handles.maxAmpBottomDevEdit, 'String', settings.maxAmpBottomDevEdit);
    set(handles.maxAmpMidErrEdit, 'string', settings.maxAmpMidErrEdit);
    set(handles.maxAmpMidDevEdit, 'String', settings.maxAmpMidDevEdit);
    set(handles.maxAmpTopErrEdit, 'string', settings.maxAmpTopErrEdit);
    set(handles.maxAmpTopDevEdit,'String', settings.maxAmpTopDevEdit);

    % Simulation Parameters panel:
    set(handles.loSimAmpEdit, 'string', settings.loSimAmpEdit);
    set(handles.LEdit, 'string', settings.LEdit);
    set(handles.tau_mEdit, 'string', settings.tau_mEdit);
    set(handles.tau_PSPmEdit, 'string', settings.tau_PSPmEdit);

    handles.options = settings.options;

    parallelCores = str2double(handles.options.parallelCores);
    pool = gcp('nocreate');
    if isempty(pool)
        openWorkers = 0;
    else
        openWorkers = pool.NumWorkers;
    end
    scheduler = parcluster('local');
    if parallelCores > scheduler.NumWorkers
      parallelCores = scheduler.NumWorkers;
    end
    if openWorkers ~= parallelCores
      delete(gcp('nocreate'));
      if parallelCores > 1
        parpool(parallelCores);
        handles.options.parallelCores = parallelCores;
      end
    end
else
    % Detection Parameters panel:
    set(handles.maxTimeToPeakEdit, 'string', 10);
    set(handles.baselineDurationEdit, 'string', 2);
    set(handles.peakIntegrationPeriodEdit, 'string', 2.5);
    set(handles.AmploboundEdit, 'string', 0.02);
    set(handles.AmpupboundEdit, 'string', 10);
    set(handles.smoothWindowEdit, 'string', 1.5);
    set(handles.RTintEdit, 'value', 1);
    set(handles.RTbinSizeEdit, 'value', 1);
    set(handles.startPulseTargetEdit, 'string', '...');
    set(handles.endPulseTargetEdit, 'string', '...');
    set(handles.startGlitchTargetEdit, 'string', '...');
    set(handles.endGlitchTargetEdit, 'string', '...');
    set(handles.startPulseNoiseEdit, 'string', '...');
    set(handles.endPulseNoiseEdit, 'string', '...');
    set(handles.startGlitchNoiseEdit, 'string', '...');
    set(handles.endGlitchNoiseEdit, 'string', '...');
    set(handles.pulseDurationEdit, 'string', '0.5');
    set(handles.downGoingCheckbox, 'value', 0);
    set(handles.voltageClampCheckbox, 'value', 0);

    % Optimisation Parameters panel:
    set(handles.distTypeEdit, 'value', 4);
    set(handles.distBaselineEdit, 'value', 1);
    set(handles.SDlobound, 'string', 0.025);
    set(handles.SDupbound, 'string', 0.035);
    set(handles.maxErrEdit, 'string', 2000);
    set(handles.maxAmpErrEdit, 'string', 1000);
    set(handles.maxRTErrEdit, 'string', 1000);
    set(handles.maxDevEdit, 'string', 150);
    set(handles.maxDevAmpEdit, 'string', 300);
    set(handles.maxDevRTEdit, 'string', 300);
    set(handles.maxAmpBottomErrEdit, 'string', 1000);
    set(handles.maxAmpBottomDevEdit, 'String', 100);
    set(handles.maxAmpMidErrEdit, 'string', 1000);
    set(handles.maxAmpMidDevEdit, 'String', 100);
    set(handles.maxAmpTopErrEdit, 'string', 400);
    set(handles.maxAmpTopDevEdit, 'String', 40);

    % Simulation Parameters panel:
    set(handles.loSimAmpEdit, 'string', 0.01);
    set(handles.LEdit, 'string', 0.6);
    set(handles.tau_mEdit, 'string', 10);
    set(handles.tau_PSPmEdit, 'string', 20);

    % Restore default options:
    %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  mu1_2  sigma1_2  scale_2  mu2_2  sigma2_2  rho_2  mu1_3  sigma1_3  scale_3  mu2_3  sigma2_3  rho_3    mu1_4  sigma1_4  scale_4  mu2_4  sigma2_4    rho_4
    loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,    -0.1,      .01,  -10000,    -1,      .01,      -1];
    upBounds = [   1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,   	   1,        1,   10000,    10,       10,       1];
    bounds = [loBounds; upBounds];

    pool = gcp('nocreate');
    if isempty(pool)
        openWorkers = 0;
    else
        openWorkers = pool.NumWorkers;
    end
    scheduler = parcluster('local');
    if openWorkers ~= scheduler.NumWorkers
        delete(gcp('nocreate'));
        pool = parpool;
    end
    parallelCores = pool.NumWorkers;

    handles.options = struct('bounds', bounds, 'nGenerations', 50, 'parallelCores', parallelCores, 'fullParallel', 0, 'tauRange', 0,...
        'figureDisplay', 1, 'cluster', 0, 'clusterProfile', 'local', 'cliff', 0, 'SDlobound', 0, 'SDupbound', 0);
end

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function clearSettings_Callback(hObject, ~, handles)
% handles    structure with handles and user data (see GUIDATA)

% Detection Parameters panel:
set(handles.maxTimeToPeakEdit, 'string', '');
set(handles.baselineDurationEdit, 'string', '');
set(handles.peakIntegrationPeriodEdit, 'string', '');
set(handles.AmploboundEdit, 'string', '');
set(handles.AmpupboundEdit, 'string', '');
set(handles.smoothWindowEdit, 'string', '');
set(handles.RTintEdit, 'value', 1);
set(handles.RTbinSizeEdit, 'value', 1);
set(handles.startPulseTargetEdit, 'string', '...');
set(handles.endPulseTargetEdit, 'string', '...');
set(handles.startGlitchTargetEdit, 'string', '...');
set(handles.endGlitchTargetEdit, 'string', '...');
set(handles.startPulseNoiseEdit, 'string', '...');
set(handles.endPulseNoiseEdit, 'string', '...');
set(handles.startGlitchNoiseEdit, 'string', '...');
set(handles.endGlitchNoiseEdit, 'string', '...');
set(handles.pulseDurationEdit, 'string', '');
set(handles.downGoingCheckbox, 'value', 0);
set(handles.voltageClampCheckbox, 'value', 0);

% Optimisation Parameters panel:
set(handles.distTypeEdit, 'value', 1);
set(handles.distBaselineEdit, 'value', 1);
set(handles.SDlobound, 'string', '');
set(handles.SDupbound, 'string', '');
set(handles.maxErrEdit, 'string', '');
set(handles.maxAmpErrEdit, 'string', '');
set(handles.maxRTErrEdit, 'string', '');
set(handles.maxDevEdit, 'string', '');
set(handles.maxDevAmpEdit, 'string', '');
set(handles.maxDevRTEdit, 'string', '');
set(handles.maxAmpBottomErrEdit, 'string', '');
set(handles.maxAmpBottomDevEdit, 'String', '');
set(handles.maxAmpMidErrEdit, 'string', '');
set(handles.maxAmpMidDevEdit, 'String', '');
set(handles.maxAmpTopErrEdit, 'string', '');
set(handles.maxAmpTopDevEdit, 'String', '');

% Simulation Parameters panel:
set(handles.loSimAmpEdit, 'string', '');
set(handles.LEdit, 'string', '');
set(handles.tau_mEdit, 'string', '');
set(handles.tau_PSPmEdit, 'string', '');

% Update handles structure
guidata(hObject, handles);










%% Options menu:
% --------------------------------------------------------------------
function optimisationOptionsMenu_Callback(~, ~, ~)


% --------------------------------------------------------------------
function changeOptions_Callback(hObject, ~, handles)
% handles    structure with handles and user data (see GUIDATA)

distributionTypeString = get(handles.distTypeEdit, 'string');
distributionTypeValue = get(handles.distTypeEdit, 'value');
distributionType = distributionTypeString(distributionTypeValue);

if ~isfield(handles.options, 'bounds')
    if strcmpi(distType,'Log Normal') || strcmpi(distType,'Bimodal Log Normal') || strcmpi(distType,'Trimodal Log Normal')
        %          mu1_1       A_1  alpha_1 scale_1     mu2_1    B_1 beta_1    rho_1    mu1_2    A_2   alpha_2 scale_2  mu2_2      B_2   beta_2  rho_2     mu1_3    A_3  alpha_3  scale_3    mu2_3      B_3   beta_3    rho_3
        loBounds = [-0.1,      .01,     .01,      0,       -1,   .01,   .01,      -1,    -0.1,   .01,      .01, -10000,    -1,     .01,     .01,    -1,     -0.1,   .01,     .01,  -10000,      -1,     .01,     .01,      -1];
        upBounds = [   1,        1,     1.5,  10000,       10,    20,   1.5,       1,       1,     1,      1.5,  10000,    10,      20,     1.5,     1,        1,     1,     1.5,   10000,      10,      20,     1.5,       1];
    elseif ~strcmpi(distType,'Skew-normal') && ~strcmpi(distType,'Bimodal skew-normal') && ~strcmpi(distType,'Trimodal skew-normal')
        %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  mu1_2  sigma1_2  scale_2  mu2_2  sigma2_2  rho_2  mu1_3  sigma1_3  scale_3  mu2_3  sigma2_3  rho_3    mu1_4  sigma1_4  scale_4  mu2_4  sigma2_4    rho_4
        loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,    -0.1,      .01,  -10000,    -1,      .01,      -1];
        upBounds = [   1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,   	   1,        1,   10000,    10,       10,       1];
    elseif strcmpi(distType,'Skew-normal')
        %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  skew1_1  skew2_1
        loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,     -10,     -10, -10000,    -1,      .01,    -1,   .01,      .01,  -10000,    -1,      .01,    -1,     -10,     -10,      -10,   -10,      -10,     -10];
        upBounds = [   1,        1,   10000,    10,       10,     1,      10,      10,  10000,    10,       10,     1,     2,        1,   10000,    10,       10,     1,      10,      10,       10,    10,       10,      10];
    elseif ~strcmpi(distType,'Trimodal skew-normal')
        %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  mu1_2  sigma1_2  scale_2  mu2_2  sigma2_2  rho_2  skew1_1  skew2_1  skew1_2  skew2_2
        loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,     -10,     -10,     -10,     -10,   .01,    -1,     -10,     -10,      -10,   -10,      -10,     -10];
        upBounds = [   1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,      10,      10,      10,      10,    10,     1,      10,      10,       10,    10,       10,      10];
    else
        %          mu1_1  sigma1_1  scale_1  mu2_1  sigma2_1  rho_1  mu1_2  sigma1_2  scale_2  mu2_2  sigma2_2  rho_2  mu1_3  sigma1_3  scale_3  mu2_3  sigma2_3  rho_3  skew1_1  skew2_1  skew1_2  skew2_2  skew1_3  skew2_3
        loBounds = [-0.1,      .01,       0,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,  -0.1,      .01,  -10000,    -1,      .01,    -1,     -10,     -10,     -10,     -10,     -10,     -10];
        upBounds = [   1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,      10,      10,      10,      10,      10,      10];
    end
    handles.options.bounds = [loBounds; upBounds];
end
if ~isfield(handles.options, 'nGenerations')
    handles.options.nGenerations = 50;
end
if ~isfield(handles.options, 'parallelCores')
    pool = gcp('nocreate');
    if isempty(pool)
        handles.options.parallelCores = 0;
    else
        handles.options.parallelCores = pool.NumWorkers;
    end
end
if ~isfield(handles.options, 'fullParallel')
    handles.options.fullParallel = 0;
end
if ~isfield(handles.options, 'tauRange')
    handles.options.tauRange = 0;
end
if ~isfield(handles.options, 'figureDisplay')
    handles.options.figureDisplay = 1;
end
if ~isfield(handles.options, 'cluster')
    handles.options.cluster = 0;
end
if ~isfield(handles.options, 'clusterProfile')
    handles.options.clusterProfile = 'local';
end
if ~isfield(handles.options, 'cliff')
    handles.options.cliff = 0;
end
if ~isfield(handles.options, 'SDupbound')
    handles.options.SDupbound = 0;
end
if ~isfield(handles.options, 'SDlobound')
    handles.options.SDlobound = 0;
end

options = minisOptions(distributionType, handles.options);
if ~isempty(options)
    handles.options = options;
end

% Update handles structure
guidata(hObject, handles);










%% Info menu:
% --------------------------------------------------------------------
function info_Callback(~, ~, ~)


% --------------------------------------------------------------------
function aboutInfo_Callback(~, ~, ~)


% --------------------------------------------------------------------
function helpInfo_Callback(~, ~, ~)


% --------------------------------------------------------------------
function symbolsInfo_Callback(~, ~, ~)










%% Play button:
% --------------------------------------------------------------------
function startTool_ClickedCallback(hObject, ~, handles)
% handles    structure with handles and user data (see GUIDATA)

% Initialise global variables
global STOP
STOP = 0;

try
  
    % Initialise general variables:
    parallelCores = handles.options.parallelCores;
    if ischar(parallelCores)
        parallelCores = str2double(parallelCores);
    end

    graphicsFormats = {'*.fig','Matlab figure (*.fig)'
        '*.ai','Adobe Illustrator `88 (*.ai)'
        '*.bmp','Windows bitmap (*.bmp)'
        '*.emf','Enhanced metafile (*.emf)'
        '*.eps','EPS Level 1 (*.eps)'
        '*.jpg','JPEG image (*.jpg)'
        '*.m','MATLAB file (*.m)'
        '*.pbm','Portable bitmap (*.pbm)'
        '*.pcx','Paintbrush 24-bit (*.pcx)'
        '*.pdf','Portable Document Format (*.pdf)'
        '*.pgm','Portable Graymap (*.pgm)'
        '*.png','Portable Network Graphics (*.png)'
        '*.ppm','Portable Pixmap (*.ppm)'
        '*.tif','TIFF image, compressed (*.tif)'};

    [initialised, detectionParameters] = initDetectParam(handles);
    if ~initialised
        return
    end

    classificationParameters = initClassParam(handles, detectionParameters.Ampupbound);

    % Determine the choice of action:
    preprocess = get(handles.preprocessButton,'value');
    detection = get(handles.detectOnly,'value');
    detectBatch = get(handles.detectBatch,'value');
    detectCompare = get(handles.detectCompare,'value');
    errorBounds = get(handles.errorBounds,'value');
    autoDistributionFit = get(handles.autoDistFit, 'value');
    simulation = get(handles.simulate, 'value');


    if preprocess
        disp('Preprocessing task initiated');
        handles.state.String = '  State: Preprocessing task initiated';
        guidata(hObject, handles);

        if fafTest
            errmsgNoFile = 'Error: Preprocess mode is unavailable in the trial version of the software';
            msgbox(errmsgNoFile,'Error','Error');
            return
        end

        % Pulses:
        [initialised, startPulse, endPulse] = initPulseTarget(handles);
        if ~initialised
            return
        end

        % Pre-process data:
        preprocessMinis(startPulse, endPulse, detectionParameters, classificationParameters, handles.ld, graphicsFormats, parallelCores, handles);
        return

    elseif detection
        disp('Detection task initiated');
        handles.state.String = '  State: Detection task initiated';
        guidata(hObject, handles);

        % Target file:
        if fafTest
            load('fafTarget.mat'); %#ok<*LOAD>
        else
            targetFilename = get(handles.loadTargetFileInput, 'string');
            if ~preprocess && (strcmpi(targetFilename, '>>> <<<') || isempty(targetFilename))
                errmsgNoFile = 'Error: no target file supplied';
                msgbox(errmsgNoFile,'Error','Error');
                return
            end
        end

        % Pulses and glitches:
        [initialised, targetExcludedTimes] = initExclTimesTarget(handles);
        if ~initialised
            return
        end

        % Estimate tau_m based on impulses for later using it as a reference time interval for tau_m estimation based on spontaneous PSPs:
        if detectionParameters.voltageClamp
            waveform.estimate = 1;
            waveform.riseTimeArrayExt = classificationParameters.riseTimeArray;
            waveform.tau_m = 10;
        else
            waveform = initWaveform(handles, detectionParameters.pulseDuration, classificationParameters.riseTimeArray, targetExcludedTimes, targetFilename);
        end

        % Filtering:
        filtering = noiseFilterDlg;
        filtering.excludedTimes = targetExcludedTimes;

        % Detect minis-like events in a file of choice:
        [minisArray, ~, waveform, ~, F] = detectMinisStandalone(true, targetFilename, detectionParameters, targetExcludedTimes, classificationParameters,...
            parallelCores, 'none', true, waveform, filtering, [1 1 1 1]);
        if ~isempty(waveform)
            minisArray = [minisArray repmat([waveform.parameters.averageAmp waveform.parameters.medianAmp waveform.parameters.tau_m], size(minisArray,1), 1)];
        end

        % Save the detected event log and figures:
        handles.ld = saveTargetEventLog(minisArray, F, waveform, filtering, detectionParameters.RTinterval, graphicsFormats, handles.ld);

    elseif detectBatch
        disp('Batch detection task initiated');
        handles.state.String = '  State: Batch detection task initiated';
        guidata(hObject, handles);

        % Target file:
        if fafTest
            load('fafTarget.mat'); %#ok<*LOAD>
        else
            targetFilename = get(handles.loadTargetFileInput, 'string');
            if ~preprocess && (strcmpi(targetFilename, '>>> <<<') || isempty(targetFilename))
                errmsgNoFile = 'Error: no target file supplied';
                msgbox(errmsgNoFile,'Error','Error');
                return
            end
        end

        % Pulses and glitches:
        [initialised, targetExcludedTimes] = initExclTimesTarget(handles);
        if ~initialised
            return
        end

        % Estimate tau_m based on impulses for later using it as a reference time interval for tau_m estimation based on spontaneous PSPs:
        if detectionParameters.voltageClamp
            waveform.estimate = 1;
            waveform.riseTimeArrayExt = classificationParameters.riseTimeArray;
            waveform.tau_m = 10;
        else
            waveform = initWaveform(handles, detectionParameters.pulseDuration, classificationParameters.riseTimeArray, targetExcludedTimes, targetFilename);
            waveform.estimate = 2;
        end

        % Filtering:
        filtering = noiseFilterDlg;
        filtering.excludedTimes = targetExcludedTimes;
        options.Resize = 'on';
        options.WindowStyle = 'normal';
        if strcmpi(filtering.state, 'on')
            filtering.filtfs = inputdlg('Choose stop band frequencies separated by commas:','FFT',1,{'50, 150'},options);
        end
        if isfield(filtering, 'filtfs') && isempty(filtering.filtfs)
          filtering.state = 'off';
        end

        % Get the target folder structure:
        targetFolder = fileparts(targetFilename);
        folderStruct = dir(targetFolder);

        % Detect minis-like events in all ABF files inside the target folder:
        for iFile = 1:numel(folderStruct)
            targetFilename = fullfile(folderStruct(iFile).folder, folderStruct(iFile).name);
            if endsWith(targetFilename, '.abf')
                [minisArray, ~, waveform] = detectMinisStandalone(false, targetFilename, detectionParameters, targetExcludedTimes, classificationParameters,...
                    parallelCores, 'none', false, waveform, filtering, [1 1 1 1], false);
                if ~isempty(waveform)
                    minisArray = [minisArray repmat([waveform.parameters.averageAmp waveform.parameters.medianAmp waveform.parameters.tau_m], size(minisArray,1), 1)]; %#ok<AGROW> 
                end
        
                % Save the detected event log and figures:
                eventLogFile = [targetFilename(1:end-3) 'txt'];
                handles.ld = saveTargetEventLogReduced(minisArray, eventLogFile, detectionParameters.RTinterval);
            end
        end

    elseif detectCompare
        disp('Detect and compare task initiated');
        handles.state.String = '  State: Detect and compare task initiated';
        guidata(hObject, handles);

        % Target file:
        if fafTest
            load('fafTarget.mat');
        else
            targetFilename = get(handles.loadTargetFileInput, 'string');
            if ~preprocess && (strcmpi(targetFilename, '>>> <<<') || isempty(targetFilename))
                errmsgNoFile = 'Error: no target file supplied';
                msgbox(errmsgNoFile,'Error','Error');
                return
            end
        end

        % Noise file:
        if fafTest
            load('fafNoise.mat');
        else
            noiseFilename = get(handles.loadNoiseFileInput,'string');
            if ~errorBounds && (strcmpi(noiseFilename, '>>> <<<') || isempty(noiseFilename))
                errmsgNoFile = 'Error: no noise file supplied';
                msgbox(errmsgNoFile,'Error','Error');
                return
            end
        end

        % Pulses and glitches:
        [initialised, targetExcludedTimes] = initExclTimesTarget(handles);
        if ~initialised
            return
        end

        [initialised, noiseExcludedTimes] = initExclTimesNoise(handles);
        if ~initialised
            return
        end

        % Estimate tau_m based on impulses for later using it as a reference time interval for tau_m estimation based on spontaneous PSPs:
        waveform = initWaveform(handles, detectionParameters.pulseDuration, classificationParameters.riseTimeArray, targetExcludedTimes, targetFilename);

        % Filtering:
        filtering = noiseFilterDlg;

        % Detect minis-like events in both target and noise files and compare them:
        [minisTarget, minisNoise, waveform, F] = compareMinis(targetFilename, noiseFilename, targetExcludedTimes, noiseExcludedTimes, detectionParameters,...
            classificationParameters, filtering, waveform, parallelCores);

        % Save the detected event logs and figures:
        handles.ld = saveEventLog(minisTarget, minisNoise, F, waveform, filtering, detectionParameters.RTinterval, graphicsFormats, handles.ld);

    elseif errorBounds
        disp('Error bound estimation task initiated');
        handles.state.String = '  State: Error bound estimation task initiated';
        guidata(hObject, handles);

        if fafTest
            errmsgNoFile = 'Error: Preprocess mode is unavailable in the trial version of the software';
            msgbox(errmsgNoFile,'Error','Error');
            return
        end

        % Pulses:
        [initialised, startPulse, endPulse] = initPulseTarget(handles);
        if ~initialised
            return
        end
        targetExcludedTimes = struct('startPulse', startPulse, 'endPulse', endPulse);

        % Estimate error margin for optimisation:
        errorMinis(targetExcludedTimes, detectionParameters, classificationParameters, handles.ld, handles.wd, parallelCores, handles);
        return

    elseif autoDistributionFit
        disp('Distribution fitting task initiated');
        handles.state.String = '  State: Distribution fitting task initiated';
        guidata(hObject, handles);

        % Target file:
        if fafTest
            load('fafTarget.mat');
        else
            targetFilename = get(handles.loadTargetFileInput, 'string');
            if ~preprocess && (strcmpi(targetFilename, '>>> <<<') || isempty(targetFilename))
                errmsgNoFile = 'Error: no target file supplied';
                msgbox(errmsgNoFile,'Error','Error');
                return
            end
        end

        % Noise file:
        if fafTest
            load('fafNoise.mat');
        else
            noiseFilename = get(handles.loadNoiseFileInput,'string');
            if ~errorBounds && (strcmpi(noiseFilename, '>>> <<<') || isempty(noiseFilename))
                errmsgNoFile = 'Error: no noise file supplied';
                msgbox(errmsgNoFile,'Error','Error');
                return
            end
        end

        % Pulses and glitches:
        [initialised, targetExcludedTimes] = initExclTimesTarget(handles);
        if ~initialised
            return
        end

        [initialised, noiseExcludedTimes] = initExclTimesNoise(handles);
        if ~initialised
            return
        end

        % Optimisation parameters:
        [initialised, optimisationParameters] = initOptimParam(handles);
        if ~initialised
            return
        end

        % Simulation parameters:
        [initialised, simulationParameters] = initSimParam(handles, detectionParameters.pulseDuration, classificationParameters.riseTimeArray,...
            targetExcludedTimes, targetFilename);
        if ~initialised
            return
        end

        % Filtering:
        filtering = noiseFilterDlg;
        filtering.targetExcludedTimes = targetExcludedTimes;
        filtering.noiseExcludedTimes = noiseExcludedTimes;

        % Start automatic distribution fitting:
        optimiseMinis(handles.wd, targetFilename, noiseFilename, targetExcludedTimes, noiseExcludedTimes, detectionParameters, simulationParameters,...
            optimisationParameters, classificationParameters, filtering, 'TwoDs', parallelCores);
    elseif simulation
        disp('Simulation task initiated');
        handles.state.String = '  State: Simulation task initiated';
        guidata(hObject, handles);

        % Noise file:
        if fafTest
            load('fafNoise.mat');
        else
            noiseFilename = get(handles.loadNoiseFileInput,'string');
            if strcmpi(noiseFilename, '>>> <<<') || strcmpi(noiseFilename, '...') || strcmpi(noiseFilename, '') || isempty(noiseFilename)
                noiseFilename = [];
            end
        end

        % Pulses and glitches:
        [initialised, noiseExcludedTimes] = initExclTimesNoise(handles);
        if ~initialised
            return
        end

        % Optimisation parameters:
        [initialised, optimisationParameters] = initOptimParam(handles);
        if ~initialised
            return
        end

        % Simulation parameters:
        [initialised, simulationParameters] = initSimParamReduced(handles);
        if ~initialised
            return
        end

        % Filtering:
        if ~isempty(noiseFilename)
            filtering = noiseFilterDlg;
            filtering.noiseExcludedTimes = noiseExcludedTimes;
        else
            filtering.state = 'off';
            filtering.noiseExcludedTimes = [];
        end

        % Simulate minis:
        simulateDetectEvaluate(noiseFilename, noiseExcludedTimes, detectionParameters, simulationParameters, optimisationParameters,...
            classificationParameters, filtering, parallelCores);
    end

    fclose all;
    disp('Task completed');
    handles.state.String = '  State: Task completed';
catch me
    if strcmp(me.message, 'Interrupted by user')
        handles.state.String = '  State: Interrupted by user';
    else
        handles.state.String = '  State: Execution failed. See minisLogFile inside the working directory.';
        guidata(hObject, handles);
        disp('Program execution failed');
        rethrow(me);
    end
end

% Update handles structure
guidata(hObject, handles);










%% Stop button:
% --------------------------------------------------------------------
function stopTool_ClickedCallback(hObject, ~, handles)
% hObject    handle to stopTool (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

global STOP                                                                 % STOP issues a global stop message
button = questdlg('Do you want to stop the current task execution?','Stop Program','Yes','No','No');
if strcmpi(button, 'Yes')
    handles.state.String = '  State: Execution termination requested';
    guidata(hObject, handles);
    disp('Program execution termination requested');
    STOP = 1;
end

% Update handles structure
guidata(hObject, handles);

%% GUI delete function
% --- Executes during object deletion, before destroying properties.
function detectionPanel_DeleteFcn(~, ~, ~)
% hObject    handle to detectionPanel (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

if (~isdeployed)
    disp('Logging is off');
    diary off
end
disp('minis is terminating');
