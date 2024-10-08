function varargout = changeParametersGUI(varargin)
% CHANGEPARAMETERSGUI MATLAB code for changeParametersGUI.fig
%      CHANGEPARAMETERSGUI, by itself, creates a new CHANGEPARAMETERSGUI or raises the existing
%      singleton*.
%
%      H = CHANGEPARAMETERSGUI returns the handle to a new CHANGEPARAMETERSGUI or the handle to
%      the existing singleton*.
%
%      CHANGEPARAMETERSGUI('CALLBACK',hObject,~,handles,...) calls the local
%      function named CALLBACK in CHANGEPARAMETERSGUI.M with the given input arguments.
%
%      CHANGEPARAMETERSGUI('Property','Value',...) creates a new CHANGEPARAMETERSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before changeParametersGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to changeParametersGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help changeParametersGUI

% Last Modified by GUIDE v2.5 16-Mar-2012 15:56:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @changeParametersGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @changeParametersGUI_OutputFcn, ...
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


% --- Executes just before changeParametersGUI is made visible.
function changeParametersGUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to changeParametersGUI (see VARARGIN)

% Initialise variables:
handles.simulationParameters = varargin{1};
set(handles.tau_mEdit,'String',handles.simulationParameters.tau_m);
set(handles.tau_sy1Edit,'String',handles.simulationParameters.tau_sy1);
set(handles.tau_sy2Edit,'String',handles.simulationParameters.tau_sy2);
set(handles.nSeriesEdit,'String',handles.simulationParameters.nSeries);
set(handles.LEdit,'String',handles.simulationParameters.basalL);

% Simulation Panel:
backColour = javaObjectEDT(java.awt.Color(.79, .82, .87));
pposSimulationPanel = getpixelposition(handles.simulationPanel);

% tau_m text:
ppostau_mText = [.0925*pposSimulationPanel(3) pposSimulationPanel(2)+.59*pposSimulationPanel(4) .25*pposSimulationPanel(3) .18*pposSimulationPanel(4)];
labelStr = '<html><font style="font-size:1.3em;" face = "MS Sans Serif">&#964<sub>m</sub> (passive membrane time constant, ms):</font></html>';
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
setBackground(jLabel, backColour);
[handles.ppostau_mText, ppostau_mTextContainer] = javacomponent(jLabel, ppostau_mText, gcf); %#ok<*JAVCM>
set(ppostau_mTextContainer,'units','normalized');

% nSeries text:
pposnSeries = [.0925*pposSimulationPanel(3) pposSimulationPanel(2)+.36*pposSimulationPanel(4) .25*pposSimulationPanel(3) .18*pposSimulationPanel(4)];
labelStr = '<html><font style="font-size:1.3em;" face = "MS Sans Serif">Number of series terms:</font></html>';
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
setBackground(jLabel, backColour);
[handles.nSeries, nSeriesContainer] = javacomponent(jLabel, pposnSeries, gcf);
set(nSeriesContainer,'units','normalized');

% L text:
pposL = [.0925*pposSimulationPanel(3) pposSimulationPanel(2)+.12*pposSimulationPanel(4) .25*pposSimulationPanel(3) .18*pposSimulationPanel(4)];
labelStr = '<html><font style="font-size:1.3em;" face = "MS Sans Serif">L (electrotonic length, &#955):</font></html>';
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
setBackground(jLabel, backColour);
[handles.L, LContainer] = javacomponent(jLabel, pposL, gcf);
set(LContainer,'units','normalized');

% tau_sy1 text:
ppostau_sy1 = [.48*pposSimulationPanel(3) pposSimulationPanel(2)+.59*pposSimulationPanel(4) .31*pposSimulationPanel(3) .21*pposSimulationPanel(4)];
labelStr = '<html><font style="font-size:1.3em;" face = "MS Sans Serif">&#964<sub>sy1</sub> (synaptic decay time constant for current input, ms):</font></html>';
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
setBackground(jLabel, backColour);
[handles.tau_sy1, tau_sy1Container] = javacomponent(jLabel, ppostau_sy1, gcf);
set(tau_sy1Container,'units','normalized');

% tau_sy2 text:
ppostau_sy2 = [.48*pposSimulationPanel(3) pposSimulationPanel(2)+.33*pposSimulationPanel(4) .31*pposSimulationPanel(3) .21*pposSimulationPanel(4)];
labelStr = '<html><font style="font-size:1.3em;" face = "MS Sans Serif">&#964<sub>sy2</sub> (synaptic rise time constant for current input, ms):</font></html>';
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
setBackground(jLabel, backColour);
[handles.tau_sy2, tau_sy2Container] = javacomponent(jLabel, ppostau_sy2, gcf);
set(tau_sy2Container,'units','normalized');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes changeParametersGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = changeParametersGUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
    varargout{1} = handles.simulationParameters;
    delete(handles.figure1);
end



function tau_mEdit_Callback(hObject, ~, handles) %#ok<*DEFNU>
% hObject    handle to tau_mEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

tau_mEdit = str2double(get(hObject,'String'));
if isnan(tau_mEdit)
    errmsg = 'Error: at least one of the entries is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.tau_mEdit,'String',tau_mEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function tau_mEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to tau_mEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function tau_sy1Edit_Callback(hObject, ~, handles)
% hObject    handle to tau_sy1Edit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

tau_sy1Edit = str2double(get(hObject,'String'));
if isnan(tau_sy1Edit)
    errmsg = 'Error: at least one of the entries is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.tau_sy1Edit,'String',tau_sy1Edit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function tau_sy1Edit_CreateFcn(hObject, ~, ~)
% hObject    handle to tau_sy1Edit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tau_sy2Edit_Callback(hObject, ~, handles)
% hObject    handle to tau_sy2Edit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

tau_sy2Edit = str2double(get(hObject,'String'));
if isnan(tau_sy2Edit)
    errmsg = 'Error: at least one of the entries is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.tau_sy2Edit,'String',tau_sy2Edit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function tau_sy2Edit_CreateFcn(hObject, ~, ~)
% hObject    handle to tau_sy2Edit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nSeriesEdit_Callback(hObject, ~, handles)
% hObject    handle to nSeriesEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

nSeriesEdit = str2double(get(hObject,'String'));
if isnan(nSeriesEdit)
    errmsg = 'Error: at least one of the entries is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.nSeriesEdit,'String',nSeriesEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function nSeriesEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to nSeriesEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function LEdit_Callback(hObject, ~, handles)
% hObject    handle to LEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

LEdit = str2double(get(hObject,'String'));
if isnan(LEdit)
    errmsg = 'Error: at least one of the entries is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.LEdit,'String',LEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function LEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to LEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in okButton.
function okButton_Callback(hObject, ~, handles)
% hObject    handle to okButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% Assign simulation parameters:
loSimAmp = handles.simulationParameters.loSimAmp;
L = str2double(get(handles.LEdit,'string'));
Q = handles.simulationParameters.Q;
d = handles.simulationParameters.basald;
nSeries = str2double(get(handles.nSeriesEdit,'string'));
tau_sy2 = str2double(get(handles.tau_sy2Edit,'string'));
tau_sy1 = str2double(get(handles.tau_sy1Edit,'string'));
R_i = handles.simulationParameters.R_i;
tau_m = str2double(get(handles.tau_mEdit,'string'));
C_m = handles.simulationParameters.C_m;
if isnan(L) || isnan(nSeries) || isnan(tau_sy2) || isnan(tau_sy1) || isnan(tau_m)
    errmsg = 'Error: at least one of the entries is NaN';
    msgbox(errmsg,'Error','Error');
    return
elseif tau_m < 1
    errmsg = 'Error: at least one of the entries is outside the range of acceptable values';
    msgbox(errmsg,'Error','Error');
    return
end
handles.simulationParameters = struct('tau_m', tau_m, 'C_m', C_m, 'R_i', R_i, 'tau_sy1', tau_sy1, 'tau_sy2', tau_sy2, 'nSeries', nSeries,...
    'basald', d, 'apicald', [], 'Q', Q, 'basalL', L, 'apicalL', [], 'sampleSize', 10, 'loSimAmp', loSimAmp);

% Update handles structure
guidata(hObject, handles);

uiresume(handles.figure1);
