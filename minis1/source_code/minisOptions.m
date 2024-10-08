function varargout = minisOptions(varargin)
%MINISOPTIONS M-file for minisOptions.fig
%      MINISOPTIONS, by itself, creates a new MINISOPTIONS or raises the existing
%      singleton*.
%
%      H = MINISOPTIONS returns the handle to a new MINISOPTIONS or the handle to
%      the existing singleton*.
%
%      MINISOPTIONS('Property','Value',...) creates a new MINISOPTIONS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to minisOptions_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MINISOPTIONS('CALLBACK') and MINISOPTIONS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MINISOPTIONS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help minisOptions

% Last Modified by GUIDE v2.5 07-Nov-2012 14:55:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @minisOptions_OpeningFcn, ...
    'gui_OutputFcn',  @minisOptions_OutputFcn, ...
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


% --- Executes just before minisOptions is made visible.
function minisOptions_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to minisOptions (see VARARGIN)

% Initialise variables:
distributionType = varargin{1};
distributionType = strtrim(distributionType);
options = varargin{2};

set(handles.nGenerationsEdit,'String',options.nGenerations);
set(handles.parallelisationEdit,'String',options.parallelCores);
set(handles.figureCheckbox,'Value',options.figureDisplay);
set(handles.parallelCheckbox,'Value',options.fullParallel);
set(handles.tauRangeCheckbox,'Value',options.tauRange);
set(handles.clusterCheckbox,'Value',options.cluster);
set(handles.cliffCheckbox,'Value',options.cliff);
set(handles.SDupboundCheckbox,'Value',options.SDupbound);
set(handles.SDloboundCheckbox,'Value',options.SDlobound);
set(handles.clusterEdit,'String',options.clusterProfile);

set(handles.loboundEdit1,'String',options.bounds(1,1));
set(handles.loboundEdit2,'String',options.bounds(1,2));
set(handles.loboundEdit3,'String',options.bounds(1,3));
set(handles.loboundEdit4,'String',options.bounds(1,4));
set(handles.loboundEdit5,'String',options.bounds(1,5));
set(handles.loboundEdit6,'String',options.bounds(1,6));
set(handles.loboundEdit7,'String',options.bounds(1,7));
set(handles.loboundEdit8,'String',options.bounds(1,8));
set(handles.loboundEdit9,'String',options.bounds(1,9));
set(handles.loboundEdit10,'String',options.bounds(1,10));
set(handles.loboundEdit11,'String',options.bounds(1,11));
set(handles.loboundEdit12,'String',options.bounds(1,12));
set(handles.loboundEdit13,'String',options.bounds(1,13));
set(handles.loboundEdit14,'String',options.bounds(1,14));
set(handles.loboundEdit15,'String',options.bounds(1,15));
set(handles.loboundEdit16,'String',options.bounds(1,16));
set(handles.loboundEdit17,'String',options.bounds(1,17));
set(handles.loboundEdit18,'String',options.bounds(1,18));
set(handles.loboundEdit19,'String',options.bounds(1,19));
set(handles.loboundEdit20,'String',options.bounds(1,20));
set(handles.loboundEdit21,'String',options.bounds(1,21));
set(handles.loboundEdit22,'String',options.bounds(1,22));
set(handles.loboundEdit23,'String',options.bounds(1,23));
set(handles.loboundEdit24,'String',options.bounds(1,24));

set(handles.upboundEdit1,'String',options.bounds(2,1));
set(handles.upboundEdit2,'String',options.bounds(2,2));
set(handles.upboundEdit3,'String',options.bounds(2,3));
set(handles.upboundEdit4,'String',options.bounds(2,4));
set(handles.upboundEdit5,'String',options.bounds(2,5));
set(handles.upboundEdit6,'String',options.bounds(2,6));
set(handles.upboundEdit7,'String',options.bounds(2,7));
set(handles.upboundEdit8,'String',options.bounds(2,8));
set(handles.upboundEdit9,'String',options.bounds(2,9));
set(handles.upboundEdit10,'String',options.bounds(2,10));
set(handles.upboundEdit11,'String',options.bounds(2,11));
set(handles.upboundEdit12,'String',options.bounds(2,12));
set(handles.upboundEdit13,'String',options.bounds(2,13));
set(handles.upboundEdit14,'String',options.bounds(2,14));
set(handles.upboundEdit15,'String',options.bounds(2,15));
set(handles.upboundEdit16,'String',options.bounds(2,16));
set(handles.upboundEdit17,'String',options.bounds(2,17));
set(handles.upboundEdit18,'String',options.bounds(2,18));
set(handles.upboundEdit19,'String',options.bounds(2,19));
set(handles.upboundEdit20,'String',options.bounds(2,20));
set(handles.upboundEdit21,'String',options.bounds(2,21));
set(handles.upboundEdit22,'String',options.bounds(2,22));
set(handles.upboundEdit23,'String',options.bounds(2,23));
set(handles.upboundEdit24,'String',options.bounds(2,24));

backColour = javaObjectEDT(java.awt.Color(.79, .82, .87));
pposParamBoundPanel = getpixelposition(handles.parameterBoundsPanel);

% Define position variables:
pposParamNames = zeros(24,4);
for position = 1:size(pposParamNames,1)
    pposParamNames(position,:) = [(.027+.04058*(position-1))*pposParamBoundPanel(3) pposParamBoundPanel(2)+.57*pposParamBoundPanel(4) .09*pposParamBoundPanel(3) .18*pposParamBoundPanel(4)];
end

% Define script variables:
mu1_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#956;1<sub>1</sub></font></html>';
sigma1_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#963;1<sub>1</sub></font></html>';
scale_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">scale1<sub>1</sub></font></html>';
mu1_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#956;1<sub>2</sub></font></html>';
sigma1_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#963;1<sub>2</sub></font></html>';
scale_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">scale1<sub>2</sub></font></html>';
mu1_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#956;1<sub>3</sub></font></html>';
sigma1_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#963;1<sub>3</sub></font></html>';
scale_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">scale1<sub>3</sub></font></html>';
mu1_4Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#956;1<sub>4</sub></font></html>';
sigma1_4Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#963;1<sub>4</sub></font></html>';
scale_4Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">scale1<sub>4</sub></font></html>';

mu2_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#956;2<sub>1</sub></font></html>';
sigma2_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#963;2<sub>1</sub></font></html>';
mu2_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#956;2<sub>2</sub></font></html>';
sigma2_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#963;2<sub>2</sub></font></html>';
mu2_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#956;2<sub>3</sub></font></html>';
sigma2_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#963;2<sub>3</sub></font></html>';
mu2_4Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#956;2<sub>4</sub></font></html>';
sigma2_4Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#963;2<sub>4</sub></font></html>';

rho_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#961;<sub>1</sub></font></html>';
rho_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#961;<sub>2</sub></font></html>';
rho_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#961;<sub>3</sub></font></html>';
rho_4Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#961;<sub>4</sub></font></html>';

A_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">A<sub>1</sub></font></html>';
A_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">A<sub>2</sub></font></html>';
A_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">A<sub>3</sub></font></html>';
alpha_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#945;<sub>1</sub></font></html>';
alpha_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#945;<sub>2</sub></font></html>';
alpha_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#945;<sub>3</sub></font></html>';
alpha1_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#945;1<sub>1</sub></font></html>';
alpha2_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#945;2<sub>1</sub></font></html>';
alpha1_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#945;1<sub>2</sub></font></html>';
alpha2_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#945;2<sub>2</sub></font></html>';
alpha1_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#945;1<sub>3</sub></font></html>';
alpha2_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#945;2<sub>3</sub></font></html>';
B_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">B<sub>1</sub></font></html>';
B_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">B<sub>2</sub></font></html>';
B_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">B<sub>3</sub></font></html>';
beta_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#946;<sub>1</sub></font></html>';
beta_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#946;<sub>2</sub></font></html>';
beta_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#946;<sub>3</sub></font></html>';

theta_1Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#952;<sub>1</sub></font></html>';
theta_2Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#952;<sub>2</sub></font></html>';
theta_3Str = '<html><font style="font-size:1.1em;" face = "MS Sans Serif">&#952;<sub>3</sub></font></html>';

% Draw:
if ~strcmpi(distributionType,'Log Normal') && ~strcmpi(distributionType,'Bimodal Log Normal') && ~strcmpi(distributionType,'Trimodal Log Normal')
    jLabel = javaObjectEDT('javax.swing.JLabel', mu1_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText1 = javacomponent(jLabel, pposParamNames(1,:), gcf); %#ok<*JAVCM>
    
    jLabel = javaObjectEDT('javax.swing.JLabel', sigma1_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText2 = javacomponent(jLabel, pposParamNames(2,:), gcf);
    
    jLabel = javaObjectEDT('javax.swing.JLabel', scale_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText3 = javacomponent(jLabel, pposParamNames(3,:), gcf);
    
    jLabel = javaObjectEDT('javax.swing.JLabel', mu2_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText4 = javacomponent(jLabel, pposParamNames(4,:), gcf);
    
    jLabel = javaObjectEDT('javax.swing.JLabel', sigma2_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText5 = javacomponent(jLabel, pposParamNames(5,:), gcf);
    
    if ~strcmpi(distributionType,'Gaussian') && ~strcmpi(distributionType,'Bimodal Gaussian') && ~strcmpi(distributionType,'Trimodal Gaussian')...
            && ~strcmpi(distributionType,'Quadrimodal Gaussian')
        jLabel = javaObjectEDT('javax.swing.JLabel', rho_1Str);
        setBackground(jLabel, backColour);
        handles.ParamNamesText6 = javacomponent(jLabel, pposParamNames(6,:), gcf);
        if ~strcmpi(distributionType,'Normal') && ~strcmpi(distributionType,'Skew-normal')
            jLabel = javaObjectEDT('javax.swing.JLabel', mu1_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText7 = javacomponent(jLabel, pposParamNames(7,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', sigma1_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText8 = javacomponent(jLabel, pposParamNames(8,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', scale_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText9 = javacomponent(jLabel, pposParamNames(9,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', mu2_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText10 = javacomponent(jLabel, pposParamNames(10,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', sigma2_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText11 = javacomponent(jLabel, pposParamNames(11,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', rho_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText12 = javacomponent(jLabel, pposParamNames(12,:), gcf);
            
            if ~strcmpi(distributionType,'Bimodal Normal') && ~strcmpi(distributionType,'Bimodal Skew-normal')
                jLabel = javaObjectEDT('javax.swing.JLabel', mu1_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText13 = javacomponent(jLabel, pposParamNames(13,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', sigma1_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText14 = javacomponent(jLabel, pposParamNames(14,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', scale_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText15 = javacomponent(jLabel, pposParamNames(15,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', mu2_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText16 = javacomponent(jLabel, pposParamNames(16,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', sigma2_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText17 = javacomponent(jLabel, pposParamNames(17,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', rho_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText18 = javacomponent(jLabel, pposParamNames(18,:), gcf);
                if strcmpi(distributionType,'Trimodal Skew-normal')
                    jLabel = javaObjectEDT('javax.swing.JLabel', alpha1_1Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText19 = javacomponent(jLabel, pposParamNames(19,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', alpha2_1Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText20 = javacomponent(jLabel, pposParamNames(20,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', alpha1_2Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText21 = javacomponent(jLabel, pposParamNames(21,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', alpha2_2Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText22 = javacomponent(jLabel, pposParamNames(22,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', alpha1_3Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText23 = javacomponent(jLabel, pposParamNames(23,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', alpha2_3Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText24 = javacomponent(jLabel, pposParamNames(24,:), gcf);
                elseif ~strcmpi(distributionType,'Trimodal Normal')
                    jLabel = javaObjectEDT('javax.swing.JLabel', mu1_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText19 = javacomponent(jLabel, pposParamNames(19,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', sigma1_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText20 = javacomponent(jLabel, pposParamNames(20,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', scale_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText21 = javacomponent(jLabel, pposParamNames(21,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', mu2_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText22 = javacomponent(jLabel, pposParamNames(22,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', sigma2_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText23 = javacomponent(jLabel, pposParamNames(23,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', rho_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText24 = javacomponent(jLabel, pposParamNames(24,:), gcf);
                end
            elseif strcmpi(distributionType,'Bimodal Skew-normal')
                jLabel = javaObjectEDT('javax.swing.JLabel', alpha1_1Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText13 = javacomponent(jLabel, pposParamNames(13,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', alpha2_1Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText14 = javacomponent(jLabel, pposParamNames(14,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', alpha1_2Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText15 = javacomponent(jLabel, pposParamNames(15,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', alpha2_2Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText16 = javacomponent(jLabel, pposParamNames(16,:), gcf);
            end
        elseif strcmpi(distributionType,'Skew-normal')
            jLabel = javaObjectEDT('javax.swing.JLabel', alpha1_1Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText7 = javacomponent(jLabel, pposParamNames(7,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', alpha2_1Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText8 = javacomponent(jLabel, pposParamNames(8,:), gcf);
        end
    else
        jLabel = javaObjectEDT('javax.swing.JLabel', theta_1Str);
        setBackground(jLabel, backColour);
        handles.ParamNamesText6 = javacomponent(jLabel, pposParamNames(6,:), gcf);
        if strcmpi(distributionType,'Bimodal Gaussian') || strcmpi(distributionType,'Trimodal Gaussian') || strcmpi(distributionType,'Quadrimodal Gaussian')
            jLabel = javaObjectEDT('javax.swing.JLabel', mu1_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText7 = javacomponent(jLabel, pposParamNames(7,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', sigma1_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText8 = javacomponent(jLabel, pposParamNames(8,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', scale_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText9 = javacomponent(jLabel, pposParamNames(9,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', mu2_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText10 = javacomponent(jLabel, pposParamNames(10,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', sigma2_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText11 = javacomponent(jLabel, pposParamNames(11,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', theta_2Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText12 = javacomponent(jLabel, pposParamNames(12,:), gcf);
            if strcmpi(distributionType,'Trimodal Gaussian') || strcmpi(distributionType,'Quadrimodal Gaussian')
                jLabel = javaObjectEDT('javax.swing.JLabel', mu1_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText7 = javacomponent(jLabel, pposParamNames(13,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', sigma1_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText8 = javacomponent(jLabel, pposParamNames(14,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', scale_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText9 = javacomponent(jLabel, pposParamNames(15,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', mu2_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText10 = javacomponent(jLabel, pposParamNames(16,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', sigma2_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText11 = javacomponent(jLabel, pposParamNames(17,:), gcf);
                
                jLabel = javaObjectEDT('javax.swing.JLabel', theta_3Str);
                setBackground(jLabel, backColour);
                handles.ParamNamesText12 = javacomponent(jLabel, pposParamNames(18,:), gcf);
                if strcmpi(distributionType,'Quadrimodal Gaussian')
                    jLabel = javaObjectEDT('javax.swing.JLabel', mu1_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText7 = javacomponent(jLabel, pposParamNames(19,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', sigma1_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText8 = javacomponent(jLabel, pposParamNames(20,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', scale_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText9 = javacomponent(jLabel, pposParamNames(21,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', mu2_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText10 = javacomponent(jLabel, pposParamNames(22,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', sigma2_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText11 = javacomponent(jLabel, pposParamNames(23,:), gcf);
                    
                    jLabel = javaObjectEDT('javax.swing.JLabel', theta_4Str);
                    setBackground(jLabel, backColour);
                    handles.ParamNamesText12 = javacomponent(jLabel, pposParamNames(24,:), gcf);
                end
            end
        end
    end
else
    jLabel = javaObjectEDT('javax.swing.JLabel', mu1_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText1 = javacomponent(jLabel, pposParamNames(1,:), gcf);
    
    jLabel = javaObjectEDT('javax.swing.JLabel', A_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText1 = javacomponent(jLabel, pposParamNames(2,:), gcf);
    
    jLabel = javaObjectEDT('javax.swing.JLabel', alpha_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText2 = javacomponent(jLabel, pposParamNames(3,:), gcf);
    
    jLabel = javaObjectEDT('javax.swing.JLabel', scale_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText3 = javacomponent(jLabel, pposParamNames(4,:), gcf);
    
    jLabel = javaObjectEDT('javax.swing.JLabel', mu2_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText1 = javacomponent(jLabel, pposParamNames(5,:), gcf);
    
    jLabel = javaObjectEDT('javax.swing.JLabel', B_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText4 = javacomponent(jLabel, pposParamNames(6,:), gcf);
    
    jLabel = javaObjectEDT('javax.swing.JLabel', beta_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText5 = javacomponent(jLabel, pposParamNames(7,:), gcf);
    
    jLabel = javaObjectEDT('javax.swing.JLabel', rho_1Str);
    setBackground(jLabel, backColour);
    handles.ParamNamesText6 = javacomponent(jLabel, pposParamNames(8,:), gcf);
    if strcmpi(distributionType,'Bimodal Log Normal') || strcmpi(distributionType,'Trimodal Log Normal')
        jLabel = javaObjectEDT('javax.swing.JLabel', mu1_2Str);
        setBackground(jLabel, backColour);
        handles.ParamNamesText1 = javacomponent(jLabel, pposParamNames(9,:), gcf);
        
        jLabel = javaObjectEDT('javax.swing.JLabel', A_2Str);
        setBackground(jLabel, backColour);
        handles.ParamNamesText7 = javacomponent(jLabel, pposParamNames(10,:), gcf);
        
        jLabel = javaObjectEDT('javax.swing.JLabel', alpha_2Str);
        setBackground(jLabel, backColour);
        handles.ParamNamesText8 = javacomponent(jLabel, pposParamNames(11,:), gcf);
        
        jLabel = javaObjectEDT('javax.swing.JLabel', scale_2Str);
        setBackground(jLabel, backColour);
        handles.ParamNamesText9 = javacomponent(jLabel, pposParamNames(12,:), gcf);
        
        jLabel = javaObjectEDT('javax.swing.JLabel', mu2_2Str);
        setBackground(jLabel, backColour);
        handles.ParamNamesText1 = javacomponent(jLabel, pposParamNames(13,:), gcf);
        
        jLabel = javaObjectEDT('javax.swing.JLabel', B_2Str);
        setBackground(jLabel, backColour);
        handles.ParamNamesText10 = javacomponent(jLabel, pposParamNames(14,:), gcf);
        
        jLabel = javaObjectEDT('javax.swing.JLabel', beta_2Str);
        setBackground(jLabel, backColour);
        handles.ParamNamesText11 = javacomponent(jLabel, pposParamNames(15,:), gcf);
        
        jLabel = javaObjectEDT('javax.swing.JLabel', rho_2Str);
        setBackground(jLabel, backColour);
        handles.ParamNamesText12 = javacomponent(jLabel, pposParamNames(16,:), gcf);
        if strcmpi(distributionType,'Trimodal Log Normal')
            jLabel = javaObjectEDT('javax.swing.JLabel', mu1_3Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText1 = javacomponent(jLabel, pposParamNames(17,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', A_3Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText7 = javacomponent(jLabel, pposParamNames(18,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', alpha_3Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText8 = javacomponent(jLabel, pposParamNames(19,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', scale_3Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText9 = javacomponent(jLabel, pposParamNames(20,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', mu2_3Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText1 = javacomponent(jLabel, pposParamNames(21,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', B_3Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText10 = javacomponent(jLabel, pposParamNames(22,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', beta_3Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText11 = javacomponent(jLabel, pposParamNames(23,:), gcf);
            
            jLabel = javaObjectEDT('javax.swing.JLabel', rho_3Str);
            setBackground(jLabel, backColour);
            handles.ParamNamesText12 = javacomponent(jLabel, pposParamNames(24,:), gcf);
        end
    end
end

% Choose default command line output for minisOptions
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes minisOptions wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = minisOptions_OutputFcn(hObject, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
    varargout{1} = handles.output;
    delete(hObject);
else
    varargout{1} = [];
end



function upboundEdit1_Callback(hObject, ~, handles) %#ok<*DEFNU>
% hObject    handle to upboundEdit1 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit1 = str2double(get(hObject,'String'));
if isnan(upboundEdit1)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit1,'String',upboundEdit1);
    
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function upboundEdit1_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit1 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit2_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit2 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit2 = str2double(get(hObject,'String'));
if isnan(upboundEdit2)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit2,'String',upboundEdit2);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit2_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit2 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit3_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit3 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit3 = str2double(get(hObject,'String'));
if isnan(upboundEdit3)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit3,'String',upboundEdit3);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit3_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit3 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit4_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit4 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit4 = str2double(get(hObject,'String'));
if isnan(upboundEdit4)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit4,'String',upboundEdit4);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit4_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit4 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit5_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit5 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit5 = str2double(get(hObject,'String'));
if isnan(upboundEdit5)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit5,'String',upboundEdit5);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit5_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit5 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit6_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit6 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit6 = str2double(get(hObject,'String'));
if isnan(upboundEdit6)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit6,'String',upboundEdit6);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit6_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit6 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit1_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit1 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit1 = str2double(get(hObject,'String'));
if isnan(loboundEdit1)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit1,'String',loboundEdit1);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit1_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit1 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit2_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit2 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit2 = str2double(get(hObject,'String'));
if isnan(loboundEdit2)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit2,'String',loboundEdit2);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit2_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit2 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit3_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit3 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit3 = str2double(get(hObject,'String'));
if isnan(loboundEdit3)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit3,'String',loboundEdit3);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit3_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit3 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit4_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit4 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit4 = str2double(get(hObject,'String'));
if isnan(loboundEdit4)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit4,'String',loboundEdit4);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit4_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit4 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit5_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit5 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit5 = str2double(get(hObject,'String'));
if isnan(loboundEdit5)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit5,'String',loboundEdit5);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit5_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit5 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit6_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit6 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit6 = str2double(get(hObject,'String'));
if isnan(loboundEdit6)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit6,'String',loboundEdit6);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit6_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit6 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit7_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit7 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit7 = str2double(get(hObject,'String'));
if isnan(upboundEdit7)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit7,'String',upboundEdit7);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit7_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit7 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit8_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit8 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit8 = str2double(get(hObject,'String'));
if isnan(upboundEdit8)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit8,'String',upboundEdit8);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit8_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit8 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit9_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit9 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit9 = str2double(get(hObject,'String'));
if isnan(upboundEdit9)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit9,'String',upboundEdit9);
    
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function upboundEdit9_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit9 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit10_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit10 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit10 = str2double(get(hObject,'String'));
if isnan(upboundEdit10)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit10,'String',upboundEdit10);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit10_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit10 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit11_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit11 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit11 = str2double(get(hObject,'String'));
if isnan(upboundEdit11)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit11,'String',upboundEdit11);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit11_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit11 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit12_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit12 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit12 = str2double(get(hObject,'String'));
if isnan(upboundEdit12)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit12,'String',upboundEdit12);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit12_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit12 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit13_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit13 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit13 = str2double(get(hObject,'String'));
if isnan(upboundEdit13)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit13,'String',upboundEdit13);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit13_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit13 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit14_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit14 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit14 = str2double(get(hObject,'String'));
if isnan(upboundEdit14)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit14,'String',upboundEdit14);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit14_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit14 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit15_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit15 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit15 = str2double(get(hObject,'String'));
if isnan(upboundEdit15)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit15,'String',upboundEdit15);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit15_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit15 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit16_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit16 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit16 = str2double(get(hObject,'String'));
if isnan(upboundEdit16)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit16,'String',upboundEdit16);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit16_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit16 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit17_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit17 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit17 = str2double(get(hObject,'String'));
if isnan(upboundEdit17)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit17,'String',upboundEdit17);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit17_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit17 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit18_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit18 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit18 = str2double(get(hObject,'String'));
if isnan(upboundEdit18)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit18,'String',upboundEdit18);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit18_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit18 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit7_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit7 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit7 = str2double(get(hObject,'String'));
if isnan(loboundEdit7)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit7,'String',loboundEdit7);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit7_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit7 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit8_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit8 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit8 = str2double(get(hObject,'String'));
if isnan(loboundEdit8)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit8,'String',loboundEdit8);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit8_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit8 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit9_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit9 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit9 = str2double(get(hObject,'String'));
if isnan(loboundEdit9)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit9,'String',loboundEdit9);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit9_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit9 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit10_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit10 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit10 = str2double(get(hObject,'String'));
if isnan(loboundEdit10)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit10,'String',loboundEdit10);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit10_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit10 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit11_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit11 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit11 = str2double(get(hObject,'String'));
if isnan(loboundEdit11)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit11,'String',loboundEdit11);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit11_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit11 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit12_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit12 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit12 = str2double(get(hObject,'String'));
if isnan(loboundEdit12)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit12,'String',loboundEdit12);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit12_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit12 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit13_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit13 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit13 = str2double(get(hObject,'String'));
if isnan(loboundEdit13)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit13,'String',loboundEdit13);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit13_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit13 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit14_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit14 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit14 = str2double(get(hObject,'String'));
if isnan(loboundEdit14)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit14,'String',loboundEdit14);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit14_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit14 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit15_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit15 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit15 = str2double(get(hObject,'String'));
if isnan(loboundEdit15)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit15,'String',loboundEdit15);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit15_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit15 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit16_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit16 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit16 = str2double(get(hObject,'String'));
if isnan(loboundEdit16)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit16,'String',loboundEdit16);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit16_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit16 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit17_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit17 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit17 = str2double(get(hObject,'String'));
if isnan(loboundEdit17)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit17,'String',loboundEdit17);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit17_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit17 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit18_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit18 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit18 = str2double(get(hObject,'String'));
if isnan(loboundEdit18)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit18,'String',loboundEdit18);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit18_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit18 (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nGenerationsEdit_Callback(hObject, ~, handles)
% hObject    handle to nGenerationsEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

nGenerationsEdit = str2double(get(hObject,'String'));
if isnan(nGenerationsEdit)
    errmsg = 'Error: Number of generations is NaN';
    msgbox(errmsg,'Error','Error');
else
    set(handles.nGenerationsEdit,'String',nGenerationsEdit);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function nGenerationsEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to nGenerationsEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function parallelisationEdit_Callback(hObject, ~, handles)
% hObject    handle to parallelisationEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

parallelisation = get(hObject,'String');
if ~strcmpi(parallelisation, 'min') && ~strcmpi(parallelisation, 'max')
    parallelisation = str2double(parallelisation);
    if isnan(parallelisation)
        errmsg = 'Error: Number of parallel cores is NaN';
        msgbox(errmsg,'Error','Error');
        return
    end
    if parallelisation > 1
        scheduler = parcluster(parallel.defaultClusterProfile);
        if parallelisation > scheduler.NumWorkers
            errmsg = 'Error: not enough processor cores to support the required level of parallelisation. Please check your system settings.';
            msgbox(errmsg,'Error','Error');
            return
        end
    elseif parallelisation < 1
        errmsg = 'Error: Number of parallel cores is smaller than 1';
        msgbox(errmsg,'Error','Error');
        return
    end
end

set(handles.parallelisationEdit,'String',parallelisation);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function parallelisationEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to parallelisationEdit (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in acceptChangesButton.
function acceptChangesButton_Callback(hObject, ~, handles)
% hObject    handle to acceptChangesButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

data = guidata(handles.figure1);
parallelisation = get(data.parallelisationEdit,'String');
if ~strcmpi(parallelisation, 'min') && ~strcmpi(parallelisation, 'max')
    parallelisation = str2double(parallelisation);
    if isnan(parallelisation)
        errmsg = 'Error: Number of parallel cores is NaN';
        msgbox(errmsg,'Error','Error');
        return
    elseif parallelisation < 1
        errmsg = 'Error: Number of parallel cores is smaller than 1';
        msgbox(errmsg,'Error','Error');
        return
    end
end

if parallelisation > 1
    parallelPool = gcp('nocreate');
    if isempty(parallelPool)
      openWorkers = 1;
    else
      openWorkers = parallelPool.NumWorkers;
    end
    scheduler = parcluster(parallel.defaultClusterProfile);
    switch parallelisation
        case {'min', 1}
            parallelisation = 1;
            if parallelisation ~= openWorkers
                delete(parallelPool)
            end
        case 'max'
            parallelisation = scheduler.NumWorkers;
            if parallelisation ~= openWorkers
                delete(parallelPool)
                parpool;
            end
        otherwise
            parallelisation = floor(parallelisation);
            if parallelisation ~= openWorkers
                if parallelisation < 1
                    parallelisation = 1;
                    delete(parallelPool)
                elseif parallelisation > scheduler.NumWorkers
                    parallelisation = scheduler.NumWorkers;
                    set(data.parallelisationEdit,'String',parallelisation);
                    delete(parallelPool)
                    parpool;
                else
                    delete(parallelPool)
                    parpool(parallelisation);
                end
            end
    end
end
parallelisation = num2str(parallelisation);

lobounds = char(get(data.loboundEdit1,'string'), get(data.loboundEdit2,'string'), get(data.loboundEdit3,'string'),...
    get(data.loboundEdit4,'string'), get(data.loboundEdit5,'string'), get(data.loboundEdit6,'string'),get(data.loboundEdit7,'string'),...
    get(data.loboundEdit8,'string'), get(data.loboundEdit9,'string'), get(data.loboundEdit10,'string'), get(data.loboundEdit11,'string'),...
    get(data.loboundEdit12,'string'), get(data.loboundEdit13,'string'), get(data.loboundEdit14,'string'), get(data.loboundEdit15,'string'),...
    get(data.loboundEdit16,'string'), get(data.loboundEdit17,'string'), get(data.loboundEdit18,'string'), get(data.loboundEdit19,'string'),...
    get(data.loboundEdit20,'string'), get(data.loboundEdit21,'string'), get(data.loboundEdit22,'string'), get(data.loboundEdit23,'string'),...
    get(data.loboundEdit24,'string'));
lobounds = strarray2numarray(lobounds);
upbounds = char(get(data.upboundEdit1,'string'), get(data.upboundEdit2,'string'), get(data.upboundEdit3,'string'),...
    get(data.upboundEdit4,'string'), get(data.upboundEdit5,'string'), get(data.upboundEdit6,'string'), get(data.upboundEdit7,'string'),...
    get(data.upboundEdit8,'string'), get(data.upboundEdit9,'string'), get(data.upboundEdit10,'string'), get(data.upboundEdit11,'string'),...
    get(data.upboundEdit12,'string'), get(data.upboundEdit13,'string'), get(data.upboundEdit14,'string'), get(data.upboundEdit15,'string'),...
    get(data.upboundEdit16,'string'), get(data.upboundEdit17,'string'), get(data.upboundEdit18,'string'), get(data.upboundEdit19,'string'),...
    get(data.upboundEdit20,'string'), get(data.upboundEdit21,'string'), get(data.upboundEdit22,'string'), get(data.upboundEdit23,'string'),...
    get(data.upboundEdit24,'string'));
upbounds = strarray2numarray(upbounds);
bounds = [lobounds; upbounds];
nGenerations = str2double(get(data.nGenerationsEdit,'string'));
fullParallel = get(data.parallelCheckbox,'Value');
if strcmp(parallelisation, '1')
    delete(gcp('nocreate'));
    fullParallel = 0;
end
tauRange = get(data.tauRangeCheckbox,'Value');
cluster = get(data.clusterCheckbox,'Value');
figureDisplay = get(data.figureCheckbox,'Value');
cliff = get(data.cliffCheckbox,'Value');
SDupbound = get(data.SDupboundCheckbox,'Value');
SDlobound = get(data.SDloboundCheckbox,'Value');
clustProfile = get(data.clusterEdit,'string');

if sum(sum(isnan(lobounds))) || isnan(nGenerations)
    errmsg = 'Error: at least one of the lower bound entries is NaN';
    msgbox(errmsg,'Error','Error');
    return
end
if sum(sum(isnan(upbounds)))
    errmsg = 'Error: at least one of the upper bound entries is NaN';
    msgbox(errmsg,'Error','Error');
    return
end
if isnan(nGenerations)
    errmsg = 'Error: Number of generations is NaN';
    msgbox(errmsg,'Error','Error');
    return
end

handles.output = struct('bounds', bounds, 'nGenerations', nGenerations, 'parallelCores', parallelisation,...
    'fullParallel', fullParallel, 'tauRange', tauRange, 'cluster', cluster, 'cliff', cliff', 'figureDisplay', figureDisplay,...
    'clusterProfile', clustProfile, 'SDlobound', SDlobound, 'SDupbound', SDupbound);

guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes on button press in cliffCheckbox.
function cliffCheckbox_Callback(hObject, ~, handles)
% hObject    handle to cliffCheckbox (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in figureCheckbox.
function figureCheckbox_Callback(hObject, ~, handles)
% hObject    handle to figureCheckbox (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);


% --- Executes on button press in parallelCheckbox.
function parallelCheckbox_Callback(hObject, ~, handles)
% hObject    handle to parallelCheckbox (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);


% --- Executes on button press in clusterCheckbox.
function clusterCheckbox_Callback(hObject, ~, handles)
% hObject    handle to clusterCheckbox (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);



% --- Executes on button press in SDupboundCheckbox.
function SDupboundCheckbox_Callback(hObject, ~, handles)
% hObject    handle to SDupboundCheckbox (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);


% --- Executes on button press in SDloboundCheckbox.
function SDloboundCheckbox_Callback(hObject, ~, handles)
% hObject    handle to SDloboundCheckbox (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);



function clusterEdit_Callback(hObject, ~, handles)
% hObject    handle to clusterEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function clusterEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to clusterEdit (see GCBO)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tauRangeCheckbox.
function tauRangeCheckbox_Callback(hObject, ~, handles)
% hObject    handle to tauRangeCheckbox (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);



function upboundEdit19_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit19 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit19 = str2double(get(hObject,'String'));
if isnan(upboundEdit19)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit19,'String',upboundEdit19);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit19_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit20_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit20 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit20 = str2double(get(hObject,'String'));
if isnan(upboundEdit20)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit20,'String',upboundEdit20);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit20_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit21_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit21 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit21 = str2double(get(hObject,'String'));
if isnan(upboundEdit21)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit21,'String',upboundEdit21);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit21_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit22_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit22 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit22 = str2double(get(hObject,'String'));
if isnan(upboundEdit22)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit22,'String',upboundEdit22);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit22_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit23_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit23 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit23 = str2double(get(hObject,'String'));
if isnan(upboundEdit23)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit23,'String',upboundEdit23);
    
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function upboundEdit23_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upboundEdit24_Callback(hObject, ~, handles)
% hObject    handle to upboundEdit24 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

upboundEdit24 = str2double(get(hObject,'String'));
if isnan(upboundEdit24)
    msgbox('Error: at least one of the upper bound entries is NaN', 'Error', 'Error');
else
    set(handles.upboundEdit24,'String',upboundEdit24);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function upboundEdit24_CreateFcn(hObject, ~, ~)
% hObject    handle to upboundEdit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit19_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit19 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit19 = str2double(get(hObject,'String'));
if isnan(loboundEdit19)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit19,'String',loboundEdit19);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit19_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit20_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit20 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit20 = str2double(get(hObject,'String'));
if isnan(loboundEdit20)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit20,'String',loboundEdit20);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit20_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit21_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit20 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit21 = str2double(get(hObject,'String'));
if isnan(loboundEdit21)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit21,'String',loboundEdit21);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit21_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit22_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit20 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit22 = str2double(get(hObject,'String'));
if isnan(loboundEdit22)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit22,'String',loboundEdit22);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit22_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit23_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit23 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit23 = str2double(get(hObject,'String'));
if isnan(loboundEdit23)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit23,'String',loboundEdit23);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit23_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loboundEdit24_Callback(hObject, ~, handles)
% hObject    handle to loboundEdit24 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

loboundEdit24 = str2double(get(hObject,'String'));
if isnan(loboundEdit24)
    msgbox('Error: at least one of the lower bound entries is NaN', 'Error', 'Error');
else
    set(handles.loboundEdit24,'String',loboundEdit24);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function loboundEdit24_CreateFcn(hObject, ~, ~)
% hObject    handle to loboundEdit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
