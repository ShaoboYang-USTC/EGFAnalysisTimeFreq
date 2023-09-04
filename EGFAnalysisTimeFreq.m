%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Wave EGF/CF Phase Velocity Dipersion Analysis Software
% The software can do individual, semiautomated or automated phase velocity
% dispersion measurements.
% by Huajian Yao, 2009 at MIT, hjyao@mit.edu, huajianyao@gmail.com
% For reference, please use: Yao, van der Hilst, de Hoop, 2006, GJI
% Acknowledgements: Hui Huang designed part of the codes for fully
%                   automated phase v dispersion pick by analyzing the
%                   continuity of the dispersion curve, signal to noise
%                   ratio, and behavior of the dispersion curve (increasing 
%                   or decreasing rate, ...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------History of modification------------------------
% Modified: 2009/11/4
% Modified: 2010/02/26
%           (1) allow for the input of CF, the previous version only for
%               EGFs. The code takes Hilbert Transform of CF to obtain EGF.
% Modified: 2010/03/03
%           (1) all for the input of multiple click points used as
%               constraints for the automatic pick of the dispersion curve.
%               This will greatly help to improve completeness of the dispersion
%               measurements when dispersion branch is not that continuous.
%               Left mouse button picks points and right mouse button picks 
%               the last point when the The cursor is shown on the software
%               interface.
% Modified: 2010/03/04
%            (1) correct for a potential bug when the waveform time length
%                is shorter than the time corresponding to set min. group 
%                velocity
% Modified: 2010/09/07
%           change the way SNR is calculated by considering the time length
%           differences of signal and noise windows when using pmtm
%           spectral energy density calculations.
% Modified: 2010/09/14
%           allow for stop the semiautomatic processing by click 'Stop' on
%           the 'Revise dispersion data' question dialog box instead of
%           clicking 'Stop Processing' radio button.
% MAJOR Modification: 2010/09/20    at Scripps Inst. of Oceanography, UCSD
%           Allow for automatic/semiautomatic/individual processing of
%           group velocity dispersion analysis, and then use the group
%           velocity measurements (group travel time) for time-variable
%           filtiering analysis of phase velocities when using 'Time
%           Domain' option in the 'Methods' panel, which improves phase
%           velocity dispersion quality. If no group velocity dispersion
%           picked, the codes also allow for the input of period-dependent
%           group velocity window (compared to fixed group velocity window)
%           . 
%           The new codes also peform new analysis of signal to noise ratio
%           estimation: SNR(f) = (maximum amplitude of envelope around 
%           frequency f in the signal window) / (mean amplitude of envelope
%           of 150 s long noise window right after signal window).       
%           This new SNR method is better than pmtm method in matlab, which
%           sometimes depends on the length of time window used for SNR
%           analysis.
%Modified:  2010/10/22,2010/11/1  Scripps, UCSD
%          fix several warning and potential bugs; add new criteria to
%          judge the quality of group and phase v dispersion curve
%Modified: 2015/06/08 USTC
%          plot the group velocity dispersion curve on the phase velocity
%          image as a reference
%Modified: 2020/08 USTC, by Ning Gu
%          The Automatic dispersion analysis can output the group and phase
%          dispersion images. These are the input of DisperPicker (Yang et
%          al., 2022, SRL) which can automatically pick the dispersion
%          curves based on CNN.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = EGFAnalysisTimeFreq(varargin)
% EGFANALYSISTIMEFREQ M-file for EGFAnalysisTimeFreq.fig
%      EGFANALYSISTIMEFREQ, by itself, creates a new EGFANALYSISTIMEFREQ or raises the existing
%      singleton*.
%
%      H = EGFANALYSISTIMEFREQ returns the handle to a new EGFANALYSISTIMEFREQ or the handle to
%      the existing singleton*.
%
%      EGFANALYSISTIMEFREQ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EGFANALYSISTIMEFREQ.M with the given input arguments.
%
%      EGFANALYSISTIMEFREQ('Property','Value',...) creates a new EGFANALYSISTIMEFREQ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EGFAnalysisTimeFreq_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EGFAnalysisTimeFreq_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EGFAnalysisTimeFreq

% Last Modified by GUIDE v2.5 19-Sep-2010 01:13:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EGFAnalysisTimeFreq_OpeningFcn, ...
                   'gui_OutputFcn',  @EGFAnalysisTimeFreq_OutputFcn, ...
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


% --- Executes just before EGFAnalysisTimeFreq is made visible.
function EGFAnalysisTimeFreq_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EGFAnalysisTimeFreq (see VARARGIN)

% Choose default command line output for EGFAnalysisTimeFreq
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EGFAnalysisTimeFreq wait for user response (see UIRESUME)
% uiwait(handles.EGFAnalysis);


% --- Outputs from this function are returned to the command line.
function varargout = EGFAnalysisTimeFreq_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%
% --- Executes during object creation, after setting all properties.
function EGFAnalysis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EGFAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% call data struct creation function
DataStructCreate;

%%
function DataStructCreate

global GreenFcnInfo filter cross
global g_DataDirectory g_DispDirectory
global g_PrevStaLonLat

g_DataDirectory = pwd;
g_DispDirectory = pwd;
g_PrevStaLonLat = [0 0;0 0];

GreenFcnInfo = struct('Lat1',0,...
    'Lon1',0,...
    'Lat2',0,...
    'Lon2',0,...
    'Name1','',...
    'Name2','',...
    'StaDist',0,...
    'GreenFcn',zeros(2,1),...
    'Time',zeros(1,1),...
    'PtNum',0,...
    'SampleF',0,...
    'SampleT',0,...
    'WinMinV',0,...
    'WinMaxV',0,...
    'DeltaV',0,...
    'PhaseImg1',zeros(1,1),...
    'PhaseImg2',zeros(1,1),...
    'PhaseImg3',zeros(1,1),...
    'GroupImg1',zeros(1,1),...
    'GroupImg2',zeros(1,1),...
    'GroupImg3',zeros(1,1),...
    'PhaseVImg',zeros(1,1),...
    'GroupVImg',zeros(1,1),...
    'PhaseVDisp',zeros(1,1),...
    'GroupVDisp',zeros(1,1),...
    'FileName', '');

% filter for time domain EGF analysis
FilterInfo = struct('Mode',0,...
    'Domain',0,...
    'Window',0,...
    'CtrT',0,...
    'LowT',0,...
    'HighT',0,...
    'CtrF',0,...
    'LowF',0,...
    'HighF',0,...
    'SampleF',0,...
    'SampleT',0,...
    'Length',0,...
    'BandWidth',0,...
    'KaiserPara',1,...
    'GaussAlfa',2.5,...
    'Data',zeros(1,100));

% Period range and interval for narrow band pass filtering
CrossInfo = struct('StartT',0,...
    'EndT',0,...
    'StartF',0,...
    'EndF',0,...
    'DeltaT',0,...
    'DeltaF',0,...
    'NumCtrT',0);

filter = struct(FilterInfo);
cross = struct(CrossInfo);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function EGFAnalysis_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to EGFAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global g_ClickPoint
g_ClickPoint = get(handles.axes2, 'CurrentPoint');

%%
% --- Executes on button press in StationFileInput.
function StationFileInput_Callback(hObject, eventdata, handles)
% hObject    handle to StationFileInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global g_allsta

[pfile1, pname1, index] = uigetfile({'*.txt';'*.dat';'*.*'},'Open file containing all stations information', pwd);
files = strcat(pname1,pfile1);

fstat = fopen(files,'r');
g_allsta = struct('name', {}, 'net', {}, 'lat', {}, 'lon', {});
% read and plotwindowedegf stations

h1 = handles.axes1;
hold(h1,'off');
set(gcf,'CurrentAxes',h1)
load topomap
x=linspace(0,360,1080);
y=linspace(-90,90,2160);
[cmap clim] = demcmap(topomap);
hold on
imagesc(x,y,topomap,clim);%colorbar('vert');
colormap(cmap);axis image; grid on; axis on;

load coast % the coast line is lat, long
kk = find(long < 0);
long(kk) = 360 - abs(long(kk));
ii = find( abs(long) < 0.5);
long(ii) = NaN;
lat(ii) = NaN;
hold(h1,'on');
plot(h1, long,lat,'k', 'LineWidth',2);
set(gca,'ydir','normal');
set(gca,'FontSize', 6, 'FontWeight','bold');
axis equal

i=1;
while ~feof(fstat)
    name = fscanf(fstat,'%s',1);
    if ~strcmp(name,'')
	   g_allsta(i).name = name; %station name
%  	   g_allsta(i).net = fscanf(fstat,'%s',1); % station network
       g_allsta(i).lon = fscanf(fstat,'%f',1); % station latitude      
	   g_allsta(i).lat = fscanf(fstat,'%f',1); % station longitude
       
       if g_allsta(i).lon < 0
           g_allsta(i).lon = 360 + g_allsta(i).lon;
       end
       Lon(i) = g_allsta(i).lon;
       Lat(i) = g_allsta(i).lat;
       
       hold(h1, 'on');
       plot(h1, g_allsta(i).lon, g_allsta(i).lat, 'k^', 'MarkerSize',6, 'MarkerFaceColor','k');
    else
        break
    end
    temp = fgetl(fstat);
    i=i+1;
end

% nn = find(Lon < 0);
% Lon(nn) = 360 - abs(Lon(nn));

ExtraLon = 0.1*(max(Lon)-min(Lon));
ExtraLat = 0.1*(max(Lat)-min(Lat));
ExtraLon = min(ExtraLon, 1);
ExtraLat = min(ExtraLat, 1);
xlim(h1,[min(Lon)-ExtraLon max(Lon)+ExtraLon]);
ylim(h1,[min(Lat)-ExtraLat max(Lat)+ExtraLat]);

% title('Station Location', 'FontSize',8,'FontWeight','bold');
stanum=i-1;
fclose(fstat);

%%
% --- Executes on button press in WaveformFileInput.
function WaveformFileInput_Callback(hObject, eventdata, handles)
% hObject    handle to WaveformFileInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global g_DataDirectory g_allsta g_EGFdata g_EGForCFIndex

h1 = handles.axes1;
hold(h1,'on');
set(gcf,'CurrentAxes',h1);

stanum = length(g_allsta);
g_EGFdata = struct('fname', '');

if g_DataDirectory == 0
    g_DataDirectory = pwd;
end

[pfile1, pname1, index] = uigetfile({'*.txt';'*.dat';'*.*'},'Open file containing all data file names', g_DataDirectory);
files = strcat(pname1,pfile1);
g_DataDirectory = pname1;

fname = fopen(files, 'r');

i = 1;
while ~feof(fname)
    name = fscanf(fname, '%s', 1);
    if ~strcmp(name,'')
        g_EGFdata(i).fname = name;
    else
        break
    end
    temp = fgetl(fname);
    i = i + 1;
end
fclose(fname);

datatype = questdlg('Please select input data type (EGF or CF):','Data Type', 'EGF', 'CF', 'EGF');
if strcmp(datatype, 'EGF')
    g_EGForCFIndex = 1;
elseif strcmp(datatype, 'CF')
    g_EGForCFIndex = 2;
end

filenum = length(g_EGFdata);
if filenum > 0
    set(handles.MsgEdit, 'String', 'Read EGF/CF file names successfully!');
end


%%
% --- Executes on button press in DisperFolder.
function DisperFolder_Callback(hObject, eventdata, handles)
% hObject    handle to DisperFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global g_DispDirectory

g_DispDirectory = uigetdir(pwd, 'Select folder to save dispersion data:');



function MsgEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MsgEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MsgEdit as text
%        str2double(get(hObject,'String')) returns contents of MsgEdit as a double


% --- Executes during object creation, after setting all properties.
function MsgEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MsgEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function GFcnWinMinV_Callback(hObject, eventdata, handles)
% hObject    handle to GFcnWinMinV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GFcnWinMinV as text
%        str2double(get(hObject,'String')) returns contents of GFcnWinMinV as a double



function new_Callback(hObject, eventdata, handles)
% hObject    handle to new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of new as text
%        str2double(get(hObject,'String')) returns contents of new as a double


% --- Executes during object creation, after setting all properties.
function new_CreateFcn(hObject, eventdata, handles)
% hObject    handle to new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StartPeriod_Callback(hObject, eventdata, handles)
% hObject    handle to StartPeriod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StartPeriod as text
%        str2double(get(hObject,'String')) returns contents of StartPeriod as a double


% --- Executes during object creation, after setting all properties.
function StartPeriod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StartPeriod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EndPeriod_Callback(hObject, eventdata, handles)
% hObject    handle to EndPeriod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EndPeriod as text
%        str2double(get(hObject,'String')) returns contents of EndPeriod as a double


% --- Executes during object creation, after setting all properties.
function EndPeriod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EndPeriod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function DeltaPeriod_Callback(hObject, eventdata, handles)
% hObject    handle to DeltaPeriod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DeltaPeriod as text
%        str2double(get(hObject,'String')) returns contents of DeltaPeriod as a double


% --- Executes during object creation, after setting all properties.
function DeltaPeriod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DeltaPeriod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function GFcnWinMaxV_Callback(hObject, eventdata, handles)
% hObject    handle to GFcnWinMaxV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GFcnWinMaxV as text
%        str2double(get(hObject,'String')) returns contents of GFcnWinMaxV as a double


% --- Executes during object creation, after setting all properties.
function GFcnWinMaxV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GFcnWinMaxV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function GFcnWinMinV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GFcnWinMinV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TimeDomain.
function TimeDomain_Callback(hObject, eventdata, handles)
% hObject    handle to TimeDomain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TimeDomain
if get(handles.TimeDomain, 'Value')
    set(handles.FreqDomain, 'Value', 0);
    set(handles.AkiSPAC, 'Value', 0);    
else
    set(handles.TimeDomain, 'Value', 1);
end
set(handles.MsgEdit,'String','Set "Time Domain Filter" box!');

% --- Executes on button press in FreqDomain.
function FreqDomain_Callback(hObject, eventdata, handles)
% hObject    handle to FreqDomain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FreqDomain
if get(handles.FreqDomain, 'Value')
    set(handles.TimeDomain, 'Value', 0);
    set(handles.AkiSPAC, 'Value', 0);    
else
    set(handles.FreqDomain, 'Value', 1);
end

% --- Executes on button press in AkiSPAC.
function AkiSPAC_Callback(hObject, eventdata, handles)
% hObject    handle to AkiSPAC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AkiSPAC
if get(handles.AkiSPAC, 'Value')
    set(handles.TimeDomain, 'Value', 0);
    set(handles.FreqDomain, 'Value', 0);    
else
    set(handles.AkiSPAC, 'Value', 1);
end


% --- Executes on button press in AutoProcessIndex.
function AutoProcessIndex_Callback(hObject, eventdata, handles)
% hObject    handle to AutoProcessIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AutoProcessIndex


% --- Executes on button press in ProcessSelectedData.
function ProcessSelectedData_Callback(hObject, eventdata, handles)
% hObject    handle to ProcessSelectedData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ProcessSelectedData
if get(handles.ProcessSelectedData, 'Value')
    set(handles.ProcessAllData, 'Value', 0);
else
    set(handles.ProcessAllData, 'Value', 1);
end

% --- Executes on button press in ProcessAllData.
function ProcessAllData_Callback(hObject, eventdata, handles)
% hObject    handle to ProcessAllData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ProcessAllData
if get(handles.ProcessAllData, 'Value')
    set(handles.ProcessSelectedData, 'Value', 0);
else
    set(handles.ProcessSelectedData, 'Value', 1);
end


function DataEndIndex_Callback(hObject, eventdata, handles)
% hObject    handle to DataEndIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataEndIndex as text
%        str2double(get(hObject,'String')) returns contents of DataEndIndex as a double


% --- Executes during object creation, after setting all properties.
function DataEndIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataEndIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DataStartIndex_Callback(hObject, eventdata, handles)
% hObject    handle to DataStartIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataStartIndex as text
%        str2double(get(hObject,'String')) returns contents of DataStartIndex as a double


% --- Executes during object creation, after setting all properties.
function DataStartIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataStartIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in StopDataProcessing.
function StopDataProcessing_Callback(hObject, eventdata, handles)
% hObject    handle to StopDataProcessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of StopDataProcessing

%%
% --- Executes on button press in StartDataProcessing.
function StartDataProcessing_Callback(hObject, eventdata, handles)
% hObject    handle to StartDataProcessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global GreenFcnInfo gfcn g_allsta g_EGFdata g_PrevStaLonLat g_FileName
global g_DataDirectory 

g_PrevStaLonLat = [0 0;0 0];

stanum = length(g_allsta);
filenum = length(g_EGFdata);

h1 = handles.axes1;
h2 = handles.axes2;

if get(handles.ProcessAllData, 'Value')    
    StartIndex = 1;
    EndIndex = filenum;
    set(handles.DataStartIndex, 'String', num2str(StartIndex));
    set(handles.DataEndIndex, 'String', num2str(EndIndex));
else
    StartIndex = str2num(get(handles.DataStartIndex, 'String'));
    EndIndex = str2num(get(handles.DataEndIndex, 'String'));
end

set(handles.StopDataProcessing, 'Value', 0);

for i = StartIndex:EndIndex
    
    gfcn = struct(GreenFcnInfo);
    if get(handles.StopDataProcessing,'Value')
       set(handles.MsgEdit,'String',['Data Processing Stopped! Path Index = ' num2str(i-1)]);
       break
    end    
        
    display(['Data Processing ... Path Index = ' num2str(i)]);
    
    % read EGF data file
    datafile = strcat(g_DataDirectory, g_EGFdata(i).fname);
    

    g_FileName = g_EGFdata(i).fname;
    maxamp = RdGreenFcn(datafile);
    
%     % reset the length of EGF for better showing the waveform
%     gfcn.WinMinV = str2num(get(handles.GFcnWinMinV,'String'));
%     maxtime = gfcn.StaDist/gfcn.WinMinV + 25;
%     ptnum = gfcn.PtNum;
%     gfcn.PtNum = min(round(maxtime*gfcn.SampleF)+1, ptnum);

    % obtain the name of two stations
    for kk = 1:stanum
        if abs(g_allsta(kk).lon - gfcn.Lon1)<1e-4 && abs(g_allsta(kk).lat - gfcn.Lat1)<1e-4
        
            gfcn.Name1 = g_allsta(kk).name;
            break
        end
    end
    for kk = 1:stanum
        if abs(g_allsta(kk).lon - gfcn.Lon2)<1e-4 && abs(g_allsta(kk).lat - gfcn.Lat2)<1e-4
        
            gfcn.Name2 = g_allsta(kk).name;
            break
        end
    end
    
    % updata message box information
    UpdataMsgBoxInfo(hObject, eventdata, handles);
    
    ProcessIndex = 0;
    if ~isnan(maxamp) && ~strcmp(gfcn.Name1, gfcn.Name2)
        if maxamp > 0
            ProcessIndex = 1;
        end        
    end   
        
    if ProcessIndex == 1
        % main function for two-station data processing
        EGFDataProcessing(hObject, eventdata, handles);
        set(handles.MsgEdit,'String',strcat('Data Processing is continuing: file index = ', num2str(i+1)));

        pause(0.05)
        
    end
    
%     clear gfcn
    
end

if ~get(handles.StopDataProcessing,'Value')
    set(handles.MsgEdit,'String','Data Processing Finished!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main function for two-station data processing
function EGFDataProcessing(hObject, eventdata, handles)

global gfcn filter cross g_HighSNRIndex
global g_PrevStaLonLat
global g_DisperType

gfcn.WinMinV = str2num(get(handles.GFcnWinMinV,'String'));
gfcn.WinMaxV = str2num(get(handles.GFcnWinMaxV,'String'));

cross.StartT = str2num(get(handles.StartPeriod,'String'));
cross.EndT = str2num(get(handles.EndPeriod,'String'));
cross.DeltaT = str2num(get(handles.DeltaPeriod,'String'));
cross.NumCtrT = round((cross.EndT - cross.StartT)/cross.DeltaT)+ 1;
g_WinMaxPtNum = round(gfcn.SampleF*gfcn.StaDist/gfcn.WinMinV) + 1;
g_WinMinPtNum = round(gfcn.SampleF*gfcn.StaDist/gfcn.WinMaxV) + 1;
if g_WinMaxPtNum > gfcn.PtNum
    g_WinMaxPtNum = gfcn.PtNum-1;
    gfcn.WinMinV = ceil(10*gfcn.StaDist/gfcn.Time(end))/10;
    set(handles.MsgEdit,'String',strcat('Min velocity reset to ',num2str(gfcn.WinMinV)));
end
time = ((g_WinMinPtNum-1):1:(g_WinMaxPtNum-1))/gfcn.SampleF;
SignalLength = length(time);

gfcn.DeltaV = 0.005;
VPoint = gfcn.WinMinV:gfcn.DeltaV:gfcn.WinMaxV;
VImgPt = length(VPoint);
TPoint = cross.StartT:cross.DeltaT:cross.EndT;
TravPtV = gfcn.StaDist./time;
gfcn.Time = (0:(gfcn.PtNum-1))*gfcn.SampleT;


%  function to create signal window function 
[Window, TaperNum] = MakeSignalWindow(hObject, eventdata, handles);
WavePtNum = min((g_WinMaxPtNum + TaperNum), gfcn.PtNum);
WinWave = zeros(1, WavePtNum);
%  window the stacked EGF (causal + acausal)
stackEGF = gfcn.GreenFcn(1,:) + gfcn.GreenFcn(2,:);
WinWave = stackEGF(1:WavePtNum).*Window(1:WavePtNum);
WinWave = WinWave/max(abs(WinWave));
stackEGF = stackEGF/max(abs(WinWave));

% plot two stations in processing
h1 = handles.axes1;
set(gcf,'CurrentAxes',h1);
hold(h1,'on');
if sum(sum(g_PrevStaLonLat)) ~= 0
    plot(h1, g_PrevStaLonLat(1,1), g_PrevStaLonLat(1,2), 'k^', 'MarkerSize',6, 'MarkerFaceColor','k');
    hold(h1,'on');
    plot(h1, g_PrevStaLonLat(2,1), g_PrevStaLonLat(2,2), 'k^', 'MarkerSize',6, 'MarkerFaceColor','k');
    hold(h1,'on');
    plot(h1, gfcn.Lon1, gfcn.Lat1, 'r^', 'MarkerSize',6, 'MarkerFaceColor','r');
    hold(h1, 'on');
    plot(h1, gfcn.Lon2, gfcn.Lat2, 'b^', 'MarkerSize',6, 'MarkerFaceColor','b');
else
    plot(h1, gfcn.Lon1, gfcn.Lat1, 'r^', 'MarkerSize',6, 'MarkerFaceColor','r');
    hold(h1, 'on');
    plot(h1, gfcn.Lon2, gfcn.Lat2, 'b^', 'MarkerSize',6, 'MarkerFaceColor','b');
end


% define 50s long noise window after the windowed surface wave
NoisePt = round(50/gfcn.SampleT);
NoiseWinWave = zeros(1, NoisePt);
if ((g_WinMaxPtNum + TaperNum) + (cross.EndT/gfcn.SampleT))< gfcn.PtNum
    nn = (g_WinMaxPtNum + TaperNum + 1):min((g_WinMaxPtNum + TaperNum + NoisePt),gfcn.PtNum);
    NoiseWinWave(1:length(nn)) = stackEGF(nn);
    NoiseIndex = 1;
    NoiseLength = length(nn);
else
    NoiseIndex = 0;
end



% plot signal (blue) and noise (red)
h3 = handles.axes3;
set(gcf,'CurrentAxes',h3);
hold(h3,'off');
timeshow = (0:min((g_WinMaxPtNum + TaperNum + NoisePt),gfcn.PtNum)-1)*gfcn.SampleT;
plot(timeshow, stackEGF(1:length(timeshow))+2,'k');
hold on; plot((1:length(WinWave))*gfcn.SampleT, WinWave);
if NoiseIndex == 1   
    hold on; plot((nn(1)-1:(nn(1)+NoisePt-2))*gfcn.SampleT, NoiseWinWave,'r');
end

% ylim([-1.5 3.5]);
set(gca, 'FontSize', 8);

% calculate envelope images for signal and noise and estimate SNR
% SNR(T) =  max(signal envelope at period T)/mean(noise envelope at period T)
h4 = handles.axes4;
set(gcf,'CurrentAxes',h4);
hold(h4,'off');
minSNR = str2num(get(handles.minSNR, 'String'));
SNR_T = zeros(1, cross.NumCtrT);

EnvelopeImageSignal = EnvelopeImageCalculation(WinWave, gfcn.SampleF, TPoint, gfcn.StaDist);
AmpS_T = max(EnvelopeImageSignal,[],2);

if NoiseIndex == 1  % noise window long enough
    EnvelopeImageNoise = EnvelopeImageCalculation(NoiseWinWave.*tukeywin(NoisePt,0.2)', gfcn.SampleF, TPoint, gfcn.StaDist);
    for i = 1:cross.NumCtrT
        SNR_T(i) = AmpS_T(i)/mean(EnvelopeImageNoise(i,1:NoiseLength));
    end
    semilogy(TPoint, SNR_T);
    HighSNRIndex = find(SNR_T > minSNR); % find SNR > MinSNR;
    g_HighSNRIndex = HighSNRIndex;
    hold on;
    semilogy(TPoint(HighSNRIndex), SNR_T(HighSNRIndex), 'r*');
    grid on
    xlim([cross.StartT cross.EndT]);
    xlabel('Period (s)','FontSize', 8);
    ylabel('SNR','FontSize', 8); 
    if mean(SNR_T) < minSNR/2 || max(SNR_T) < minSNR || length(HighSNRIndex) < 0.1*cross.NumCtrT
        IsDispGood = 0; % overall SNR too low, do not pick dispersion
    else
        IsDispGood = 1; % pick dispersion
    end 
else  % noise window too short
    semilogy(TPoint, AmpS_T);
    xlim([cross.StartT cross.EndT]);
    xlabel('Period (s)','FontSize', 8);
    ylabel('PSD of Signal', 'Color','r','FontSize', 8);
    title('Noise window does not exit!', 'Color','r', 'FontSize',8);    
    HighSNRIndex = 1:cross.NumCtrT; % set high SNR for all periods
    IsDispGood = 1;  % pick dispersion
end
set(gca, 'FontSize', 8);


% % estimate Signal to Noise Ratio using multi-taper method in matlab
% fftNumPt = 2^nextpow2(2^8*gfcn.SampleF);
% [PxxS,f] = pmtm(WinWave,5/2,fftNumPt,gfcn.SampleF);
% if NoiseIndex == 1
%     [PxxN,f] = pmtm(NoiseWinWave,5/2,fftNumPt,gfcn.SampleF);
%     mm = find(PxxN == 0); nn = find(PxxN > 0);
%     PxxN(mm) = min(PxxN(nn));
%     SNR_f = PxxS./PxxN*NoiseLength/SignalLength;
%     SNR_T = interp1(f, SNR_f, 1./TPoint, 'cubic');   
% end
% 
% h4 = handles.axes4;
% set(gcf,'CurrentAxes',h4);
% hold(h4,'off');
% minSNR = str2num(get(handles.minSNR, 'String'));
% if NoiseIndex == 1
%     semilogy(TPoint, SNR_T);
%     HighSNRIndex = find(SNR_T > minSNR); % find SNR > MinSNR;
%     g_HighSNRIndex = HighSNRIndex;
%     hold on;
%     semilogy(TPoint(HighSNRIndex), SNR_T(HighSNRIndex), 'r*');
%     grid on
%     xlabel('Period (s)','FontSize', 8);
%     ylabel('SNR','FontSize', 8);
%     if mean(SNR_T) < minSNR/2
%         IsDispGood = 2; % overall SNR too low, do not pick dispersion
%     else
%         IsDispGood = 1; % pick dispersion
%     end
% else
%     PxxS_T = interp1(f, PxxS, 1./TPoint, 'v4');
%     semilogy(TPoint, PxxS_T);
% %     meanPxxS = mean(PxxS_T);
% %     HighSNRIndex = find(PxxS_T > meanPxxS);
% %     hold on; semilogy(TPoint(HighSNRIndex), PxxS_T(HighSNRIndex), 'r*');
%     HighSNRIndex = 1:cross.NumCtrT;
%     
%     xlabel('Period (s)','FontSize', 8);
%     ylabel('PSD of Signal', 'Color','r','FontSize', 8);
%     title('Noise window does not exit!', 'Color','r', 'FontSize',8);
%     IsDispGood = 1;  % pick dispersion
% end
% set(gca, 'FontSize', 8);


SNRIndex = zeros(1, cross.NumCtrT);
SNRIndex(HighSNRIndex) = 1;
SNR_T = zeros(1, cross.NumCtrT);

% for those not so bad pts, if they are in the middle of good pts, accept them
for ii = 2:length(SNRIndex)-1
    if SNRIndex(ii) == 0
        if SNR_T(ii) > minSNR/2 && SNRIndex(ii-1) == 1 && SNRIndex(ii+1) == 1
            SNRIndex(ii)=1;
        end
    else
        SNRIndex(ii)=1;
    end
end


if IsDispGood == 1   % good data quality
    
    
    if get(handles.GroupVIndex, 'Value')  % also obtain group velocity dispersion curve
        timeptnum = g_WinMinPtNum:1:g_WinMaxPtNum;
        for i = 1:cross.NumCtrT
            gfcn.GroupVImg(1:VImgPt, i) = interp1(TravPtV, EnvelopeImageSignal(i,timeptnum)'/AmpS_T(i), VPoint, 'spline');
        end
        minamp = min(min(gfcn.GroupVImg(1:VImgPt,1:cross.NumCtrT)));
        
        
        g_DisperType = 2; % for group velocity dispersion
        
        while IsDispGood == 1
            h2 = handles.axes2;
            set(gcf,'CurrentAxes',h2);
            hold(h2,'off');
            imagesc(TPoint, VPoint, gfcn.GroupVImg(1:VImgPt,1:cross.NumCtrT), [minamp, 1]);
            
            set(gca,'YDir','normal');
            xlabel('Period (s)', 'FontSize', 10, 'FontWeight', 'bold');
            ylabel('Group Velocity (km/s)', 'FontSize', 10, 'FontWeight', 'bold');
            set(gca, 'FontSize', 10, 'FontWeight', 'bold');
            
           
                
%             if cross.EndT < 4
%                 set(gca, 'XTick',0:0.5:4,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%             elseif cross.EndT >= 4 && cross.EndT <= 20
%                 set(gca, 'XTick',0:2:20,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%             elseif cross.EndT > 20
%                 set(gca, 'XTick',0:5:cross.EndT,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%             end  

%%%%%%%%%%%%%% -- modified by Ning Gu,2019.04.17 --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  global  g_FileName g_DispDirectory 
               
                 AI_file = strcat(g_DispDirectory, '/','G.',g_FileName);   
                 dlmwrite(AI_file,gfcn.GroupVImg,'delimiter',' ');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

            if get(handles.AutoProcessIndex, 'Value')==1  % for automatic pick of dispersion curve
                IsDispGood = AutoGroupDisperProcessing(hObject, handles, TPoint, VPoint, gfcn.GroupVImg, SNRIndex, SNR_T, minSNR);
           
%              global  g_FileName g_DispDirectory 
               
                AI_file = strcat(g_DispDirectory, '/','G.',g_FileName);   
                dlmwrite(AI_file,gfcn.GroupVImg,'delimiter',' ');
            else  % for semiautomatic pick of dispersion curve; user need to pick dispersion curve by user           
                [IsDispGood, StartTIndex, EndTIndex] = DisperPick(hObject, handles); % obtain the disperion curve from user mouse input points
                

            end
        end

        % if the group velocity dispersion curve has been picked
        if IsDispGood == 2
            GroupTimeIndex = 1;
            GroupTime = gfcn.StaDist./gfcn.GroupVDisp(1:cross.NumCtrT);
        else
            GroupTimeIndex = 0; % no group dispersion/travel time picked
            GroupTime = NaN*zeros(1,cross.NumCtrT);
        end

        hold off
        pause(0.05);
    else
        GroupTimeIndex = 0; % no group dispersion/travel time picked
        GroupTime = NaN*zeros(1,cross.NumCtrT);
    end
    
    if IsDispGood ~= 0  % IsDispGood = 0 is for the case that quality of
                        % group velocity dispersion is too bad and no
                        % group v disperion saved, or the processing is
                        % stopped
        g_DisperType = 1; % for phase velocity dispersion
        IsDispGood = 1; % reset IsDispGood = 1 for next step phase v dispersion analysis
        % function to calculate phase v (c-T) image matrix
        PhaseVImageCalculation(hObject, eventdata, handles, WinWave, GroupTimeIndex, GroupTime);
    end
end

% phase velocity dispersion pick
while IsDispGood == 1
    % plot c-T image
    h2 = handles.axes2;
    set(gcf,'CurrentAxes',h2);
    hold(h2,'off');
    imagesc(TPoint, VPoint, gfcn.PhaseVImg); 
    % cmap = colormap('gray');
    % colormap(cmap(end:-1:1,:)); colorbar;

    set(gca,'YDir','normal');
    set(gca,'YDir','normal','FontSize', 8, 'FontWeight', 'bold','FontName','Arial');
    xlabel('Period (s)', 'FontSize', 8, 'FontWeight', 'bold','FontName','Arial');
    ylabel('Phase Velocity (km/s)', 'FontSize', 8, 'FontWeight', 'bold','FontName','Arial');
%     if cross.EndT < 4
%         set(gca, 'XTick',0:0.5:4,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%     elseif cross.EndT >= 4 && cross.EndT <= 20
%         set(gca, 'XTick',0:2:20,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%     elseif cross.EndT > 20
%         set(gca, 'XTick',0:5:cross.EndT,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%     end
    
%%%%%%%%%%%%%% -- modified by Ning Gu,2019.04.17 --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 global  g_FileName g_DispDirectory 
               
                AI_file = strcat(g_DispDirectory, '/','C.',g_FileName);   
                dlmwrite(AI_file,gfcn.PhaseVImg,'delimiter',' ');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % if group disper picked, then plot group disper on phase v image as a reference
    if get(handles.GroupVIndex, 'Value') && GroupTimeIndex == 1  
        hold on; plot(TPoint, gfcn.GroupVDisp, 'm-', 'LineWidth',2);
        
   %%%%%%%%%% guning %%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%
%         IIGood = find(SNRIndex == 1);       
%          hold on; plot(TPoint(StartTIndex:EndTIndex), gfcn.GroupVDisp(StartTIndex:EndTIndex), 'co', 'MarkerSize', 8);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
    
    if get(handles.AutoProcessIndex, 'Value')==1  % for automatic pick of dispersion curve
        
    %%%% guning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
         IsDispGood = AutoPhaseDisperProcessing(hObject, handles, TPoint, VPoint, gfcn.PhaseVImg, SNRIndex, SNR_T, minSNR);
%          IsDispGood = 2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    else  % for semiautomatic pick of dispersion curve; user need to pick dispersion curve by user
        [IsDispGood, StartTIndex, EndTIndex] = DisperPick(hObject, handles); % obtain the disperion curve
    end
end

hold off


%% --- calculate envelope image, i.e., to obtain envelope at each T
function EnvelopeImage = EnvelopeImageCalculation(WinWave, fs, TPoint, StaDist)
%% new code for group velocity analysis using frequency domain Gaussian filter
alfa = [0 100 250 500 1000 2000 4000 20000; 5  8  12  20  25  35  50 75];
guassalfa = interp1(alfa(1,:), alfa(2, :), StaDist);

NumCtrT = length(TPoint);
PtNum = length(WinWave);

nfft = 2^nextpow2(max(PtNum,1024*fs));
xxfft = fft(WinWave, nfft);
fxx = (0:(nfft/2))/nfft*fs; 
IIf = 1:(nfft/2+1);
JJf = (nfft/2+2):nfft;

EnvelopeImage = zeros(NumCtrT, PtNum);
for i = 1:NumCtrT
    CtrT = TPoint(i);
    fc = 1/CtrT;                
    Hf = exp(-guassalfa*(fxx - fc).^2/fc^2);
    yyfft = zeros(1,nfft);
    yyfft(IIf) = xxfft(IIf).*Hf;
    yyfft(JJf) = conj(yyfft((nfft/2):-1:2));
    yy = real(ifft(yyfft, nfft));
    filtwave = abs(hilbert(yy(1:nfft)));
    EnvelopeImage(i, 1:PtNum) = filtwave(1:PtNum);
end


%% --- calculate the T-C (phase v) image from 3 different methods
function PhaseVImageCalculation(hObject, eventdata, handles, WinWave, GroupTimeIndex, GroupTime)
global filter gfcn cross
global g_WinMaxPtNum g_WinMinPtNum
global g_WindnumT g_Winmintime g_refGVwinIndex g_refWinV

df_freq = str2num(get(handles.df_Spectral, 'String')); % frequency interval of the spectrum
fftNumPt = 2^nextpow2(gfcn.SampleF/df_freq); % number of fft point when using freq. domain or Aki's method

TPoint = cross.StartT:cross.DeltaT:cross.EndT;
VPoint = gfcn.WinMinV:gfcn.DeltaV:gfcn.WinMaxV;
VImgPt = length(VPoint); 

WaveNumPt = length(WinWave);

% if no information of reference group velocity window is given for
% time-variable filtering analysis
if isempty(g_refGVwinIndex)
    g_refGVwinIndex = 0;
end
       
% if time domain method, set time-variable filtering parameters
if get(handles.TimeDomain, 'Value')==1   % time domain method 

    if get(handles.AutoProcessIndex, 'Value')==0 % for semi-automatic dispersion analysis
        % select the method for time-variable filtering analysis
        if g_refGVwinIndex == 1 || get(handles.GroupVIndex, 'Value');
            output = questdlg('Please select window type for time-variable filtering analysis:',...
                'Time-varialbe filtering type', 'obs', 'ref', 'no', 'obs');
            if strcmp(output, 'ref')
                if g_refGVwinIndex == 0
                    set(handles.MsgEdit,'String','No ref. group v window, change to uniform group v window');
                    output = 'no';
                else
                    % interpolation to obtain the group v window at each period 
                    GroupVwinMin = interp1(g_refWinV(:,1), g_refWinV(:,3),TPoint, 'pchip');
                    GroupVwinMax = interp1(g_refWinV(:,1), g_refWinV(:,4),TPoint, 'pchip');
                end
            elseif strcmp(output, 'obs')
                if GroupTimeIndex == 0
                    set(handles.MsgEdit,'String','No obs. group v, change to uniform group v window');
                    output = 'no';
                else
                    set(handles.TimeVariableFiltering, 'Value',1); % time-variable filtering
                    g_WindnumT = str2double(get(handles.windowPeriodNum, 'String'));
                    g_Winmintime = str2double(get(handles.minWindowTime, 'String'));
                    GroupVwinMin = gfcn.StaDist./(GroupTime + max(g_WindnumT/2*TPoint,g_Winmintime));
                    GroupVwinMax = gfcn.StaDist./(GroupTime - max(g_WindnumT/2*TPoint,g_Winmintime));
                    III = find(GroupVwinMax <= 0);
                    GroupVwinMax(III) = gfcn.WinMaxV*2;
                end
            else
                output = 'no'; % use uniform group v window
            end
        else
            output = 'no'; % use uniform group v window
        end

    elseif get(handles.AutoProcessIndex, 'Value')==1 % for automatic dispersion analysis
        if (GroupTimeIndex == 1) && get(handles.TimeVariableFiltering, 'Value') % using picked (observed) period-dependent group velocity window
            g_WindnumT = str2double(get(handles.windowPeriodNum, 'String'));
            g_Winmintime = str2double(get(handles.minWindowTime, 'String'));
            GroupVwinMin = gfcn.StaDist./(GroupTime + max(g_WindnumT/2*TPoint,g_Winmintime));
            GroupVwinMax = gfcn.StaDist./(GroupTime - max(g_WindnumT/2*TPoint,g_Winmintime));
            III = find(GroupVwinMax <= 0);
            GroupVwinMax(III) = gfcn.WinMaxV*2; % reset velocity of negative travel time 
            output = 'obs';
        else
            if g_refGVwinIndex == 1 % using input reference period-dependent group velocity window
                % interpolation to obtain the group v window at each period 
                GroupVwinMin = interp1(g_refWinV(:,1), g_refWinV(:,3),TPoint, 'pchip');
                GroupVwinMax = interp1(g_refWinV(:,1), g_refWinV(:,4),TPoint, 'pchip'); 
                output = 'ref';
            else  % use the uniform group velocity window in Settings Panel, i.e., no MFT analysis
                output = 'no';
            end
        end
    end

    % plot time-variable group velocity window at different periods
    if strcmp(output, 'obs') || strcmp(output, 'ref')   
        h2 = handles.axes2;
        set(gcf,'CurrentAxes',h2);
        % reset period-dependent group v winodw: has to be less than
        % 0.98*gfcn.WinMaxV or larger than 1.02*gfcn.WinMinV
        pWinMinV = max(1.02*gfcn.WinMinV*ones(1, cross.NumCtrT), GroupVwinMin);
        pWinMaxV = min(0.98*gfcn.WinMaxV*ones(1, cross.NumCtrT), GroupVwinMax);
        GroupVwinMin = pWinMinV;
        GroupVwinMax = pWinMaxV;

        hold(h2, 'on');plot(TPoint, pWinMinV, 'w--', 'LineWidth', 2);
        hold(h2, 'on'); plot(TPoint, pWinMaxV,'w--', 'LineWidth', 2);        
        pause(0.05);
    end
end

        
% get the method for phase velocity dispersion
if get(handles.TimeDomain, 'Value')==1   % time domain method   
    
    % if do not set time domain filter 
    if filter.SampleF == 0 || filter.BandWidth == 0
        filter.SampleF = gfcn.SampleF;
        filter.BandWidth = cross.DeltaT;
        set(handles.FilterBandWidthT, 'String', num2str(filter.BandWidth));
%         dfmin = 1/(cross.EndT - 0.5*filter.BandWidth) - 1/(cross.EndT + 0.5*filter.BandWidth);
%         filter.Length = 2^nextpow2(gfcn.SampleF/dfmin);
        filter.Length = 2^nextpow2(1024*gfcn.SampleF);
        filter.KaiserPara = 6;
    elseif filter.SampleF ~= gfcn.SampleF
        filter.SampleF = gfcn.SampleF;
    end    

    HalfFilterNum =  round(filter.Length/2);
    WinWave((WaveNumPt + 1):(WaveNumPt + HalfFilterNum)) = 0;
    % band pass filtering
    for numt = 1:cross.NumCtrT
        filter.CtrT = cross.StartT + (numt - 1)*cross.DeltaT;
        filter.CtrF = (2/filter.SampleF)/filter.CtrT;
        
        filter.LowF = (2/filter.SampleF)/(filter.CtrT + 0.5*filter.BandWidth);
        filter.HighF = (2/filter.SampleF)/(filter.CtrT - 0.5*filter.BandWidth);
%         filter.LowF = (2/filter.SampleF)/(filter.CtrT + 0.5*df_freq);
%         filter.HighF = (2/filter.SampleF)/(filter.CtrT - 0.5*df_freq);      
        filter.Data = fir1(filter.Length, [filter.LowF, filter.HighF], kaiser(filter.Length + 1,filter.KaiserPara));
        FilteredWave = zeros(1,WaveNumPt + HalfFilterNum);
        
        %two-pass filtering (time and time-reverse) in order to remove phase shift
        if strcmp(output,'no')   % using fixed wave window within [gfcn.WinMinV gfcn.WinMaxV]
            FilteredWave = fftfilt(filter.Data, WinWave(1:(WaveNumPt + HalfFilterNum)));
        else % using frequency-time variable MFT to window the original waveform           
            if strcmp(output, 'obs')
                winpt = round(max(g_WindnumT*filter.CtrT,g_Winmintime)*filter.SampleF);
                if mod(winpt,2) == 1  % to ensure winpt is even number
                    winpt = winpt + 1;
                end                
                wintukey = tukeywin(winpt, 0.2);      
                grouppt = winpt + round(GroupTime(numt)*filter.SampleF + 1);
            elseif strcmp(output, 'ref')
                winpt = round(gfcn.StaDist/(GroupVwinMax(numt) - GroupVwinMin(numt))*filter.SampleF + 1);
                if mod(winpt,2) == 1  % to ensure winpt is even number
                    winpt = winpt + 1;
                end
                wintukey = tukeywin(winpt, 0.2);
                GroupTime(numt) = (gfcn.StaDist/GroupVwinMax(numt) + gfcn.StaDist/GroupVwinMin(numt))/2;
                if GroupTime(numt) >= WaveNumPt
                    GroupTime(numt) = WaveNumPt - 1;
                end                
                grouppt = winpt + round(GroupTime(numt)*filter.SampleF + 1);
            end
            tmpWave = [zeros(1,winpt) WinWave(1:WaveNumPt) zeros(1,winpt)];
            tmpWave((grouppt-winpt/2):(grouppt+winpt/2-1)) = tmpWave((grouppt-winpt/2):(grouppt+winpt/2-1)).*wintukey';
            tmpWave(1:(grouppt-winpt/2)) = 0;
            tmpWave((grouppt+winpt/2-1):end) = 0;
            NewWinWave = zeros(1, WaveNumPt + HalfFilterNum);
            NewWinWave(1:WaveNumPt) = tmpWave((winpt+1):(winpt+WaveNumPt));
            FilteredWave = fftfilt(filter.Data, NewWinWave(1:(WaveNumPt + HalfFilterNum)));
%            figure; plot(WinWave(1:WaveNumPt),'r'); hold on; plot(NewWinWave); 
        end
              
        FilteredWave = FilteredWave((WaveNumPt + HalfFilterNum):-1:1);
        FilteredWave = fftfilt(filter.Data, FilteredWave(1:(WaveNumPt + HalfFilterNum)));
        FilteredWave = FilteredWave((WaveNumPt + HalfFilterNum):-1:1);
        gfcn.PhaseImg(1:WaveNumPt,numt) = FilteredWave(1:WaveNumPt);
        gfcn.PhaseImg(1:WaveNumPt,numt) = gfcn.PhaseImg(1:WaveNumPt,numt)/max(abs(gfcn.PhaseImg(1:WaveNumPt,numt)));
    end
        
%     hold off
%     time = ((g_WinMinPtNum-1):(g_WinMaxPtNum-1))/gfcn.SampleF; 
%     imagesc(TPoint, time, gfcn.PhaseImg(g_WinMinPtNum:g_WinMaxPtNum,1:cross.NumCtrT),[-1,1]);
%     xlabel('Period (s)', 'FontSize', 8, 'FontWeight', 'bold','FontName','Arial');
%     ylabel('t (s)', 'FontSize', 8, 'FontWeight', 'bold','FontName','Arial');
%     set(gca, 'FontSize', 8, 'FontWeight', 'bold','FontName','Arial');   
%     pause(1);

    hold off
    timeptnum = g_WinMinPtNum:1:g_WinMaxPtNum;
    time = ((g_WinMinPtNum-1):1:(g_WinMaxPtNum-1))/gfcn.SampleF;
    gfcn.PhaseVImg = zeros(VImgPt, cross.NumCtrT);
    for i = 1:cross.NumCtrT
        CenterT = cross.StartT + (i - 1)*cross.DeltaT;
        TravPtV = gfcn.StaDist./(time - CenterT/8); 
 %%%%%%%%%%%%%%%%    guning  %%%%%%%%%%%%%%%%%%%
        [~,aaa] = find(TravPtV == inf);
        if ~isempty(aaa)
            TravPtV(aaa) = 100;
           
        end
      
 %%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        gfcn.PhaseVImg(1:VImgPt, i) = interp1(TravPtV, gfcn.PhaseImg(timeptnum,i), VPoint, 'spline');
%          figure; plot(TravPtV, gfcn.PhaseImg(timeptnum,i)); hold on; plot(VPoint, gfcn.PhaseVImg(1:VImgPt, i), 'r--');
        gfcn.PhaseVImg(1:VImgPt, i) = gfcn.PhaseVImg(1:VImgPt, i)/max(abs(gfcn.PhaseVImg(1:VImgPt, i)));
    end
 
elseif get(handles.FreqDomain, 'Value')==1   % frequency domain phase method
    numfft = fftNumPt;
    halfnfft = numfft/2;
    EGFcnfft = fft(WinWave, numfft);
    fs = gfcn.SampleF/2;
    fdata = (0:halfnfft)/halfnfft*fs;
    df = fdata(2);
    phaseEGF = angle(EGFcnfft);
    nl = floor(1/cross.EndT/df)+1;
    nh = ceil(1/cross.StartT/df)+1;
    nfpt = nh - nl + 1;
    fvImage = zeros(VImgPt, nfpt);
    
    for i = nl:nh
        n = i - nl + 1;
        kr = 2*pi*fdata(i)*gfcn.StaDist./VPoint;

        % exact phase
        % phaseGF = angle(complex(-bessely(0,kr), -besselj(0,kr)));

        % far field phase: correction will be made to compensate the phase
        % difference between far field phase and exact phase later
    %     phaseGF = -(2*pi*fdata(i)*gfcn.StaDist./cgrid + pi/4);
        phaseGF = angle(complex(cos(kr+pi/4), -sin(kr+pi/4)));

        fvImage(:,n) = cos(phaseGF-phaseEGF(i));
    end
    gfcn.PhaseVImg = zeros(VImgPt, cross.NumCtrT);
%     for i = 1:VImgPt
%         gfcn.PhaseVImg(i,:) = interp1(1./fdata(nl:nh), fvImage(i,:), TPoint, 'cubic');
%     end
    [X, Y] = meshgrid(1./fdata(nl:nh), VPoint);
    [XI, YI] = meshgrid(TPoint, VPoint);
    gfcn.PhaseVImg = interp2(X, Y, fvImage, XI, YI, 'pchip');
    for i = 1:cross.NumCtrT
        gfcn.PhaseVImg(:,i) = gfcn.PhaseVImg(:,i)/max(gfcn.PhaseVImg(:,i));
    end
    
elseif get(handles.AkiSPAC, 'Value')==1   % Aki's SPAC method
    
    numfft = fftNumPt;
    halfnfft = numfft/2;
    CFcn = -imag(hilbert(WinWave));
    CFcnfft = fft(CFcn, numfft);
    fs = gfcn.SampleF/2;
    fdata = (0:halfnfft)/halfnfft*fs;
%     df = fdata(2);    
%     LowF = 1/cross.EndT;
%     HighF = 1/cross.StartT;
%     nl = floor(LowF/df)+1;
%     nh = ceil(HighF/df)+1;
    nl = 2;
    nh = halfnfft -1;
    wfftReal = real(CFcnfft);
%     wfftImag = imag(CFcnfft);  % using imagenary part of CFcnfft
    wfftImag = imag(hilbert(real(CFcnfft))); % using the hilbert transform of the real part of CFcnfft
    
    % estimate phase velocity from zero crossing points of J0 and Y0
    mm = find(wfftReal == 0);
    wfftReal(mm) = 1.0e-10;
    mm = find(wfftImag == 0);
    wfftImag(mm) = 1.0e-10;

    AA = wfftReal(nl:nh).*wfftReal(nl+1:nh+1);
    BB = wfftImag(nl:nh).*wfftImag(nl+1:nh+1);
    II = find(AA <= 0);
    JJ = find(BB <= 0);
    II = II + nl - 1;
    JJ = JJ + nl - 1;
    NN = sort([II JJ]);
    nzeros = length(NN);
    zerof = zeros(1,nzeros);
    nII = length(II);
    nJJ = length(JJ);
    for i = 1:nII
        zerof(i) = interp1(wfftReal(II(i):II(i)+1), fdata(II(i):II(i)+1),0);
    end
    for j = 1:nJJ
        zerof(j+nII) = interp1(wfftImag(JJ(j):JJ(j)+1), fdata(JJ(j):JJ(j)+1),0);
    end
    [zerof, index] = sort(zerof); 
    
    fcimage = zeros(VImgPt, nzeros);
    for i = 1:nzeros
        kr = 2*pi*zerof(i)*gfcn.StaDist./VPoint;
        if index(i) <= nII
            fcimage(:,i) = besselj(0, kr).^2;
        else
            fcimage(:,i) = bessely(0, kr).^2;
        end
    end

%     minT = 1./zerof(end);
%     maxT = 1./zerof(1);
%     Tpt = ceil(minT):floor(maxT);
%     nTpt = length(Tpt);
%     Tcimage = zeros(cross.NumCtrT, VImgPt);    
%     for i = 1:nc
%         Tcimage(:,i) = interp1(1./zerof, fcimage(:,i), Tpt, 'cubic');
%     end

    gfcn.PhaseVImg = zeros(VImgPt, cross.NumCtrT);
    [X, Y] = meshgrid(1./zerof, VPoint);
    [XI, YI] = meshgrid(TPoint, VPoint);
    gfcn.PhaseVImg = interp2(X, Y, fcimage, XI, YI, 'pchip');
    for i = 1:cross.NumCtrT
        gfcn.PhaseVImg(:,i) = cos(gfcn.PhaseVImg(:,i)/max(abs(gfcn.PhaseVImg(:,i)))*pi);
    end
    
end

%% --- function for automatic pick of group velocity dispersion curve
function IsDispGood = AutoGroupDisperProcessing(hObject, handles, TPoint, VPoint, GroupVImg, SNRIndex, SNR_T, minSNR)
global g_refgroupdisp gfcn

% GroupVImg(1:VImgPt,1:NumCtrT)

VImgPt = length(VPoint);
NumCtrT = length(TPoint);
dc = VPoint(2) - VPoint(1);
dT = TPoint(2) - TPoint(1);

gRef_low = interp1(g_refgroupdisp(:,1), g_refgroupdisp(:,2), TPoint, 'pchip');
gRef_high = interp1(g_refgroupdisp(:,1), g_refgroupdisp(:,3), TPoint, 'pchip');
gRef = (gRef_low + gRef_high)/2;

% find the approximate maximum period for dispersion analysis which
% satisfies interstation_distance > minlamdaRatio*c*T
minlamdaRatio = str2num(get(handles.minWavelength, 'String'));
lamda = gRef.*TPoint;
II = find((lamda*minlamdaRatio) >= gfcn.StaDist);
if length(II) > 0
    nMaxT = II(1);
else
    nMaxT = NumCtrT;
end

% hold on; plot([1 1]*TPoint(nMaxT), [VPoint(1) VPoint(end)],'g--', 'LineWidth', 2);

gmax_ref = max(gRef_high); % max g in the ref. disper
gmin_ref = min(gRef_low); % min g in the ref. disper

% set tial T and c for searching dispersion curves
T_try_index = round((0:10)/10*(TPoint(nMaxT) - TPoint(1))/dT)+1;  % trial T index for searching dispersion curve
T_try = TPoint(1) + (T_try_index - 1)*dT;  % trial T for searching dispersion index
[maxamp, g_try_index] = max(GroupVImg(:,T_try_index),[],1);
g_try = VPoint(1) + (g_try_index - 1)*dc;

% search disperison curves
GroupVDisp_try = zeros(length(T_try), NumCtrT);
k = 0;
for i = 1:length(T_try)
%     hold on; plot(T_try(i), g_try(i), 'w*');
    Initialg = g_try_index(i);
    InitialT = T_try_index(i);
    DispPt = AutoSearch(Initialg, InitialT, GroupVImg);

    if k == 0
        k =  k + 1;
        GroupVDisp_try(k,:) = VPoint(1) + (DispPt - 1)*dc;
%         hold on; plot(TPoint, GroupVDisp_try(k,:), 'g');
    else
        tempDisp = VPoint(1) + (DispPt - 1)*dc;
        SameDisperIndex = 0;
        for nn = 1:k
            if sum(abs(GroupVDisp_try(nn,:) - tempDisp)) < 1.0e-4
                SameDisperIndex = 1;
            end
        end
        if SameDisperIndex == 0
            k = k + 1;
            GroupVDisp_try(k,:) = tempDisp;
%             hold on; plot(TPoint, GroupVDisp_try(k,:), 'g');
        end
    end
end

% determine the quality of the dispersion curve by looking at how many
% disperion points fall in the reference dispersion range
NumDispCurve = k;
InRangePt = zeros(1,NumDispCurve); % 
for i = 1:NumDispCurve
    GoodIndex = sign(GroupVDisp_try(i,1:nMaxT) - gRef_low(1:nMaxT)) + sign(gRef_high(1:nMaxT) - GroupVDisp_try(i,1:nMaxT));
    II = find(GoodIndex == 2);
    InRangePt(i) = length(II);  
end
maxpt = max(InRangePt);
meanpt = mean(InRangePt);
II = find(InRangePt >= (2*maxpt+meanpt)/3); % find the best several dispersion curves within reference range
if length(II) == 1
    GroupVDisp = GroupVDisp_try(II,:);
else
    RefObsDispDiff = zeros(1, length(II));
    ObsSumAbsDiff = zeros(1, length(II));
    for i = 1:length(II)
        RefObsDispDiff(i) = sum(abs(GroupVDisp_try(II(i),1:nMaxT) - (gRef_low(1:nMaxT) + gRef_high(1:nMaxT))/2));
        ObsSumAbsDiff(i) = sum(abs(diff(GroupVDisp_try(II(i),1:nMaxT))));
    end
    [mindiff, index1] = min(RefObsDispDiff); % find lowest difference dispersion curve with respect to reference
    [minabs, index2] = min(ObsSumAbsDiff); % find smoothest dispersion curve
    if index1 == index2
        GroupVDisp = GroupVDisp_try(II(index1),:);
    else
        BestTwoDiff = abs(GroupVDisp_try(II(index1),1:nMaxT) - GroupVDisp_try(II(index2),1:nMaxT));
        if length(find(BestTwoDiff < 1.0e-3)) > 0.67*nMaxT  % 2/3 of two best dispersion curves are overlapping
            GroupVDisp = GroupVDisp_try(II(index2),:); % choose the smoothest one
        else
            GroupVDisp = GroupVDisp_try(II(index1),:); % choose the smaller difference one if 2/3 dispersion are different
        end
    end    
end
hold on; plot(TPoint, GroupVDisp,'k', 'LineWidth', 2); % plot the final selected dispersion curve

RawDisper = [TPoint' GroupVDisp'];

% MinTWidth = dT*3;   % set minimum length (sec) of dispersion data to be saved

NewDisper = [RawDisper ones(NumCtrT,1)];

% whether NewDisper falls in Group V range or not
for ii = 1:size(NewDisper,1)
    if NewDisper(ii,2) > gRef_high(ii) || NewDisper(ii,2) < gRef_low(ii)
        NewDisper(ii,3)=0;
    end
end

% find the group velocity corresponding to the maximum amplitude of the
% envelope at each period
[GroupVMaxAmp, maxIndex] = max(GroupVImg,[],1);
JJ = abs(GroupVMaxAmp - GroupVDisp) < 0.01;
% if the picked dispersion is very different from the dispersion
% corresponding to the maximum amplitude of the envelope (only 1/10 points
% are overlapping), not save the dispersion curve and will not pick phase
% velocity dispersion curve by setting IsDispGood = 0
if length(JJ) < 0.1*NumCtrT  
    IsDispGood = 0;
else   
    % find reasonable dispersion points with high SNR
    GoodIndex = NewDisper(:,3) + SNRIndex';
    II = find((NewDisper(:,3) + SNRIndex') == 2);
    hold on; plot(NewDisper(:,1), NewDisper(:,2), 'co','MarkerSize',8,'MarkerFaceColor','c');
    hold on; plot(NewDisper(II,1), NewDisper(II,2), 'ko','MarkerSize',6,'MarkerFaceColor','k');

    % save dispersion data
    SaveAutoPickDisper(hObject, handles, GroupVDisp, GoodIndex);

    IsDispGood = 2;
end


%% --- function for automatic pick of phase velocity dispersion curve
function IsDispGood = AutoPhaseDisperProcessing(hObject, handles, TPoint, VPoint, PhaseVImg, SNRIndex, SNR_T, minSNR)
global g_refphasedisp gfcn

% PhaseVImg(1:VImgPt,1:NumCtrT)
VImgPt = length(VPoint);  % samples of V, size = 176
NumCtrT = length(TPoint);  % samples of T, size = 49
dc = VPoint(2) - VPoint(1);
dT = TPoint(2) - TPoint(1);
cRef_low = interp1(g_refphasedisp(:,1), g_refphasedisp(:,2), TPoint, 'pchip');
cRef_high = interp1(g_refphasedisp(:,1), g_refphasedisp(:,3), TPoint, 'pchip');
cRef = (cRef_low + cRef_high)/2;

% find the approximate maximum period for dispersion analysis which
% satisfies interstation_distance > minlamdaRatio*c*T
minlamdaRatio = str2num(get(handles.minWavelength, 'String'));
lamda = cRef.*TPoint;
II = find((lamda*minlamdaRatio) >= gfcn.StaDist);
if length(II) > 0
    nMaxT = II(1);
else
    nMaxT = NumCtrT;
end

% hold on; plot([1 1]*TPoint(nMaxT), [VPoint(1) VPoint(end)],'g--', 'LineWidth', 2);

cmax_ref = max(cRef_high); % max c in the ref. disper
cmin_ref = min(cRef_low); % min c in the ref. disper

% set tial T and c for searching dispersion curves
T_try_index = round([1 5]/10*(TPoint(nMaxT) - TPoint(1))/dT)+1;  % trial T index for searching dispersion curve [6 25]
T_try = TPoint(1) + (T_try_index - 1)*dT;  % trial T for searching dispersion index [0.7 2.6]
c_try_index = round(((cmin_ref - VPoint(1)) + (cmax_ref - cmin_ref)*(1:9)/10)/dc) + 1;
c_try = VPoint(1) + (c_try_index - 1)*dc;

% search disperison curves
PhaseVDisp_try = zeros(length(T_try)*length(c_try), NumCtrT);
k = 0;
for i = 1:length(T_try)
%     hold on; plot(T_try(i)*ones(1, length(c_try)), c_try, 'w*');
    for j = 1:length(c_try)
        Initialc = c_try_index(j);
        InitialT = T_try_index(i);
        DispPt = AutoSearch(Initialc, InitialT, PhaseVImg);
        
        % ysb: delete the repeated disper curves.
        if k == 0
            k =  k + 1;
            PhaseVDisp_try(k,:) = VPoint(1) + (DispPt - 1)*dc;
            % hold on; plot(TPoint, PhaseVDisp_try(k,:), 'g');
        else
            tempDisp = VPoint(1) + (DispPt - 1)*dc;
            SameDisperIndex = 0;
            for nn = 1:k
                if sum(abs(PhaseVDisp_try(nn,:) - tempDisp)) < 1.0e-4
                    SameDisperIndex = 1;
                end
            end
            if SameDisperIndex == 0
                k = k + 1;
                PhaseVDisp_try(k,:) = tempDisp;
%                 hold on; plot(TPoint, PhaseVDisp_try(k,:), 'g');
            end
        end
                    
    end
end

% determine the quality of the dispersion curve by looking at how many
% disperion points fall in the reference dispersion range
NumDispCurve = k;
InRangePt = zeros(1,NumDispCurve); % ysb: number of disperion points fall in the reference dispersion range
for i = 1:NumDispCurve
    GoodIndex = sign(PhaseVDisp_try(i,1:nMaxT) - cRef_low(1:nMaxT)) + sign(cRef_high(1:nMaxT) - PhaseVDisp_try(i,1:nMaxT));
    II = find(GoodIndex == 2);
    InRangePt(i) = length(II);  
end
maxpt = max(InRangePt);
meanpt = mean(InRangePt);
II = find(InRangePt >= (2*maxpt+meanpt)/3); % find the best several dispersion curves within reference range
if length(II) == 1
    PhaseVDisp = PhaseVDisp_try(II,:);
    hold on; plot(TPoint, PhaseVDisp,'k', 'LineWidth', 2);
else
    RefObsDispDiff = zeros(1, length(II));
    ObsSumAbsDiff = zeros(1, length(II));
    for i = 1:length(II)
        RefObsDispDiff(i) = sum(abs(PhaseVDisp_try(II(i),1:nMaxT) - (cRef_low(1:nMaxT) + cRef_high(1:nMaxT))/2));
        ObsSumAbsDiff(i) = sum(abs(diff(PhaseVDisp_try(II(i),1:nMaxT)))); % ysb: the smoothness of the curve
    end
    [mindiff, index1] = min(RefObsDispDiff); % find lowest difference dispersion curve with respect to reference
    [minabs, index2] = min(ObsSumAbsDiff); % find smoothest dispersion curve
    
    if index1 == index2
        PhaseVDisp = PhaseVDisp_try(II(index1),:);
    else
        BestTwoDiff = abs(PhaseVDisp_try(II(index1),1:nMaxT) - PhaseVDisp_try(II(index2),1:nMaxT));
        if length(find(BestTwoDiff < 1.0e-3)) > 0.67*nMaxT  % 2/3 of two best dispersion curves are overlapping
            PhaseVDisp = PhaseVDisp_try(II(index2),:); % choose the smoothest one
        else
            PhaseVDisp = PhaseVDisp_try(II(index1),:); % choose the smaller difference one if 2/3 dispersion are different
        end
    end        
end

hold on; plot(TPoint, PhaseVDisp,'k', 'LineWidth', 2); % plot the final selected dispersion curve

RawDisper = [TPoint' PhaseVDisp'];
MinTWidth = dT*3;   % set minimum length (sec) of dispersion data to be saved
LoveRayIndex = get(handles.RayleighWave, 'Value'); % 1 for Rayleigh, 0 for Love

if LoveRayIndex == 0
    NewDisper = AutomaticDisperLove(RawDisper,1/dT,MinTWidth, SNR_T, minSNR, TPoint(nMaxT));
else
    NewDisper = AutomaticDisperRayl(RawDisper,1/dT,MinTWidth, SNR_T, minSNR, TPoint(nMaxT));
end

% whether NewDisper falls in PhaseV range or not
for ii = 1:size(NewDisper,1)
    if NewDisper(ii,2) > cRef_high(ii) || NewDisper(ii,2) < cRef_low(ii)
        NewDisper(ii,2)=0;
        NewDisper(ii,3)=0;
    end
end

% NewDisper = [RawDisper, ones(size(RawDisper,1),1)];

% find reasonable dispersion points
% II = find(NewDisper(:,3) > 0);
% hold on; plot(NewDisper(II,1), NewDisper(II,2), 'k*');
% find reasonable dispersion points with high SNR
GoodIndex = NewDisper(:,3) + SNRIndex';
II = find((NewDisper(:,3) + SNRIndex') == 2);

% save dispersion data when having at least 4 or 0.1*NumCtrT good points
if length(II) >= 4 || length(II) > 0.1*NumCtrT
    hold on; plot(NewDisper(:,1), NewDisper(:,2), 'co','MarkerSize',8,'MarkerFaceColor','c');
    hold on; plot(NewDisper(II,1), NewDisper(II,2), 'ko','MarkerSize',6,'MarkerFaceColor','k');
    % save dispersion data
    SaveAutoPickDisper(hObject, handles, PhaseVDisp, GoodIndex);
    IsDispGood = 2;
else
    IsDispGood = 0;
end




% save and plot the reasonable dispersion points
function SaveAutoPickDisper(hObject, handles, VDisp, GoodIndex)

global cross gfcn g_PrevStaLonLat g_FileName 
global g_DispDirectory
global g_DisperType

h1 = handles.axes1;
h2 = handles.axes2;

StaDistance = gfcn.StaDist;
TPoint = cross.StartT:cross.DeltaT:cross.EndT;
wavelength = TPoint.*VDisp(1:cross.NumCtrT);
I1 = find(wavelength > gfcn.StaDist);
I2 = find(wavelength > gfcn.StaDist/2);
if length(I1) > 0
    hold on; plot([TPoint(I1(1)) TPoint(I1(1))], [gfcn.WinMinV gfcn.WinMaxV], 'r--', 'LineWidth', 2);
end
if length(I2) > 0
    hold on; plot([TPoint(I2(1)) TPoint(I2(1))], [gfcn.WinMinV gfcn.WinMaxV], 'y--', 'LineWidth', 2);
end  

minlamdaRatio = str2num(get(handles.minWavelength, 'String'));
I3 = find((wavelength*minlamdaRatio) > gfcn.StaDist);
if length(I3) > 0
    hold on; plot([TPoint(I3(1)) TPoint(I3(1))], [gfcn.WinMinV gfcn.WinMaxV], 'g--', 'LineWidth', 2);
end  

if g_DisperType == 1
    if get(handles.TimeDomain, 'Value')==1
        dispfilename = strcat('CDisp.T.', g_FileName);
        set(handles.MsgEdit, 'String', strcat('Output Disper File Name: ', dispfilename));
    elseif get(handles.FreqDomain, 'Value')==1
        dispfilename = strcat('CDisp.F.', g_FileName);
        set(handles.MsgEdit, 'String', strcat('Output Disper File Name: ', dispfilename));
    elseif get(handles.AkiSPAC, 'Value')==1 
        dispfilename = strcat('CDisp.A.', g_FileName);
        set(handles.MsgEdit, 'String', strcat('Output Disper File Name: ', dispfilename));
    end
elseif g_DisperType == 2
    dispfilename = strcat('GDisp.', g_FileName);
    set(handles.MsgEdit, 'String', strcat('Output Disper File Name: ', dispfilename));
end

%write T-V to file 
if g_DispDirectory == 0
    g_DispDirectory = pwd;
end

if sum(sum(g_PrevStaLonLat)) ~= 0
    set(gcf,'CurrentAxes',h1);
    hold(h1,'on');
    plot(g_PrevStaLonLat(1:2,1), g_PrevStaLonLat(1:2,2), 'b-');
%         plot(h1, g_PrevStaLonLat(1:2,1), g_PrevStaLonLat(1:2,2), 'k^', 'MarkerSize',6, 'MarkerFaceColor','k');
end


DataExistIndex = 0;
minlamdaRatio = str2num(get(handles.minWavelength, 'String'));
for i = 1:cross.NumCtrT
   wavelength = VDisp(i)*TPoint(i);
   if (StaDistance >= minlamdaRatio*wavelength) && GoodIndex(i) == 2
        DataExistIndex = 1;
        
        if g_DisperType == 1
            % make correction of phase velocity at near field for time domain dispersion measurements
            if get(handles.TimeDomain, 'Value')==1 || get(handles.FreqDomain, 'Value')==1  % time/freq domain method   
                kr = 2*pi*StaDistance/wavelength; 
                G2d = complex(-bessely(0,kr),-besselj(0,kr));  % exact Green's function
                G2d_ff = complex(cos(kr+pi/4), -sin(kr+pi/4));  % Green's function at far field
                phiG2d = -angle(G2d)*57.3; 
                phiG2d_ff = -angle(G2d_ff)*57.3;
                phiDiff = phiG2d_ff - phiG2d; % phase difference in degree for G2d
                if phiDiff >= 180
                    phiDiff = phiDiff - 360;
                elseif phiDiff < -180
                    phiDiff = phiDiff + 360;
                end
                t_corrected = (StaDistance/VDisp(i) + TPoint(i)/8) - (phiDiff/360)*TPoint(i);
                VDisp(i) = StaDistance/(t_corrected - TPoint(i)/8);           
            end
        end

   else
        GoodIndex(i) = 0;
   end
end

II  = find(GoodIndex == 2);
JJ = find(GoodIndex < 2);
GoodIndex(II) = 1;
GoodIndex(JJ) = 0;

if DataExistIndex == 1
    TVfile = strcat(g_DispDirectory, '/', dispfilename);    
    ftv = fopen(TVfile,'w');
    fprintf(ftv,'%f     ', gfcn.Lon1);
    fprintf(ftv,'%f\n', gfcn.Lat1);
    fprintf(ftv,'%f     ', gfcn.Lon2);
    fprintf(ftv,'%f\n', gfcn.Lat2);
    for i = 1:cross.NumCtrT
        fprintf(ftv,'%4.3f    %4.3f    %4.3f   %4d\n',[TPoint(i) VDisp(i)*GoodIndex(i) 0.0  GoodIndex(i)]);
    end
    fclose(ftv); 
    set(gcf,'CurrentAxes',h2);
    hold(h2,'on');
    plot(TPoint(II), VDisp(II), 'ro', 'MarkerSize', 4, 'MarkerFaceColor','r');
  
    
    set(gcf,'CurrentAxes',h1);
    hold(h1,'on');
    plot([gfcn.Lon1, gfcn.Lon2], [gfcn.Lat1, gfcn.Lat2], 'r-');
%    plot(h1, [gfcn.Lon1, gfcn.Lon2], [gfcn.Lat1, gfcn.Lat2], 'k^', 'MarkerSize',6, 'MarkerFaceColor','g');
    g_PrevStaLonLat = [gfcn.Lon1, gfcn.Lat1;gfcn.Lon2, gfcn.Lat2];
end

if g_DisperType == 1
    gfcn.PhaseVDisp = VDisp;
elseif g_DisperType == 2
    gfcn.GroupVDisp = VDisp;
end

%% --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NewDisper = AutomaticDisperLove(RawDisper,SampleF,MinTWidth, SNR, minSNR, MaxTCal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fove Love waves
% This function extracts the reasonable part of a whole disper curve
% between two stations.

% MinTWidth=5;  % choose segments containing at least 5s data

% Rewrite RawDisper
if size(SNR,1) ~= size(RawDisper,1)
    SNR = SNR';
end    
RawDisper = [RawDisper, SNR];
% RawDisper
% 1     2     3
% T   PhaseV  SNR

% Range of dcdt: if needed, this part can be put in input parameters.
% 1       2          3     
% T   upper-limit  lower-limit
r_dcdt=[  0  0.30  -0.025
       10  0.10  -0.025
       30  0.06  -0.010
       45  0.03  -0.006
      200  0.02  -0.003];  
    
r_dcdt=r_dcdt/SampleF;

% dcdt: time derivative of phase velocity
% 1     2     3
% T    dcdt  g/b 
%
% if lower-limit <= dcdt < upper-limit : good = 1
% if 2 * lower-limit < dcdt < 1.5 * upper-limit  : half-good =1/2
% other cases : bad = 0

% Create smooth range of dcdt
Smooth_r_dcdt = zeros(size(RawDisper,1)-1,3);
Smooth_r_dcdt(:,1) = RawDisper(1:end-1,1);
Smooth_r_dcdt(:,2) = interp1(r_dcdt(:,1),r_dcdt(:,2),Smooth_r_dcdt(:,1),'pchip');
Smooth_r_dcdt(:,3) = interp1(r_dcdt(:,1),r_dcdt(:,3),Smooth_r_dcdt(:,1),'pchip');

% dcdt: time derivative of phase velocity
% 1     2     3
% T    dcdt  g/b 
%
% if lower-limit <= dcdt < upper-limit : good = 1
% if 2 * lower-limit < dcdt < 1.5 * upper-limit  : half-good =1/2
% other cases : bad = 0

dcdt=zeros(size(RawDisper,1)-1,3);
dcdt(:,1)=RawDisper(1:end-1,1);
dcdt(:,2)=RawDisper(2:end,2)-RawDisper(1:end-1,2);

for ii=1:size(dcdt,1)
    upper_limit=Smooth_r_dcdt(ii,2);
    lower_limit=Smooth_r_dcdt(ii,3);
    
    if dcdt(ii,2)>0 && dcdt(ii,2)<=upper_limit
        dcdt(ii,3)=1;
    elseif dcdt(ii,2)>= 1.5 * lower_limit && dcdt(ii,2)<=upper_limit*1.5
        dcdt(ii,3)=0.5;
    end
end

% Recheck data quality
% if d(dc/dt)/dt changes fast => bad 
for ii = 2:size(dcdt,1)
    ddcdtt = abs(dcdt(ii,2)-dcdt(ii-1,2));
    if ddcdtt > Smooth_r_dcdt(ii,2)-Smooth_r_dcdt(ii,3)
        dcdt(ii,3)=0;
    end
end
    
% if there is less than 2 half-good pts in the middle of good pts, just
% accepted them
if SampleF == 1
    for ii=2:size(dcdt,1)-2
        if dcdt(ii,3)==0.5 && dcdt(ii-1,3)==1 && dcdt(ii+1,3)==1
            dcdt(ii,3)=1;
        elseif dcdt(ii,3)==0.5 && dcdt(ii-1,3)==1 && dcdt(ii+1,3)>0 && dcdt(ii+2,3)==1
            dcdt(ii,3)=1;
        end
    end
elseif SampleF == 2
    for ii = 2:size(dcdt,1)-SampleF
        if dcdt(ii,3)==0 && dcdt(ii-1,3)>0 && dcdt(ii+1,3)>0
            dcdt(ii,3)=0.5;
        end
    end
    for ii = 2:size(dcdt,1)-SampleF
        if dcdt(ii,3) == 0.5 && dcdt(ii-1,3)==1
            QualityAfter=dcdt(ii+1:end,3);
            NearestOne=min(QaulityAfter==1);
            if sum(QaulityAfter(1:NearestOne))/NearestOne > 0.5
                dcdt(ii,3)=1;
            end
        end
    end
end

% Calculate the length of each good-data segment: dcdt(ii,3)==1
% Nonzero:
%   1      2       3 
% Tstart  Tend   length
Nonzero=zeros(size(dcdt,1),3);
NonzeroIndex=0;
for ii=1:size(dcdt,1)
    if dcdt(ii,3)==1
        if ii==1
            NonzeroIndex=NonzeroIndex+1;
            Nonzero(NonzeroIndex,1)=dcdt(ii,1);
        elseif dcdt(ii-1,3)<1
            NonzeroIndex=NonzeroIndex+1;
            Nonzero(NonzeroIndex,1)=dcdt(ii,1);
        end
        Nonzero(NonzeroIndex,3)=Nonzero(NonzeroIndex,3)+1;
    end
end
Nonzero(:,2)=Nonzero(:,1)+Nonzero(:,3)/SampleF;

MinSegLength=MinTWidth*SampleF;  % choose segments containing at least MinTWidth data 
AcceptedSeg=find(Nonzero(:,3)>=MinSegLength);
SegNum=length(AcceptedSeg);

if SegNum<1
    NewDisper=zeros(size(RawDisper,1),3);
    NewDisper(:,1) = RawDisper(:,1);
    return
end

TstartIX=find( RawDisper(:,1)==Nonzero(AcceptedSeg(1),1) );
TendIX=find( RawDisper(:,1)==Nonzero(AcceptedSeg(1),2) );

NewDisper=RawDisper(TstartIX:TendIX,:);

for ii=2:SegNum
    if Nonzero(AcceptedSeg(ii),1) > MaxTCal - MinTWidth
        break;
    end
    PrevTstartIX = TstartIX;
    PrevTendIX = TendIX;    
    TstartIX = find( RawDisper(:,1)==Nonzero(AcceptedSeg(ii),1) );
    TendIX = find( RawDisper(:,1)==Nonzero(AcceptedSeg(ii),2) );
    
    % whether this segment and its previous one both belong to the same
    % dispersion curve
    % The dispersion curve should increases in a reasonable range
    PhaseVGap = RawDisper(TstartIX,2) - RawDisper(PrevTendIX,2);
    TGap = RawDisper(TstartIX,1) - RawDisper(PrevTendIX,1);
    upper_PhaseVGap1 = mean(dcdt(max(PrevTendIX-4,PrevTstartIX):PrevTendIX-1,2))*TGap*2;
    upper_PhaseVGap2 = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(PrevTendIX,1),RawDisper(TstartIX,1)); 
    upper_PhaseVGap = max(upper_PhaseVGap1,upper_PhaseVGap2);
    lower_PhaseVGap = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(PrevTendIX,1),RawDisper(TstartIX,1))/5+min(r_dcdt(:,3))/2*TGap;
    
    
    if PhaseVGap > upper_PhaseVGap || PhaseVGap < lower_PhaseVGap
        meanSNR1 = mean(NewDisper(:,3));
        meanSNR2 = mean(RawDisper(TstartIX:TendIX,3));
        IIII = find(NewDisper(:,3) >= minSNR);
        JJJJ = find(RawDisper(TstartIX:TendIX,3) >= minSNR);
        
        if length(IIII) > 1.5*length(JJJJ)
            continue
        end
        
%         if  meanSNR1 < 2*minSNR || meanSNR2 < 2*minSNR
            if meanSNR1 >= meanSNR2 && NewDisper(end,1) > 10
                break;
            else
                NewDisper=[NaN,NaN,NaN];
            end
%         end
    end
    
    NewDisper=[NewDisper;RawDisper(TstartIX:TendIX,:)];
end
    
% Rewrite NewDisper
FullNewDisper = zeros(size(RawDisper,1),3);
FullNewDisper(:,1) = RawDisper(:,1);
for ii=1:size(FullNewDisper,1)
    T = FullNewDisper(ii,1);
    TIX = find( NewDisper(:,1) == T );
    if length(TIX) == 1
        FullNewDisper(ii,2) = NewDisper(TIX,2);
        FullNewDisper(ii,3) = 1;
    end
end
clear NewDisper
NewDisper = FullNewDisper;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NewDisper = AutomaticDisperRayl(RawDisper,SampleF,MinTWidth, SNR, minSNR, MaxTCal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for rayleigh wave 
% This function extracts the reasonable part of a whole disper curve
% between two stations.

% MinTWidth=5;  % choose segments containing at least 5s data

% Rewrite RawDisper
if size(SNR,1) ~= size(RawDisper,1)
    SNR = SNR';
end    
RawDisper = [RawDisper, SNR];
% RawDisper
% 1     2     3
% T   PhaseV  SNR

% Range of dcdt: if needed, this part can be put in input parameters.
% 1       2          3     
% T   upper-limit  lower-limit
r_dcdt=[  0  0.20  -0.20
         10  0.30  -0.15
         20  0.15  -0.10
         30  0.10  -0.01
         40  0.06  -0.005
         50  0.04  -0.003
        200  0.02  -0.003];   

% r_dcdt=r_dcdt/SampleF;
% ysb:
 r_dcdt=r_dcdt/2;

% Create smooth range of dcdt
% ysb: format: 
%      1: T   2: lower-limit    3: up-limit
Smooth_r_dcdt = zeros(size(RawDisper,1)-1,3);
Smooth_r_dcdt(:,1) = RawDisper(1:end-1,1);
Smooth_r_dcdt(:,2) = interp1(r_dcdt(:,1),r_dcdt(:,2),Smooth_r_dcdt(:,1),'pchip');
Smooth_r_dcdt(:,3) = interp1(r_dcdt(:,1),r_dcdt(:,3),Smooth_r_dcdt(:,1),'pchip');

% dcdt: time derivative of phase velocity
% 1     2     3
% T    dcdt  g/b 
%
% if lower-limit <= dcdt < upper-limit : good = 1
% if 2 * lower-limit < dcdt < 1.5 * upper-limit  : half-good =1/2
% other cases : bad = 0

dcdt=zeros(size(RawDisper,1)-1,3);
dcdt(:,1)=RawDisper(1:end-1,1);
dcdt(:,2)=RawDisper(2:end,2)-RawDisper(1:end-1,2);

for ii=1:size(dcdt,1)
    upper_limit=Smooth_r_dcdt(ii,2);
    lower_limit=Smooth_r_dcdt(ii,3);
    
    if dcdt(ii,2)>lower_limit && dcdt(ii,2)<=upper_limit
        if dcdt(ii,1) >= 30  % for period less than 30 s, phase velocity of Rayleigh wave usually increases as T increases
            if dcdt(ii,2)<=0
                dcdt(ii,3)=0.5;  % ysb: the weight of decrease V is lower
            else
                dcdt(ii,3)=1;
            end
        else  % for period < 30s
            dcdt(ii,3)=1;
        end     
    elseif dcdt(ii,2)>= 1.5*lower_limit && dcdt(ii,2)<=upper_limit*1.5
        dcdt(ii,3)=0.5;
    end
end

% Recheck data quality
% if d(dc/dt)/dt changes fast => bad 
for ii = 2:size(dcdt,1)
    ddcdtt = abs(dcdt(ii,2)-dcdt(ii-1,2));
    if ddcdtt > Smooth_r_dcdt(ii,2)-Smooth_r_dcdt(ii,3)
        dcdt(ii,3)=0;
    end
end

% if there is less than 2 half-good pts in the middle of good pts, just
% accepted them
if SampleF == 1
    for ii=2:size(dcdt,1)-2
        if dcdt(ii,3)==0.5 && dcdt(ii-1,3)==1 && dcdt(ii+1,3)==1
            dcdt(ii,3)=1;
        elseif dcdt(ii,3)==0.5 && dcdt(ii-1,3)==1 && dcdt(ii+1,3)>0 && dcdt(ii+2,3)==1
            dcdt(ii,3)=1;
        end
    end
elseif SampleF == 2
    for ii = 2:size(dcdt,1)-SampleF
        if dcdt(ii,3)==0 && dcdt(ii-1,3)>0 && dcdt(ii+1,3)>0
            dcdt(ii,3)=0.5;
        end
    end
    for ii = 2:size(dcdt,1)-SampleF
        if dcdt(ii,3) == 0.5 && dcdt(ii-1,3)==1
            QaulityAfter=dcdt(ii+1:end,3);
            NearestOne=min(QaulityAfter==1);
            if sum(QaulityAfter(1:NearestOne))/NearestOne > 0.5
                dcdt(ii,3)=1;
            end
        end
    end
end

% Calculate the length of each good-data segment: dcdt(ii,3)==1
% Nonzero:
%   1      2       3 
% Tstart  Tend   length
Nonzero=zeros(size(dcdt,1),3);
NonzeroIndex=0;
for ii=1:size(dcdt,1)
    if dcdt(ii,3)==1
        if ii==1
            NonzeroIndex=NonzeroIndex+1;
            Nonzero(NonzeroIndex,1)=dcdt(ii,1);
        elseif dcdt(ii-1,3)<1
            NonzeroIndex=NonzeroIndex+1;
            Nonzero(NonzeroIndex,1)=dcdt(ii,1);
        end
        Nonzero(NonzeroIndex,3)=Nonzero(NonzeroIndex,3)+1;
    end
end
Nonzero(:,2)=Nonzero(:,1)+Nonzero(:,3)/SampleF;

MinSegLength=MinTWidth*SampleF;  % choose segments containing at least MinTWidth data 
% AcceptedSeg=find(Nonzero(:,3)>=MinSegLength);
AcceptedSeg=find(Nonzero(:,3)>=MinSegLength);
SegNum=length(AcceptedSeg);

if SegNum<1
    NewDisper=zeros(size(RawDisper,1),3);
    NewDisper(:,1) = RawDisper(:,1);
    return
end

% ysb: seletc the first(short period) segment?
% ysb: find return a empty matrix because of the precision of system
% TstartIX=find( abs(RawDisper(:,1)-Nonzero(AcceptedSeg(1),1))<1e-5 );
% TendIX=find( abs(RawDisper(:,1)-Nonzero(AcceptedSeg(1),2))<1e-5 );  
max_len_index = find(Nonzero(:,3)==max(Nonzero(:,3)));
if length(max_len_index) ~= 1
    max_len_index = max_len_index(1);
end

TstartIX=find( abs(RawDisper(:,1)-Nonzero(max_len_index,1))<1e-5 );
TendIX=find( abs(RawDisper(:,1)-Nonzero(max_len_index,2))<1e-5 );  

NewDisper=RawDisper(TstartIX:TendIX,:);

% % ysb:
% for ii=2:SegNum
%     if Nonzero(AcceptedSeg(ii),1) > MaxTCal - MinTWidth
%         break;
%     end
%     PrevTstartIX = TstartIX;
%     PrevTendIX = TendIX;    
%     TstartIX = find( abs(RawDisper(:,1)-Nonzero(AcceptedSeg(ii),1))<1e-5 );
%     TendIX = find( abs(RawDisper(:,1)-Nonzero(AcceptedSeg(ii),2))<1e-5 );
%     
%     % whether this segment and its previous one both belong to the same
%     % dispersion curve
%     % The dispersion curve should increases in a reasonable range
%     PhaseVGap = RawDisper(TstartIX,2) - RawDisper(PrevTendIX,2);
%     TGap = RawDisper(TstartIX,1) - RawDisper(PrevTendIX,1);
%     upper_PhaseVGap1 = mean(dcdt(max(PrevTendIX-4,PrevTstartIX):PrevTendIX-1,2))*TGap*2;
%     upper_PhaseVGap2 = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(PrevTendIX,1),RawDisper(TstartIX,1)); 
%     upper_PhaseVGap = max(upper_PhaseVGap1,upper_PhaseVGap2);
%     lower_PhaseVGap = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(PrevTendIX,1),RawDisper(TstartIX,1))/5;
%     
%     % PhaseVGap
%     % upper_PhaseVGap
%     if PhaseVGap > upper_PhaseVGap || PhaseVGap < lower_PhaseVGap
%         meanSNR1 = mean(NewDisper(:,3));
%         meanSNR2 = mean(RawDisper(TstartIX:TendIX,3));
%         IIII = find(NewDisper(:,3) >= minSNR);
%         JJJJ = find(RawDisper(TstartIX:TendIX,3) >= minSNR);
%         if length(IIII) > 1.5*length(JJJJ)
%             continue
%         end
%             
%          % if  meanSNR1 < 2*minSNR || meanSNR2 < 2*minSNR
%             if meanSNR1 >= meanSNR2 && NewDisper(end,1) > 10
%                 break;
%             else
%                 NewDisper=[NaN,NaN,NaN];
%             end
%          % end
%     end
% end

% Rewrite NewDisper
FullNewDisper = zeros(size(RawDisper,1),3);
FullNewDisper(:,1) = RawDisper(:,1);
for ii=1:size(FullNewDisper,1)
    T = FullNewDisper(ii,1);
    TIX = find( NewDisper(:,1) == T );
    if length(TIX) == 1
        FullNewDisper(ii,2) = NewDisper(TIX,2);
        FullNewDisper(ii,3) = 1;
    end
end
clear NewDisper
NewDisper = FullNewDisper;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r_PhaseVGap=PhaseVGapIntegral(Smooth_r_dcdt,Tstart,Tend)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_PhaseVGap = 0;
%%% gunign %%%%%%%%
if ~isempty(Tstart)
    TstartIX = find(Smooth_r_dcdt(:,1)==Tstart);
    TendIX = find(Smooth_r_dcdt(:,1)==Tend);
    r_PhaseVGap = sum(Smooth_r_dcdt(TstartIX:TendIX-1,2));
else
   TendIX = find(Smooth_r_dcdt(:,1)==Tend);
   r_PhaseVGap = sum(Smooth_r_dcdt(1:TendIX-1,2)); 
end

% This function calculte the integral of a known function
% WhichBand1=min(find(r_dcdt(:,1)>=Tstart));
% WhichBand2=min(find(r_dcdt(:,1)>=Tend));
% SegNum=WhichBand2-WhichBand1+1;
% 
% if SegNum==1
%     r_PhaseVGap=r_dcdt(WhichBand1,2)*(Tend-Tstart);
% else
%     r_PhaseVGap=r_dcdt(WhichBand1,2)*(r_dcdt(WhichBand1,1)-Tstart);
%     for ii=2:SegNum-1
%         r_PhaseVGap=r_PhaseVGap+r_dcdt(WhichBand1+ii-1,2)*(r_dcdt(WhichBand1+ii-1,1)-r_dcdt(WhichBand1+ii-2,1));
%     end
%     r_PhaseVGap=r_PhaseVGap+r_dcdt(WhichBand2,2)*(Tend-r_dcdt(WhichBand2-1,1));
% end


function FilterBandWidthT_Callback(hObject, eventdata, handles)
% hObject    handle to FilterBandWidthT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FilterBandWidthT as text
%        str2double(get(hObject,'String')) returns contents of FilterBandWidthT as a double


% --- Executes during object creation, after setting all properties.
function FilterBandWidthT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilterBandWidthT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FilterCenterT_Callback(hObject, eventdata, handles)
% hObject    handle to FilterCenterT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FilterCenterT as text
%        str2double(get(hObject,'String')) returns contents of FilterCenterT as a double


% --- Executes during object creation, after setting all properties.
function FilterCenterT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilterCenterT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%
% --- Executes on button press in PhaseVFilterDesign.
function PhaseVFilterDesign_Callback(hObject, eventdata, handles)
% hObject    handle to PhaseVFilterDesign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filter

%     errordlg('No data sampling frequency information! Please check it!');
prompt = {'Please enter your data sampling frequency (Hz):'};
title = ['Data Sampling Frequency'];
line = 2;
DataSPF = inputdlg(prompt,title,line);
filter.SampleF = str2num(DataSPF{1});
filter.SampleT = 1/filter.SampleF;

filter.Length = 2^nextpow2(2^10*filter.SampleF);
% EndT = str2num(get(handles.EndPeriod, 'String'));
% dfmin = 1/(EndT - 0.5*filter.BandWidth) - 1/(EndT + 0.5*filter.BandWidth);
% filter.Length = 2^nextpow2(filter.SampleF/dfmin);
filter.CtrT = str2num(get(handles.FilterCenterT,'String'));
filter.BandWidth = str2num(get(handles.FilterBandWidthT,'String'));
filter.LowF = (2/filter.SampleF)/(filter.CtrT + 0.5*filter.BandWidth);
filter.HighF = (2/filter.SampleF)/(filter.CtrT - 0.5*filter.BandWidth);    
% Kaiser Window
prompt = {'Please input beta value of Kaiser window (e.g.: ...,6, 7, 8, 9...) (Increasing beta widens the main lobe and decreases the amplitude of the sidelobes'};
title = ['Kaiser Window Parameter for Phase V Dispersion'];
line = 2;
KaiserBeta = inputdlg(prompt,title,line);
filter.KaiserPara = str2num(KaiserBeta{1});   
filter.Data = fir1(filter.Length, [filter.LowF, filter.HighF], kaiser(filter.Length + 1,filter.KaiserPara));
wvtool(filter.Data);


% --- Executes on button press in GroupVIndex.
function GroupVIndex_Callback(hObject, eventdata, handles)
% hObject    handle to GroupVIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GroupVIndex


% --- Executes on button press in TimeVariableFiltering.
function TimeVariableFiltering_Callback(hObject, eventdata, handles)
% hObject    handle to TimeVariableFiltering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TimeVariableFiltering

if get(handles.TimeVariableFiltering, 'Value')
    set(handles.GroupVIndex, 'Value', 1);
else
    set(handles.GroupVIndex, 'Value', 0);
end


function windowPeriodNum_Callback(hObject, eventdata, handles)
% hObject    handle to windowPeriodNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of windowPeriodNum as text
%        str2double(get(hObject,'String')) returns contents of windowPeriodNum as a double


% --- Executes during object creation, after setting all properties.
function windowPeriodNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowPeriodNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minWindowTime_Callback(hObject, eventdata, handles)
% hObject    handle to minWindowTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minWindowTime as text
%        str2double(get(hObject,'String')) returns contents of minWindowTime as a double


% --- Executes during object creation, after setting all properties.
function minWindowTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minWindowTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PeriodVariableWindow.
function PeriodVariableWindow_Callback(hObject, eventdata, handles)
% hObject    handle to PeriodVariableWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global g_refGVwinIndex g_refWinV

% set multiple filtering technique parameters for phase v analysis using
% input reference group v windown

string = 'Input certain group v window for Time-variable Filtering analysis if no group v dispersion picked?';
button = questdlg(string,'Time-variable Filtering analysis using period-variable group v window?', 'Yes', 'No', 'Yes');

if strcmp(button, 'Yes')
    g_refGVwinIndex = 1;
    g_refWinV = RefGVwin(hObject, eventdata, handles);
else
    g_refGVwinIndex = 0;
end

%% read reference group velocity window (period dependent)
function g_refWinV = RefGVwin(hObject, eventdata, handles)
% hObject    handle to RefGDisper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[pfile1, pname1, index] = uigetfile({'*.txt';'*.dat';'*.*'},'Open period-variable group velocity window file', pwd);
files = strcat(pname1,pfile1);

g_refWinV = load(files);
[nT, nc] = size(g_refWinV);
if nc ~= 4
    errordlg('Group v window format error! (need to be T(i) g(i) gmin(i) gmax(i))');
else
    set(handles.MsgEdit,'String','Reference group v window read successfully!');
end

h2 = handles.axes2;
set(gcf,'CurrentAxes',h2);
hold(h2,'off');

plot(g_refWinV(:,1), g_refWinV(:,3),'r--');
hold on; plot(g_refWinV(:,1), g_refWinV(:,4),'r--');
hold on; plot(g_refWinV(:,1), g_refWinV(:,2), 'k', 'LineWidth',2);
xlabel('Period (s)', 'FontSize',8);
ylabel('Group Vel. (km/s)', 'FontSize',8);


%%
% --- Executes on button press in RdCFcn.
function RdCFcn_Callback(hObject, eventdata, handles)
% hObject    handle to RdCFcn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gfcn g_PrevStaLonLat g_FileName g_EGForCFIndex
global g_DataDirectory 

if g_DataDirectory == 0
    g_DataDirectory = pwd;
end

[pfile1, pname1] = uigetfile({'*.dat';'*.*'},'Open Green"s Function/CFs between two stations', g_DataDirectory);
greenfcnfile = strcat(pname1,pfile1);
g_DataDirectory = pname1;
gfcn.FileName = pfile1;
g_FileName = pfile1;

set(handles.MsgEdit, 'String', gfcn.FileName);

datatype = questdlg('Please select input data type (EGF or CF):','Data Type', 'EGF', 'CF', 'EGF');
if strcmp(datatype, 'EGF')
    g_EGForCFIndex = 1;
elseif strcmp(datatype, 'CF')
    g_EGForCFIndex = 2;
end

maxamp = RdGreenFcn(greenfcnfile);

% reset the length of EGF for better showing the waveform
gfcn.WinMinV = str2num(get(handles.GFcnWinMinV,'String'));
maxtime = gfcn.StaDist/gfcn.WinMinV + 250;
ptnum = min(round(maxtime*gfcn.SampleF)+1, gfcn.PtNum);

% ptnum = gfcn.PtNum;
% gfcn.PtNum = min(round(maxtime*gfcn.SampleF)+1, ptnum);
% gfcn.Time = (0:(gfcn.PtNum-1))*gfcn.SampleT;

h2 = handles.axes2;
if ~isnan(maxamp)
    set(gcf,'CurrentAxes',h2);
    hold(h2,'off');
    plot(gfcn.Time(1:ptnum), gfcn.GreenFcn(1,(1:ptnum)),'r-', 'LineWidth', 2);
    hold on
    plot(-gfcn.Time(1:ptnum), gfcn.GreenFcn(2,(1:ptnum)),'b-', 'LineWidth', 2);
    hold off
    xlabel('t (sec)', 'FontSize', 8, 'FontWeight', 'bold');
    set(gca, 'FontSize', 8, 'FontWeight', 'bold');
end   

% UpdataMsgBoxInfo(hObject, eventdata, handles);
set(handles.EditNameSta1,'String','');
set(handles.EditNameSta2,'String','');
set(handles.EditLonSta1,'String',num2str(gfcn.Lon1));
set(handles.EditLonSta2,'String',num2str(gfcn.Lon2));
set(handles.EditLatSta1,'String',num2str(gfcn.Lat1));
set(handles.EditLatSta2,'String',num2str(gfcn.Lat2));
set(handles.EditSrcStaDist,'String',num2str((gfcn.StaDist)));

h1 = handles.axes1;
set(gcf,'CurrentAxes',h1);
hold(h1,'on');
if sum(sum(g_PrevStaLonLat)) ~= 0
    plot(h1, g_PrevStaLonLat(1,1), g_PrevStaLonLat(1,2), 'k^', 'MarkerSize',6, 'MarkerFaceColor','k');
    hold(h1,'on');
    plot(h1, g_PrevStaLonLat(2,1), g_PrevStaLonLat(2,2), 'k^', 'MarkerSize',6, 'MarkerFaceColor','k');
    hold(h1,'on');
    plot(h1, gfcn.Lon1, gfcn.Lat1, 'r^', 'MarkerSize',6, 'MarkerFaceColor','r');
    hold(h1, 'on');
    plot(h1, gfcn.Lon2, gfcn.Lat2, 'b^', 'MarkerSize',6, 'MarkerFaceColor','b');
else
    plot(h1, gfcn.Lon1, gfcn.Lat1, 'r^', 'MarkerSize',6, 'MarkerFaceColor','r');
    hold(h1, 'on');
    plot(h1, gfcn.Lon2, gfcn.Lat2, 'b^', 'MarkerSize',6, 'MarkerFaceColor','b');
end
g_PrevStaLonLat = [gfcn.Lon1, gfcn.Lat1;gfcn.Lon2, gfcn.Lat2];

%% 
% --- Executes on button press in PlotWindowedEGF.
function PlotWindowedEGF_Callback(hObject, eventdata, handles)
% hObject    handle to PlotWindowedEGF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gfcn 

%  function to create signal window function 
[Window, TaperNum] = MakeSignalWindow(hObject, eventdata, handles);
WinWave1(1:gfcn.PtNum) = gfcn.GreenFcn(1,1:gfcn.PtNum).*Window(1:gfcn.PtNum);
WinWave2(1:gfcn.PtNum) = gfcn.GreenFcn(2,1:gfcn.PtNum).*Window(1:gfcn.PtNum);
WinWave3(1:gfcn.PtNum) = WinWave1(1:gfcn.PtNum) + WinWave2(1:gfcn.PtNum);
MaxAmpWinWave3 = max(abs(WinWave3(1:gfcn.PtNum)));
if MaxAmpWinWave3 > 0
    WinWave3(1:gfcn.PtNum) = WinWave3(1:gfcn.PtNum)/MaxAmpWinWave3;
end

h2 = handles.axes2;
set(gcf,'CurrentAxes',h2);
hold(h2,'off');

plot(gfcn.Time, gfcn.GreenFcn(1,1:gfcn.PtNum),'r--', 'LineWidth', 0.1);
hold on
plot(-gfcn.Time, gfcn.GreenFcn(2,1:gfcn.PtNum),'b--', 'LineWidth',0.1);
hold on

plot(gfcn.Time, WinWave1,'r-', 'LineWidth', 2);
hold on
plot(-gfcn.Time, WinWave2,'b-', 'LineWidth', 2);
hold on

plot(gfcn.Time, WinWave3 - 2,'k', 'LineWidth', 2);
hold on
plot(-gfcn.Time, WinWave3 - 2,'k', 'LineWidth', 2);
hold off

xlabel('t (sec)', 'FontSize', 8, 'FontWeight', 'bold');
set(gca, 'FontSize', 8, 'FontWeight', 'bold');

% --- function to create signal window function 
function [Window, TaperNum] = MakeSignalWindow(hObject, eventdata, handles)
global gfcn
global g_WinMaxPtNum g_WinMinPtNum

gfcn.WinMinV = str2num(get(handles.GFcnWinMinV,'String'));
gfcn.WinMaxV = str2num(get(handles.GFcnWinMaxV,'String'));
g_WinMaxPtNum = round(gfcn.SampleF*gfcn.StaDist/gfcn.WinMinV) + 1;
g_WinMinPtNum = round(gfcn.SampleF*gfcn.StaDist/gfcn.WinMaxV) + 1;
if g_WinMaxPtNum > gfcn.PtNum
    g_WinMaxPtNum = gfcn.PtNum-1;
    gfcn.WinMinV = ceil(10*gfcn.StaDist/gfcn.Time(end))/10;
    set(handles.MsgEdit,'String',strcat('Min velocity reset to ',num2str(gfcn.WinMinV)));
end

Window = zeros(1,gfcn.PtNum);
Window(g_WinMinPtNum:g_WinMaxPtNum) = 1;
TaperNum = round(20/gfcn.SampleT);
if TaperNum > (g_WinMinPtNum-1)
    TaperNum1 = g_WinMinPtNum-1;
    Window(1:(g_WinMinPtNum-1)) = sin(0.5*pi*(0:(TaperNum1-1))/TaperNum1);
else
    Window((g_WinMinPtNum-TaperNum):(g_WinMinPtNum-1)) = sin(0.5*pi*(0:(TaperNum-1))/TaperNum);
end
if (g_WinMaxPtNum + TaperNum) < gfcn.PtNum
    Window((g_WinMaxPtNum+1):(g_WinMaxPtNum+TaperNum)) = sin(0.5*pi*((TaperNum-1):-1:0)/TaperNum);
else
    TaperNum2 = gfcn.PtNum - (g_WinMaxPtNum+1) +1;
    Window((g_WinMaxPtNum+1):gfcn.PtNum) = sin(0.5*pi*((TaperNum2-1):-1:0)/TaperNum2);
end

%%
% --- Executes on selection change in PhaseGroupVImg.
function PhaseGroupVImg_Callback(hObject, eventdata, handles)
% hObject    handle to PhaseGroupVImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns PhaseGroupVImg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PhaseGroupVImg
global gfcn cross
global g_WinMaxPtNum g_WinMinPtNum g_HighSNRIndex
global g_GroupTimeIndex g_GroupTime g_DisperType

SignalLength = g_WinMaxPtNum - g_WinMinPtNum + 1; % signal window points

% get which part of the GFcn/CFcn for analysis (causal, acausal, causal + acausal
GreenFcnType = get(handles.PhaseGroupVImg, 'Value');

cross.StartT = str2num(get(handles.StartPeriod,'String'));
cross.EndT = str2num(get(handles.EndPeriod,'String'));
cross.DeltaT = str2num(get(handles.DeltaPeriod,'String'));
cross.NumCtrT = (cross.EndT - cross.StartT)/cross.DeltaT + 1;
TPoint = cross.StartT:cross.DeltaT:cross.EndT;
gfcn.DeltaV = 0.005;
VPoint = gfcn.WinMinV:gfcn.DeltaV:gfcn.WinMaxV;
VImgPt = length(VPoint);
time = ((g_WinMinPtNum-1):1:(g_WinMaxPtNum-1))/gfcn.SampleF;
TravPtV = gfcn.StaDist./time;

%  function to create signal window function 
[Window, TaperNum] = MakeSignalWindow(hObject, eventdata, handles);
%  window EGF
WinWave1(1:gfcn.PtNum) = gfcn.GreenFcn(1,1:gfcn.PtNum).*Window(1:gfcn.PtNum);
WinWave2(1:gfcn.PtNum) = gfcn.GreenFcn(2,1:gfcn.PtNum).*Window(1:gfcn.PtNum);
WinWave3(1:gfcn.PtNum) = (WinWave1(1:gfcn.PtNum) + WinWave2(1:gfcn.PtNum))/2;

% define 150s long noise window after the windowed surface wave
NoisePt = round(150/gfcn.SampleT);
NoiseWinWave = zeros(1, NoisePt);
if ((g_WinMaxPtNum + TaperNum) + (cross.EndT/gfcn.SampleT))< gfcn.PtNum
    nn = (g_WinMaxPtNum + TaperNum + 1):min((g_WinMaxPtNum + TaperNum + NoisePt),gfcn.PtNum);
    if GreenFcnType == 1
        NoiseWinWave(1:length(nn)) = gfcn.GreenFcn(1,nn);
    elseif GreenFcnType == 2
        NoiseWinWave(1:length(nn)) = gfcn.GreenFcn(2,nn);
    elseif GreenFcnType == 3
        NoiseWinWave(1:length(nn)) = (gfcn.GreenFcn(1,nn)+ gfcn.GreenFcn(2,nn))/2;
    end
    NoiseIndex = 1;
    NoiseLength = length(nn); % noise window points
else
    NoiseIndex = 0;
end

WavePtNum = min((g_WinMaxPtNum + TaperNum), gfcn.PtNum);
WinWave = zeros(1, WavePtNum);
switch GreenFcnType
    case {1,4}
        WinWave = WinWave1(1:WavePtNum);  %%%%%%% for WinWave1 A->B
    case {2,5}
        WinWave = WinWave2(1:WavePtNum);   %%%%%%% for WinWave2 B->A
    case {3,6}
        WinWave = WinWave3(1:WavePtNum);   %%%%%%% for WinWave3 A+B
end 

% plot signal (blue) and noise (red)
h3 = handles.axes3;
set(gcf,'CurrentAxes',h3);
hold(h3,'off');
plot((0:length(WinWave)-1)*gfcn.SampleT, WinWave);
if NoiseIndex == 1
    hold on; plot((nn(1)-1:(nn(1)+NoisePt-2))*gfcn.SampleT, NoiseWinWave,'r');
end

hold(h3, 'off');
set(gca, 'FontSize', 8);

% calculate envelope images for signal and noise and estimate SNR
% SNR(T) =  max(signal envelope at period T)/mean(noise envelope at period T)
h4 = handles.axes4;
set(gcf,'CurrentAxes',h4);
hold(h4,'off');
minSNR = str2num(get(handles.minSNR, 'String'));
SNR_T = zeros(1, cross.NumCtrT);

EnvelopeImageSignal = EnvelopeImageCalculation(WinWave, gfcn.SampleF, TPoint, gfcn.StaDist);
AmpS_T = max(EnvelopeImageSignal,[],2);

if NoiseIndex == 1  % noise window long enough
    EnvelopeImageNoise = EnvelopeImageCalculation(NoiseWinWave.*tukeywin(NoisePt,0.2)', gfcn.SampleF, TPoint, gfcn.StaDist);
    for i = 1:cross.NumCtrT
        SNR_T(i) = AmpS_T(i)/mean(EnvelopeImageNoise(i,1:NoiseLength));
    end
    semilogy(TPoint, SNR_T);
    HighSNRIndex = find(SNR_T > minSNR); % find SNR > MinSNR;
    g_HighSNRIndex = HighSNRIndex;
    hold on;
    semilogy(TPoint(HighSNRIndex), SNR_T(HighSNRIndex), 'r*');
    grid on
    xlim([cross.StartT cross.EndT]);
    xlabel('Period (s)','FontSize', 8);
    ylabel('SNR','FontSize', 8); 
    if mean(SNR_T) < 1 || max(SNR_T) < minSNR
        IsDispGood = 0; % overall SNR too low, do not pick dispersion
    else
        IsDispGood = 1; % pick dispersion
    end 
else  % noise window too short
    semilogy(TPoint, AmpS_T);
    xlim([cross.StartT cross.EndT]);
    xlabel('Period (s)','FontSize', 8);
    ylabel('PSD of Signal', 'Color','r','FontSize', 8);
    title('Noise window does not exit!', 'Color','r', 'FontSize',8);    
    HighSNRIndex = 1:cross.NumCtrT; % set high SNR for all periods
    IsDispGood = 1;  % pick dispersion
end
set(gca, 'FontSize', 8);

% fftNumPt = 2^nextpow2(2^8*gfcn.SampleF);
% % estimate Signal to Noise Ratio
% [PxxS,f] = pmtm(WinWave,5/2,fftNumPt,gfcn.SampleF);
% if NoiseIndex == 1
%     [PxxN,f] = pmtm(NoiseWinWave,5/2,fftNumPt,gfcn.SampleF);
%     mm = find(PxxN == 0); nn = find(PxxN > 0);
%     PxxN(mm) = min(PxxN(nn));
%     SNR_f = PxxS./PxxN*NoiseLength/SignalLength;
%     SNR_T = interp1(f, SNR_f, 1./TPoint, 'cubic');   
% end
% 
% h4 = handles.axes4;
% set(gcf,'CurrentAxes',h4);
% hold(h4,'off');
% minSNR = str2num(get(handles.minSNR, 'String'));
% if NoiseIndex == 1
%     semilogy(TPoint, SNR_T);
%     HighSNRIndex = find(SNR_T > minSNR); % find SNR > minSNR;
%     hold on;
%     semilogy(TPoint(HighSNRIndex), SNR_T(HighSNRIndex), 'r*');
%     xlabel('Period (s)','FontSize', 8);
%     ylabel('SNR','FontSize', 8);
% else
%     PxxS_T = interp1(f, PxxS, 1./TPoint, 'v4');
%     semilogy(TPoint, PxxS_T);
%     HighSNRIndex = 1:cross.NumCtrT;
%     xlabel('Period (s)', 'FontSize', 8);
%     ylabel('PSD of Signal', 'Color','r', 'FontSize', 8);
%     title('Noise window too short!', 'Color','r', 'FontSize',8);
% end
% set(gca, 'FontSize', 8);

SNRIndex = zeros(1, cross.NumCtrT);
SNRIndex(HighSNRIndex) = 1;
SNR_T = zeros(1, cross.NumCtrT);

% for those not so bad pts, if they are in the middle of good pts, accept them
for ii = 2:length(SNRIndex)-1
    if SNRIndex(ii) == 0
        if SNR_T(ii) > minSNR/2 && SNRIndex(ii-1) == 1 && SNRIndex(ii+1) == 1
            SNRIndex(ii)=1;
        end
    else
        SNRIndex(ii)=1;
    end
end

g_HighSNRIndex = find(SNRIndex == 1);


switch GreenFcnType
    case {1,2,3}
        
        g_DisperType = 2;  % for group velocity dispersion
        timeptnum = g_WinMinPtNum:1:g_WinMaxPtNum;
        for i = 1:cross.NumCtrT
            gfcn.GroupVImg(1:VImgPt, i) = interp1(TravPtV, EnvelopeImageSignal(i,timeptnum)'/AmpS_T(i), VPoint, 'spline');
%            gfcn.GroupVImg(1:VImgPt, i) = interp1(TravPtV, EnvelopeImageSignal(i,timeptnum)', VPoint, 'spline');

        end
        minamp = min(min(gfcn.GroupVImg(1:VImgPt,1:cross.NumCtrT)));
        maxamp = max(max(gfcn.GroupVImg(1:VImgPt,1:cross.NumCtrT)));
                
        h2 = handles.axes2;
        set(gcf,'CurrentAxes',h2);
        hold(h2,'off');
        imagesc(TPoint, VPoint, gfcn.GroupVImg(1:VImgPt,1:cross.NumCtrT), [minamp, maxamp]);
        set(gca,'YDir','normal');
        xlabel('Period (s)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Group Velocity (km/s)', 'FontSize', 10, 'FontWeight', 'bold');
        set(gca, 'FontSize', 10, 'FontWeight', 'bold');
%         if cross.EndT < 4
%             set(gca, 'XTick',0:0.5:4,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%         elseif cross.EndT >= 4 && cross.EndT <= 20
%             set(gca, 'XTick',0:2:20,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%         elseif cross.EndT > 20
%             set(gca, 'XTick',0:5:cross.EndT,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%         end

    case {4, 5, 6}
        
        g_DisperType = 1;  % for phase velocity dispersion
        
        % function to calculate phase v (c-T) image matrix
        PhaseVImageCalculation(hObject, eventdata, handles, WinWave, g_GroupTimeIndex, g_GroupTime);
        % plot c-T image
        h2 = handles.axes2;
        set(gcf,'CurrentAxes',h2);
        hold(h2,'off');
        imagesc(TPoint, VPoint, gfcn.PhaseVImg, [-1 1]); 
        % cmap = colormap('gray');
        % colormap(cmap(end:-1:1,:)); colorbar;
        set(gca,'YDir','normal');
        set(gca,'YDir','normal','FontSize', 8, 'FontWeight', 'bold','FontName','Arial');
        xlabel('Period (s)', 'FontSize', 8, 'FontWeight', 'bold','FontName','Arial');
        ylabel('Phase Velocity (km/s)', 'FontSize', 8, 'FontWeight', 'bold','FontName','Arial');
%         if cross.EndT < 4
%             set(gca, 'XTick',0:0.5:4,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%         elseif cross.EndT >= 4 && cross.EndT <= 20
%             set(gca, 'XTick',0:2:20,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%         elseif cross.EndT > 20
%             set(gca, 'XTick',0:5:cross.EndT,'YTick',0:0.25:6,'XGrid','on','YGrid','on','TickDir','out');
%         end        
end      


% --- Executes during object creation, after setting all properties.
function PhaseGroupVImg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PhaseGroupVImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PickDisperCurve.
function PickDisperCurve_Callback(hObject, eventdata, handles)
% hObject    handle to PickDisperCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global cross gfcn g_DisperType g_GroupTimeIndex g_GroupTime

h2 = handles.axes2;

VPoint = gfcn.WinMinV:gfcn.DeltaV:gfcn.WinMaxV;
VImgPt = round((gfcn.WinMaxV - gfcn.WinMinV)/gfcn.DeltaV + 1);
TPoint = cross.StartT:cross.DeltaT:cross.EndT;

% for phase or group velocity dispersion
IsDispGood = 1;  % ... = 1: calculate dispersion; ... = other: do not caculate
while IsDispGood == 1
    set(gcf,'CurrentAxes',h2);
    hold(h2,'off');
    if g_DisperType == 1
        imagesc(TPoint, VPoint, gfcn.PhaseVImg(1:VImgPt,1:cross.NumCtrT), [-1, 1]);
        ylabel('Phase Velocity (km/s)', 'FontSize', 8, 'FontWeight', 'bold');
        
    elseif g_DisperType == 2
        minamp = min(min(gfcn.GroupVImg(1:VImgPt,1:cross.NumCtrT)));
        imagesc(TPoint, VPoint, gfcn.GroupVImg(1:VImgPt,1:cross.NumCtrT), [minamp, 1]);
        ylabel('Group Velocity (km/s)', 'FontSize', 8, 'FontWeight', 'bold');
    end
    set(gca,'YDir','normal');
    xlabel('Period (s)', 'FontSize', 8, 'FontWeight', 'bold');
    
    set(gca, 'FontSize', 8, 'FontWeight', 'bold');
    set(gca, 'XTick',0:5:300,'YTick',1:0.25:6,'XGrid','on','YGrid','on');
    [IsDispGood, StartTIndex, EndTIndex] = DisperPick(hObject, handles); % obtain the disperion curve
    
    % if the group velocity dispersion curve has been picked
    if IsDispGood == 2
        g_GroupTimeIndex = 1;
        g_GroupTime = gfcn.StaDist./gfcn.GroupVDisp(1:cross.NumCtrT);
    else
        g_GroupTimeIndex = 0; % no group dispersion/travel time picked
        g_GroupTime = NaN*zeros(1,cross.NumCtrT);
    end
        
    hold off
end
hold off

%% function for pick dispersion curve by clicking on the c-T image
function [IsDispGood, StartTIndex, EndTIndex] = DisperPick(hObject, handles)

global cross g_ClickPoint gfcn g_PrevStaLonLat g_FileName g_HighSNRIndex
global g_DispDirectory 
global g_DisperType

h1 = handles.axes1;
h2 = handles.axes2;
StartTIndex = 1; % corresponds to index of Start Period for saving disper
EndTIndex = cross.NumCtrT;  % corresponds to index of End Period for saving disper

% set(gcf,'CurrentAxes',h2);
axes(h2);

SaveOrNot = questdlg('Do you want to get dispersion curve?','Dispersion Curve','Yes','No','Stop','Yes');


if strcmp(SaveOrNot, 'Yes') == 1
    
    k = 1;
%  old codes for allowing one point click
%     g_ClickPoint = zeros(2,3);
%     while k ~= 0
%         set(handles.MsgEdit,'String','Please click the left button of the mouse on lower figure to select disperion curve!');
%         k = waitforbuttonpress;
%         x = get(h2,'XLim');
%         y = get(h2,'YLim');
%         if g_ClickPoint(1,1) >= x(1,1) && g_ClickPoint(1,1) <= x(1,2) && g_ClickPoint(1,2) >= y(1,1) && g_ClickPoint(1,2) <= y(1,2)
%             hold on
%             plot(h2, g_ClickPoint(1,1), g_ClickPoint(1,2),'*g');
%             k = 0;
%             set(handles.MsgEdit, 'String', 'Successful click!');
%         else      
%             k = 1;
%             set(handles.MsgEdit, 'String', 'Wrong click! Please click on the lower figure!');
%         end
%     end
    
    while k ~= 0
        set(handles.MsgEdit,'String','Left button click for first point and right for last point to select disperion curve!');
        % Initially, the list of points is empty.
        xy = [];
        n = 0;
        % Loop, picking up the points.
        disp('Left mouse button picks points.')
        disp('Right mouse button picks last point.')
        button = 1;
        x = get(h2,'XLim');
        y = get(h2,'YLim');        
        while button == 1
            [xi,yi,button] = ginput(1);
            if xi > x(2) || xi < x(1) || yi > y(2) || yi < y(1)
                set(handles.MsgEdit, 'String', 'Wrong click! Please click on the lower image!');
            else
                hold on; plot(xi,yi,'g*')
                n = n+1;
                xy(:,n) = [xi;yi];  
            end
        end

        if n == 0
            k = 1;
            set(handles.MsgEdit, 'String', 'Wrong click! Please click on the lower image!');
        else
            k = 0;
        end
    end
    
    DisperStartT = cross.StartT;
    DisperEndT = cross.EndT;

    StaDistance = gfcn.StaDist;
    TPoint = cross.StartT:cross.DeltaT:cross.EndT;
%     InitialT = floor((g_ClickPoint(1,1) - cross.StartT)/cross.DeltaT + 1);
%     InitialV = g_ClickPoint(1,2);
%     InitialY = floor((InitialV - gfcn.WinMinV)/gfcn.DeltaV + 1);
    InitialT = round((xy(1,:) - cross.StartT)/cross.DeltaT + 1);
    InitialY = round((xy(2,:) - gfcn.WinMinV)/gfcn.DeltaV + 1);
    nSelectPt = length(InitialT);
    % sort picked points according to increasing periods
    [InitialT, II] = sort(InitialT);
    InitialY = InitialY(II);
    
    set(gcf,'CurrentAxes',h2);
    if g_DisperType == 1  % phase velocity dispersion
        gfcn.PhaseVDisp = zeros(1, cross.NumCtrT);  
        VImgPt = round((gfcn.WinMaxV - gfcn.WinMinV)/gfcn.DeltaV + 1);
        if nSelectPt == 1
            DispPt(1:cross.NumCtrT) = AutoSearch(InitialY, InitialT, gfcn.PhaseVImg(1:VImgPt,1:cross.NumCtrT));
        elseif nSelectPt > 1
            DispPt(1:cross.NumCtrT) = AutoSearchMultiplePoints(InitialY, InitialT, gfcn.PhaseVImg(1:VImgPt,1:cross.NumCtrT));
        end
        gfcn.PhaseVDisp(1:cross.NumCtrT) = gfcn.WinMinV + (DispPt - 1)*gfcn.DeltaV;
        hold on
        plot(TPoint, gfcn.PhaseVDisp(1:cross.NumCtrT), 'k', 'LineWidth', 2);

        wavelength = TPoint.*gfcn.PhaseVDisp(1:cross.NumCtrT);
        DisperV = gfcn.PhaseVDisp(1:cross.NumCtrT);
        hold on;
        plot(TPoint(g_HighSNRIndex), gfcn.PhaseVDisp(g_HighSNRIndex), 'co', 'MarkerSize', 8);
        
        [Tnew, gv] = c2g(TPoint, DisperV); % obtain group vel from phase vel
        hold on; 
        plot(Tnew, gv, 'g-', 'LineWidth', 2);
           
    elseif g_DisperType == 2  % group velocity dispersion
        gfcn.GroupVDisp = zeros(1, cross.NumCtrT);  
        VImgPt = round((gfcn.WinMaxV - gfcn.WinMinV)/gfcn.DeltaV + 1);
        if nSelectPt == 1
            DispPt(1:cross.NumCtrT) = AutoSearch(InitialY, InitialT, gfcn.GroupVImg(1:VImgPt,1:cross.NumCtrT));
        elseif nSelectPt > 1
            DispPt(1:cross.NumCtrT) = AutoSearchMultiplePoints(InitialY, InitialT, gfcn.GroupVImg(1:VImgPt,1:cross.NumCtrT));
        end
        gfcn.GroupVDisp(1:cross.NumCtrT) = gfcn.WinMinV + (DispPt - 1)*gfcn.DeltaV;
        hold on
        plot(TPoint, gfcn.GroupVDisp(1:cross.NumCtrT), 'k', 'LineWidth', 2);

        wavelength = TPoint.*gfcn.GroupVDisp(1:cross.NumCtrT);
        DisperV = gfcn.GroupVDisp(1:cross.NumCtrT);
        hold on;
        plot(TPoint(g_HighSNRIndex), gfcn.GroupVDisp(g_HighSNRIndex), 'co', 'MarkerSize', 8);
        
    end
        
    I1 = find(wavelength > gfcn.StaDist);
    I2 = find(wavelength > gfcn.StaDist/2);
    if length(I1) > 0
        hold on; plot([TPoint(I1(1)) TPoint(I1(1))], [gfcn.WinMinV gfcn.WinMaxV], 'c--', 'LineWidth', 2);
    end
    if length(I2) > 0
        hold on; plot([TPoint(I2(1)) TPoint(I2(1))], [gfcn.WinMinV gfcn.WinMaxV], 'y--', 'LineWidth', 2);
    end  
    minlamdaRatio = str2num(get(handles.minWavelength, 'String'));
    I3 = find((wavelength*minlamdaRatio) > gfcn.StaDist);
    if length(I3) > 0
        hold on; plot([TPoint(I3(1)) TPoint(I3(1))], [gfcn.WinMinV gfcn.WinMaxV], 'r-', 'LineWidth', 2);
    end      

    if g_DisperType == 1
        if get(handles.TimeDomain, 'Value')==1
            dispfilename = strcat('CDisp.T.', g_FileName);
            set(handles.MsgEdit, 'String', strcat('Output Disper File Name: ', dispfilename));
        elseif get(handles.FreqDomain, 'Value')==1
            dispfilename = strcat('CDisp.F.', g_FileName);
            set(handles.MsgEdit, 'String', strcat('Output Disper File Name: ', dispfilename));
        elseif get(handles.AkiSPAC, 'Value')==1 
            dispfilename = strcat('CDisp.A.', g_FileName);
            set(handles.MsgEdit, 'String', strcat('Output Disper File Name: ', dispfilename));
        end
    else
        dispfilename = strcat('GDisp.', g_FileName);
        set(handles.MsgEdit, 'String', strcat('Output Disper File Name: ', dispfilename));
    end
    
    
    
    IsDispGood = 2;
    
    InputOrNotIndex = 1;
    while InputOrNotIndex == 1    
        prompt = {'Enter Start Period:','Enter End Period:'};
        title = ['Set start and end period (s) for saving dispersion data'];
        line = 2; 
        def = {num2str(cross.StartT), num2str(cross.EndT)};
        DisperPeriod = inputdlg(prompt,title,line, def);
        if size(DisperPeriod,1)~=line
            InputOrNotIndex = 1;
        else
            DisperStartT = str2num(DisperPeriod{1});
            DisperEndT = str2num(DisperPeriod{2});
            if isempty(DisperStartT) || isempty(DisperEndT)
                InputOrNotIndex = 1;
            else
                StartTIndex = round((DisperStartT - cross.StartT)/cross.DeltaT) + 1;
                EndTIndex = round((DisperEndT - cross.StartT)/cross.DeltaT) + 1;    
                if DisperStartT >= cross.StartT && DisperEndT <= cross.EndT
                    InputOrNotIndex = 0;
                else
                    InputOrNotIndex = 1;
                end
            end
        end
    end
       

    %write T-V to file 
    if g_DispDirectory == 0
        g_DispDirectory = pwd;
    end
    
    if sum(sum(g_PrevStaLonLat)) ~= 0
        set(gcf,'CurrentAxes',h1);
        hold(h1,'on');
        plot(g_PrevStaLonLat(1:2,1), g_PrevStaLonLat(1:2,2), 'b-');
%         plot(h1, g_PrevStaLonLat(1:2,1), g_PrevStaLonLat(1:2,2), 'k^', 'MarkerSize',6, 'MarkerFaceColor','k');
    end

    TVfile = strcat(g_DispDirectory, '/', dispfilename);    
    ftv = fopen(TVfile,'w');
    fprintf(ftv,'%f     ', gfcn.Lon1);
    fprintf(ftv,'%f\n', gfcn.Lat1);
    fprintf(ftv,'%f     ', gfcn.Lon2);
    fprintf(ftv,'%f\n', gfcn.Lat2);
    DataExistIndex = 0;
    minlamdaRatio = str2num(get(handles.minWavelength, 'String'));
    set(gcf,'CurrentAxes',h2);
    hold(h2,'on');
    
    NewEndTIndex = EndTIndex;
    for i = 1:cross.NumCtrT
       if (i >= StartTIndex) & (i <=EndTIndex)
           wavelength = DisperV(i)*TPoint(i);
           if StaDistance >= minlamdaRatio*wavelength 
                DataExistIndex = 1;
                if g_DisperType == 1
                    % make correction of phase velocity at near field for time/frequency 
                    % domain dispersion measurements due to the use of far-field
                    % approximation of surface wave propagation
                    if get(handles.TimeDomain, 'Value')==1  ||  get(handles.FreqDomain, 'Value')==1   
                        kr = 2*pi*StaDistance/wavelength; 
                        G2d = complex(-bessely(0,kr),-besselj(0,kr));  % exact Green's function
                        G2d_ff = complex(cos(kr+pi/4), -sin(kr+pi/4));  % Green's function at far field
                        phiG2d = -angle(G2d)*57.3; 
                        phiG2d_ff = -angle(G2d_ff)*57.3;
                        phiDiff = phiG2d_ff - phiG2d; % phase difference in degree for G2d
                        if phiDiff >= 180
                            phiDiff = phiDiff - 360;
                        elseif phiDiff < -180
                            phiDiff = phiDiff + 360;
                        end
                        t_corrected = (StaDistance/DisperV(i) + TPoint(i)/8) - (phiDiff/360)*TPoint(i);
                        DisperV(i) = StaDistance/(t_corrected - TPoint(i)/8);           
                    end

                    hold(h2,'on');
                    plot(h2, TPoint(i), DisperV(i), 'ro','MarkerSize', 5, 'MarkerFaceColor','r');
                else
                    hold(h2,'on');
                    plot(h2, TPoint(i), DisperV(i), 'ro','MarkerSize', 5, 'MarkerFaceColor','r');
                end

                fprintf(ftv,'%4.3f    %4.3f    %4.3f   %4d\n',[TPoint(i) DisperV(i) 0.0  1]);

           else
               NewEndTIndex = min(NewEndTIndex, i);
    %             % plot period high limit for the dispersion measurements
    %             plot([TPoint(i) TPoint(i)], [gfcn.WinMaxV gfcn.WinMinV],
    %             'g-', 'LineWidth',2);
                fprintf(ftv,'%4.3f    %4.3f    %4.3f   %4d\n',[TPoint(i) 0.0  0.0  0]);
           end
       else
           fprintf(ftv,'%4.3f    %4.3f    %4.3f   %4d\n',[TPoint(i) 0.0  0.0  0]);
       end
           
    end
    fclose(ftv);
    EndTIndex = NewEndTIndex;
    
    ReviseOrNot = questdlg('Want to revise the saved dispersion?','Revise dispersion data','Yes','No','Stop','No');
    if strcmp(ReviseOrNot, 'Yes') == 1
        IsDispGood = 1;
        delete(TVfile);
        set(handles.MsgEdit, 'String', 'No Disperion Data Written! Disper File Deleted!');
    elseif strcmp(ReviseOrNot, 'No') == 1
        if DataExistIndex == 1
            set(gcf,'CurrentAxes',h1);
            hold(h1,'on');
            plot([gfcn.Lon1, gfcn.Lon2], [gfcn.Lat1, gfcn.Lat2], 'r-');
%             plot(h1, [gfcn.Lon1, gfcn.Lon2], [gfcn.Lat1, gfcn.Lat2], 'k^', 'MarkerSize',6, 'MarkerFaceColor','g');
            g_PrevStaLonLat = [gfcn.Lon1, gfcn.Lat1;gfcn.Lon2, gfcn.Lat2];
            IsDispGood = 2;  % dispersion picked
        else
            IsDispGood = 0;  % no dispersion picked
        end
    elseif strcmp(ReviseOrNot, 'Stop') == 1
        IsDispGood = 0;
        set(handles.StopDataProcessing, 'Value', 1); 
        set(handles.MsgEdit, 'String', 'Processing stopped!');        
    end

elseif strcmp(SaveOrNot, 'No') == 1
    IsDispGood = 0;
elseif strcmp(SaveOrNot, 'Stop') == 1
    IsDispGood = 0;
    set(handles.StopDataProcessing, 'Value', 1); 
    set(handles.MsgEdit, 'String', 'Processing stopped!');   
end

% ---------------------------
function [Tnew, gv] = c2g(T, c)
% function of computing group velocity dispersion from phase velocity dispersion
% using the equation g = c + k(dc/dk), where k = w/c = 2*pi/(cT)
% T: vector of period (monotonically increasing or decreasing)
% c: vector of phase velocity corresponding to T

% 3 point average of phase v disper in order to obtain smooth group v
c(2:end-1) = (c(1:end-2) + c(2:end-1) + c(3:end))/3;  

k = 2*pi./(c.*T); % wavenumber
% dc = c(3:end) - c(1:end-2);
% dk = k(3:end) - k(1:end-2);
% gv = c(2:end-1) + k(2:end-1).*(dc./dk);

dc1 = c(2:end-1) - c(1:end-2);
dk1 = k(2:end-1) - k(1:end-2);
gv1 = c(2:end-1) + k(2:end-1).*(dc1./dk1); % from left derivative 

dc2 = c(3:end) - c(2:end-1);
dk2 = k(3:end) - k(2:end-1);
gv2 = c(2:end-1) + k(2:end-1).*(dc2./dk2); % from right derivative
gv = (gv1 + gv2)/2;  % mean group velocity between T(2) and T(end-1)

Tnew = T(2:end-1);

% ---------------------------
function EditNameSta1_Callback(hObject, eventdata, handles)
% hObject    handle to EditNameSta1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditNameSta1 as text
%        str2double(get(hObject,'String')) returns contents of EditNameSta1 as a double


% --- Executes during object creation, after setting all properties.
function EditNameSta1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditNameSta1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditNameSta2_Callback(hObject, eventdata, handles)
% hObject    handle to EditNameSta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditNameSta2 as text
%        str2double(get(hObject,'String')) returns contents of EditNameSta2 as a double


% --- Executes during object creation, after setting all properties.
function EditNameSta2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditNameSta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditLonSta1_Callback(hObject, eventdata, handles)
% hObject    handle to EditLonSta1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditLonSta1 as text
%        str2double(get(hObject,'String')) returns contents of EditLonSta1 as a double


% --- Executes during object creation, after setting all properties.
function EditLonSta1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditLonSta1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditLonSta2_Callback(hObject, eventdata, handles)
% hObject    handle to EditLonSta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditLonSta2 as text
%        str2double(get(hObject,'String')) returns contents of EditLonSta2 as a double


% --- Executes during object creation, after setting all properties.
function EditLonSta2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditLonSta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditLatSta1_Callback(hObject, eventdata, handles)
% hObject    handle to EditLatSta1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditLatSta1 as text
%        str2double(get(hObject,'String')) returns contents of EditLatSta1 as a double


% --- Executes during object creation, after setting all properties.
function EditLatSta1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditLatSta1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditLatSta2_Callback(hObject, eventdata, handles)
% hObject    handle to EditLatSta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditLatSta2 as text
%        str2double(get(hObject,'String')) returns contents of EditLatSta2 as a double


% --- Executes during object creation, after setting all properties.
function EditLatSta2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditLatSta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditSrcStaDist_Callback(hObject, eventdata, handles)
% hObject    handle to EditSrcStaDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSrcStaDist as text
%        str2double(get(hObject,'String')) returns contents of EditSrcStaDist as a double


% --- Executes during object creation, after setting all properties.
function EditSrcStaDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSrcStaDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to read CFcn data file
function maxamp = RdGreenFcn(greenfcnfile)
global gfcn g_EGForCFIndex
fgfcn = fopen(greenfcnfile, 'r');
gfcn.Lon1 = fscanf(fgfcn, '%f', 1);
gfcn.Lat1 = fscanf(fgfcn, '%f', 1);
temp = fgetl(fgfcn);
if strcmp(temp,'')==1
    sta1_elev = 0.0;
else
    sta1_elev = str2num(temp);
end
gfcn.Lon2 = fscanf(fgfcn, '%f', 1);
gfcn.Lat2 = fscanf(fgfcn, '%f', 1);
temp = fgetl(fgfcn);
if strcmp(temp,'')==1
    sta2_elev = 0.0;
else
    sta2_elev = str2num(temp);
end
GreenFcn = fscanf(fgfcn,'%f',[3,inf]);
gfcn.PtNum = size(GreenFcn, 2);
gfcn.Time = GreenFcn(1,:);
maxamp = max(max(GreenFcn(2,:)), max(GreenFcn(3,:)));
if maxamp > 0
    GreenFcn(2,:) = GreenFcn(2,:)/maxamp;
    GreenFcn(3,:) = GreenFcn(3,:)/maxamp;
end

% using hilbert tranform to obtain EGF from CF if reading CF
if g_EGForCFIndex == 2
    GreenFcn(2,:) = imag(hilbert(GreenFcn(2,:)));
    GreenFcn(3,:) = imag(hilbert(GreenFcn(3,:)));
end

gfcn.GreenFcn = GreenFcn(2:3,:);
fclose(fgfcn);
gfcn.SampleF = 1/(gfcn.Time(2) - gfcn.Time(1));
gfcn.SampleT = gfcn.Time(2) - gfcn.Time(1);

if gfcn.Lon1 < 0
    gfcn.Lon1 = 360 + gfcn.Lon1;
end
if gfcn.Lon2 < 0
    gfcn.Lon2 = 360 + gfcn.Lon2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified by Ning Gu
% [x,y,~] = deg2utm([gfcn.Lat1;gfcn.Lat2],[gfcn.Lon1;gfcn.Lon2]);
% tempDist = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2)/1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tempDist = deg2km(distance(gfcn.Lat1, gfcn.Lon1, gfcn.Lat2, gfcn.Lon2));
staElevDiff = abs(sta2_elev - sta1_elev)/1000;
% correct station distance due to elevation difference
gfcn.StaDist = sqrt(tempDist*tempDist + staElevDiff*staElevDiff);

% for lunar data considering the radius difference between earth and moon
% gfcn.StaDist = gfcn.StaDist*1738.1/6371;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % updata message box information
function UpdataMsgBoxInfo(hObject, eventdata, handles)
global gfcn
set(handles.EditNameSta1,'String',gfcn.Name1);
set(handles.EditNameSta2,'String',gfcn.Name2);
set(handles.EditLonSta1,'String',num2str(gfcn.Lon1));
set(handles.EditLonSta2,'String',num2str(gfcn.Lon2));
set(handles.EditLatSta1,'String',num2str(gfcn.Lat1));
set(handles.EditLatSta2,'String',num2str(gfcn.Lat2));
set(handles.EditSrcStaDist,'String',num2str(gfcn.StaDist,'%.3f'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Automatically search arrival time line on a image ---
function ArrPt = AutoSearch(InitialY, InitialX, ImageData)
% Input: InitialY  InitlaX   ImageData
% OutPut: ArrPt

YSize = size(ImageData, 1);
XSize = size(ImageData, 2);

ArrPt = zeros(1,XSize);

% Center_T search up
step = 3;
point_left = 0;
point_right = 0;

for i = InitialX:XSize
	index1 = 0;
	index2 = 0; 
	point_left = InitialY;
	point_right = InitialY;
	while index1 == 0
        point_left_new = max(1, point_left - step);
	    if ImageData(point_left,i) < ImageData(point_left_new,i)
		  point_left = point_left_new;
%          point_right = point_right - step;
		else
		   index1 = 1;
		   point_left = point_left_new;
        end
    end
    while index2 == 0
        point_right_new = min(point_right + step, YSize);
		if ImageData(point_right,i) < ImageData(point_right_new,i)
		   point_right = point_right_new;
%           point_left = point_left + step;
		else
           index2=1;
		   point_right = point_right_new;
		end
    end

    [MaxAmp, index_max] = max(ImageData(point_left:point_right,i));
    ArrPt(i) = index_max + point_left - 1;
    InitialY = ArrPt(i);
        
end  %end for

% Center_T search down

InitialY = ArrPt(InitialX);
for i = (InitialX - 1):(-1):1
	index1 = 0;
	index2 = 0; 
    point_left = InitialY;
	point_right = InitialY;

	while index1 == 0
        point_left_new = max(1, point_left - step);
	    if ImageData(point_left,i) < ImageData(point_left_new,i)
		  point_left = point_left_new;
%          point_right = point_right - step;
		else
		   index1 = 1;
		   point_left = point_left_new;
        end
    end
    while index2 == 0
        point_right_new = min(point_right + step, YSize);
		if ImageData(point_right,i) < ImageData(point_right_new,i)
		   point_right = point_right_new;
%           point_left = point_left + step;
		else
           index2=1;
		   point_right = point_right_new;
		end
    end
    
    [MaxAmp, index_max] = max(ImageData(point_left:point_right,i));
    ArrPt(i) = index_max + point_left - 1;
    InitialY = ArrPt(i);
    
end  %end for

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Automatically search arrival time line on a image with multiple points constraints---
function ArrPt = AutoSearchMultiplePoints(ptY, ptX, ImageData)
% Input: InitialY  InitlaX   ImageData
% OutPut: ArrPt

% sort ptX and ptY according to the increasing of ptX
[ptX, II] = sort(ptX);
ptY = ptY(II);

YSize = size(ImageData, 1);
XSize = size(ImageData, 2);
nPt = length(ptX);  % number of input searching points

ArrPt = zeros(1,XSize);

% X searching up for the point with maximum X
step = 3;
point_left = 0;
point_right = 0;

InitialX = ptX(nPt);
InitialY = ptY(nPt);
for i = ptX(nPt):XSize
	index1 = 0;
	index2 = 0; 
	point_left = InitialY;
	point_right = InitialY;
	while index1 == 0
        point_left_new = max(1, point_left - step);
	    if ImageData(point_left,i) < ImageData(point_left_new,i)
		  point_left = point_left_new;
%          point_right = point_right - step;
		else
		   index1 = 1;
		   point_left = point_left_new;
        end
    end
    while index2 == 0
        point_right_new = min(point_right + step, YSize);
		if ImageData(point_right,i) < ImageData(point_right_new,i)
		   point_right = point_right_new;
%           point_left = point_left + step;
		else
           index2=1;
		   point_right = point_right_new;
		end
    end

    [MaxAmp, index_max] = max(ImageData(point_left:point_right,i));
    ArrPt(i) = index_max + point_left - 1;
    InitialY = ArrPt(i);
        
end  %end for

% X searching down for the point with maximum X. There will other points
% with smaller X which will act as internal constraints for the searching
% process

InitialX = ptX(nPt);
InitialY = ArrPt(ptX(nPt));
midX = ptX(nPt-1);
midY = ptY(nPt-1);
kk = 0;
for i = ptX(nPt):(-1):1
	index1 = 0;
	index2 = 0;
    
    if i == midX
        InitialY = midY;
        kk = kk + 1;
        if (nPt - kk) > 1
            midX = ptX(nPt - kk - 1);
            midY = ptY(nPt - kk - 1);
        end
    end
            
    point_left = InitialY;
	point_right = InitialY;

	while index1 == 0
        point_left_new = max(1, point_left - step);
	    if ImageData(point_left,i) < ImageData(point_left_new,i)
		  point_left = point_left_new;
%          point_right = point_right - step;
		else
		   index1 = 1;
		   point_left = point_left_new;
        end
    end
    while index2 == 0
        point_right_new = min(point_right + step, YSize);
		if ImageData(point_right,i) < ImageData(point_right_new,i)
		   point_right = point_right_new;
%           point_left = point_left + step;
		else
           index2=1;
		   point_right = point_right_new;
		end
    end
    
    [MaxAmp, index_max] = max(ImageData(point_left:point_right,i));
    ArrPt(i) = index_max + point_left - 1;
    InitialY = ArrPt(i);
    
end  %end for


% --- Executes on button press in RayleighWave.
function RayleighWave_Callback(hObject, eventdata, handles)
% hObject    handle to RayleighWave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RayleighWave
if get(handles.RayleighWave, 'Value')
    set(handles.LoveWave, 'Value', 0);
else
    set(handles.RayleighWave, 'Value', 1);
end


% --- Executes on button press in LoveWave.
function LoveWave_Callback(hObject, eventdata, handles)
% hObject    handle to LoveWave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LoveWave
if get(handles.LoveWave, 'Value')
    set(handles.RayleighWave, 'Value', 0);
else
    set(handles.LoveWave, 'Value', 1);
end


% --- Executes on button press in RefCGDisper.
function RefCGDisper_Callback(hObject, eventdata, handles)
% hObject    handle to RefCGDisper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global g_refphasedisp g_refgroupdisp

% input phase velocity reference dispersion and ranges
[pfile1, pname1, index] = uigetfile({'*.txt';'*.dat';'*.*'},'Open reference phase velocity file', pwd);
files = strcat(pname1,pfile1);

refcdisp = load(files);
g_refphasedisp = refcdisp;
[nT, nc] = size(g_refphasedisp);
if nc == 3
    g_refphasedisp(:,2) = refcdisp(:,2) - refcdisp(:,3);
    g_refphasedisp(:,3) = refcdisp(:,2) + refcdisp(:,3);
    set(handles.MsgEdit,'String','Reference dispersion read successfully!');
elseif nc == 2
    dc = str2num(get(handles.deltaPhaseV, 'String'));
    g_refphasedisp(:,2) = refcdisp(:,2) - dc;
    g_refphasedisp(:,3) = refcdisp(:,2) + dc;
    set(handles.MsgEdit,'String','Reference dispersion read successfully!');
end

h2 = handles.axes2;
set(gcf,'CurrentAxes',h2);
hold(h2,'off');

plot(g_refphasedisp(:,1), g_refphasedisp(:,2),'b--');
hold on; plot(g_refphasedisp(:,1), g_refphasedisp(:,3),'b--');
hold on; plot(g_refphasedisp(:,1), refcdisp(:,2), 'b', 'LineWidth',2);
xlabel('Period (s)', 'FontSize',8);
ylabel('Ref Velocity (km/s)', 'FontSize',8);

% input group velocity reference dispersion and ranges
if get(handles.GroupVIndex,'Value')  % also measure group v dispersion
    [pfile1, pname1, index] = uigetfile({'*.txt';'*.dat';'*.*'},'Open reference group velocity file', pwd);
    files = strcat(pname1,pfile1);

    refgdisp = load(files);
    g_refgroupdisp = refgdisp;
    [nT, nc] = size(g_refgroupdisp);
    if nc == 3
        g_refgroupdisp(:,2) = refgdisp(:,2) - refgdisp(:,3);
        g_refgroupdisp(:,3) = refgdisp(:,2) + refgdisp(:,3);
        set(handles.MsgEdit,'String','Reference dispersion read successfully!');
    elseif nc == 2
        dc = str2num(get(handles.deltaPhaseV, 'String'));
        g_refgroupdisp(:,2) = refgdisp(:,2) - dc;
        g_refgroupdisp(:,3) = refgdisp(:,2) + dc;
        set(handles.MsgEdit,'String','Reference dispersion read successfully!');
    end

    hold(h2,'on');

    plot(g_refgroupdisp(:,1), g_refgroupdisp(:,2),'g--');
    hold on; plot(g_refgroupdisp(:,1), g_refgroupdisp(:,3),'g--');
    hold on; plot(g_refgroupdisp(:,1), refgdisp(:,2), 'g', 'LineWidth',2);
    xlabel('Period (s)', 'FontSize',8);
    ylabel('Ref Velocity (km/s)', 'FontSize',8);
end
    


function deltaPhaseV_Callback(hObject, eventdata, handles)
% hObject    handle to deltaPhaseV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deltaPhaseV as text
%        str2double(get(hObject,'String')) returns contents of deltaPhaseV as a double


% --- Executes during object creation, after setting all properties.
function deltaPhaseV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltaPhaseV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function minWavelength_Callback(hObject, eventdata, handles)
% hObject    handle to minWavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minWavelength as text
%        str2double(get(hObject,'String')) returns contents of minWavelength as a double


% --- Executes during object creation, after setting all properties.
function minWavelength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minWavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% caculate distance, modified by Ning Gu
function  [x,y,utmzone] = deg2utm(Lat,Lon)
% -------------------------------------------------------------------------
% [x,y,utmzone] = deg2utm(Lat,Lon)
%
% Description: Function to convert lat/lon vectors into UTM coordinates (WGS84).
% Some code has been extracted from UTM.m function by Gabriel Ruiz Martinez.
%
% Inputs:
%    Lat: Latitude vector.   Degrees.  +ddd.ddddd  WGS84
%    Lon: Longitude vector.  Degrees.  +ddd.ddddd  WGS84
%
% Outputs:
%    x, y , utmzone.   See example
%
% Example 1:
%    Lat=[40.3154333; 46.283900; 37.577833; 28.645650; 38.855550; 25.061783];
%    Lon=[-3.4857166; 7.8012333; -119.95525; -17.759533; -94.7990166; 121.640266];
%    [x,y,utmzone] = deg2utm(Lat,Lon);
%    fprintf('%7.0f ',x)
%       458731  407653  239027  230253  343898  362850
%    fprintf('%7.0f ',y)
%      4462881 5126290 4163083 3171843 4302285 2772478
%    utmzone =
%       30 T
%       32 T
%       11 S
%       28 R
%       15 S
%       51 R
%
% Example 2: If you have Lat/Lon coordinates in Degrees, Minutes and Seconds
%    LatDMS=[40 18 55.56; 46 17 2.04];
%    LonDMS=[-3 29  8.58;  7 48 4.44];
%    Lat=dms2deg(mat2dms(LatDMS)); %convert into degrees
%    Lon=dms2deg(mat2dms(LonDMS)); %convert into degrees
%    [x,y,utmzone] = deg2utm(Lat,Lon)
%
% Author: 
%   Rafael Palacios
%   Universidad Pontificia Comillas
%   Madrid, Spain
% Version: Apr/06, Jun/06, Aug/06, Aug/06
% Aug/06: fixed a problem (found by Rodolphe Dewarrat) related to southern 
%    hemisphere coordinates. 
% Aug/06: corrected m-Lint warnings
%-------------------------------------------------------------------------

% Argument checking
%
error(nargoutchk(2, 2, nargin));  %2 arguments required
n1=length(Lat);
n2=length(Lon);
if (n1~=n2)
   error('Lat and Lon vectors should have the same length');
end


% Memory pre-allocation
%
x=zeros(n1,1);
y=zeros(n1,1);
utmzone(n1,:)='60 X';

% Main Loop
%
for i=1:n1
   la=Lat(i);
   lo=Lon(i);

   sa = 6378137.000000 ; sb = 6356752.314245;
         
   %e = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sa;
   e2 = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sb;
   e2cuadrada = e2 ^ 2;
   c = ( sa ^ 2 ) / sb;
   %alpha = ( sa - sb ) / sa;             %f
   %ablandamiento = 1 / alpha;   % 1/f

   lat = la * ( pi / 180 );
   lon = lo * ( pi / 180 );

   Huso = fix( ( lo / 6 ) + 31);
   S = ( ( Huso * 6 ) - 183 );
   deltaS = lon - ( S * ( pi / 180 ) );

   if (la<-72), Letra='C';
   elseif (la<-64), Letra='D';
   elseif (la<-56), Letra='E';
   elseif (la<-48), Letra='F';
   elseif (la<-40), Letra='G';
   elseif (la<-32), Letra='H';
   elseif (la<-24), Letra='J';
   elseif (la<-16), Letra='K';
   elseif (la<-8), Letra='L';
   elseif (la<0), Letra='M';
   elseif (la<8), Letra='N';
   elseif (la<16), Letra='P';
   elseif (la<24), Letra='Q';
   elseif (la<32), Letra='R';
   elseif (la<40), Letra='S';
   elseif (la<48), Letra='T';
   elseif (la<56), Letra='U';
   elseif (la<64), Letra='V';
   elseif (la<72), Letra='W';
   else Letra='X';
   end

   a = cos(lat) * sin(deltaS);
   epsilon = 0.5 * log( ( 1 +  a) / ( 1 - a ) );
   nu = atan( tan(lat) / cos(deltaS) ) - lat;
   v = ( c / ( ( 1 + ( e2cuadrada * ( cos(lat) ) ^ 2 ) ) ) ^ 0.5 ) * 0.9996;
   ta = ( e2cuadrada / 2 ) * epsilon ^ 2 * ( cos(lat) ) ^ 2;
   a1 = sin( 2 * lat );
   a2 = a1 * ( cos(lat) ) ^ 2;
   j2 = lat + ( a1 / 2 );
   j4 = ( ( 3 * j2 ) + a2 ) / 4;
   j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat) ) ^ 2) ) / 3;
   alfa = ( 3 / 4 ) * e2cuadrada;
   beta = ( 5 / 3 ) * alfa ^ 2;
   gama = ( 35 / 27 ) * alfa ^ 3;
   Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
   xx = epsilon * v * ( 1 + ( ta / 3 ) ) + 500000;
   yy = nu * v * ( 1 + ta ) + Bm;

   if (yy<0)
       yy=9999999+yy;
   end

   x(i)=xx;
   y(i)=yy;
   utmzone(i,:)=sprintf('%02d %c',Huso,Letra);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

