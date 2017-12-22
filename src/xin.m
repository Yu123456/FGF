function varargout = xin(varargin)
%XIN M-file for xin.fig
%      XIN, by itself, creates a new XIN or raises the existing
%      singleton*.
%
%      H = XIN returns the handle to a new XIN or the handle to
%      the existing singleton*.
%
%      XIN('Property','Value',...) creates a new XIN using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to xin_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      XIN('CALLBACK') and XIN('CALLBACK',hObject,...) call the
%      local function named CALLBACK in XIN.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help xin

% Last Modified by GUIDE v2.5 17-Apr-2010 10:16:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @xin_OpeningFcn, ...
                   'gui_OutputFcn',  @xin_OutputFcn, ...
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


% --- Executes just before xin is made visible.
function xin_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for xin
handles.output = hObject;
global num
num=0;
set(handles.start,'visible','on');
set(handles.pause,'visible','off');
set(handles.close,'visible','on');
set(handles.stop,'visible','off');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes xin wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = xin_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.

% --- Executes during object creation, after setting all properties.


% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global num
t=0:0.1:20;
m=201;
x=sin(t);
y=x;
n=num+1;
set(handles.start,'visible','off');
set(handles.pause,'visible','on');
for i=n:m
plot(i,y(i));
num=i;
hold on
grid on
pause(0.01);
end



function shuzi_Callback(hObject, eventdata, handles)
% hObject    handle to shuzi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shuzi as text
%        str2double(get(hObject,'String')) returns contents of shuzi as a double


% --- Executes during object creation, after setting all properties.
function shuzi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shuzi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.


% --- Executes on button press in pause.
function pause_Callback(hObject, eventdata, handles)
% hObject    handle to pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
t=0:0.1:20;
m=201;
x=sin(t);
y=x;
set(handles.start,'visible','on');
set(handles.pause,'visible','off');
%pause
uiwait
% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global num 
num=0;
close();
