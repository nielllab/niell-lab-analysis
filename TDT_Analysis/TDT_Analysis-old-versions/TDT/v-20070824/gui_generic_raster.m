function varargout = gui_generic_raster(varargin)
% GUI_GENERIC_RASTER M-file for gui_generic_raster.fig
%      GUI_GENERIC_RASTER, by itself, creates a new GUI_GENERIC_RASTER or raises the existing
%      singleton*.
%
%      H = GUI_GENERIC_RASTER returns the handle to a new GUI_GENERIC_RASTER or the handle to
%      the existing singleton*.
%
%      GUI_GENERIC_RASTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GENERIC_RASTER.M with the given input arguments.
%
%      GUI_GENERIC_RASTER('Property','Value',...) creates a new GUI_GENERIC_RASTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_generic_raster_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_generic_raster_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_generic_raster

% Last Modified by GUIDE v2.5 24-Aug-2007 13:52:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_generic_raster_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_generic_raster_OutputFcn, ...
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


% --- Executes just before gui_generic_raster is made visible.
function gui_generic_raster_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_generic_raster (see VARARGIN)

% Choose default command line output for gui_generic_raster
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_generic_raster wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.figure1,'Units','Pixels')
set(handles.figure1,'Position',[1280 770 400 260]);

% --- Outputs from this function are returned to the command line.
function varargout = gui_generic_raster_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function Tank_Name_Callback(hObject, eventdata, handles)
% hObject    handle to Tank_Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tank_Name as text
%        str2double(get(hObject,'String')) returns contents of Tank_Name as a double


% --- Executes during object creation, after setting all properties.
function Tank_Name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tank_Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Block_Name_Callback(hObject, eventdata, handles)
% hObject    handle to Block_Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Block_Name as text
%        str2double(get(hObject,'String')) returns contents of Block_Name as a double


% --- Executes during object creation, after setting all properties.
function Block_Name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Block_Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Channel_List_Callback(hObject, eventdata, handles)
% hObject    handle to Channel_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Channel_List as text
%        str2double(get(hObject,'String')) returns contents of Channel_List as a double


% --- Executes during object creation, after setting all properties.
function Channel_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Channel_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bin_size_Callback(hObject, eventdata, handles)
% hObject    handle to bin_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bin_size as text
%        str2double(get(hObject,'String')) returns contents of bin_size as a double


% --- Executes during object creation, after setting all properties.
function bin_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bin_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_rate_Callback(hObject, eventdata, handles)
% hObject    handle to max_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_rate as text
%        str2double(get(hObject,'String')) returns contents of max_rate as a double


% --- Executes during object creation, after setting all properties.
function max_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nrows_Callback(hObject, eventdata, handles)
% hObject    handle to nrows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nrows as text
%        str2double(get(hObject,'String')) returns contents of nrows as a double


% --- Executes during object creation, after setting all properties.
function nrows_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nrows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ncols_Callback(hObject, eventdata, handles)
% hObject    handle to ncols (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ncols as text
%        str2double(get(hObject,'String')) returns contents of ncols as a double


% --- Executes during object creation, after setting all properties.
function ncols_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ncols (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunBtn.
function RunBtn_Callback(hObject, eventdata, handles)
% hObject    handle to RunBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Block_Name = get(handles.Block_Name,'String')
Tank_Name = get(handles.Tank_Name,'String')
channel_list = str2num(get(handles.Channel_List,'String'));
bin_size = str2num(get(handles.bin_size,'String'));
max_rate = str2num(get(handles.max_rate,'String'));
nrows = str2num(get(handles.nrows,'String'));
ncols = str2num(get(handles.ncols,'String'));

generic_raster(Tank_Name, Block_Name,channel_list, nrows, ncols, bin_size,max_rate)



% --- Executes on button press in BlockBrowse.
function TankBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to BlockBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'lastblock') %% already set directory
    base = handles.lastblock;
else
    base = 'C:\data';
end

pname = uigetdir(base,'Select Block:');   %%% start location of tanks for search

if (pname ~= 0)
    delims = strfind(pname,'\');
    Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1);
    Block_Name = pname(delims(length(delims))+1 :length(pname));

    set(handles.Tank_Name','String',Tank_Name);
    set(handles.Block_Name','String',Block_Name);

    handles.lastblock = pname; %%pname(1:delims(length(delims)));
    guidata(hObject, handles);
end


% --- Executes on button press in clearfigures.
function clearfigures_Callback(hObject, eventdata, handles)
% hObject    handle to clearfigures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global all_generic_raster_figures;

channel_list = str2num(get(handles.Channel_List,'String'));
numchannels = length(channel_list);

for f = 1:2*numchannels
    close(f);
end
