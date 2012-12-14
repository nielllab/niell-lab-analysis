function varargout = mainCVAP(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mainCVAP_OpeningFcn, ...
                   'gui_OutputFcn',  @mainCVAP_OutputFcn, ...
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

% --- Executes just before mainCVAP is made visible.
function mainCVAP_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for mainCVAP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
axes(handles.axes1);
plot(0.5,0.5,'.');
cla;
set(gca,'XTickLabel','');
set(gca,'YTickLabel','');
text(0.03,0.94, 'Plotting Window','FontSize',12);
handles.data = [];
handles.ialg = [];
handles.ialgvalue = [];
handles.indexname = [];
handles.labels = [];
handles.truelabels = [];
handles.distype = [];
handles.runstate = [];
handles.datastate =[];
handles.validstate = [];
guidata(hObject,handles);

% --- Outputs from this function are returned to the command line.
function varargout = mainCVAP_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function pushClustering_Callback(hObject, eventdata, handles)
if isempty(handles.datastate)
   Q = ['Clustering state: No data for clustering. Please load data first !'];
   set(handles.Outext1, 'String', Q);
   return;
end
single = get(handles.popClutertimes, 'Value');
if single == 1
  k2 = get(handles.edit1, 'String');
else
  k2 = get(handles.edit2, 'String');
end
valid_runclustering;

function pushClearCluster_Callback(hObject, eventdata, handles)
valid_clear_clustering;
Q = ['Clustering state: Clustering results are cleared !'];
set(handles.Outext1, 'String', Q);

function pushValidation_Callback(hObject, eventdata, handles)
valid_validation_check;

function pushClearValid_Callback(hObject, eventdata, handles)
valid_clear_validation;
Q = ['Validation state: ' 'Validity index values are cleared !'];
set(handles.Outext2, 'String', Q);

function pushClearPlot_Callback(hObject, eventdata, handles)
valid_clear_plotting;

function pushRedraw_Callback(hObject, eventdata, handles)
Q = handles.validstate;
if isempty(Q)
   return;
end
id = get(handles.Plotnewfigure, 'Checked');
if strcmp(id, 'on')
   figure('color','white');
else
   axes(handles.axes1);
   cla;
end
valid_redraw;

% ============================================ %
function popDataset_Callback(hObject, eventdata, handles)
Q = get(handles.Loadsolution, 'Checked');
if strcmp(Q, 'on')
  return;
end
id = get(handles.popDataset, 'Value');
dname = get(handles.popDataset, 'String');
if strcmp(dname{id}, 'iris.txt')
  set(handles.popMeasure, 'Value', 1);
  set(handles.Normaliz, 'Checked', 'off');
else
  set(handles.popMeasure, 'Value', 2);
end
file = dname{id};
data = load(['datafile\' file]);
R = ['Running state: Data file ' file ' is loaded.'];
valid_data_load;

function popDataset_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function popPosition_Callback(hObject, eventdata, handles)
id = get(handles.popDataset, 'Value');
dname = get(handles.popDataset, 'String');
demo = get(handles.Loademo, 'Checked');
if strcmp(demo, 'on')
   data = load(['datafile\' dname{id}]);
else
   data = load(dname);
end
valid_data_split

function popPosition_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
dname = {'Class labels in first column/row', 'Class labels in end column/row', ...
    'All columns/rows are data'};
set(hObject, 'String', dname);

function popAlgorithm_Callback(hObject, eventdata, handles)
function popAlgorithm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
dname = {'PAM', 'K-means', 'Hierarchical-single', 'Hierarchical-complete', ...
     'Hierarchical-centroid', 'SOM', 'Affinity Propagation', 'New algorithm'};
set(hObject, 'String', dname);

function popClutertimes_Callback(hObject, eventdata, handles)
function popClutertimes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
dname = {'Single Clustering', 'Multiple Clustering'};
set(hObject, 'String', dname);
set(hObject, 'Value', 2);

function popIndex_Callback(hObject, eventdata, handles)
function popIndex_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'));
    set(hObject,'BackgroundColor','white');
end
dname = {'Rand', 'Adjusted Rand', 'Jaccard', 'FM', 'Silhouette', ...
                'Davies-Bouldin', 'Calinski-Harabasz', 'Dunn', 'R-Squared', ...
                'Homogeneity-Separation', 'All validity indices'};
set(hObject, 'String', dname);

function popMeasure_Callback(hObject, eventdata, handles)
function popMeasure_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
dname = {'Euclidean distance', 'Pearson correlation'};
set(hObject, 'String', dname);
set(hObject, 'Value', 1);

% ============================================ %
function Outext1_Callback(hObject, eventdata, handles)
function Outext1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Outext2_Callback(hObject, eventdata, handles)
function Outext2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ============================================ %
function FileMenu_Callback(hObject, eventdata, handles)
function Loadata_Callback(hObject, eventdata, handles)
[file,filePath] = uigetfile('*.txt'); 
if isequal(file, 0)
   return;
end
dname = [filePath file];
try
  data = load(dname);
catch
  set(handles.Outext1, 'String', 'Running state: incorrect data file !');
  return
end
R = ['Running state: Data file ' file ' is loaded.'];
valid_data_load;

set(hObject, 'Checked', 'on');         % 'uimenu' class
set(handles.Loademo, 'Checked', 'off');
set(handles.Loadsolution, 'Checked', 'off');
set(handles.popDataset, 'Value', 1);
set(handles.popDataset, 'String', dname);
set(handles.popMeasure, 'Value', 1);
%set(handles.LoadSmatrix, 'Checked', 'off');

function Loademo_Callback(hObject, eventdata, handles)
dname = {'iris.txt', 'leuk72_3k.txt'};
file = dname{1}; 
data = load(['datafile\' file]);
R = ['Running state: Default data set ' file ' is loaded.'];
valid_data_load; % Go to pop-up menu "Data set" to load other data.

set(handles.popDataset, 'String', dname);
set(hObject, 'Checked', 'on'); 
set(handles.Loadata, 'Checked', 'off');
set(handles.Loadsolution, 'Checked', 'off');
set(handles.popDataset, 'Value', 1);
%set(handles.LoadSmatrix, 'Checked', 'off');

function Loadsolution_Callback(hObject, eventdata, handles)
valid_solution_load;

function Save_Callback(hObject, eventdata, handles)
valid_save_solution;

function CloseMenuItem_Callback(hObject, eventdata, handles)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end
delete(handles.figure1)

% --------------------------------------------------------------------
function Validation_Callback(hObject, eventdata, handles)
function Validone_Callback(hObject, eventdata, handles)
set(hObject, 'Checked', 'on'); 
set(handles.Validmutiple, 'Checked', 'off');

function Validmutiple_Callback(hObject, eventdata, handles)
valid_mutiple_validation;

function Validityanalys_Callback(hObject, eventdata, handles)
set(hObject, 'Checked', 'on'); 
set(handles.EstimateK, 'Checked', 'off');
dname = {'Rand', 'Adjusted Rand', 'Jaccard', 'FM', 'Silhouette', ...
                'Davies-Bouldin', 'Calinski-Harabasz', 'Dunn', 'R-Squared', ...
                 'Homogeneity-Separation', 'All validity indices'};
set(handles.popIndex, 'String', dname);
set(handles.popIndex, 'Value', 1);

function EstimateK_Callback(hObject, eventdata, handles)
set(hObject, 'Checked', 'on'); 
set(handles.Validityanalys, 'Checked', 'off');
dname = {'Rand', 'Adjusted Rand', 'Jaccard', 'FM', 'Silhouette', ...
                'Davies-Bouldin', 'Calinski-Harabasz', 'Dunn', 'R-Squared', ...
                'C-index', 'Krzanowski-Lai', 'Hartigan', 'RMSSTD group', ...
                'weighted inter/intra', 'Homogeneity-Separation',...
                '-  -  -  -  -  -  -  -  -  -', 'All validity indices'};
set(handles.popIndex, 'String', dname);
set(handles.popIndex, 'Value', 5);
set(handles.popClutertimes, 'Value', 2);
set(handles.Validone, 'Checked', 'on');
set(handles.Validmutiple, 'Checked', 'off');

function ClusterRow_Callback(hObject, eventdata, handles)
set(hObject, 'Checked', 'on'); 
set(handles.ClusterColumn, 'Checked', 'off');

function ClusterColumn_Callback(hObject, eventdata, handles)
set(hObject, 'Checked', 'on'); 
set(handles.ClusterRow, 'Checked', 'off');

% --------------------------------------------------------------------
function Options_Callback(hObject, eventdata, handles)
function Standz_Callback(hObject, eventdata, handles)
id = get(handles.Standz, 'Checked');
if strcmp(id, 'on')
   set(hObject, 'Checked', 'off');
else
   set(hObject, 'Checked', 'on');
end
set(handles.Normaliz, 'Checked', 'off');

function Normaliz_Callback(hObject, eventdata, handles)
id = get(handles.Normaliz, 'Checked');
if strcmp(id, 'on')
   set(hObject, 'Checked', 'off');
else
   set(hObject, 'Checked', 'on');
end
set(handles.Standz, 'Checked', 'off'); 

% --------------------------------------------------------------------
function View_Callback(hObject, eventdata, handles)
function Plotdata_Callback(hObject, eventdata, handles)
valid_data_plot_check;
if id_return
   return
end
cla;
labels = handles.truelabels;
data = handles.data;
id ='co';
plot_data_bylabels;

function PlotdataLabel_Callback(hObject, eventdata, handles)
valid_data_plot_check;
if id_return
   return
end
valid_data_plot;

function PlotHierarchy_Callback(hObject, eventdata, handles)
if isempty(handles.datastate)
   R = 'Running state: Please load data first for using this function !';
   set(handles.Outext1, 'String', R);
   return
end
valid_plot_hierarchytree;

function Plotnewfigure_Callback(hObject, eventdata, handles)
id = get(handles.Plotnewfigure, 'Checked');
if strcmp(id, 'on')
   set(hObject, 'Checked', 'off');
else
   set(hObject, 'Checked', 'on');
end

% --------------------------------------------------------------------
function Tool_Callback(hObject, eventdata, handles)
function ErrorRate_Callback(hObject, eventdata, handles)
valid_errorate_check;

function PrintLabel_Callback(hObject, eventdata, handles)
valid_print_label;

% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
function helpf_Callback(hObject, eventdata, handles)
winopen('Readme.txt'); 

function Copyright_Callback(hObject, eventdata, handles)
winopen('Copyright.txt'); 
