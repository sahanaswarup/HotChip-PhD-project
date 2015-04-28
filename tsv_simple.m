function varargout = tsv_simple(varargin)
% TSV_SIMPLE MATLAB code for tsv_simple.fig
%      TSV_SIMPLE, by itself, creates a new TSV_SIMPLE or raises the existing
%      singleton*.
%
%      H = TSV_SIMPLE returns the handle to a new TSV_SIMPLE or the handle to
%      the existing singleton*.
%
%      TSV_SIMPLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TSV_SIMPLE.M with the given input arguments.
%
%      TSV_SIMPLE('Property','Value',...) creates a new TSV_SIMPLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tsv_simple_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tsv_simple_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tsv_simple

% Last Modified by GUIDE v2.5 29-Apr-2014 13:32:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tsv_simple_OpeningFcn, ...
                   'gui_OutputFcn',  @tsv_simple_OutputFcn, ...
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


% --- Executes just before tsv_simple is made visible.
function tsv_simple_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tsv_simple (see VARARGIN)

% Choose default command line output for tsv_simple
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tsv_simple wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tsv_simple_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    x_y_dim       = str2double(get(handles.x_y_dim,'String'));
    num_layers    = str2double(get(handles.num_layers,'String'));
    z_dim         = str2double(get(handles.z_dim,'String'));
    lc            = 0.1; %str2double(get(handles.lc,'String'));
    tsv_r         = str2double(get(handles.tsv_r,'String'));
    tsv_t         = str2double(get(handles.tsv_thickness,'String'));
    num_cells     = floor(x_y_dim / (2 * (tsv_r + tsv_t)));
    d_value       = x_y_dim / num_cells;
    dz_value      = z_dim / 10;
    tsv_pad_t     = dz_value * 0.2;
    micro_bump_r  = dz_value * 0.4;
    tsv_option    = get(handles.tsv_placement_automatic,'Value');
    contents      = cellstr(get(handles.popupmenu2,'String'));
    tsv_file      = contents{get(handles.popupmenu2,'Value')};
    gmsh_filename = get(handles.gmsh_filename,'String');
    gmsh_path     = get(handles.gmsh_path,'String');
    contents      = cellstr(get(handles.power_loc_file,'String'));
    power_file    = contents{get(handles.power_loc_file,'Value')};
    gap_bw_layers = micro_bump_r*2 + tsv_pad_t;
    
    if tsv_option == 0
        try
            csvread(tsv_file);
        catch
            err_str = ['Manual TSV placement Chosen but invalid csv file specified.' ...
                       'Please select a valid csv file or choose automatic placement'];
            errordlg(err_str,'Invalid Input');
            return;
        end
    end
    if exist(gmsh_path,'file') ~= 2
        err_str = ['Incorrect Path to gmsh executable.' ...
                   'Launching gmsh will be skipped.' ...
                   'Please run the generated gmsh file "' gmsh_filename ...
                   '" manually'];
        errordlg(err_str, 'Invalid Input');
        return;
    end
    try
        csvread(power_file);
    catch
        err_str = ['Please select a valid file containing power locations'];
        errordlg(err_str,'Invalid Input');
        return;
    end
    
    tic
    generate_plot(x_y_dim,z_dim,num_layers, d_value, dz_value, ...
                  gap_bw_layers, tsv_option, tsv_file, tsv_r, tsv_t, ...
                  'tsv_simple_struct_fdd_input');
    [heat_map, min_temp, max_temp] = struct_fdd('tsv_simple_struct_fdd_input', power_file);
    generate_msh(x_y_dim,z_dim,num_layers, d_value, dz_value,  ...
             tsv_option, tsv_file, tsv_r,tsv_t, tsv_pad_t, ...
             micro_bump_r, gmsh_filename,lc, heat_map, min_temp, max_temp);
    enable_2d = get(handles.enable_2d_mesh,'Value');
    if enable_2d == 1
        gmsh_path1 = [gmsh_path ' -2 ' gmsh_filename];
        system(gmsh_path1);
        gmsh_path2 = [gmsh_path ' ' gmsh_filename ' custom.msh' '&'];
        system(gmsh_path2);
    else
        gmsh_path1 = [gmsh_path ' ' gmsh_filename '&'];
        system(gmsh_path1);
    end 
    total_time=toc


function x_y_dim_Callback(hObject, eventdata, handles)
% hObject    handle to x_y_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_y_dim as text
%        str2double(get(hObject,'String')) returns contents of x_y_dim as a double



% --- Executes during object creation, after setting all properties.
function x_y_dim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_y_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_tsv_y_Callback(hObject, eventdata, handles)
% hObject    handle to num_tsv_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_tsv_y as text
%        str2double(get(hObject,'String')) returns contents of num_tsv_y as a double


% --- Executes during object creation, after setting all properties.
function num_tsv_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_tsv_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tsv_r_Callback(hObject, eventdata, handles)
% hObject    handle to tsv_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tsv_r as text
%        str2double(get(hObject,'String')) returns contents of tsv_r as a double


% --- Executes during object creation, after setting all properties.
function tsv_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tsv_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_dim_Callback(hObject, eventdata, handles)
% hObject    handle to y_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_dim as text
%        str2double(get(hObject,'String')) returns contents of y_dim as a double


% --- Executes during object creation, after setting all properties.
function y_dim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_dim_Callback(hObject, eventdata, handles)
% hObject    handle to z_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_dim as text
%        str2double(get(hObject,'String')) returns contents of z_dim as a double


% --- Executes during object creation, after setting all properties.
function z_dim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_tsv_x_Callback(hObject, eventdata, handles)
% hObject    handle to num_tsv_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_tsv_x as text
%        str2double(get(hObject,'String')) returns contents of num_tsv_x as a double


% --- Executes during object creation, after setting all properties.
function num_tsv_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_tsv_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gmsh_filename_Callback(hObject, eventdata, handles)
% hObject    handle to gmsh_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gmsh_filename as text
%        str2double(get(hObject,'String')) returns contents of gmsh_filename as a double


% --- Executes during object creation, after setting all properties.
function gmsh_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gmsh_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gmsh_path_Callback(hObject, eventdata, handles)
% hObject    handle to gmsh_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gmsh_path as text
%        str2double(get(hObject,'String')) returns contents of gmsh_path as a double


% --- Executes during object creation, after setting all properties.
function gmsh_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gmsh_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tsv_placement_automatic.
function tsv_placement_automatic_Callback(hObject, eventdata, handles)
% hObject    handle to tsv_placement_automatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tsv_placement_automatic
set(handles.tsv_placement_manual, 'Value', 0);


% --- Executes on button press in tsv_placement_manual.
function tsv_placement_manual_Callback(hObject, eventdata, handles)
% hObject    handle to tsv_placement_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tsv_placement_manual
set(handles.tsv_placement_automatic, 'Value', 0);
x = set(handles.num_tsv_x, 'enable','off');


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
[filename, pathname] = uigetfile('*.csv', 'Pick a file containing TSV coordinates');
if filename ~= 0
    val = get(hObject, 'Value');
    set(hObject, 'Value', val);
    set(hObject, 'String', [pathname filename]);
end



% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_Callback(hObject, eventdata, handles)
% hObject    handle to lc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc as text
%        str2double(get(hObject,'String')) returns contents of lc as a double


% --- Executes during object creation, after setting all properties.
function lc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in generate_plot.
function generate_plot_Callback(hObject, eventdata, handles)
% hObject    handle to generate_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    x_dim = str2double(get(handles.x_y_dim,'String'));
    y_dim = str2double(get(handles.y_dim,'String'));
    z_dim = str2double(get(handles.z_dim,'String'));
    num_layers = str2double(get(handles.num_layers,'String'));
    num_tsv_x = str2double(get(handles.num_tsv_x,'String'));
    num_tsv_y = str2double(get(handles.num_tsv_y,'String'));
    tsv_r = str2double(get(handles.tsv_r,'String'));
    tsv_t = str2double(get(handles.tsv_thickness,'String'));
    tsv_option = get(handles.tsv_placement_automatic,'Value');
    contents = cellstr(get(handles.popupmenu2,'String'));
    tsv_file = contents{get(handles.popupmenu2,'Value')};
    gmsh_filename = get(handles.gmsh_filename,'String');
    if tsv_option == 0
        try
            csvread(tsv_file);
        catch
            err_str = ['Manual TSV placement Chosen but invalid csv file specified.' ...
                       'Please select a valid csv file or choose automatic placement'];
            errordlg(err_str,'Invalid Input');
            return;
        end
    end
    contents = cellstr(get(handles.power_loc_file,'String'));
    power_file = contents{get(handles.power_loc_file,'Value')};
    try
        csvread(power_file);
    catch
        err_str = ['Please select a valid file containing power locations'];
        errordlg(err_str,'Invalid Input');
        return;
    end
    tsv_pad_t = dz_value * 0.2;
    micro_bump_r = dz_value * 0.4;
    generate_plot(x_dim,y_dim,z_dim,num_layers, micro_bump_r*2 + tsv_pad_t, tsv_option, tsv_file, num_tsv_x,num_tsv_y,tsv_r, tsv_t, 'tsv_simple_struct_fdd_input');
    struct_fdd('tsv_simple_struct_fdd_input', power_file);



function num_layers_Callback(hObject, eventdata, handles)
% hObject    handle to num_layers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_layers as text
%        str2double(get(hObject,'String')) returns contents of num_layers as a double


% --- Executes during object creation, after setting all properties.
function num_layers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_layers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gap_bw_layers_Callback(hObject, eventdata, handles)
% hObject    handle to gap_bw_layers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gap_bw_layers as text
%        str2double(get(hObject,'String')) returns contents of gap_bw_layers as a double


% --- Executes during object creation, after setting all properties.
function gap_bw_layers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gap_bw_layers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tsv_thickness_Callback(hObject, eventdata, handles)
% hObject    handle to tsv_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tsv_thickness as text
%        str2double(get(hObject,'String')) returns contents of tsv_thickness as a double


% --- Executes during object creation, after setting all properties.
function tsv_thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tsv_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in enable_2d_mesh.
function enable_2d_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to enable_2d_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of enable_2d_mesh


% --- Executes on button press in enable_3d_mesh.
function enable_3d_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to enable_3d_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of enable_3d_mesh


% --- Executes on button press in enable_heat_sink.
function enable_heat_sink_Callback(hObject, eventdata, handles)
% hObject    handle to enable_heat_sink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of enable_heat_sink


% --- Executes on button press in enable_fr4_layer.
function enable_fr4_layer_Callback(hObject, eventdata, handles)
% hObject    handle to enable_fr4_layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of enable_fr4_layer


% --- Executes on selection change in power_loc_file.
function power_loc_file_Callback(hObject, eventdata, handles)
% hObject    handle to power_loc_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns power_loc_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from power_loc_file
[filename, pathname] = uigetfile('*.csv', 'Pick a file containing Power Locations');
if filename ~= 0
    val = get(hObject, 'Value');
    set(hObject, 'Value', val);
    set(hObject, 'String', [pathname filename]);
end


% --- Executes during object creation, after setting all properties.
function power_loc_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to power_loc_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
