function varargout = OneColorFIONA(varargin)
% ONECOLORFIONA MATLAB code for OneColorFIONA.fig
%      ONECOLORFIONA, by itself, creates a new ONECOLORFIONA or raises the existing
%      singleton*.
%
%      H = ONECOLORFIONA returns the handle to a new ONECOLORFIONA or the handle to
%      the existing singleton*.
%
%      ONECOLORFIONA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ONECOLORFIONA.M with the given input arguments.
%
%      ONECOLORFIONA('Property','Value',...) creates a new ONECOLORFIONA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OneColorFIONA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OneColorFIONA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Shorcuts
% s for Left
% f for Right
% e for Up
% d for Multiple
% c for Down
% r for Reorient
% z for Save

% Edit the above text to modify the response to help OneColorFIONA

% Last Modified by GUIDE v2.5 25-Mar-2015 14:43:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OneColorFIONA_OpeningFcn, ...
                   'gui_OutputFcn',  @OneColorFIONA_OutputFcn, ...
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

% --- Executes just before OneColorFIONA is made visible.
function OneColorFIONA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OneColorFIONA (see VARARGIN)
% UIWAIT makes OneColorFIONA wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Choose default command line output for OneColorFIONA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Add listener for keypress function
hListener=addlistener(handles.figure1,'WindowKeyPress', @keyPress);

% Set the folder directory to be current directory
directory = pwd;
set(handles.FileLocation,'String',directory);
set(handles.listbox1,'String',{''});
% Note, in order for selection of listbox to be disabled, Value needs to be
% set to an empty matrix, and the max needs to be 2 or more (max-min >1)
set(handles.listbox1,'Value',[]);
set(handles.listbox2,'String',{''});
set(handles.listbox2,'Value',[]);
set(handles.listbox3,'String',{''});
set(handles.listbox3,'Value',[]);
set(handles.listbox4,'String',{''});
set(handles.listbox4,'Value',[]);
setappdata(handles.SliderFrames,'Modes','PostProcess');

% Adjust settings for PostProcess mode
%PostProcessVal = get(handles.PostProcessing,'Value');
set(handles.SelectRegion,'Visible','off');
set(handles.CreateRegion,'Visible','off');
set(handles.DeletePoint,'Visible','off');
set(handles.listbox1,'Visible','off');
set(handles.DeleteFrame,'Visible','off');
%set(handles.TextName,'Visible','off');
%set(handles.FolderName,'Visible','off');
%set(handles.FileNameButton,'Visible','off');
set(handles.listbox2,'Visible','off');
%set(handles.Save,'Visible','off');
%set(handles.PreviewButton,'Visible','off');
%set(handles.DeleteCroppedFile,'Visible','off');
set(handles.SliderFrames,'Visible','off');
set(handles.TextFrames,'Visible','off');
%set(handles.CheckMapping,'Visible','off');

%set(handles.listbox1,'String',{'Drawing mode selected'});
set(handles.OpenFile,'Visible','off');
%set(handles.Screen,'Visible','off');
%set(handles.OpenReference,'Visible','off');
%set(handles.StartDrawing,'Visible','off');
setappdata(handles.SliderFrames,'Modes','PostProcess');
%set(handles.CreateROI,'Visible','off');
%set(handles.StartFrame,'Visible','off');
%set(handles.EndFrame,'Visible','off');
set(handles.FindPoints,'Visible','off');
set(handles.listbox3,'Visible','off');
set(handles.uipanel6,'Visible','off');
set(handles.Hist,'Visible','off');
set(handles.uipanel2,'Visible','off');
%set(handles.Modes,'Visible','off');
set(handles.StartFIONA,'Visible','off');
set(handles.BrowseFolder,'String','Code Path');
set(handles.BrowseFile,'String','Data Folder');
%set(handles.PostProcessing,'Visible','on');
set(handles.uipanel7,'Visible','on');
set(handles.SearchText,'Visible','on');
set(handles.listbox4,'Visible','on');
set(handles.OpenCleaned,'Visible','on');
set(handles.OpenTransformed,'Visible','on');
set(handles.SingleFile,'Visible','off');
set(handles.AllFiles,'Visible','off');
set(handles.FinalFrameNo,'Visible','off');
set(handles.TextFinalFrame,'Visible','off');
set(handles.DriftPanel,'Visible','off');

% --- Outputs from this function are returned to the command line.
function varargout = OneColorFIONA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in BrowseFolder.
function BrowseFolder_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
directory = uigetdir;
set(handles.FileLocation,'String',directory);

function FileLocation_Callback(hObject, eventdata, handles)
% hObject    handle to FileLocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of FileLocation as text
%        str2double(get(hObject,'String')) returns contents of FileLocation as a double

% --- Executes during object creation, after setting all properties.
function FileLocation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileLocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FileName_Callback(hObject, eventdata, handles)
% hObject    handle to FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of FileName as text
%        str2double(get(hObject,'String')) returns contents of FileName as a double

% --- Executes during object creation, after setting all properties.
function FileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in BrowseFile.
function BrowseFile_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Browse for file name and set the string to the FileName textbox

Tag = getappdata(handles.SliderFrames,'Modes');
%assignin('base','Tag',Tag);

switch Tag
    case 'FIONA'
        FullFolderPath = uigetdir;
        CodePath = get(handles.FileLocation,'string');
        name=strrep(FullFolderPath,[CodePath '\'],'');
        set(handles.FileName,'String',name);
    case 'PostProcess'
        FullFolderPath = uigetdir;
        CodePath = get(handles.FileLocation,'string');
        name=strrep(FullFolderPath,[CodePath '\'],'');
        set(handles.FileName,'String',name);
    case 'Split'
        name = uigetfile({'*.tif';'*.tiff';'*.avi'});
        set(handles.FileName,'String',name);
    case 'DriftCorrect'
        FullFolderPath = uigetdir;
        CodePath = get(handles.FileLocation,'string');
        name=strrep(FullFolderPath,[CodePath '\'],'');
        set(handles.FileName,'String',name);
    case 'Map'
        name = uigetfile({'*.tif';'*.tiff';'*.avi'});
        set(handles.FileName,'String',name);
    
    
end

% --- Executes on button press in OpenFile.
function OpenFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tag = getappdata(handles.SliderFrames,'Modes');

if strcmp(Tag,'DriftCorrect')
    % Get CodePath
    directory = get(handles.FileLocation,'string');
    
    % Get DataPath
    Folder = get(handles.FileName,'string');
    % Check if Folder is empty. if it not, check if the directory exists
    if isempty(Folder)
        % Ask for data path if Folder is empty
        DataPath = uigetdir;
    else
        DataPath = [directory '\' Folder];
        % If directory does not exist, get directory
        if ~exist(DataPath,'dir'); DataPath = uigetdir; end
    end
    name=strrep(DataPath,[directory '\'],'');
    set(handles.FileName,'String',name);
    cd(DataPath);
    
    % Find files
    FileType = '.tif';
    FileIn = dir(['*' FileType]);
    FileInput = cell(length(FileIn),1);
    for ind=1:length(FileIn)
        FileInput{ind}=FileIn(ind).name;
    end
    setappdata(handles.SliderFrames,'FileInput',FileInput);
    FileName = FileIn(1).name;
    FullPath = [directory '\' Folder '\' FileName];
    
    % Get Tiff info and Sum Factor
    info = imfinfo(FullPath);
    num_images = numel(info); 
    FinalFrameNo = str2double(get(handles.TextFinalFrame,'string'));
    SumFactor = floor(num_images/FinalFrameNo);
    
    % Sum images to FinalFrameNo
    dim = [info(1).Height info(1).Width FinalFrameNo];
    data = double(zeros(dim(1),dim(2),FinalFrameNo));
    data_adj = uint8(zeros(dim));
    for k = 0:FinalFrameNo-1
        for j = 1:SumFactor
            data(:,:,k+1) = data(:,:,k+1)+double(imread(FileInput{1}, (SumFactor*k)+j, 'Info', info));
        end

        if mod(k,1)==0
            set(handles.listbox2,'String',{['Summing images: ' num2str(k+1) ' / ' num2str(FinalFrameNo)]});
            drawnow();
        end
        %imwrite(FinalFrame, [FileInputName{ind} '(' num2str(SumFactor) 'x).tif'], 'WriteMode', 'append',  'Compression','none');
    end
    setappdata(handles.SliderFrames,'data_DriftBig',data);
    data = uint8(data/(SumFactor*256));
    setappdata(handles.SliderFrames,'data_Drift',data);
    data1 = data(:,:,1);
else
    directory = get(handles.FileLocation,'string');
    FileName = get(handles.FileName,'string');
    FileType = '.tif';
    FileInput = regexprep(FileName, FileType, '');
    FullPath = [directory '\' FileName];
    info = imfinfo(FullPath);
    num_images = numel(info);       
    dim = [info(1).Height info(1).Width num_images];
    data_adj = uint8(zeros(dim));

    % Read image
    data1=imread(FullPath, 1, 'Info', info);
end

% Save directory, FileInput, FullPath and info                  
setappdata(handles.SliderFrames,'directory',directory);
setappdata(handles.SliderFrames,'FileInput',FileInput);
setappdata(handles.SliderFrames,'FullPath',FullPath);
setappdata(handles.SliderFrames,'info',info);

% Read first frame
[counts,x]=imhist(data1);
[C,I]=max(counts);
range = max(x)-min(x);
x = (x-min(x))/range;
set(handles.SliderBrightMin,'UserData',[C;x;counts]);
BrightMinVal=ceil(x(I)*2*100)/100;
set(handles.SliderBrightMin,'Value',BrightMinVal);
set(handles.TextBrightMin,'String',BrightMinVal);
CumulativeSum=cumsum(counts);
BrightMaxVal=floor(0.99*sum(counts));
BrightMaxVal=floor(x((find(CumulativeSum>BrightMaxVal, 1 ))*2)*100)/100;
if BrightMinVal>=BrightMaxVal; BrightMaxVal=BrightMinVal+0.01; end;
set(handles.SliderBrightMax,'Value',BrightMaxVal);
set(handles.TextBrightMax,'String',BrightMaxVal);
axes(handles.Hist); hold off;
bar(x,counts,'blue'); hold on;
plot([0 1],[0 C],'black'); axis off;
set(handles.Hist,'yticklabel',[]);
xlim([0 1]); ylim([min(counts) max(counts)]);

if strcmp(Tag,'DriftCorrect')
    for i=1:dim(3)
        data_adj(:,:,i)=round(imadjust(data(:,:,i),[BrightMinVal BrightMaxVal],[]));
    end
    
    % Set the FileNo on SliderDrift and TextDrift to be 1
    set(handles.SliderDrift,'Value',1);
    set(handles.TextDrift,'String',1);

    % Set properties of SliderFileNo
    set(handles.SliderDrift,'Max',length(FileIn));
    set(handles.SliderDrift,'SliderStep',[1/(length(FileIn)-1) 1/(length(FileIn)-1)]);
else
    for i=1:dim(3)
        data = imread(FullPath, i, 'Info', info);
        data_adj(:,:,i)=round(imadjust(data,[BrightMinVal BrightMaxVal],[])*(255/65535));
        if mod(i,10)==0
            set(handles.listbox1,'String',{['Loading images: ' num2str(i) ' / ' num2str(dim(3))]});
            drawnow();
        end
    end
    set(handles.listbox1,'String',{['Loading images: ' num2str(dim(3)) ' / ' num2str(dim(3))]});
    drawnow();
end

axes(handles.Preview);
setappdata(handles.SliderFrames,'data_adj',data_adj);
imshow(data_adj(:,:,1));

% Set the value of slider.Frames and text.Frames to be 1
set(handles.SliderFrames,'Value',1);
set(handles.TextFrames,'String',1);

% Set properties of SliderFrames
set(handles.SliderFrames,'Max',dim(3));
set(handles.SliderFrames,'SliderStep',[1/(dim(3)-1) 1/(dim(3)-1)]);

% Initialize value of FrameNo
setappdata(handles.SliderFrames,'FrameNoPre',1);

% Initialize values of BrightMaxPre and BrightMinPre
setappdata(handles.SliderFrames,'BrightMaxPre',BrightMaxVal);
setappdata(handles.SliderFrames,'BrightMinPre',BrightMinVal);

%assignin('base','data_adj',data_adj);

% Initialize values of LimitPre
xlimit = floor(xlim);
ylimit = floor(ylim);
setappdata(handles.SliderFrames,'LimitPre',[xlimit ylimit]);

% Initialize the ImageSaved counter
ImageSaved = 0;
setappdata(handles.SliderFrames,'ImageSaved',ImageSaved);

% Initialize/ Reset values of FileName
setappdata(handles.SliderFrames,'FileName','');
setappdata(handles.SliderFrames,'ROI',[]);
setappdata(handles.SliderFrames,'ROIPos',[]);
setappdata(handles.SliderFrames,'HROIred',[]);

% Reset listboxes
%set(handles.listbox1,'String',{''});
%FolderName = getappdata(handles.SliderFrames,'FolderName');
%FileName = getappdata(handles.SliderFrames,'FileName');
%ImageSaved = getappdata(handles.SliderFrames,'ImageSaved');
%String1 = ['Folder Name   : ' FolderName];
%String2 = ['File Name        : ' FileName];
%String3 = ['Image Saved   : ' num2str(ImageSaved)];
%String4 = ['Files saved as : '];
%set(handles.listbox2,'String',{String1;String2;String3;String4});

% Clear data_Preview and data_PreviewAdj
%setappdata(handles.SliderFrames,'data_Preview',[]);
%setappdata(handles.SliderFrames,'data_PreviewAdj',[]);

cd(directory);

function TextBrightMax_Callback(hObject, eventdata, handles)
% hObject    handle to TextBrightMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TextBrightMax as text
%        str2double(get(hObject,'String')) returns contents of TextBrightMax as a double

BrightMaxVal = str2num(get(handles.TextBrightMax,'String'));
BrightMaxVal = floor(BrightMaxVal*100)/100;
BrightMinVal = get(handles.SliderBrightMin,'Value');
BrightMinVal = floor(BrightMinVal*100)/100;
FrameNo = floor(get(handles.SliderFrames,'Value'));
FullPath = getappdata(handles.SliderFrames,'FullPath');
info = imfinfo(FullPath);
data_Preview = getappdata(handles.SliderFrames,'data_Preview');
data_Drift = getappdata(handles.SliderFrames,'data_Drift');

% If data_Preview is empty, obtain data from usual file. If it is not
% empty, get data from data_Preview
if isempty(data_Preview)
    if isempty(data_Drift)
        data = imread(FullPath, FrameNo, 'Info', info);
    else
        data = data_Drift(:,:,FrameNo);
    end
else
    data = data_Preview(:,:,FrameNo);
end

% Issue warning dialog if BrightMaxVal is less than 0 or more than 1
if BrightMaxVal < 0
    h=warndlg('The minimum number is 0');
    uiwait(h);
    BrightMaxVal = 0;
elseif BrightMaxVal > 1
    h=warndlg('The maximum number is 1');
    uiwait(h);
    BrightMaxVal = 1;
end

if BrightMinVal>BrightMaxVal
    BrightMinVal=BrightMaxVal-0.01;
    if BrightMinVal < 0
        BrightMinVal = 0;
    end
    set(handles.TextBrightMin,'String',BrightMinVal);
    set(handles.SliderBrightMin,'Value',BrightMinVal);
end

% Set values of SliderBrigthMax and TextBrightMax
set(handles.SliderBrightMax,'Value',BrightMaxVal);
set(handles.TextBrightMax,'String',BrightMaxVal);

axes(handles.Preview);
xlimit = xlim; ylimit = ylim;
data_adj=imadjust(data,[BrightMinVal BrightMaxVal],[]);
imshow(data_adj);
xlim(xlimit); ylim(ylimit);

% Plot histogram
histData=get(handles.SliderBrightMin,'UserData');
datalength=(length(histData)-1)/2;
C=histData(1);
x=histData(2:(datalength+1));
counts=histData((datalength+2):end);
axes(handles.Hist);
hold off;
bar(x,counts,'blue');
xlim([0 1]);
ylim([min(counts) max(counts)]);
hold on;
plot([BrightMinVal BrightMaxVal],[0 C],'black');
set(handles.Hist,'yticklabel',[]);

% --- Executes during object creation, after setting all properties.
function TextBrightMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TextBrightMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TextBrightMin_Callback(hObject, eventdata, handles)
% hObject    handle to TextBrightMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TextBrightMin as text
%        str2double(get(hObject,'String')) returns contents of TextBrightMin as a double

BrightMaxVal = get(handles.SliderBrightMax,'Value');
BrightMaxVal = floor(BrightMaxVal*100)/100;
BrightMinVal = str2num(get(handles.TextBrightMin,'String'));
BrightMinVal = floor(BrightMinVal*100)/100;
FrameNo = floor(get(handles.SliderFrames,'Value'));
FullPath = getappdata(handles.SliderFrames,'FullPath');
info = imfinfo(FullPath);
data_Preview = getappdata(handles.SliderFrames,'data_Preview');
data_Drift = getappdata(handles.SliderFrames,'data_Drift');

% If data_Preview is empty, obtain data from usual file. If it is not
% empty, get data from data_Preview
if isempty(data_Preview)
    if isempty(data_Drift)
        data = imread(FullPath, FrameNo, 'Info', info);
    else
        data = data_Drift(:,:,FrameNo);
    end
else
    data = data_Preview(:,:,FrameNo);
end

% Issue warning dialog if BrightMinVal is less than 0 or more than 1
if BrightMinVal < 0
    h=warndlg('The minimum number is 0');
    uiwait(h);
    BrightMinVal = 0;
elseif BrightMinVal > 1
    h=warndlg('The maximum number is 1');
    uiwait(h);
    BrightMinVal = 1;
end

if BrightMinVal>BrightMaxVal
    BrightMaxVal=BrightMinVal+0.01;
    if BrightMaxVal > 1
        BrightMaxVal = 1;
    end
    set(handles.TextBrightMax,'String',BrightMaxVal);
    set(handles.SliderBrightMax,'Value',BrightMaxVal);
end

% Set values of SliderBrigthMin and TextBrightMin
set(handles.SliderBrightMin,'Value',BrightMinVal);
set(handles.TextBrightMin,'String',BrightMinVal);

axes(handles.Preview);
xlimit = xlim; ylimit = ylim;
data_adj=imadjust(data,[BrightMinVal BrightMaxVal],[]);
imshow(data_adj);
xlim(xlimit); ylim(ylimit);

% Plot histogram
histData=get(handles.SliderBrightMin,'UserData');
datalength=(length(histData)-1)/2;
C=histData(1);
x=histData(2:(datalength+1));
counts=histData((datalength+2):end);
axes(handles.Hist);
hold off;
bar(x,counts,'blue');
xlim([0 1]);
ylim([min(counts) max(counts)]);
hold on;
plot([BrightMinVal BrightMaxVal],[0 C],'black');
set(handles.Hist,'yticklabel',[]);

% --- Executes during object creation, after setting all properties.
function TextBrightMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TextBrightMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TextRad_Callback(~, eventdata, handles)
% hObject    handle to TextRad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TextRad as text
%        str2double(get(hObject,'String')) returns contents of TextRad as a double

Radius = get(handles.TextRad,'String');
Radius = floor(str2num(Radius));
set(handles.SliderRadius,'Value',Radius);
set(handles.TextRad,'String',Radius);

% --- Executes during object creation, after setting all properties.
function TextRad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TextRad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function SliderBrightMax_Callback(hObject, eventdata, handles)
% hObject    handle to SliderBrightMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

BrightMaxVal = get(handles.SliderBrightMax,'Value');
BrightMaxVal = floor(BrightMaxVal*100)/100;
BrightMinVal = get(handles.SliderBrightMin,'Value');
BrightMinVal = floor(BrightMinVal*100)/100;
FrameNo = floor(get(handles.SliderFrames,'Value'));
FullPath = getappdata(handles.SliderFrames,'FullPath');
info = imfinfo(FullPath);
data_Preview = getappdata(handles.SliderFrames,'data_Preview');
data_Drift = getappdata(handles.SliderFrames,'data_Drift');

% If data_Preview is empty, obtain data from usual file. If it is not
% empty, get data from data_Preview
if isempty(data_Preview)
    if isempty(data_Drift)
        data = imread(FullPath, FrameNo, 'Info', info);
    else
        data = data_Drift(:,:,FrameNo);
    end
else
    data = data_Preview(:,:,FrameNo);
end

if BrightMinVal>BrightMaxVal
    BrightMinVal=BrightMaxVal-0.01;
    if BrightMinVal < 0
        BrightMinVal = 0;
    end
    set(handles.TextBrightMin,'String',BrightMinVal);
    set(handles.SliderBrightMin,'Value',BrightMinVal);
end

set(handles.TextBrightMax,'String',BrightMaxVal);
axes(handles.Preview);
xlimit = xlim; ylimit = ylim;
data_adj=imadjust(data,[BrightMinVal BrightMaxVal],[]);
imshow(data_adj);
xlim(xlimit); ylim(ylimit);

% Plot histogram
histData=get(handles.SliderBrightMin,'UserData');
datalength=(length(histData)-1)/2;
C=histData(1);
x=histData(2:(datalength+1));
counts=histData((datalength+2):end);
axes(handles.Hist);
hold off;
bar(x,counts,'blue');
xlim([0 1]);
ylim([min(counts) max(counts)]);
hold on;
plot([BrightMinVal BrightMaxVal],[0 C],'black');
set(handles.Hist,'yticklabel',[]);

% --- Executes during object creation, after setting all properties.
function SliderBrightMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderBrightMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function SliderBrightMin_Callback(hObject, eventdata, handles)
% hObject    handle to SliderBrightMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

BrightMaxVal = get(handles.SliderBrightMax,'Value');
BrightMaxVal = floor(BrightMaxVal*100)/100;
BrightMinVal = get(handles.SliderBrightMin,'Value');
BrightMinVal = floor(BrightMinVal*100)/100;
FrameNo = floor(get(handles.SliderFrames,'Value'));
FullPath = getappdata(handles.SliderFrames,'FullPath');
info = imfinfo(FullPath);
data_Preview = getappdata(handles.SliderFrames,'data_Preview');
data_Drift = getappdata(handles.SliderFrames,'data_Drift');

% If data_Preview is empty, obtain data from usual file. If it is not
% empty, get data from data_Preview
if isempty(data_Preview)
    if isempty(data_Drift)
        data = imread(FullPath, FrameNo, 'Info', info);
    else
        data = data_Drift(:,:,FrameNo);
    end
else
    data = data_Preview(:,:,FrameNo);
end

if BrightMinVal>BrightMaxVal
    BrightMaxVal=BrightMinVal+0.01;
    if BrightMaxVal > 1
        BrightMaxVal = 1;
    end
    set(handles.TextBrightMax,'String',BrightMaxVal);
    set(handles.SliderBrightMax,'Value',BrightMaxVal);
end

set(handles.TextBrightMin,'String',BrightMinVal);
axes(handles.Preview);
xlimit = xlim; ylimit = ylim;
data_adj=imadjust(data,[BrightMinVal BrightMaxVal],[]);
imshow(data_adj);
xlim(xlimit); ylim(ylimit);

% Plot histogram
histData=get(handles.SliderBrightMin,'UserData');
datalength=(length(histData)-1)/2;
C=histData(1);
x=histData(2:(datalength+1));
counts=histData((datalength+2):end);
axes(handles.Hist);
hold off;
bar(x,counts,'blue');
xlim([0 1]);
ylim([min(counts) max(counts)]);
hold on;
plot([BrightMinVal BrightMaxVal],[0 C],'black');
set(handles.Hist,'yticklabel',[]);

% --- Executes during object creation, after setting all properties.
function SliderBrightMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderBrightMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function SliderRadius_Callback(hObject, eventdata, handles)
% hObject    handle to SliderRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

Radius = get(handles.SliderRadius,'Value');
Radius = floor(Radius);
set(handles.SliderRadius,'Value',Radius);
set(handles.TextRad,'String',Radius);

% --- Executes during object creation, after setting all properties.
function SliderRadius_CreateFcn(hObject, ~, handles)
% hObject    handle to SliderRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in Right.
function Right_Callback(hObject, eventdata, handles)
% hObject    handle to Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');

axes(handles.Preview);
dcm_obj=datacursormode;
Cursor=getCursorInfo(dcm_obj);
NoOfPoints = size(dataInput,1);

x=Cursor(1).Position(1);

j=1;
for l = 1:NoOfPoints
    if dataInput(j,1) >= x
        dataInput(j,:) = [];
        j=j-1;
    end
    j=j+1;
end

axes(handles.Preview); hold off;
plot(dataInput(:,1),dataInput(:,2),'g'); hold on;
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10);
hold off;
datacursormode on;    

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',dataInput);

% --- Executes on button press in Left.
function Left_Callback(hObject, eventdata, handles)
% hObject    handle to Left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');

axes(handles.Preview);
dcm_obj=datacursormode;
Cursor=getCursorInfo(dcm_obj);
NoOfPoints = size(dataInput,1);

x=Cursor(1).Position(1);

j=1;
for l = 1:NoOfPoints
    if dataInput(j,1) <= x
        dataInput(j,:) = [];
        j=j-1;
    end
    j=j+1;
end

axes(handles.Preview); hold off;
plot(dataInput(:,1),dataInput(:,2),'g'); hold on;
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10);
hold off;
datacursormode on;    

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',dataInput);

% --- Executes on button press in Up.
function Up_Callback(hObject, eventdata, handles)
% hObject    handle to Up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');

axes(handles.Preview);
dcm_obj=datacursormode;
Cursor=getCursorInfo(dcm_obj);
NoOfPoints = size(dataInput,1);

y=Cursor(1).Position(2);

j=1;
for l = 1:NoOfPoints
    if dataInput(j,2) >= y
        dataInput(j,:) = [];
        j=j-1;
    end
    j=j+1;
end

axes(handles.Preview); hold off;
plot(dataInput(:,1),dataInput(:,2),'g'); hold on; 
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
hold off;
datacursormode on;    

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',dataInput);

% --- Executes on button press in Down.
function Down_Callback(hObject, eventdata, handles)
% hObject    handle to Down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');

axes(handles.Preview);
dcm_obj=datacursormode;
Cursor=getCursorInfo(dcm_obj);
NoOfPoints = size(dataInput,1);

y=Cursor(1).Position(2);

j=1;
for l = 1:NoOfPoints
    if dataInput(j,2) <= y
        dataInput(j,:) = [];
        j=j-1;
    end
    j=j+1;
end

axes(handles.Preview); hold off;
plot(dataInput(:,1),dataInput(:,2),'g'); hold on;
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10);
hold off;
datacursormode on;    

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',dataInput);

% --- Executes on button press in Multiple.
function Multiple_Callback(hObject, eventdata, handles)
% hObject    handle to Multiple (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');

axes(handles.Preview);
dcm_obj=datacursormode;
Cursor=getCursorInfo(dcm_obj);
%assignin('base','Cursor',Cursor);
%assignin('base','dataInput',dataInput);
NoOfPoints = size(dataInput,1);  

for j=1:numel(Cursor);
    indexXGreen = dataInput(:,1)==Cursor(j).Position(1);
    indexYGreen = dataInput(:,2)==Cursor(j).Position(2);
    dataInput(indexXGreen & indexYGreen == 1,:)=[];
end

axes(handles.Preview); hold off;
plot(dataInput(:,1),dataInput(:,2),'g'); hold on; 
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
hold off;
datacursormode on; 

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',dataInput);

% --- Executes on slider movement.
function SliderFileNo_Callback(hObject, eventdata, handles)
% hObject    handle to SliderFileNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Get CodePath and DataPath
CodePath = getappdata(handles.SliderFrames,'CodePath');
DataPath = getappdata(handles.SliderFrames,'DataPath');
cd(DataPath);

% Delete existing figures
h=getappdata(handles.SliderFrames,'h');
try close(h); catch; end

% Get the FileNo from SliderFileNo
FileIn = dir('*.txt');
if length(FileIn) > 1
    FileNo = floor(get(handles.SliderFileNo,'Value'));
else
    FileNo = 1;
end

% Adjust the value of TextFileNo
set(handles.TextFileNo,'String',FileNo);

% Get FileInput and FileInputName
FileInput=getappdata(handles.SliderFrames,'FileInput');

% Open and plot the first file
fid=fopen(FileInput{FileNo});
OpenCleanedFile = get(handles.OpenCleaned,'Value');
OpenTransformedFile = get(handles.OpenTransformed,'Value');
if OpenCleanedFile || OpenTransformedFile ; 
    Input = textscan(fid,'%f%f%f','CommentStyle','##'); 
else
    Input = textscan(fid,'%*f%*f%*f%f%f%*f%*f','CommentStyle','##');
end
Input = [Input{1} Input{2} (1:length(Input{1}))'];
axes(handles.Preview); hold off;
plot(Input(:,1),Input(:,2),'g'); hold on;
plot(Input(1,1),Input(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(Input(end,1),Input(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
hold off;
datacursormode on;

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',Input);

% Switch to CodePath
cd(CodePath);

% --- Executes during object creation, after setting all properties.
function SliderFileNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderFileNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function TextFileNo_Callback(hObject, eventdata, handles)
% hObject    handle to TextFileNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TextFileNo as text
%        str2double(get(hObject,'String')) returns contents of TextFileNo as a double

% Get CodePath and DataPath
CodePath = getappdata(handles.SliderFrames,'CodePath');
DataPath = getappdata(handles.SliderFrames,'DataPath');
cd(DataPath);

% Get FileInput and FileInputName
FileInput=getappdata(handles.SliderFrames,'FileInput');
FileInputName=getappdata(handles.SliderFrames,'FileInputName');
num_images = length(FileInput);

% Get the FileNo from TextFileNo
FileNo = floor(str2num(get(handles.TextFileNo,'String')));

% Issue warning dialog if value of TextfileNo is less than 1 or more than
% the number of files
if FileNo < 1
    h=warndlg('The minimum number is 1');
    uiwait(h);
    FileNo = 1;
elseif FileNo > num_images
    h=warndlg(['The maximum number is ' num2str(num_images)]);
    uiwait(h);
    FileNo = num_images;
end

% Adjust the values of TextFileNo and SliderFileNo
set(handles.TextFileNo,'String',FileNo);
set(handles.SliderFileNo,'Value',FileNo);

% Open and plot the first file
fid=fopen(FileInput{FileNo});
Input = textscan(fid,'%f%f%f%f','CommentStyle','##');
Input = [Input{1} Input{2} Input{3} Input{4}];
axes(handles.Preview); hold off;
plot(Input(:,1),Input(:,2),'g'); hold on;
plot(Input(1,1),Input(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(Input(end,1),Input(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
plot(Input(:,3),Input(:,4),'r'); hold off;
hold off;
datacursormode on;

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',Input);

% Switch to CodePath
cd(CodePath);

% --- Executes during object creation, after setting all properties.
function TextFileNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TextFileNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in DeleteFile.
function DeleteFile_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get CodePath, DataPath and FileInput
CodePath = getappdata(handles.SliderFrames,'CodePath');
DataPath = getappdata(handles.SliderFrames,'DataPath');
FolderName=strrep(DataPath,[CodePath '\'],'');
cd(DataPath);

% Get FileInput and FileInputName
FileInput=getappdata(handles.SliderFrames,'FileInput');
FileInputName=getappdata(handles.SliderFrames,'FileInputName');

% Get the FileNo from SliderFileNo
FileNo = floor(get(handles.SliderFileNo,'Value'));

% Get OpenCleanedFile and OpenTransformedFile
OpenCleanedFile = get(handles.OpenCleaned,'Value');
OpenTransformedFile = get(handles.OpenTransformed,'Value');

% Delete the file associated with the FileInput
if exist(FileInput{FileNo},'file')
    button1 = questdlg(['Are you sure you want to delete "' FileInput{FileNo} '"'], 'Alert', 'No');
    switch button1
        case 'Yes'
            fclose('all');
            delete(FileInput{FileNo});
            if OpenCleanedFile == 1 || OpenTransformedFile == 1
                JPGname=[FileInputName{FileNo} '.jpg'];
                if exist(JPGname,'file'); delete(JPGname); end;
            end
            
            % Delete FileInput and FileInputName for the particular FileNo
            FileInput(FileNo) = [];
            FileInputName(FileNo) = [];
            FileInputString = cell(length(FileInput),1);

            % Get the number of images left
            num_images = length(FileInput);

            % Open and plot the next file or the previous file
            if FileNo <= num_images
                fid=fopen(FileInput{FileNo});
            else
                FileNo = num_images;
                fid=fopen(FileInput{FileNo});
            end
            if OpenTransformedFile || OpenCleanedFile; 
                Input = textscan(fid,'%f%f%f','CommentStyle','##');
            else
                Input = textscan(fid,'%*f%*f%*f%f%f%*f%*f','CommentStyle','##');
            end
            Input = [Input{1} Input{2} (1:length(Input{1}))'];
            axes(handles.Preview); hold off;
            plot(Input(:,1),Input(:,2),'g'); hold on;
            plot(Input(1,1),Input(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
            plot(Input(end,1),Input(end,2),'o','MarkerFaceColor','black','MarkerSize',10);
            hold off;
            datacursormode on;

            % Store the FileInput and FileInputName
            setappdata(handles.SliderFrames,'FileInput',FileInput);
            setappdata(handles.SliderFrames,'FileInputName',FileInputName);

            % Set the FileNo on SliderFileNo and TextFileNo to be 1
            set(handles.SliderFileNo,'Value',FileNo);
            set(handles.TextFileNo,'String',FileNo);

            % Set properties of SliderFileNo
            if num_images > 1
                set(handles.SliderFileNo,'Max',num_images);
                set(handles.SliderFileNo,'SliderStep',[1/(num_images-1) 1/(num_images-1)]);
            end
            
            % Store dataInput
            setappdata(handles.SliderFrames,'dataInput',Input);

            % Display the number of files and file name on listbox
            if OpenCleanedFile == 1 || OpenTransformedFile == 1
                num_images = length(FileInput);
                String1 = ['No of files : ' num2str(num_images)];
                set(handles.listbox4,'String',[String1; ' '; FileInput]);
                set(handles.listbox4,'Value',[]);
            else
                num_images = length(FileInput);
                for i = 1:num_images
                    if exist(['Cleaned\' FileInputName{i} '-Cleaned.txt'],'file') == 2
                        FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"     (saved)'];
                    else
                        FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"'];
                    end        
                    if exist(['Transformed\' FileInputName{i} '-Transformed.txt'],'file') == 2
                        FileInputString{i} = [FileInputString{i} '     (transformed)'];
                    end        
                end
                String1 = ['No of files : ' num2str(num_images)];
                set(handles.listbox4,'String',[String1; ' '; FileInputString]);
                set(handles.listbox4,'Value',[]);
            end
    end
else
    h=warndlg(['"' FileInput{FileNo} '" does not exist in folder "' FolderName '" . Try another folder or file name.']);
    uiwait(h);
end  

% Switch to CodePath
cd(CodePath);

% --- Executes on button press in ScreenOutlier.
function ScreenOutlier_Callback(hObject, eventdata, handles)
% hObject    handle to ScreenOutlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in Reorient.
function Reorient_Callback(hObject, eventdata, handles)
% hObject    handle to Reorient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get CodePath, DataPath and FileInput
CodePath = getappdata(handles.SliderFrames,'CodePath');
DataPath = getappdata(handles.SliderFrames,'DataPath');
cd(DataPath);

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');

% Get OpenCleanedFile and OpenTransformedFile
OpenCleanedFile = get(handles.OpenCleaned,'Value');
OpenTransformedFile = get(handles.OpenTransformed,'Value');

% Get FileInput and FileInputName
FileInputName=getappdata(handles.SliderFrames,'FileInputName');

% Get the FileNo from SliderFileNo
FileNo = floor(get(handles.SliderFileNo,'Value'));

% Open cleaned data and assign data to xFin and yFin
xInput=dataInput(:,1); yInput=dataInput(:,2); tInput=dataInput(:,3);

% Plot raw trace
xline=[(min(xInput)-((max(xInput)-min(xInput))/10)) max(xInput)];
dcm_obj = datacursormode;
Cursor = getCursorInfo(dcm_obj);
%assignin('base','Cursor',Cursor);
if length(Cursor)==2
    xPoints = [Cursor(1).Position(1) Cursor(2).Position(1)];
    yPoints = [Cursor(1).Position(2) Cursor(2).Position(2)];
    %fit=polyfit(xInput,yInput,1);
    gradient = (yPoints(2)-yPoints(1))/(xPoints(2)-xPoints(1));
    intercept = yPoints(1)-gradient*xPoints(1);
    yline=gradient*xline+intercept;
    axes(handles.Preview); hold off;
    plot(xInput,yInput,'g'); hold on;
    plot(xline,yline,'red','LineWidth',2);
    title(['Raw Trace '  FileInputName{FileNo}],'FontName','Palatino Linotype');
    xlabel('x (nm)','FontName','Palatino Linotype');
    ylabel('y (nm)','FontName','Palatino Linotype');
    hold off;

    % Find out rotated trace
    x=xInput-xline(1);
    y=yInput-yline(1);
    radius=x.*x+y.*y;
    theta = atan(gradient);
    tanNew = tan(atan(y./x)-theta);
    xNew = sqrt(radius./(1+tanNew.*tanNew));
    yNew = tanNew.*xNew;

    % Invert negative sloped traces
    if xNew(end)-xNew(1)<0
        yNew = flipud(yNew);
        xNew = flipud(xNew);
        tInput = abs(flipud(tInput-(max(tInput)+1)));
    end
    
    % Plot rotated trace
    axes(handles.Preview); hold off;
    plot(tInput, xNew);
    title('On-Axis Trace','FontName','Palatino Linotype');
    xlabel('Frame','FontName','Palatino Linotype');
    ylabel('On-Axis Distance (nm)','FontName','Palatino Linotype');
    hold off;
    
    % Save rotated trace variable
    setappdata(handles.SliderFrames,'OnAxisTrace',[xNew tInput]);
    
    h = figure;
    plot(yNew,xNew,'g'); hold on;
    plot([0 0],[min(xNew),max(xNew)],'red','LineWidth',2);
    title(['Rotated Trace Perpendicular to Microtubule ' FileInputName{FileNo}],'FontName','Palatino Linotype');
    xlabel('x (nm)','FontName','Palatino Linotype');
    ylabel('y (nm)','FontName','Palatino Linotype');
    setappdata(handles.SliderFrames,'h',h);

    % Find files
    FileType = '.txt';
    FileIn = dir(['*' FileType]);
    FileInput = cell(length(FileIn),1);
    FileInputName = cell(length(FileIn),1);
    for ind=1:length(FileIn)
        FileInput{ind}=FileIn(ind).name;
        FileInputName{ind}=strrep(FileIn(ind).name,FileType,'');
    end
    FileInputString = cell(length(FileInput),1);

    if OpenCleanedFile == 1
        DataPath = strrep(DataPath,'\Cleaned','');
        FileName = strrep(FileInputName{FileNo},'-Cleaned','');
        cd(DataPath);
        %assignin('base','DataPath',DataPath);
        % Save picture of transformed trace
        print('-djpeg','-r300',['Transformed\' FileName '-TransformedLin']);

        % Save data
        dlmwrite(['Transformed\' FileName '-TransformedLin.txt'],[xNew yNew tInput],'\t');

        % Save data for t-test input
        dlmwrite(['T-test\' FileName 'TtestLin.txt'],[xNew tInput],'\t');
    elseif OpenTransformedFile == 1
        DataPath = strrep(DataPath,'\Transformed','');
        FileName = strrep(FileInputName{FileNo},'-Transformed','');
        cd(DataPath);
        %assignin('base','DataPath',DataPath);
        % Save picture of transformed trace
        print('-djpeg','-r300',['Transformed\' FileName '-TransformedLin']);

        % Save data
        dlmwrite(['Transformed\' FileName '-TransformedLin.txt'],[xNew yNew tInput],'\t');

        % Save data for t-test input
        dlmwrite(['T-test\' FileName 'TtestLin.txt'],[xNew tInput],'\t');
    else
        % Save picture of transformed trace
        print('-djpeg','-r300',['Transformed\' FileInputName{FileNo} '-TransformedLin']);

        % Save data
        dlmwrite(['Transformed\' FileInputName{FileNo} '-TransformedLin.txt'],[xNew yNew tInput],'\t');

        % Save data for t-test input
        dlmwrite(['T-test\' FileInputName{FileNo} 'TtestLin.txt'],[xNew tInput],'\t');

        % Update listbox
        num_images = length(FileInput);
        for i = 1:num_images
            if exist(['Cleaned\' FileInputName{i} '-Cleaned.txt'],'file') == 2
                FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"     (saved)'];
            else
                FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"'];
            end
            if exist(['Transformed\' FileInputName{i} '-Transformed.txt'],'file') == 2
                FileInputString{i} = [FileInputString{i} '     (transformed)'];
            end     
        end
        String1 = ['No of files : ' num2str(num_images)];
        set(handles.listbox4,'String',[String1; ' '; FileInputString]);
    end
else
    h=warndlg('You need to select two datapoints to reorient. Use data cursor tab to select a point. Use alt to select a second point. ');
    uiwait(h);
end

% Switch to CodePath
cd(CodePath);

axes(handles.Preview);
datacursormode on;

% --- Executes on button press in DeleteBefore.
function DeleteBefore_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteBefore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');

axes(handles.Preview);
dcm_obj=datacursormode;
Cursor=getCursorInfo(dcm_obj);

StartPoint=Cursor(1).DataIndex;

dataInput(1:StartPoint,:)=[];

axes(handles.Preview); hold off;
plot(dataInput(:,1),dataInput(:,2),'g'); hold on; 
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
hold off;
datacursormode on;    

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',dataInput);

% --- Executes on button press in DeleteAfter.
function DeleteAfter_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteAfter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');

axes(handles.Preview);
dcm_obj=datacursormode;
Cursor=getCursorInfo(dcm_obj);

StartPoint=Cursor(1).DataIndex;

dataInput(StartPoint:end,:)=[];

axes(handles.Preview); hold off;
plot(dataInput(:,1),dataInput(:,2),'g'); hold on;
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
hold off;
datacursormode on;    

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',dataInput);

% --- Executes on button press in SaveText.
function SaveText_Callback(hObject, eventdata, handles)
% hObject    handle to SaveText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get CodePath, DataPath and FileInput
CodePath = getappdata(handles.SliderFrames,'CodePath');
DataPath = getappdata(handles.SliderFrames,'DataPath');
cd(DataPath);

% Delete existing figures outside GUI
h=getappdata(handles.SliderFrames,'h');
try close(h); catch; end

% Get FileInput and FileInputName
FileInput=getappdata(handles.SliderFrames,'FileInput');
FileInputName=getappdata(handles.SliderFrames,'FileInputName');
FileInputString = cell(length(FileInput),1);
%assignin('base','FileInput',FileInput);
%assignin('base','FileInputName',FileInputName);

% Get the FileNo from SliderFileNo
FileNo = floor(get(handles.SliderFileNo,'Value'));

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');

% Get OpenCleanedFile and OpenTransformedFile
OpenCleanedFile = get(handles.OpenCleaned,'Value');
OpenTransformedFile = get(handles.OpenTransformed,'Value');

% Save cleaned-up data in the 'Cleaned' folder
if OpenCleanedFile == 1
    dlmwrite(FileInputName{FileNo},dataInput,'\t');
    
    % Display the number of files and file name on listbox
    num_images = length(FileInput);
    String1 = ['No of files : ' num2str(num_images) '         "' FileInputName{FileNo} '  saved.'];
    set(handles.listbox4,'String',[String1; ' '; FileInput]);
    
    % Save picture of raw and fitted trace
    print('-djpeg','-r300',FileInputName{FileNo});
elseif OpenTransformedFile == 1
    h=warndlg('Use "Reorient" button to save transformed images');
    uiwait(h);
else
    dlmwrite(['Cleaned\' FileInputName{FileNo} '-Cleaned.txt'],dataInput,'\t');
    
    % Display the number of files and file name on listbox
    num_images = length(FileInput);
    for i = 1:num_images
        if exist(['Cleaned\' FileInputName{i} '-Cleaned.txt'],'file') == 2
            FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"     (saved)'];
        else
            FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"'];
        end
        if exist(['Transformed\' FileInputName{i} '-Transformed.txt'],'file') == 2
            FileInputString{i} = [FileInputString{i} '     (transformed)'];
        end
    end
    String1 = ['No of files : ' num2str(num_images) '         "' [FileInputName{FileNo} '-Cleaned.txt"'] '  saved.'];
    set(handles.listbox4,'String',[String1; ' '; FileInputString]);
    
    % Save picture of raw and fitted trace
    print('-djpeg','-r300',['Cleaned\' FileInputName{FileNo} '-Cleaned']);
end

% Switch to CodePath
cd(CodePath);

% --- Executes on slider movement.
function SliderFrames_Callback(hObject, eventdata, handles)
% hObject    handle to SliderFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if 0
    BrightMaxVal = get(handles.SliderBrightMax,'Value');
    BrightMaxVal = floor(BrightMaxVal*100)/100;
    BrightMinVal = get(handles.SliderBrightMin,'Value');
    BrightMinVal = floor(BrightMinVal*100)/100;

    BrightMaxPre = getappdata(handles.SliderFrames,'BrightMaxPre');
    BrightMinPre = getappdata(handles.SliderFrames,'BrightMinPre');
    %assignin('base','BrightMinPre',BrightMinPre);
    %assignin('base','BrightMaxPre',BrightMaxPre);

    if BrightMaxVal~=BrightMaxPre || BrightMinVal~=BrightMinPre
        FullPath = getappdata(handles.SliderFrames,'FullPath');
        info = imfinfo(FullPath);
        dim = [info(1).Height info(1).Width numel(info)];
        data_adj = uint8(zeros(dim));
        for i = 1:dim(3)
            data = imread(FullPath, i, 'Info', info);
            data_adj(:,:,i)=round(imadjust(data,[BrightMinVal BrightMaxVal],[])*(255/65535));
            %data_adj(:,:,i)=imadjust(data(:,:,i),[BrightMinVal BrightMaxVal],[]);
            if mod(i,10)==0
                set(handles.listbox1,'String',{['Adjusting brightness: ' num2str(i) ' / ' num2str(dim(3))]});
                drawnow();
            end
        end
        set(handles.listbox1,'String',{['Adjusting brightness: ' num2str(dim(3)) ' / ' num2str(dim(3))]});
        drawnow();
        setappdata(handles.SliderFrames,'data_adj',data_adj);
        FrameNo = floor(get(handles.SliderFrames,'Value'));
        set(handles.TextFrames,'String',FrameNo);
        set(handles.SliderFrames,'Value',FrameNo);
        axes(handles.Preview);
        xlimit = xlim; ylimit = ylim;
        imshow(data_adj(:,:,FrameNo));
        xlim(xlimit); ylim(ylimit);
        setappdata(handles.SliderFrames,'BrightMaxPre',BrightMaxVal);
        setappdata(handles.SliderFrames,'BrightMinPre',BrightMinVal);
        %assignin('base','data_adj',data_adj);
    end
end

hListener=addlistener(handles.SliderFrames,'ContinuousValueChange', @FrameChange);

% --- Executes during object creation, after setting all properties.
function SliderFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function TextFrames_Callback(hObject, eventdata, handles)
% hObject    handle to TextFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TextFrames as text
%        str2double(get(hObject,'String')) returns contents of TextFrames as a double

FrameNo = get(handles.TextFrames,'String');
FrameNo = floor(str2num(FrameNo));

ROI=getappdata(handles.SliderFrames,'ROI');
    
% Get data_Preview
data_Preview = getappdata(handles.SliderFrames,'data_Preview');
BrightMaxVal = get(handles.SliderBrightMax,'Value');
BrightMaxVal = floor(BrightMaxVal*100)/100;
BrightMinVal = get(handles.SliderBrightMin,'Value');
BrightMinVal = floor(BrightMinVal*100)/100;

BrightMaxPre = getappdata(handles.SliderFrames,'BrightMaxPre');
BrightMinPre = getappdata(handles.SliderFrames,'BrightMinPre');

if isempty(data_Preview)
    if BrightMaxVal~=BrightMaxPre || BrightMinVal~=BrightMinPre
        FullPath = getappdata(handles.SliderFrames,'FullPath');
        info = imfinfo(FullPath);
        dim = [info(1).Height info(1).Width numel(info)];
        data_adj = uint8(zeros(dim));
        for i = 1:dim(3)
            data = imread(FullPath, i, 'Info', info);
            data_adj(:,:,i)=round(imadjust(data,[BrightMinVal BrightMaxVal],[])*(255/65535));
            %data_adj(:,:,i)=imadjust(data(:,:,i),[BrightMinVal BrightMaxVal],[]);
            if mod(i,10)==0
                set(handles.listbox1,'String',{['Adjusting brightness: ' num2str(i) ' / ' num2str(dim(3))]});
                drawnow();
            end
        end
        set(handles.listbox1,'String',{['Adjusting brightness: ' num2str(dim(3)) ' / ' num2str(dim(3))]});
        drawnow();
        setappdata(handles.SliderFrames,'data_adj',data_adj);
        setappdata(handles.SliderFrames,'BrightMaxPre',BrightMaxVal);
        setappdata(handles.SliderFrames,'BrightMinPre',BrightMinVal);
    end

    %data_adj=imadjust(data,[BrightMinVal BrightMaxVal],[]);
    axes(handles.Preview);
    xlimit = xlim; ylimit = ylim;
    data_adj=getappdata(handles.SliderFrames,'data_adj');
    num_images = size(data_adj,3);
    
    % Issue warning dialog if value of FrameNo is less than 1 or more than
    % the number of files
    if FrameNo < 1
        h=warndlg('The minimum number is 1');
        uiwait(h);
        FrameNo = 1;
    elseif FrameNo > num_images
        h=warndlg(['The maximum number is ' num2str(num_images)]);
        uiwait(h);
        FrameNo = num_images;
    end
    
    data_adj=data_adj(:,:,FrameNo);
    if ~isempty(ROI); data_adj(ROI)=255; end;
    imshow(data_adj);
    xlim(xlimit); ylim(ylimit);
    setappdata(handles.SliderFrames,'FrameNoPre',FrameNo);
    set(handles.TextFrames,'String',FrameNo);
    set(handles.SliderFrames,'Value',FrameNo);
else
    if BrightMaxVal~=BrightMaxPre || BrightMinVal~=BrightMinPre
        data_Preview = getappdata(handles.SliderFrames,'data_Preview');
        dim = size(data_Preview);
        data_adj = uint8(zeros(dim));
        for i = 1:dim(3)
            data_adj(:,:,i) = round(imadjust(data_Preview(:,:,i),[BrightMinVal BrightMaxVal],[])*(255/65535));
            %data_adj(:,:,i)=imadjust(data(:,:,i),[BrightMinVal BrightMaxVal],[]);
        end
        setappdata(handles.SliderFrames,'data_PreviewAdj',data_adj);
        setappdata(handles.SliderFrames,'BrightMaxPre',BrightMaxVal);
        setappdata(handles.SliderFrames,'BrightMinPre',BrightMinVal);
    end

    %data_adj=imadjust(data,[BrightMinVal BrightMaxVal],[]);
    axes(handles.Preview);
    xlimit = xlim; ylimit = ylim;
    data_adj=getappdata(handles.SliderFrames,'data_PreviewAdj');
    num_images = size(data_adj,3);
    
    % Issue warning dialog if value of FrameNo is less than 1 or more than
    % the number of files
    if FrameNo < 1
        h=warndlg('The minimum number is 1');
        uiwait(h);
        FrameNo = 1;
    elseif FrameNo > num_images
        h=warndlg(['The maximum number is ' num2str(num_images)]);
        uiwait(h);
        FrameNo = num_images;
    end
    
    data_adj=data_adj(:,:,FrameNo);
    if ~isempty(ROI); data_adj(ROI)=255; end;
    imshow(data_adj);
    xlim(xlimit); ylim(ylimit);
    setappdata(handles.SliderFrames,'FrameNoPre',FrameNo);
    set(handles.TextFrames,'String',FrameNo);
    set(handles.SliderFrames,'Value',FrameNo);
end

% Store values of SliderFrames and TextFrames
set(handles.SliderFrames,'Value',FrameNo);
set(handles.TextFrames,'String',FrameNo);

% --- Executes during object creation, after setting all properties.
function TextFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TextFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FrameChange(hObject,eventdata,handles)

if ~(exist('handles','var'))
     handles=guidata(hObject);
end

% To make sure that the Preview figure is updated, we need to tweak the GUI
% option a bit (from:
% http://www.mathworks.com/support/solutions/en/data/1-192W8/index.html)
% You will also need to make sure that the figure's and axes' handles are 'visible' to the Callback functions. This will make them available for plotting commands. Set their HandleVisibility properties to 'callback'. You should also check the GUI's Command-line accessibility option:
% 1. Choose GUI Options (R13) or Application Options (R12.1/R12.0) from the Tools menu of GUIDE.
% 2. In the Command-line accessiblity pulldown menu, select 'On'.
% 3. Click OK.

FrameNo = floor(get(handles.SliderFrames,'Value'));
FrameNoPre = getappdata(handles.SliderFrames,'FrameNoPre');
assignin('base','FrameNo',FrameNo);
assignin('base','FrameNoPre',FrameNoPre);

% If FrameNo is the same as FrameNoPre, don't need to implement code
if FrameNo ~= FrameNoPre
    % Get ROI for red and green rectangles
    ROI=getappdata(handles.SliderFrames,'ROI');
    
    % Get data_Preview
    data_Preview = getappdata(handles.SliderFrames,'data_Preview');
    data_Drift = getappdata(handles.SliderFrames,'data_Drift');
    assignin('base','data_Preview',data_Preview);
    BrightMaxVal = get(handles.SliderBrightMax,'Value');
    BrightMaxVal = floor(BrightMaxVal*100)/100;
    BrightMinVal = get(handles.SliderBrightMin,'Value');
    BrightMinVal = floor(BrightMinVal*100)/100;

    BrightMaxPre = getappdata(handles.SliderFrames,'BrightMaxPre');
    BrightMinPre = getappdata(handles.SliderFrames,'BrightMinPre');
    
    if BrightMaxVal~=BrightMaxPre || BrightMinVal~=BrightMinPre
        if isempty(data_Drift)
            FullPath = getappdata(handles.SliderFrames,'FullPath');
            info = imfinfo(FullPath);
            dim = [info(1).Height info(1).Width numel(info)];
            data_adj = uint8(zeros(dim));
            for i = 1:dim(3)
                data = imread(FullPath, i, 'Info', info);
                data_adj(:,:,i)=round(imadjust(data,[BrightMinVal BrightMaxVal],[])*(255/65535));
                %data_adj(:,:,i)=imadjust(data(:,:,i),[BrightMinVal BrightMaxVal],[]);
                if mod(i,10)==0
                    set(handles.listbox1,'String',{['Adjusting brightness: ' num2str(i) ' / ' num2str(dim(3))]});
                    drawnow();
                end
            end
        else
            dim = size(data_Drift);
            data_adj = uint8(zeros(dim));
            for i = 1:dim(3)
                data_adj(:,:,i)=round(imadjust(data_Drift(:,:,i),[BrightMinVal BrightMaxVal],[]));
            end
        end
        set(handles.listbox1,'String',{['Adjusting brightness: ' num2str(dim(3)) ' / ' num2str(dim(3))]});
        drawnow();
        setappdata(handles.SliderFrames,'data_adj',data_adj);
        setappdata(handles.SliderFrames,'BrightMaxPre',BrightMaxVal);
        setappdata(handles.SliderFrames,'BrightMinPre',BrightMinVal);
    end

    %data_adj=imadjust(data,[BrightMinVal BrightMaxVal],[]);
    axes(handles.Preview); hold off;
    xlimit = xlim; ylimit = ylim;
    data_adj=getappdata(handles.SliderFrames,'data_adj');
    data_adj=data_adj(:,:,FrameNo);
    if ~isempty(ROI); data_adj(ROI)=255; end;
    assignin('base','data_adj',data_adj);
    imshow(data_adj);
    xlim(xlimit); ylim(ylimit);
    setappdata(handles.SliderFrames,'FrameNoPre',FrameNo);
    set(handles.TextFrames,'String',FrameNo);
    set(handles.SliderFrames,'Value',FrameNo);

end

% --- Executes on button press in DeleteFrame.
function DeleteFrame_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes when selected object is changed in Modes.
function Modes_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Modes 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

newButton=get(eventdata.NewValue,'tag');
assignin('base','newButton',newButton);
switch newButton
    case 'Correct'
        set(handles.listbox1,'String',{'Correct mode selected'});
        set(handles.OpenFile,'Visible','on');
        %set(handles.Screen,'Visible','off');
        %set(handles.OpenReference,'Visible','on');
        %set(handles.StartDrawing,'Visible','off');
        setappdata(handles.SliderFrames,'Modes','Correct');
        %set(handles.CreateROI,'Visible','off');
        %set(handles.StartFrame,'Visible','off');
        %set(handles.EndFrame,'Visible','off');
        set(handles.SelectRegion,'Visible','off');
        set(handles.CreateRegion,'Visible','off');
        set(handles.DeletePoint,'Visible','off');
        %set(handles.FindPoints,'Visible','off');
        set(handles.listbox3,'Visible','off');
        set(handles.uipanel6,'Visible','off');
        set(handles.Hist,'Visible','on');
        set(handles.uipanel2,'Visible','on');
        set(handles.Modes,'Visible','on');
        set(handles.StartFIONA,'Visible','off');
        set(handles.BrowseFolder,'String','Browse Folder');
        set(handles.BrowseFile,'String','Browse File');
        %set(handles.PostProcessing,'Visible','off');
        set(handles.uipanel7,'Visible','off');
        set(handles.listbox1,'Visible','on');
        %set(handles.DeleteFrame,'Visible','on');
        set(handles.SearchText,'Visible','off');
        %set(handles.TextName,'Visible','off');
        %set(handles.FolderName,'Visible','off');
        %set(handles.FileNameButton,'Visible','off');
        set(handles.listbox2,'Visible','off');
        %set(handles.Save,'Visible','off');
        %set(handles.PreviewButton,'Visible','off');
        %set(handles.DeleteCroppedFile,'Visible','off');
        set(handles.SliderFrames,'Visible','on');
        set(handles.TextFrames,'Visible','on');
        set(handles.listbox4,'Visible','off');
        %set(handles.CheckMapping,'Visible','off');
        set(handles.OpenCleaned,'Visible','off');
        set(handles.OpenTransformed,'Visible','off');
        set(handles.SingleFile,'Visible','off');
        set(handles.AllFiles,'Visible','off');
        set(handles.FinalFrameNo,'Visible','off');
        set(handles.TextFinalFrame,'Visible','off');
        set(handles.DriftPanel,'Visible','off');
    case 'Split'
        set(handles.listbox1,'String',{'Split mode selected'});
        set(handles.OpenFile,'Visible','on');
        %set(handles.Screen,'Visible','off');
        %set(handles.OpenReference,'Visible','on');
        %set(handles.StartDrawing,'Visible','off');
        setappdata(handles.SliderFrames,'Modes','Split');
        %set(handles.CreateROI,'Visible','off');
        %set(handles.StartFrame,'Visible','off');
        %set(handles.EndFrame,'Visible','off');
        set(handles.SelectRegion,'Visible','on');
        set(handles.CreateRegion,'Visible','on');
        set(handles.DeletePoint,'Visible','on');
        set(handles.FindPoints,'Visible','on');
        set(handles.listbox3,'Visible','off');
        set(handles.uipanel6,'Visible','off');
        set(handles.Hist,'Visible','on');
        set(handles.uipanel2,'Visible','on');
        set(handles.Modes,'Visible','on');
        set(handles.StartFIONA,'Visible','off');
        set(handles.BrowseFolder,'String','Browse Folder');
        set(handles.BrowseFile,'String','Browse File');
        %set(handles.PostProcessing,'Visible','off');
        set(handles.uipanel7,'Visible','off');
        set(handles.listbox1,'Visible','on');
        %set(handles.DeleteFrame,'Visible','on');
        set(handles.SearchText,'Visible','off');
        %set(handles.TextName,'Visible','off');
        %set(handles.FolderName,'Visible','off');
        %set(handles.FileNameButton,'Visible','off');
        set(handles.listbox2,'Visible','off');
        %set(handles.Save,'Visible','off');
        %set(handles.PreviewButton,'Visible','off');
        %set(handles.DeleteCroppedFile,'Visible','off');
        set(handles.SliderFrames,'Visible','on');
        set(handles.TextFrames,'Visible','on');
        set(handles.listbox4,'Visible','off');
        %set(handles.CheckMapping,'Visible','off');
        set(handles.OpenCleaned,'Visible','off');
        set(handles.OpenTransformed,'Visible','off');
        set(handles.SingleFile,'Visible','off');
        set(handles.AllFiles,'Visible','off');
        set(handles.FinalFrameNo,'Visible','off');
        set(handles.TextFinalFrame,'Visible','off');
        set(handles.DriftPanel,'Visible','off');
    case 'DriftCorrect'
        set(handles.listbox2,'String',{'Drift Correct mode selected'});
        set(handles.OpenFile,'Visible','on');
        %set(handles.Screen,'Visible','off');
        %set(handles.OpenReference,'Visible','on');
        %set(handles.StartDrawing,'Visible','off');
        setappdata(handles.SliderFrames,'Modes','DriftCorrect');
        %set(handles.CreateROI,'Visible','off');
        %set(handles.StartFrame,'Visible','off');
        %set(handles.EndFrame,'Visible','off');
        set(handles.SelectRegion,'Visible','off');
        set(handles.CreateRegion,'Visible','off');
        set(handles.DeletePoint,'Visible','off');
        set(handles.FindPoints,'Visible','on');
        set(handles.listbox3,'Visible','off');
        set(handles.uipanel6,'Visible','off');
        set(handles.Hist,'Visible','on');
        set(handles.uipanel2,'Visible','on');
        set(handles.Modes,'Visible','on');
        set(handles.StartFIONA,'Visible','off');
        set(handles.BrowseFolder,'String','Code Path');
        set(handles.BrowseFile,'String','Data Folder');
        %set(handles.PostProcessing,'Visible','off');
        set(handles.uipanel7,'Visible','off');
        set(handles.listbox1,'Visible','off');
        %set(handles.DeleteFrame,'Visible','on');
        set(handles.SearchText,'Visible','off');
        %set(handles.TextName,'Visible','off');
        %set(handles.FolderName,'Visible','off');
        %set(handles.FileNameButton,'Visible','off');
        set(handles.listbox2,'Visible','on');
        %set(handles.Save,'Visible','off');
        %set(handles.PreviewButton,'Visible','off');
        %set(handles.DeleteCroppedFile,'Visible','off');
        set(handles.SliderFrames,'Visible','on');
        set(handles.TextFrames,'Visible','on');
        set(handles.listbox4,'Visible','off');
        %set(handles.CheckMapping,'Visible','off');
        set(handles.OpenCleaned,'Visible','off');
        set(handles.OpenTransformed,'Visible','off');
        set(handles.SingleFile,'Visible','on');
        set(handles.AllFiles,'Visible','on');
        set(handles.FinalFrameNo,'Visible','on');
        set(handles.TextFinalFrame,'Visible','on');
        set(handles.DriftPanel,'Visible','on');
    case 'PostProcess'
        %PostProcessVal = get(handles.PostProcessing,'Value');
        set(handles.SelectRegion,'Visible','off');
        set(handles.CreateRegion,'Visible','off');
        set(handles.DeletePoint,'Visible','off');
        set(handles.listbox1,'Visible','off');
        set(handles.DeleteFrame,'Visible','off');
        %set(handles.TextName,'Visible','off');
        %set(handles.FolderName,'Visible','off');
        %set(handles.FileNameButton,'Visible','off');
        set(handles.listbox2,'Visible','off');
        %set(handles.Save,'Visible','off');
        %set(handles.PreviewButton,'Visible','off');
        %set(handles.DeleteCroppedFile,'Visible','off');
        set(handles.SliderFrames,'Visible','off');
        set(handles.TextFrames,'Visible','off');
        %set(handles.CheckMapping,'Visible','off');
        
        %set(handles.listbox1,'String',{'Drawing mode selected'});
        set(handles.OpenFile,'Visible','off');
        %set(handles.Screen,'Visible','off');
        %set(handles.OpenReference,'Visible','off');
        %set(handles.StartDrawing,'Visible','off');
        setappdata(handles.SliderFrames,'Modes','PostProcess');
        %set(handles.CreateROI,'Visible','off');
        %set(handles.StartFrame,'Visible','off');
        %set(handles.EndFrame,'Visible','off');
        set(handles.FindPoints,'Visible','off');
        set(handles.listbox3,'Visible','off');
        set(handles.uipanel6,'Visible','off');
        set(handles.Hist,'Visible','off');
        set(handles.uipanel2,'Visible','off');
        %set(handles.Modes,'Visible','off');
        set(handles.StartFIONA,'Visible','off');
        set(handles.BrowseFolder,'String','Code Path');
        set(handles.BrowseFile,'String','Data Folder');
        %set(handles.PostProcessing,'Visible','on');
        set(handles.uipanel7,'Visible','on');
        set(handles.SearchText,'Visible','on');
        set(handles.listbox4,'Visible','on');
        set(handles.OpenCleaned,'Visible','on');
        set(handles.OpenTransformed,'Visible','on');
        set(handles.SingleFile,'Visible','off');
        set(handles.AllFiles,'Visible','off');
        set(handles.FinalFrameNo,'Visible','off');
        set(handles.TextFinalFrame,'Visible','off');
        set(handles.DriftPanel,'Visible','off');
    case 'FIONA'
        %PostProcessVal = get(handles.PostProcessing,'Value');
        set(handles.SelectRegion,'Visible','off');
        set(handles.CreateRegion,'Visible','off');
        set(handles.DeletePoint,'Visible','off');
        set(handles.listbox1,'Visible','off');
        set(handles.DeleteFrame,'Visible','off');
        %set(handles.TextName,'Visible','off');
        %set(handles.FolderName,'Visible','off');
        %set(handles.FileNameButton,'Visible','off');
        set(handles.listbox2,'Visible','off');
        %set(handles.Save,'Visible','off');
        %set(handles.PreviewButton,'Visible','off');
        %set(handles.DeleteCroppedFile,'Visible','off');
        set(handles.SliderFrames,'Visible','off');
        set(handles.TextFrames,'Visible','off');
        %set(handles.CheckMapping,'Visible','off');

        set(handles.listbox1,'String',{'Drawing mode selected'});
        set(handles.OpenFile,'Visible','off');
        %set(handles.Screen,'Visible','off');
        %set(handles.OpenReference,'Visible','off');
        %set(handles.StartDrawing,'Visible','off');
        setappdata(handles.SliderFrames,'Modes','FIONA');
        %set(handles.CreateROI,'Visible','off');
        %set(handles.StartFrame,'Visible','off');
        %set(handles.EndFrame,'Visible','off');
        set(handles.FindPoints,'Visible','off');
        set(handles.listbox3,'Visible','on');
        set(handles.uipanel6,'Visible','on');
        set(handles.Hist,'Visible','off');
        set(handles.uipanel2,'Visible','off');
        set(handles.Modes,'Visible','on');
        set(handles.StartFIONA,'Visible','on');
        set(handles.BrowseFolder,'String','Code Path');
        set(handles.BrowseFile,'String','Data Folder');
        %set(handles.PostProcessing,'Visible','on');
        set(handles.uipanel7,'Visible','off');
        set(handles.SearchText,'Visible','off');
        set(handles.listbox4,'Visible','off');
        set(handles.OpenCleaned,'Visible','off');
        set(handles.OpenTransformed,'Visible','off');
        set(handles.SingleFile,'Visible','off');
        set(handles.AllFiles,'Visible','off');
        set(handles.FinalFrameNo,'Visible','off');
        set(handles.TextFinalFrame,'Visible','off');
        set(handles.DriftPanel,'Visible','off');
        
        FIONASetting = getappdata(handles.SliderFrames,'FIONASetting');
        if isempty(FIONASetting)
            if exist('FIONASetting.mat','file')
                load('FIONASetting.mat','FIONASetting');
            else
                Opt1 = 'ActPix'; Val1 = '16000';            % Actual pixel size in nanometer
                Opt2 = 'ObjMag'; Val2 = '100';              % Objective magnification
                Opt3 = 'AddMag'; Val3 = '1.5';              % Additional magnification
                Opt4 = 'Camera'; Val4 = '  Scorpion';
                Opt5 = 'Readout'; Val5 = '  10 MHz 14 bit';
                Opt6 = 'Preamp'; Val6 = '  5.2x';
                Opt7 = 'EMgain'; Val7 = '200';              % Electron multiplying (EM) gain setting of camera during acquisition
                Opt8 = 'Parallel'; Val8 = '  No';
                Opt9 = 'CCDSens'; Val9 = '12.13';
                Opt10 = 'CameraVal'; Val10 = 3;
                Opt11 = 'ReadoutVal'; Val11 = 4;
                Opt12 = 'PreampVal'; Val12 = 3;
                Opt13 = 'ParallelVal'; Val13 = 2;
                FIONASetting = struct(Opt1,Val1,Opt2,Val2,Opt3,Val3,Opt4,Val4,Opt5,Val5,Opt6,Val6,Opt7,Val7,Opt8,Val8,Opt9,Val9,Opt10,Val10,Opt11,Val11,Opt12,Val12,Opt13,Val13);
                save('FIONASetting.mat','FIONASetting');
            end
            setappdata(handles.SliderFrames,'FIONASetting',FIONASetting);
        end
        
        String1 = ['Actual Pixel :                        ' FIONASetting.ActPix ' nm'];
        String2 = ['Objective Magnification :      ' FIONASetting.ObjMag];
        String3 = ['Additional Magnification :     ' FIONASetting.AddMag];
        String4 = ['Camera :                            ' FIONASetting.Camera];
        String5 = ['Readout Rate :                   ' FIONASetting.Readout];
        String6 = ['Preamp Setting :                 ' FIONASetting.Preamp];
        String7 = ['EM gain :                               ' FIONASetting.EMgain];
        String8 = ['Parallel Processing:            ' FIONASetting.Parallel];
        String9 = ['CCD Sensitivity :                   ' FIONASetting.CCDSens];

        set(handles.listbox3,'String',{String1;String2;String3;String4;String5;String6;String7;String8;String9});

        % Add FolderName to the string of FileName
        FolderName = getappdata(handles.SliderFrames,'FolderName');
        % If FolderName exist, place FolderName on the string, otherwise,
        % leave it as it is
        if ~isempty(FolderName)
            set(handles.FileName,'string',FolderName);
        end
        
        %assignin('base','FIONASetting',FIONASetting);
end

% --- Executes on button press in SearchText.
function SearchText_Callback(hObject, eventdata, handles)
% hObject    handle to SearchText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get CodePath
CodePath = get(handles.FileLocation,'string');

% Get DataPath
Folder = get(handles.FileName,'string');
% Check if Folder is empty. if it not, check if the directory exists
if isempty(Folder)
    % Ask for data path if Folder is empty
    DataPath = uigetdir;
else
    DataPath = [CodePath '\' Folder];
    % If directory does not exist, get directory
    if ~exist(DataPath,'dir'); DataPath = uigetdir; end
end
name=strrep(DataPath,[CodePath '\'],'');
set(handles.FileName,'String',name);
cd(DataPath);

% Create Cleaned, Transformed and T-test folder if they have not existed
if exist('Cleaned','dir')~=7; mkdir('Cleaned'); end
if exist('Transformed','dir')~=7; mkdir('Transformed'); end
if exist('T-test','dir')~=7; mkdir('T-test'); end

% Get OpenCleanedFile and OpenTransformedFile and change DataPath if any is switched on
OpenCleanedFile = get(handles.OpenCleaned,'Value');
OpenTransformedFile = get(handles.OpenTransformed,'Value');
if OpenCleanedFile; DataPath = [DataPath '\Cleaned']; end
if OpenTransformedFile; DataPath = [DataPath '\Transformed']; end
cd(DataPath);

% Find files
FileType = '.txt';
FileIn = dir(['*' FileType]);
FileInput = cell(length(FileIn),1);
FileInputName = cell(length(FileIn),1);
for ind=1:length(FileIn)
    FileInput{ind}=FileIn(ind).name;
    FileInputName{ind}=strrep(FileIn(ind).name,FileType,'');
end
FileInputString = cell(length(FileInput),1);
%assignin('base','FileInput',FileInput);

% Display the number of files and file name on listbox
if OpenCleanedFile == 1 || OpenTransformedFile == 1
    num_images = length(FileInput);
    String1 = ['No of files : ' num2str(num_images)];
    set(handles.listbox4,'String',[String1; ' '; FileInput]);
else
    num_images = length(FileInput);
    for i = 1:num_images
        if exist(['Cleaned\' FileInputName{i} '-Cleaned.txt'],'file') == 2
            FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"     (saved)'];
        else
            FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"'];
        end        
        if exist(['Transformed\' FileInputName{i} '-Transformed.txt'],'file') == 2
            FileInputString{i} = [FileInputString{i} '     (transformed)'];
        end        
    end
    String1 = ['No of files : ' num2str(num_images)];
    set(handles.listbox4,'String',[String1; ' '; FileInputString]);
    assignin('base','FileInputName',FileInputName);
    assignin('base','FileInputString',FileInputString);
end

% Open and plot the first file
fid=fopen(FileInput{1});
Input = textscan(fid,'%*f%*f%*f%f%f%*f%*f','CommentStyle','##');
Input = [Input{1} Input{2} (1:length(Input{1}))'];
axes(handles.Preview); hold off;
plot(Input(:,1),Input(:,2),'g'); hold on;
plot(Input(1,1),Input(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(Input(end,1),Input(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
hold off;
datacursormode on;

% Store the FileInput and FileInputName
setappdata(handles.SliderFrames,'FileInput',FileInput);
setappdata(handles.SliderFrames,'FileInputName',FileInputName);

% Set the FileNo on SliderFileNo and TextFileNo to be 1
set(handles.SliderFileNo,'Value',1);
set(handles.TextFileNo,'String',1);

% Set properties of SliderFileNo
set(handles.SliderFileNo,'Max',num_images);
set(handles.SliderFileNo,'SliderStep',[1/(num_images-1) 1/(num_images-1)]);

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',Input);

% Store CodePath and DataPath
setappdata(handles.SliderFrames,'CodePath',CodePath);
setappdata(handles.SliderFrames,'DataPath',DataPath);

% Switch to CodePath
cd(CodePath);

% --- Executes on button press in OpenCleaned.
function OpenCleaned_Callback(hObject, eventdata, handles)
% hObject    handle to OpenCleaned (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OpenCleaned

% Switch OpenTransformed value off
set(handles.OpenTransformed,'Value',0);

% Get CodePath
CodePath = get(handles.FileLocation,'string');

% Get DataPath
Folder = get(handles.FileName,'string');
% Check if Folder is empty. if it not, check if the directory exists
if isempty(Folder)
    DataPath = uigetdir;
else
    DataPath = [CodePath '\' Folder];
    % If directory does not exist, get directory
    if ~exist(DataPath,'dir'); DataPath = uigetdir; end
end
name=strrep(DataPath,[CodePath '\'],'');
set(handles.FileName,'String',name);
cd(DataPath);

% Get the FileNo from TextFileNo
FileIn = dir('*.txt');
if length(FileIn) > 1
    FileNo = floor(str2double(get(handles.TextFileNo,'String')));
else
    FileNo = 1;
end

% Create Cleaned and Transformed folder if they have not existed
if exist('Cleaned','dir')~=7; mkdir('Cleaned'); end
if exist('Transformed','dir')~=7; mkdir('Transformed'); end
if exist('T-test','dir')~=7; mkdir('T-test'); end

% Get OpenCleanedFile and change DataPath if it is switched on
OpenCleanedFile = get(handles.OpenCleaned,'Value');
if OpenCleanedFile; DataPath = [DataPath '\Cleaned']; end
cd(DataPath);

% Find files
FileType = '.txt';
FileIn = dir(['*' FileType]);
FileInput = cell(length(FileIn),1);
FileInputName = cell(length(FileIn),1);
for ind=1:length(FileIn)
    FileInput{ind}=FileIn(ind).name;
    FileInputName{ind}=strrep(FileIn(ind).name,FileType,'');
end
FileInputString = cell(length(FileInput),1);
%assignin('base','FileInput',FileInput);

% Display the number of files and file name on listbox
if OpenCleanedFile
    num_images = length(FileInput);
    String1 = ['No of files : ' num2str(num_images)];
    set(handles.listbox4,'String',[String1; ' '; FileInput]);
else
    num_images = length(FileInput);
    for i = 1:num_images
        if exist(['Cleaned\' FileInputName{i} '-Cleaned.txt'],'file') == 2
            FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"     (saved)'];
        else
            FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"'];
        end
        if exist(['Transformed\' FileInputName{i} '-Transformed.txt'],'file') == 2
            FileInputString{i} = [FileInputString{i} '     (transformed)'];
        end     
    end
    String1 = ['No of files : ' num2str(num_images)];
    set(handles.listbox4,'String',[String1; ' '; FileInputString]);
end

% Open and plot the first file
if FileNo > num_images; FileNo = num_images; end;
fid=fopen(FileInput{FileNo});
if OpenCleanedFile; 
    Input = textscan(fid,'%f%f%f','CommentStyle','##');
else
    Input = textscan(fid,'%*f%*f%*f%f%f%*f%*f','CommentStyle','##');
end
Input = [Input{1} Input{2} (1:length(Input{1}))'];
axes(handles.Preview); hold off;
plot(Input(:,1),Input(:,2),'g'); hold on;
plot(Input(1,1),Input(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(Input(end,1),Input(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
hold off;
datacursormode on;

% Store the FileInput and FileInputName
setappdata(handles.SliderFrames,'FileInput',FileInput);
setappdata(handles.SliderFrames,'FileInputName',FileInputName);

% Set the FileNo on SliderFileNo and TextFileNo to be 1
set(handles.SliderFileNo,'Value',FileNo);
set(handles.TextFileNo,'String',FileNo);

% Set properties of SliderFileNo
if num_images > 1
    set(handles.SliderFileNo,'Max',num_images);
    set(handles.SliderFileNo,'SliderStep',[1/(num_images-1) 1/(num_images-1)]);
end

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',Input);

% Store CodePath and DataPath
setappdata(handles.SliderFrames,'CodePath',CodePath);
setappdata(handles.SliderFrames,'DataPath',DataPath);

% Switch to CodePath
cd(CodePath);

% --- Executes on button press in OpenTransformed.
function OpenTransformed_Callback(hObject, eventdata, handles)
% hObject    handle to OpenTransformed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OpenTransformed

% Switch OpenCleaned value off
set(handles.OpenCleaned,'Value',0);

% Get CodePath
CodePath = get(handles.FileLocation,'string');

% Get DataPath
Folder = get(handles.FileName,'string');
% Check if Folder is empty. if it not, check if the directory exists
if isempty(Folder)
    DataPath = uigetdir;
else
    DataPath = [CodePath '\' Folder];
    % If directory does not exist, get directory
    if ~exist(DataPath,'dir'); DataPath = uigetdir; end
end
name=strrep(DataPath,[CodePath '\'],'');
set(handles.FileName,'String',name);
cd(DataPath);

% Get the FileNo from TextFileNo
FileIn = dir('*.txt');
if length(FileIn) > 1
    FileNo = floor(str2double(get(handles.TextFileNo,'String')));
else
    FileNo = 1;
end

% Create Cleaned and Transformed folder if they have not existed
if exist('Cleaned','dir')~=7; mkdir('Cleaned'); end
if exist('Transformed','dir')~=7; mkdir('Transformed'); end
if exist('T-test','dir')~=7; mkdir('T-test'); end

% Get OpenCleanedFile and change DataPath if it is switched on
OpenTransformedFile = get(handles.OpenTransformed,'Value');
if OpenTransformedFile; DataPath = [DataPath '\Transformed']; end
cd(DataPath);

% Find the files without 'Green' or 'Red' on the file name
FileType = '.txt';
FileIn = dir(['*' FileType]);
FileInput = cell(length(FileIn),1);
FileInputName = cell(length(FileIn),1);
for ind=1:length(FileIn)
    FileInput{ind}=FileIn(ind).name;
    FileInputName{ind}=strrep(FileIn(ind).name,FileType,'');
end
FileInputString = cell(length(FileInput),1);
%assignin('base','FileInput',FileInput);

% Display the number of files and file name on listbox
if OpenTransformedFile
    num_images = length(FileInput);
    String1 = ['No of files : ' num2str(num_images)];
    set(handles.listbox4,'String',[String1; ' '; FileInput]);
else
    num_images = length(FileInput);
    for i = 1:num_images
        if exist(['Cleaned\' FileInputName{i} '-Cleaned.txt'],'file') == 2
            FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"     (saved)'];
        else
            FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"'];
        end
        if exist(['Transformed\' FileInputName{i} '-Transformed.txt'],'file') == 2
            FileInputString{i} = [FileInputString{i} '     (transformed)'];
        end     
    end
    String1 = ['No of files : ' num2str(num_images)];
    set(handles.listbox4,'String',[String1; ' '; FileInputString]);
end

% Open and plot the first file
if FileNo > num_images; FileNo = num_images; end;
fid=fopen(FileInput{FileNo});
if OpenTransformedFile; 
    Input = textscan(fid,'%f%f%f','CommentStyle','##');
else
    Input = textscan(fid,'%*f%*f%*f%f%f%*f%*f','CommentStyle','##');
end
Input = [Input{1} Input{2} (1:length(Input{1}))'];
axes(handles.Preview); hold off;
plot(Input(:,1),Input(:,2),'g'); hold on;
plot(Input(1,1),Input(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(Input(end,1),Input(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
hold off;
datacursormode on;

% Store the FileInput and FileInputName
setappdata(handles.SliderFrames,'FileInput',FileInput);
setappdata(handles.SliderFrames,'FileInputName',FileInputName);

% Set the FileNo on SliderFileNo and TextFileNo to be FileNo
set(handles.SliderFileNo,'Value',FileNo);
set(handles.TextFileNo,'String',FileNo);

% Set properties of SliderFileNo
if num_images > 1
    set(handles.SliderFileNo,'Max',num_images);
    set(handles.SliderFileNo,'SliderStep',[1/(num_images-1) 1/(num_images-1)]);
end

% Store dataInput
setappdata(handles.SliderFrames,'dataInput',Input);

% Store CodePath and DataPath
setappdata(handles.SliderFrames,'CodePath',CodePath);
setappdata(handles.SliderFrames,'DataPath',DataPath);

% Switch to CodePath
cd(CodePath);

% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4

% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in PolyTransform.
function PolyTransform_Callback(hObject, eventdata, handles)
% hObject    handle to PolyTransform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get CodePath, DataPath and FileInput
CodePath = getappdata(handles.SliderFrames,'CodePath');
DataPath = getappdata(handles.SliderFrames,'DataPath');
cd(DataPath);

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');

% Get OpenCleanedFile and OpenTransformedFile
OpenCleanedFile = get(handles.OpenCleaned,'Value');
OpenTransformedFile = get(handles.OpenTransformed,'Value');

% Get FileInput and FileInputName
FileInputName=getappdata(handles.SliderFrames,'FileInputName');

% Get the FileNo from SliderFileNo
FileNo = floor(get(handles.SliderFileNo,'Value'));

% Open cleaned data and assign data to xFin and yFin
xFin=dataInput(:,1); yFin=dataInput(:,2); tInput=dataInput(:,3);

% Make sure that the input has a slope below 1. The more vertical it
% is, the less accurate is the transformation
[xmin,indmin]=min(xFin);
[xmax,indmax]=max(xFin);
ymin=yFin(indmin);
ymax=yFin(indmax);
Slope=abs((ymax-ymin)/(xmax-xmin));
% Force rotate
if 0
    Slope = 2;
end
if Slope>1
    xFin=dataInput(:,2);
    yFin=dataInput(:,1);
end

% Determine whether to fit data to fifth or sixth degree polynomial
xFin=xFin-xFin(1);
yFin=yFin-yFin(1);
%xFin=xFin-xFin(floor(length(xFin)/2));
%yFin=yFin-yFin(floor(length(yFin)/2));
FinePrecision=0.01;             % Use numbers in power of 10
MidPrecision=0.1;               % Use numbers in power of 10
RoughPrecision=1;               % Use numbers in power of 10
xlineFine=(floor(min(xFin)-((max(xFin)-min(xFin))/10)):FinePrecision:floor(max(xFin)+((max(xFin)-min(xFin))/10)))';
xlineMid=(floor(min(xFin)-((max(xFin)-min(xFin))/10)):MidPrecision:floor(max(xFin)+((max(xFin)-min(xFin))/10)))';
xlineRough=(floor(min(xFin)-((max(xFin)-min(xFin))/10)):RoughPrecision:floor(max(xFin)+((max(xFin)-min(xFin))/10)))';
[fit5,S5]=polyfit(xFin,yFin,5);
[fit6,S6]=polyfit(xFin,yFin,6);
if S5.normr<S6.normr
    ylineFine=fit5(1)*xlineFine.^5+fit5(2)*xlineFine.^4+fit5(3)*xlineFine.^3+fit5(4)*xlineFine.^2+fit5(5)*xlineFine+fit5(6);
    ylineMid=fit5(1)*xlineMid.^5+fit5(2)*xlineMid.^4+fit5(3)*xlineMid.^3+fit5(4)*xlineMid.^2+fit5(5)*xlineMid+fit5(6);
    ylineRough=fit5(1)*xlineRough.^5+fit5(2)*xlineRough.^4+fit5(3)*xlineRough.^3+fit5(4)*xlineRough.^2+fit5(5)*xlineRough+fit5(6);
else
    ylineFine=fit6(1)*xlineFine.^6+fit6(2)*xlineFine.^5+fit6(3)*xlineFine.^4+fit6(4)*xlineFine.^3+fit6(5)*xlineFine.^2+fit6(6)*xlineFine+fit6(7);
    ylineMid=fit6(1)*xlineMid.^6+fit6(2)*xlineMid.^5+fit6(3)*xlineMid.^4+fit6(4)*xlineMid.^3+fit6(5)*xlineMid.^2+fit6(6)*xlineMid+fit6(7);
    ylineRough=fit6(1)*xlineRough.^6+fit6(2)*xlineRough.^5+fit6(3)*xlineRough.^4+fit6(4)*xlineRough.^3+fit6(5)*xlineRough.^2+fit6(6)*xlineRough+fit6(7);
end

% Plot raw and fitted trace
plot(xFin,yFin,'g');
hold on;
plot(xlineRough,ylineRough,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',1);
title(['Raw Trace '  FileInputName{FileNo}],'FontName','Palatino Linotype');
xlabel('x (nm)','FontName','Palatino Linotype');
ylabel('y (nm)','FontName','Palatino Linotype');
hold off

% Find distance along axis for fine curve
OnAxisDistance=sqrt((xlineFine(2:end)-xlineFine(1:end-1)).^2+(ylineFine(2:end)-ylineFine(1:end-1)).^2);
OnAxisDistance=[0;OnAxisDistance];
OnAxisDistance=cumsum(OnAxisDistance);
% Need to check that the length for OnAxisDistance before and after
% cumsum are the same
%figure;
%plot(OnAxisDistance);

% Initiate vectors
OnAxis=zeros(length(xFin),1);
OffAxis=zeros(length(xFin),1);

% Find closest point
for k=1:length(xFin)
    % Take the point on the same x-coordinate and determine direction
    indTest=find(xlineRough==floor(xFin(k)));
    DistanceOnSite=(xlineRough(indTest)-xFin(k))*(xlineRough(indTest)-xFin(k))+(ylineRough(indTest)-yFin(k))*(ylineRough(indTest)-yFin(k));
    DistanceForward=(xlineRough(indTest+1)-xFin(k))*(xlineRough(indTest+1)-xFin(k))+(ylineRough(indTest+1)-yFin(k))*(ylineRough(indTest+1)-yFin(k));
    DistanceBackward=(xlineRough(indTest-1)-xFin(k))*(xlineRough(indTest-1)-xFin(k))+(ylineRough(indTest-1)-yFin(k))*(ylineRough(indTest-1)-yFin(k));

    % Finding out the direction to search, left (-1) or right (+1) of the indTest
    [~,Direction] = min([DistanceBackward DistanceOnSite DistanceForward]);
    Direction=Direction-2;

    % Loop and search for the minimum distance
    loop=1;
    if Direction == 0
        loop=0;
        indPre=indTest-1;   % Indeces used to refine search for minimum in between indPre and indPost
        indPost=indTest+1;        
    else
        while loop==1
            indTest=indTest+Direction;
            % Compute distance difference and stop looping when it is positive
            xVec=xlineRough(indTest)-xFin(k);
            yVec=ylineRough(indTest)-yFin(k);
            DistanceNew=xVec*xVec+yVec*yVec;
            DistanceDiff=DistanceNew-DistanceOnSite;
            DistanceOnSite=DistanceNew;
            if DistanceDiff>0
                if Direction==1
                    loop=0;
                    indPre=indTest-2;
                    indPost=indTest;
                else
                    loop=0;
                    indPre=indTest;
                    indPost=indTest+2;
                end
            end
        end
    end

    % Refine search for minimum in between indPre and indPost
    indPre=((indPre-1)*(RoughPrecision/MidPrecision))+1;
    indPost=((indPost-1)*(RoughPrecision/MidPrecision))+1;
    xRefine=xlineMid(indPre:indPost);
    yRefine=ylineMid(indPre:indPost);
    RefineDistance=(xRefine-xFin(k)).^2+(yRefine-yFin(k)).^2;
    [~,indTest] = min(RefineDistance);
    indTest=indPre+indTest-1;
    indPre=((indPre-1)*(MidPrecision/FinePrecision))+1;
    indPost=((indPost-1)*(MidPrecision/FinePrecision))+1;
    xRefine=xlineFine(indPre:indPost);
    yRefine=ylineFine(indPre:indPost);
    RefineDistance=(xRefine-xFin(k)).^2+(yRefine-yFin(k)).^2;
    [~,indTest] = min(RefineDistance);
    indTest=indPre+indTest-1;
    OnAxis(k)=OnAxisDistance(indTest);
    OffAxisDirection=sign(ylineFine(indTest)-yFin(k));
    OffAxis(k)=OffAxisDirection*sqrt((xlineFine(indTest)-xFin(k))^2+(ylineFine(indTest)-yFin(k))^2);
end

% Invert negative sloped traces
if OnAxis(end)-OnAxis(1)<0
    OnAxisG = flipud(OnAxis);
    OffAxisG = flipud(OffAxis);
else
    OnAxisG = OnAxis;
    OffAxisG = OffAxis;
end

% Plot rotated trace
axes(handles.Preview); hold off;
plot(tInput, OnAxisG)
title(['Transformed Trace ' FileInputName{FileNo}],'FontName','Palatino Linotype');
xlabel('Frame','FontName','Palatino Linotype');
ylabel('On-Axis Distance (nm)','FontName','Palatino Linotype');
hold off;

% Save rotated trace variable
setappdata(handles.SliderFrames,'OnAxisTrace',[OnAxisG tInput]);

% Find files
FileType = '.txt';
FileIn = dir(['*' FileType]);
FileInput = cell(length(FileIn),1);
FileInputName = cell(length(FileIn),1);
for ind=1:length(FileIn)
    FileInput{ind}=FileIn(ind).name;
    FileInputName{ind}=strrep(FileIn(ind).name,FileType,'');
end
FileInputString = cell(length(FileInput),1);

if OpenCleanedFile == 1
    DataPath = strrep(DataPath,'\Cleaned','');
    FileName = strrep(FileInputName{FileNo},'-Cleaned','');
    cd(DataPath);
    %assignin('base','DataPath',DataPath);
    
    % Save data
    dlmwrite(['Transformed\' FileName '-Transformed.txt'],[OnAxisG OffAxisG tInput],'\t');

    % Save data for t-test input
    dlmwrite(['T-test\' FileName 'Ttest.txt'],[OnAxisG tInput],'\t');
elseif OpenTransformedFile == 1
    DataPath = strrep(DataPath,'\Transformed','');
    FileName = strrep(FileInputName{FileNo},'-Transformed','');
    cd(DataPath);
    %assignin('base','DataPath',DataPath);

    % Save data
    dlmwrite(['Transformed\' FileName '-Transformed.txt'],[OnAxisG OffAxisG tInput],'\t');

    % Save data for t-test input
    dlmwrite(['T-test\' FileName 'Ttest.txt'],[OnAxisG tInput],'\t');
else
    % Save data
    dlmwrite(['Transformed\' FileInputName{FileNo} '-Transformed.txt'],[OnAxisG OffAxisG tInput],'\t');

    % Save data for t-test input
    dlmwrite(['T-test\' FileInputName{FileNo} 'Ttest.txt'],[OnAxisG tInput],'\t');
    
    % Update listbox
    num_images = length(FileInput);
    for i = 1:num_images
        if exist(['Cleaned\' FileInputName{i} '-Cleaned.txt'],'file') == 2
            FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"     (saved)'];
        else
            FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"'];
        end
        if exist(['Transformed\' FileInputName{i} '-Transformed.txt'],'file') == 2
            FileInputString{i} = [FileInputString{i} '     (transformed)'];
        end     
    end
    String1 = ['No of files : ' num2str(num_images)];
    set(handles.listbox4,'String',[String1; ' '; FileInputString]);
end

% Switch to CodePath
cd(CodePath);

% Plot transformed trace and the fitted curve
h = figure;
plot(OnAxis,OffAxis,'g'); hold on;
plot([min(OnAxis) max(OnAxis)],[0 0],'red','LineWidth',2);
title(['Transformed Trace ' FileInputName{FileNo}],'FontName','Palatino Linotype');
assignin('base','File',FileInputName{FileNo});
xlabel('x (nm)','FontName','Palatino Linotype');
ylabel('y (nm)','FontName','Palatino Linotype');
hold off
setappdata(handles.SliderFrames,'h',h);

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in SelectRegion.
function SelectRegion_Callback(hObject, eventdata, handles)
% hObject    handle to SelectRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Adjust intensity when brightness values are changed
BrightMaxVal = get(handles.SliderBrightMax,'Value');
BrightMaxVal = floor(BrightMaxVal*100)/100;
BrightMinVal = get(handles.SliderBrightMin,'Value');
BrightMinVal = floor(BrightMinVal*100)/100;
BrightMaxPre = getappdata(handles.SliderFrames,'BrightMaxPre');
BrightMinPre = getappdata(handles.SliderFrames,'BrightMinPre');
assignin('base','BrightMaxVal',BrightMaxVal);
assignin('base','BrightMinVal',BrightMinVal);
assignin('base','BrightMaxPre',BrightMaxPre);
assignin('base','BrightMinPre',BrightMinPre);

if BrightMaxVal~=BrightMaxPre || BrightMinVal~=BrightMinPre
    FullPath = getappdata(handles.SliderFrames,'FullPath');
    info = imfinfo(FullPath);
    dim = [info(1).Height info(1).Width numel(info)];
    data_adj = uint8(zeros(dim));
    for i = 1:dim(3)
        data = imread(FullPath, i, 'Info', info);
        data_adj(:,:,i)=round(imadjust(data,[BrightMinVal BrightMaxVal],[])*(255/65535));
        %data_adj(:,:,i)=imadjust(data(:,:,i),[BrightMinVal BrightMaxVal],[]);
        if mod(i,10)==0
            set(handles.listbox1,'String',{['Adjusting brightness: ' num2str(i) ' / ' num2str(dim(3))]});
            drawnow();
        end
    end
    set(handles.listbox1,'String',{['Adjusting brightness: ' num2str(dim(3)) ' / ' num2str(dim(3))]});
    drawnow();
    setappdata(handles.SliderFrames,'data_adj',data_adj);
    setappdata(handles.SliderFrames,'BrightMaxPre',BrightMaxVal);
    setappdata(handles.SliderFrames,'BrightMinPre',BrightMinVal);
else
    data_adj=getappdata(handles.SliderFrames,'data_adj');
    dim = size(data_adj);
end

FrameNo = floor(get(handles.SliderFrames,'Value'));
Pos1 = [0 0 dim(2) (dim(1)/2-floor(dim(1)/32))];
Pos2 = [0 (dim(1)/2+floor(dim(1)/32)) dim(2) (dim(1)/2-floor(dim(1)/32))];
axes(handles.Preview);
hold off;
imshow(data_adj(:,:,FrameNo));
H1 = imrect(handles.Preview,Pos1);
H2 = imrect(handles.Preview,Pos2);
%assignin('base','H1',H1);
setappdata(handles.SliderFrames,'H1',H1);
setappdata(handles.SliderFrames,'H2',H2);
set(handles.listbox1,'String',{'Select two rectangular regions'});
%set(handles.listbox1,'Value',[]);

% --- Executes on button press in CreateRegion.
function CreateRegion_Callback(hObject, eventdata, handles)
% hObject    handle to CreateRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get info on tiff file
FullPath = getappdata(handles.SliderFrames,'FullPath');
info = imfinfo(FullPath);
FileName = get(handles.FileName,'string');
FileType = '.tif';
FileInput = regexprep(FileName, FileType, '');

% Adjust intensity when brightness values are changed
BrightMaxVal = get(handles.SliderBrightMax,'Value');
BrightMaxVal = floor(BrightMaxVal*100)/100;
BrightMinVal = get(handles.SliderBrightMin,'Value');
BrightMinVal = floor(BrightMinVal*100)/100;
BrightMaxPre = getappdata(handles.SliderFrames,'BrightMaxPre');
BrightMinPre = getappdata(handles.SliderFrames,'BrightMinPre');

if BrightMaxVal~=BrightMaxPre || BrightMinVal~=BrightMinPre
    FullPath = getappdata(handles.SliderFrames,'FullPath');
    info = imfinfo(FullPath);
    dim = [info(1).Height info(1).Width numel(info)];
    data_adj = uint8(zeros(dim));
    for i = 1:dim(3)
        data = imread(FullPath, i, 'Info', info);
        data_adj(:,:,i)=round(imadjust(data,[BrightMinVal BrightMaxVal],[])*(255/65535));
        %data_adj(:,:,i)=imadjust(data(:,:,i),[BrightMinVal BrightMaxVal],[]);
        if mod(i,10)==0
            set(handles.listbox1,'String',{['Adjusting brightness: ' num2str(i) ' / ' num2str(dim(3))]});
            drawnow();
        end
    end
    set(handles.listbox1,'String',{['Adjusting brightness: ' num2str(dim(3)) ' / ' num2str(dim(3))]});
    drawnow();
    setappdata(handles.SliderFrames,'data_adj',data_adj);
    setappdata(handles.SliderFrames,'BrightMaxPre',BrightMaxVal);
    setappdata(handles.SliderFrames,'BrightMinPre',BrightMinVal);
else
    data_adj=getappdata(handles.SliderFrames,'data_adj');
    dim = size(data_adj);
end

FrameNo = floor(get(handles.SliderFrames,'Value'));
data_adj=data_adj(:,:,FrameNo);
axes(handles.Preview)
H1=getappdata(handles.SliderFrames,'H1');
H2=getappdata(handles.SliderFrames,'H2');
% The first two values of positions define the first point. The last two
% values define the third point of the rectangle
Pos1=floor(getPosition(H1)); Pos1(3:4)=Pos1(3:4)+Pos1(1:2); Pos1(Pos1<1)=1;
Pos2=floor(getPosition(H2)); Pos2(3:4)=Pos2(3:4)+Pos2(1:2); Pos2(Pos2<1)=1;
if Pos1(3)>dim(2); Pos1(3) = dim(2); end; if Pos1(4)>dim(1); Pos1(4) = dim(1); end
if Pos2(3)>dim(2); Pos2(3) = dim(2); end; if Pos2(4)>dim(1); Pos2(4) = dim(1); end
setappdata(handles.SliderFrames,'MaskPositions',[Pos1; Pos2]);
BW1=uint8(zeros(dim(1),dim(2),3));
BW2=uint8(zeros(dim(1),dim(2),3));
BW1(:,:,2)=uint8(createMask(H1)).*data_adj;
BW2(:,:,1)=uint8(createMask(H2)).*data_adj;

Mask = zeros(dim(1),dim(2),2);
Mask(:,:,1)=createMask(H1);
Mask(:,:,2)=createMask(H2);
setappdata(handles.SliderFrames,'Mask',Mask);

setColor(H1,'green');
setColor(H2,'red');
%assignin('base','BW1',BW1);
%assignin('base','BW2',BW2);

axes(handles.Preview);
%data_adj=imadjust(BW1+BW2,[BrightMinVal BrightMaxVal],[]);
data_adj=BW1+BW2;
H3=imshow(data_adj);
%setappdata(handles.SliderFrames,'ColoredImage',data_adj);

% Plot histogram
histData=get(handles.SliderBrightMin,'UserData');
datalength=(length(histData)-1)/2;
C=histData(1);
x=histData(2:(datalength+1));
counts=histData((datalength+2):end);
axes(handles.Hist);
hold off;
bar(x,counts,'blue');
xlim([0 1]);
ylim([min(counts) max(counts)]);
hold on;
plot([BrightMinVal BrightMaxVal],[0 C],'black');
set(handles.Hist,'yticklabel',[]);

assignin('base','Pos1',Pos1);
assignin('base','Pos2',Pos2);

if exist([FileInput 'Green.tif'],'file') || exist([FileInput 'Red.tif'],'file')
    button1 = questdlg([FileInput 'Green.tif or ' FileInput 'Red.tif or both already exist. Overwrite?']);
    if strcmp(button1,'Yes')
        delete([FileInput 'Green.tif']);
        delete([FileInput 'Red.tif']);
    end
end

% Import tiff and save
for i=1:numel(info)
    data = imread(FullPath, i, 'Info', info);
    assignin('base','data',data);
    dataGreen = data(Pos1(2):Pos1(4),Pos1(1):Pos1(3));
    dataRed = data(Pos2(2):Pos2(4),Pos2(1):Pos2(3));

    imwrite(dataGreen, [FileInput 'Green.tif'], 'WriteMode', 'append',  'Compression','none');
    imwrite(dataRed, [FileInput 'Red.tif'], 'WriteMode', 'append',  'Compression','none');
    if mod(i,10)==0
        set(handles.listbox1,'String',{['Saving split images: ' num2str(i) ' / ' num2str(numel(info))]});
        drawnow();
    end
end
set(handles.listbox1,'String',{['Saving split images: ' num2str(i) ' / ' num2str(numel(info))],['Images saved as ' FileInput 'Green.tif and ' FileInput 'Red.tif']});
drawnow();

% --- Executes on button press in DeletePoint.
function DeletePoint_Callback(hObject, eventdata, handles)
% hObject    handle to DeletePoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in FindPoints.
function FindPoints_Callback(hObject, eventdata, handles)
% hObject    handle to FindPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get FilePath and FileNo
FileInput=getappdata(handles.SliderFrames,'FileInput');
FileNo = floor(get(handles.SliderDrift,'Value'));

% List of parameters
CCDsens = 10.82;            % CCD sensitivity of the camera at specific readout rate and pre-amp setting. See note below for values
EMgain = 20;                % Electron multiplying (EM) gain setting of camera during acquisition
spacing = 20;               % Spacing between grid points, in pixels
CountToPhoton = CCDsens/EMgain;         % Count to photon conversion
CenterOnly = 1;                         % If 1, FIONA will only return center of points
DistanceThreshold = 3.5;                % Distance threshold in pixel to determine if two points in different channels matches
OverlapThreshold = 1;                   % The number of pixels allowed for overlaps. Sometimes FIONA analysis may be more accurate with bigger region.

BrightMaxVal = get(handles.SliderBrightMax,'Value');
BrightMaxVal = floor(BrightMaxVal*100)/100;
BrightMinVal = get(handles.SliderBrightMin,'Value');
BrightMinVal = floor(BrightMinVal*100)/100;
BrightMaxPre = getappdata(handles.SliderFrames,'BrightMaxPre');
BrightMinPre = getappdata(handles.SliderFrames,'BrightMinPre');
data = getappdata(handles.SliderFrames,'data_DriftBig');
data_Drift = getappdata(handles.SliderFrames,'data_Drift');
data_adj = getappdata(handles.SliderFrames,'data_adj');
%assignin('base','data',data);
dim = size(data); if numel(dim)==2; dim = [dim 1]; end

% Change data_adj if the brightness is changed
if BrightMaxVal~=BrightMaxPre || BrightMinVal~=BrightMinPre
    if isempty(data_Drift)
        FullPath = getappdata(handles.SliderFrames,'FullPath');
        info = imfinfo(FullPath);
        dim = [info(1).Height info(1).Width numel(info)];
        data_adj = uint8(zeros(dim));
        for i = 1:dim(3)
            data = imread(FullPath, i, 'Info', info);
            data_adj(:,:,i)=round(imadjust(data,[BrightMinVal BrightMaxVal],[])*(255/65535));
            %data_adj(:,:,i)=imadjust(data(:,:,i),[BrightMinVal BrightMaxVal],[]);
            if mod(i,10)==0
                set(handles.listbox2,'String',{['Adjusting brightness: ' num2str(i) ' / ' num2str(dim(3))]});
                drawnow();
            end
        end
    else
        dim = size(data_Drift);
        data_adj = uint8(zeros(dim));
        for i = 1:dim(3)
            data_adj(:,:,i)=round(imadjust(data_Drift(:,:,i),[BrightMinVal BrightMaxVal],[]));
        end
    end
    set(handles.listbox2,'String',{['Adjusting brightness: ' num2str(dim(3)) ' / ' num2str(dim(3))]});
    drawnow();
    setappdata(handles.SliderFrames,'data_adj',data_adj);
    setappdata(handles.SliderFrames,'BrightMaxPre',BrightMaxVal);
    setappdata(handles.SliderFrames,'BrightMinPre',BrightMinVal);
end

% Get radius
Radius=get(handles.SliderRadius,'Value');
RadiusThresholded = Radius - OverlapThreshold;
setappdata(handles.SliderFrames,'Radius',Radius);

% Initialize TotalFrameAnalyzed and set ProgressStep size
MaxOutput = 0;
TotalFrameAnalyzed = 0;
ProgressStep = 10;                      % How frequent progress is updated (in frames)

% Initialize UnmappedOutput
%UnmappedOutput=zeros(dim(3));

tic;
% Loop for each frame
UnmappedOutput = cell(dim(3)+2,1);
for loop = 0:dim(3)+1    
    if loop == 0; i=1;
    elseif loop == dim(3)+1; i=loop-1;
    else i=loop; end
    
    % Convert data into binary image (ones and zeros)
    assignin('base','data_adj',data_adj);
    Channel1 = im2bw(data_adj(:,:,i));

    % Draw Channel
    if 1
        axes(handles.Preview); hold off;
        imshow(double(Channel1)); hold on; drawnow();
    end

    % Label connected components
    [Labels1, num1] = bwlabel(Channel1, 8);
    Labels1=uint16(Labels1);

    % Find the centroid of each connected region. Round and report result
    % in matrix form
    Cent1 = regionprops(Channel1,'Centroid'); Cent1 = round(cell2mat({Cent1.Centroid}'));
    assignin('base','Cent1',Cent1);

    % Find box ranges for centroids within acceptable overlap defined
    % by the OverlapThreshold
    Centroid1 = [Cent1(:,1)-RadiusThresholded Cent1(:,1)+RadiusThresholded Cent1(:,2)-RadiusThresholded Cent1(:,2)+RadiusThresholded];
    Centroid1(Centroid1<1)=1; Centroid1(Centroid1(:,2)>dim(2),2)=dim(2); Centroid1(Centroid1(:,4)>dim(1),4)=dim(1);

    % Initialize SquareBoxes and LabeledBoxes matrix. SquareBoxes will
    % contain 0, 1, and more than ones for overlapping boxes.
    % LabeledBoxes will contain boxes labeled from 1 to the length of
    % the Centroids.
    SquareBoxes1=uint8(zeros(dim(1),dim(2),3)); LabeledBoxes1=uint16(zeros(dim(1),dim(2)));
    assignin('base','SquareBoxes1',SquareBoxes1); assignin('base','Centroid1',Centroid1);

    for l=1:length(Centroid1);
        assignin('base','l',l);
        SquareBoxes1(Centroid1(l,3):Centroid1(l,4),Centroid1(l,1):Centroid1(l,2),3)=SquareBoxes1(Centroid1(l,3):Centroid1(l,4),Centroid1(l,1):Centroid1(l,2),3)+1;
        LabeledBoxes1(Centroid1(l,3):Centroid1(l,4),Centroid1(l,1):Centroid1(l,2))=l;
    end
    %assignin('base','LabeledBoxes1',LabeledBoxes1); %assignin('base','SquareBoxes1',SquareBoxes1);

    % Check for overlaps between boxes to discount boxes that overlap.
    % This is done by looking at SquareBoxes. Pixels that have values
    % more than 1 means there is overlap, and this will be omitted. By
    % comparing BWBoxes (only 1 and 0) and GrayBoxes (can be more than 1),
    % we can find out which boxes overlap. FinalBoxes is a 2D SquareBoxes
    % without overlapping boxes
    GrayBoxes1=SquareBoxes1(:,:,3);
    BWBoxes1=im2bw(GrayBoxes1,0);
    FinalBoxes1=single(zeros(dim(1),dim(2)));
    [LBWBoxes1, numBWBoxes1] = bwlabel(BWBoxes1, 4);

    for j=1:numBWBoxes1
        Index=find(LBWBoxes1==j);
        FinalBoxes1(Index)=mean(GrayBoxes1(Index));
    end

    % Increase box size of Centroids so that they are the size defined
    % by the Radius. This simply increases the size by the overlap threshold
    Centroid1 = [Cent1(:,1)-Radius Cent1(:,1)+Radius Cent1(:,2)-Radius Cent1(:,2)+Radius];
    Centroid1(Centroid1<1)=1; Centroid1(Centroid1(:,2)>dim(2),2)=dim(2); Centroid1(Centroid1(:,4)>dim(1),4)=dim(1);

    % Omit out Centroids that are overlapping (the index of the labeled
    % boxes that has its corresponding FinalBoxes value having a mean
    % of more than 1 will tell us which Centroids are overlapping)
    Centroid1(unique(LabeledBoxes1(FinalBoxes1>1)),:)=[];

    % Draw SquareBoxes with bigger sized square, reflecting the value of Radius
    SquareBoxes1=uint8(zeros(dim(1),dim(2),3));
    for l=1:length(Centroid1);
        SquareBoxes1(Centroid1(l,3):Centroid1(l,4),Centroid1(l,1):Centroid1(l,2),3)=1;
    end
    SquareBoxes1(:,:,1)=SquareBoxes1(:,:,3);

    % Draw SquareBoxes
    if 1; h=imshow(double(SquareBoxes1)); set(h,'AlphaData',0.3); end

    % Find the number of points that do not overlap with one another
    numFin1 = length(Centroid1);

    % Calculate the number of overlapping points, then display on listbox
    OverlapPoints1=num1-numFin1;
    String1=['Frame Analyzed: ' num2str(i)]; 
    String2=['Points found: ' num2str(numFin1)];
    String3=['Overlapping points: ' num2str(OverlapPoints1)];
    set(handles.listbox2,'String',{String1;String2;String3});

    % Do FIONA on each point
    nCentroid1 = length(Centroid1);

    % Initialize Output
    if CenterOnly == 1; Output1 = zeros(nCentroid1,2);
    else Output1 = zeros(nCentroid1,7); end
    %assignin('base','data',data); %assignin('base','Centroid1',Centroid1); %assignin('base','Centroid2',Centroid2);      

    % Looping for FIONA 
    for l=1:length(Centroid1)
        assignin('base','data',data);
        assignin('base','l',l);
        assignin('base','i',i);
        assignin('base','Centroid1',Centroid1);

        % Get the information of the first and last data point when
        % i=dim(3)+1 and i=dim(3)+2
        if loop==0;
            CodePath = get(handles.FileLocation,'string');
            Folder = get(handles.FileName,'string');
            DataPath = [CodePath '\' Folder];
            cd(DataPath);
            infoTemp = imfinfo(FileInput{FileNo});
            setappdata(handles.SliderFrames,'TotalFrame',numel(infoTemp));
            dataTemp = imread(FileInput{FileNo}, 1, 'Info', infoTemp);
            FIONAInput = double(dataTemp(Centroid1(l,3):Centroid1(l,4),Centroid1(l,1):Centroid1(l,2)));
            cd(CodePath);
        elseif loop == dim(3)+1;
            cd(DataPath);
            dataTemp = imread(FileInput{FileNo}, numel(infoTemp), 'Info', infoTemp);
            FIONAInput = double(dataTemp(Centroid1(l,3):Centroid1(l,4),Centroid1(l,1):Centroid1(l,2)));
            cd(CodePath);
        else
            FIONAInput = double(data(Centroid1(l,3):Centroid1(l,4),Centroid1(l,1):Centroid1(l,2),i));
        end

        [ny,nx]=size(FIONAInput);               % Find out the dimensions
        grid = [nx ny 1:nx 1:ny];         % Gridding input for gauss2dfunct and gauss2dfit
        tilt=0;
        % Parameters: p(1): z-offset, p(2): amplitude, p(3): xStdev, p(4): yStdev, p(5): xCenter, p(6): yCenter, p(7): tilt.
        try
            [popt,resnorm,residual,ret]=gauss2dfit(FIONAInput,grid,tilt);
            
            % Check if center is out of the box dimension
            if popt(5)>0 && popt(5)<nx && popt(6)>0 && popt(6)<ny
                % Getting center and precision
                xCenter = (Centroid1(l,1)-1+popt(5));        % Center of x in pixels
                yCenter = (Centroid1(l,3)-1+popt(6));        % Center of y in pixels
                sx2 = popt(3)*popt(3);              % Square of xStdev
                sx4 = sx2*sx2;                      % xStdev^4
                sy2 = popt(4)*popt(4);              % Square of yStdev
                sy4 = sy2*sy2;                      % yStdev^4

                % To estimate b, which is the standard deviation of the background,
                % we'll look at the z-offset (popt(1)) and calculate the standard
                % deviation based on anything below the z-offset. Before that we
                % would want to make everything 3 standard deviations away from
                % our spot to be (z-offset + 1) or NaN.

                % Find the limits of the data
                xmin = floor(popt(5)-4*popt(3)); xmax = ceil(popt(5)+4*popt(3));
                ymin = floor(popt(6)-4*popt(4)); ymax = ceil(popt(6)+4*popt(4));

                % Correct limits when they are out of boundaries
                if xmin < 1; xmin = 1; elseif xmax > nx; xmax = nx;
                elseif ymin < 1; ymin = 1; elseif ymax > ny; ymax = ny; end

                % Calculate PhotonNo and precisions if more output is needed
                if CenterOnly == 1
                    Output1(l,:) = [xCenter yCenter];
                else
                    z = reshape(FIONAInput,ny*nx,1); 
                    SumSquaresTotal = sum(z.*z);
                    RSquare = 1-(resnorm/SumSquaresTotal);
                    residual = reshape(residual,ny,nx);     % Extract residual
                    residual = -residual;                   % Invert such that residual = FIONAInputl - zfit
                    residual(xmin:xmax,ymin:ymax)=0;        % Set values of residual around center to be NaN
                    residual = residual(residual<0);
                    b = sqrt(sum(residual.*residual)/(length(residual)-1));
                    b = b*CountToPhoton;                    % Standard deviation of the background
                    b2 = b*b;                               % Square of background
                    PhotonNo = abs(2*pi*popt(2)*popt(3)*popt(4)*CountToPhoton);   % Number of Photons calculated using volume under gaussian, which is 2*pi*A*stdev(x)*stdev(y)
                    PhotonNo2 = PhotonNo*PhotonNo;          % Square of PhotonNo
                    xPrecision = sqrt((sx2/PhotonNo) + (a2/(12*PhotonNo)) + (8*pi*sx4*b2/(a2*PhotonNo2)));
                    yPrecision = sqrt((sy2/PhotonNo) + (a2/(12*PhotonNo)) + (8*pi*sy4*b2/(a2*PhotonNo2)));
                    Output1(l,:) = [PhotonNo xPrecision yPrecision xCenter yCenter RSquare ErrorStatus];
                end
                
                % Define no error status
                ErrorStatus = 0;
            else
                % Define error status 1 if center is beyond box dimension
                ErrorStatus = 1;
                if CenterOnly == 1            
                    if OutputNo == 1, Output1(l,:) = [0 0];
                    else Output2(l,:) = [0 0];
                    end
                else
                    if OutputNo == 1, Output1(l,:) = [0 0 0 0 0 0 ErrorStatus];
                    else Output2(l,:) = [0 0 0 0 0 0 ErrorStatus];
                    end
                end
            end
        catch
            % Define error status 2 if FIONA did not complete
            ErrorStatus = 2;
            if CenterOnly == 1            
                Output1(l,:) = [0 0];
            else
                Output1(l,:) = [0 0 0 0 0 0 ErrorStatus];
            end
        end

        % Display progress
        if mod(l,ProgressStep)==0
            TotalFrameAnalyzed = TotalFrameAnalyzed + ProgressStep;
            String1 = [' Point Analyzed :             ' num2str(l) ' / ' num2str(length(Centroid1))];
            String2 = [' Frame Analyzed :           ' num2str(i) '/ ' num2str(dim(3))];
            String3 = [' Elapsed time:                   ' num2str(floor(toc)) ' s'];
            String4 = [' Total points analyzed:    ' num2str(TotalFrameAnalyzed)];
            set(handles.listbox2,'String',{' Analyzing...', String1, String2, String3, String4});
            drawnow();
        end
    end
    TotalFrameAnalyzed = TotalFrameAnalyzed + mod(l,ProgressStep);
    if length(Centroid1)>MaxOutput
        MaxOutput = length(Centroid1);
    end
    
    % Save centroid information in UnmappedOutput
    if loop==0; UnmappedOutput(1,1)={Output1};
    elseif loop == dim(3)+1; UnmappedOutput(dim(3)+2,1)={Output1};
    else UnmappedOutput(loop+1,1)={Output1};
    end

    % Plotting Centroids
    if 1
        plot(Output1(:,1), Output1(:,2), 'b*'); hold off; drawnow();
        %assignin('base','Output1',Output1); %assignin('base','Output2',Output2);
    end

end
assignin('base','UnmappedOutput',UnmappedOutput);

% Update listbox and show that analysis is complete
String1 = [' Point Analyzed :             ' num2str(l) ' / ' num2str(length(Centroid1))];
String2 = [' Frame Analyzed :           ' num2str(i) '/ ' num2str(dim(3))];
String3 = [' Elapsed time:                   ' num2str(floor(toc)) ' s'];
String4 = [' Total points analyzed:    ' num2str(TotalFrameAnalyzed)];
set(handles.listbox2,'String',{' Analyzing...', String1, String2, String3, String4, ' ', ' Analysis Complete!'});
drawnow();

% Find out the nearest neighbor between adjacent frames
TrackingPath = zeros(MaxOutput,dim(3)+1);
for i=1:dim(3)+1
    [NeighborInd,Distance] = knnsearch(UnmappedOutput{i+1}, UnmappedOutput{i});
    NeighborInd(Distance>DistanceThreshold)=0;
    TrackingPath(1:length(NeighborInd),i)=NeighborInd;
end
assignin('base','TrackingPath',TrackingPath);

% Initialize StartPath. This matrix determines whether we will start
% finding a track path from one particular point. If that point has been
% used in another track, the StartPath value for that point will be zero,
% and a search will not start from that point. That point can still be
% part of another track, just that it can no longer be a start to a track
StartPath = ones(size(TrackingPath));
StartPath(TrackingPath==0)=0;

if 1
TrackCount = 1;
for Col = 1:size(TrackingPath,2);
    for Row = 1:size(TrackingPath,1);
        if StartPath(Row,Col) > 0
            TrackCol = Col;
            TrackRow = Row;
            Track = zeros(dim(3),3);
            while TrackCol > 0
                if TrackRow > 0 && TrackCol <= size(TrackingPath,2)
                    Track(TrackCol,:)=[UnmappedOutput{TrackCol}(TrackRow,:) TrackCol];
                    StartPath(TrackRow,TrackCol) = 0;
                    TrackRow = TrackingPath(TrackRow,TrackCol);
                    TrackCol = TrackCol + 1;
                elseif TrackRow > 0 && TrackCol == size(TrackingPath,2)+1
                    Track(TrackCol,:)=[UnmappedOutput{TrackCol}(TrackRow,:) TrackCol];
                    TrackCol = 0; Track(Track(:,1)==0,:)=[];
                    if ~isempty(Track)
                        FinalTrack(TrackCount,1)={Track};
                        TrackCount = TrackCount + 1;
                    end
                else
                    TrackCol = 0; Track(Track(:,1)==0,:)=[];
                    if ~isempty(Track)
                        FinalTrack(TrackCount,1)={Track};
                        TrackCount = TrackCount + 1;
                    end
                end
            end
        end
    end
end
end

assignin('base','FinalTrack',FinalTrack);
% Plot tracks
TrackCount = length(FinalTrack);
ColorSet = hsv(TrackCount);
axes(handles.Preview); hold off;
plot(FinalTrack{1}(:,1),FinalTrack{1}(:,2),'Color', ColorSet(1,:)); hold on;
for i = 2:TrackCount
    plot(FinalTrack{i}(:,1),FinalTrack{i}(:,2),'Color', ColorSet(i,:));
end

% FinalTrack2 substract the first point of every track to the rest
% of the track. It also exclude those tracks shorter or equal to LengthThreshold
LengthThreshold = dim(3)+1;
assignin('base','LengthThreshold',LengthThreshold);
TrackNo = 1;
axes(handles.Preview); hold off;
x=FinalTrack{1}(:,1)-FinalTrack{1}(1,1); y=FinalTrack{1}(:,2)-FinalTrack{1}(1,2);
FinalTrack2{1}=[x y FinalTrack{1}(:,3)];
plot(x,y,'Color', ColorSet(1,:)); hold on;
for i = 2:TrackCount
    if length(FinalTrack{i})>LengthThreshold
        x=FinalTrack{i}(:,1)-FinalTrack{i}(1,1); y=FinalTrack{i}(:,2)-FinalTrack{i}(1,2);
        FinalTrack2(TrackNo,1)={[x y FinalTrack{i}(:,3)]};
        TrackNo = TrackNo + 1;
        plot(x,y,'Color', ColorSet(i,:));
    end
end
assignin('base','FinalTrack2',FinalTrack2);

% Calculate distance difference between sets of data
TrackCount = length(FinalTrack2);
RSquareMatrix = zeros(TrackCount,TrackCount);
x0=0;
for i = 1:TrackCount
    for j = i:TrackCount
        StartInd=max([FinalTrack2{i}(1,3) FinalTrack2{j}(1,3)]);
        EndInd=min([FinalTrack2{i}(end,3) FinalTrack2{j}(end,3)]);
        if EndInd-StartInd >= LengthThreshold
            data1=FinalTrack2{i}(find(FinalTrack2{i}(:,3)==StartInd):find(FinalTrack2{i}(:,3)==EndInd),1:2); 
            data2=FinalTrack2{j}(find(FinalTrack2{j}(:,3)==StartInd):find(FinalTrack2{j}(:,3)==EndInd),1:2);
            %[popt,resnorm,residual,~,output] = lsqcurvefit(@gauss2dfunct,p0,grid,z,[],[],options);
            [popt,resnormx,residual,~,output] = lsqcurvefit(@AlignData,x0,data1(:,1),data2(:,1));
            [popt,resnormy,residual,~,output] = lsqcurvefit(@AlignData,x0,data1(:,2),data2(:,2));
            SumSquaresTotalx = sum(data1(:,1).*data1(:,1));
            SumSquaresTotaly = sum(data1(:,2).*data1(:,2));
            RSquare = 1-((resnormx/SumSquaresTotalx)+(resnormy/SumSquaresTotaly));
            RSquareMatrix(i,j)=RSquare;
        else
            RSquareMatrix(i,j)=0;
        end
    end
end
RSquareMatrix(isnan(RSquareMatrix))=0;
RSquareMatrix(RSquareMatrix<0)=0;
RSquareMatrix=RSquareMatrix+RSquareMatrix';
RSquareMatrix(RSquareMatrix>=1)=0;
RSquareMax=max(RSquareMatrix,[],2);
%figure;plot(sort(RSquareMax));
%RSquareMax=mean(RSquareMatrix);
%HistRSquare = reshape(RSquareMatrix,TrackCount^2,1);
%HistRSquare(HistRSquare==0)=[];
%[HistRSquare,xout]=hist(HistRSquare,100);
%HistRSquare=sort(HistRSquare);

% Plot tracks that have strong correlation with the others
RSquareThreshold = get(handles.sliderThreshold,'Value');;
ColorSet = hsv(TrackCount);
axes(handles.Preview); hold off;
ChosenTrack = zeros(length(TrackCount),1);
%figure;
for i = 1:TrackCount
    if RSquareMax(i)>RSquareThreshold
        x=FinalTrack2{i}(:,1); y=FinalTrack2{i}(:,2); z=ones(length(x),1)*i;
        plot3(x,y,z, 'Color', ColorSet(i,:)); hold on; view(2);
        ChosenTrack(i) = 1;
    end
end
datacursormode on;
setappdata(handles.SliderFrames,'ChosenTrack',ChosenTrack);
setappdata(handles.SliderFrames,'RSquareMax',RSquareMax);
setappdata(handles.SliderFrames,'FinalTrack2',FinalTrack2);
%setappdata(handles.SliderFrames,'RSquareMatrix',RSquareMatrix);

% --- Executes on button press in SingleFile.
function SingleFile_Callback(hObject, eventdata, handles)
% hObject    handle to SingleFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of SingleFile

set(handles.SingleFile,'Value',1);
set(handles.AllFiles,'Value',0);

% --- Executes on button press in AllFiles.
function AllFiles_Callback(hObject, eventdata, handles)
% hObject    handle to AllFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of AllFiles

set(handles.SingleFile,'Value',0);
set(handles.AllFiles,'Value',1);

function TextDrift_Callback(hObject, eventdata, handles)
% hObject    handle to TextDrift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of TextDrift as text
%        str2double(get(hObject,'String')) returns contents of TextDrift as a double

% --- Executes during object creation, after setting all properties.
function TextDrift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TextDrift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function SliderDrift_Callback(hObject, eventdata, handles)
% hObject    handle to SliderDrift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Get CodePath
CodePath = get(handles.FileLocation,'string');

% Get DataPath
Folder = get(handles.FileName,'string');
% Check if Folder is empty. if it not, check if the directory exists
if isempty(Folder)
    % Ask for data path if Folder is empty
    DataPath = uigetdir;
else
    DataPath = [CodePath '\' Folder];
    % If directory does not exist, get directory
    if ~exist(DataPath,'dir'); DataPath = uigetdir; end
end
name=strrep(DataPath,[CodePath '\'],'');
set(handles.FileName,'String',name);
cd(DataPath);

% Get the FileNo from SliderDrift
FileIn = dir('*.tif');
if length(FileIn) > 1
    FileNo = floor(get(handles.SliderDrift,'Value'));
else
    FileNo = 1;
end

% Adjust the value of TextDrift
set(handles.TextDrift,'String',FileNo);

% Get FileInput
FileInput=getappdata(handles.SliderFrames,'FileInput');
assignin('base','FileInput',FileInput);

% Get Tiff info and Sum Factor
info = imfinfo(FileInput{FileNo});
num_images = numel(info); 
FinalFrameNo = str2double(get(handles.TextFinalFrame,'string'));
SumFactor = floor(num_images/FinalFrameNo);
%SumRemain = num_images - SumFactor*FinalFrameNo;

% Sum images to FinalFrameNo
dim = [info(1).Height info(1).Width FinalFrameNo];
data = double(zeros(dim(1),dim(2),FinalFrameNo));
data_adj = uint8(zeros(dim));
for k = 0:FinalFrameNo-1
    for j = 1:SumFactor
        data(:,:,k+1) = data(:,:,k+1)+double(imread(FileInput{FileNo}, (SumFactor*k)+j, 'Info', info));
    end
    
    if mod(k,1)==0
        set(handles.listbox2,'String',{['Summing images: ' num2str(k+1) ' / ' num2str(FinalFrameNo)]});
        drawnow();
    end
    %imwrite(FinalFrame, [FileInputName{ind} '(' num2str(SumFactor) 'x).tif'], 'WriteMode', 'append',  'Compression','none');
end
setappdata(handles.SliderFrames,'data_DriftBig',data);
data = uint8(data/(SumFactor*256));
%assignin('base','data',data);

% Read first frame
[counts,x]=imhist(data(:,:,1));
[C,I]=max(counts);
range = max(x)-min(x);
x = (x-min(x))/range;
set(handles.SliderBrightMin,'UserData',[C;x;counts]);
BrightMinVal=ceil(x(I)*2*100)/100;
set(handles.SliderBrightMin,'Value',BrightMinVal);
set(handles.TextBrightMin,'String',BrightMinVal);
CumulativeSum=cumsum(counts);
BrightMaxVal=floor(0.99*sum(counts));
BrightMaxVal=floor(x((find(CumulativeSum>BrightMaxVal, 1 ))*2)*100)/100;
if BrightMinVal>=BrightMaxVal; BrightMaxVal=BrightMinVal+0.01; end;
set(handles.SliderBrightMax,'Value',BrightMaxVal);
set(handles.TextBrightMax,'String',BrightMaxVal);
axes(handles.Hist); hold off;
bar(x,counts,'blue'); hold on;
plot([0 1],[0 C],'black'); axis off;
set(handles.Hist,'yticklabel',[]);
xlim([0 1]); ylim([min(counts) max(counts)]);

for i=1:FinalFrameNo
    data_adj(:,:,i)=round(imadjust(data(:,:,i),[BrightMinVal BrightMaxVal],[]));
end

% Plot Image
axes(handles.Preview);
imshow(data_adj(:,:,1));
setappdata(handles.SliderFrames,'data_adj',data_adj);
setappdata(handles.SliderFrames,'data_Drift',data);

% Set the value of slider.Frames and text.Frames to be 1
set(handles.SliderFrames,'Value',1);
set(handles.TextFrames,'String',1);

% Set properties of SliderFrames
set(handles.SliderFrames,'Max',FinalFrameNo);
set(handles.SliderFrames,'SliderStep',[1/(FinalFrameNo-1) 1/(FinalFrameNo-1)]);

% Initialize value of FrameNo
setappdata(handles.SliderFrames,'FrameNoPre',1);

% Initialize values of BrightMaxPre and BrightMinPre
setappdata(handles.SliderFrames,'BrightMaxPre',BrightMaxVal);
setappdata(handles.SliderFrames,'BrightMinPre',BrightMinVal);

%assignin('base','data_adj',data_adj);

% Initialize values of LimitPre
xlimit = floor(xlim);
ylimit = floor(ylim);
setappdata(handles.SliderFrames,'LimitPre',[xlimit ylimit]);

% Initialize the ImageSaved counter
ImageSaved = 0;
setappdata(handles.SliderFrames,'ImageSaved',ImageSaved);

% Initialize/ Reset values of FileName
setappdata(handles.SliderFrames,'FileName','');
setappdata(handles.SliderFrames,'ROI',[]);
setappdata(handles.SliderFrames,'ROIPos',[]);
setappdata(handles.SliderFrames,'HROIred',[]);

% Reset listboxes
%set(handles.listbox1,'String',{''});
%FolderName = getappdata(handles.SliderFrames,'FolderName');
%FileName = getappdata(handles.SliderFrames,'FileName');
%ImageSaved = getappdata(handles.SliderFrames,'ImageSaved');
%String1 = ['Folder Name   : ' FolderName];
%String2 = ['File Name        : ' FileName];
%String3 = ['Image Saved   : ' num2str(ImageSaved)];
%String4 = ['Files saved as : '];
%set(handles.listbox2,'String',{String1;String2;String3;String4});

% Clear data_Preview and data_PreviewAdj
setappdata(handles.SliderFrames,'data_Preview',[]);
setappdata(handles.SliderFrames,'data_PreviewAdj',[]);

% Switch to CodePath
cd(CodePath);

% --- Executes during object creation, after setting all properties.
function SliderDrift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderDrift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function TextFinalFrame_Callback(hObject, ~, handles)
% hObject    handle to TextFinalFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of TextFinalFrame as text
%        str2double(get(hObject,'String')) returns contents of TextFinalFrame as a double

% --- Executes during object creation, after setting all properties.
function TextFinalFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TextFinalFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in DeleteTrack.
function DeleteTrack_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ChosenTrack = getappdata(handles.SliderFrames,'ChosenTrack');
RSquareMax = getappdata(handles.SliderFrames,'RSquareMax');
FinalTrack2 = getappdata(handles.SliderFrames,'FinalTrack2');

dcm_obj=datacursormode;
Cursor=getCursorInfo(dcm_obj);
ChosenTrack(Cursor.Position(3))=0;

TrackCount = length(FinalTrack2);
RSquareThreshold = get(handles.sliderThreshold,'Value');
ColorSet = hsv(TrackCount);
axes(handles.Preview); hold off;
%figure;
for i = 1:TrackCount
    if ChosenTrack(i)==1
        x=FinalTrack2{i}(:,1); y=FinalTrack2{i}(:,2); z=ones(length(x),1)*i;
        plot3(x,y,z, 'Color', ColorSet(i,:)); hold on; view(2);
    end
end
datacursormode on;

assignin('base','ChosenTrack',ChosenTrack);
setappdata(handles.SliderFrames,'ChosenTrack',ChosenTrack);

% --- Executes on button press in ComputeDrift.
function ComputeDrift_Callback(hObject, eventdata, handles)
% hObject    handle to ComputeDrift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ChosenTrack = getappdata(handles.SliderFrames,'ChosenTrack');
assignin('base','ChosenTrack',ChosenTrack);
RSquareMax = getappdata(handles.SliderFrames,'RSquareMax');
FinalTrack2 = getappdata(handles.SliderFrames,'FinalTrack2');
TotalFrame = getappdata(handles.SliderFrames,'TotalFrame');

% Get CodePath
CodePath = get(handles.FileLocation,'string');

% Get DataPath
Folder = get(handles.FileName,'string');
% Check if Folder is empty. if it not, check if the directory exists
if isempty(Folder)
    % Ask for data path if Folder is empty
    DataPath = uigetdir;
else
    DataPath = [CodePath '\' Folder];
    % If directory does not exist, get directory
    if ~exist(DataPath,'dir'); DataPath = uigetdir; end
end
cd(DataPath);

% Get FileInput and FileNo
FileInput=getappdata(handles.SliderFrames,'FileInput');
FileNo = floor(get(handles.SliderDrift,'Value'));

SumTrack = zeros(size(FinalTrack2{1}));
TrackCount = length(FinalTrack2);
axes(handles.Preview); hold off;
for i = 1:TrackCount
    if ChosenTrack(i)==1
        SumTrack = SumTrack + FinalTrack2{i};
    end
end
SumTrack = SumTrack/length(find(ChosenTrack==1));
assignin('base','SumTrack',SumTrack);
plot3(SumTrack(:,1),SumTrack(:,2),SumTrack(:,3)); view(2);hold on;

% Interpolate data
Step = floor(TotalFrame/(length(SumTrack)-2));
xOriginal = [0 (Step/2):Step:(TotalFrame-Step/2) TotalFrame];
assignin('base','xOriginal',xOriginal);
xInterpolation = interp1(xOriginal,SumTrack(:,1),0:TotalFrame,'spline')';
yInterpolation = interp1(xOriginal,SumTrack(:,2),0:TotalFrame,'spline')';
plot(xInterpolation,yInterpolation,'r');

% Save data
FileName=strrep(FileInput{FileNo},'.tif','');
dlmwrite([FileName '-Drift.txt'],[xInterpolation yInterpolation],'\t');

% Switch back to CodePath
cd(CodePath);

function TextThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to TextThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of TextThreshold as text
%        str2double(get(hObject,'String')) returns contents of TextThreshold as a double

% --- Executes during object creation, after setting all properties.
function TextThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TextThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function sliderThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to sliderThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

RSquareMax = getappdata(handles.SliderFrames,'RSquareMax');
FinalTrack2 = getappdata(handles.SliderFrames,'FinalTrack2');

TrackCount = length(FinalTrack2);
ChosenTrack = zeros(TrackCount,1);
RSquareThreshold = get(handles.sliderThreshold,'Value');
ColorSet = hsv(TrackCount);
axes(handles.Preview); hold off;
%figure;
for i = 1:TrackCount
    if RSquareMax(i)>RSquareThreshold
        x=FinalTrack2{i}(:,1); y=FinalTrack2{i}(:,2); z=ones(length(x),1)*i;
        plot3(x,y,z, 'Color', ColorSet(i,:)); hold on; view(2);
        ChosenTrack(i) = 1;
    end
end
datacursormode on;

setappdata(handles.SliderFrames,'ChosenTrack',ChosenTrack);
assignin('base','ChosenTrack',ChosenTrack);
assignin('base','FinalTrack2',FinalTrack2);
set(handles.TextThreshold,'String',num2str(RSquareThreshold));

% --- Executes during object creation, after setting all properties.
function sliderThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function keyPress(src,eventdata,handles)
    
% This is to allow us to use keypress function with datacursormode.
% Eventually settle with add.listener because we cannot seem to damp the
% windows listener. There is an error associated with this, which is that
% when we press the alt key needed to select another cursor, there is an
% error generated. We therefore use the try and catch function to silent
% this error
% We will no longer need the figure1_WindowsKeyPressFcn, so we can delete
% that

% Since hObject is not defined, we need to define it here, so that we can
% get the handles to access all the gui variables (getappdata, get, etc)
hObject=findobj('Type','figure','Tag','figure1');
delete(setdiff(findall(0, 'type', 'figure'), hObject));

if ~(exist('handles','var'))
     handles=guidata(hObject);
end

% Change focus to SliderFileNo. Otherwise the focus is on SearchText and it will
% be activated if we press spacebar
% uicontrol(handles.SliderFileNo);  % Now this is included in the if
% statement to prevent all keys to be directed to SliderFileNo

key = get(src,'CurrentCharacter');
if 1
    try
        if key == ' ' || key == 'a'    
            uicontrol(handles.SliderFileNo);
            FileNo = get(handles.SliderFileNo,'Value');
            num_images = get(handles.SliderFileNo,'Max');
            
            if key == ' '
                if FileNo~=num_images; FileNo = FileNo+1; end
            elseif key == 'a'        
                if FileNo~=1; FileNo = FileNo-1; end 
            end
            set(handles.SliderFileNo,'Value',FileNo);
            set(handles.TextFileNo,'String',FileNo);
            SliderFileNo_Callback(hObject, eventdata, handles);
            uicontrol(handles.SliderFileNo);
            figure(hObject);
        end       
        if key == 's'; uicontrol(handles.SliderFileNo); Left_Callback(hObject, eventdata, handles); uicontrol(handles.SliderFileNo); figure(hObject); end
        if key == 'f'; uicontrol(handles.SliderFileNo); Right_Callback(hObject, eventdata, handles); uicontrol(handles.SliderFileNo); figure(hObject); end
        if key == 'e'; uicontrol(handles.SliderFileNo); Up_Callback(hObject, eventdata, handles); uicontrol(handles.SliderFileNo); figure(hObject); end
        if key == 'd'; uicontrol(handles.SliderFileNo); Multiple_Callback(hObject, eventdata, handles); uicontrol(handles.SliderFileNo); figure(hObject); end
        if key == 'c'; uicontrol(handles.SliderFileNo); Down_Callback(hObject, eventdata, handles); uicontrol(handles.SliderFileNo); figure(hObject); end
        if key == 'r'; uicontrol(handles.SliderFileNo); Reorient_Callback(hObject, eventdata, handles); uicontrol(handles.SliderFileNo); figure(hObject); end
        %if key == 'r'; uicontrol(handles.SliderFileNo); NewVal=str2num(get(handles.TextNoiseRemove,'String'));NewVal=NewVal+1;set(handles.TextNoiseRemove,'String',num2str(NewVal)); uicontrol(handles.SliderFileNo); figure(hObject); end
        %if key == 'v'; uicontrol(handles.SliderFileNo); NewVal=str2num(get(handles.TextNoiseRemove,'String'));if NewVal>1;NewVal=NewVal-1;end;set(handles.TextNoiseRemove,'String',num2str(NewVal)); uicontrol(handles.SliderFileNo); figure(hObject); end
        %if key == 'w'; uicontrol(handles.SliderFileNo); NewVal=str2num(get(handles.TextFitWeight,'String'));NewVal=NewVal+5;set(handles.TextFitWeight,'String',num2str(NewVal)); uicontrol(handles.SliderFileNo); figure(hObject); end
        %if key == 'x'; uicontrol(handles.SliderFileNo); NewVal=str2num(get(handles.TextFitWeight,'String'));if NewVal>0;NewVal=NewVal-5;end;set(handles.TextFitWeight,'String',num2str(NewVal)); uicontrol(handles.SliderFileNo); end
        if key == 'z'; uicontrol(handles.SliderFileNo); SaveText_Callback(hObject, eventdata, handles); uicontrol(handles.SliderFileNo); figure(hObject); end
           
    catch
        %disp('There is an error');
        %rethrow(err);
    end
end
% Change focus to SliderFileNo. Otherwise the focus is on SearchText and it will
% be activated if we press spacebar

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in popupmenu4.

function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function textFIONA_Callback(hObject, eventdata, handles)
% hObject    handle to textFIONA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of textFIONA as text
%        str2double(get(hObject,'String')) returns contents of textFIONA as a double

% --- Executes during object creation, after setting all properties.
function textFIONA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textFIONA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ChangeSetting.
function ChangeSetting_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeSetting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get FIONASetting
FIONASetting = getappdata(handles.SliderFrames,'FIONASetting');

% Get the selection of the popup menu
PopupValue = get(handles.popupmenu1,'Value');
switch PopupValue
    case 1
        FIONASetting.ActPix = get(handles.textFIONA,'String');
    case 2
        FIONASetting.ObjMag = get(handles.textFIONA,'String');
    case 3
        FIONASetting.AddMag = get(handles.textFIONA,'String');
    case 4
        Selection = get(handles.popupmenu2,'String');
        FIONASetting.CameraVal = get(handles.popupmenu2,'Value');
        FIONASetting.Camera = Selection{FIONASetting.CameraVal};
        switch FIONASetting.CameraVal
            case 1
                switch FIONASetting.PreampVal
                    case 1
                        FIONASetting.Preamp = '  1x';
                    case 2
                        FIONASetting.Preamp = '  2.3x';
                    case 3
                        FIONASetting.Preamp = '  4.5x';
                end
            case 2
                switch FIONASetting.PreampVal
                    case 1
                        FIONASetting.Preamp = '  1x';
                    case 2
                        FIONASetting.Preamp = '  2.4x';
                    case 3
                        FIONASetting.Preamp = '  4.9x';
                end
            case 3
                switch FIONASetting.PreampVal
                    case 1
                        FIONASetting.Preamp = '  1x';
                    case 2
                        FIONASetting.Preamp = '  2.5x';
                    case 3
                        FIONASetting.Preamp = '  5.2x';
                end
        end
    case 5
        Selection = get(handles.popupmenu2,'String');
        FIONASetting.ReadoutVal = get(handles.popupmenu2,'Value');
        FIONASetting.Readout = Selection{FIONASetting.ReadoutVal};
    case 6
        Selection = get(handles.popupmenu2,'String');
        FIONASetting.PreampVal = get(handles.popupmenu2,'Value');
        FIONASetting.Preamp = Selection{FIONASetting.PreampVal};
    case 7
        FIONASetting.EMgain = get(handles.textFIONA,'String');
    case 8
        Selection = get(handles.popupmenu2,'String');
        FIONASetting.ParallelVal = get(handles.popupmenu2,'Value');
        FIONASetting.Parallel = Selection{FIONASetting.ParallelVal};
end

switch FIONASetting.CameraVal
    case 1
        switch FIONASetting.ReadoutVal
            case 1
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '24.33';
                    case 2; FIONASetting.CCDSens = '9.74';
                    case 3; FIONASetting.CCDSens = '4.26';
                end
            case 2
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '59.84';
                    case 2; FIONASetting.CCDSens = '24.17';
                    case 3; FIONASetting.CCDSens = '10.72';
                end
            case 3
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '59.86';
                    case 2; FIONASetting.CCDSens = '24.36';
                    case 3; FIONASetting.CCDSens = '10.82';
                end
            case 4
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '67.42';
                    case 2; FIONASetting.CCDSens = '26.85';
                    case 3; FIONASetting.CCDSens = '12.27';
                end
        end
    case 2
        switch FIONASetting.ReadoutVal
            case 1
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '23.3';
                    case 2; FIONASetting.CCDSens = '9.1';
                    case 3; FIONASetting.CCDSens = '4.2';
                end
            case 2
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '55.8';
                    case 2; FIONASetting.CCDSens = '23.6';
                    case 3; FIONASetting.CCDSens = '10.3';
                end
            case 3
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '55.0';
                    case 2; FIONASetting.CCDSens = '23.9';
                    case 3; FIONASetting.CCDSens = '10.4';
                end
            case 4
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '64.5';
                    case 2; FIONASetting.CCDSens = '26.3';
                    case 3; FIONASetting.CCDSens = '11.9';
                end
        end
    case 3
        switch FIONASetting.ReadoutVal
            case 1
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '23.02';
                    case 2; FIONASetting.CCDSens = '9.29';
                    case 3; FIONASetting.CCDSens = '4.18';
                end
            case 2
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '56.17';
                    case 2; FIONASetting.CCDSens = '23.07';
                    case 3; FIONASetting.CCDSens = '10.05';
                end
            case 3
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '57.20';
                    case 2; FIONASetting.CCDSens = '23.20';
                    case 3; FIONASetting.CCDSens = '10.33';
                end
            case 4
                switch FIONASetting.PreampVal
                    case 1; FIONASetting.CCDSens = '67.75';
                    case 2; FIONASetting.CCDSens = '27.00';
                    case 3; FIONASetting.CCDSens = '12.13';
                end
        end
end
            
String1 = ['Actual Pixel :                        ' FIONASetting.ActPix ' nm'];
String2 = ['Objective Magnification :      ' FIONASetting.ObjMag];
String3 = ['Additional Magnification :     ' FIONASetting.AddMag];
String4 = ['Camera :                            ' FIONASetting.Camera];
String5 = ['Readout Rate :                   ' FIONASetting.Readout];
String6 = ['Preamp Setting :                 ' FIONASetting.Preamp];
String7 = ['EM gain :                               ' FIONASetting.EMgain];
String8 = ['Parallel Processing:            ' FIONASetting.Parallel];
String9 = ['CCD Sensitivity :                   ' FIONASetting.CCDSens];

set(handles.listbox3,'String',{String1;String2;String3;String4;String5;String6;String7;String8;String9});

% Store FIONASetting variable in GUI and directory
setappdata(handles.SliderFrames,'FIONASetting',FIONASetting);
save('FIONASetting.mat','FIONASetting');

% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3

% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in StartFIONA.
function StartFIONA_Callback(hObject, eventdata, handles)
% hObject    handle to StartFIONA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get FIONASetting
FIONASetting = getappdata(handles.SliderFrames,'FIONASetting');

% Start FIONA analysis on all files
% Output is a matrix of [PhotonNo xPrecision yPrecision xCenter yCenter RSquare ErrorStatus]
% Still need to work on allowing it to do parallel processing
% Get folder for files and assign filenames to FileInputName

% Sum frames to lengthen exposure time
SumFactor = 1;              % Number of frames to be added together

% List of parameters
ActPix = str2double(FIONASetting.ActPix);               % Actual pixel size in nanometer
ObjMag = str2double(FIONASetting.ObjMag);               % Objective magnification
AddMag = str2double(FIONASetting.AddMag);               % Additional magnification

CCDsens = str2double(FIONASetting.CCDSens);             % CCD sensitivity of the camera at specific readout rate and pre-amp setting. See note below for values
EMgain = str2double(FIONASetting.EMgain);               % Electron multiplying (EM) gain setting of camera during acquisition

FileType = '.tif';                       % Input file type
Parallel = FIONASetting.Parallel;           % (No or Yes) for no and with parallel computing using matlabpool

% Get CodePath
CodePath = get(handles.FileLocation,'string');

% Get DataPath
Folder = get(handles.FileName,'string');
% Check if Folder is empty. if it not, check if the directory exists
if isempty(Folder)
    % Check if FolderName exists
    FolderName = getappdata(handles.SliderFrames,'FolderName');
    % If FolderName exist, use the FolderName as DataPath, otherwise ask
    % to get directory
    if isempty(FolderName)
        DataPath = uigetdir;
    else
        DataPath = [CodePath '\' FolderName];
    end
else
    DataPath = [CodePath '\' Folder];
    % If directory does not exist, get directory
    if ~exist(DataPath,'dir')
        DataPath = uigetdir;
    end
end

name=strrep(DataPath,[CodePath '\'],'');
set(handles.FileName,'String',name);
cd(DataPath);

% Get all the file names and assign them under FileInput
FileIn = dir(['*' FileType]);
FileInput = cell(length(FileIn),1);
FileInputName = cell(length(FileIn),1);
for ind=1:length(FileIn)
    FileInput{ind}=FileIn(ind).name;
    FileInputName{ind}=strrep(FileIn(ind).name,FileType,'');
end
    
% Applying FIONA to each file and generate output file
PixelSize = ActPix/(ObjMag*AddMag);     % Pixel size in nanometer
CountToPhoton = CCDsens/EMgain;         % Count to photon conversion
a2 = PixelSize*PixelSize;               % Square of PixelSize

% Initialize TotalFrameAnalyzed and set ProgressStep size
TotalFrameAnalyzed = 0;
ProgressStep = 10;                      % How frequent progress is updated (in frames)

% Allow use of multiple cores if the 'Parallel' option is 'Yes'
if strcmp(Parallel,'Yes')
    % Check if existing matlabpool is opened. Matlabpool allows the use
    % of more than one core of the computer
    IsOpen=matlabpool('size')>0;
    if IsOpen == 0
        matlabpool open;
    end
    % Find out the number of cores in the computer. If the number of cores
    % returned is less than expected, it can be increased by going to
    % Parallel>Manage Configurations, right click on 'local', click 
    % 'Properties', under the 'Scheduler', you can find 'Number of workers
    % available to scheduler (ClusterSize). You can increase this number. You
    % can find out if the computer is making use of all of its cores using task
    % manager in Windows or Activity Monitor in Mac. In my experience, 
    % Matlabpool will only use half of the available cores by default.
    CoreNo = matlabpool('size');
else
    CoreNo = 1;
end

tic;
% Looping over all files
for ind=1:length(FileIn)
    % Getting info from image to make it faster to upload multiple frames
    cd(DataPath);
    info = imfinfo(FileInput{ind});
    num_images = numel(info);
    
    % Allocating memory to variables
    Output = zeros(num_images,7);
    %BothOutputs = zeros(num_images,4);
    
    % Loop over all frames
    for k = 0:num_images-1
        data = double(zeros(info(1).Height,info(1).Width));
        for j = 1:SumFactor
            data = data+double(imread(FileInput{ind}, (SumFactor*k)+j, 'Info', info));
        end
        assignin('base','data',data);

        %data=double(data);
        [ny,nx]=size(data);                 % Find out the dimensions
        grid = [nx ny 1:nx 1:ny];         % Gridding input for gauss2dfunct and gauss2dfit
        tilt=0;

        % Change to code directory
        cd(CodePath);

        % Parameters: p(1): z-offset, p(2): amplitude, p(3): xStdev, p(4): yStdev, p(5): xCenter, p(6): yCenter, p(7): tilt.
        try
            [popt,resnorm,residual,ret]=gauss2dfit(data,grid,tilt);

            if popt(5)>0 && popt(5)<nx && popt(6)>0 && popt(6)<ny

                % Getting center and precision
                xCenter = popt(5);                  % Center of x in pixels
                yCenter = popt(6);                  % Center of y in pixels

                PixelSize2 = PixelSize*PixelSize;
                sx2 = popt(3)*popt(3);              % Square of xStdev
                sx4 = sx2*sx2;                      % xStdev^4
                sy2 = popt(4)*popt(4);              % Square of yStdev
                sy4 = sy2*sy2;                      % yStdev^4

                % To estimate b, which is the standard deviation of the background,
                % we'll look at the z-offset (popt(1)) and calculate the standard
                % deviation based on anything below the z-offset. Before that we
                % would want to make everything 3 standard deviations away from
                % our spot to be (z-offset + 1) or NaN.

                % Find the limits of the data
                xmin = floor(popt(5)-4*popt(3));
                xmax = ceil(popt(5)+4*popt(3));
                ymin = floor(popt(6)-4*popt(4));
                ymax = ceil(popt(6)+4*popt(4));

                if xmin < 1
                    xmin = 1;
                elseif xmax > nx
                    xmax = nx;
                elseif ymin < 1
                    ymin = 1;
                elseif ymax > ny
                    ymax = ny;
                end

                z = reshape(data,ny*nx,1); 
                SumSquaresTotal = sum(z.*z);
                RSquare = 1-(resnorm/SumSquaresTotal);

                residual = reshape(residual,ny,nx);     % Extract residual
                residual = -residual;                   % Invert such that residual = data - zfit
                residual(xmin:xmax,ymin:ymax)=0;        % Set values of residual around center to be NaN
                residual = residual(residual<0);
                b = sqrt(sum(residual.*residual)/(length(residual)-1));
                b = b*CountToPhoton;                    % Standard deviation of the background
                b2 = b*b;                               % Square of background
                PhotonNo = abs(2*pi*popt(2)*popt(3)*popt(4)*CountToPhoton);   % Number of Photons calculated using volume under gaussian, which is 2*pi*A*stdev(x)*stdev(y)
                PhotonNo2 = PhotonNo*PhotonNo;          % Square of PhotonNo
                xPrecision = sqrt((sx2/PhotonNo) + (a2/(12*PhotonNo)) + (8*pi*sx4*b2/(a2*PhotonNo2)));
                yPrecision = sqrt((sy2/PhotonNo) + (a2/(12*PhotonNo)) + (8*pi*sy4*b2/(a2*PhotonNo2)));

                ErrorStatus = 0;

                cd(DataPath);
                Output(k+1,:) = [PhotonNo xPrecision yPrecision xCenter yCenter RSquare ErrorStatus];
            else
                cd(DataPath);
                ErrorStatus = 1;
                Output(k+1,:) = [0 0 0 0 0 0 ErrorStatus];
            end
        catch
            cd(DataPath);
            ErrorStatus = 2;
            Output(k+1,:) = [0 0 0 0 0 0 ErrorStatus];
        end

        % Display progress
        if mod(k,ProgressStep)==0
            TotalFrameAnalyzed = TotalFrameAnalyzed + ProgressStep;
            String1 = [' Frames analyzed:           ' num2str(k) ' / ' num2str(num_images)];
            String2 = [' File analyzed:                  ' num2str(ind) ' / ' num2str(length(FileIn))];
            String3 = [' Elapsed time:                   ' num2str(floor(toc)) ' s'];
            String4 = [' Total frames analyzed:    ' num2str(TotalFrameAnalyzed)];
            set(handles.listbox3,'String',{String1;String2;String3;String4});
            drawnow();
        end
        
        % Make points with errors equal previous values
        if 0
            ErrorRow=find(Output(:,7)>0);   % Row at which there are errors
            if isempty(ErrorRow)==0
                % Take out xOut at the beginning
                loop=1;
                while loop == 1
                    if ErrorRow(1)==1
                        Output(1,:)=[];
                        ErrorRow(1)=[];
                        ErrorRow=ErrorRow-1;
                    else
                        loop=0;
                    end
                end
                % Make points on ErrorRow equal previous values
                for l=1:length(ErrorRow)
                    Output(ErrorRow(l),1:6)=Output((ErrorRow(l)-1),1:6);
                end
            end
        end

        % Save Output
        dlmwrite([FileInputName{ind} '(' num2str(SumFactor) 'x).txt'],Output,'\t')
    end
    
    axes(handles.Preview);
    hold off;
    plot(Output(:,4), Output(:,5),'g');
    drawnow();
end

% Notes on CCD sensitivity for different cameras
% Sub-zero (Readout rate on left (e.g. 10 MHz 14 bit), Preamp setting in
% bracket (e.g. 1x, 2.3x, 4.5x), CCD sensitivity after colon)
%   10 MHz 14 bit: 67.42 (1x), 26.85 (2.3x), 12.27 (4.5x)
%    5 MHz 14 bit: 59.86 (1x), 24.36 (2.3x), 10.82 (4.5x)
%    3 MHz 14 bit: 59.84 (1x), 24.17 (2.3x), 10.72 (4.5x)
%    1 MHz 16 bit: 24.33 (1x),  9.74 (2.3x),  4.26 (4.5x)
% 
% Johnny Cage 
% Readout rate on left (e.g. 10 MHz 14 bit), Preamp setting in
% bracket (e.g. 1x, 2.4x, 4.9x), CCD sensitivity after colon
%   10 MHz 14 bit: 64.5 (1x), 26.3 (2.4x), 11.9 (4.9x)
%    5 MHz 14 bit: 55.0 (1x), 23.9 (2.4x), 10.4 (4.9x)
%    3 MHz 14 bit: 55.8 (1x), 23.6 (2.4x), 10.3 (4.9x)
%    1 MHz 16 bit: 23.3 (1x),  9.1 (2.4x),  4.2 (4.9x)
% 
% Scorpion
% Readout rate on left (e.g. 10 MHz 14 bit), Preamp setting in
% bracket (e.g. 1x, 2.3x, 4.5x), CCD sensitivity after colon
%   10 MHz 14 bit: 67.75 (1x), 27.00 (2.5x), 12.13 (5.2x)
%    5 MHz 14 bit: 57.20 (1x), 23.20 (2.5x), 10.33 (5.2x)
%    3 MHz 14 bit: 56.17 (1x), 23.07 (2.5x), 10.05 (5.2x)
%    1 MHz 16 bit: 23.02 (1x),  9.29 (2.5x),  4.18 (5.2x)

clear i;

% --- Executes on slider movement.
function sliderOutlier_Callback(hObject, eventdata, handles)
% hObject    handle to sliderOutlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Adjust the values of TextOutlier and SliderOutlier
StdThreshold = get(handles.sliderOutlier,'Value');
assignin('base','StdThreshold',StdThreshold);
StdThreshold = floor(StdThreshold*10)/10;
%assignin('base','StdThreshold',StdThreshold);
set(handles.textOutlier,'String',StdThreshold);
set(handles.sliderOutlier,'Value',StdThreshold);

% --- Executes during object creation, after setting all properties.
function sliderOutlier_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderOutlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function textOutlier_Callback(hObject, eventdata, handles)
% hObject    handle to textOutlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of textOutlier as text
%        str2double(get(hObject,'String')) returns contents of textOutlier as a double

StdThreshold = str2num(get(handles.textOutlier,'String'));
StdThreshold = floor(StdThreshold*10)/10;
set(handles.textOutlier,'String',StdThreshold);
set(handles.sliderOutlier,'Value',StdThreshold);

% --- Executes during object creation, after setting all properties.
function textOutlier_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textOutlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in CommandScreenOutlier.
function CommandScreenOutlier_Callback(hObject, eventdata, handles)
% hObject    handle to CommandScreenOutlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set pointer to hourglass
set(gcf, 'pointer', 'watch')
drawnow;

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');
StdThreshold = get(handles.sliderOutlier,'Value');
%assignin('base','dataInput',dataInput);
%axes(handles.Preview);
NoOfPoints = size(dataInput,1);

if NoOfPoints <= 500; span = 0.1;
elseif NoOfPoints <= 5000; span = 0.01;
elseif NoOfPoints <= 50000; span = 0.001;
else span = 0.0001;
end

xSmooth = smooth(dataInput(:,3),dataInput(:,1),span,'rlowess');
ySmooth = smooth(dataInput(:,3),dataInput(:,2),span,'rlowess');
%figure;
%plot(dataInput(:,3),dataInput(:,1),'.b'); hold on; plot(dataInput(:,3),xSmooth,'r');

xResidual = dataInput(:,1)-xSmooth; yResidual = dataInput(:,2)-ySmooth;
[xHistFreq,xHistCenter]=hist(xResidual,100); [yHistFreq,yHistCenter]=hist(yResidual,100);
xMax=max(xHistFreq); yMax=max(yHistFreq);

f = @(p,x)p(1)*exp(-((x-p(2))/p(3)).^2);
px = [xMax 0 std(xResidual)]; py = [yMax 0 std(yResidual)];
pfitx = lsqcurvefit(f,px,xHistCenter',xHistFreq'); pfity = lsqcurvefit(f,py,yHistCenter',yHistFreq');

xOutliers = find(abs(xResidual) > StdThreshold*pfitx(3));
yOutliers = find(abs(yResidual) > StdThreshold*pfity(3));
Outliers = unique([xOutliers; yOutliers]);

dataOutlier = dataInput(Outliers,:);
dataGood = dataInput; dataGood(Outliers,:)=[];

axes(handles.Preview); hold off;
plot(dataInput(:,1),dataInput(:,2),'g'); hold on;
plot(dataOutlier(:,1),dataOutlier(:,2),'or'); hold on; 
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
hold off;

datacursormode on;    

% Store dataInput and dataGood
setappdata(handles.SliderFrames,'dataInput',dataInput);
setappdata(handles.SliderFrames,'dataGood',dataGood);

% Set pointer to arrow
set(gcf, 'pointer', 'arrow')

h=figure;
subplot(5,2,1:2); plot(dataInput(:,3),dataInput(:,1),'g'); hold on;
plot(dataInput(xOutliers,3),dataInput(xOutliers,1),'or', 'MarkerSize', 3);
plot(dataInput(:,3),xSmooth(:,1),'b');
ylabel('x vs t','fontweight','b','fontsize',10,'FontName','Palatino Linotype');
subplot(5,2,3:4); plot(dataInput(:,3),dataInput(:,2),'g'); hold on;
plot(dataInput(yOutliers,3),dataInput(yOutliers,2),'or', 'MarkerSize', 3);
plot(dataInput(:,3),ySmooth(:,1),'b');
ylabel('y vs t','fontweight','b','fontsize',10,'FontName','Palatino Linotype');
subplot(5,2,[5,7,9]); plot(dataInput(:,1),dataInput(:,2),'g'); hold on;
plot(dataOutlier(:,1),dataOutlier(:,2),'or', 'MarkerSize', 3); hold on; 
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
xlabel('Outliers Plot','fontweight','b','fontsize',10,'FontName','Palatino Linotype');
hold off;
subplot(5,2,[6,8,10]); plot(dataGood(:,1),dataGood(:,2),'g'); hold on;
plot(dataGood(1,1),dataGood(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataGood(end,1),dataGood(end,2),'o','MarkerFaceColor','black','MarkerSize',10); 
xlabel('Cleaned Plot','fontweight','b','fontsize',10,'FontName','Palatino Linotype');
hold off;
setappdata(handles.SliderFrames,'h',h);

% --- Executes on button press in CommandDeleteOutlier.
function CommandDeleteOutlier_Callback(hObject, eventdata, handles)
% hObject    handle to CommandDeleteOutlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataInput = getappdata(handles.SliderFrames,'dataGood');
setappdata(handles.SliderFrames,'dataInput',dataInput);
h=getappdata(handles.SliderFrames,'h');
%assignin('base','h',h);
try close(h); catch; end

axes(handles.Preview); hold off;
plot(dataInput(:,1),dataInput(:,2),'g'); hold on;
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10);
hold off;

datacursormode on;  

% --- Executes on button press in CommandResaveTrace.
function CommandResaveTrace_Callback(hObject, eventdata, handles)
% hObject    handle to CommandResaveTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get CodePath, DataPath and FileInput
CodePath = getappdata(handles.SliderFrames,'CodePath');
DataPath = getappdata(handles.SliderFrames,'DataPath');
cd(DataPath);

% Get rotated trace variable
OnAxisTrace = getappdata(handles.SliderFrames,'OnAxisTrace');
assignin('base','OnAxisTrace',OnAxisTrace);    

% Get OpenCleanedFile and OpenTransformedFile
OpenCleanedFile = get(handles.OpenCleaned,'Value');
OpenTransformedFile = get(handles.OpenTransformed,'Value');

% Get FileInput and FileInputName
FileInputName=getappdata(handles.SliderFrames,'FileInputName');

% Get the FileNo from SliderFileNo
FileNo = floor(get(handles.SliderFileNo,'Value'));

% Delete existing figures outside GUI
h=getappdata(handles.SliderFrames,'h');
try close(h); catch; end

% Get cursors
dcm_obj = datacursormode;
Cursor = getCursorInfo(dcm_obj);
assignin('base','Cursor',Cursor);

% Issue warning dialog if the number of cursors is not even and not more than 2
if length(Cursor) < 2 || mod(length(Cursor),2)==1
    h=warndlg('You need to have at least 2 cursors and even number of cursors');
else
    % Sort Cursors
    CursorPositions = zeros(length(Cursor),3);
    for i = 1:length(Cursor)
        CursorPositions(i,1:2) = Cursor(i).Position;
        CursorPositions(i,3) = Cursor(i).DataIndex;
    end
    CursorPositions = sortrows(CursorPositions,3);
    
    % Select new data and draw figures with selected parts
    DataNew = cell(length(Cursor)/2,1);
    axes(handles.Preview); hold off;
    for i = 1:length(Cursor)/2
        DataNew{i} = OnAxisTrace(CursorPositions(2*i-1,3):CursorPositions(2*i,3),:);
        plot(DataNew{i}(:,2), DataNew{i}(:,1)); hold on;
    end

    % Find files
    FileType = '.txt';
    FileIn = dir(['*' FileType]);
    FileInput = cell(length(FileIn),1);
    FileInputName = cell(length(FileIn),1);
    for ind=1:length(FileIn)
        FileInput{ind}=FileIn(ind).name;
        FileInputName{ind}=strrep(FileIn(ind).name,FileType,'');
    end
    FileInputString = cell(length(FileInput),1);

    if OpenCleanedFile == 1
        DataPath = strrep(DataPath,'\Cleaned','');
        FileName = strrep(FileInputName{FileNo},'-Cleaned','');
        cd([DataPath '\T-test']);
        
        % Delete previous data
        fclose all;
        delete([FileName '*.txt']);

        % Save data for t-test input
        for i = 1:length(Cursor)/2
            dlmwrite([FileName 'TtestLin-' num2str(i) '.txt'],DataNew{i},'\t');
        end
    elseif OpenTransformedFile == 1
        DataPath = strrep(DataPath,'\Transformed','');
        FileName = strrep(FileInputName{FileNo},'-Transformed','');
        cd([DataPath '\T-test']);
        
        % Delete previous data
        fclose all;
        delete([FileName '*.txt']);

        % Save data for t-test input
        for i = 1:length(Cursor)/2
            dlmwrite([FileName 'TtestLin-' num2str(i) '.txt'],DataNew{i},'\t');
        end
    else
        cd([DataPath '\T-test']);
        
        % Delete previous data
        fclose all;
        delete([FileInputName{FileNo} '*.txt']);

        % Save data for t-test input
        for i = 1:length(Cursor)/2
            dlmwrite([FileInputName{FileNo} 'TtestSegment-' num2str(i) '.txt'],DataNew{i},'\t');
        end
        
        % Update listbox
        num_images = length(FileInput);
        for i = 1:num_images
            if exist(['Cleaned\' FileInputName{i} '-Cleaned.txt'],'file') == 2
                FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"     (saved)'];
            else
                FileInputString{i} = [ '(' num2str(i) ')   "' FileInput{i} '"'];
            end
            if exist(['Transformed\' FileInputName{i} '-Transformed.txt'],'file') == 2
                FileInputString{i} = [FileInputString{i} '     (transformed)'];
            end     
        end
        String1 = ['No of files : ' num2str(num_images)];
        set(handles.listbox4,'String',[String1; ' '; FileInputString]);
    end

end

datacursormode on;

cd(CodePath);
    
% --- Executes on button press in CommandDrawBorder.
function CommandDrawBorder_Callback(hObject, eventdata, handles)
% hObject    handle to CommandDrawBorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get dataInput
dataInput = getappdata(handles.SliderFrames,'dataInput');

% Plot trace
axes(handles.Preview);
plot(dataInput(:,1),dataInput(:,2),'g'); hold on;
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10);
drawnow();

% Get Polyline
[xBorder,yBorder] = getline(gca);

% Draw Polyline
axes(handles.Preview);
plot([xBorder; xBorder(1)],[yBorder; yBorder(1)],'r'); 

% Find points within polyline to delete
BoundPoints = inpolygon(dataInput(:,1),dataInput(:,2),xBorder,yBorder);
axes(handles.Preview);
plot(dataInput(BoundPoints,1),dataInput(BoundPoints,2),'or'); hold off;

setappdata(handles.SliderFrames,'BoundPoints', BoundPoints);

datacursormode on;

% --- Executes on button press in CommandDeleteBorder.
function CommandDeleteBorder_Callback(hObject, eventdata, handles)
% hObject    handle to CommandDeleteBorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataInput = getappdata(handles.SliderFrames,'dataInput');
%assignin('base','dataInput',dataInput);
BoundPoints = getappdata(handles.SliderFrames,'BoundPoints');
%assignin('base','BoundPoints',BoundPoints);
dataInput(BoundPoints,:) = [];

axes(handles.Preview); hold off;
plot(dataInput(:,1),dataInput(:,2),'g'); hold on;
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10);

setappdata(handles.SliderFrames,'dataInput',dataInput);

datacursormode on;  


% --- Executes on button press in CommandDeleteBorderOut.
function CommandDeleteBorderOut_Callback(hObject, eventdata, handles)
% hObject    handle to CommandDeleteBorderOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataInput = getappdata(handles.SliderFrames,'dataInput');
%assignin('base','dataInput',dataInput);
BoundPoints = getappdata(handles.SliderFrames,'BoundPoints');
%assignin('base','BoundPoints',BoundPoints);
dataInput = dataInput(BoundPoints,:);

axes(handles.Preview); hold off;
plot(dataInput(:,1),dataInput(:,2),'g'); hold on;
plot(dataInput(1,1),dataInput(1,2),'p','MarkerFaceColor','black','MarkerSize',15); 
plot(dataInput(end,1),dataInput(end,2),'o','MarkerFaceColor','black','MarkerSize',10);

setappdata(handles.SliderFrames,'dataInput',dataInput);

datacursormode on; 
