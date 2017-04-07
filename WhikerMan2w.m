function varargout = WhikerMan2w(varargin)
% WHIKERMAN2L M-file for WhikerMan2l.fig
%      WHIKERMAN2L, by itself, creates a new WHIKERMAN2L or raises the
%      existing
%      singleton*.
%
%      H = WHIKERMAN2L returns the handle to a new WHIKERMAN2L or the
%      handle to
%      the existing singleton*.
%4
%      WHIKERMAN2L('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WHIKERMAN2L.M with the given input arguments.
%
%      WHIKERMAN2L('Property','Value',...) creates a new WHIKERMAN2L or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before WhikerMan2l_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WhikerMan2l_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WhikerMan2l

% Last Modified by GUIDE v2.5 20-Feb-2017 14:13:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @WhikerMan2l_OpeningFcn, ...
    'gui_OutputFcn',  @WhikerMan2l_OutputFcn, ...
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


% --- Executes just before WhikerMan2l is made visible.
function WhikerMan2l_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WhikerMan2l (see VARARGIN)

% Choose default command line output for WhikerMan2l
handles.output = hObject;

clc

% default pmtrs
handles.red = [1 0 0]; handles.green = [0 1 0];

handles.sigma_prior = .1;    % ditto
set(handles.edit_sigma_prior,'String',handles.sigma_prior);

handles.default_energy_threshold_subtract_image_mean = -20;
handles.default_energy_threshold_do_not_subtract_image_mean = 75;

handles.frame_interval = 1;
set(handles.edit_frame_interval,'String',handles.frame_interval)

handles.continuous_tracking = false;
set(handles.radiobutton_continuous_tracking,'value',handles.continuous_tracking)

handles.subtract_image_mean = false;
set(handles.radiobutton_subtract_image_mean,'value',handles.subtract_image_mean)

if get(handles.radiobutton_subtract_image_mean,'Value')
    handles.energy_threshold = handles.default_energy_threshold_subtract_image_mean;
else
    handles.energy_threshold = handles.default_energy_threshold_do_not_subtract_image_mean;
end
set(handles.edit_energy_threshold,'String',handles.energy_threshold);

handles.energy_threshold_tolerance = 5;% percent
set(handles.edit_energy_threshold_tolerance,'String',handles.energy_threshold_tolerance);

handles.adaptive_energy_threshold = false;
set(handles.radiobutton_adaptive_energy_threshold,'value',handles.adaptive_energy_threshold);

handles.miumax = 0.2;
set(handles.edit_miumax,'String',handles.miumax);

handles.snout_sigma = 8;
set(handles.edit_snout_sigma,'String',handles.snout_sigma);





% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WhikerMan2l wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WhikerMan2l_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in popupmenu_choose_file.
function popupmenu_choose_file_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_choose_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_choose_file

% % Reset parameters to default values

val = get(hObject,'Value');
if(val~=1)
    string_list = get(hObject,'String');
    handles.fname = string_list{val};   % avi file or dat file.
else
    error('')
end
clear val string_list

handles = initialise_new_video(handles);
handles.trfname = '';

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = initialise_new_video(handles)

% Fixed pmtrs (and therefore don't need to be loaded from tr file)
handles.dt = 0.01;
handles.current_pt = 0;
hsize = [3 3]; % defaul
sigma = .5;
handles.gaussian = fspecial('gaussian', hsize, sigma);
clear hsize sigma

% handles.fname(end-2:end),pause


switch handles.fname(end-2:end)
    case 'avi'
        video.type = 'avi';
        video.vObj = mmreader(handles.fname);
        handles.nframes = video.vObj.NumberOfFrames;
        video.width = video.vObj.Width;
        video.height = video.vObj.Height;
        % specify first frame and last frame:
        handles.startframe = 1;
        handles.stopframe = handles.nframes;
        handles.frameidx = handles.startframe;  % current frame
    case 'dat'
        video.type = 'dat';
        video.fid = fopen(handles.fname,'r');
        video.header = read_mikrotron_datfile_header(video.fid);
        handles.nframes = video.header.nframes;
        video.width = video.header.width;
        video.height = video.header.height;
        video.offset = 8192;
        % set file position to start of first frame
        fseek(video.fid,8192,-1);
        % specify first frame
        handles.startframe = video.header.startframe+1;
        if handles.startframe > 1
            handles.stopframe = handles.startframe - 1;
        else
            handles.stopframe = handles.nframes;          
        end
        handles.frameidx = handles.startframe;  % current frame
        
    otherwise
        error(fprintf('Unhandled video file type: %s\n',handles.fname(end-2:end)))
end
handles.video = video;

axes(handles.video_frame)
p = get(gca,'position');
set(gca,'position',[p(1:2) video.width video.height])
clear p

% load the first frame:
handles.frame = load_frame(video,handles.frameidx);
% %%%%michaela
% handles.frame1(:,:,1) = rot90(handles.frame(:,:,1));
% handles.frame1(:,:,2) = rot90(handles.frame(:,:,2));
% handles.frame1(:,:,3) = rot90(handles.frame(:,:,3));
% handles.frame=[];
% handles.frame1(120:end,:,1) = 0;
% handles.frame1(120:end,:,2) = 0;
% handles.frame1(120:end,:,3) = 0;
% handles.frame=handles.frame1;
% 
% %%%%%

image(handles.frame)
clear frameidx dt t r b

set(handles.frameidx_display,'string',num2str(handles.frameidx))

fn = [handles.fname(1:end-3) 'meanframe'];
if handles.subtract_image_mean
    if exist(fn,'file')
        load(fn,'meanframe','-mat')
    else
        %     get(handles.text_status,'BackgroundColor')
        set(handles.text_status,'BackgroundColor',handles.red)
        set(handles.text_status,'String','BUSY')
        %     get(handles.text_status,'BackgroundColor')
                  
        pause(.01)
        meanframe = zeros(video.height,video.width,3);
        space = 10;
        for fidx = 2:space:handles.nframes
            meanframe = meanframe + double(load_frame(handles, fidx))/(handles.nframes/space);
        end
        save(fn,'meanframe')
        clear space fidx
        set(handles.text_status,'backgroundcolor',handles.green)
        set(handles.text_status,'String','idle')
    end
else
    
    meanframe = zeros(video.height,video.width,3);
end

axes(handles.aux_frame)
imagesc(double(handles.frame(:,:,1))-meanframe(:,:,1))
colormap gray
hold on

handles.meanframe = meanframe;
clear meanframe

% clear tracking variables in workspace:
handles.n_control_pts = 3;
handles.rall = zeros(handles.nframes,2,handles.n_control_pts);
handles.kappa_all = zeros(1,handles.nframes);
handles.theta_all = zeros(1,handles.nframes);
handles.fp_all = zeros(2,handles.nframes);
handles.tracked = zeros(1,handles.nframes);
handles.Emin_all = zeros(1,handles.nframes);
handles.fpidx = zeros(1,handles.nframes);
handles.folliclemask = cell(1,handles.nframes);

handles.current_pt = 0;
handles.s0 = [];

% determine default tracking file
% handles.trfname = '';
set(handles.edit_specify_tr_file,'String',[handles.fname(1:end-4) '.tr']);
clear fn



% guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_choose_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

ff = [dir('*.avi');dir('*.dat')];
clear ffavi ffdat

nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose File';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(hObject,'String',string_list);

%% Pushbutton callbacks

% --- Executes on button press in pushbutton_fit_snakes.
function pushbutton_fit_snakes_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fit_snakes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles, fnoutput] = fit_snakes(handles);
clear fnoutput
guidata(hObject, handles);

function [handles, fnoutput] = fit_snakes(handles)


if ~isfield(handles,'trfname') || isempty(handles.trfname)
    figure, beep
    text(.1,.5,'Specify tr file before starting to track.')
    track = false;
    fnoutput = -1;    
else    
    % initialise
    %! not sure if need the following:
    
    if ~isfield(handles,'tfp')
        handles.tfp = zeros(1,handles.nframes);
    end
    
    % new v2r:
    lastframe = modsubtract(handles.frameidx,1,handles.nframes);
    % lastbutoneframe = modsubtract(handles.frameidx,2,handles.nframes);
    
    if (handles.frameidx==handles.startframe) || ~handles.tracked(lastframe)
        first_frame = true;
    else
        first_frame = false;
    end
    
    track = true;
    nattempts = 1;
    usepredictor = true;
    fnoutput = -1;    
end

% pmtrs = handles.pmtrs;  ! nb to get rid of need for this line

while (track)
    
    lastframe = modsubtract(handles.frameidx,1,handles.nframes);
    lastbutoneframe = modsubtract(handles.frameidx,2,handles.nframes);

    if first_frame || (rem(handles.frameidx,handles.frame_interval)==0)
        doplot = true;
        axes(handles.aux_frame)
        cla reset
    else
        doplot = false;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preprocess:
    im.raw = double(handles.frame(:,:,1)) - handles.meanframe(:,:,1);
    % Smooth image a bit to facilitate whisker tracking
    im.s = imfilter(im.raw,handles.gaussian);
    [im.dx,im.dy] = gradient(im.s);
    
    if doplot
        fname = handles.fname(1:end-4);
        idx = regexp(fname,'\_');
        fname(idx) = '-';
        clear idx
        imagesc(im.s), colormap gray
        title(sprintf('%s: Frame %d',fname,handles.frameidx))
        hold on
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Whisker tracking
    % Set initial conditions for contour parameters
    
    D = handles.n_control_pts;
    rall = handles.rall;
    
    if first_frame==1
        % no tracking solution for previous frame:
        % must set control pts. can be done manually or automatically,
        % depending on setting of the auto init radio button.
        rold = zeros(2,D);
        goodsolution = false;
        if get(handles.radiobutton_auto_initialise,'value')
            % set control pts automatically, using the tr data of the file
            % currently selected in the auto init popupmenu:
            menu = cellstr(get(handles.popupmenu_choose_auto_initialise_file,'String'));
            fn = menu{get(handles.popupmenu_choose_auto_initialise_file,'Value')};
            [rold, goodsolution] = auto_initialise(fn,im.s, handles);
            clear menu fn
        end
        if ~goodsolution
            % if either the autoinit failed or user has chosen not to use
            % it, manually set the control pts.
            % it's necessary that the first point is the one nearest the follicle and that the points thereafter move outwards
            % monotonically along the whisker.
            for i = 1:D
                [rold(1,i),rold(2,i)] = ginput(1);
                plot(rold(1,i),rold(2,i),'wo')
            end
            clear i
        end
    else
        % seed the tracker using solution for previous frame
        rold = reshape(rall(lastframe,:,:),2,D);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the initial conditions using information from previous frames
    % (Quadratic) Bezier curve terminology:
    %   P0 and P2 are the endpoints
    %   P1 is the shape-determining point
    % For P0 and P1, predict location by estimating velocity:
    ddr = zeros(2,D);
    
    % new code in v2r:
    z = modsubtract(handles.frameidx,handles.startframe,handles.nframes);
    if z==handles.nframes
        z = 0;
    end
    if z >= 2
        if handles.tracked(lastframe) && handles.tracked(lastbutoneframe)
            ddr = reshape(rall(lastframe,:,:),2,D) - reshape(rall(lastbutoneframe,:,:),2,D);
        end
    end
    
    % Sometimes, eg due to a sudden deceleration, the velocity-based
    % initial conditions are not helpful.  The 'usepredictor' flag is used
    % to track using simply the solution from the previous frame as the
    % initial condition:
    if usepredictor
        r = predict_contour(rold,ddr,handles.miumax);
    else
        % ie prediction is that contour is same as previous frame
        r = rold;
    end
    clear ddr miumax
    
    if doplot
        % plot solution from previous frame, if there is one:
        plot(rold(1,:),rold(2,:),'b+'),
        % plot initial conditions:
        plot(r(1,:),r(2,:),'c+')
    end
    
    clear rold
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimise Bezier curve to best fit image and prior constraints
    % Constrain P0 and P2 to move normal to contour.  So, first compute unit
    % normal vectors at B(t=0) = P0 and B(t=1) = P2:
    tvec = bezierdtval(r,[0 1]);    % tangent vectors to bezier at s=[0,1] (ie at P0,P2)
    n = sqrt(sum(tvec.^2,1));
    tvec = tvec ./ (ones(2,1)*n);   % unit tangent vectors
    clear n
    R = [0 -1;1 0]; % 90deg rotation matrix
    nvec = R*tvec;       % normal vectors at s=[0,1]
    n = sqrt(sum(nvec.^2,1));
    nvec = nvec ./ (ones(2,1)*n);    % unit normals
    clear tvec R n
    if doplot
        quiver(r(1,[1 3]),r(2,[1 3]),10*nvec(1,:),10*nvec(2,:),0,'c')
    end
    
    % Optimisation parameters:
    % z(1) - component of normal to curve at P0=r(:,1)
    % z(2:3) - P1=r(:,2)
    % z(4) - component of normal to ourve at P2=r(:,3)
    
    % Set initial values
    z0 = [0 r(:,2)' 0];
    % Find 'z' that minimises the cost Bimfun(z)
    [z,Emin,exitflag,output] = fminunc(@(z) Bimfun(z, r, nvec, im, handles.dt, handles.sigma_prior),z0, ....
        optimset('GradObj','On','TolX',.05,'Display','off'));%,'TolFun',.01));
    %     fprintf('Iterations %d\n',output.iterations)
    
    % Bezier curve that locally minimises the cost:
    rnew = zeros(size(r));
    rnew(:,1) = r(:,1) + z(1)*nvec(:,1);
    rnew(:,2) = z(2:3)';
    rnew(:,3) = r(:,3) + z(4)*nvec(:,2);
    clear r z0 z exitflag output
    %     axes(handles.aux_frame)
    clear nvec
    
    b = bezierval(rnew,0:.01:1);
    if doplot
        plot(b(1,:),b(2,:),'r',rnew(1,:),rnew(2,:),'r.')
        title(sprintf('%s: %d (%.1f)', fname, handles.frameidx, Emin))
        %         plot(polyval(p(1,:),handles.sstar),polyval(p(2,:),handles.sstar),'go')
    end
    
    % Snout outline detection:
    % if it's the first_frame, process within large Region Of Interest
    % (ROI).  else, pick an ROI centred around follicle position of
    % previous frame:
    if handles.frameidx == handles.startframe
        if ~isempty(handles.folliclemask{lastframe})
            roi.x = round(handles.fp_all(1,lastframe)) + [-40:40];
            roi.y = round(handles.fp_all(2,lastframe)) + [-40:40];
            %roi.y = round(handles.fp_all(2,lastframe)) + [-40:40];
            
            folliclemask_old = handles.folliclemask{lastframe};
            fpidx_old = handles.fpidx(lastframe);
        else
            if isfield(handles,'snout_roi')
                roi = handles.snout_roi;
            else
                roi.x = [50:size(handles.frame,2)];%50
                roi.y = [50:size(handles.frame,1)];%50
            end
            folliclemask_old = [];
            fpidx_old = [];
        end
    else
        if isfield(handles,'snout_roi')
            roi = handles.snout_roi;
        else
            roi.x = [50:size(handles.frame,2)];
            roi.y = [50:size(handles.frame,1)];%50
        end
        folliclemask_old = [];
        fpidx_old = [];
    end
    
    roi.x = roi.x(roi.x<=size(handles.frame,2));
    roi.y = roi.y(roi.y<=size(handles.frame,1));
    
    %     [p,~] = fit_curve(b,handles.polyorder);
    theta = base_angle(rnew,0);
    [folliclemask] = protract(double(handles.frame(:,:,1)),roi,folliclemask_old, fpidx_old, theta, handles.snout_sigma);
    clear folliclemask_old fpidx_old roi theta
    handles.folliclemask{handles.frameidx} = folliclemask;
    if doplot
        plot(folliclemask(1,:),folliclemask(2,:),'y-')
    end
    
    % Post-process:
    % 1. Find follicle position (intersection of polyfit to bezier with
    % folliclemask).  Do this by fitting a standard 'polyfit' polynomial to
    % the bezier points (since this extrapolates better).
    % 2. Renormalise: shift control points such that arc length from fp to
    % P0 and to P1 equals that of the initial frame.
    if ~isfield(handles,'s0')
        handles.s0 = [];
    end
    [rnew, fp, fpidx, tfp, sfp, s, tstar, sp2] = postproc(rnew, folliclemask, handles.s0, handles.dt, handles.sstar);
    if isempty(handles.s0)
        handles.s0(1) = -sfp;
        handles.s0(2) = sp2;
    end
    
    if doplot
        b = zeros(2,length(s)+1);
        %         b(1,:) = polyval(p(1,:),[sfp s]);
        %         b(2,:) = polyval(p(2,:),[sfp s]);
        t = [0:handles.dt:1];
        b = bezierval(rnew,[tfp t]);
        %         pause
        plot(b(1,:),b(2,:),'r',fp(1),fp(2),'r.')
        title(sprintf('%s: %d (%.1f)', fname, handles.frameidx, Emin))
        plot(rnew(1,:),rnew(2,:),'msq'),
        %         plot(polyval(p(1,:),handles.sstar),polyval(p(2,:),handles.sstar),'go')
        clear b t
    end
    
    
    rall(handles.frameidx,:,:) = rnew;
    
    % Test if solution is good  (energy is low enough).
    %   If yes, go on to next frame
    %   If no, retrack this frame using predictor=false initial conditions.
    %   If this is the second time through, accept the least bad of the two
    %   solutions but stop tracking and return control to the user.
    
    % store current solution, in case need to retrack:
    attempts.Emin(nattempts) = Emin;
    attempts.fp(:,nattempts) = fp;
    attempts.fpidx(:,nattempts) = fpidx;
    attempts.rnew(nattempts,:,:) = rnew;
    %     attempts.p(nattempts,:,:) = p;
    clear Emin fp fpidx p
    
    track = get(handles.radiobutton_continuous_tracking,'Value');
    
    if handles.adaptive_energy_threshold
        % adaptive energy threshold (new 200217 v2w)
        theta = base_angle(squeeze(rnew),tstar);
        clear rnew
        idx = find(abs(handles.theta_all-theta)<3);
        if ~isempty(idx)
%             handles.energy_threshold_tolerance
%             handles.Emin_all(idx)
%             pause
            energy_threshold = (1+handles.energy_threshold_tolerance/100)*median(handles.Emin_all(idx));
            set(handles.text_energy_threshold_adaptive,'String',round(energy_threshold));
        else
            % no relevant data
            energy_threshold = handles.energy_threshold;
        end
        clear theta idx
    else
        energy_threshold = handles.energy_threshold;
    end
    
    if (attempts.Emin(nattempts) > energy_threshold) && ~first_frame
        % contour solution may be bad
        beep, disp('here')
        if nattempts==1
            % try retracking with different initial conditions
            track = 1;
            usepredictor = false;   % this helps with sudden deccelerations
            gotonextframe = false;
            title(sprintf('%s: %d: %.1f BAD SOLUTION: Attempting rescue...',fname, handles.frameidx,attempts.Emin(1)))
        elseif nattempts==2
            % retracking with different initial conditions didn't help.
            % So, stop tracking and return control to user.
            gotonextframe = true;
            track = 0;
            title(sprintf('%s: %d: %.1f %.1f',fname, handles.frameidx,attempts.Emin(1:2)))
            fnoutput = 2;
        else
            error('')
        end
        nattempts = nattempts+1;
    else
        % contour solution is good
        % step on to next frame
        if nattempts==2
            title(sprintf('%s: %d: %.1f %.1f Rescue!',fname, handles.frameidx,attempts.Emin(1:2)))
        end
        gotonextframe = true;
    end
    
    if gotonextframe
        % Take the best (or least-bad) solution amongst those attempted:
        % old (replaced 230114)
        %         [~,idx] = max(attempts.Emin);
        % new (230114):
        [~,idx] = min(attempts.Emin);
        r = attempts.rnew(idx,:,:);
        rall(handles.frameidx,:,:) = r;
        %         p = squeeze(attempts.p(idx,:,:));
        fp = attempts.fp(:,idx);
        handles.Emin_all(handles.frameidx) = attempts.Emin(idx);
        handles.energy_threshold_adaptive(handles.frameidx) = energy_threshold;
        handles.fpidx(handles.frameidx) = attempts.fpidx(idx);
        clear attempts idx
        % Use this solution to extract kinematic parameters:
        %         kappa = curvature(p,handles.sstar);
        %         theta = base_angle(p,handles.sstar);
        kappa = curvature(squeeze(r),tstar);
        theta = base_angle(squeeze(r),tstar);
        clear r
        % Save data in handles:
        handles.rall = rall;
        handles.fp_all(:,handles.frameidx) = fp;
        handles.kappa_all(handles.frameidx) = kappa;
        handles.theta_all(handles.frameidx) = theta;
        handles.tracked(handles.frameidx) = true;
        %         handles.pmtrs = pmtrs;  % 's0' can change
        % Save to disk:
        fp_all = handles.fp_all;
        kappa_all = handles.kappa_all;
        theta_all = handles.theta_all;
        Emin_all = handles.Emin_all;
        pmtrs.sstar = handles.sstar;
        pmtrs.s0 = handles.s0;
        pmtrs.sigma_prior = handles.sigma_prior;
        pmtrs.subtract_image_mean = handles.subtract_image_mean;
        pmtrs.energy_threshold = handles.energy_threshold;
        pmtrs.energy_threshold_tolerance = handles.energy_threshold_tolerance;
        pmtrs.miumax = handles.miumax;
        pmtrs.snout_roi = handles.snout_roi;
        pmtrs.snout_sigma = handles.snout_sigma;
        pmtrs.adaptive_energy_threshold = handles.adaptive_energy_threshold;
%         handles.pmtrs = pmtrs;  % this is needed for the edit_whiskertag function, in case case user changes whiskertag and trf needs to be resaved.
        save(handles.trfname,'rall','kappa_all','theta_all','fp_all','pmtrs','Emin_all')
        clear pmtrs 
        if get(handles.radiobutton_kinematics_test,'Value')
            % test curvature
            % plot circle with given curvature
            rad = abs(1/kappa);
            ttest = 0.5;
            cc = bezierval(squeeze(rall(handles.frameidx,:,:)),ttest);
            tmp = bezierval(squeeze(rall(handles.frameidx,:,:)),ttest+[-.01 .01]);
            tv = diff(tmp,1,2); tv = tv/norm(tv);
            nv = [0 1;-1 0] * tv;
            cc = cc + sign(-kappa) * rad * nv;
            tall = 0:.01:2*pi;
            c = zeros(2,length(tall));
            for tidx = 1:length(tall)
                t = tall(tidx);
                %                 rad*[cos(t);sin(t)],cc
                c(:,tidx) = cc + rad*[cos(t);sin(t)];
            end
            plot(c(1,:),c(2,:),'g:')
            clear rad cc tall c tidx t ttest
        end
        
        clear fp kappa theta
        clear rall kappa_all theta_all fp_all Emin_all
        % Unless it's last frame, proceed:
        if handles.frameidx ~= handles.stopframe;
            handles.frameidx = modadd(handles.frameidx,1,handles.nframes);
            handles.frame = load_frame(handles.video,handles.frameidx);
            set(handles.frameidx_display,'string',num2str(handles.frameidx))
            if doplot
                cla(handles.video_frame,'reset')
                 image(handles.frame,'Parent',handles.video_frame)
            end
            %     colormap gray
        else
            % end of video
            track = 0;
            fnoutput = 1;
            beep,pause(.5),beep,pause(.5),beep,pause(.5),beep,pause(.5),beep
        end
        nattempts = 1;  % new frame, first attempt
        usepredictor = true;    % default
        first_frame = 0;
        
    end
    
    
    pause(.01)
    
end % while loop
clear track


%% windowbutton callbacks

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vpos = get(handles.video_frame,'position'); % origin is upper left
pt = get(hObject,'CurrentPoint');    % mouse click location wrt figure
% set(handles.text_mouse,'String',num2str(pt))

% axes(handles.video_frame)
% hold on

% location wrt video_frame axes:
tpt(1) = pt(1)-vpos(1);
tpt(2) = vpos(4)-(pt(2)-vpos(2));
% plot(tpt(1),tpt(2),'bx')
% set(handles.text_axes,'String',num2str(tpt))
handles.CurrentPoint = tpt;




%% subfunctions

function [rinit, goodsolution] = auto_initialise(trfile, image, handles)

load(trfile,'-mat','theta_all','rall')

% choose a representative subset of frames:
[theta,frame] = sort(theta_all);
theta = theta(2:end-1);
frame = frame(2:end-1);
n = length(theta);
% Crude way to pick the subset, since biased to the most common angle.
% But, so long as the subset is big enough, it should satisfice.
interval = 2;%max([floor(n/100),1]);
frame_subset = frame(1:interval:end);

% for each selected frame, compute line integral (E) of contour on
% current image (high E is bad):
dt = 0.01;
t = 0:dt:1;
E = zeros(1,length(frame_subset));
figure,
% set(gca,'ydir','reverse')
imagesc(image), colormap gray, hold on
for fidx = 1:length(frame_subset)
    frame = frame_subset(fidx);
    bez = squeeze(rall(frame,:,:));
    r = bezierval(bez,t);
    plot(r(1,:),r(2,:),'r')
    E(fidx) = dt*sum(interp2(image,r(1,:),r(2,:)));
end
clear fidx frame bez r

% find the contour/frame with the lowest E:
[Emin,Eidx] = min(E);
xlabel(sprintf('min(E)=%f',Emin))
if Emin < handles.energy_threshold
    % we've found a decent initial condition.  Use it:
    rinit = squeeze(rall(frame_subset(Eidx),:,:));
    goodsolution = true;
else
    % we've not.  Flag this and force manual approach.
    rinit = zeros(size(rall,2),size(rall,3));
    goodsolution = false;
    title('Press any key to continue with manual initialisation')
    pause
end
clear E Eidx
close



function frame = load_frame(videopmtrs,frameidx)

switch videopmtrs.type
    case 'avi'
        video = read(videopmtrs.vObj, frameidx);
        frame = video(:,:,:,1);
    case 'dat'
        offset = videopmtrs.header.imagesize * (frameidx-1) + videopmtrs.offset;
        fseek(videopmtrs.fid,offset,-1);
        tmp = fread(videopmtrs.fid,videopmtrs.header.imagesize-24,'uint8=>uint8');
        tmp = reshape([tmp; zeros(24,1)],videopmtrs.width,videopmtrs.height)';
        frame = uint8(zeros(videopmtrs.height,videopmtrs.width,3));
        frame(:,:,1) = tmp;
        frame(:,:,2) = tmp;
        frame(:,:,3) = tmp;
        clear tmp
    otherwise
        error('Unhandled video file type')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = predict_contour(rold, ddr, miumax)
% The predictor is a linear combination of two terms.
% Term 1: P0 & P2 come from r = rold + ddr
%         'P1prev' is set so that the contour has the same shape
%         (curvature) as in the previous frame
% Term 2: P0 & P2 as above.
%         'P1linear' is the midpoint of P0 and P2.  The point of this is
%         that, when the contour is linear, the P1 coordinate is
%         ill-defined (any point along the P0-P2 line works).
%         Using the midpoint gives more stable extrapolation
%         beyond the range t=[0,1] and avoids instability.
% P1prev and P1linear are combined according to how linear the contour
% was in the previous frame.  If it was near-linear, P1linear is
% weighted most; if it was curved, P1prev is weighted most.
% Computation of term 1:
P0old = rold(:,1);
P1old = rold(:,2);
P2old = rold(:,3);
% P0 and P2 are easy:
P0 = P0old + ddr(:,1);
P2 = P2old + ddr(:,3);
% To set P1prev, the idea is that the shape of the quadratic Bezier
% curve is defined by the triangle P0-P1-P2.  So, once we've moved
% P0 and P2, we need to move P1 such that the triangle
% P0-P1-P2 is the same size and shape (but not orientation or
% position) as the triangle P0old-P1old-P2old.
% So, compute angle of line P0-P1 wrt the line P0-P2 so that we can
% find P1prev by rotation wrt the new P0-P2 line:
tmp = ((P2old-P0old)'*(P1old-P0old))/(norm(P2old-P0old)*norm(P1old-P0old));
% dirty hack 300114 copied 04/03/2014
tmp = min([tmp 1]);
theta = acos(tmp);

n = [0 -1; 1 0]*(P2old-P0old);  %normal to P0-P2
if (P1old-P0old)'*n > 0
    % sign of angle is correct
else
    theta = -theta;
end
clear tmp n
% Set P1 by rotation wrt the (new) P0-P2 and rescaling wrt the
% length of the line P0old-P1old:
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
P1prev = P0 + norm(P1old-P0old)*R*(P2-P0)/norm(P2-P0);
clear theta R
% Computation of term 2:
P1linear = (P0+P2)/2;
% Linearly combine P1prev and P1linear.  To determine the balance
% between them, quantify how curved the curve of the previous frame
% was.
% 'miu' is length of the normal from P1 to P0-P2, divided by the
% length of P0-P2.  mui = 0 for a straight line.  miu>0 for a
% curve:
miu = norm(P0old-2*P1old+P2old)/norm(P2old-P0old);
% Introduce a sensitivity parameter.  When miu>=miumax, P1=P1prev.
P1 = min([miu/miumax,1])*P1prev + (1-min([miu/miumax,1]))*P1linear;
clear P1prev P1linear miu miumax
r = [P0 P1 P2];
clear P0old P1old P2old P0 P1 P2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnew, fp, fpidx, tfp, sfp, s, tstar, sp2] = postproc(r, folliclemask, s0, dt, sstar)

% compute arc length parameter 's':
t = -1:dt:1.5;
idx = find(t==0);
b = bezierval(r,t);
db = diff(b,1,2);
ds = sqrt(sum(db.^2,1));
s = [0 cumsum(ds)]; % s(t=-1) = 0
s = s - s(idx); % now s(t=0) = 0
idx = find(t==1);
sp2 = s(idx);
clear db ds idx

% Find intersection between whisker polynomial and folliclemask
[fp, fpidx, tfp] = find_follicle(r,dt,folliclemask);
% if B(s) is the (x,y) point at 's', then fp = B(s(fpidx)).  this is used in snout detection.
[~,idx] = min((t-tfp).^2);
if length(idx)~=1
    error('')
end
sfp = s(idx);
clear idx

% Renormalise contour using polynomial (non-Bezier) approach:
% Relocate P0 and P2 so that there are at standardised distances from fp
% (s0) and adjust P1 to preserve shape.
if ~isempty(s0)
    % s0(1) is location of P0 wrt fp on the reference contour (first_frame)
    rnew = zeros(size(r));
    %             rnew(:,1) = [polyval(p(1,:),sfp+s0(1));polyval(p(2,:),sfp+s0(1))];
    %             rnew(:,3) = [polyval(p(1,:),sfp+s0(2));polyval(p(2,:),sfp+s0(2))];
    %     idx = find(t==tfp);
    sprime = s - sfp;   % now sprime(t=tfp) = 0
    d = (sprime-s0(1)).^2;
    [~,idx] = min(d);
    rnew(:,1) = b(:,idx);
    clear d idx
    d = (sprime-s0(2)).^2;
    [~,idx] = min(d);
    rnew(:,3) = b(:,idx);
    clear d idx
    % find P1:
    % trying out another way 170414:
    % old way:
    %     rnew(:,2) = 2*(bezierval(r,.5) - .25*rnew(:,1) - .25*rnew(:,3));
    % new way
    rnew(:,2) = 0.5*(rnew(:,1)-r(:,1)) + 0.5*(rnew(:,3)-r(:,3)) + r(:,2);
    
else
    % it's the first frame, so using this to define s0
    rnew = r;
end

[~,idx] = min((s-sfp-sstar).^2);
tstar = t(idx);
clear idx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function theta = base_angle(r,t)

% order = size(p,2)-1;
% pds = zeros(2,order);
% pds(1,:) = p(1,1:order).*(order:-1:1);
% pds(2,:) = p(2,1:order).*(order:-1:1);
%
% dxds = polyval(pds(1,:),s);
% dyds = polyval(pds(2,:),s);
% theta = atan2(-dyds,dxds)*(180/pi);
% clear order pds dxds dyds

dBdt = bezierdtval(r,t);
theta = atan2(-dBdt(2,:),dBdt(1,:))*(180/pi);
clear dBdt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kappa = curvature(r,t)

% order = size(p,2)-1;
% pds = zeros(2,order);
% pds(1,:) = p(1,1:order).*(order:-1:1);
% pds(2,:) = p(2,1:order).*(order:-1:1);
% pds2(1,:) = pds(1,1:order-1).*(order-1:-1:1);
% pds2(2,:) = pds(2,1:order-1).*(order-1:-1:1);
%
% dxds = polyval(pds(1,:),s);
% dyds = polyval(pds(2,:),s);
% d2xds2 = polyval(pds2(1,:),s);
% d2yds2 = polyval(pds2(2,:),s);
%
% kappa = dxds.*d2yds2 - dyds.*d2xds2;
% kappa = kappa ./ (dxds.^2+dyds.^2).^(3/2);
%
% clear order pds pds2 dxds dyds d2xds2 d2yds2

dBdt = bezierdtval(r,t);
d2Bdt2 = bezierdt2val(r,t);
kappa = (dBdt(1,:).*d2Bdt2(2,:)-dBdt(2,:).*d2Bdt2(1,:)) ./ (dBdt(1,:).^2+dBdt(2,:).^2).^(3/2);
clear dBdt d2Bt2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function bez = bezierfit(r,order)
% % Fit bezier curve of order 'order' to control points r
% % Initial conditions:
% switch order
%     case 1
%         bez0 = zeros(2,2);
%         bez0(:,1) = r(:,1);
%         bez0(:,2) = r(:,end);
%     case 2
%         bez0 = zeros(2,3);
%         bez0(:,1) = r(:,1);
%         bez0(:,3) = r(:,end);
%         m = ceil(size(r,2)/2);
%         bez0(:,2) = r(:,m);
%     otherwise
%         error('write more code!')
% end
% clear m
% % Minimise regression error:
% [bez,Emin,exitflag,output] = fminunc(@(bez) Bfun(bez,r),bez0,....
%                 optimset('Display','off'));%'TolX',.1,'TolFun',.01));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [E,gradE] = Bimfun(z, r, n, im, dt, sigma)
% Assume order 2 bezier curve
bez = zeros(2,3);
bez(:,1) = r(:,1) + z(1)*n(:,1);
bez(:,2) = z(2:3)';
bez(:,3) = r(:,3) + z(4)*n(:,2);
Ereg = 0.5 * sigma * sum(sum((bez-r).^2));
gradEreg = zeros(4,1);
gradEreg(1) = sigma * (bez(:,1)-r(:,1))'*n(:,1);
gradEreg(2:3) = sigma * (bez(:,2)-r(:,2));
gradEreg(4) = sigma * (bez(:,3)-r(:,3))'*n(:,2);
clear r
t = 0:dt:1;
r = bezierval(bez,t);
% dr = diff(r,1,2);
% ds = sqrt(sum(dr.^2,1));
% s = [0 cumsum(ds)];
% snew = linspace(s(2),s(end-1),length(s));
% clear dr ds
% % figure, plot(s,t,'.',snew,0,'.'),pause
% tnew = interp1(s,t,snew);
% r = bezierval(bez,tnew);
% t = tnew;
% clear s snew tnew

% find the points on the bezier curve that are internal to the image:
w = size(im.s,2);
h = size(im.s,1);
gdpt = find((r(1,:)>=1)&(r(1,:)<=w)&(r(2,:)>=1)&(r(2,:)<=h));
clear w h

%im.s
E = dt*sum(interp2(im.s,r(1,gdpt),r(2,gdpt)));

% fprintf('Eim=%.2f Ereg=%.2f Etot=%.2f\n',E,Ereg,E+Ereg)
E = E + Ereg;
gradE = zeros(4,1);
idx = interp2(im.dx,r(1,gdpt),r(2,gdpt));
idy = interp2(im.dy,r(1,gdpt),r(2,gdpt));
gradE(1) = dt*sum(idx.*((1-t(gdpt)).^2)*n(1,1)+idy.*((1-t(gdpt)).^2)*n(2,1));
gradE(2) = dt*sum(idx.*(2*(1-t(gdpt)).*t(gdpt)));
gradE(3) = dt*sum(idy.*(2*(1-t(gdpt)).*t(gdpt)));
gradE(4) = dt*sum(idx.*(t(gdpt).^2)*n(1,2)+idy.*(t(gdpt).^2)*n(2,2));
clear gdpt idx idy
gradE = gradE + gradEreg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = bezierval(bez,t)
% Evaluate bezier curve with parameters bez(:,1),bez(:,2),... at points t
% size(t) = [1,N], where N is number of points at which to evaluate
% function.  Typically, t = [0,1]
% size(bez) = [2,order+1]
order = size(bez,2)-1;
switch order
    case 1
        p0 = bez(:,1);
        p1 = bez(:,2);
        b = p0 + (p1-p0)*t;
    case 2
        p0 = bez(:,1);
        p1 = bez(:,2);
        p2 = bez(:,3);
        b = p0*(1-t).^2 + 2*p1*((1-t).*t) + p2*t.^2;
    otherwise
        error('Write more code!')
end
clear p0 p1 p2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dBdt = bezierdtval(bez,t)
% Evaluate 1st deriv of bezier curve wrt t, with parameters bez(:,1),bez(:,2),... at points t
% size(t) = [1,N], where N is number of points at which to evaluate function.  Typically, t = [0,1]
% size(bez) = [2,order+1]
order = size(bez,2)-1;
switch order
    case 2
        p0 = bez(:,1);
        p1 = bez(:,2);
        p2 = bez(:,3);
        %         dBdt = -p0*2*(1-t) + 2*p1*(1-2*t) + 2*p2*t;
        dBdt = 2*(p0-2*p1+p2)*t + 2*(-p0+p1)*ones(size(t));
    otherwise
        error('Write more code!')
end
clear p0 p1 p2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dB2dt2 = bezierdt2val(bez,t)
% Evaluate 2nd deriv of bezier curve wrt t, with parameters bez(:,1),bez(:,2),... at points t
% size(t) = [1,N], where N is number of points at which to evaluate function.  Typically, t = [0,1]
% size(bez) = [2,order+1]
order = size(bez,2)-1;
switch order
    case 2
        p0 = bez(:,1);
        p1 = bez(:,2);
        p2 = bez(:,3);
        dB2dt2 = 2*(p0-2*p1+p2) * ones(size(t));
    otherwise
        error('Write more code!')
end
clear p0 p1 p2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,s] = fit_curve(r,order)
% Fit polynomial curve (parameterised by 's') to the control points:
dr = sqrt(sum(diff(r,1,2).^2,1));
s = [0 cumsum(dr)];
p = zeros(2,order+1);
p(1,:) = polyfit(s,r(1,:),order);
p(2,:) = polyfit(s,r(2,:),order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fp,fidx,tfp] = find_follicle(r,dt,folliclemask)

% compute arc length parameter 's':
tneg = -1:dt:0;
w = bezierval(r,tneg);
dw = diff(w,1,2);
ds = sqrt(sum(dw.^2,1));
% sneg = [0 cumsum(ds)] - sum(ds); % s(1)=0 ~ P0
clear dw ds

% w = zeros(2,length(tneg));
% w(1,:) = polyval(p(1,:),sneg);
% w(2,:) = polyval(p(2,:),sneg);

dist = zeros(length(tneg),size(folliclemask,2));
for i = 1:length(tneg)
    for j = 1:size(folliclemask,2)
        dist(i,j) = sum((w(:,i)-folliclemask(:,j)).^2,1);
    end
end
clear i j
[m,fidx] = min(dist');  % min over folliclemask for fixed s
% locate first minimum (in order of decreasing t)
% sidx = find(diff(m)>0,1);
[~,sidx] = min(m);      % min of the min (over s)
%     fidx = fidx(sidx);      % pt on f closest to w
% now we have Order(0) follicle position estimate:
fp = [w(1,sidx); w(2,sidx)];
clear dist m tmp

% use linear interpolation to increase accuracy of fp estimate.
% find line that locally fits folliclemask (p2) and line that locally fits
% whisker contour (p1):
sidx = min([size(w,2)-1 sidx]);
p1 = polyfit(w(1,sidx-1:sidx+1),w(2,sidx-1:sidx+1),1);
fidx = fidx(sidx);
p2 = polyfit(folliclemask(1,fidx-1:fidx+1),folliclemask(2,fidx-1:fidx+1),1);
% now solve for the point where p1 and p2 intersect:
M = [p1(1) -1; p2(1) -1];
fp = -inv(M)*[p1(2);p2(2)];
clear p1 p2 M

% find the 't' value of fp.  Do this by minimising ||f(t)-fp||^2 wrt t,
% starting from initial condition 't0':
t0 = tneg(sidx);
options = optimset('display','off');
tfp = fminunc(@(t) tmpfun(t, r, fp), t0, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = tmpfun(t,r,fp)
fphat = zeros(2,1);
% fphat(1) = polyval(p(1,:),s);
% fphat(2) = polyval(p(2,:),s);
fphat = bezierval(r,t);
E = sum((fp-fphat).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [folliclemask] = protract(image, roi, folliclemask_old, fpidx_old, theta, sigma)

% Filter and differentiate:

ims.raw = image(roi.y,roi.x,1);
ims.w = size(ims.raw,2);
ims.h = size(ims.raw,1);
ims.mfilt = medfilt2(ims.raw,[5 5]);
hsize = 3*[sigma sigma];
filter = fspecial('gaussian', hsize, sigma);
ims.gfilt = imfilter(ims.mfilt,filter,'replicate');
clear sigma hsize filter

% keyboard

% ims.Del1.x = [zeros(size(ims.gfilt,1),1) diff(ims.gfilt,1,2)];
% ims.Del1.y = [zeros(1,size(ims.gfilt,2)); diff(ims.gfilt,1,1)];
% ims.Del2.x = [zeros(size(ims.gfilt,1),1) diff(ims.gfilt,2,2) zeros(size(ims.gfilt,1),1)];
% ims.Del2.x = [zeros(1,size(ims.gfilt,2)); diff(ims.gfilt,2,1); zeros(1,size(ims.gfilt,2))];
[ims.Del1.x,ims.Del1.y] = gradient(ims.gfilt);

% Specify the desired gradient measure:
if ~isempty(folliclemask_old)
    % Revised code (290414)
    % Go back to using normal to snout mask instead of tangent to whisker
    % solution:
    %     % use normal to snout mask in previous frame to define gradient
    %     % direction:
    tvec = folliclemask_old(:,fpidx_old) - folliclemask_old(:,fpidx_old+1);
    tvec = tvec / norm(tvec);
    % So following lines (290414) commented out:
    %     % use tangent to whisker solution in previous frame to define gradient
    %     % direction
    %     nvec = [-cosd(theta); sind(theta)];
    %     nvec = nvec/norm(nvec);
    
else
    % no previous frame, so use default:
%     nvec = -[1; -1]/sqrt(2);
    % 200217 automatic default
    bw =  edge(ims.gfilt);  %edge detect smooth image to find snout
    [y,x] = find(bw);   % points along the edge
    p = polyfit(x,y,1); % linear fit to the snout edge
    tvec = [200 polyval(p,200)]' - [100 polyval(p,100)]';
    tvec = tvec / norm(tvec);
    clear bw y x p
end
R = [0 1;-1 0];
nvec = -R*tvec;
clear tvec R
ims.grad = nvec(1)*ims.Del1.x+nvec(2)*ims.Del1.y;

snout.x = 1:length(roi.x);
[~,snout.y] = min(ims.grad);
% Take every Nth pixel:
N = 5;
snout.x = snout.x(1:N:end);
snout.y = snout.y(1:N:end);
Npts = length(snout.x);
clear N

% Quality control on snout detection, debugging:
%     figure
%     subplot(2,2,1), imagesc(ims.raw), colormap gray
%     hold on, plot(snout.x,snout.y,'r.-')
%     subplot(2,2,2), imagesc(ims.gfilt)
%     subplot(2,2,3), imagesc(ims.mfilt)
%     subplot(2,2,4), imagesc(ims.grad)
%     pause

% The above produces a snout estimate that is outside the face
% slightly.  Warp the contour in a bit in direction of local normals

z = 5;
folliclemask(1,:) = snout.x + z*nvec(1,:);
folliclemask(2,:) = snout.y + z*nvec(2,:);
clear nvec

%  Debug (may not work 200217):
%     figure
%     imagesc(ims.raw), colormap gray, colorbar
%     hold on,
%     plot(folliclemask(1,:),folliclemask(2,:),'m-',trackmask(1,:),trackmask(2,:),'g-')

% Align masks to complete image:
folliclemask(1,:) = folliclemask(1,:) + roi.x(1);
folliclemask(2,:) = folliclemask(2,:) + roi.y(1);
clear roi

% The extreme points are probably junk.  Cut them:
L = 1;
folliclemask(:,1:L) = [];
folliclemask(:,end-L+1:end) = [];
clear L

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

frameidx = handles.frameidx;

change_frame = 0;
change_cpt = 0;
nudge = 0;
dnudge = 0.5;

if length(eventdata.Modifier) == 0
    switch eventdata.Key
        case 'leftarrow'
            frameidx = modsubtract(frameidx,1,handles.nframes); change_frame = 1;
        case 'rightarrow'
            frameidx = modadd(frameidx,1,handles.nframes); change_frame = 1;
        case 'uparrow'
            cpt = handles.current_pt;
            cpt = rem(cpt+1,handles.n_control_pts+1); change_cpt = 1;
        case 'downarrow'
            cpt = handles.current_pt;
            cpt = (cpt-1) + (cpt==0)*(handles.n_control_pts+1); change_cpt = 1;
        case 'm'
            handles = move_control_point_Callback(handles);
    end
elseif strmatch(eventdata.Modifier{1},'shift')
    switch eventdata.Key
        case 'leftarrow'
            frameidx = modsubtract(frameidx,10,handles.nframes); change_frame = 1;
        case 'rightarrow'
            frameidx = modadd(frameidx,10,handles.nframes); change_frame = 1;
    end
elseif strmatch(eventdata.Modifier{1},'control')
    switch eventdata.Key
        case 'leftarrow'
            frameidx = modsubtract(frameidx,100,handles.nframes); change_frame = 1;
        case 'rightarrow'
            frameidx = modadd(frameidx,100,handles.nframes); change_frame = 1;
    end
elseif strmatch(eventdata.Modifier{1},'alt')
    switch eventdata.Key
        case 'leftarrow'
            nudge = -1;
        case 'rightarrow'
            nudge = 1;
        case 'uparrow'
            nudge = -2;
        case 'downarrow'
            nudge = 2;
    end
    
end

if change_frame
    % plot new image in the baseframe.  if it exists, superimpose tracking solution.
    %     display_frame(frameidx)
    handles.frame = load_frame(handles.video,frameidx);
    titlestring = sprintf('%d',frameidx);
    %     if isfield(handles,'goodsolution')
    %         if handles.goodsolution(frameidx)
    %             titlestring = sprintf('%d (%.1f)',frameidx,handles.Emin_all(frameidx));
    %         end
    %     end
    if isfield(handles,'tracked')
        if handles.tracked(frameidx) && isfield(handles,'Emin')
            titlestring = sprintf('%d (%.1f)',frameidx,handles.Emin_all(frameidx));
        end
    end
    set(handles.frameidx_display,'string',titlestring)
    clear titlestring
    %     im = double(handles.frame(:,:,1)) - handles.meanframe(:,:,1);
    handles.frameidx = frameidx;
    %     figure(handles.figure1)
    axes(handles.video_frame)
    cla
    %     imagesc(im), colormap gray
    image(handles.frame)
    if isfield(handles,'rall')
        hold on
        rall = handles.rall;
        dt = 0.05;
        t = 0:dt:1;
        r = squeeze(rall(frameidx,:,:));
        b = bezierval(r,t);
        plot(b(1,:),b(2,:),'r-',r(1,:),r(2,:),'r.',handles.fp_all(1,handles.frameidx),handles.fp_all(2,handles.frameidx),'y.')
        % how do manual editing of control points?
        % could use up/downarrows to cycle between control points.  default is no
        % selection.  return key -> ginput -> updated control point. update display
        clear dt t r b
    end
    clear frameidx
end

if change_cpt
    axes(handles.aux_frame)
    cla
    imagesc(double(handles.frame(:,:,1))-handles.meanframe(:,:,1)),
    colormap gray
    title(handles.frameidx)
    hold on
    %     size(handles.rall),handles.frameidx,pause
    r = squeeze(handles.rall(handles.frameidx,:,:));
    dt = 0.02;
    t = 0:dt:1;
    b = bezierval(r,t);
    plot(b(1,:),b(2,:),'r-',r(1,:),r(2,:),'r.')
    if cpt~=0
        plot(r(1,cpt),r(2,cpt),'ow')
    end
    handles.current_pt  = cpt;
    clear r t b
end

if nudge ~= 0
    cpt = handles.current_pt;
    axes(handles.aux_frame)
    cla
    imagesc(double(handles.frame(:,:,1))-handles.meanframe(:,:,1)),
    colormap gray
    title(handles.frameidx)
    hold on
    r = squeeze(handles.rall(handles.frameidx,:,:));
    plot(r(1,:),r(2,:),'r.',r(1,cpt),r(2,cpt),'wo')
    %     order = handles.pmtrs.order;
    rnew = r;
    switch abs(nudge)
        case 1
            % horizontal
            rnew(1,cpt) = r(1,cpt) + nudge*dnudge;
        case 2
            % vertical
            rnew(2,cpt) = r(2,cpt) + sign(nudge)*dnudge;
        otherwise
            error
    end
    clear r
    plot(rnew(1,cpt),rnew(2,cpt),'.w')
    
    % Compute energy of the solution:
    z = [0 rnew(:,2)' 0];
    % So, first compute unit normal vectors at B(t=0) = P0 and B(t=1) = P2:
    tvec = bezierdtval(rnew,[0 1]);
    n = sqrt(sum(tvec.^2,1));
    tvec = tvec ./ (ones(2,1)*n);
    clear n
    R = [0 -1;1 0]; % 90deg rotation matrix
    nvec = R*tvec;       % Normal to start and end of contour
    n = sqrt(sum(nvec.^2,1));
    nvec = nvec ./ (ones(2,1)*n);    % unit normals
    clear tvec R n
    im.raw = double(handles.frame(:,:,1)) - handles.meanframe(:,:,1);
    im.s = imfilter(im.raw,handles.gaussian);
    [im.dx,im.dy] = gradient(im.s);
    [E, ~] = Bimfun(z, rnew, nvec, im, handles.dt, handles.sigma_prior);
    handles.Emin_all(handles.frameidx) = E;
    title(sprintf('%d:E=%.2f',handles.frameidx,E))
    clear  z n im E
    
    % Find intersection between whisker polynomial and folliclemask
    roi.x = round(handles.fp_all(1,handles.frameidx-1)) + [-40:40];
    roi.y = round(handles.fp_all(2,handles.frameidx-1)) + [-40:40];
    roi.x = roi.x((roi.x>=1)&(roi.x<=size(handles.frame,2)));
    roi.y = roi.y((roi.y>=1)&(roi.y<=size(handles.frame,1)));
    
    [p,~] = fit_curve(rnew,handles.polyorder);
    theta = base_angle(p,0);
    clear p
    folliclemask = protract(double(handles.frame(:,:,1)),roi,...
        handles.folliclemask{handles.frameidx-1},handles.fpidx(handles.frameidx-1),theta,handles.snout_sigma);
    clear roi theta
    
    % Postprocess:
    %     [~, p, fp, fpidx, sfp, s] = postproc(rnew, folliclemask, handles.s0, handles.dt, handles.polyorder);
    r = rnew;
    [rnew, fp, fpidx, tfp, sfp, s, tstar, sp2] = postproc(rnew, folliclemask, handles.s0, handles.dt, handles.sstar);
    
    
    doplot = true;  % should really be made consistent with 'doplot' in the main loop
    if doplot
        b = zeros(2,length(s)+1);
        %         b(1,:) = polyval(p(1,:),[sfp s]);
        %         b(2,:) = polyval(p(2,:),[sfp s]);
        t = [-tfp:handles.dt:1];
        b = bezierval(rnew,[tfp t]);
        clear t
        plot(b(1,:),b(2,:),'r',r(1,:),r(2,:),'r.',fp(1),fp(2),'y.')
        plot(rnew(1,:),rnew(2,:),'msq'),
        %         plot(polyval(p(1,:),handles.sstar),polyval(p(2,:),handles.sstar),'go')
        clear b
    end
    clear sfp s doplot
    
    handles.rall(handles.frameidx,:,:) = rnew;
    handles.fp_all(:,handles.frameidx) = fp;
    handles.folliclemask{handles.frameidx} = folliclemask;
    handles.fpidx_all(:,handles.frameidx) = fpidx;
    handles.tracked(handles.frameidx) = true;
    handles.kappa_all(handles.frameidx) = curvature(rnew,tstar);
    handles.theta_all(handles.frameidx) = base_angle(rnew,tstar);
    clear fp fpidx
    
    % Save to disk:
    rall = handles.rall;
    fp_all = handles.fp_all;
    kappa_all = handles.kappa_all;
    theta_all = handles.theta_all;
    pmtrs.sstar = handles.sstar;
    pmtrs.s0 = handles.s0;
    pmtrs.sigma_prior = handles.sigma_prior;
    pmtrs.subtract_image_mean = handles.subtract_image_mean;
    save(handles.trfname,'rall','kappa_all','theta_all','fp_all','pmtrs')
    clear rall kappa_all theta_all fp_all pmtrs
    
end


guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_move_control_point.
function pushbutton_move_control_point_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_move_control_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

error('Update the curve to postprocess contour')

handles = move_control_point_Callback(handles);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = move_control_point_Callback(handles)

cpt = handles.current_pt;
if cpt==0
    return
end

axes(handles.aux_frame)
cla
imagesc(double(handles.frame(:,:,1))-handles.meanframe(:,:,1)),
colormap gray
title(handles.frameidx)
hold on

r = squeeze(handles.rall(handles.frameidx,:,:));
plot(r(1,:),r(2,:),'r.-')
plot(r(1,cpt),r(2,cpt),'wo')
[x,y] = ginput(1);
plot(x,y,'w+')
r(:,cpt) = [x y]';
plot(r(1,:),r(2,:),'r--o')

D = size(r,2);
handles.rall(handles.frameidx,:,:) = r;

clear cpt r ddr

% guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in radiobutton_kinematics_test.
function radiobutton_kinematics_test_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_kinematics_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_kinematics_test

% --- Executes on button press in radiobutton_continuous_tracking.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function radiobutton_continuous_tracking_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_continuous_tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_continuous_tracking

% dummy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in radiobutton_subtract_image_mean.
function radiobutton_subtract_image_mean_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_subtract_image_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_subtract_image_mean

% Compute meanframe
meanframe = zeros(handles.video.height,handles.video.width,3);
if get(hObject,'Value')==true
    set(handles.text_status,'BackgroundColor',handles.red)
    set(handles.text_status,'String','BUSY')
    pause(.01)
    space = 10;
    for frameidx = 1:space:handles.nframes
        meanframe = meanframe + double(load_frame(handles.video,frameidx))/(handles.nframes/space);
    end
    clear space frameidx
    set(handles.text_status,'BackgroundColor',handles.green)
    set(handles.text_status,'String','idle')
end
fn = [handles.fname(1:end-3) 'meanframe'];
save(fn,'meanframe')
handles.meanframe = meanframe;
clear meanframe fn

axes(handles.aux_frame)
imagesc(double(handles.frame(:,:,1))-handles.meanframe(:,:,1))
colormap gray
hold on

% Update the energy threshold:
if get(handles.radiobutton_subtract_image_mean,'Value')
    handles.energy_threshold = handles.default_energy_threshold_subtract_image_mean;
else
    handles.energy_threshold = handles.default_energy_threshold_do_not_subtract_image_mean;
end
set(handles.edit_energy_threshold,'String', handles.energy_threshold);

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.trfname)
handles.trfname

handles.rall = zeros(handles.nframes,2,handles.n_control_pts);
handles.kappa_all = zeros(1,handles.nframes);
handles.theta_all = zeros(1,handles.nframes);
handles.fp_all = zeros(2,handles.nframes);
handles.fpidx_all = zeros(1,handles.nframes);
handles.Emin_all = zeros(1,handles.nframes);
handles.tracked = zeros(1,handles.nframes);

handles.current_pt = 0;
handles.s0 = [];

handles.frameidx = handles.startframe;
handles.frame = load_frame(handles.video,handles.frameidx);
set(handles.frameidx_display,'string',num2str(handles.frameidx))

axes(handles.video_frame)
cla
image(handles.frame)

axes(handles.aux_frame)
cla
imagesc(double(handles.frame(:,:,1))-handles.meanframe(:,:,1))
% title(frameidx)
colormap gray
hold on

handles.continuous_tracking = false;
set(handles.radiobutton_continuous_tracking,'value',handles.continuous_tracking)
% handles.subtract_image_mean = true;
% set(handles.radiobutton_subtract_image_mean,'value',handles.subtract_image_mean)

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_plot.
function pushbutton_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rall = handles.rall;

nframes = size(rall,1);
frames = 1:nframes;
idx = find(handles.tracked);
% D = size(rall,3);

% idx
% size(handles.tracked)
% size(handles.theta_all)

% % compute base angle
% theta = zeros(length(idx),D-1);
% for d = 1:D-1
%     dy = rall(idx,2,d+1)-rall(idx,2,1);
%     dx = rall(idx,1,d+1)-rall(idx,1,1);
%     theta(idx,d) = atan2(-dy,dx)*(180/pi);
% end
% clear d dx dy

% % compute curvature: kappa=(x^'y^('')-y^'x^(''))/((x^('2)+y^('2))^(3/2)).
% tmp = diff(rall,1,3);
% dt = reshape(sum(tmp.^2,2),nframes,D-1);    % differential lengths along the contours
% clear tmp
% % t = cumsum(dt,2);   % distance along the contours
% dxdt = reshape(diff(rall(:,1,:),1,3),nframes,D-1) ./ dt;
% dxdt = dxdt(:,2:end);
% d2xdt2 = reshape(diff(rall(:,1,:),2,3),nframes,D-2) ./ dt(:,2:end).^2;
% dydt = reshape(diff(rall(:,2,:),1,3),nframes,D-1) ./ dt;
% dydt = dydt(:,2:end);
% d2ydt2 = reshape(diff(rall(:,2,:),2,3),nframes,D-2) ./ dt(:,2:end).^2;
% kappa = dxdt.*d2ydt2 - dydt.*d2xdt2;
% kappa = kappa ./ (dxdt.^2+dydt.^2).^(3/2);

figure
spr = 4;
h(1) = subplot(spr,1,1);
% plot(frames(idx),theta(idx,:),'k',frames(idx),handles.theta_all(idx),'r')
plot(frames(idx),handles.theta_all(idx),'.')
title('Base angle'), ylabel('Deg')%, xlabel('Frame')
h(2) = subplot(spr,1,2);
% plot(frames(idx),medfilt1(kappa(idx,:),5))
% plot(frames(idx),kappa(idx,:),'k',frames(idx),handles.kappa_all(idx),'r')
plot(frames(idx),handles.kappa_all(idx),'.')
title('Curvature'), ylabel('pixels^{-1}')
h(3) = subplot(spr,1,3);
plot(frames(idx),handles.fp_all(1,idx),'.')
title('Follicle X position'), xlabel('Frame'), ylabel('pixels')
h(4) = subplot(spr,1,4);
plot(frames(idx),handles.fp_all(2,idx),'.')
title('Follicle Y position'), xlabel('Frame'), ylabel('pixels')

linkaxes(h,'x')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_sigma_prior_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sigma_prior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sigma_prior as text
%        str2double(get(hObject,'String')) returns contents of edit_sigma_prior as a double

handles.sigma_prior = str2double(get(hObject,'String'));
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function edit_sigma_prior_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sigma_prior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_frame_interval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frame_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_frame_interval as text
%        str2double(get(hObject,'String')) returns contents of edit_frame_interval as a double

handles.frame_interval = str2double(get(hObject,'String'));
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function edit_frame_interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_frame_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_energy_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_energy_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_energy_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_energy_threshold as a double

handles.energy_threshold = str2double(get(hObject,'String'));
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function edit_energy_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_energy_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_auto_initialise.
function radiobutton_auto_initialise_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_auto_initialise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_auto_initialise


% --- Executes on selection change in popupmenu_choose_auto_initialise_file.
function popupmenu_choose_auto_initialise_file_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_auto_initialise_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_choose_auto_initialise_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_choose_auto_initialise_file


% --- Executes during object creation, after setting all properties.
function popupmenu_choose_auto_initialise_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_auto_initialise_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

pwd, disp('asd')
ff = dir('*.tr')
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose init file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(hObject,'String',string_list);
clear ff nr_files string_list i


% --- Executes on button press in pushbutton_start_batch.
function pushbutton_start_batch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.radiobutton_auto_initialise,'Value',true)
set(handles.radiobutton_continuous_tracking,'Value',true);
menu = cellstr(get(handles.popupmenu_choose_batch_file,'String'));
batchfn = menu{get(handles.popupmenu_choose_batch_file,'Value')};
batch_fid = fopen(batchfn,'rt');
clear menu
reportfn = 'report.rep';
report_fid = fopen(reportfn,'at');
time = fix(clock);
time = sprintf('%d:%d:%d',time(4:6));
fprintf(report_fid,'\nStarting batch %s at %s on %s...\n',batchfn,time,date);
clear time batch_fn
stopbatch = false;
while (batch_fid && ~stopbatch)
    handles.fname = fgetl(batch_fid);
    fprintf(report_fid,' Opening video %s\n',handles.fname);
    % stuff
    % choose the most recent, appropriate tr file in working directory, if there is one, as the auto init file
    ff = dir(['*.tr']);
    if length(ff)>0
        nfiles=size(ff,1);
        datenum = zeros(1,nfiles);
        for i=1:nfiles
            datenum(i) = ff(i).datenum;
        end
        clear i nfiles
        [~,idx] = max(datenum);
        fn = ff(idx).name;
        set(handles.popupmenu_choose_auto_initialise_file,'String',fn);
        clear datenum idx ff fn
        clear ff nr_files string_list i
    else
        % there are no suitable tr files in current directory.  use default.
    end
    
    handles = initialise_new_video(handles);
    handles = initialise_new_track(handles);
    [handles,fnoutput] = fit_snakes(handles);
    switch fnoutput
        case -1
            % should only happen if user switched off the
            % continuous_tracking radio button
            stopbatch = true;
            fprintf(report_fid,' Track failure: output=%d\n',fnoutput);
        case 1
            % fit_snakes successfully reached end of video
            fprintf(report_fid,' Tracker reached end of video without error\n');
        case 2
            % Tracking failed somewhere in the video
            stopbatch = true;
            fprintf(report_fid,' Track failure: output=%d\n',fnoutput);
        otherwise
            error('Unrecognised output from fit_snakes %d', fnoutput)
    end
    if ~get(handles.radiobutton_continuous_tracking,'Value')
        % ie user interrupted the track
        stopbatch = true;
        fprintf(report_fid,' Track failure: user abort\n');
    end
end
fclose(batch_fid);
fclose(report_fid);
clear batch_fid report_fid

guidata(hObject, handles);


% --- Executes on selection change in popupmenu_choose_batch_file.
function popupmenu_choose_batch_file_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_batch_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_choose_batch_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_choose_batch_file

ff = dir('*.bat');
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose batch file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(hObject,'String',string_list);
clear ff nr_files string_list i


% --- Executes during object creation, after setting all properties.
function popupmenu_choose_batch_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_batch_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

ff = dir('*.bat');
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose batch file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(hObject,'String',string_list);
clear ff nr_files string_list i



function edit_miumax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_miumax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_miumax as text
%        str2double(get(hObject,'String')) returns contents of edit_miumax as a double

handles.miumax = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_miumax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_miumax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z = modsubtract(x,y,n)
% modular subtraction.
% eg if x=5,y=2,n=9, z = 3 (as normal)
% eg  if x=2,y=2,n=9, z = 9 (wrap-around)

z = x-y;
if z <= 0
    z = z+n;
end


function z = modadd(x,y,n)
% modular addition
% eg if x=5,y=2,n=9, z = 7 (as normal)
% eg  if x=8,y=2,n=9, z = 1 (wrap-around)

z = x+y;
if z > n
    z = z - n;
end



function handles = initialise_new_track(handles)

% if tracking data file exists, load both data and pmtrs:
if exist(handles.trfname,'file')
    load(handles.trfname,'rall','kappa_all','theta_all','fp_all','Emin_all','pmtrs','-mat')
    handles.rall = rall;
    handles.n_control_pts = size(handles.rall,3);
    if handles.n_control_pts~=3
        error('Bad tr file: N control points must be 3')
    end
    handles.tracked = squeeze(sum(abs(rall(:,:,1)),2))>0;
    handles.kappa_all = kappa_all;
    handles.theta_all = theta_all;
    handles.fp_all = fp_all;
    handles.Emin_all = Emin_all;
    clear rall kappa_all theta_all fp_all Emin_all
    handles.sstar = pmtrs.sstar;
    handles.s0 = pmtrs.s0;
    handles.sigma_prior = pmtrs.sigma_prior;
    set(handles.edit_sigma_prior,'String',handles.sigma_prior);
    handles.subtract_image_mean = pmtrs.subtract_image_mean;
    set(handles.radiobutton_subtract_image_mean,'value',handles.subtract_image_mean)
    handles.energy_threshold = 140%pmtrs.energy_threshold;
    handles.adaptive_energy_threshold = pmtrs.adaptive_energy_threshold;
    set(handles.radiobutton_adaptive_energy_threshold,'value',handles.adaptive_energy_threshold)
    set(handles.edit_energy_threshold,'String',handles.energy_threshold);
    if isfield(pmtrs,'miumax')
        handles.miumax =  pmtrs.miumax;
        set(handles.edit_miumax,'String',handles.miumax);
    else
        handles.miumax = str2double(get(handles.edit_miumax,'String'));
    end
    if isfield(pmtrs,'snout_sigma')
        handles.snout_sigma = pmtrs.snout_sigma;
        set(handles.edit_snout_sigma,'String',handles.snout_sigma);
    end
    if isfield(pmtrs,'snout_roi')
        handles.snout_roi = pmtrs.snout_roi;
    end   
    clear pmtrs
else    % if no tracking data file exists, initialise.
    handles.n_control_pts = 3;
    handles.tracked = zeros(1,handles.nframes);
    handles.rall = zeros(handles.nframes,2,handles.n_control_pts);
    handles.kappa_all = zeros(1,handles.nframes);
    handles.theta_all = zeros(1,handles.nframes);
    handles.fp_all = zeros(2,handles.nframes);
    handles.fpidx_all = zeros(1,handles.nframes);
    handles.Emin_all = zeros(1,handles.nframes);
    handles.sstar = 20; % distance from follicle at which to evaluate angle, curvature
    handles.sigma_prior = str2double(get(handles.edit_sigma_prior,'String'));
    handles.subtract_image_mean = get(handles.radiobutton_subtract_image_mean,'value');
    handles.energy_threshold = str2double(get(handles.edit_energy_threshold,'String'));
    handles.energy_threshold_tolerance = str2double(get(handles.edit_energy_threshold_tolerance,'String'));
    handles.miumax = str2double(get(handles.edit_miumax,'String'));
    handles.snout_sigma = str2double(get(handles.edit_snout_sigma,'String'));
    handles.adaptive_energy_threshold = get(handles.radiobutton_adaptive_energy_threshold,'value');
    handles.current_pt = 0;
    handles.s0 = [];


end
handles.polyorder = handles.n_control_pts-1;

handles.frameidx = handles.startframe;
handles.frame = load_frame(handles.video,handles.frameidx);
set(handles.frameidx_display,'string',num2str(handles.frameidx))

% plot tracked contour for current frame
axes(handles.video_frame)
cla
image(handles.frame)
hold on

rall = handles.rall;
dt = 0.05;
t = 0:dt:1;
r = squeeze(rall(handles.frameidx,:,:));
b = bezierval(r,t);
plot(b(1,:),b(2,:),'r-')%,r(1,:),r(2,:),'r.',handles.fp_all(1,handles.frameidx),handles.fp_all(2,handles.frameidx),'y.')
clear dt t r b

% refresh the aux frame too
axes(handles.aux_frame)
cla
if handles.subtract_image_mean
    load([handles.fname(1:end-3) 'meanframe'],'meanframe','-mat')
else
    meanframe = zeros(handles.video.height,handles.video.width,3);
end

imagesc(double(handles.frame(:,:,1))-meanframe(:,:,1))
colormap gray
% hold on
clear meanframe

% handles.continuous_tracking = false;
% set(handles.radiobutton_continuous_tracking,'value',handles.continuous_tracking)

% guidata(hObject, handles);



function edit_specify_tr_file_Callback(hObject, eventdata, handles)
% hObject    handle to edit_specify_tr_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_specify_tr_file as text
%        str2double(get(hObject,'String')) returns contents of edit_specify_tr_file as a double

handles.trfname = get(hObject,'String');
handles = initialise_new_track(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_specify_tr_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_specify_tr_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_set_snout_roi.
function pushbutton_set_snout_roi_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set_snout_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'snout_roi')
    roi_default = [handles.snout_roi.x([1 end]) handles.snout_roi.y([1 end])];
else
    roi_default = [1 size(handles.frame,2) 1 size(handles.frame,1)];
end

roi = setroi_gui(handles.frame,roi_default);
clear roi_default

handles.snout_roi.x = roi(1):roi(2);
handles.snout_roi.y = roi(3):roi(4);
clear roi

guidata(hObject, handles);



function edit_snout_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_snout_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_snout_sigma as text
%        str2double(get(hObject,'String')) returns contents of edit_snout_sigma as a double

handles.snout_sigma = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_snout_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_snout_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_energy_threshold_tolerance_Callback(hObject, eventdata, handles)
% hObject    handle to edit_energy_threshold_tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_energy_threshold_tolerance as text
handles.energy_threshold_tolerance =  str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_energy_threshold_tolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_energy_threshold_tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_adaptive_energy_threshold.
function radiobutton_adaptive_energy_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_adaptive_energy_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.adaptive_energy_threshold = get(hObject,'Value');% returns toggle state of radiobutton_adaptive_energy_threshold
guidata(hObject, handles);
