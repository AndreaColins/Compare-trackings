function [barPos, image] = poletracker(fname, radius,drawing,interval)
tic
% helper function for extracting pole location
% Created RSP 040914
%
% Developed by Mathew Evans 14.10.2014
%
% AIM: use this script to develop pole tracking code and to find
% parameters, eg for the pole radius

%close all


% the video:
handles.fname = fname;

% assumed pole radius:
handles.pole.radius = radius;   % units of pixels
%handles.pole.roi = [1 300 1 180];
handles.pole.roi = [1 300 1 300];   % x,ys

% get header information from the video file:
video.fid = fopen(handles.fname,'r');
video.header = read_mikrotron_datfile_header(video.fid);
handles.nframes = video.header.nframes;
video.width = video.header.width;
video.height = video.header.height;
video.offset = 8192;
video.triggerframe=video.header.triggerframe;
video.triggertick=video.header.triggertick;
video.startframe=video.header.startframe
% set file position to start of first frame
fseek(video.fid,8192,-1);
% keyboard
% figure;
sframes=video.startframe+200;
if sframes>video.header.nframes
    sframes=sframes-video.header.nframes;
end
%tr4=load(strcat(fname(1:(end-3)),'tr4'),'-mat');

%folcenter=tr4.whisker.fp3_all(:,1:2);
frames=[sframes:interval:sframes+1000];
%% loop to find pole
handles.radius2=50;
barPos = zeros(2,length(frames));
i = 0;
for frameidx = frames; %video.header.nframes
    i = i+1;
    
    %frameidx
    if frameidx>video.header.nframes
    frameidx=frameidx-video.header.nframes;
    end
    % Load a frame of data
    handles.frame = load_frame(video,frameidx);
    handles.framemean = mean(handles.frame(:));
    image = double(handles.frame(:,:,1));
    %image = rotate(double(handles.frame(:,:,1)),-60);
    if drawing;
        clf;
        imagesc(image)
        axis image
        colormap gray
      
        %[point]=findpoint(image);
        %imshowpair(image,BW,'montage')
        %     keyboard
        
        % Re-scale figure brightness to 0:0.995 of range
        n = hist(image(:),0:255);
        ncut = find(cumsum(n)/sum(n)>0.995,1);
        caxis([0 ncut])
        clear n ncut
    end
    
    % Find pole centre based on convolution of a pole-sized disk across the
    % image
    % gof = goodness of fit
    [polecentre, gof] = polefit(handles.frame,handles.pole.radius,handles.pole.roi);
    %[polecentre2, gof] = polefit(handles.frame,handles.radius2,[ 1 100 301 470]);
    if drawing
        hold on
        
        plot_circle(polecentre,handles.pole.radius)
        %plot_circle(polecentre2,handles.radius2)
        %plot_circle(folcenter(frameidx,:),10)
        title(sprintf('frame %d: pc = (%d,%d) gof=%.1f', frameidx, polecentre,gof))
        drawnow
        pause(0.01)
    end
    
    %barPos(:,i) = polecentre-folcenter(frameidx,:);
    barPos(:,i) = polecentre;

    %     clf
    
end
fclose('all')
toc
clear video
end

%% Subfunctions
function [point]=findpoint(image)
BW = edge(image);
%imshow(BW)
BW=rot90(BW,3);
point=[0 0];
[rows,cols]=find(BW(40:100,1:100));
scatter(rows,cols)
hold on
x=rows;
s0=[-10 100];
fun= @(p,x)(x*p(1)+p(2));
[p]= lsqcurvefit(fun,s0,x,cols)
ysim=p(1)*x+p(2);
plot(rows,ysim,'r');
hold off
end

function frame = load_frame(videopmtrs,frameidx)

offset = videopmtrs.header.imagesize * (frameidx-1) + videopmtrs.offset;
fseek(videopmtrs.fid,offset,-1);
tmp = fread(videopmtrs.fid,videopmtrs.header.imagesize-24,'uint8=>uint8');
tmp = reshape([tmp; zeros(24,1)],videopmtrs.width,videopmtrs.height)';
frame = uint8(zeros(videopmtrs.height,videopmtrs.width,3));
frame(:,:,1) = tmp;
frame(:,:,2) = tmp;
frame(:,:,3) = tmp;
clear tmp
%fclose('all') ; % Trying to fix 'too many open files' error 
end

%

function plot_circle(centre,radius)

dtheta = 2*pi/100;
theta = 0:dtheta:2*pi-dtheta;

x = centre(1) + radius*cos(theta);
y = centre(2) + radius*sin(theta);

plot(x,y,'r:',centre(1),centre(2),'rx')
% plot(centre(1),centre(2),'rx')
end

%

function [polecentre, gof] = polefit(frame,radius,roi)
% TO DO:
%       find pole of arbitrary size - maybe use a template with a gradiant?
%       otherwise, try a range of sizes and output the best

% ALSO. Take whisker contact into account. Use whisker theta to find
% contact side, also previous pole position to know whether contact is
% expected.

frame = double(frame(:,:,1));
frame = frame - mean(frame(:));                         % zero mean
frame = frame([roi(3):roi(4)],[roi(1):roi(2)]);         % Restrict search to ROI


[X,Y] = meshgrid(-radius-2:radius+2,-radius-2:radius+2);
kernel = ones(2*radius+5);
% kernel(1:2*radius,1:2*radius) = 0;    % pacman
idx = X.^2 + Y.^2 <= radius^2;
alpha = (sum(kernel(:))-sum(idx(:)))/sum(idx(:));   % pacman
% alpha = (sum(kernel(:))-.75*sum(idx(:)))/sum(idx(:));   % pacman
kernel(idx) = -alpha;

% clear idx
z = conv2(frame,kernel,'valid');
[gof, midx] = max(z(:)/sum(idx(:)));
[X,Y] = meshgrid([1:size(z,2)],[1:size(z,1)]);
X = X(:); Y = Y(:);
polecentre(1) = X(midx) + (length(kernel)-1)/2;
polecentre(2) = Y(midx) + (length(kernel)-1)/2;

end