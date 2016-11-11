function [idx,result]=comparetr4(files,nm)
%%%%%%%%%%%para cambiar la posicion de inicio del video
txt=importdata('listofvideos.txt');
textdata=string(char(txt.textdata));

%%%%%sacar las comas
textdata=textdata(:,2:end-1)
%%%%%
for f=1:size(files,1)   
v(f)=load(files(f,:),'-mat');

%encontrar la primera frame y gardarla en la estructura v

%dar vuelta los vectores desde la primera frame a la ultima
 % idx2=circshift(idx,size(idx,1)-315,2)
end
size(files,1)  

result=zeros(length(v(1).whisker(1).tracked),4,size(files,1));
if nm==1
f1=figure;
f2=figure;
end
colours={'b.-'  'g.-'  'r.-'  'c.-'  'y.-'  'm.-'};
spr = 4;
for f=1:length(v)
for w = 1:length(v(f).whisker)
    wmod=f;
    %wmod = rem(w,length(v(f).whisker)) + length(v(f).whisker)*(w==length(v(f).whisker)); % to avoid running out of colours
    idx = find(v(f).whisker(w).tracked);
    
    theta = zeros(3,length(v(f).whisker(w).tracked));
    azimuth = zeros(1,length(v(f).whisker(w).tracked));
    elevation = zeros(1,length(v(f).whisker(w).tracked));
    twist = zeros(1,length(v(f).whisker(w).tracked));
    curv_hc = zeros(3,length(v(f).whisker(w).tracked));  % head-centred
    curv_fc = zeros(3,length(v(f).whisker(w).tracked));  % follicle-centred - later
    curv3 = zeros(1,length(v(f).whisker(w).tracked));
    np = zeros(2,length(v(f).whisker(w).tracked));
    r3 = squeeze(v(f).whisker(w).r3all);
    
    for i = 1:length(idx)
        fr = idx(i);
        r = squeeze(v(f).whisker(w).r3all(fr,:,:));
        theta(:,fr) = base_angle3(r,0);
        curv_hc(1,fr) = curvature(r([2 3],:),0);
        curv_hc(2,fr) = curvature(r([1 3],:),0);
        curv_hc(3,fr) = curvature(r([1 2],:),0);
        curv3(fr) = curvature3(r,0);
        tv = bezierdtval(r,0);  % tangent to whisker at base
        tv = tv/norm(tv);
        cv = bezierdt2val(r,0); % 2nd derivative to whisker at base
        NPO = eye(3)-tv*tv'; % projection operator, onto plane normal to 'tv'
        %         xpr = NPO*[1 0 0]'; % projection of x axis
        %         zpr = NPO*[0 0 1]'; % projection of z axis
        np(1,fr) = [1 0 0]*NPO*cv;   % projection of cv onto NP, projected onto x axis
        np(2,fr) = [0 0 1]*NPO*cv;   % projection of cv onto NP, projected onto z axis
        %         curv(1,fr) = cv'*xpr; % projection of 2nd deriv vector
        %         curv(2,fr) = cv'*zpr;
    end
    clear i fr tv cv
    azimuth(idx) = theta(3,idx);
    elevation(idx) = theta(1,idx);
    twist(idx) = atan2(-np(2,idx),np(1,idx))*(180/pi);
    
  

    if nm==1
    figure(f1)
    h1(1) = subplot(spr,2,1); hold on
   % idx2=circshift(idx,size(idx,1)-315,2)
    plot(idx,azimuth(idx),colours{wmod}), ylabel('azimuth')
    title('azimuth: 90deg is normal to ant-post axis; >90deg for whisker protracted')
    h1(2) = subplot(spr,2,3); hold on
    plot(idx,elevation(idx),colours{wmod}), ylabel('elevation')
   
    title('elevation: 90deg is horizontal; >90deg for whisker tip oriented up')
    h1(3) = subplot(spr,2,5); hold on
    plot(idx,twist(idx),colours{wmod}), ylabel('twist')
    title('twist: 90deg is concave down; <90deg for concave posterior')
    
    h2(1) = subplot(spr,2,2); hold on
    plot(idx,curv_hc(1,idx),colours{wmod}), ylabel('kappa coronal')
    h2(2) = subplot(spr,2,4); hold on
    plot(idx,curv_hc(3,idx),colours{wmod}), ylabel('kappa horizontal')
    h2(3) = subplot(spr,2,6); hold on
    plot(idx,curv3(idx),colours{wmod}), ylabel('kappa3')
    
if f==length(v)
    legendCell =cellstr(num2str([1:length(v)]', 'N=%-d'));
    l1=legend(legendCell,'Location', 'BestOutside');
    set(l1,'position',[0.45 0.1 0.1 0.1])
end
a = axis;
    tmp = curv_hc(:,idx);
    tmp = [tmp(:)' curv3(idx)];
%     axis([a(1:2) min(tmp) max(tmp)])
    clear tmp
    linkaxes(h2,'x')
    linkaxes(h1,'x')
   
    figure(f2)
    subplot(2,2,1);    hold on
    plot(azimuth(idx),elevation(idx),colours{wmod})
    xlabel('azimuth'), ylabel('elevation')
    title('tangent direction')
    axis square
    %     linkaxes(h,'x')
    %     subplot(2,1,2); hold on
    % %     plot(curv(1,idx),curv(2,idx),handles.colours{w})
    %     compass(curv(1,idx),curv(2,idx),handles.lines{wmod})
    %     xlabel('horizontal projection of normal plane'), ylabel('vertical projection of normal plane')
    %     axis square
    %     a = axis;
    %     xlim([-max(abs(a(1:2))),max(abs(a(1:2)))])
    %     ylim([-max(abs(a(3:4))),max(abs(a(3:4)))])
    
    subplot(2,2,2), hold on
%     plot(np(1,idx),np(2,idx),handles.points{wmod})
%     xlabel('x'), ylabel('z')
%     a = axis;
%     xlim(max(abs(a(1:2)))*[-1 1])
%     ylim(max(abs(a(3:4)))*[-1 1])
%     title('curvature direction')
    plot(azimuth(idx),twist(idx),colours{wmod})
    xlabel('azimuth'), ylabel('twist')
    title('tangent direction')
    axis square
    
    subplot(2,2,3), hold on
%     plot(frames(idx),sqrt(sum(np(:,idx).^2+np(:,idx).^2,1)),handles.lines{wmod})
%     ylabel('|B''''(0)| (pixels)'), title('Intrinsic curvature index')
    
%     subplot(2,2,4), hold on
%     plot(frames(idx),twist(idx),handles.lines{wmod})
    plot(twist(idx),elevation(idx),colours{wmod})
    xlabel('twist'), ylabel('elevation')
    title('tangent direction')
    axis square
    end
end
result(:,1,f)=azimuth;
result(:,2,f)=elevation;
result(:,3,f)=curv_hc(1,idx);
result(:,4,f)=curv_hc(2,idx);
end
idx=idx';
clear w

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

function theta = base_angle3(r,t)
% angle of tangent vector to Bezier defined by control points {r} at point
% t
% theta(1) is in the x-y plane
% theta(2) is in the y-z plane
% theta(3) is in the x-z plane (not sure meaningful)

dBdt = bezierdtval(r,t);    % tangent vector
theta = zeros(3,length(t));
theta(3) = atan2(-dBdt(2,:),dBdt(1,:))*(180/pi);        % azimuth - whisker pointing caudal is 0'.
theta(1) = 180-atan2(-dBdt(2,:),dBdt(3,:))*(180/pi);    % elevation - whisker pointing ventral is 0'.
theta(2) = atan2(-dBdt(3,:),dBdt(1,:))*(180/pi);
clear dBdt

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

function kappa = curvature3(r,t)
% extension of the kappa formula to 3 dimensions (wikipedia curvature page)
dBdt = bezierdtval(r,t);
d2Bdt2 = bezierdt2val(r,t);

if numel(t)>1
    error('Generalise the code!')
end

kappa = norm(cross(dBdt,d2Bdt2))/norm(dBdt).^3;

