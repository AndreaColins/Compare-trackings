function selected=selector
videosselected=[];
positionselected=[];
polenameselected=[];
radiusselected=[];
xlsfile='preselectedair.xlsx';
condition=char('Smooth pole','Open coil','Closed coil','Bamboo','Black Sandpaper','Carbon Pole','Cardboard','Toothpick','Wood','air');
[xls_info,txt] = xlsread(xlsfile);
totaln=12;
%%%%%seleccion del pole a estudiar
for pole=10:10
  
    xls_idx=find(xls_info(:,4)==pole);
    xls_pole=xls_info(xls_idx,:);
    txt_pole=txt(xls_idx,1:2);
    if size(xls_idx)<=totaln
        size(xls_idx)
        condition(pole,:)
        selected=txt_pole;
        selectedp=[xls_pole(:,7),xls_pole(:,8)];
        videosselected=[videosselected;selected];
        positionselected=[positionselected;selectedp];
        
        radius=xls_info(xls_idx(1),6)*ones(size(xls_idx));
        polename(1:size(xls_idx,1),1)={condition(pole,:)}
        polenameselected=[polenameselected;polename];
        radiusselected=[radiusselected;radius];
        
        
    else
        radius=xls_info(xls_idx(1),6)*ones(totaln,1);
        polename(1:totaln,1)={condition(pole,:)}
        polenameselected=[polenameselected;polename];
        radiusselected=[radiusselected;radius];
        days=xls_pole(end,5);
        nbyday=totaln/days;
        selected(1:totaln/days,2)={''};
        selectedpos=zeros(floor(totaln/days),2);
        %%%%formar los 4 grupos diferenciados por posicion del pole
        x=xls_pole(:,7);
        y=xls_pole(:,8);
        
        
        %%%%%form groups
        for i=1:days
            figure
            iselected=0;
            I=find(xls_pole(:,5)==i);
            position=[x(I),y(I)];
            Igroup=zeros(size(position,1),4);
            txt_position=txt_pole(I,:);
            scatter(x(I),y(I))
            hold on
            ylimit=(max(y(I))-min(y(I)))/2+min(y(I));
            xlimit=(max(x(I))-min(x(I)))/2+min(x(I));
            plot(xlimit*ones(1,2),[min(y(I)),max(y(I))],'g');
            plot([min(x(I)),max(x(I))],ylimit*ones(1,2),'r');
            
            pause(1)
            Igroup(:,1)=(position(:,1)<=xlimit)&(position(:,2)<=ylimit);
            Igroup(:,2)=(position(:,1)>xlimit)&(position(:,2)<=ylimit);
            Igroup(:,3)=(position(:,1)<=xlimit)&(position(:,2)>ylimit);
            Igroup(:,4)=(position(:,1)>xlimit)&(position(:,2)>ylimit);
            
            emptygroups=isempty(find(Igroup(:,1),1))+isempty(find(Igroup(:,2),1))+isempty(find(Igroup(:,3),1))+isempty(find(Igroup(:,4),1));
            ngroups=4-emptygroups;
            nbygroup=floor((totaln/days)/ngroups);
            for k=1:4
                I1=find(Igroup(:,k));
                
                for j=1:nbygroup
                    if ~isempty(I1)
                        idx=round(rand(1)*(size(I1,1)-1)+1);
                        iselected=iselected+1;
                        fileidx(iselected)=I1(idx);
                        selected(iselected,:)=txt_position(fileidx(iselected),:);
                        selectedp(iselected,:)=position(fileidx(iselected),:);
                        scatter(position(fileidx(iselected),1),position(fileidx(iselected),2),'r')
                        pause(2)
                        I1(idx)=[];
                        
                    end
                end
                
            end
            txt_position(fileidx,:)=[];
            position(fileidx,:)=[];
            
            if size(selectedp,1)<nbyday
                newselected=nbyday-size(selectedp,1);
                for j=1:newselected
                    idx=round(rand(1)*(size(txt_position,1)-1)+1);
                    iselected=iselected+1;
                    selected(iselected,:)=txt_position(idx,:);
                    selectedp(iselected,:)=position(idx,:);
                    scatter(position(idx,1),position(idx,2),'r')
                    pause(2)
                    txt_position(idx,:)=[];
                    position(idx,:)=[];
                end
            end
            title(strcat('Smooth pole, Day ',int2str(i)))
            xlabel('x position [px]')
            ylabel('y position [px]')
            hold off
            fileidx=[];
            videosselected=[videosselected;selected];
            positionselected=[positionselected;selectedp];
            selectedp=[];
            
            
        end
    end
    clear selected selectedp
    polename={};
    radius=[];
end
A={'Horizontal view','Vertical view','Pole','radius','X position','Y position'};

xlswrite('Selectedvideos.xlsx',A,2,'A1')
xlswrite('Selectedvideos.xlsx',videosselected,2,'A2')
xlswrite('Selectedvideos.xlsx',polenameselected,2,'C2')
xlswrite('Selectedvideos.xlsx',radiusselected,2,'D2')
xlswrite('Selectedvideos.xlsx',positionselected,2,'E2')
end