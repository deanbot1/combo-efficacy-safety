clear all
close all
%%
rng('default')

load Loewe_MTE_allPatients.mat
load Preclinical_surfacefits
%col_id=['r','b','m'];
%col_id=hsv(length(Cellline_name));
%col_id=lines(length(Cellline_name));
col_id=varycolor(length(Cellline_name));

X=C(:,2:end)';
X=X/1000;
B= X(:,2);
A= X(:,1);
AD= A./0.004560;%Covert exposure to dose
BD=B./0.0007291; %Covert exposure to dose
A1=A./max(A);
B1=B./max(B);
[theta,rho]=cart2pol(A1,B1);
ax1=subplot(40,1,1:32)
hold on
set(gca, 'ColorOrder', col_id);
%theta2 = datasample(theta,length(Cellline_name),'Replace',false)
theta2=repmat([max(theta), min(theta)],1, ceil(length(Cellline_name)/2));
allign1=repmat({ 'left','right'}, 1 , ceil(length(Cellline_name)/2));
allign2=repmat({'bottom', 'top'}, 1 , ceil(length(Cellline_name)/2));
%Cellline_name={'Xenograft-1','Xenograft-2','Xenograft-3','Xenograft-4','Xenograft-5','Xenograft-6','Xenograft-7','Xenograft-8','Xenograft-9','Xenograft-10','Xenograft-11'}

%for i=1:length(Cellline_name)
for i=2  % for MDA-DB cell line only!!
beta=Beta_est{i};
TGI=simTGI(Model_strid{i},[A,B],beta)
ax1(i)=plot(theta,TGI,'LineWidth',2,'color',col_id(i,:))
%ax1=plot(theta,TGI,'LineWidth',2,'color','b')
igood=find(TGI==max(TGI));
foo=get(ax1(i),'color')
ax11(i)=plot(theta(igood),TGI(igood),'o','MarkerSize',8,'MarkerFaceColor',foo,'MarkerEdgeColor','k')
ax2=text(theta2(i),TGI(find(theta==theta2(i))), Cellline_name{i},'Color',foo,'FontSize',11,'horizontalAlignment', allign1{i},'VerticalAlignment','middle')
set(ax2, 'rotation', 0,'FontWeight','Bold')
legendInfo{i} = Cellline_name{i}; % or whatever is appropriate
end
ylabel('%GRI')
set(gca, 'XTickLabel',[])
box on
set(gca,'units','normalized','xlim',[0 max(theta)],'color','none','XTick',[min(theta):max(theta)/6:max(theta)],'ylim',[0 140])
%set(gca,'FontName','Arial','FontWeight','Bold','FontSize',12)
set(gca,'FontName','Arial','FontSize',16)
%legend(ax1,legendInfo, 'Location','northeast','Orientation','vertical')
%columnlegend(3, legendInfo, 'Location', 'NorthWest', 'boxon'); 

%% now add suplo axis
ax1 = ax1(end);
ax2=subplot(40,1,33:34);
copyobj(allchild(ax1),ax2);
L=get(gca,'children');
for i=1:length(L)
       set(L(i),'visible','off');
       % Use the following to check mark all the lines : set(L(i),'visible','on');
end

%patch([min(theta),max(theta),max(theta),min(theta)],[0,0,1,0],[1 0 0 1])
%colormap(pink)
set(gca,'units','normalized','xlim',[0 max(theta)],'color','none','XTick',[min(theta):max(theta)/6:max(theta)],'Ylim',[0 1],'YTick',[],'YColor','none','XColor','r','Ticklength',[0.001 0.001])
set(get(gca,'YLabel'),'visible','off')
% set(gca,'xticklabel',num2str(get(gca,'xtick')','%.3f'))
h=get(gca,'xticklabel');
for i=1:length(h);
    hh(i)=str2num(h{i});
end
p=polyfitZero(theta,BD,5);
%BDpred=polyval(p,theta);
%plot(theta,BD,'o',theta,BDpred,'-r')
hBDpred=polyval(p,hh)*4.5;
patch([min(hh) hh max(hh)],[0 hBDpred/max(hBDpred) 0],[0 hBDpred/max(hBDpred) 1])
set(gca,'xticklabel',num2str(hBDpred,'%0.0f\n'))
colormap(flipud(pink))
%xlabel('TD-1 Dose(mg) for QW','Color','r')
xlabel('TAK-228 Dose(mg) for QDX3','Color','r')
set(gca,'FontName','Arial','FontSize',16)
box off
%% last subplot

ax3=subplot(40,1,39:40);
copyobj(allchild(ax1),ax3);
L=get(gca,'children');
for i=1:length(L)
       set(L(i),'visible','off');
       % Use the following to check mark all the lines : set(L(i),'visible','on');
end

% patch([min(theta),max(theta),max(theta),min(theta)],[0,0,1,0],[1 0 0 1])
% colormap(pink)
set(gca,'units','normalized','xlim',[0 max(theta)],'color','none','XTick',[min(theta):max(theta)/6:max(theta)],'Ylim',[0 1],'YTick',[],'YColor','none','XColor','b','Ticklength',[0.001 0.001])
%xlabel('TD-2 Dose(mg) for QWX3','color','b')
xlabel('TAK-117 Dose(mg) for QDX3','color','b')
set(ax3,'xdir','reverse','Ticklength',[0.001 0.001])
set(get(gca,'YLabel'),'visible','off')
h=get(gca,'xticklabel');
hh=[];
for i=1:length(h);
    hh(i)=str2num(h{i});
end
p1=polyfitZero(theta,fliplr(AD')',5);
%ADpred=polyval(p1,theta);
%plot(theta,fliplr(AD')','o',theta,ADpred,'-r')
hADpred=polyval(p1,hh)*.69;
hADpred=fliplr(hADpred')';
set(gca,'xticklabel',num2str(hADpred,'%0.0f\n'))
patch([min(hh) hh max(hh)],[0 hADpred/max(hADpred) 0],[0 hADpred/max(hADpred) 1])
colormap(flipud(pink))
set(gca,'FontName','Arial','FontSize',16)
box off

%%
set(gcf,'units','inches',...
            'paperposition',[1.5 1 8 6],'position',[7.0625    4.0208    5.9792    5.6458]);

set(findall(gca,'type','all'),...
        'linewidth',2,...
        'fontsize',20,...
        'fontname','arial',...
        'box','off',...
        'tickdir','out');

print -dtiff -r300 test.tif