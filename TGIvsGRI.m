%% makes a graph for supplementary materials showing TGI/GRI relationship
gvehMDA = 0.0693;
tf = 21;
tgi = @(GRI,psi)100*((exp(psi)./(exp(psi)-1))*(1-exp(-psi.*GRI/100)));
figure; hold on; 
        text(30,100,'~70%RR','Horizontalalignment','right');
        plot([30 220],[100 100],'k-');
        plot([30 220],[60 60],'k-');
        text(30,60,'~00%RR','Horizontalalignment','right');

GRI = [0:1:200]; 
for gveh = union(10.^[-3:1:-1],gvehMDA)
    psi = tf*gveh;
    TGI = tgi(GRI,psi); 
    h=plot(GRI,TGI,'-','Linewidth',2); 
    if abs(gveh-gvehMDA)>eps
        text(GRI(end),TGI(end),['g_{veh}=' num2str(gveh)],'color',get(h,'color'));
    else
        set(h,'LineWidth',3)
        text(GRI(end),TGI(end),['g_{veh}^{361}=' num2str(gveh)],'color',get(h,'color'),'fontweight','bold');
    end
end
grid on;
xlabel('GRI');ylabel('TGI')
set(gca,'XLim',[0 220]);
title(sprintf('TGI/GRI relationship (%d day follow up)',tf));

