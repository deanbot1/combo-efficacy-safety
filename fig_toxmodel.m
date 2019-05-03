%% supplementary figure showing level curves of bliss-like-model that we are using

clear all; close all

logit = @(p)log(p./(1-p));
pcrit = 0.25; % probablity at which we have maximum tolerability of DLTs
b = logit(0.01); % spontaneous frequency of DLT in absence of drugs must be < logit(0.25) but > logit(0)=-Inf
mx =(-log(3)-b)/3; % logistic regression slope of drug x
my =(-log(3)-b)/2; % logistic regression slope of drug y
Cx = 3*[0:.01:1.2];
Cy = 2*[0:.01:1.2];


prob = @(x,y,alpha)1./(1+exp(-(b + mx*x + my*y + alpha*x.*y)));
[XX,YY] = meshgrid(Cx,Cy);

alist =  (mx*my/(b-logit(pcrit)))*[1 0.5 0 -10] ; 
acolors = {'b','m','r','g'};
figure;
astrings = {'\alpha = \alpha*','\alpha* < \alpha < 0','\alpha = 0','\alpha > 0'};
%astrings = {'$$ \alpha = {m_xm_y}/{log(p/(1-p))+b} $$','\alpha < 0','\alpha = 0','\alpha > 0'}
for j = 1:length(alist)
    alpha = alist(j);
   
%     [C,h] = contourf(XX,YY,prob(XX,YY,alpha),[0:.05:1]); hold on; colorbar; 
%     clabel(C,h);
%     caxis([0 1]);
%     title(num2str(alpha));
    C = contourc(Cx,Cy,prob(XX,YY,alpha),pcrit*[1 1]);
    kk = find(C(1,:)==pcrit);
    for k = kk
        
    h=plot(C(1,k+1:k+C(2,k)),C(2,k+1:k+C(2,k)),[acolors{j} '-'],'LineWidth',2); hold on
    NC = size(C,2)-1;
    if j>1
    text(C(1,round(NC/2)),C(2,round(NC/2)),astrings{j},'Color',acolors{j},'BackgroundColor','w','Rotation',-45);%,'horizontalalignment','center');
    else
    text(2.4,1.8,astrings{j},'Color',acolors{j},'BackgroundColor','w','Rotation',0);%,'horizontalalignment','center');
    end
    end
end
%legend(h,num2str(alist'));
set(gca,'XLim',[0 max(Cx)],'Ylim',[0 max(Cy)]);
set(gca,'Xtick',[0 3],'XtickLabel',{'0','MTE_X'});
xlabel('C_X');
set(gca,'Ytick',[0 2],'YtickLabel',{'0','MTE_Y'});
ylabel('C_Y');
