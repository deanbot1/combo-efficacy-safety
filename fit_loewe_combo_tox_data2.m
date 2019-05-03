function [bfit,boot,figs] = fit_loewe_combo_tox_data(Xvar,Yvar,btox,Ptox,Nboot,Xname,Yname,aref)

if nargin < 8
    aref = -1;
end

figs(1) = figure;

% vtox = @(v,b)(glmval(b,log(v),'logit')); % monotherapy tox in 1 variable
% ctox = @(x,y,bx,by,a)(max(0.0001,min(0.9999,vtox(x,bx)+vtox(y,by)+a*vtox(x,bx).*vtox(y,by)))); % combo tox function

rs = @(x)(reshape(x,[],1));
xax=subplot(2,2,4);
if ~isempty(find(Yvar==0))
 [bfit.bx,mtex] = unilogifit(Xvar(Yvar==0),btox(Yvar==0),'linear',Ptox); xlabel(Xname);% fit X monotherapy data
end
yax=subplot(2,2,1);
if ~isempty(find(Xvar==0))
 [bfit.by,mtey] = unilogifit(Yvar(Xvar==0),btox(Xvar==0),'linear',Ptox); xlabel(Yname);% fit Y monotherapy data
set(yax,'view',[-90 90],'Xaxislocation','top');
end

% objf = @(X,Y,DLT,BX,BY,alpha)(-sum(DLT.*log(ctox(X,Y,BX,BY,alpha)) + (1-DLT).*log(1-ctox(X,Y,BX,BY,alpha)))); % LLF
% objfun = @(alpha)(objf(Xvar,Yvar,btox,bx,by,alpha));
% 
% [abest,fbest] = fminsearch(objfun,aref);%,optimset('Display','iter'));

xf = @(x)(x);

mtex=1; mtey=1;
X = Xvar/mtex; Y = Yvar/mtey;
[bbest, dev]= glmfit(xf([X Y X.*Y]),btox,'binomial');
abest = bbest(4);

[XX,YY] = meshgrid([0:.01:1]*max(Xvar),[0:.01:1]*max(Yvar));
PP = reshape(glmval(bbest,xf([rs(XX/mtex) rs(YY/mtey) rs((XX/mtex).*(YY/mtey))]),'logit'),size(XX));
xyax = subplot(2,2,2);



% calculate AICC of loewe model; nLL value is dev/2; 
 aicc_loewe = aicc(dev/2, length(Xvar), length(bbest))

contourf(XX,YY,PP); hold on; colormap pink; caxis([0 1]); shading flat;
plot(Xvar(btox==0),Yvar(btox==0),'go','MarkerFaceColor','g');
plot(Xvar(btox==1),Yvar(btox==1),'ro','MarkerFaceColor','r');
bref = bbest; bref(4)=aref;
PPref = reshape(glmval(bref,xf([rs(XX/mtex) rs(YY/mtey) rs((XX/mtex).*(YY/mtey))]),'logit'),size(XX));
[C,h] = contour(XX,YY,PPref,[Ptox Ptox]);
set(h,'Color','c','LineWidth',2);
[C,h] = contour(XX,YY,PP,[Ptox Ptox]);
set(h,'Color','b','LineWidth',3);
xlabel(Xname);
ylabel(Yname);
MTEcurve = C;

title(['\alpha_{best}=' num2str(abest) ',\alpha_{ref}=' num2str(aref)]);

set(gca,'Xtick',get(xax,'Xtick'));
set(gca,'Ytick',get(yax,'Xtick'));
axes(xax); plot(XX(1,:),PP(1,:),'m-');
axes(yax); plot(YY(:,1),PP(:,1),'m-');

subplot(2,2,3); 

hwait = waitbar(0,'bootstrapping...');
for j = 1:Nboot
	jboot = ceil(length(Xvar)*rand(size(Xvar))); % sample with replacement
	Xboot = Xvar(jboot)/mtex; Yboot = Yvar(jboot)/mtey; Tboot = btox(jboot);
	bboot(:,j) = glmfit(xf([Xboot Yboot Xboot.*Yboot]),Tboot,'binomial');
	waitbar(j/Nboot);
end
close(hwait);
if Nboot > 0
	Nbad = length(find(isnan(bboot(4,:)))); % number of 'bad' bootstraps
end
if Nboot>0 & Nbad/Nboot < 0.1
	aboot = bboot(4,:);
	xlabel('\alpha estimates');
	%plot(aref*[1 1],Ylim,'r-');text(aref,max(Ylim),'\alpha_{true}','Color','r');
	pvec = [0.05 .5 .95];
	aquant = quantile(aboot,pvec);
    %hist(aboot,[aquant(1):(aquant(3)-aquant(1))/20:aquant(3)]); hold on;
    hist(aboot); hold on;
    Ylim = get(gca,'Ylim');
	vline = @(x,color)(plot(x*[1 1],Ylim,[color '-']));
	vlabel = @(x,label,color)(text(x,max(Ylim),label,'Color',color,'VerticalAlignment','bottom'));
	vline(aref,'r');vlabel(aref,'\alpha_{ref}','r');
	vline(abest,'m');vlabel(abest,'\alpha_{best}','m');    
	for k = 1:length(aquant)
		vline(aquant(k),'g');
		vlabel(aquant(k),['\alpha_{' num2str(pvec(k)) '}'],'g');
	end
	%set(gca,'Xlim',[min(union(aquant,aref)) max(union(aquant,aref))]);
    colorbar; caxis([0 1]); 
    axes(xyax);
    for j = 1:Nboot
        PPP(:,j) = glmval(bboot(:,j),xf([rs(XX/mtex) rs(YY/mtey) rs((XX/mtex).*(YY/mtey))]),'logit');
    end
    PPP = PPP(:,~isnan(aboot));
    PPPP = quantile(PPP,pvec,2);
    for k = [1 3]
        [C,h] = contour(XX,YY,reshape(PPPP(:,k),size(XX)),[Ptox Ptox]);
        set(h,'Color','b','LineWidth',1);
        if k ==1
            CI_region_5= C(:,2:length(C));
        elseif k==3
            CI_region_95=C(:,2:length(C));
        end
    end
    
    if isempty(CI_region_5) ==true
        CI_region_5 = [];
    elseif (max(CI_region_5(1,:))-max(Xvar))>0.001 ||...
            (max(CI_region_5(2,:))-max(Yvar))>0.001 ||...
            (abs(max(Yvar)-max(CI_region_5(2,length(CI_region_5))))<0.001) ||...
            (abs(max(Xvar)-CI_region_5(1,1))<0.001) 
            
        CI_region_5 = [];
    end
    
    if (max(CI_region_95(1,:))-max(Xvar))>0.001 || (max(CI_region_95(2,:))-max(Yvar))>0.001
        CI_region_95 = [];
    end

    
else
	if Nboot > 0 % then it must have been too many bad boots
		warning(sprintf('TOO MANY (%d/%d) BAD BOOTSTRAP RUNS; BOOTSTRAP RESULTS WILL NOT BE SHOWN',Nbad,Nboot));
        CI_region_5 =[];
        CI_region_95=[];
        PPPP = [];
  	end
end


bfit.abest = abest;
bfit.bbest = bbest;
%bfit = struct('abest',abest,'bx',bx,'by',by,'bbest', bbest);
if Nboot > 0
	boot = struct ('bboot',bboot,'CI_5',CI_region_5,'CI_95', CI_region_95, ...
        'PPPP',PPPP,'XX',XX,'YY',YY,'PP',PP,'MTE_curve',MTEcurve);
else
	boot = struct('MTE_curve',MTEcurve);
end

newfig = figure;
newh = copyobj(get(xyax,'Children'),gca);
axis tight; colormap pink; caxis([0 1]); xlabel(Xname); ylabel(Yname); colorbar;
title(['\alpha_{best} = ' num2str(abest) ', \alpha_{ref} = ' num2str(aref)]);
figs(2) = newfig;

function aicc = aicc(nLL,n,k)
% corrected Aikake Information Criterion
% nLL = negative log likelihood 
% n = number of data points
% k = number of parameters

aic = 2*k + 2*nLL;
aicc = aic + 2*k*(k+1)/(n-k-1);
