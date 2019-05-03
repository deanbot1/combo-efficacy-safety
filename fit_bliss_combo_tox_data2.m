function [bfit,boot,figs] = fit_bliss_combo_tox_data(Xvar,Yvar,btox,Ptox,Nboot,Xname,Yname,aref)

if nargin < 8
    aref = -1;
end

figs(1)=figure;

vtox = @(v,b)(glmval(b,log(v),'logit')); % monotherapy tox in 1 variable
ctox = @(x,y,bx,by,a)(max(0.0001,min(0.9999,vtox(x,bx)+vtox(y,by)+a*vtox(x,bx).*vtox(y,by)))); % combo tox function
rs = @(x)(reshape(x,[],1));

%only return b from unilogifit function, mtx is not returned
xax=subplot(2,2,4); bx = unilogifit(Xvar(Yvar==0),btox(Yvar==0),'loglin',Ptox); xlabel(Xname);% fit X monotherapy data
yax=subplot(2,2,1); by = unilogifit(Yvar(Xvar==0),btox(Xvar==0),'loglin',Ptox); xlabel(Yname);% fit Y monotherapy data
set(yax,'view',[-90 90],'Xaxislocation','top');

%cost function of logistic regression f=objf
objf = @(X,Y,DLT,BX,BY,alpha)(-sum(DLT.*log(ctox(X,Y,BX,BY,alpha)) + (1-DLT).*log(1-ctox(X,Y,BX,BY,alpha)))); % LLF(log likelihood function)
objfun = @(alpha)(0*(alpha-aref).^2 + objf(Xvar,Yvar,btox,bx,by,alpha));

[abest,fbest] = fminsearch(objfun,aref);%,optimset('Display','iter'));

% calculate AICC of bliss model, number of parameters in bliss is 5 bx,by,a
aicc_bliss = aicc(fbest, length(Xvar),5)


[XX,YY] = meshgrid([0:.01:1]*max(Xvar),[0:.01:1]*max(Yvar));
PP = reshape(ctox(rs(XX),rs(YY),bx,by,abest),size(XX));
xyax = subplot(2,2,2);

contourf(XX,YY,PP); hold on; colormap pink; caxis([0 1]); shading flat;
%  surfc(XX,YY,PP)
plot(Xvar(btox==0),Yvar(btox==0),'go','MarkerFaceColor','g');
plot(Xvar(btox==1),Yvar(btox==1),'ro','MarkerFaceColor','r');
PPref = reshape(ctox(rs(XX),rs(YY),bx,by,aref),size(XX));

%plot the MTEtrue line / original param [-1;5] [-10;4]
% PPtrue= reshape(ctox(rs(XX),rs(YY),[-1;5],[1;8],10)  ,size(XX));
% [C,h] = contour(XX,YY,PPtrue,[Ptox Ptox]);
% set(h,'Color','y','LineWidth',3);

[C,h] = contour(XX,YY,PPref,[Ptox Ptox]);
set(h,'Color','c','LineWidth',2);
[C,h] = contour(XX,YY,PP,[Ptox Ptox]);
set(h,'Color','b','LineWidth',2);
% [C,h] = contour(XX,YY,PPtrue, [Ptox Ptox]);
% set(h, 'Color', 'g','LineWidth',2);
MTE_curve = C;
xlabel(Xname);
ylabel(Yname);

title(['\alpha_{best}=' num2str(abest) ',\alpha_{ref}=' num2str(aref)]);

set(gca,'Xtick',get(xax,'Xtick'));
set(gca,'Ytick',get(yax,'Xtick'));
axes(xax); plot(XX(1,:),PP(1,:),'m-');
axes(yax); plot(YY(:,1),PP(:,1),'m-');

subplot(2,2,3); 

hwait = waitbar(0,'bootstrapping...');
for j = 1:Nboot
	jboot = ceil(length(Xvar)*rand(size(Xvar))); % sample with replacement
	Xboot = Xvar(jboot); Yboot = Yvar(jboot); Tboot = btox(jboot);  %get a subpop and randomize the order of Xvar
	bxboot(:,j) = glmfit(log(Xboot(Yboot==0)),Tboot(Yboot==0),'binomial');
	byboot(:,j) = glmfit(log(Yboot(Xboot==0)),Tboot(Xboot==0),'binomial');
	objboot = @(alpha)(objf(Xboot,Yboot,Tboot,bxboot(:,j),byboot(:,j),alpha));
	[aboot(j),fval(j),exitflag] = fminsearch(objboot,abest);
    if isnan(fval(j)) == true
        warning(['throwing out aboot = ' num2str(aboot(:,j)')]);
        %keyboard
%         %keyboardrowing out
        aboot(j) = NaN; 
    end
	waitbar(j/Nboot,hwait);
end
close(hwait);
if Nboot > 0
	Nbad = length(find(isnan(aboot))); % number of 'bad' bootstraps
end
if Nboot>0 && Nbad/Nboot < 0.1
	xlabel('\alpha estimates');
	%plot(aref*[1 1],Ylim,'r-');text(aref,max(Ylim),'\alpha_{true}','Color','r');
	pvec = [0.05 0.95];
	aquant = quantile(aboot,pvec);
    histogram(aboot,[aquant(1):(aquant(2)-aquant(1))/20:aquant(2)]); hold on;
%     hist(aboot); hold on;
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
        PPP(:,j) = ctox(rs(XX),rs(YY),bxboot(:,j),byboot(:,j),aboot(j));
    end
    PPP = PPP(:,~isnan(aboot));
    PPPP = quantile(PPP,pvec,2);
    for k = 1:length(pvec) % simply change and use 0.05 and 0.95 quantile 
        % this contour C are the 'outline' of the 90% confidence 
        % region for MTEbest
        [C,h] = contour(XX,YY,reshape(PPPP(:,k),size(XX)),[Ptox Ptox]);
        set(h,'Color','y','LineWidth',3);
        if k ==1
            CI_region_5= C(:,2:length(C));
%             CI_region_5= C(:,find(C(1,:)==0):find(C(2,:)==0));
        elseif k==2
            CI_region_95= C(:,2:length(C));
%             CI_region_95= C(:,find(C(1,:)==0):find(C(2,:)==0));
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
        PPPP=[];
	end
end

bfit = struct('abest',abest,'bx',bx,'by',by);
if Nboot > 0
	boot = struct('abest',aboot,'bx',bxboot,'by',byboot,'CI_5',CI_region_5,'CI_95', CI_region_95,...
        'PPPP',PPPP,'XX',XX,'YY',YY,'PP',PP,'MTE_curve',MTE_curve);
else
	boot = struct('MTE_curve',MTE_curve);
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