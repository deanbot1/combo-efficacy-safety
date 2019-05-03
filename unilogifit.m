function [b,mtx]=unilogifit(X,Y,scale,Ptox)

switch scale
	case 'linear'
		xf = @(x)(x);
        mtxf = @(p,b1,b2)((log(p/(1-p))-b1)/b2);
	case 'log'
		xf = @(x)(log(x));
        mtxf = @(p,b1,b2)(exp((log(p/(1-p))-b1)/b2));
	case 'loglin'
		xf = @(x)(log(x));
        mtxf = @(p,b1,b2)(exp((log(p/(1-p))-b1)/b2));
		scale = 'linear';
end


violin(X(Y==0),0,'go','markerfacecolor','g','markersize',3); hold on;
violin(X(Y==1),1,'ro','markerfacecolor','r','markersize',3); hold on;

ntiles = 3; % number of X-tiles to plot little whiskers for
qtiles = [0:1/ntiles:1];
cutoffs = quantile(X,qtiles);
for j = 1:ntiles
	igood = find(X >= cutoffs(j) & X < cutoffs(j+1));
	xmean = exp(mean(log(X(igood)))); % geometric mean
	[phat,pci] = binofit(sum(Y(igood)),length(igood));
	plot(xmean,phat,'mo','MarkerFaceColor','m');
	plot(xmean*[1 1],pci,'m-','LineWidth',2);
end

%plot(X(Y==0),Y(Y==0),'go','MarkerFaceColor','g'); hold on;
%plot(X(Y==1),Y(Y==1),'ro','MarkerFaceColor','r'); hold on;
[b,dev,stats] = glmfit(xf(X),Y,'binomial');
xx = [min(X):(max(X)-min(X))/100:max(X)];

%compute the 95% confidence bounds for the predicted values 
[yhat,dylo,dyhi] = glmval(b,xf(xx),'logit',stats);
plot(xx,yhat,'k-','linewidth',3);
plot(xx,yhat-dylo,'k-');plot(xx,yhat+dyhi,'k-');
mtx = mtxf(Ptox, b(1), b(2)); % max tolerated exposure
plot(mtx*[1 1],[0 Ptox],'b-','LineWidth',2);
plot(mtx,0,'bv','MarkerFaceColor','b');
set(gca,'Xscale',scale)
ylabel('p(DLT)');
axis tight
%title(scale);
plot(xx([1 end]),Ptox*[1 1],'k:');
%text(xx(1),Ptox,'Ptox ','HorizontalAlignment', 'right');
set(gca,'Xlim',[min(X) max(X)]);
set(gca,'Ylim',[-.05 1.05],'Ytick',[0:.25:1]); 
%text(min(xx),0.7,{' p values:',[' Yint \neq 0: ' num2str(stats.p(1))],[' slope \neq 0: ' num2str(stats.p(2))]});
title(['N=' num2str(length(X)) ',p=' num2str(stats.p(2))],'HorizontalAlignment','center');


function violin(X,Y,varargin)
pd = fitdist(X,'kernel');
maxpdf = max(pd.pdf(X));
wscale = (.25/maxpdf)*sign(0.5-Y);
yjit = wscale*rand(size(X)).*pd.pdf(X); % jitter scaled by PDF fit from above line
plot(X,Y+yjit,varargin{:});
