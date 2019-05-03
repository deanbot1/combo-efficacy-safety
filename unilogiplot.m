function [b,dev,stats]=unilogiplot(X,Y,scale)

switch scale
	case 'linear'
		xf = @(x)(x);
	case 'log'
		xf = @(x)(log(x));
end

plot(X,Y,'ko','MarkerFaceColor','k'); hold on;
[b,dev,stats] = glmfit(xf(X),Y,'binomial');
xx = [min(X):(max(X)-min(X))/100:max(X)];
[yhat,dylo,dyhi] = glmval(b,xf(xx),'logit',stats);
plot(xx,yhat,'r-','linewidth',3);
plot(xx,yhat-dylo,'r-');plot(xx,yhat+dyhi,'r-');
set(gca,'Xscale',scale)
ylabel('p(DLT)');
axis tight
%title(scale);
plot(xx([1 end]),.25*[1 1],'k:');
%text(xx(1),0.25,'0.25 ','HorizontalAlignment', 'right');
set(gca,'Ytick',[0:.25:1]);
%text(min(xx),0.7,{' p values:',[' Yint \neq 0: ' num2str(stats.p(1))],[' slope \neq 0: ' num2str(stats.p(2))]});
title([' p(obs|slope=0)=' num2str(stats.p(2))],'Color','r');
