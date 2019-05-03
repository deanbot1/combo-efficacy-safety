clear all 
close all

%% efficacy surface fit
cellline= 'Cell line MDA-MB-361';
load MDA-MB-361.mat
%MLN0128 % Plasma protein binding human=70.5 Mice=47.8
%MLN117 % Plasma protein binding human=82 Mice=53
cor0128=(100-47.8)/(100-70.5);
cor1117=(100-53)/(100-82);
TGI = TGI*100;
AUC1117 = A*cor1117;%converted to human equivalent cavg
AUC0128 = B*cor0128;%converted to human equivalent cavg
popAUC1117 = [];
popAUC0128 = [];
dose1117 = [];
dose0128 = [];
dosegrp = [];
%MLN0128 % Plasma protein binding human=70.5 Mice=47.8
%MLN117 % Plasma protein binding human=82 Mice=53



%% visulize INTERPOLATED surface
a = linspace(min(AUC1117),max(AUC1117),30);
b = linspace(min(AUC0128),max(AUC0128),30);

[X,Y] = meshgrid(a,b);
%Z=griddata(A,B,TGI,X,Y,'cubic');
Z=griddata(AUC1117,AUC0128,TGI,X,Y,'cubic');
surf(X,Y,Z,'EdgeColor','none') %interpolated
%surf(X,Y,Z) %interpolated
axis tight; hold on
plot3(AUC1117,AUC0128,TGI,'.','MarkerSize',15);
colorbar
view([-45, 30])
xlabel('MLN1117 C_{avg} (mg/L)')
h=get(gca,'xlabel');
set(h,'rotation',25)
ylabel('MLN0128 C_{avg} (mg/L)')
h=get(gca,'ylabel');
set(h,'rotation',-20)
zlabel('%GRI')
%surf(...,...,...,'EdgeColor','none')
title(cellline)

%% plot # of mice per dose (Cavg) group
figure
UU = unique([AUC1117 AUC0128],'rows');
for j = 1:size(UU,1)
    ii = find(AUC1117==UU(j,1) & AUC0128==UU(j,2));
    plot(UU(j,1),UU(j,2),'r+');
    text(UU(j,1),UU(j,2),[' ' num2str(length(ii))],'fontsize',12,'color','r'); hold on
end
set(gca,'Xtick',unique(AUC1117),'Xlim',[-1 7]);
set(gca,'Ytick',unique(AUC0128),'Ylim',[-.01 0.065]);
grid on
xlabel(plotanno.xlab);
ylabel(plotanno.ylab);
title('# mice at each exposure combo');

%% now attempt surface fitting with uncertainty

comboeff = @(mx,ey50,emaxy,beta,x,y)...
    mx*x + emaxy*y./(ey50+y) + beta*mx*x.*(emaxy*y./(ey50+y));

fito = fit([AUC1117,AUC0128],TGI,comboeff,...
    'StartPoint',[31 0.02 166 9e-5],...
    'Lower',[0 0 0 -Inf]);

figure;
plot(fito,[AUC1117,AUC0128],TGI);
xlabel('AUC1117'); ylabel('AUC0128'); zlabel('GRI');
colormap jet;
colorbar

%% now bootstrap...
Nboot = 2500;
Ndat = length(TGI);
fits = {};
params = zeros(Nboot,4);
for j = 1:Nboot
    ifit = ceil(Ndat*rand(size(TGI))); % bootstrap index
    fits = fit([AUC1117(ifit),AUC0128(ifit)],TGI(ifit),comboeff,...
    'StartPoint',[31 0.02 166 9e-5],...
    'Lower',[0 0 0 -Inf]);
    params(j,:) = coeffvalues(fits);
end

%% plots bootstrapping results
pnames = coeffnames(fito);
quants = quantile(params,[0.025 .5 .975]);
results = table(coeffvalues(fito)',quants(2,:)',quants(1,:)',quants(3,:)','Variablenames',{'best','Boot_med','CI_95_lo','CI_95_hi'},'Rownames',pnames)
writetable(results,[cellline '_fitparams.xlsx'],'writerownames',1);


%% not sure what the rest of this file does, or did...
stop
%clear all
%Reload data
%% setup
unidosegrp = unique(dosegrp);

for i = 1:length(unidosegrp)
    grp = unidosegrp(i);
    grpidx = dosegrp == grp;
    
    grp1117AUC(i,1) = mean(popAUC1117(grpidx));
    grp0128AUC(i,1) = mean(popAUC0128(grpidx));

    grpTGImean(i,1) = mean(TGI(grpidx));
    
end
    
%% fit surface to all data

% % linear
%   eff = @(b,x) b(1).*x(:,1) + b(2).*x(:,2) + b(3).*b(1).*b(2).*x(:,1).*x(:,2);
%   beta0 = [.02 1 0];

% saturating A (MLN1117)
%  eff = @(b,x) b(1).*x(:,1)./(b(3) + x(:,1)) + b(2).*x(:,2) + b(4).*b(1).*b(2).*x(:,1).*x(:,2)./(b(3) + x(:,1));
%  beta0 = [200 1 10 0];

% % saturating B (MLN0128)
  eff = @(b,x) b(1).*x(:,1) + b(2).*x(:,2)./(b(3) + x(:,2)) + b(4).*b(1).*b(2).*x(:,1).*x(:,2)./(b(3) + x(:,2));
  beta0 = [.003 .4 20 -1];

%    eff = @(b,x)exp(b(1)).*x(:,1) + exp(b(2)).*x(:,2)./(exp(b(3)) + x(:,2)) + b(4).*exp(b(1)).*exp(b(2)).*x(:,1).*x(:,2)./(exp(b(3)) + x(:,2));
%    beta0 = [log(5) log(25) log(0.01) 0];

 %  beta0 = [log(5 10 20 0];

   % % saturating A&B
% eff = @(b,x) b(1).*x(:,1)./(b(3) + x(:,1)) + b(2).*x(:,2)./(b(4) + x(:,2)) + b(5).*b(1).*b(2).*x(:,1).*x(:,2)./(b(3) + x(:,1))./(b(4) + x(:,2));
% beta0 = [100  100 0.003 0.003 0];


[beta,R,J,CovB,MSE] = nlinfit([AUC1117,AUC0128],TGI,@(b,x) eff(b,x),beta0);
%beta=[exp(beta(1))
betaci = nlparci(beta,R,'covar',CovB)

yfit = eff(beta,[AUC1117,AUC0128]);

R2 = 1 - sum(R.^2)/sum((yfit - mean(TGI)).^2)
AIC = length(TGI)* log(MSE)+ 2*length(beta)

figure
plot3(AUC1117,AUC0128,TGI,'k.','MarkerSize',20)
% scatter3(AUC1117,AUC0128,TGI,'filled')
 hold on
%xerr = [AUC1117 AUC1117]';
%yerr = [AUC0128 AUC0128]';
%zerr = [TGI-TGIstd, TGI+TGIstd]';

%plot3(xerr,yerr,zerr,'k-','linewidth',2)

% generate interpolation grid
a = linspace(min(AUC1117),max(AUC1117),10);
b = linspace(min(AUC0128),max(AUC0128),10);

[X,Y] = meshgrid(a,b);
iso = eff(beta,[X(:),Y(:)]);
iso = reshape(iso,10,10);

mesh(X,Y,iso,'FaceColor','interp','EdgeColor','none')
alpha(0.9)     % slightly transparent to see data points underneath surface
view([-45, 30])

xlabel('MLN1117 C_{avg} (mg/L)')
h=get(gca,'xlabel');
set(h,'rotation',25)
ylabel('MLN0128 C_{avg} (mg/L)')
h=get(gca,'ylabel');
set(h,'rotation',-20)
zlabel('%TGI')
title(cellline)

%%
beta(1)
beta(2)
beta(3)
beta(4)
beta(5)