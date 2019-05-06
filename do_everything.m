%% generates synthetic tox/efficacy data then fits it and generates graphs

%% clear all, close all
clear all
close all

Xname = 'Drug X (mg/L)';  % for graph labeling
Yname = 'Drug Y (mg/L)';  % for graph labeling

%% generate synthetic mouse exposure/efficacy data from 'truth model(s)'
% all exposures are assumed to be HUMAN free fraction for interspecies
% translation. Please see manuscript for how to do this correction.

XPK_mouse = [0 0 0.5 0.5 1  1]; % Y PK exposure for mouse groups (all the same because we can't collect PK from individual mice)
YPK_mouse = [0 0  0.75 0.75 1.5 1.5]; % Y PK exposure for mouse groups (all the same because we can't collect PK from individual mice)

[Xm,Ym] = meshgrid(XPK_mouse,YPK_mouse);
Xm = reshape(Xm,[],1);
Ym = reshape(Ym,[],1);

X_PK_eff = @(slope,X)slope*X;  % PK/efficacy function for drug X in terms of GRI
Y_PK_eff = @(Emax,EC50,Y)Emax*Y./(EC50+Y); % PK/efficacy function for drug Y in terms of GRI

Xslope = 100; % slope of GRI in drug X as single agent
YEmax = 300; % max GRI for drug Y
YEC50 = 1; % EC50 of efficacy of drug Y as single agent
GRInoise = 20; % sd of error in GRI observations
beta = .05; % synergy term

XY_PK_eff = @(Xslope,YEmax,YEC50,beta,x,y)...  % the combo efficacy equals...
	X_PK_eff(Xslope,x) + Y_PK_eff(YEmax,YEC50,y) + ... % sum of single agent effects
	beta*X_PK_eff(Xslope,x).*Y_PK_eff(YEmax,YEC50,y); % interaction term

GRIdata = XY_PK_eff(Xslope,YEmax,YEC50,beta,Xm,Ym)+GRInoise*randn(size(Xm));
figure; plot3(Xm,Ym,GRIdata,'ro'); xlabel(Xname);ylabel(Yname);zlabel('GRI');

%% now fit a surface to the efficacy data

fitobj = fit([Xm,Ym],GRIdata,XY_PK_eff,...
	'Lower',[0 0 0 -Inf],...  % make sure you let beta go negative if it wants to!
	'StartPoint',[50 100 3 0])

[XXm,YYm] = meshgrid(linspace(0,max(Xm),10),linspace(0,max(Ym),10));
figure; surf(XXm,YYm,fitobj(XXm,YYm));
colorbar
colormap jet

%% generate exposure/safety data in humans
% once again all concentrations are assumed to be in free fraction in human
% blood
Nh = 50; % number of humans in each of the 3 trials (single agent X, SA Y, and combox3)
% note we are assumign uniform sampling in concentration space, but this is not realistic if there is any overlapping efficacy
XAUC = [zeros(Nh,1);2*rand(Nh,1);2*rand(3*Nh,1)]; % exposure of drug X in each patient
YAUC = [3*rand(Nh,1);zeros(Nh,1);3*rand(3*Nh,1)]; % exposure of Y in each patient

igood = find(XAUC+2*YAUC < 3);
XAUC = XAUC(igood); YAUC = YAUC(igood);

b = -4; % Y intercept in logit (p(DLT)) space
mx = 3; % slope of X in logit(p(DLT))
my = 4; % slope of Y in logit(p(DLT))
alpa = 1; % intereaction term for tox

logit_tox = b + mx*XAUC + my*YAUC + alpa*XAUC.*YAUC;
p_tox = 1./(1+exp(-logit_tox)); % probability of tox of each patient
DLT = rand(size(p_tox)) < p_tox; % whether or not each patient had a DLT


%% plot tox unidimensionally ALL SUBJECTS
% if matlab throws an error below, that means either all subjects had DLT's or all subjects didn't
% have DLT's. Just run the section above again, and hopefully the random
% data will include both DLT and non-DLT patients. 

scale = 'loglin';
figure
PTOX = 0.25; % probability threshold to consider as max tolerated exposure
subplot(1,2,1);unilogifit(XAUC(YAUC==0),DLT(YAUC==0),scale,PTOX); xlabel(Xname);
subplot(1,2,2);unilogifit(YAUC(XAUC==0),DLT(XAUC==0),scale,PTOX); xlabel(Yname);

%% "clean" the data set of NaNs
Igood = ~isnan(DLT.*XAUC.*YAUC);
XX = XAUC(Igood); 
YY = YAUC(Igood); 
DLT = DLT(Igood);

%% Plot  clean data set in 2d

figure
clear h
h(1) = plot(XX(DLT==0),YY(DLT==0),'go','MarkerFaceColor','g'); hold on;
h(2) =  plot(XX(DLT==1),YY(DLT==1),'ro','MarkerFaceColor','r');

xlabel(Xname);
ylabel(Yname);
legend(h,{'no DLT','DLT'});
grid on

%% now we fit a Loewe-like model to the data and boostrap
Nboot = 2500; % number of bootstraps to run

[abest_loewe,aboot_loewe,figs_loewe] = fit_loewe_combo_tox_data2(YY,XX,DLT,0.25,Nboot,Yname,Xname,0);
C = aboot_loewe.MTE_curve;
%save(sprintf('Loewe_MTE_%s.mat',titlestr),'C');

%set(gca,'Ylim',[-.5 18]/1000,'Xlim',[-100 4500]/1000)
grid on

% make a table summarising the bootstrapping results
if Nboot > 0
    quants = [0.05 0.5 0.95];
    bootq = quantile(aboot_loewe.bboot',quants)';
    bootresults = table(abest_loewe.bbest,bootq(:,2),bootq(:,1),bootq(:,3),...
        'VariableNames',{'bestfit','boot_median','CI90_lo','CI90_hi'},...
        'RowNames',{'b','mx','my','alpha'})
    %writetable(bootresults,'loewe_params.xlsx','WriteRowNames',1);
end


%% now we overlay efficacy onto MTE curve derived above and slice through using polar coordinates
XXX = C(1,2:end)'; % X coordinates of MTE curve
YYY = C(2,2:end)'; % Y coordinates of MTE curve

figure;
surf(XXm,YYm,fitobj(XXm,YYm));
colorbar
colormap jet
hold on; 
ZZZ = fitobj(XXX,YYY); 
ibest = find(ZZZ==max(ZZZ));  Zbest = fitobj(XXX(ibest),YYY(ibest));
plot3(XXX,YYY,ZZZ,'m-','LineWidth',3); xlabel(Xname); ylabel(Yname); zlabel('GRI');
plot(XXX,YYY,'m-');
plot(XXX(ibest),YYY(ibest),'mp','MarkerFaceColor','m');
stem3(XXX(ibest),YYY(ibest),ZZZ(ibest),'mp','MarkerFaceColor','m');
alpha 0.5 % make surface semi-transparent


% now the slice thing
theta = atan(YYY./XXX);
figure; 
subplot(4,1,1:2);
plot(theta,ZZZ,'m-','LineWidth',3); ylabel('GRI');  grid on; hold on
plot(theta(ibest),ZZZ(ibest),'mp','MarkerFaceColor','m','MarkerSize',10);
title('optimal exposure combination');
subplot(4,1,3); 
area(theta,XXX,'FaceColor','c'); ylabel(Xname); hold on
plot(theta(ibest),XXX(ibest),'mp','MarkerFaceColor','m'); 
text(theta(ibest),XXX(ibest),[' ' num2str(XXX(ibest)) Xname],'color','m');
subplot(4,1,4);
area(theta,YYY,'FaceColor','c'); ylabel(Yname); hold on
plot(theta(ibest),YYY(ibest),'mp','MarkerFaceColor','m'); 
text(theta(ibest),YYY(ibest),[' ' num2str(YYY(ibest)) Yname],'color','m');
xlabel('\theta=tan^{-1}([Y]/[X])');
% the magenta star is optimal tolerated exposure combination


