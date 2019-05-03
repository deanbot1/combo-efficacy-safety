%% loads tox data, and who knows, maybe it plots it too... and also models it

%% clear all, close all
clear all
close all

milledonlyflag = 0 % ********** set to 1 if you want MILLED ONLY, -1 unmilled only, 0 all patients ***********

Xname = 'TAK-228 C_{ave} (mg/L)';
Yname = 'TAK-117 C_{ave} (mg/L)';

%% Tmike is incorrect, fix it:

% Cave = (AUC24_raw * number of doses in a cycle / number of days in cycle)/24 
% ng*hr/ml * (ndoses/cycle) / (ndays/cycle) / hours = ng/ml
% mayank's modeling is in mg/L = ug/ml = 1000ng/ml


%% read in PK data using readtable and fix it
Tmike = readtable('../data/Dean_Worksheet.xls','sheet','Sheet1'); % mike's data
Tkey = readtable('../data/Dean_Worksheet.xls','sheet','sched_key'); % key
ndoseperweek = @(sched)Tkey{find(strcmp(Tkey.Schedule_name,sched)),'ndoseperweek'};
Nrows = size(Tmike,1);
%[AUC_1117,AUC_0128] = deal(NaN(Nrows,1));
AUC_1117_raw = Tmike.AUC24_MLN1117_UNMODIFIED;
zeros = find(AUC_1117_raw < eps); 
AUC_1117_raw(zeros) = Tmike.AUCLST_MLN1117(zeros);
AUC_0128_raw = Tmike.AUC24_MLN0128_UNMOIDIFIED;
zeros = find(AUC_0128_raw < eps); 
AUC_0128_raw(zeros) = Tmike.AUCLST_MLN0128(zeros);
[ndose_1117,ndose_0128] = deal(NaN(Nrows,1));
for j = 1:size(Tmike,1) 
    ndose_1117(j) = ndoseperweek(char(Tmike{j,'SCHED1117'}));
    ndose_0128(j) = ndoseperweek(char(Tmike{j,'SCHED128'}));
end

Tmike.AUC_0128 = AUC_0128_raw.*ndose_0128/7/24;
Tmike.AUC_1117 = AUC_1117_raw.*ndose_1117/7/24;

weird228 = find(Tmike.AUC_0128_BAD>eps & Tmike.AUC_0128< eps);

figure; plot(Tmike.AUC_0128_BAD/24,Tmike.AUC_0128,'ro'); xlabel('old'); ylabel('new'); title('TAK-228'); hold on; plot([0 90],[0 90],'k-');
figure; plot(Tmike.AUC_1117_BAD/24,Tmike.AUC_1117,'ro'); xlabel('old'); ylabel('new'); title('TAK-117'); hold on; plot([0 7000],[0 7000],'k-');

%% read in more data (tox data)

Trach = readtable('../data/C32001DoseEscForClinPharm.csv'); % rachel's data
length(find(strcmp(Tmike.STUDY,'C32001')))
studies = unique(Tmike.STUDY);

ikeep = [];
for j = 1:length(studies)
    study = studies{j};
    istud = find(strcmp(Tmike.STUDY,study));
    idlist = Tmike(istud,:).SUBJID;
    [~,uidi] = unique(idlist);
    ikeep = union(ikeep,istud(uidi));
end

Tmike = Tmike(ikeep,:);

%Tcombo = join(Trach,Tmike,'keys','SUBJID')
foo = cell(height(Trach),1);
[foo{:}] = deal('C32001');
Trach.STUDY = foo;
Trm = join(Trach,Tmike, 'keys',{'STUDY','SUBJID'});
IDrm = Trm.SUBJID;
T = Tmike((strcmp(Tmike.STUDY,'C32001') & ismember(Tmike.SUBJID,IDrm))|~strcmp(Tmike.STUDY,'C32001'),:);
% check to make sure they have the same info
T.AUC_0128(isnan(T.AUC_0128)) = 0;
T.AUC_1117(isnan(T.AUC_1117)) = 0;
Tmike.AUC_0128(isnan(Tmike.AUC_0128)) = 0;
Tmike.AUC_1117(isnan(Tmike.AUC_1117)) = 0;

% make a table of how many patients in which set

clear summtab
for j = 1:length(studies)
    clear N_auc N_safety N_both
    study = studies{j};
    switch study
        case 'C32001'
            auctest = @(tab)(tab.AUC_0128.*tab.AUC_1117>0);
        otherwise
            auctest = @(tab)(tab.AUC_0128+tab.AUC_1117>0);
    end
    studynice = strrep(study,'-','_');
    rnames = {'Milled','Unmilled','Total'};
    N_auc(1,1) = length(find(Tmike.Milled==1 & strcmp(Tmike.STUDY,study) & auctest(Tmike)));
    N_auc(2,1) = length(find(Tmike.Milled==0 & strcmp(Tmike.STUDY,study) & auctest(Tmike)));
    N_auc(3,1) = length(find(strcmp(Tmike.STUDY,study) & auctest(Tmike)));    
    N_safety(1,1) = length(find(Tmike.Milled==1 & strcmp(Tmike.STUDY,study)));
    N_safety(2,1) = length(find(Tmike.Milled==0 & strcmp(Tmike.STUDY,study)));
    N_safety(3,1) = length(find(strcmp(Tmike.STUDY,study)));
    N_dlt(1,1) = length(find(Tmike.Milled==1 & strcmp(Tmike.STUDY,study) & Tmike.DLT==1));
    N_dlt(2,1) = length(find(Tmike.Milled==0 & strcmp(Tmike.STUDY,study)& Tmike.DLT==1));
    N_dlt(3,1) = length(find(strcmp(Tmike.STUDY,study)& Tmike.DLT==1));
    N_both(1,1) = length(find(T.Milled==1 & strcmp(T.STUDY,study) & auctest(T)));
    N_both(2,1) = length(find(T.Milled==0 & strcmp(T.STUDY,study) & auctest(T)));
    N_both(3,1) = length(find(strcmp(T.STUDY,study) & auctest(T)));
    N_bdlt(1,1) = length(find(T.Milled==1 & strcmp(T.STUDY,study) & auctest(T) & T.DLT));
    N_bdlt(2,1) = length(find(T.Milled==0 & strcmp(T.STUDY,study) & auctest(T) & T.DLT));
    N_bdlt(3,1) = length(find(strcmp(T.STUDY,study) & auctest(T) & T.DLT));
    summtab.(studynice) = table(N_auc,N_safety,N_both,N_dlt,N_bdlt,'Rownames',rnames);
    disp(['************** ' study ' ***************']);
    disp(summtab.(studynice));
end

%% cleaner summary table ignoring milled/unmilled and counting dlts
Tmega = table('Size',[length(studies) 6],'VariableTypes',{'string','double','double','double','double','double'},...
    'VariableNames',{'STUDY','N_auc','N_safety','N_both','N_dlt','N_bdlt'});

Tmega.STUDY = studies;
for j = 1:length(studies)
     studynice = strrep(studies{j},'-','_');
     Tmega{j,2:end} = summtab.(studynice){'Total',:};
end
Tmega.PCT_sDLT = round(100*Tmega.N_dlt./Tmega.N_safety);
Tmega.PCT_bDLT = round(100*Tmega.N_bdlt./Tmega.N_both);


disp(Tmega);

%% process the data

clear pat 

for j = 1:size(T,1);
    pat(j).SUBJID = T.SUBJID{j};
    pat(j).STUDY = T.STUDY{j};
    pat(j).DLTNAME = T.DLT1{j};
    pat(j).DLT = T.DLT(j);
    pat(j).X_Cmax = T.CMAX_MLN0128(j)/1000;
    pat(j).Y_Cmax = T.CMAX_MLN1117(j)/1000;
    pat(j).X_AUClast = T.AUC_0128(j)/1000;
    pat(j).Y_AUClast = T.AUC_1117(j)/1000;
    pat(j).MILLED = T.Milled(j);
    pat(j).X_DOSE = T.dose128(j);
    pat(j).Y_DOSE = T.DOSE1117(j);
end

% keep only the milled or no 128 patients
for j = 1:length(pat)
    if isnan(pat(j).X_Cmax), pat(j).X_Cmax=0;end
    if isnan(pat(j).Y_Cmax), pat(j).Y_Cmax=0;end
    if isnan(pat(j).X_AUClast), pat(j).X_AUClast=0;end
    if isnan(pat(j).Y_AUClast), pat(j).Y_AUClast=0;end
    if isnan(pat(j).X_DOSE),pat(j).X_DOSE=0;end
    if isnan(pat(j).Y_DOSE),pat(j).Y_DOSE=0;end
end



XAUC = cat(1,pat.X_AUClast);
YAUC = cat(1,pat.Y_AUClast);
XC = cat(1,pat.X_Cmax);
YC = cat(1,pat.Y_Cmax);

allp = length(pat);
pat = pat((XAUC+YAUC>0)); % filter out all pats that have no AUC's!
aucp = length(pat);

disp(sprintf('keep the %d of %d patients who have nonzero X or Y AUC values', aucp,allp));

XAUC = cat(1,pat.X_AUClast);
YAUC = cat(1,pat.Y_AUClast);

iY = find(cat(1,pat.X_DOSE)==0);
iX = find(cat(1,pat.Y_DOSE)==0);
iXY = find(cat(1,pat.Y_DOSE)>0 & cat(1,pat.X_DOSE)>0);

MILLED = cat(1,pat.MILLED);
%N117 = length(find(XAUC==0 & YAUC>0));
%N128 = length(find(YAUC==0 & XAUC>0));
%Nbot = length(find(XAUC>0 & YAUC>0)); 
N117 = length(iY);
N128 = length(iX);
Nbot = length(iXY);

disp(sprintf('starting with %d patients',length(pat)));
%disp(sprintf('of these, %d are 117 only, %d are 128 only, and %d are combinations',N117,N128,Nbot));

switch milledonlyflag
case 1 
ikeep = find(MILLED | XAUC==0); 
npatall = length(pat);
pat = pat(ikeep);
npatmil = length(pat);
disp(sprintf('milled 128 only: keeping %d of %d patients',npatmil,npatall));
case -1 
ikeep = find(~MILLED | XAUC==0); 
npatall = length(pat);
pat = pat(ikeep);
npatmil = length(pat);
disp(sprintf('unmilled 128 only: keeping %d of %d patients',npatmil,npatall));
end

XAUC = cat(1,pat.X_AUClast);
YAUC = cat(1,pat.Y_AUClast);
N117 = length(find(XAUC==0 & YAUC>0));
N128 = length(find(YAUC==0 & XAUC>0));
Nbot = length(find(XAUC>0 & YAUC>0));
disp(sprintf('of these, %d are 117 only, %d are 128 only, and %d are combinations',N117,N128,Nbot));


%% plot things


%XAUC=XAUC/(28);%converted to cavg %% DEBUG WHY 28 and not 24??
%YAUC=YAUC/(28);%converted to cavg %% DEBUG WHY 28 and not 24??


XC = cat(1,pat.X_Cmax);
YC = cat(1,pat.Y_Cmax);
PTOX = 0.25;

DLT = cat(1,pat.DLT);
Ino = find(DLT==0);
Iyes = find(DLT==1);

%% plot statistics of exposures

XCAVENORM = XAUC./cat(1,pat.X_DOSE);
XCMAXNORM = XC./cat(1,pat.X_DOSE);

figure('Position',[274   445   383   420]);
switch milledonlyflag
    case 0
        titlestr = 'allPatients';
    case 1
        titlestr = 'milledOnly';
    case -1
        titlestr = 'unmilledOnly';
end

subplot(2,1,2); 
histogram(log10(XCAVENORM(XAUC>0)),[-2:.2:2]);
xlabel('C_{ave}/dose (ng/ml/mg)');
set(gca,'Xticklabel',num2str(10.^get(gca,'Xtick')','%4.4g'));

subplot(2,1,1);
histogram(log10(XCMAXNORM(XC>0)),[-2:.2:2]);
xlabel('C_{max}/dose (ng/ml/mg)');
set(gca,'Xticklabel',num2str(10.^get(gca,'Xtick')','%4.4g'));
title(titlestr);


%% plot unidimensionally ALL SUBJECTS
scale = 'loglin';
figure
subplot(2,2,1);igood=YC==0&~isnan(XC)&XC>0;unilogifit(XC(igood),DLT(igood),scale,PTOX); xlabel('TAK-228 Cmax(mg/L)');
subplot(2,2,2);unilogifit(YC(XC==0 & YC>0),DLT(XC==0 & YC>0),scale,PTOX); xlabel('TAK-117 Cmax (mg/L)');
%subplot(2,2,3);unilogifit(XAUC(YAUC==0),DLT(YAUC==0),scale,PTOX); xlabel('0128 AUClast');
%subplot(2,2,4);unilogifit(YAUC(XAUC==0),DLT(XAUC==0),scale,PTOX); xlabel('1117 AUClast');
subplot(2,2,3);unilogifit(XAUC(YAUC==0),DLT(YAUC==0),scale,PTOX); xlabel('TAK-228 Cavg (mg/L)');
subplot(2,2,4);unilogifit(YAUC(XAUC==0),DLT(XAUC==0),scale,PTOX); xlabel('TAK-117 Cavg (mg/L)');

%% "clean" the data set
%for j = 1:length(pat); stud{j} = pat(j).STUDY; end
%Imilled = strcmp(stud','INK1117-001')|strcmp(stud','C32001');
Igood = ~isnan(DLT.*XAUC.*YAUC);
%Igood = Imilled & Igood;
 

XX = XAUC(Igood); 
YY = YAUC(Igood); 
DLT = DLT(Igood);

%% Plot  clean data set in 2d

figure
clear h
h(1) = plot(XX(DLT==0),YY(DLT==0),'go','MarkerFaceColor','g'); hold on;
h(2) =  plot(XX(DLT==1),YY(DLT==1),'ro','MarkerFaceColor','r');

xlabel('TAK-228 C_{avg} (mg/L)');
ylabel('TAK-117 C_{avg} (mg/L)');

title(titlestr);
legend(h,{'no DLT','DLT'});
set(gca,'Xlim',[-.5 30]/1000,'Ylim',[-100 4500]/1000)
grid on
colorbar
colormap pink
caxis([0 1])
%% how many bootstraps?
Nboot = 2500;

%% Now try a simple 'bliss boosting' combo model

% [abest_bliss,aboot_bliss,figs_bliss] = fit_bliss_combo_tox_data2(XX,YY,DLT,0.25,Nboot,Xname,Yname,-1);
% C = aboot_bliss.MTE_curve;
% save(sprintf('Bliss_MTE_%s.mat',titlestr),'C');

%% lowewe model

[abest_loewe,aboot_loewe,figs_loewe] = fit_loewe_combo_tox_data2(YY,XX,DLT,0.25,Nboot,Yname,Xname,0);
C = aboot_loewe.MTE_curve;
save(sprintf('Loewe_MTE_%s.mat',titlestr),'C');

set(gca,'Ylim',[-.5 18]/1000,'Xlim',[-100 4500]/1000)
grid on

%% make a table summarising the bootstrapping results
if Nboot > 0
    quants = [0.05 0.5 0.95];
    bootq = quantile(aboot_loewe.bboot',quants)';
    bootresults = table(abest_loewe.bbest,bootq(:,2),bootq(:,1),bootq(:,3),...
        'VariableNames',{'bestfit','boot_median','CI90_lo','CI90_hi'},...
        'RowNames',{'b','mx','my','alpha'})
    writetable(bootresults,'loewe_params.xlsx','WriteRowNames',1);
end



%% Thall model

%[abest_thall,aboot_thall,figs_thall] = fit_thall_combo_tox_data(XX,YY,DLT,0.25,Nboot,Xname,Yname,[0;1]);

%% to overlay efficacy, run FinalCLINICALTGIsimulation.m separately, editing it to point to the mat file created above



%% Devide surface in 9 blocks and plot binomial mean and cofidenece interval for both observed and predicted
% Observed data plotted as open circle and pred is shown as filled circle
% if Pred is outside observed 95%CI circle will be filled with red unless green

%Bliss
%     if isempty(aboot_bliss.PPPP) == false
%         [ comarea_mx, comarea_my, compFig , pci] = compPlot_data( XX, YY, DLT );
%         [ptox_best]=compPlot_model(comarea_mx, comarea_my, aboot_bliss, compFig, pci); 
%     end
%     
% %% Loewe
% 
%    if isempty(aboot_loewe.PPPP) == false
%         [ comarea_mx, comarea_my, compFig , pci] = compPlot_data( XX, YY, DLT );
%         [ptox_best]=compPlot_model(comarea_mx, comarea_my, aboot_loewe, compFig, pci); 
%     end