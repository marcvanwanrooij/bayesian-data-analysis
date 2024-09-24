%% Setup
% - Create repository folders, explain what these folders are for. DAC, RDC, 
% DSC. - Store data in DAC - Create matlab folder - create a tmp.m file (temporary, 
% at the end of the session this should be renamed
%% Initialization
% when running the script, choose change folder

close all;
clear;

% running the script will prompt the message: tmp.m cannot be found on
% path. It will ask to do something. Choose: change folder.
%% Load Data
% Start with absolute path, fname = '/Users/marcw/Library/CloudStorage/Dropbox/Onderwijs/course 
% - Analysis Workshop/DAC/NH-0115-22-10-07-data.mat'; % absolute and then explain 
% relative path. Now the folder structure and the change folder make sense: everyone 
% with the same folder structure can run this script, and it will work! Not so 
% for absolute paths or when you add this tmp.m folder to Matlab's path.

fname = fullfile('..','datasets','tactile_data_one_participant.xlsx'); % file name, relative path
T = readtable(fname);

% Explain the experiment, explain data (rows individual trials, columns
% different parameters). Note: this is tidydata (not raw).
T(1,:)
return
%% Graphics
% response azimuth vs stimulus azimuth

y = Data(:,5); % response azimuth
y = y+50*randn(size(y));
x = Data(:,1); % stimulus  azimuth

% Build this up: first just plot, and discuss what would make this figure
% better
figure(1);
clf;
plot(x,y,'ko','MarkerFaceColor','k','MarkerEdgeColor','w','Markersize',13,'LineWidth',2);
xlabel('stimulus azimuth (deg)');
ylabel('response azimuth (deg)');
% set(gca,'FontSize',20);
set(gca,'Xtick',-90:30:90,'Ytick',-90:30:90);
xlim([-100 100]);
ylim([-100 100]);
% axis square;
% box off;
nicegraph;

return
%% Doing Frequentist linear regression

stats	= regstats(y,x,'linear',{'beta','tstat','rsquare','fstat'});
beta	= stats.beta; % beta(1) = offset, beta(2) = slope
pval = stats.tstat.pval;
% what do these parameters mean? Are they meaningful?
% And continue with 'prediction'
xi		= -90:90;
ypred	= beta(2)*xi+beta(1);

figure(1)
hold on;
plot(xi,ypred,'-','LineWidth',3);
unityline;
horline;
verline;
%% Doing Bayesian analysis
% Show that Bayesian analysis produces the same point estmate but does not produce 
% a single number, it produces a distribution

mcmc = regjags(y,x);

betamcmc = [mean(mcmc.beta0) mean(mcmc.beta1)]';
% beta


figure(42)
clf
plotpost(mcmc.beta1,'xlim',[0.8 1.2]);
nicegraph;