%% Setup
% - Create repository folders, explain what these folders are for. DAC, RDC,
% DSC. 
% - Store data in DAC
% - Create matlab folder
% - create a tmp.m file (temporary, at the end of the session this should
% be renamed


%% Initialization
% when running the script, choose change folder
close all;
clear;

% running the script will prompt the message: tmp.m cannot be found on
% path. It will ask to do something. Choose: change folder.

%% Load Data
% Start with absolute path,
% fname = '/Users/marcw/Library/CloudStorage/Dropbox/Onderwijs/course -
% Analysis Workshop/DAC/NH-0115-22-10-07-data.mat'; % absolute
% and then explain relative path. Now the folder structure and the change
% folder make sense: everyone with the same folder structure can run this
% script, and it will work! Not so for absolute paths or when you add this
% tmp.m folder to Matlab's path.
fname = fullfile('..','..','DAC','NH-0115-22-10-07-data.mat'); % file name, relative path
load(fname);

% Explain the experiment, explain data (rows individual trials, columns
% different parameters). Note: this is tidydata (not raw):
% https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html
%
% The data consists of a N x 6 matrix whose columns contain:
% stimulus azimuth | stimulus elevation | stimulus frequency-band (1=BB, 2=HP, 3=LP) | stimulus level | response azimuth | response elevation, and each row indicates trial number.


%% Graphics
F	= Data(:,3);
uF	= unique(F);

% uF = [4 6 9]
% sel = F==3; % selection vector, selecting BB sounds
% sel
% whos sel
% N = sum(sel);

% F(sel)
% find(sel)

figure(1);
clf;
titlestr = {'BB','HP','LP'};
for ii = 1:numel(uF)
	[ii uF(ii)]
	sel = F==uF(ii);
	% response azimuth vs stimulus azimuth
	y = Data(sel,6); % response azimuth
	% y = y+50*randn(size(y));
	x = Data(sel,2); % stimulus  azimuth
	
	%% graphics
	plotstimres(x,y,ii);
	title(titlestr{ii});
	
	
	%% Statistics
	stats	= regstats(y,x,'linear',{'beta','tstat','rsquare','fstat'});
beta	= stats.beta; % beta(1) = offset, beta(2) = slope
str = ['y = ' num2str(beta(2),2) 'x + ' num2str(beta(1),2)];

	bf_text(0.4,0.9,str,'FontSize',24);
pval = stats.tstat.pval(2);
str = ['p = ' num2str(pval) ];

	bf_text(0.2,0.8,str,'FontSize',24);
	
	%% 
	mcmc(ii) = regjags(y,x);

end



%% 

figure(42)
clf
for ii = 1:3
	b = mcmc(ii);
	subplot(2,3,ii)
	plotpost(b.beta1,'xlim',[-0.2 1.2] );
end

subplot(234)
plotpost(mcmc(1).beta1,'xlim',[-0.2 1.2],'showCurve',true);
hold on
plotpost(mcmc(2).beta1,'xlim',[-0.2 1.2],'showCurve',true,'Color','r');

subplot(235);
plotpost(mcmc(1).beta1-mcmc(2).beta1,'xlim',[-1 1]);

%%

	y = Data(:,6); % response elevation
	% y = y+50*randn(size(y));
	x = Data(:,2); % stimulus  elevation
	z = F;
samples = hregjags(x,y,z);

%%
%% 

figure(43)
clf
for ii = 1:3
	b = samples.beta1(:,ii);
	subplot(2,3,ii)
	plotpost(b,'xlim',[-0.2 1.2] );
end

subplot(234)
plotpost(samples.beta1(:,1),'xlim',[-0.2 1.2],'showCurve',true);
hold on
plotpost(samples.beta1(:,2),'xlim',[-0.2 1.2],'showCurve',true,'Color','r');

subplot(235);
plotpost(samples.beta1(:,1)-samples.beta1(:,2),'xlim',[-1 1]);

