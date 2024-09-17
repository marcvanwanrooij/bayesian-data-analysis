close all;
clearvars;
clc;

% Bayesian linear regression with JAGS
addpath(fullfile('..','DAC'))

% Load data for all participants
load LocalizationPoolData.mat

% For this excercise, we will use broadband only and azimuth
% target/response

% stimulus azimuth | stimulus elevation | 
% stimulus frequency-band (1=BB, 2=HP, 3=LP) | 
% stimulus level | response azimuth | 
% response elevation, and each row indicates trial number.
sel = data(:,3) == 1;

x = data(sel,1); % stimulus azimuth
y = data(sel,5); % response azimuth
s = data(sel,8); % subject number

uS = unique(s); % unique subjects

%% A) Simple model for linear regression 
% Step 1) Identify the data relevant to the research questions 
% Look at the data, identify what is the independent/dependent
% variables

figure(1)
for ns = 1:length(uS)

    subplot(2,7,ns)
    sel = s == uS(ns);

    plot(x(sel),y(sel),'ko')

    xlim([-100,100]);
    ylim([-100,100]);
    axis square
    title(['S' num2str(uS(ns))])
    horline;
    verline;
    unityline;

    xlabel('target (deg)')
    ylabel('response (deg)')
    box off;

end


% Step 2) Define a descriptive model for the relevant data.
%
% y~dnorm(mu,sigma)
% mu = beta0 + beta1*x

% Step 3) Specify a prior distribution on the parameters.
%

% Step 4) Use Bayesian inference to re-allocate credibility across 
% parameter values. Interpret the posterior distribution on parameters
% Inference
% (Un)comment to choose the excercise 
type = 'simple';
% type = 'robust';
% type = 'hierarchical';
samples = getsamples_jags(x,y,s,type);

% plot posterior draws
figure(2)
subplot(121)
plotpost(samples.beta0);
xlabel('\beta_0: bias')
xlim([-5,5])
nicegraph;

subplot(122)
plotpost(samples.beta1);
xlabel('\beta_1: gain')
xlim([0.6,1.1])
nicegraph;

% Step 5) Check that the posterior predictions mimic the data with 
% reasonable accuracy
figure(3)
sel = s == 1;

hold on

% plot 80 credible lines
vv = floor(linspace(1,length(samples.beta0),100));
for nv = 1:length(vv)
    plot(-90:90,...
        samples.beta0(vv(nv))+...
        samples.beta1(vv(nv))*(-90:90),...
        'color',[.5,.5,.5,.6],'LineWidth',1)
end
% plot best fit line
plot(-90:90,...
    mean(samples.beta0)+mean(samples.beta1)*(-90:90),...
    'color',[.4,.1,.1],'LineWidth',2)
% plot data
plot(x(sel),y(sel),'ko','MarkerSize',8)
% making this look nicer
xlim([-100,100]);
ylim([-100,100]);
axis square
horline;
verline;
unityline;
title('simple fit')

xlabel('target (deg)')
ylabel('response (deg)')
box off;

%% B) Robust model for linear regression 
% Here, we change two things:
% 1) the noise distribution is now a t-student, to accomodate for outliers
% 2) To improve sampling efficiency we standardized the data.

% Sampling
type = 'robust';
% type = 'hierarchical';
samples = getsamples_jags(x,y,s,type);

% plot posterior draws
figure(4)
subplot(121)
plotpost(samples.beta0);
xlabel('\beta_0: bias')
xlim([-5,5])
nicegraph;

subplot(122)
plotpost(samples.beta1);
xlabel('\beta_1: gain')
xlim([0.6,1.1])
nicegraph;

% plot the fit..although for this data, there aren't many outliers
figure(5)
sel = s == 1;

hold on

% plot 80 credible lines
vv = floor(linspace(1,length(samples.beta0),100));
for nv = 1:length(vv)
    plot(-90:90,...
        samples.beta0(vv(nv))+...
        samples.beta1(vv(nv))*(-90:90),...
        'color',[.5,.5,.5,.6],'LineWidth',1)
end
% plot best fit line
plot(-90:90,...
    mean(samples.beta0)+mean(samples.beta1)*(-90:90),...
    'color',[.4,.1,.1],'LineWidth',2)
% plot data
plot(x(sel),y(sel),'ko','MarkerSize',8)
% making this look nicer
xlim([-100,100]);
ylim([-100,100]);
axis square
horline;
verline;
unityline;
title('robust fit')

xlabel('target (deg)')
ylabel('response (deg)')
box off;

%% C) Hierarchical model for linear regression 
% Here, we extend the previous model to include subjects (or participants)
% and we make the model hierarchical (we set priors on priors), to capture
% population parameters

% Sampling
type = 'hierarchical';
samples = getsamples_jags(x,y,s,type);

% plot posterior draws
figure(5)
subplot(221)
plotpost(samples.beta0mu,'color','k');
xlabel('\mu_{\beta_{0}}: bias')
title('population mean bias')
xlim([-7,7])
nicegraph;
set(gca,'fontsize',14)

subplot(222)
plotpost(samples.beta1mu,'color','k');
xlabel('\mu_{\beta_{1}}: gain')
title('population mean gain')
xlim([0.6,1.2])
nicegraph;
set(gca,'fontsize',14)

subplot(223)
for ns = 1:length(uS)
    [yy,xx]= ksdensity(samples.beta0(:,ns));
    plot(xx,yy,'color',[ns./(length(uS)+1),.1,.1,.6],'linew',2);
    hold on
end
xlabel('\beta_{0,s}: bias|s')
xlim([-7,7])
nicegraph;
set(gca,'fontsize',14)

subplot(224)
for ns = 1:length(uS)
    [yy,xx]= ksdensity(samples.beta1(:,ns));
    plot(xx,yy,'color',[ns./(1+length(uS)),.1,.1,.6],'linew',2);
    hold on
end
xlabel('\beta_{1,s}: gain|s')
xlim([0.6,1.2])
nicegraph;
set(gca,'fontsize',14)

figure(6)
uS = unique(s);
for ns = 1:length(uS)

    subplot(2,7,ns)
    sel = s == uS(ns);
    
    plot(x(sel),y(sel),'ko')

    xlim([-100,100]);
    ylim([-100,100]);
    axis square
    title(['S' num2str(uS(ns))])
    horline;
    verline;
    unityline;

    xlabel('target (deg)')
    ylabel('response (deg)')
    box off;


    hold on

    vv = floor(linspace(1,1000,80));
    for nv = 1:length(vv)
        plot(-90:90,...
            mean(samples.beta0(vv(nv),ns))+...
            mean(samples.beta1(vv(nv),ns))*(-90:90),...
            'color',[.1,.1,.1,.5],'LineWidth',1)
    end

    plot(-90:90,...
        mean(samples.beta0(:,ns))+...
        mean(samples.beta1(:,ns))*(-90:90),...
        'LineWidth',1.5)

end

%% D) BONUS: Hierarchical models 
% What if we add a participant (S14) with a single response? What would the
% model predict? Will it be able to fit? 
% Answer is yes, using the information from other subjects!

% We first add a target at 30 deg. with a response at 45 deg. for S14
x = [x;30]; % add target
y = [y;45]; % add response
s = [s;14]; % add subject
% Sampling
type = 'hierarchical';
samples = getsamples_jags(x,y,s,type);

% plot posterior draws
figure(7)
subplot(221)
plotpost(samples.beta0mu,'color','k');
xlabel('\mu_{\beta_{0}}: bias')
title('population mean bias')
xlim([-7,7])
nicegraph;
set(gca,'fontsize',14)

subplot(222)
plotpost(samples.beta1mu,'color','k');
xlabel('\mu_{\beta_{1}}: gain')
title('population mean gain')
xlim([0.6,1.2])
nicegraph;
set(gca,'fontsize',14)

subplot(223)
for ns = 1:length(uS)
    [yy,xx]= ksdensity(samples.beta0(:,ns));
    plot(xx,yy,'color',[ns./(length(uS)+1),.1,.1,.6],'linew',2);
    hold on
end
xlabel('\beta_{0,s}: bias|s')
xlim([-7,7])
nicegraph;
set(gca,'fontsize',14)

subplot(224)
for ns = 1:length(uS)
    [yy,xx]= ksdensity(samples.beta1(:,ns));
    plot(xx,yy,'color',[ns./(1+length(uS)),.1,.1,.6],'linew',2);
    hold on
end
xlabel('\beta_{1,s}: gain|s')
xlim([0.6,1.2])
nicegraph;
set(gca,'fontsize',14)


figure(8)
uS = unique(s);
for ns = 1:length(uS)

    subplot(2,7,ns)
    sel = s == uS(ns);
    
    plot(x(sel),y(sel),'ko')

    xlim([-100,100]);
    ylim([-100,100]);
    axis square
    title(['S' num2str(uS(ns))])
    horline;
    verline;
    unityline;

    xlabel('target (deg)')
    ylabel('response (deg)')
    box off;


    hold on

    vv = floor(linspace(1,1000,80));
    for nv = 1:length(vv)
        plot(-90:90,...
            mean(samples.beta0(vv(nv),ns))+...
            mean(samples.beta1(vv(nv),ns))*(-90:90),...
            'color',[.1,.1,.1,.5],'LineWidth',1)
    end

    plot(-90:90,...
        mean(samples.beta0(:,ns))+...
        mean(samples.beta1(:,ns))*(-90:90),...
        'LineWidth',1.5)

end


%% %%%%%%%%%%%%%%%%%% LOCAL FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function samples = getsamples_jags(x,y,s,type)

% Initialization
chain_globals;
runjagsMethod = 'rjags';
numSavedSteps	= 1e3;
thinSteps		= 1;
burnInSteps		= 1e3;
nChains			= 2;
dic             = false;

% Actual regression
if strcmpi(type,'simple')

    % DATA
    sel = s == 1;
    x       = x(sel);
    y       = y(sel);
    Ntotal  = length(y);

    % Which model to use?
    modelname       = fullfile(pwd, 'Week4_Ex1Mod1.txt');
    parameters		= {'beta0','beta1','sigma'};		% The parameter(s) to be monitored.

    dataStruct      = struct('x',x,'y',y,'Ntotal',Ntotal);

elseif strcmpi(type,'robust')

    sel = s == 1;
    x = x(sel);
    y = y(sel);
    Ntotal  = length(y);

    modelname       = fullfile(pwd, 'Week4_Ex2Mod2.txt');
    parameters		= {'beta0','beta1','sigma','nu'};		% The parameter(s) to be monitored.

    % Specify data, as a structure
    dataStruct = struct('x',x,'y',y,'Ntotal',Ntotal);


elseif strcmpi(type,'hierarchical')

    x       = x;
    y       = y;
    s       = s;
    Ntotal  = length(y);
    Nsubj   = length(unique(s));

    modelname       = fullfile(pwd, 'Week4_Ex3Mod3.txt');
    parameters		= {'beta0','beta1','sigma','nu','beta1mu','beta0mu','zbeta0sigma','zbeta1sigma'};		% The parameter(s) to be monitored.

    % Specify data, as a structure
    dataStruct = struct('x',x,'y',y,'s',s,'Nsubj',Nsubj,'Ntotal',Ntotal);

end



% Chains
initsStruct = struct([]);
for ii = 1:nChains
    if strcmpi(type,'simple')
        % Starting values based on some expectation
        initsStruct(ii).beta0		= 0; % bias starting point
        initsStruct(ii).beta1		= 1; % gain starting point
        initsStruct(ii).sigma		= 1; % standard deviation

    elseif strcmpi(type,'robust')
        % Here data are standardized!
        initsStruct(ii).zbeta0		= 0;
        initsStruct(ii).zbeta1		= corr(x,y);
        initsStruct(ii).zsigma		= 1-corr(x,y)^2;
        initsStruct(ii).nuMinusOne  = 80;

    elseif strcmpi(type,'hierarchical')
        initsStruct(ii).zbeta0		= zeros(Nsubj,1);
        initsStruct(ii).zbeta1		= ones(Nsubj,1)*corr(x,y);
        initsStruct(ii).zsigma		= 1-corr(x,y)^2;
        initsStruct(ii).nuMinusOne  = 80;
    end

end

% run chains
nIter			= ceil((numSavedSteps*thinSteps )/nChains); % Steps per chain.

if strcmp(runjagsMethod,'parallel')
    doparallel		= 1; % do use parallelization
else
    doparallel		= 0; % do not use parallelization
end
fprintf( 'Running JAGS...\n' );
% [samples, stats, structArray] = matjags( ...
[samples, stats] = matjags( ...
    dataStruct, ...                     % Observed data
    modelname, ...    % File that contains model definition
    initsStruct, ...                          % Initial values for latent variables
    'doparallel' , true, ...      % Parallelization flag
    'nchains', nChains,...              % Number of MCMC chains
    'nburnin', burnInSteps,...              % Number of burnin steps
    'nsamples', nIter, ...           % Number of samples to extract
    'thin', thinSteps, ...                      % Thinning parameter
    'dic',dic, ...                       % Do the DIC?
    'monitorparams', parameters, ...     % List of latent variables to monitor
    'savejagsoutput',0, ...          % Save command line output produced by JAGS?
    'verbosity',1, ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup',1);                    % clean up of temporary files?

% Extract chain values:
samples = extractchain(samples); % from multiple -dimension matrix to 1- or 2-D

end