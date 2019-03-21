% Code adapted from Allen et al., "Tracking whole-brain connectivity dynamics in the resting state." Cereb Cortex. 2014 24(3):663-76. 
% Dependency: SimTB toolbox (http://mialab.mrn.org/software/simtb/)
%
% Objective: create a toy simulation with time varying connectivity
% In this example we model 10 components whose connectivity structure varies through 4 discrete states.
% Variables of interest are the state vector ('STATE') and the final TCs ('TC').
% The order and duration (dwell time) of states are controlled by variables 'Sorder' and 'Sdwell'
% The amplitude of unique events (aU) relative to the shared events will  make the modular structure more or less apparent.
% 
%
%Here we simulate 10 "neural" components, to which we add the "BOLD in
%vessels" component, reaching the other regions at different times, and the
%Global Signal, as average of all the above.
%%
clear all; clc
p = which('simtb_create_sP');
if isempty(p)
    INSTALL_DIR = fileparts(which('simtb'));
    addpath(genpath(INSTALL_DIR));
end
%rng(100) % to enable repeated generation of the same simulation
for isub=1:100

% number of components
nC = 10;

% number of time points
nT = 1210;

% TR
TR = .72;

% number of different connectivity states
nStates = 1;

% probability of unique events
pU = 0.5;

% amplitude of unique events (relative to module-specific events)
aU = .5;

% probability of state specific events 
pState = .5;

%Module membership for each state
ModMem = zeros(nC,nStates);

% Number of event types (i.e., number of different modules)
nE = 5;

% Modules are formed by sharing of events.  In each state there are up to
% nE different modules. The matrix ModMem stores the membership of each
% component to a module.Note tjat module 2 in state 2 has nothing to do with 
% module 2 in the other states, it's just an index.
% Negative numbers indicate activity is negatively related to the events in
% the given module.
ModMem(1,:) = [1];
ModMem(2,:) = [2];
ModMem(3,:) = [2];
ModMem(4,:) = [3];
ModMem(5,:) = [3];
ModMem(6,:) = [3];
ModMem(7,:) = [4];
ModMem(8,:) = [4];
ModMem(9,:) = [5];
ModMem(10,:)= [5];


%% Create the event time courses

% random aspects (different for each component)
eT = rand(nT, nC) < pU;
eT = eT.*sign(rand(nT, nC)-0.5);
eT = eT*aU;

% define the order and time in each state
Sorder = [1 ]; % state order
Sdwell = nT;
%NOTE: the Sdwell should sum to nT, check here and amend the last partition:
if sum(Sdwell) ~= nT
    Sdwell(end) = nT - sum(Sdwell(1:end-1));
end
Cdwell = cumsum(Sdwell);
Cdwell = [0 Cdwell];
STATE = zeros(1,nT); % state vector
for ii = 1:length(Sorder)
    sIND = Cdwell(ii)+1:Cdwell(ii+1);
    % events related to each module
    e = rand(length(sIND),nE) < pState;
    e = e.*sign(rand(length(sIND), nE)-0.5);
    for cc = 1:nC
        eT(sIND,cc) = eT(sIND,cc) + sign(ModMem(cc,Sorder(ii)))*e(:,abs(ModMem(cc,Sorder(ii))));
    end
    STATE(sIND) = Sorder(ii);
end

% event time series are stored in eT
%% Convolve event TCs
[tc, MDESC, P, PDESC] = simtb_TCsource(eT(:,1), TR, 1);

P(1) = 6;     % delay of response (relative to onset)
P(2) = 15;    % delay of undershoot (relative to onset)
P(3) = 1;     % dispersion of response
P(4) = 1;     % dispersion of undershoot
P(5) = 3;     % ratio of response to undershoot
P(6) = 0;     % onset (seconds)
P(7) = 32;    % length of kernel (seconds)

PB(1) = 20;     % delay of response (relative to onset)
PB(2) = 30;    % delay of undershoot (relative to onset)
PB(3) = 1;     % dispersion of response
PB(4) = 1;     % dispersion of undershoot
PB(5) = 3;     % ratio of response to undershoot
PB(6) = 0;     % onset (seconds)
PB(7) = 32;    % length of kernel (seconds)

TC  = zeros(nT,nC);
TC(:,1) = simtb_TCsource(eT(:,1), TR, 1, PB); % blood component
for cc = 2:nC
    P(1)=ceil(.5*cc); % here we assign a given increasing fixed delay to each component
    TC(:,cc) = simtb_TCsource(eT(:,cc), TR, 1, P);%, P); % all use same HRF
end

% Add a little gaussian noise
% TC = TC + .1*randn(nT,nC); %0.1
%TC(:,1) = TC(:,1) + 1.1*randn(nT,1); %0.1

blood=TC(:,1);
TC=TC(:,2:10);
TC = TC + .5*randn(size(TC)); %0.1

% the blood signal reaches the sources with increasing delay, so each
% source has its own behavior (but with HRF onset proportional to blood
% arrival time, plus the blood

%%% delay %%%
for i=1:9;TC(1:nT-i,i)=.8*TC(1:nT-i,i)+(2*rand/.2)*blood(i+1:end);end
%%%
%%% no blood %%%
%for i=1:9;TC(1:nT,i)=.5*TC(1:nT,i);end
%%%

TC=TC(1:nT-10,:); %here cut ten points to account for temporal shift
%TCtemp=[TC blood(1:nT-10)];
TCtemp=[10*TC repmat(blood(1:nT-10),1,1005)];
GS=mean(TCtemp,2);
GS=GS+.001*randn(size(GS));
TC=[TC GS]; % add the GS as 10th component
TC=[TC blood(1:nT-10)]; % add blood as 11th component
timeseries=TC;
save(['C:\Users\dmarinaz\Dropbox\code\gsr_iss\simtb\sub' num2str(isub) '.mat'], 'timeseries');
end
plot_results_simtb;