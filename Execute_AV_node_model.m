%% Example code for running the AV node model

%% Create a MEX-file
% The MEX file is generated from the C++ file that Matlab can use

% Uncomment for Mac
% mex -setup C
% mex('AV_Node_model.cpp', '-compatibleArrayDims');
% Uncomment for Windows
% mex AV_Node_model.cpp -largeArrayDims

%%
clear all;
close all;
clc;

%% Generate AA series
P4D_mu = 150; % Mean of Pearson Type IV Distribution
P4D_sigma = 10; % Standard Deviation of Pearson Type IV Distribution
P4D_gamma = 1; % Skewness of Pearson Type IV Distribution
P4D_kappa = 6; % Kurtosis of Pearson Type IV Distribution

L_aa = 10000; % Length of series with atrial arrival times

impulses = pearsrnd(P4D_mu,P4D_sigma,P4D_gamma,P4D_kappa,L_aa-1,1);
aa_impulses = [0;cumsum(impulses)];

%% Set AV node model parameters
R_FP = [300 400 250]; % Refractory period of fast pathway
R_SP = [200 300 250]; % Refractory period of slow pathway
D_FP = [5 7 250]; % Conduction delay of fast pathway
D_SP = [15 7 250]; % Conduction delay of slow pathway

% The refractory period and conduction delay is hard-coded in the model file and would have
% following parameters
%     R_CN = [250 0 1]; % Refractory period of coupling node
%     D_CN = [0 0 1]; % Conduction delay of coupling node

A_ampl = 0.2; % peak-to-peak amplitude of the respiratory modulation
A_freq = 0.2; % respiration frequency in Hz
A = [A_ampl,A_freq]; % time-varying ANS factor that affects the refractory period and conduction delay

%% Run AV Node Model
[hisEntries,AHPathway,AHRT,AHDL,AHET,~,~,~] = AV_Node_model(...
    aa_impulses,R_FP,R_SP,D_FP,D_SP,A);

%% Calculate the RR series from the ventricular activation times
% HisEntries has the same amount of rows as aa_impulses, as the ventricular
% activations are tracked for each atrial impulse individually.
% HisEntries has two columns, as the ventricular activations are tracked
% for the slow and fast pathway individually. We do this because an atrial
% activation can cause two ventricular activations, one over the slow 
% pathway and one over the fast pathway for some parameter model sets.
%
% First, all the NaNs from impulses that didn't result in ventricular
% activations are excluded.
% Then, the resulting array is sorted chronologically
% Finally, the RR series is computed from the ventricular activation times.
rr_seq = diff(sort(hisEntries(~isnan(hisEntries))));

% AHPathway:
% For each atrial impulse, this array stores whether this impulse resulted
% in a ventricular excitation and over which pathway the impulse traveled.
% 0: Concealed conduction
% 1: The atrial impulse resulted in a ventricular excitation and travelled over the slow pathway
% 2: The atrial impulse resulted in a ventricular excitation and travelled over the fast pathway
% 3: The atrial impulse resulted in two ventricular excitation, where one excitation resulted from the impulse travelling over the slow pathway and the other from an impulse from the slow pathway

% AHRT:
% For each atrial impulse that travels from node to node, the refractory
% periods that are computed are stored.

% AHDL:
% For each atrial impulse that travels from node to node, the conduction
% delays that are computed are stored.

% AHET:
% For each atrial impulse that travels from node to node, the time at which
% the impulse arrives at the nodes are stored.

%% Plottinf the simulated RR interval series in a histogram
histogram(rr_seq, 200:50:2000)
xlabel('RR (ms)')
