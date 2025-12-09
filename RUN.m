%% Run this file section by section to get an overview over the outcome of this project
%% Load a given Setup

%Initialize
MM = Macro_model();

%EITHER Load demo setup
MM.load_setup('Data\config_Reg25Cyl1.45id5308.json');

%OR Select Setup file of json format
%MM.load_setup(0);

%Load observed data
MM.load_data('Data\positionsVelReg25Cyl1.45id5308.csv');

%OR Select suitable file of csv format
%MM.load_data(0);

%Visualize
MM.plot_obstacle();

%% Visualize Data
MM.plot_data();

%% Example Computations

%Computation with standard kernel model
MM.compute();
MM.plot()

%% Pretrained Nets

%Load nets
load('TrainedNets.mat')

%% Compute and Plot Trained Convolution model
NT = MM.NT;
MM.compute_Macro_NN(trainedNetConvolution,NT);
MM.plot_net();

%% Compute and Plot Trained Nonlinear Neural network model
NT = MM.NT;
MM.compute_Macro_NN(trainedNetNN,NT);
MM.plot_net();


%% Sampling and Neural-Network Training
% Will train new NetworkstrainedNetConvolution and trainedNetNN

%Sample datapoints from the flow data of different geometries
DataCollection

% Train Networks based on these samples
Machine_Learning




