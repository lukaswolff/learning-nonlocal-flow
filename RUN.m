%% Load a given Setup

%Initialize
MM = Macro_model();

%Load demo setup
MM.load_setup('Data\config_Reg30Cyl1.45id5398.json');

%OR Select Setup file of json format
%MM.load_setup(0);

%Load observed data
MM.load_data('Data\positionsVelReg30Cyl1.45id5398.csv');

%OR Select suitable file of csv format
MM.load_data(0);

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
NT = MM.NT;
MM.compute_Macro_NN(trainedNetConvolution,NT);
MM.plot_net();

%%To be continued


