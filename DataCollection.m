%% From a velocity field and a densitiy field over time extract snapshots ofthe flow dependent on the surroundings

%% Select Files
%[files_data,path_data] = uigetfile('*.csv', 'Select Multiple Files', 'MultiSelect', 'on' );
%[files_setup,path_setup] = uigetfile('*.json', 'Select Multiple Files', 'MultiSelect', 'on' );

%% Preselected Files (alternative to above)
path_data = 'Data\';
path_setup = 'Data\';
load('Data\ListOfFiles.mat');

%% Read Data

N = length(files_data);
data = cell(1,N);

% Using parallel computing for speed. Alternatively replace parfor with for 
parfor i = 1:N
    
    MM_Data = Macro_model();
    MM_Data.load_setup(fullfile(path_setup,files_setup{i}));
    
    MM_Data.load_data(fullfile(path_data,files_data{i}));
    
    data{i} = MM_Data.extract_fluxdata();
    
    %data{i} = MM.compute_flux(data{i});
    fprintf('Iteration Number: = %g\n', i);
    fprintf('Number of Datapoints sampled: = %g\n', length(data{i}));
    
end

%% Save
%save('Data\selection',"data","files_data","-v7.3","-nocompression")

%% Concatenate
samples = {};
N = length(data);
for i = 1:N
    r = ceil(length(data{i})/10000);%If too many samples per experiment, reduce.
    samples = [samples data{i}(1:r:end)];
end

%% Load 

% for i = 1:length(data)
%     for j = 1:length(data{i})
%         %data{i}{j}.density = single(data{i}{j}.density);
%         %data{i}{j}.obstacle = boolean(data{i}{j}.obstacle);
%         data{i}{j}.velocity = data{i}{j}.velocity/0.1; 
%     end
%     i
% end

%% Plot histogram of the velocities (relative to the belt velocity)
% Note: Spikes at 0 and -1
velocities = zeros(2,length(samples));

for i = 1:length(samples)
    velocities(:,i) = samples{i}.velocity;
end

%%
figure(6)
hist(velocities(1,:),100)
title('Velocity histogram perpendicular to flow direction')
% Spike at 0 (no change), mostly positive values <1 due to regulators at
% the right side of the belt

figure(7)
hist(velocities(2,:),100)
title('Velocity histogram in flow direction')
% Spike at 0 (no change) and at -1 (full blockage), all values in between
% are obtained.
