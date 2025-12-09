%% Definition of Networks

% Network for Kernel Method

layersConv = layerGraph();
l_temp = [
    imageInputLayer([21 21 2],"Name","imageinput")
    ];
layersConv = addLayers(layersConv,l_temp);

l_temp = [
    convolution2dLayer([21 21],1,"Name","conv1",'WeightsInitializer','zeros','BiasLearnRateFactor',0)
    ];
layersConv = addLayers(layersConv,l_temp);

l_temp = [
    convolution2dLayer([21 21],1,"Name","conv2",'WeightsInitializer','zeros','BiasLearnRateFactor',0)
    ];
layersConv = addLayers(layersConv,l_temp);

l_temp = [
    depthConcatenationLayer(2,"Name","depthcat")
    regressionLayer("Name","regressionoutput")];
layersConv = addLayers(layersConv,l_temp);

layersConv = connectLayers(layersConv,"imageinput","conv1");
layersConv = connectLayers(layersConv,"imageinput","conv2");
layersConv = connectLayers(layersConv,"conv1","depthcat/in1");
layersConv = connectLayers(layersConv,"conv2","depthcat/in2");

% Network for Nonlinear Method
layersNN = [
    imageInputLayer([21 21 2],"Name","imageinput")
    averagePooling2dLayer([3 3],"Name","avgpool2d","Stride",[3 3])
    %convolution2dLayer([7 7],20,"Name","conv","Stride",[7 7])
    % fullyConnectedLayer(32,"Name","fc")
    reluLayer("Name","relu1")
    fullyConnectedLayer(64,"Name","fc1")
    reluLayer("Name","relu2")
    fullyConnectedLayer(32,"Name","fc2")
    reluLayer("Name","relu3")
    fullyConnectedLayer(16,"Name","fc_1")
    reluLayer("Name","relu4")
    fullyConnectedLayer(8,"Name","fc_11")
    reluLayer("Name","relu5")
    fullyConnectedLayer(4,"Name","fc_12")
    reluLayer("Name","relu6")
    fullyConnectedLayer(2,"Name","fc_13")
    regressionLayer("Name","regressionoutput")];

%% Assemble data



N = length(samples);

%Input into Network: "Image" of density
Images = [];
Images2 = [];

%Output: Velocity in x and y dimension
Labels = [];
Labels2 = [];

j = 0;

for i = 1:N
    vel = samples{i}.velocity(1);
    if abs(vel) > 1.2
        j = j+1;
        continue %Velocity above 1 corresponds to unphysical behaviour or measurement-errors.
    end
    if abs(vel) < 0.05
        if rand < 0.9
            j = j+1;
            continue %reduce number of samples with velocity 0.
        end
    end
    Images(:,:,1,i-j) =  samples{i}.density;
    Images2(:,:,1,i-j) =  flipud(samples{i}.density); %mirror by x Axis
    Labels(1,1,:,i-j) = [samples{i}.velocity(1),samples{i}.velocity(2)];
    Labels2(1,1,:,i-j) = [-samples{i}.velocity(1),samples{i}.velocity(2)]; %x-Velocity is not influenced by mirroring, y_velocity reversed
    Images(:,:,2,i-j) =  samples{i}.obstacle;
    Images2(:,:,2,i-j) =  flipud(samples{i}.obstacle);
end
Images = cat(4,Images, Images2);
% Obstacles = cat(3,Obstacles, Obstacles2);
Labels = cat(4,Labels, Labels2);
clear Labels2 Images2 N j vel i



%% Define options for the training
opts = trainingOptions("sgdm", ...
    'InitialLearnRate',0.00001,...   %0.000001 works fine for layer1
    'MaxEpochs',6, ...
    'Shuffle','every-epoch', ...
    'Plots','training-progress', ...
    'MiniBatchSize',128, ...
    'L2Regularization',0.0000001, ...
    'Verbose',true,...
    'LearnRateDropPeriod',1,...
    'GradientThreshold',inf,...
    'LearnRateSchedule', "piecewise");
%'InitialLearnRate',0.01,...

%% Training of Kernel-network
opts.InitialLearnRate = 0.00001;
[trainedNetConvolution, info] = trainNetwork(Images,Labels,layersConv, opts); %layerGraph(trainedNet2)

%% Training of Nonlinear-network
opts.InitialLearnRate = 0.05; %For different data one might have to fine-tune this and the other options
[trainedNetNN, info] = trainNetwork(Images,Labels,layersNN, opts);

%% Output
%Select net for which to show the Output
trainedNet = trainedNetNN;

factor = 2;
YTest = predict(trainedNet, Images(:,:,:,1:factor:end));
YTrue = squeeze(Labels(:,:,:,1:factor:end))';
%% Compare the expected and the predicted solution of the nonlinear function
compare_x = [YTest(:,1) YTrue(:,1)];
compare_y = [YTest(:,2) YTrue(:,2)];
% figure(17)
% scatter(YTest, YTrue);

figure(18)
tiledlayout(1,2)
nexttile
[values, centers] = hist3(compare_x,[51 51]);
imagesc(centers{:}, values.')
hold on
xlabel('Prediction')
ylabel('Measurement')
colorbar
plot(centers{1},centers{1},'r')
hold off
nexttile
[values, centers] = hist3(compare_y,[51 51]);
imagesc(centers{:}, values.')
hold on
xlabel('Prediction')
ylabel('Measurement')
colorbar
plot(centers{1},centers{1},'r')
hold off
%% Errorbar Plot
velocities = -1:0.05:1;
binEdges = conv(velocities, [0.5, 0.5], 'valid');
velocities = velocities(2:end-1);
for xy = 1:2
    YTrue_ = YTrue(:,xy);
    YTest_ = YTest(:,xy);
    [N,edges,bin] = histcounts(YTrue_,binEdges);

    N_bins = length(N);
    variances = zeros(1,N_bins);
    means = zeros(1,N_bins);
    for i = 1:length(N)
        predictionerror = YTest_(bin == i)-YTrue_(bin == i);
        variances(i) = var(predictionerror);
        quantiles(i,:) = quantile(predictionerror,[0.2,0.8]);
        means(i) = median(predictionerror);
        binsize(xy,i) = length(predictionerror);
    end

    fig_errorbar = figure(xy+5);
    errorbar(velocities,velocities+means,quantiles(:,1)-means',quantiles(:,2)-means')
    title('Estimated vs True velocity')
    hold on
    %errorbar(velocities(1:end-1),velocities(1:end-1)+means,variances.^.5)
    plot(velocities,velocities,'r')
    switch xy
        case 1
            xlim([-0.5,0.5])
            ylim([-0.5,0.5])     
        case 2
            xlim([-1.1,0])
            ylim([-1.1,0])
    end

    legend('Estimation: Mean and 10%/90% Quantiles','True value', 'Location','northwest')
    xlabel('Observed real-world velocity of samples')
    ylabel('Estimation')
    grid on
    hold off
end


%%
trainedNet = trainedNetConvolution;
figure(19)
title('Estimated kernels')
tiledlayout(2,3)
nexttile
imagesc(trainedNet.Layers(2, 1).Weights(:,:,1));
title('x-Direction, Density')
nexttile
imagesc(trainedNet.Layers(2, 1).Weights(:,:,2));
title('x-Direction, Obstacle')
nexttile
imagesc(trainedNet.Layers(2, 1).Weights(:,:,1)+trainedNet.Layers(2, 1).Weights(:,:,2));
title('x-Direction, Combined')
nexttile
imagesc(trainedNet.Layers(3, 1).Weights(:,:,1));
title('y-Direction, Density')
nexttile
imagesc(trainedNet.Layers(3, 1).Weights(:,:,2));
title('y-Direction, Obstacle')
nexttile
imagesc(trainedNet.Layers(3, 1).Weights(:,:,1)+trainedNet.Layers(3, 1).Weights(:,:,2));
title('y-Direction, Combined')
colorbar


%% *Output Macromodel*
% 
% YTrue = cell2mat(Labels)';
% YComp = cell2mat(MM_computed)'
% compare2 = [YComp YTrue];
% figure(20)
% [values, centers] = hist3(compare2,[51 51]);
% imagesc(centers{:}, values.')
% colorbar