clear;
addpath(genpath('../functions'))
addpath(genpath('../datasets'))

%% load data
file = 'Libras';

Parameter.num_view = 2
for temp_num_view = 1:Parameter.num_view
    path_train = ['datasets/' file '/' file 'Dimension' num2str(temp_num_view) '_TRAIN.arff']
    Data_train = read_file(path_train);
    N_train = size(Data_train,1);
    path_test = ['datasets/' file '/' file 'Dimension' num2str(temp_num_view) '_TEST.arff']
    Data_test = read_file(path_test);
    max_value = max(max([Data_train(:,1:end-1);Data_test(:,1:end-1)]));
    min_value = min(min([Data_train(:,1:end-1);Data_test(:,1:end-1)]));
    Data = [Data_train; Data_test];
    [mD,nD]=size(Data);
    Y_true=Data(:,end);
%     Y_true=Data(N_train+1:end,end); % test
    DT=Data(:,1:end-1);
    DT=(DT-min_value)/(max_value-min_value); % regularize time series data;
    T{temp_num_view,1}=[(nD-1)*ones(mD,1),DT]; % Each element in first column of time series matrix is the length of the time series in that row;
end
N=mD; % Number of time series

%% Parameters initialization
clear length
% fixed
Parameter.C=length(unique(Y_true));% the number of clusters
Parameter.alpha=-100; % parameter in Soft Minimum Function
Parameter.routing_iter_num = 3;
S_tp1_init_regu = 1; % S_tp1_init_regu = 1 on libras otherwise = 0
Parameter.Imax=50; % the number of internal iterations
Parameter.eta=0.01; % learning rate
Parameter.epsilon=0.1; % internal convergence parameter
if N/Parameter.C < 20
    Parameter.k_neigh = 5;  % number of neighbors to determine the initial graph of A
else
    Parameter.k_neigh = 15;  % number of neighbors to determine the initial graph of A
end
Parameter.lambda_0_determined = 0;
Parameter.islocal = 1; % islocalµÄÖµÎª1: only update the similarities of the k neighbor pairs, faster; ÖµÎª0: update all the similarities
Num_iter = 30;
Parameter.sigma1 = 1; % parameter in spectral kernel
Parameter.sigma2 = 1; % parameter in shapelets' similarity
Max_iter = 15; % 15 30
Parameter.k=5;% the number of shapelets in equal length 
Parameter.Wv_determined = 0;
Parameter.Y_is_initialization = 1;
Parameter.Lmin_space = [0.05, 0.01, 0.15, 0.2, 0.25,0.3, 0.35, 0.4, 0.45, 0.5];
Parameter.R_space = [1,2,3];
Parameter.lambda_space = [10^-4,10^-2,10^0,10^2,10^4];
% selected
Parameter.fixed_attetions = 0; % 1 0.5 0
Lmin_j=3; Parameter.Lmin=floor(0.05*Lmin_j*nD);
Parameter.R=3;
Parameter.lambda_0=Parameter.lambda_space(5);
for Rj_i=1:Parameter.R
    Parameter.lambda_1(Rj_i) = Parameter.lambda_space(2);
end
Parameter.lambda_2=Parameter.lambda_space(1);

%% run MUSLA

% load initialization_shapelets
load([ 'out/'  file '/' file '_initialization1_all_views.mat']) % load initialization_shapelets
% initialization_shapelets % or run initialization_shapelets.m
S_tp0=S_0_space{Lmin_j,Parameter.R,Parameter.k};
if S_tp1_init_regu == 1
    S_tp0 = [S_tp0(:,1), z_regularization(S_tp0(:,2:end))]; 
end
for Rj_i=1:Parameter.R
    S_tp1{Rj_i,1} = S_tp0(Rj_i*Parameter.k-Parameter.k+1:Rj_i*Parameter.k,1:S_tp0(Rj_i*Parameter.k,1)+1);
end

gap=100; F_tp1=10000; F_t=F_tp1+10^10; wh_time=0;
while gap>Parameter.epsilon && wh_time<Max_iter
    %%%%%%Calculation matrix
    for Rj_i=1:Parameter.R
        Parameter.Att_views_b = zeros(Parameter.num_view, Parameter.k); 
        Parameter.Att_views_fixed = Parameter.fixed_attetions*ones(Parameter.num_view, Parameter.k)/Parameter.num_view;
        [X_tp1{Rj_i,1},Xkj_tp1_skl{Rj_i,1},Parameter] = distance_timeseries_Cap_shapelets(T,S_tp1{Rj_i,1},Parameter);
    end
    % update L_A_tp1;
    [L_A_tp1,A_tp1,temp_A0_tp1,Y_tp1,distX_tp1,Wv,Parameter] = update_MUSLA_AY(X_tp1,Parameter,Num_iter); % update L_G_tp1, the Laplacian matrix of similarity matrix G of time series. 
    for Rj_i=1:Parameter.R
        [SS_tp1{Rj_i,1},XS_tp1{Rj_i,1},SSij_tp1_sil{Rj_i,1},Parameter] = shapelet_similarity(S_tp1{Rj_i,1},Parameter); % update shapelets similarity matrix H_tp1;
    end
    F_tp1 = function_value_MUSLA_Simplify(X_tp1,distX_tp1,Y_tp1,L_A_tp1,A_tp1,SS_tp1,Parameter); % calculate the value of objective function;
    gap = F_t-F_tp1; 

    if isnan(F_tp1)
        wh_time = Max_iter;
    else
        S_tp1 = update_MUSLA_S(X_tp1,A_tp1,temp_A0_tp1,S_tp1,Xkj_tp1_skl,SSij_tp1_sil,SS_tp1,Wv,Parameter); % update S_tp1;
        F_t=F_tp1;
        wh_time=wh_time+1;
    end
end
[L_A_tp1,A_tp1,temp_A0_tp1,Y_tp1,distX_tp1, Wv,Parameter] = update_MUSLA_AY(X_tp1,Parameter,Num_iter);% update L_G_tp1, the Laplacian matrix of similarity matrix G of time series. 

[RI_test_combined,RI_test_alone,NMI_test_combined,NMI_test_alone] = compute_cluster_index(Y_true,Parameter,A_tp1,N_train)
