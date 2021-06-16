function Data_z=z_normlization(Data)
% 每行一个时间序列
D_max=repmat(max(Data,[],2),1,size(Data,2));
D_min=repmat(min(Data,[],2),1,size(Data,2));
Data_z=(Data-D_min+10^-10)./(D_max-D_min+10^-10);