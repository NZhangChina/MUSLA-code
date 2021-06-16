for Lmin_j=1:5
    Lmin_value = floor(0.05*Lmin_j*nD);
    for k_j = 1:5
        temp_S_0{k_j,1}=zeros(k_j * 3, 1+3*Lmin_value*Parameter.num_view );
    end
    for R_j=1
        length = R_j * Lmin_value*Parameter.num_view;
        temp_Segment_matrix = [];
        for vv=1:Parameter.num_view
            temp_Segment_matrix = [temp_Segment_matrix segment_obtain(T{vv,1},length/Parameter.num_view)];
        end
        for k_j =1:4
            [temp_S,~]=my_kmeans(temp_Segment_matrix',k_j);
            temp_S = temp_S';
            temp_S_0{k_j,1}((R_j-1)*k_j+1:R_j*k_j,1)=length*ones(k_j,1);
            temp_S_0{k_j,1}((R_j-1)*k_j+1:R_j*k_j,2:length+1)=temp_S;
        end
    end
    for R_j=1
        for k_j = 1:5
            temp_S0 = temp_S_0{k_j,1}(1:k_j*R_j,1:1+R_j*Lmin_value*Parameter.num_view);

            % temp_W_0 = update_W(temp_X_0,temp_Y_0,Parameter);
            % 这里的初始化在后面求出W_0
            S_0_space{Lmin_j,R_j,k_j} = temp_S0;  % initialize S_0;
        end
    end
end
for Lmin_j=6:10
    Lmin_value = floor(0.05*Lmin_j*nD);
    for k_j = 1:5
        temp_S_0{k_j,1}=zeros(k_j * 3, 1+3*Lmin_value*Parameter.num_view );
    end
    for R_j=1:1
        length = R_j * Lmin_value*Parameter.num_view;
        temp_Segment_matrix = [];
        for vv=1:Parameter.num_view
            temp_Segment_matrix = [temp_Segment_matrix segment_obtain(T{vv,1},length/Parameter.num_view)];
        end
        for k_j =1:4
            [temp_S,~]=my_kmeans(temp_Segment_matrix',k_j);
            temp_S = temp_S';
            temp_S_0{k_j,1}((R_j-1)*k_j+1:R_j*k_j,1)=length*ones(k_j,1);
            temp_S_0{k_j,1}((R_j-1)*k_j+1:R_j*k_j,2:length+1)=temp_S;
        end
    end
    for R_j=1:1
        for k_j = 1:5
            temp_S0 = temp_S_0{k_j,1}(1:k_j*R_j,1:1+R_j*Lmin_value*Parameter.num_view);

            % temp_W_0 = update_W(temp_X_0,temp_Y_0,Parameter);
            % 这里的初始化在后面求出W_0
            S_0_space{Lmin_j,R_j,k_j} = temp_S0;  % initialize S_0;
        end
    end
end
% save(['out/' file '/' file '_initialization1_all_views'], 'S_0_space')