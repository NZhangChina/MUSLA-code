function [D_T_S, Xkj_skl,Parameter]=distance_timeseries_Cap_shapelets(T,S,Parameter)

alpha = Parameter.alpha;
num_view = Parameter.num_view;
% Parameter.fixed_attetions 
[mT,nT]=size(T{1,1});
for vv=1:num_view
    DT{vv,1}=T{vv,1}(:,2:end);
end
[mS0,nS]=size(S);
DS=S(:,2:end);

Xkj_skl=zeros(mS0,mT,nS-1);
for routing_iter_ii = 1:Parameter.routing_iter_num
    Parameter.Att_views_c = (1-Parameter.fixed_attetions)*exp(Parameter.Att_views_b) ./ repmat(sum(exp(Parameter.Att_views_b),1),Parameter.num_view,1);
    Parameter.Att_views_c =  Parameter.Att_views_c + Parameter.Att_views_fixed;
    
    for j=1:mT % the j-th time series
        for k=1:mS0
            len=S(k,1);% the length of k-th shapelet
            shapelet=DS(k,1:len); % the k-th shapelet
            Q_j=T{1,1}(j,1); %the length of time series
            for vv=1:num_view
                time_series(vv,:)=DT{vv,1}(j,1:Q_j); % the j-th time series
            end
            % distance_shapelet_views_style = 1
           
            series_long = time_series;
            series_short = shapelet;
            series_shorts = [];
            [num_view mL]=size(series_long); % num_view,mL
            mS=length(series_short);
            for vv = 1:num_view
                series_shorts(vv,:) = series_short((vv-1)*mS/num_view+1:vv*mS/num_view);
            end
            num_seg=mL-mS/num_view+1;
            for vv = 1:num_view
                D1{vv,1}=[];
            end
            D2_views =[];
            D2 = [];
            for q=1:num_seg
                segs=series_long(:,q:q+mS/num_view-1);% the q-th segment of the long series
                seg = [];
                seg_short =[];
                D2_views(q,1) = 0;
                for vv = 1:num_view
                    seg = segs(vv,:);
                    seg_short = series_shorts(vv,:);

                    D2{vv,1}(q,1)=num_view/mS*norm(seg_short-seg)^2;
                    D2_views(q,1) = D2_views(q,1) + D2{vv,1}(q,1) * Parameter.Att_views_c(vv,k);
                    D1_p=num_view/mS*(seg_short-seg);
                    D1{vv,1}=[D1{vv,1};D1_p];
                end
            end

    %         X_forward = [];
            for vv = 1:num_view
                X_up{vv,1} = sum(D2{vv,1}.*exp(alpha*D2_views));
                X_down{vv,1} = sum(exp(alpha*D2_views));
                X_forward(j,k,vv) = X_up{vv,1}/X_down{vv,1}; % the distance between the series_long and series_short;
            end
        end
    end
    temp_Att_views_c(1,:,:) = Parameter.Att_views_c';
%     X_forward_SqNorm = repmat(sum(X_forward.^2,1),mT,1,1);
    X_s = sum(repmat(temp_Att_views_c,mT,1,1) .* X_forward,3); 
    X_s_SqNorm = repmat(sum(X_s.^2,1),mT,1);
    X_v = X_s .* sqrt(X_s_SqNorm) ./ (1 + X_s_SqNorm);
%     Att_views_b_add(:,:) = sum(X_forward .* repmat(X_v,1,1,vv) ./ X_forward_SqNorm,1);
    Att_views_b_add(:,:) = sum(X_forward .* repmat(X_v,1,1,vv),1);
    Parameter.Att_views_b = Parameter.Att_views_b + Att_views_b_add';
end


Parameter.Att_views_c = (1-Parameter.fixed_attetions)*exp(Parameter.Att_views_b) ./ repmat(sum(exp(Parameter.Att_views_b),1),Parameter.num_view,1);
Parameter.Att_views_c =  Parameter.Att_views_c + Parameter.Att_views_fixed;
X = [];
X_up = [];
X_down = [];
for j=1:mT % the j-th time series
    for k=1:mS0
        len=S(k,1);% the length of k-th shapelet
        shapelet=DS(k,1:len); % the k-th shapelet
        Q_j=T{1,1}(j,1); %the length of time series
        for vv=1:num_view
            time_series(vv,:)=DT{vv,1}(j,1:Q_j); % the j-th time series
        end
        % distance_shapelet_views_style = 1

        series_long = time_series;
        series_short = shapelet;
        series_shorts = [];
        [num_view mL]=size(series_long); % num_view,mL
        mS=length(series_short);
        for vv = 1:num_view
            series_shorts(vv,:) = series_short((vv-1)*mS/num_view+1:vv*mS/num_view);
        end
        num_seg=mL-mS/num_view+1;
%         for vv = 1:num_view
%             D1{vv,j,k}=[];
%         end
        D1=[];
        D2_views =[];
        for q=1:num_seg
            segs=series_long(:,q:q+mS/num_view-1);% the q-th segment of the long series
            seg = [];
            seg_short =[];
            D2_views(q,1) = 0;
            D1_p=[];
            for vv = 1:num_view
                seg = segs(vv,:);
                seg_short = series_shorts(vv,:);
                D2{vv,1}(q,1)=num_view/mS*norm(seg_short-seg)^2;
                D2_views(q,1) = D2_views(q,1) + D2{vv,1}(q,1) * Parameter.Att_views_c(vv,k);
                D1_p_vv=Parameter.Att_views_c(vv,k)*num_view/mS*(seg_short-seg);
                D1_p=[D1_p D1_p_vv];
%                 D1{vv,j,k}=[D1{vv,j,k};D1_p];
            end
            D1 = [D1; D1_p];
        end

%         X(j,k) = 0;
%         for vv = 1:num_view
%             X_up{vv,j,k} = sum(D2{vv,j,k}.*exp(alpha*D2_views{j,k}));
%             X_down{vv,j,k} = sum(exp(alpha*D2_views{j,k}));
%             X(j,k) = X(j,k) + Parameter.Att_views_c(vv,k)*X_up{vv,j,k}/X_down{vv,j,k}; % the distance between the series_long and series_short;
%         end
        X_up=sum( D2_views .* exp(alpha*D2_views) );
        X_down=sum(exp(alpha*D2_views));
        X(j,k)=X_up/X_down;
        
        Xkj_sk =[];
        len=S(k,1);
        for l=1:len 
            part1=1/X_down^2; 
            part2=D1(:,l);
            part3=exp(alpha*D2_views).*(X_down*(1+alpha*D2_views)-alpha*X_up);
            Xkj_sk(l)=part1*sum(part2.*part3); % the derivative of X_(kj) on S_(kl)
        end
        Xkj_skl(k,j,1:len)=Xkj_sk;
%         for vv = 1:num_view
%             for l=1:len/num_view
%                 ll = ll + 1;
%                 part1=1/X_down{vv,j,k}^2;
%                 part2=D1{vv,j,k}(:,l);
%                 part3=exp(alpha*D2_views{j,k}).*(X_down{vv,j,k}*(1+alpha*D2_views{j,k})-alpha*X_up{vv,j,k});
%                 Xkj_sk(ll)=part1*sum(part2.*part3); % the derivative of X_(kj) on S_(kl)
%             end
%         end
%         Xkj_skl(k,j,1:len)=Xkj_sk;
    end
end

D_T_S=X';
            
