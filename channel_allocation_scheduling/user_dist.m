function [X_pu ,...
    Y_pu,...
    X_su_bs,...
    Y_su_bs,...
    COMP_compatiblity_su,...
    cell_association_state,...
    p_rx_su_bs,...
    p_rx_su_pu,...
    p_rx_su_bs_COMP,...
    p_rx_bs_pu] =user_dist(X_max,...
    Y_max,...
    N_pu,...
    N_su,...
    N_bs,...
    X_org,...
    Y_org,...
    R2,...
    p_tx_cbs,...
    p_tx_pu,...
    d0,...
    n)

count_su=0;
N_temp_su=N_su*10; % we will take N_su SUs from a bigger number of SUs


%% pu distribution

X_pu =X_max*rand(1,N_pu);
Y_pu =Y_max*rand(1,N_pu);
%% su distribution
X_su =X_max*rand(1,N_temp_su);% available SU users
Y_su =Y_max*rand(1,N_temp_su);


%% Associating SU with BS
X_su_bs=cell(N_bs,1);% SU's within coverage of BS
Y_su_bs=cell(N_bs,1);

for i=1:N_temp_su
    
    if count_su==N_su
        break
    end
    dist_su_bs=sqrt((X_org-X_su(i)).^2+(Y_org-Y_su(i)).^2);
    [temp,ind_su_bs]=min(dist_su_bs);
    
    if temp<=R2
        
        X_su_bs{ind_su_bs}=[X_su_bs{ind_su_bs} X_su(i)];
        Y_su_bs{ind_su_bs}=[Y_su_bs{ind_su_bs} Y_su(i)];
        count_su=count_su+1; 
    end
end

for i=1:N_bs
    N_su_bs(i)=length(X_su_bs{i});
end

%% transforming X_su_bs cell array to a vector
X_temp=X_su_bs;
Y_temp=Y_su_bs;
X_su_bs=[];
Y_su_bs=[];
cell_association_state=[];
p_rx_su_bs=[];
p_rx_su_pu=[];
for i=1:N_bs
    X_su_bs=[ X_su_bs X_temp{i}];
    Y_su_bs=[ Y_su_bs Y_temp{i}];
    cell_association_state=[cell_association_state i*ones(1,length(X_temp{i}))];
    
    dist_su_bs=(sqrt((X_org(i)-X_temp{i}).^2+(Y_org(i)-Y_temp{i}).^2));
    p_rx_su_bs=[p_rx_su_bs ...
        p_tx_cbs*(dist_su_bs/d0).^(-1*n)];
    % received power at SU from BS
    X_temp_temp=X_temp{i};
    Y_temp_temp=Y_temp{i};
    for j=1:length(X_temp{i})
        dist_su_pu=sqrt((X_pu-X_temp_temp(j)).^2+(Y_pu-Y_temp_temp(j)).^2);% distance vector
        ... between a particular SU and all PU
            p_rx_su_pu=[p_rx_su_pu ; p_tx_pu*(dist_su_pu/d0).^(-1*n)];% received power vector
        ...at SU from all PU
    end
end

%% finding COMP compatible SU's (for 2 BS only)
COMP_compatiblity_su=zeros(1,length(X_su_bs));
p_rx_su_bs_COMP=zeros(1,length(X_su_bs));
for i=1:length(X_su_bs)
    temp1=sqrt((X_su_bs(i)-X_org(1)).^2+(Y_su_bs(i)-Y_org(1)).^2);
    temp2=sqrt((X_su_bs(i)-X_org(2)).^2+(Y_su_bs(i)-Y_org(2)).^2);
    
    if temp1<R2&&temp2<R2
        COMP_compatiblity_su(i)=1;
        
    end
    if cell_association_state(i)==1
        p_rx_su_bs_COMP(i)=p_tx_cbs*(temp2/d0).^(-1*n);
    else
        p_rx_su_bs_COMP(i)=p_tx_cbs*(temp1/d0).^(-1*n);
        
    end
end
p_rx_bs_pu=zeros(N_bs,N_pu);
for i=1:length(X_pu)
    temp1=sqrt((X_pu(i)-X_org(1)).^2+(Y_pu(i)-Y_org(1)).^2);
    temp2=sqrt((X_pu(i)-X_org(2)).^2+(X_pu(i)-Y_org(2)).^2);
    
    p_rx_bs_pu(1,i)=p_tx_pu*(temp1/d0).^(-1*n);
    p_rx_bs_pu(2,i)=p_tx_pu*(temp2/d0).^(-1*n); 
end
