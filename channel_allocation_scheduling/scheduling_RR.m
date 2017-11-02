function [bits]=...
    scheduling_RR(cell_association_state,...
    F_i_list,...
    F_j_list,...
    F_ij_list,...
    F_ij_JT_list,...
    N_su,...
    N_ch,...
    N_slot,...
    T_slot,...
    bits,...
    p_rx_su_bs,...
    p_n,...
    B)


su_timeslot=zeros(N_su,N_slot);
...transmit only once at each time sample
    % bs_state=zeros(1,N_bs) ;% variable to switch between BS and transmit to a SU
channel_assigned_su=cell(1,N_su);
ch_state=zeros(N_ch,N_slot);
su_state=ones(1,N_su);
count_su=zeros(1,N_su);

%% allocating channels to i and j th BSs

F_i_list=[F_i_list F_ij_list(1:ceil(length(F_ij_list)/2))...
    F_ij_JT_list(1:floor(ceil(length(F_ij_JT_list)/2)))]

F_j_list=[F_j_list F_ij_list(ceil(length(F_ij_list)/2+1):end)...
    F_ij_JT_list(ceil(length(F_ij_JT_list)/2)+1:end)]

%% updating user queue after each slot assignment
% ind_su contains all user's index
[dummy,ind_su] = sort(count_su);

rem_slot=(length(F_ij_JT_list)+length(F_ij_list)+length(F_i_list)+length(F_j_list))*N_slot;

while rem_slot>0 && sum(su_state)>0
    rem_slot;
    rem_slot_previous=rem_slot;
    %% get the index of user that is currently being scheduled
    for i=1:length(ind_su)
        if su_state(ind_su(i))==1;
            scheduling_su_ind=ind_su(i);
            break
        end
    end
    
    scheduled=0;
    
    
    %% First looking at F_i/F_j depending on cell association
    if cell_association_state(scheduling_su_ind)==1
        for i_i=1:length(F_i_list)
            for i_slot=1:N_slot
                if  su_timeslot(scheduling_su_ind,i_slot)==0 && ...
                        ch_state(F_i_list(i_i),i_slot)==0
                    
                    channel_assigned_su{scheduling_su_ind}=...
                        [channel_assigned_su{scheduling_su_ind} F_i_list(i_i) ];
                    su_timeslot(scheduling_su_ind,i_slot)=1;
                    ch_state(F_i_list(i_i),i_slot)=1;
                    rem_slot=rem_slot-1;
                    %% updating bits count
                    p_temp=p_rx_su_bs(scheduling_su_ind);
                    C=B*log2(1+p_temp/p_n);
                    bits(scheduling_su_ind)=bits(scheduling_su_ind)+C*T_slot;
                    count_su(scheduling_su_ind)=count_su(scheduling_su_ind)+1;
                    %% updating user queue after each slot assignment
                    [dummy,ind_su] = sort(count_su);
                    
                    scheduled=1;
                    break
                end
                if scheduled
                    break
                end
            end
            if scheduled
                break
            end
        end
    else
        for i_j=1:length(F_j_list)
            for i_slot=1:N_slot
                if  su_timeslot(scheduling_su_ind,i_slot)==0 && ...
                        ch_state(F_j_list(i_j),i_slot)==0
                    
                    channel_assigned_su{scheduling_su_ind}=...
                        [channel_assigned_su{scheduling_su_ind} F_j_list(i_j) ];
                    su_timeslot(scheduling_su_ind,i_slot)=1;
                    ch_state(F_j_list(i_j),i_slot)=1;
                    rem_slot=rem_slot-1;
                    %% updating bits count
                    p_temp=p_rx_su_bs(scheduling_su_ind);
                    C=B*log2(1+p_temp/p_n);
                    bits(scheduling_su_ind)=bits(scheduling_su_ind)+C*T_slot;
                    count_su(scheduling_su_ind)=count_su(scheduling_su_ind)+1;
                    
                    %% updating user queue after each slot assignment
                    [dummy,ind_su] = sort(count_su);
                    
                    scheduled=1;
                    break
                end
                if scheduled
                    break
                end
            end
            if scheduled
                break
            end
        end
    end
    
    if rem_slot_previous==rem_slot
        su_state(scheduling_su_ind)=0;
    end
end

