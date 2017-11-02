function [bits]=...
    scheduling_COMP_user_POV(cell_association_state,...
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
    p_rx_su_bs_COMP,...
    p_n,...
    B,...
    COMP_compatibility_su_list)


su_timeslot=zeros(N_su,N_slot);
...transmit only once at each time sample
    % bs_state=zeros(1,N_bs) ;% variable to switch between BS and transmit to a SU
channel_assigned_su=cell(1,N_su);
ch_state=zeros(N_ch,N_slot);
su_state=ones(1,N_su);

%% updating user queue after each slot assignment
% ind_su contains all user's index
[dummy,ind_su] = sort(bits);

rem_slot=(length(F_ij_JT_list)+length(F_ij_list)+length(F_i_list)+length(F_j_list))*N_slot;

while rem_slot>0 && sum(su_state)>0
    
    rem_slot_previous=rem_slot;
    
    for i=1:length(ind_su)
        if su_state(ind_su(i))==1;
            scheduling_su_ind=ind_su(i);
            break
        end
    end
    
    scheduled=0;
    %% serving the user if it is JT  
     
    if find(COMP_compatibility_su_list==scheduling_su_ind)
        %% looking at F_ij_JT first
        for i_ij_JT=1:length(F_ij_JT_list)
            for i_slot=1:N_slot
                if  su_timeslot(scheduling_su_ind,i_slot)==0 && ...
                        ch_state(F_ij_JT_list(i_ij_JT),i_slot)==0
                    
                    channel_assigned_su{scheduling_su_ind}=...
                        [channel_assigned_su{scheduling_su_ind} F_ij_JT_list(i_ij_JT) ];
                    su_timeslot(scheduling_su_ind,i_slot)=1;
                    ch_state(F_ij_JT_list(i_ij_JT),i_slot)=1;
                    rem_slot=rem_slot-1;
                    %% updating bits count
                    p_temp=p_rx_su_bs(scheduling_su_ind)+p_rx_su_bs_COMP(scheduling_su_ind);
                    C=B*log2(1+p_temp/p_n);
                    bits(scheduling_su_ind)=bits(scheduling_su_ind)+C*T_slot;
                    
                    %% updating user queue after each slot assignment
                    [dummy,ind_su] = sort(bits);
                    
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
        %% If user is not served, looking at F_i/F_j depending on cell association
        if scheduled==0
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
                            
                            %% updating user queue after each slot assignment
                            [dummy,ind_su] = sort(bits);
                            
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
                            
                            %% updating user queue after each slot assignment
                            [dummy,ind_su] = sort(bits);
                            
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
        end
        
        %% Finally, If still not served looking at F_ij
        
        if scheduled==0
            for i_ij=1:length(F_ij_list)
                for i_slot=1:N_slot
                    if  su_timeslot(scheduling_su_ind,i_slot)==0 && ...
                            ch_state(F_ij_list(i_ij),i_slot)==0
                        
                        channel_assigned_su{scheduling_su_ind}=...
                            [channel_assigned_su{scheduling_su_ind} F_ij_list(i_ij) ];
                        su_timeslot(scheduling_su_ind,i_slot)=1;
                        ch_state(F_ij_list(i_ij),i_slot)=1;
                        rem_slot=rem_slot-1;
                        %% updating bits count
                        p_temp=p_rx_su_bs(scheduling_su_ind);
                        C=B*log2(1+p_temp/p_n);
                        bits(scheduling_su_ind)=bits(scheduling_su_ind)+C*T_slot;
                        
                        %% updating user queue after each slot assignment
                        [dummy,ind_su] = sort(bits);
                        
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
        
        %% if the user is not JT
    else
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
                        
                        %% updating user queue after each slot assignment
                        [dummy,ind_su] = sort(bits);
                        
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
                        
                        %% updating user queue after each slot assignment
                        [dummy,ind_su] = sort(bits);
                        
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
    end
    
    %% Next, If still not served looking at F_ij
    
    if scheduled==0
        for i_ij=1:length(F_ij_list)
            for i_slot=1:N_slot
                if  su_timeslot(scheduling_su_ind,i_slot)==0 && ...
                        ch_state(F_ij_list(i_ij),i_slot)==0
                    
                    channel_assigned_su{scheduling_su_ind}=...
                        [channel_assigned_su{scheduling_su_ind} F_ij_list(i_ij) ];
                    su_timeslot(scheduling_su_ind,i_slot)=1;
                    ch_state(F_ij_list(i_ij),i_slot)=1;
                    rem_slot=rem_slot-1;
                    %% updating bits count
                    p_temp=p_rx_su_bs(scheduling_su_ind);
                    C=B*log2(1+p_temp/p_n);
                    bits(scheduling_su_ind)=bits(scheduling_su_ind)+C*T_slot;
                    
                    %% updating user queue after each slot assignment
                    [dummy,ind_su] = sort(bits);
                    
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
    
    %% Finally, If still not served looking at F_ij_JT

    if scheduled==0
        for i_ij_JT=1:length(F_ij_JT_list)
            for i_slot=1:N_slot
                if  su_timeslot(scheduling_su_ind,i_slot)==0 && ...
                        ch_state(F_ij_JT_list(i_ij_JT),i_slot)==0
                    
                    channel_assigned_su{scheduling_su_ind}=...
                        [channel_assigned_su{scheduling_su_ind} F_ij_JT_list(i_ij_JT) ];
                    su_timeslot(scheduling_su_ind,i_slot)=1;
                    ch_state(F_ij_JT_list(i_ij_JT),i_slot)=1;
                    rem_slot=rem_slot-1;
                    %% updating bits count
                    p_temp=p_rx_su_bs(scheduling_su_ind);
                    C=B*log2(1+p_temp/p_n);
                    bits(scheduling_su_ind)=bits(scheduling_su_ind)+C*T_slot;
                    
                    %% updating user queue after each slot assignment
                    [dummy,ind_su] = sort(bits);
                    
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


