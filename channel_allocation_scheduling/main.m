% This code creates a system model consisting of primary users and
% secondary users. Next, it creates primary user activity and performs the
% energy detection in each cognitive base stations. The central server
% takes the decisions about each channel and classifies them to
% 4 categories of channels as described in section 3.4. Finally, channel
% allocation and scheduling is done using both the benchmark and proposed
% model and figures 7 to 12 are generated
clear all
close all
clc
tic
%% Vary this
T=5;% total observation time in seconds

N_pu=20;
N_su=300;
alpha=1; % average busy time
N_network=20;
OL_fac=.5;

T_ss=5e-3;% sensing time
N_slot=100;% no. of slots in a time frame
t_d=30e-3;% running energy detector every t_d seconds
T_a=5e-3;
T_slot=(t_d-T_ss-T_a)/N_slot;
N=1e3; % no. of samples
N_COMP=2e3;
Fs=N_COMP/T_ss;
Ts=1/Fs;

p_thresh_db=-129.7970;
p_thresh=10^(p_thresh_db/10);
p_thresh_COMP_db=-129.87;
p_thresh_COMP=10^(p_thresh_COMP_db/10);


beta=1; % avg idle time
R_pu=500;
R_su=384.9326;


M = 4;                  % modulation index for psk
hpsk = comm.PSKModulator('ModulationOrder',M,...
    'BitInput',false,...
    'PhaseOffset',0);   % M-psk modulator

%% Creating system model
d0=1;
n=3; % path loss exponent
p_n_dbm=-100;% noise power in dbm
p_n_db=p_n_dbm-30;
p_n=10^((p_n_dbm-30)/10);% noise power in watt

snr_su_lim_db=3;
snr_su_lim=10^(snr_su_lim_db/10);
p_tx_cbs=(R_su/d0)^n*snr_su_lim*p_n;
p_tx_cbs_db=10*log10(p_tx_cbs);
p_tx_cbs_dbm=p_tx_cbs_db+30;% transmission power of CBS in dbm

snr_pu_lim_db=3;
snr_pu_lim=10^(snr_pu_lim_db/10);
p_tx_pu=(R_pu/d0)^n*snr_su_lim*p_n;
p_tx_pu_db=10*log10(p_tx_pu);
p_tx_pu_dbm=10^((p_tx_pu-30)/10); % power in watt

% system model parameters
N_bs=2;
X_max=2500;
Y_max=1500;
OL=R_su*OL_fac;
X_fac=.4;
Y_fac=.5;
temp=X_max/2-(R_su-OL/2);
X_org=[temp temp+2*R_su-R_su*OL_fac] ;
Y_org=[Y_max*Y_fac Y_max*Y_fac];

throughput_mat_v2=zeros(N_network,N_su);
throughput_mat_COMP_v2=zeros(N_network,N_su);



%% Channel allocation constants
B=10e6;

%% PU constants
N_ch=N_pu;
lambda= [1./beta'  1./alpha'];
Fs_ed=Fs/2; % Sampling rate of ED
Fs_ed_COMP=Fs; % Sampling rate of ED

%% directory to save output
ind=1;
keep_going=1;
while keep_going
    folder_name=['Run ' num2str(ind)];
    ind=ind+1;
    if exist(folder_name)==7
        continue
    end
    mkdir(folder_name)
    keep_going=0;
end


for i_network=1:N_network
    i_network
    Tic_3=tic;
    % create a system model using the defined parameters
    [ X_pu ,...
        Y_pu,...
        X_su_bs,...
        Y_su_bs,...
        COMP_compatibility_su,...
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
        R_su,...
        p_tx_cbs,...
        p_tx_pu,...
        d0,...
        n);
    
    COMP_compatibility_su_list=find(COMP_compatibility_su);
    
    snr_su_pu_db=10*log10(p_rx_su_pu/p_n);
    dist_bs_pu=zeros(N_bs,N_pu);
    
    for i=1:N_bs
        dist_bs_pu(i,:)=sqrt((X_org(i)-X_pu).^2+(Y_org(i)-Y_pu).^2);
    end
    
    
    Tic_2=tic;
    
    %% Channel allocation variables
    
    bits_v2=zeros(1,N_su);
    bits_COMP_v2=zeros(1,N_su);
    C=B*log2(1+p_rx_su_bs/p_n);% SU capacity
    
    for k=1:N_ch
        %% Creating PU signal
        [ch_tx,t_ch,pu_state]=gen_PU_PSK(T,lambda,Fs,M,hpsk);
        ch_tx_mat(k,:)=ch_tx;
        pu_state_mat(k,:)=pu_state;
        pu_state=[];
    end
    
    %% Setting up variables for energy detection
    T_d=t_ch(1:round(t_d/(1/Fs_ed)):end);% time array for running ED
    T_d_COMP=t_ch(1:round(t_d/(1/Fs_ed_COMP)):end);% time array for running ED
    d=[];% detection result for all channels from each SU of a parular BS
    EOF_sig=0;
    pu_state_observed=zeros(1,N_ch);% does the
    % signal really exist for a substantial period of T_ss
    
    for i=1:T/t_d
        % ==>
        ind=find(t_ch>=T_d(i),1);
        i
        Tss_end=ind+Fs/Fs_ed*N_COMP-1;
        if (ind+Fs/Fs_ed*N_COMP-1)>T*Fs% checking if we have reached the end of signal
            break
        end
        
        for i_N_ch=1:N_ch
            pu_state_observed(i_N_ch)=sum(pu_state_mat(i_N_ch,ind:ind+Fs/Fs_ed*N_COMP-1))...
                /length(pu_state_mat(i_N_ch,ind:ind+Fs/Fs_ed*N_COMP-1))>.5;
        end
        
        %% Energy Detection by base stations
        % In the following loop each station senses every
        % channel and stores their decision about each channel
        % regarding idle/busy and COMP/non-COMP
        
        d_bs=zeros(N_bs,N_ch);% detection result of each base
        d_bs_COMP=zeros(N_bs,N_ch);% detection result of each base
        
        for i_N_bs=1:N_bs
            d=zeros(1,N_ch);
            d_COMP=zeros(1,N_ch);
            
            for i_N_ch=1:N_ch
                
                temp1=sqrt(p_rx_bs_pu(i_N_bs,i_N_ch));
                
                ch_rx=temp1*ch_tx_mat(i_N_ch,:);
                
                noise=(randn(1,N)+1i*randn(1,N))./(sqrt(2));
                ch_rx=ch_rx(ind:Fs/Fs_ed:round(ind+Fs/Fs_ed*N-1))...
                    + sqrt(p_n)*noise;
                
                p_r=sum((abs(ch_rx)).^2)/(length(ch_rx)-1);%recieved energy
                
                if p_r > p_thresh % if energy is greater than threshold then signal is present
                    d(i_N_ch) = 1;
                    
                else d(i_N_ch)=0;
                    
                end
                if p_r > p_thresh_COMP % if energy is greater than threshold then signal is present
                    d_COMP(i_N_ch) = 1;
                    
                else d_COMP(i_N_ch)=0;
                end
            end
            d_bs(i_N_bs,:)=d_bs(i_N_bs,:)|d;
            d_bs_COMP(i_N_bs,:)=d_bs_COMP(i_N_bs,:)|d;
            
        end
        
        
        % The decision are sent to central server which performs channel classification
        F_ij=zeros(1,N_ch);% JT incompatible free channels in both BS
        % for a particular T_d_COMP(i)
        for i_N_bs=1:N_bs
            F_ij=F_ij|(d_bs(i_N_bs,:)==0 & d_bs_COMP(i_N_bs,:)==1);
        end
        F_ij_list=find(F_ij);
        
        F_ij_JT=zeros(1,N_ch);% JT compatible free channels in both BS
        % for a particular T_d_COMP(i)
        for i_N_bs=1:N_bs
            F_ij_JT=F_ij_JT|d_bs_COMP(i_N_bs,:);
        end
        F_ij_JT=~F_ij_JT;
        F_ij_JT_list=find(F_ij_JT);
        
        F_i=zeros(1,N_ch);% free channel in 'i' BS only
        for i_N_ch=1:N_ch
            
            if d_bs(1,i_N_ch)==0 && d_bs(2,i_N_ch)==1
                F_i(i_N_ch)=1;
            end
        end
        F_i_list=find(F_i);
        
        
        F_j=zeros(1,N_ch);% free channel in 'j' BS only
        
        for i_N_ch=1:N_ch
            
            if d_bs(2,i_N_ch)==0 && d_bs(1,i_N_ch)==1
                F_j(i_N_ch)=1;
            end
        end
        F_j_list=find(F_j);
        
        
        
        %% Scheduling and channel allocation
        %         v2 ==> benchmark (without COMP)
        %         COMP_v2 ==> proposed algorithm (with COMP)
        
        
        [bits_v2]=...
            scheduling_RR(cell_association_state,...
            F_i_list,...
            F_j_list,...
            F_ij_list,...
            F_ij_JT_list,...
            N_su,...
            N_ch,...
            N_slot,...
            T_slot,...
            bits_v2,...
            p_rx_su_bs,...
            p_n,...
            B);
        
        
        [bits_COMP_v2]=...
            scheduling_proposed(cell_association_state,...
            F_i_list,...
            F_j_list,...
            F_ij_list,...
            F_ij_JT_list,...
            N_su,...
            N_ch,...
            N_slot,...
            T_slot,...
            bits_COMP_v2,...
            p_rx_su_bs,...
            p_rx_su_bs_COMP,...
            p_n,...
            B,...
            COMP_compatibility_su_list);
        
    end % end of T
    
    throughput_v2=bits_v2/T;
    throughput_mat_v2(i_network,:)=throughput_mat_v2(i_network,:)+throughput_v2;
    throughput_COMP_v2=bits_COMP_v2/T;
    throughput_mat_COMP_v2(i_network,:)=throughput_mat_COMP_v2(i_network,:)+throughput_COMP_v2;
    
    time_iter(i_iter)=toc(Tic_2)/60;
    disp([num2str(time_iter(i_iter)) ' minutes elapsed for each iteration'])
    
    
    throughput_mat_v2(i_network,:)= throughput_mat_v2(i_network,:);
    throughput_mat_COMP_v2(i_network,:)= throughput_mat_COMP_v2(i_network,:);
    
    %% COMP compatible SU throughput
    %% Total
    COMP_SU=find(COMP_compatibility_su);
    throughput_v2_mid_BS(i_network)=sum(throughput_mat_v2(i_network,COMP_compatibility_su>0));
    throughput_COMP_v2_mid_BS(i_network)=sum(throughput_mat_COMP_v2(i_network,(COMP_compatibility_su>0)));
    
    
    time_network(i_network)=toc(Tic_3)/60;
    disp([num2str(time_network(i_network)) ' minutes elapsed for each network'])
    
    
    N_bs_i(i_network)=sum(cell_association_state==1);
    N_bs_j(i_network)=sum(cell_association_state==2);
    
    save([folder_name '\Workspace_Network_' num2str(i_network)])
    
    
end


%% System Total throughput
throughput_v2_total=sum(throughput_mat_v2,2);
throughput_COMP_v2_total=sum(throughput_mat_COMP_v2,2);

%% Throughput gain
gain_throughput_total=(throughput_COMP_v2_total-throughput_v2_total)./throughput_v2_total*100;
gain_throughput_total_mid_BS=(throughput_COMP_v2_mid_BS-throughput_v2_mid_BS)./throughput_v2_mid_BS*100;

%% mean
gain_throughput_total_mean=mean(gain_throughput_total);
gain_throughput_total_mid_BS_mean=mean(gain_throughput_total_mid_BS);
%% max/min
gain_throughput_total_max=max(gain_throughput_total);
gain_throughput_total_min=min(gain_throughput_total);
gain_throughput_total_mid_BS_max=max(gain_throughput_total_mid_BS);
gain_throughput_total_mid_BS_min=min(gain_throughput_total_mid_BS);

marker=['o'	,'+','*','x','s','d', '^','v','>','<','p','h'];

%% minimum throughput of CRNs
for i=1:N_network
    throughput_v2_min(i)=min(throughput_mat_v2(i,:));
    throughput_COMP_v2_min(i)=min(throughput_mat_COMP_v2(i,:));
end
throughput_inc_avg=mean((throughput_COMP_v2_min-throughput_v2_min)/throughput_v2_min)*100;
%% Calculating Fairness of CRN-1
FI_v2=sum(throughput_mat_v2(1,:))^2/(N_su*sum(throughput_mat_v2(1,:).^2));
FI_COMP_v2=sum(throughput_mat_COMP_v2(1,:))^2/(N_su*sum(throughput_mat_COMP_v2(1,:).^2));

for i=1:N_network
    FI_mat_v2(i)=sum(throughput_mat_v2(i,:))^2/(N_su*sum(throughput_mat_v2(i,:).^2));
    FI_mat_COMP_v2(i)=sum(throughput_mat_COMP_v2(i,:))^2/(N_su*sum(throughput_mat_COMP_v2(i,:).^2));
end


%% plotting base stations, PU and SU (Figure 7)

figure(1)

h= plot(X_org,Y_org,'r^','markersize', 8,'linewidth',1.5);
hold on

plot(X_pu,Y_pu,'ks','markersize', 8,'linewidth',1.5);
labels=cellstr(num2str((1:N_pu)'));

hold on

plot(X_su_bs,Y_su_bs,'k.','markersize', 8,'linewidth',1.5);
legend(cellstr(['CBS';
    'PU ';
    'SU ']))

plot(R_su*cos(2*pi*(0:.001:1))+X_org(1),R_su*sin(2*pi*(0:.001:1))+Y_org(1),'linewidth',1.5)
plot(R_su*cos(2*pi*(0:.001:1))+X_org(2),R_su*sin(2*pi*(0:.001:1))+Y_org(2),'linewidth',1.5)

uistack(h,'top')
axis equal
axis([0 X_max+20 0 Y_max])
hold off

xlabel('X-coordinate (m)')
ylabel('Y-coordinate (m)')



%% plotting user throughput of CRN-1 (Figure 8 and 9)
figure()
bar(throughput_mat_v2(1,:)/10^6)
xlim([0 N_su+1])
ylabel('Throughput (Mbps)','FontSize',12)
xlabel('User Index','FontSize',12)

figure()
bar(throughput_mat_COMP_v2(1,:)/10^6)
xlim([0 N_su+1])
ylabel('Throughput (Mbps)','FontSize',12)
xlabel('User Index','FontSize',12)

%% plotting fairness (Figure 10)
figure(2)
bar(FI_mat_COMP_v2,'FaceColor',[195 224 246]/255)
hold on
bar(FI_mat_v2,0.4,'b')
ylim([0 1.2])
xlim([0 N_network+1])
ylabel('FI','FontSize',12)
xlabel('Network Index','FontSize',12)

legend(cellstr(['CoMP   ';
    'No CoMP';]))


%% plotting min throughput comparison of all networks (Figure 11)

figure(103)
hold on

bar(throughput_COMP_v2_min/10^6,'FaceColor',[195 224 246]/255)
bar(throughput_v2_min/10^6,0.4,'b')

set(gca,'XTick',index)
xlim([0 N_network+1])

ylabel('Throughput (Mbps)','FontSize',12)
xlabel('Network Index','FontSize',12)

legend(cellstr(['CoMP   ';
    'No CoMP';]))

ylim([0 0.8])
box on

%% plotting total throughput comparison of all networks (Figure 12)

index=1:N_network;
% Plotting Throughput
figure(102)

bar(throughput_v2_total/10^6,'FaceColor',[195 224 246]/255)
hold on
bar(throughput_COMP_v2_total/10^6,0.4)

set(gca,'XTick',index)
xlim([0 N_network+1])
ylim([0 290])
ylabel('Throughput (Mbps)','FontSize',12)
xlabel('Network Index','FontSize',12)

legend(cellstr(['No CoMP';
    'CoMP   ';]))



%% Plotting gain

figure(104)

bar(gain_throughput_total)

set(gca,'XTick',index)
xlim([0 N_network+1])

ylabel('Throughput Gain ($\%$)','FontSize',12)
xlabel('Network Index','FontSize',12)

time_total=toc/3600;
disp([num2str(time_total) ' hours elapsed for complete simulation '])


