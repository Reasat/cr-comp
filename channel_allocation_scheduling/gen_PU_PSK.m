% initial sample of transmitted signal is assumed zero but pu tx state is
% assumed busy, fix it later
function [ch_tx,t_ch, pu_state]=gen_PU_PSK(T,lambda,Fs,M,hpsk)
%% Creating PU activity
Ts=1/Fs;
old = digits(3); % specify number of significant decimal digits
Ts=double(vpa(Ts));
digits(old);%Restore the default accuracy setting for further computations

state=randi(2,1,1)-1;     % state= 1==>busy
% state= 0==>idle
seq=state;   % list of transition between states busy/idle
emission=-1*log(1-rand)/lambda(2);% amount of time the channel will be busy
emission=floor(emission/Ts)*Ts;
while sum(emission)<T
    if state==1% if busy toggle
        state=0;% make idle
        seq=[seq,state];
        temp_busy=-1*log(1-rand)/lambda(1);    
        temp_busy=floor(temp_busy/Ts)*Ts;
        emission=[emission temp_busy];
        
        
    else if state==0% if idle toggle
            state=1;% make busy
            seq=[seq,state];
            temp_idle=-1*log(1-rand)/lambda(2);
            temp_idle=floor(temp_idle/Ts)*Ts;
            emission=[emission temp_idle];
            
        end
    end
end

if sum(emission)>T
 emission(end)=emission(end)-(sum(emission)-T);
end



%% creating PU signal

%% PU signal variables

ch_tx=[];% ch is signal in channel
t_ch=[];% time of channel samples
pu_state=[];

for i=1:length(seq)
    state_ch=seq(i);
    t = Ts:Ts:emission(i); % Sampling times  
    if isempty(t_ch)
        temp=0;
    else temp=t_ch(end);
    end
    t=t+temp;
    t_ch=[t_ch t];
    
    %% creating and modulating signal at PU transmitter
    
    N_t=length(t);
    
    if state_ch==1 % busy
        infoSignal = randi(M,N_t,1)-1;  % random binary signal  (bits = log2(M))
        signal_ch_tx = step(hpsk,infoSignal)';   % M-psk signal
        % y is the transmitted signal
        
    else if state_ch==0 %idle
            
            signal_ch_tx = zeros(1,N_t); % no signal
            
            
        end
    end
    ch_tx=[ch_tx signal_ch_tx];
    pu_state=[pu_state seq(i)*ones(1,length(signal_ch_tx))];
    
end

t_ch=t_ch-Ts;

N=T*Fs;
if(length(t_ch)<N)
    temp=t_ch(end)+Ts*(1:(N-length(t_ch)));
    t_ch=[ t_ch temp ];
    ch_tx=[ch_tx zeros(1,length(temp))];
    pu_state=[pu_state zeros(1,length(temp))];
end
if(length(t_ch)>N)
    t_ch= t_ch(1:N);
    ch_tx=ch_tx(1:N);
    pu_state=pu_state(1:N);
end

% if(length(t_ch)~=N)
%     temp=t_ch(end)+Ts*(1:(N-length(t_ch)));
%     t_ch=[ t_ch temp ];
%     ch_tx=[ch_tx zeros(1,length(temp))];
%     pu_state=[pu_state zeros(1,length(temp))];
% end
% if(length(t_ch)>N)
%     t_ch= t_ch(1:N);
%     ch_tx=ch_tx(1:N);
%     pu_state=pu_state(1:N);
% end


