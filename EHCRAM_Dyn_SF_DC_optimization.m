clc
close all;
clear all;

%%% 1 - DEFINITION OF LoRa PARAMETERS
T = 300;           % Total simulation time - 60 mins
counter = 1;        % counting variable
T0 = 0;             % initial value of T0
Window = 1;         % size of window observation
nrofdevices = 100;
nrofchannels = 170; %no. of FHSS channels %ceil(freqspan /freqinterval);
nrofpackets = 1;    %[packet] no. of packet from each LoRa ED
BW = 125e3;         %[Hz] Bandwidth - 125khz
start_SF = 1;       % SF values used
end_SF = 1;         % representing for the no of SF values used
nrofSF = end_SF;    % no. of SF values used SF: 7 - 12 
CR = 1;             % Coding rate
Header = 0;         % [byte] LoRaWAN Header: at least 13 for regular ULs
IH = 0; %IH = 0 --> Expelicit header %IH = 1 --> Implicit header
CRC = 0;            % Correction Rate Consequence
L = 51;             % [bytes] the length of the Payload
Pr_sym = 8;         % [symbols] Preamable symbols
DE_ini = [1; 1; 0; 0; 0; 0]; %DE is set to 1 if T_sym exceeds 16ms
SF_ini = [12; 11; 10; 9; 8; 7];
DC_ini = [0.1; 0.01; 0.001]; % Set of duty cycle: 10%, 10%, 0.1% 
packetduration_ini = [2;1;0.6;0.3;0.2;0.1];
Tx_power = 14;      %[dBm] LoRa Moderm Tool
U = 3.3;            % [V] Voltage is applied for LORA32U4II
I_tx = 44e-3;       %[A] LoRa Moderm Tool /Transmit current :128mA for 70mS
I_rx = 10.8e-3; %[A] % LoRa Moderm Tool /Receive current without sleep: 14mA
I_sl = 100e-9; %[A] % LoRa Moderm Tool /Current receive + sleep : 1mA
I_deep_sl = 300e-6; % [A] Current super sleep : 300uA
E_ini = 10*ones(nrofdevices,1); %[J]
time_offset_tot = [];

%%%%%%%%%%%% 2 - DATA TRAFFIC AND HARVESTED ENERGY INFORMATION
%------ Thien adds on 29/4/2020 and 3/5/2020
%--- 2.1 - Data traffic parameters
%arr_time1 =1/0.1*ones(nrofdevices,1); % (s) this value is the arrival time of packet. 
%lambda_data1 = 1./arr_time1; % (pkt/s) the number of packets per second

%--- 2.2 - Energy traffic parameters
arr_time_ener =1/0.5*ones(nrofdevices,1); % (s) this value is the arrival time of packet. 
lambda_ener = 1./arr_time_ener; % (pkt/s) the number of packets per second

%--- 2.3 - Kalman Filter parameters
Q = 1.1; % variance of state noise
R = 0.95;  % variance of observation noise   
UseTrueVariances = 'true';
lookbackWindow = 10;

%--- 2.3.1 - For data
x0 = 0;
z0 = x0 + sqrt(Q)*randn;
smoothed_z0 = z0;
M0 = 1; 
Mp0 = 0;
State0 = 0; 
StateP = 0; 
K0 =0;      
nP0 = 0;
nP02 = 0;
aa=1; % state parameter
sigma2_u=1.2; % variance of state noise
sigma2_w=0.95; % variance of observation noise 
S0=1; % initial state
X0=S0+sqrt(sigma2_u)*randn;

%--- 2.3.2 - For solar-based EH - Tested on 3 May 2020
x_ener0 = 0;
z_ener0 = x_ener0 + sqrt(Q)*randn;
smoothed_ener_z0 = z_ener0;
M_ener0 = 1; 
Mp_ener0 = 0;
State_ener0 = 0; 
State_enerP = 0; 
K_ener0 =0;      
nP_ener0 = 0;
nP_ener02 = 0;
aa_ener=1; % state parameter
sigma_ener2_u=1.2; % variance of state noise
sigma_ener2_w=0.95; % variance of observation noise 
S_ener0=1; % initial state
X_ener0=S_ener0+sqrt(sigma_ener2_u)*randn;

%%%%%%%%%%%%%%% 3 - INITIALISING PARAMETERS FOR PROPOSED ALGORITHM
%--- 3.1 - Initilising the value of SF
sf = randi([start_SF end_SF]); % the index of SF
SF0 = SF_ini(sf,1);     % this is initial value of SF --> random value
DE0 = DE_ini(sf,1);     % this is initial value of DE
DC0 = DC_ini(2,1);      % select the initial value of DC: 1%
packetduration0 = packetduration_ini(find(SF_ini== SF0),1);
timespan = 3*Window*60;   %Window*Trans_Int = observation time
timeinterval = 10e-3;
nrofslots = timespan/timeinterval;
interval = 1;
devicestepsize = nrofdevices/interval;
results = zeros(nrofdevices/devicestepsize, 2);
nr_of_dev = devicestepsize:devicestepsize:nrofdevices;
arr_time0 = 0.01; % arrival time of packet 0.1 = 6pkts/60s

while (T0 < T),
%--- 3.2 - Calculating LoRa parameters
    Rate = SF0*(BW/(2^SF0))*(4/(4+CR)); % [bps]
    T_sym = ((2^SF0)/BW);               %[s];
    T_Pr_sym = (Pr_sym + 4.25)*T_sym;   %[s]
    pay_sym = 8 + max(ceil((8*(L+Header)-4*SF0 + 28 + 16*CRC - 20*IH)/(4*(SF0 - 2*DE0)))*(CR+4), 0);
    T_pay_sym = pay_sym*T_sym;          %[s]
    ToA = T_pay_sym + T_Pr_sym;         %[s] packet duration/or ToA
    T_sl = ToA*((1 - DC0)/DC0);         % [s] sleeping period
    Ts = L*8/Rate;                      % [s] trans. time of a pkt in LoRa
    Trans_Int = ToA/DC0;                %[s] Transmission Interval/Cycle
    lora_duration = [SF0 Rate ToA DE0]; 
    packetduration = packetduration0;
    %packetduration = round(ToA);     %[s] packet duration/or ToA
    %packetduration = [2.3020; 1.2329; 0.5755; 0.3287; 0.1746; 0.0975]%[s]
    lambda_data = arr_time0*ones(nrofdevices,1);
    
%--- 3.3 - Calculating information in the current period - the i_{th}
%--- 3.3.1 - Initialising the incoming data and energy packets
%--- the data and energy packets are observed in every observation period
%--- observation period = Window*Trans_Int or Window*Transmission_Cycle
    for nn=1:1:nrofdevices
        %--- For data traffic
        np(nn) = (randi([1,2],1)*lambda_data(nn)*Window*Trans_Int);
        %--- For solar-based EH - Tested on 3 May 2020
        if (T0>=0 && T0 <=7200)                 %5-7am
            np_ener(nn) = ceil(0.05*randi([1,3],1)*lambda_ener(nn)*Window*Trans_Int);
            eh_unit(nn) = 25e-5;%1e-5*randi([10,15],1)*rand(1);%([32,76],1);
        elseif (T0>7200 && T0 <=14400)          %7-9am
            np_ener(nn) = ceil(0.08*randi([3,6],1)*lambda_ener(nn)*Window*Trans_Int);
            eh_unit(nn) = 15e-5;%1e-5*randi([26, 40],1)*rand(1);%([76, 228],1)
        elseif (T0>14400 && T0 <=25200)         %9-12pm
            np_ener = ceil(0.1*randi([6,10],1)*lambda_ener(nn)*Window*Trans_Int);
            eh_unit(nn) = 15e-5;%1e-5*randi([40, 80],1)*rand(1); % ([228, 675],1)
        elseif (T0>25200 && T0 <=32400)         %12-2pm
            np_ener(nn) = ceil(0.12*randi([7,10],1)*lambda_ener(nn)*Window*Trans_Int);
            eh_unit(nn) = 15e-5;%1e-5*randi([80, 110],1)*rand(1); %([385, 675],1)
        elseif (T0>32400 && T0 <=43200)         %2-5pm
            np_ener (nn) = ceil(0.11*randi([7,10],1)*lambda_ener(nn)*Window*Trans_Int);
            eh_unit(nn) = 15e-5;%1e-5*randi([45, 70],1)*rand(1); %([60, 385],1)
        elseif (T0>43200 && T0 <=54000)         %5-8pm
            np_ener(nn) = 0.0001;
            eh_unit(nn) = 15e-5;%1e-5*randi([10, 20],1)*rand(1); %([1, 60],1)
        elseif (T0>54000 && T0 <=64800)         %8-11pm
            np_ener = 0.0001;
            eh_unit(nn) = 1e-5*randi([1, 2],1)*rand(1);
        elseif (T0>64800 && T0 <=86400)         %11-5am
            np_ener(nn) = 0.0001;
            eh_unit(nn) = 1e-5*randi([1, 2],1);
        elseif (T0>86400 && T0 <=93600)         %2-5am
            np_ener(nn) = ceil(0.05*randi([1,3],1)*lambda_ener(nn)*Window*Trans_Int);
            eh_unit(nn) = 1e-5*randi([1, 40],1);
        elseif (T0>93600 && T0 <=100800)        %5-7am
            np_ener(nn) = ceil(0.08*randi([3,6],1)*lambda_ener(nn)*Window*Trans_Int);
            eh_unit(nn) = 1e-5*randi([32, 70],1);
        end
    end
    np_ED = np; % Number of data packet for each ED
    np_ED_count(:,counter) = np_ED.';

    ep = np_ener;
    EP(:,counter) = ep.';   % incoming energy packets in current period
    ener_unit = eh_unit;
    ener_unit_count(:,counter) = eh_unit.';

%---- 3.3.2 - Generating the arrival data and energy packets 
% This step can be replaced by using the random access of EDs (CRAM - FHSS)
%--- For data traffic in the current period 
    nP= sum(np);
    nP_count(counter) = nP;   
    nP0_count(counter) = nP0;
    nP0 = nP;
    nP_new = nP_count(counter) - nP0_count(counter);
    nP_new_count(counter) = nP_new;

%--- For energy traffic in the current period 
    nP_ener= sum(np_ener);
    nP_ener_count(counter) = nP_ener;   
    nP_ener0_count(counter) = nP_ener0;
    nP_ener0 = nP_ener;
    nP_ener_new = nP_ener_count(counter) - nP_ener0_count(counter);
    nP_ener_new_count(counter) = nP_ener_new;

%--- Real measured incoming energy packets in the current period
    D1_ener(counter) = nP_ener; 
    
%--- Generating the random access of EDs to transmit data on FHSS channels
%--- Calculate exact number of ED can be approved to access FHSS channels 
%--- Generate random FHSS Channel from 1 - 170
    for nr = nr_of_dev   
        nr;
        ft = zeros(nrofslots, nrofchannels);
        ft2 = zeros(nrofslots, nrofchannels);
        colission = zeros(nr, nrofpackets);
        freq = randi([1 nrofchannels], [nr nrofpackets]);
        time = randi([1 nrofslots - (packetduration*nrofpackets)/timeinterval], [nr 1]);
        for i = 1:nr
            time_offset = time(i, 1);
            for p = 1:nrofpackets    
                for k = -1:1:1
                    if (freq(i, p)+k<1) || (freq(i, p)+k>nrofchannels)
                        continue
                    end
                    for j = 1:packetduration/timeinterval
                        if j+time_offset > nrofslots
                            continue
                        end
                        if  ft(j+time_offset, freq(i, p)+k) == 0
                            ft(j+time_offset, freq(i, p)+k) = 1;
                            ft2(j+time_offset,freq(i, p)+k) = i; 
                        else
                            ft(j+time_offset, freq(i, p)+k) = ft(j+time_offset, freq(i, p)+k) + 1;
                            colission(i, p) = 1;
                            colission(ft2(j+time_offset,freq(i, p)+k), p) = 1;
                        end
                    end
                end
                time_offset = time_offset + packetduration/timeinterval;
            end
        E_tx(i) = U*I_tx*ToA*p;     % Consumed energy for TX
        E_sl(i) = U*I_sl*T_sl;      % Consumed energy for SL
        E_c(i) = E_tx(i) + E_sl(i); % Total consumed energy
        end
        results(floor(nr/devicestepsize), 1) = sum(sum(colission));
        fail = sum(colission, 2) == nrofpackets;
        results(floor(nr/devicestepsize), 2) = sum(fail);
        results(floor(nr/devicestepsize), 3) = 100*sum(sum(fail))/nr;
        fail_1 = sum(fail);
        success_1 = nrofpackets*nr - sum(fail);
        fail_count(counter) = fail_1;
        success_count(counter) = (nrofpackets*nr) - fail_1;
        PLR = (fail_1/(nrofpackets*nr))*100;
        PLR_count(counter) = (fail_count(counter)./(nrofpackets*nr))*100;
        PDR_count(counter) = (success_count(counter)./(nrofpackets*nr))*100;
        PLR_count_2(counter) = (fail_count(counter)./nP_count(counter))*100;
        PDR_count_2(counter) = (success_count(counter)./nP_count(counter))*100;
    end
    
%--- Determing the harvested energy
    E_h = T_sl.*(ep.').*(eh_unit.');
    E_r = E_ini - E_c.' + E_h;

%--- In reality: some EDs cannot access/transmit data on FHSS channels
    nP_theory = sum(np);   % In theory
    
%--- The incoming packet in the current period D1 = Di;
    % In theory, NP = nrofdevices*nrofpackets
    D1(counter) = nP; %% this code is used on 2 May 2020  
    
%--- The successfultransmitted packets in the current period D11 = Di'
    D11 = success_1; %nP - pkt_drop_theory;
    D11_count(counter) = D11; 
    
%--- 3.3.2 - Estimating data/energy in the next period with Kalman Filter
%--- For data traffic
    EWMALength_KF = Window*Trans_Int;
    lambda_KF = 1-2/(EWMALength_KF+1);
    x1 = nP_count(counter);%nP_new_count(counter);
    z1 = x1 + sqrt(R)*randn;
    smoothed_z1 =  lambda_KF*smoothed_z0 + (1-lambda_KF)*z1;   
    x1_count(counter) = x1;
    smoothed_z1_count(counter) = smoothed_z1;
    smoothed_z0 = smoothed_z1; 

    S1=aa*nP;
    X1=S1+sqrt(sigma2_w^(counter))*randn;
    S(counter)=aa*nP_count(counter);%S0+sqrt(sigma2_u)*randn;
    X(counter)=S(counter)+sqrt(sigma2_w^(counter))*randn;
    StateP(counter)=aa*nP_count(counter);%State0; % predic
    Mp1 = (aa^2)*M0 + sigma2_u;%Q;
    Mp(counter)= (aa^2)*M0 + sigma2_u;%Q; % covariance of prediction
    KF(counter)= Mp(counter)/(sigma2_w^counter+Mp(counter));% Kalman gain
    M(counter)=(1-KF(counter))*Mp(counter); 
    State(counter)=StateP(counter)+KF(counter)*(X(counter)-StateP(counter)); % update 
    State01=StateP(counter)+KF(counter)*(X(counter)-StateP(counter)); % update 
    M0 = Mp1;
    State0 = State01;

%--- Real measured incoming data packets in the next period
    nP2 = smoothed_z1;
    D2(counter) = smoothed_z1;
    
%--- Compare the incoming packets btw 2 adjacent period: nP2 and nP
    % If nP2 > nP --> there is an increase in the no. of incoming pkts
    % --> reduce the value of DC: (DC_new = DC_old - X --> considering)
    % If nP2 < nP --> there is not an increase in the no. of incoming pkts
    % --> remain the value of DC: (DC_new = DC_old - X --> considering)
    delta1(counter) = smoothed_z1 - nP; % nP2 = smoothed_z1;
    
    %--- the value of X: used to change the values of SF
    if (nP2 >= nP)
        X_value(counter) = ceil(log2(nP2/nP)); %nP2/nP; %
    else
        X_value(counter) = floor(log2(nP2/nP)) -1; %0;
    end
    
    %--- Compare the succesfully transmitted data packets in prev. period 
    %--- and the estimated incoming packets in the next period
    delta2(counter) = smoothed_z1 - D11; %nP2 = smoothed_z1;

%--- For energy traffic
    %EWMALength_KF = Window*Trans_Int; % equal to the case of data traffic
    %lambda_KF = 1-2/(EWMALength_KF+1);% equal to the case of data traffic
    x_ener1 = nP_ener_count(counter);
    z_ener1 = x_ener1 + sqrt(R)*randn;
    smoothed_ener_z1 =  lambda_KF*smoothed_ener_z0 + (1-lambda_KF)*z_ener1;   
    x_ener1_count(counter) = x_ener1;
    smoothed_ener_z1_count(counter) = smoothed_ener_z1;
    smoothed_ener_z0 = smoothed_ener_z1; 

    S_ener1=aa_ener*nP_ener;
    X_ener1=S_ener1+sqrt(sigma_ener2_w^(counter))*randn;
    S_ener(counter)=aa_ener*nP_ener_count(counter);
    X_ener(counter)=S_ener(counter)+sqrt(sigma_ener2_w^(counter))*randn;
    State_enerP(counter)=aa_ener*nP_ener_count(counter);
    Mp_ener1 = (aa_ener^2)*M_ener0 + sigma_ener2_u;%Q;
    Mp_ener(counter)= (aa_ener^2)*M_ener0 + sigma_ener2_u;%Q; % covariance of prediction
    KF_ener(counter)= Mp_ener(counter)/(sigma_ener2_w^counter+Mp_ener(counter));% Kalman gain
    M_ener(counter)=(1-KF_ener(counter))*Mp_ener(counter); 
    State_ener(counter)=State_enerP(counter)+KF_ener(counter)*(X_ener(counter)-State_enerP(counter)); % update 
    State_ener01=State_enerP(counter)+KF_ener(counter)*(X_ener(counter)-State_enerP(counter)); % update 
    M_ener0 = Mp_ener1;
    State_ener0 = State_ener01;
    
    %--- Estimate the incoming energy packets in the next period
    D2_ener(counter) = smoothed_ener_z1;%smoothed_z1

%--- Determining the energy of EDs in the current period
    Etx(:,counter) = E_tx.';
    Esl(:,counter) = E_sl.';
    Ec(:,counter) = E_c.';
    Ec2 = E_c.';
    Eh(:,counter) = E_h;
    Er(:,counter) = E_r;
    EA = sum(E_c)/D11;      % Aver. energy consumed per a successful packet
                            % in a observation period
    EA_count(counter) = EA;
    Life(:,counter) = (Trans_Int/3600).*(Er(:, counter)./Ec(:, counter));
    T_sl_count(counter) = T_sl;
    Trans_Int_count(counter) = Trans_Int;
    ToA_count (counter) = ToA;

%--- 3.4 - PROPOSED ALGORITHM
    delta_E = E_h - Ec2;
    delta_E_count(:,counter) = delta_E;
    min_E = min(E_h - Ec2);
    min_E_ind = find((E_h - Ec2) <= min_E);
    A = ((pay_sym + Pr_sym + 4.25)/BW);%*((1-DC0)/DC0);
    B = (sum(Ec2)*(nP2))/(sum(E_h)*D11);
    tau = 2^SF0;
    dc = ((1-DC0)/DC0);
    if (min_E < 0)
        fprintf('Case 1 - DC Optimisation \n');
        cvx_begin 
            cvx_precision high
            expression dc_opt
            variables a1 a2    % this variable allows to be feasible
            obj = dc_opt*A*tau - A*B*tau*dc + 1e-7*a1 + 1e-7*a2;
            maximize obj
            subject to
                0 <= obj;
                99 <= dc_opt - 1e-7*a1;
                dc_opt + 1e-7*a2 <= 999;
        cvx_end
        DC_opt = 1/(1+dc_opt);
        dc_opt_count(counter) = dc_opt;
        DC_opt_initial(counter) = DC_opt;
        DC_new = DC_opt;
                
        %--- Adjusting the value of DC to meet the variation of data traffic
        if (delta2(counter) >=0)
            DC_new = min(DC_opt,(DC0 - (DC0*(ceil(log2(nP2/nP)))/100)));
            fprintf('Case 1.1 (D2 > D11) - Reduce DC \n');
        else
            DC_new = DC0;
            fprintf('Case 1.2 (D2 < D1) - Remain DC  \n');
        end
    elseif (min_E > 0)||(PLR >= 5)
        fprintf('Case 2 - Adjusting Duty Cycle \n');
        if (delta1(counter)>0)
            DC_new = DC0 - (DC0*(ceil(log2(nP2/nP)))/100);
            fprintf('Case 2.1 (D2 > D1) - Reduce DC  \n');
        else
            DC_new = DC0 - (DC0*(floor(log2(nP2/nP)))/100);
            fprintf('Case 2.2 (D2 < D11) - Increase DC  \n');
        end
    end

%--- update the value of SF and DC
    DC_new_count(counter) = DC_new;
    DC(counter) = DC0;
    if (DC_new > 0.01)
        DC_new = 0.01;
        fprintf('Re-updating DC value \n');
    end
    arr_time0_count(counter) = arr_time0;
    arr_time = min(arr_time0,(nrofchannels*DC_new*ToA)); % suitable for SF>7
    arr_time_count(counter) = arr_time;             % suitable for SF>7
    arr_time0 = arr_time;
    %arr_time = nrofchannels*DC_new*ToA; % suitable for SF7
    %arr_time_count(counter) = arr_time; % suitable for SF7
    E_ini = E_r;
    SF(counter) = SF0;
    DC0 = DC_new;
    T0 = T0 + Window*Trans_Int; % the length of the observation windows.
    T0_count(counter) = T0;  % The time to start calculating arrival pkts
    counter = counter +1
end

%%%%%%%%%%%%%%% Recording the data
%---4. RECORDING DATA
%---4.1 The real data is genereated by sensor nodes
T0_mat(1,1:counter-1)= T0_count;
T0_mat_2 = [0, T0_mat(1:length(T0_mat)-1)];
nP_mat(1,1:counter-1) = nP_count;
nP_start = [0,nP_mat(1:length(nP_mat)-1)];
nP_end = nP_mat; 
nP_arr = nP_end - nP_start;
nP_arr_sum = sum(nP_arr);
nP_new_mat(1,1:counter-1) = nP_new_count;

%---4.2 Measured and Estimated Data packets
nP_count; %assigned to D1;
smoothed_z1_count; %estimated from nP_count
D1; % Incoming packets in the current period D1 = nP_new_mat = x1_count
D11_count; %Succesfully transmitted packets in the current period
D2; % Estimated incoming packets in the next period D2 = x12
delta1; %The comparison btw D2 and D1: delta1 = D2 - D1;
delta2; %The comparison btw D2 and D11: delta2 = D2 - D11;
X_value;
delta_E; %comparison btw consumed and harvested energy: delta_E = EH - Ec1
D1_ener; % for solar-based EH - Tested on 3 May 2020
D2_ener; % for solar-based EH - Tested on 3 May 2020

%%%%%%--- Performance results
tt = 1:1:counter-1;
Obser_time = T0_mat - T0_mat_2;

%--- Packet Loss Ratio (PLR) and Throughput
tot_trans_pkt = sum(nP_count);          % Tot no. of generated pkts in 34200s
tot_drop_pkt = sum(fail_count);         % Tot no. of dropped pkts in 34200s
tot_succ_pkt = sum(success_count);      % Tot no. of pkts in 34200 seconds
Th_count = (success_count)./Obser_time; % Take account for N EDs
aver_PLR = mean(PLR_count);             % Take account for 34200 seconds
aver_PDR = mean(PDR_count);             % Take account for 34200 seconds
aver_Th = mean(Th_count);               % Take account for 34200 seconds
aver_Th_per_ED = mean(Th_count./nrofdevices); % Aver_Throughput per each ED
PLR_count_2 = (fail_count./nP_count)*100;
PDR_count_2 = (success_count./nP_count)*100;
aver_PLR_2 = mean(PLR_count_2);         % Take account for 34200 seconds
aver_PDR_2 = mean(PDR_count_2);

%--- Calculating the harvested/energy efficiency
aver_T_sl = mean(T_sl_count);
aver_EA = mean(EA_count);
aver_EH = mean(sum(Eh));% aver. total energy harvested in an observation time 
                        % account for N end-devices in an observation time
aver_EC = mean(sum(Ec));% aver. total energy consumed in an observation time 
                        % account for N end-devices in an observation time
aver_EH_ED = mean(mean(Eh)); % average energy harvested by each ED
aver_EC_ED = mean(mean(Ec)); % average energy consumed by each ED
aver_Eff_Eh = mean(sum(Eh))/tot_succ_pkt;
aver_Eff_Ec = mean(sum(Ec))/tot_succ_pkt;
aver_Eff_Eh_suc_pkt = sum(sum(Eh))/tot_succ_pkt; %EH per a successful pkt
aver_Eff_Ec_suc_pkt = sum(sum(Ec))/tot_succ_pkt; %EC per a successful pkt
tot_Eff_Eh = sum(sum(Eh))/tot_trans_pkt; %EH per a pkt
tot_Eff_Ec = sum(sum(Ec))/tot_trans_pkt; %EC per a pkt 
Eff_Ec_Th = (aver_Th)/((sum(sum(Ec))*lambda_data(1,1))/tot_trans_pkt);
Eff_Eh_Th = (aver_Th)/((sum(sum(Eh))*lambda_data(1,1))/tot_trans_pkt);
Eff_Ec_Th_succ = (aver_Th)/((sum(sum(Ec))*lambda_data(1,1))/tot_succ_pkt);
Eff_Eh_Th_succ = (aver_Th)/((sum(sum(Eh))*lambda_data(1,1))/tot_succ_pkt);
Result_QoS = [aver_PDR, aver_Th];
Results_EH = [aver_T_sl, aver_EH, aver_Eff_Eh, tot_Eff_Eh, Eff_Eh_Th, Eff_Eh_Th_succ];
Results_EC = [aver_EA, aver_EC, aver_Eff_Ec, tot_Eff_Ec, Eff_Ec_Th, Eff_Ec_Th_succ];
Er_1 = Er(1,:);
Eh_1 = Eh(1,:);
Ec_1 = Ec(1,:);
Life_1 = (Trans_Int/3600).*(Er_1./Ec_1);
alpha_10 = aver_EH/tot_Eff_Ec;

figure
subplot(1,3,1)
grid on
plot(tt, SF,'-r','linewidth',1.5)
ylabel({'SF Values'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
subplot(1,3,2)
grid on
plot(tt, DC,'-r','linewidth',1.5)
hold on
plot(DC_opt_initial,'-.b','linewidth',1.5)
hold on 
plot(DC_new_count,'-g','linewidth',1.5)
ylabel({'DC Values'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
legend({'DC', 'DC-opt-ini', 'DC-opt-edited'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
subplot(1,3,3)
grid on
plot(arr_time_count,'-r','linewidth',1.5)
hold on
plot(arr_time0_count,'-.b','linewidth',1.5)
ylabel({'Arrival Data Rate (packet/s)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')

figure 
grid on
plot(T0_mat_2, D1,'-b', T0_mat_2, State,'r','linewidth', 1.2)
ylabel({'Incoming Data Packets'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
legend({'Actual Values', 'Estimated Values'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')


%--- Figure 2a - Kalman Filter for Data Estimation: KF Gain and MSE
figure
grid on
plot(T0_mat_2,KF,'-r','linewidth', 1.2)
ylabel({'Kalman Gain - Data Estimation'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')

figure 
grid on
plot(T0_mat_2,M,'-r','linewidth', 1.2)
ylabel({'MSE - Data Estimation'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')


%--- Figure 3a - Kalman Filter for Energy Estimation: KF Gain and MSE
figure
grid on
plot(T0_mat_2,KF_ener,'-r','linewidth', 1.2)
ylabel({'Kalman Gain - Energy Estimation'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')

figure
grid on
plot(T0_mat_2,M_ener,'-r','linewidth', 1.2)
ylabel({'MSE - Energy Estimation'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')

figure
grid on
plot(T0_mat_2, D1_ener,'-b', T0_mat_2, State_ener,'r','linewidth', 1.2)
ylabel({'Incoming Energy Packets'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
legend({'Actual Values', 'Estimated Values'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')

%--- Figure 4 - SF and DC values with optimisation
figure
plot (Life_1, '-r','linewidth',1.5)
ylabel({'Time'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
xlabel({'Lifetime of an ED'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
grid on

t = 0;

