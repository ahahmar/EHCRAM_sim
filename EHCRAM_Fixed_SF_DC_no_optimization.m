clc
close all;
clear all;

%%% 1 - DEFINITION OF LoRa PARAMETERS
T = 60;           % Total simulation time - 60 mins
counter = 1;        % counting variable
T0 = 0;             % initial value of T0
Window = 1;         % size of window observation
nrofdevices = 100;
nrofchannels = 170; %no. of FHSS channels %ceil(freqspan /freqinterval);
nrofpackets = 5;    %[packet] no. of packet from each LoRa ED
BW = 125e3;         %[Hz] Bandwidth - 125khz
start_SF = 1;       % SF values used
end_SF = 1;         % representing for the no of SF values used
nrofSF = end_SF;    % no. of SF values used SF: 7 - 12 
CR = 1;             % Coding rate
Header = 0;         % [byte] LoRaWAN Header: at least 13 for regular ULs
IH = 0; %IH = 0 --> Expelicit header %IH = 1 --> Implicit header
CRC = 0;            % Correction Rate Consequence
L = 4.375;             % [bytes] the length of the Payload
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
% arr_time =1/0.1*ones(nrofdevices,1); % (s) this value is the arrival time of packet. 
% lambda_data = 1./arr_time; % (pkt/s) the number of packets per second

%--- 2.2 - Energy traffic parameters
arr_time_ener =1/0.5*ones(nrofdevices,1); % (s) this value is the arrival time of packet. 
lambda_ener = 1./arr_time_ener; % (pkt/s) the number of packets per second

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
        %np(nn) = ceil(0.01*randi([1,6],1)*lambda_data(nn)*Window*Trans_Int);
        np(nn) = lambda_data(nn)*Window*Trans_Int;
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
%---- Thien adds on 29/4/2020
%--- For data traffic in the current period 
    nP= sum(np);
    nP_count(counter) = nP;   
    
%--- For energy traffic in the current period 
    nP_ener= sum(np_ener);
    nP_ener_count(counter) = nP_ener;   
    
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
    % In theory
    D11 = success_1; %nP - pkt_drop_theory;
    D11_count(counter) = D11; 
    
%--- Determining the energy of EDs in the current period
    Etx(:,counter) = E_tx.';
    Esl(:,counter) = E_sl.';
    Ec(:,counter) = E_c.';
    Ec2 = E_c.';
    Eh(:,counter) = E_h;
    Er(:,counter) = E_r;
    EA = sum(E_c)/D11;   % Aver. energy consumed for a successful packet
    EA_count(counter) = EA;
    Life(:,counter) = (Trans_Int/3600).*(Er(:, counter)./Ec(:, counter));
    T_sl_count(counter) = T_sl;
    Trans_Int_count(counter) = Trans_Int;
    ToA_count (counter) = ToA;
    DC(counter) = DC0;
    E_ini = E_r;
    SF(counter) = SF0;
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

%---4.2 Measured and Estimated Data packets
nP_count; %assigned to D1;
D1; % Incoming packets in the current period D1 = nP_new_mat = x1_count
D11_count; %Succesfully transmitted packets in the current period
D1_ener; % for solar-based EH - Tested on 3 May 2020

%%%%%%--- Performance results
tt = 1:1:counter-1;
Obser_time = T0_mat - T0_mat_2;

%--- Packet Loss Ratio (PLR) and Throughput
tot_trans_pkt = sum(nP_count); 
tot_drop_pkt = sum(fail_count);
tot_succ_pkt = sum(success_count);
Th_count = (success_count)./Obser_time;
aver_PLR = mean(PLR_count);
aver_PDR = mean(PDR_count);
aver_Th = mean(Th_count);
aver_Th_per_ED = mean(Th_count./nrofdevices); % Aver_Throughput per each ED

%--- Calculating the harvested/energy efficiency
aver_T_sl = mean(T_sl_count);
aver_EA = mean(EA_count);
aver_EH = mean(sum(Eh));
aver_EC = mean(sum(Ec));
aver_Eff_Eh = mean(sum(Eh))/tot_succ_pkt;
aver_Eff_Ec = mean(sum(Ec))/tot_succ_pkt;
tot_Eff_Eh = sum(sum(Eh))/tot_trans_pkt;
tot_Eff_Ec = sum(sum(Ec))/tot_trans_pkt;
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


Harvested_energy_efficiency= (aver_EH/tot_Eff_Ec);


fprintf('EH= %d, \n',aver_EH);
fprintf('EC= %d, \n',tot_Eff_Ec);
fprintf('Harvested Energy efficiency = %d, \n',Harvested_energy_efficiency);
fprintf('Succ.Pkt = %d, \n',tot_succ_pkt);
fprintf('Drop.Pkt = %d, \n',tot_drop_pkt);
fprintf('Throuput = %d, \n',aver_Th_per_ED);
fprintf('Average Harvested Energy = %d, \n',aver_EH);
fprintf('Average Consumed Energy = %d, \n',aver_EC);

%pdr = (tot_drop_pkt*100)/(tot_succ_pkt+tot_drop_pkt);
fprintf('PDR = %2f, \n',aver_PDR);

fprintf('Useful # of bits = %2f, \n',tot_succ_pkt*L*8);
fprintf('Useful energy per pkt J = %2f, \n',((aver_EC/(tot_succ_pkt*L*8))/nrofdevices)); % micro joule

ener1=(aver_EC/(tot_succ_pkt*nrofdevices*8*L))*1000000;
fprintf('Useful energy per pkt microJ = %2f, \n',ener1);

%LoRaWAN,
enerlorawan=(aver_EC/(205*nrofdevices*8*L))*1000000;
fprintf('Useful energy per pkt microJ LoRaWAN = %2f, \n',enerlorawan);


%--- Figure 4 - SF and DC values with optimisation
%figure
%plot (Life_1, '-r','linewidth',1.5)
%ylabel({'Time'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
%xlabel({'Lifetime of an ED'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
%grid on

%figure
%subplot(1,2,1)
%grid on
%plot(tt, SF,'-r','linewidth',1.5)
%ylabel({'SF Values'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
%xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
%subplot(1,2,2)
%grid on
%plot(tt, DC,'-r','linewidth',1.5)
%ylabel({'DC Values'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
%xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')

t = 0;
%end
%--- Figure - FHSS Channels Distribution
% figure; 
% pcolor(ft);
% map2 = [0 0 0;  0 1 0 ;0 0.5 0; 0.6 0 0;];
% colormap(map2);
% shading flat;
% ylabel({'Observation Time (ms)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
% xlabel({'Cryptography FHSS Channels'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
% ax = gca;

%--- Figure 2b - Incoming data Traffic with Kalman Filter
% figure
% grid on
% plot(T0_mat_2, D1,'-b', T0_mat_2, State,'r','linewidth', 1.2)
% ylabel({'Incoming Data Packets'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
% xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
% legend({'Actual Values', 'Estimated Values'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
% 

%--- Figure 3b - Incoming energy Traffic with Kalman Filter
% figure
% grid on
% plot(T0_mat_2, D1_ener,'-b', T0_mat_2, State_ener,'r','linewidth', 1.2)
% ylabel({'Incoming Energy Packets'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
% xlabel({'Simulation Time (sec)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')
% legend({'Actual Values', 'Estimated Values'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Helvetica')

