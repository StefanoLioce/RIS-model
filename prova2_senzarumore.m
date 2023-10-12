clear all
close all
clc
%%
numTx = 1; numRx = 1;
c = physconst('LightSpeed'); %light speed
f = 28e9; % [GHz] - reference frequency
Frequency = 28e9; %[GHz] - antenna element resonance frequncy 
B= .5e9; %bandwidth
freq = f-B/2:.1e9:f+B/2;
lambda = c/f; %wavelength
k0 = 2*pi/lambda; %wave number
%% Design antenna element
Zopt = 50; %Value impedance matching
eps_r = 3.3;
h = 0.5e-3; % high substrate[m]
W = c/(2*f)*sqrt(2/(eps_r+1)); %Patch width
eps_reff_r = (eps_r+1)/2 + (((eps_r-1)/2)*(1+12*h/W)^(-0.5)); %Effective dielectric constant
DL = 0.412*h*((eps_reff_r+0.3)*(W/h+0.264))/((eps_reff_r-0.258)*(W/h+0.8));%Extension length
L = c/(2*f*sqrt(eps_reff_r))-2*DL; %Patch length
Leff = L +2*DL; %Effective length of the patch
LG = 6*h+L; %Ground length
WG = 6*h+W; %Ground width
fr = c/(2*L*sqrt(eps_r)); %Resonant frequency of dominant TM010 mode
frc = c/(2*Leff*sqrt(eps_r)); %Resonant frequency of dominant TM010 mode to include edge effects
q = frc/fr; %Fringe factor

ant = patchMicrostrip;
ant.Length = L;              %patch length along x-axis
ant.Width = W;             %patch width along y-axis
ant.GroundPlaneLength = LG; %ground plane length along x-axis
ant.GroundPlaneWidth = WG;  %ground plane width along y-axis
ant.Height = h; %height of substrate
ant.Substrate = dielectric('Name','METEORWAVE8000','EpsilonR',eps_reff_r,'LossTangent',0.0016,...
    'Thickness',h);
ant.TiltAxis = 'Z'; 

% FeedOffset for real antenna element
ant.FeedOffset = [-L/2+ant.FeedWidth,0];
Rin0 = real(impedance(ant,f));

y = linspace(ant.FeedWidth,L-ant.FeedWidth,1e6);

clear j;
for j = 1:length(y)
    Rin(j) = Rin0*(cos(pi/L*y(j))^2);
end

% Find point near to 50 Ohm for matching

Diff = abs(Rin - ones(1,length(y)).*Zopt);
[m1 index_m1] = min(Diff(1:floor(length(Diff)/2)));
[m2 index_m2] = min(Diff(floor(length(Diff)/2):end));
Z50_1 = Rin(index_m1);
y50_1 = y(index_m1);
Z50_2 = Rin(index_m2 + round(length(Diff)/2)-1);
y50_2 = y(index_m2 + round(length(Diff)/2)-1);
Z50 = [Z50_1 Z50_2];
y50 = [y50_1 y50_2];

%theoretical model
ant.FeedOffset = [y50(1)- L/2 , 0];

p = PatternPlotOptions;
p.Transparency = 0.5;

ant.Tilt=0;
G0 = pattern(ant,f,'Type','realizedgain');
maxG0 = max(G0(:));
disp(['Max antenna TX/RX Gain: ' num2str(maxG0) '[dBi]']);
ant.Tilt = 90; ant.TiltAxis = [0 1 0];

%% 
dx = lambda/2; %spacing along x-direction
dy = lambda/2; %spacing along y-direction
N_ele = 20; %sqrt of number of elemets - even number

RX_loc = [5; 1; 3]; %receiver position
TX_loc = [4; 4; 3 ]; %transmitter position
RIS_loc = [0.0332; 2; 2]; % RIS position

d_T_RIS = norm(TX_loc-RIS_loc); %distance between transmitter and RIS
d_RIS_R=norm(RIS_loc-RX_loc); %distance between receiver and RIS

I_phi=sign(TX_loc(2)-RIS_loc(2));
phi_TX = I_phi* atand ( abs( RIS_loc(2)-TX_loc(2)) / abs(RIS_loc(1)-TX_loc(1)) ); %azimuth angle between normal of RIS and transmitter
% phi_TX = 270 + phi_TX; %reference to x-axis
   
I_theta=sign(TX_loc(3)-RIS_loc(3));
theta_TX=  I_theta * asind ( abs (RIS_loc(3)-TX_loc(3) ) / d_T_RIS ); %elevation angle between normal of RIS and transmitter 

I_theta=sign(RX_loc(3) - RIS_loc(3));
theta_RX=  I_theta * asind( abs(RX_loc(3)-RIS_loc(3))/d_RIS_R ); %elevation angle between normal of RIS and receiver

I_phi=sign(RX_loc(2) - RIS_loc(2)); 
phi_RX=I_phi * atand( abs(RX_loc(2)-RIS_loc(2))/ abs(RX_loc(1)-RIS_loc(1)) );%azimuth angle between normal of RIS and receiver
% phi_RX = 270 + phi_RX; %reference to x-axis

%% Compute reflected electric field 
Ei = 1; %incident electric field
Er = zeros(181,361); %initial value of reflected electric field
Gamma = 1; %reflection amplitude
Phi = zeros(N_ele,N_ele); %initialization of reflection phase 
A = 1; %illuminating amplitude
a = pi/6; %illuminating phase
PHI = linspace(-180,180,361);
THETA = linspace(-90,90,181);
clear z, clear y;
for yy=1:N_ele
    for zz=1:N_ele
        Tx_loc(zz,yy) = {[TX_loc(1); -(yy-1)*dx+TX_loc(2); (zz-1)*dy+TX_loc(3)]};
        Rx_loc(zz,yy) = {[RX_loc(1); -(yy-1)*dx+RX_loc(2); (zz-1)*dy+RX_loc(3)]};
        
        phi_Txx(zz,yy) = sign((-(yy-1)*dx+TX_loc(2))-RIS_loc(2))*atand(abs(RIS_loc(2)-(-(yy-1)*dx+TX_loc(2)))/abs(RIS_loc(1)-TX_loc(1)));
        theta_Txx(zz,yy) =  sign(((zz-1)*dy+TX_loc(3))-RIS_loc(3))*asind(abs(RIS_loc(3)-((zz-1)*dy+TX_loc(3)))/norm([TX_loc(1); -(yy-1)*dx+TX_loc(2); (zz-1)*dy+TX_loc(3)]-RIS_loc ));

        phi_Rxx(zz,yy) = sign((-(yy-1)*dx+RX_loc(2))-RIS_loc(2))*atand(abs((-(yy-1)*dx+RX_loc(2))-RIS_loc(2))/abs(RX_loc(1)-RIS_loc(1)));
        theta_Rxx(zz,yy) =  sign(((zz-1)*dy+RX_loc(3))-RIS_loc(3))*asind(abs(((zz-1)*dy+RX_loc(3))-RIS_loc(3))/norm(RIS_loc-[RX_loc(1); -(yy-1)*dx+RX_loc(2); (zz-1)*dy+RX_loc(3)]));
    end
end
%%
for z = 1:N_ele
    for y = 1:N_ele
        zn = (z-1)*dy; ym = (y-1)*dy;
        Phi(z,y) = k0*(ym*sind(90-theta_RX)*sind(-phi_RX)+zn*cosd(90-theta_RX));
        PhiTx(z,y) = k0*(ym*sind(90-theta_TX)*sind(-phi_TX)+zn*cosd(90-theta_TX));

        while Phi(z,y)>=2*pi
            Phi(z,y)=Phi(z,y)-2*pi;
        end
        while PhiTx(z,y)>=2*pi
            PhiTx(z,y)=PhiTx(z,y)-2*pi;
        end
        if Phi(z,y)>=0 && Phi(z,y)<=pi
            Phiq(z,y) = 0;
        else
            Phiq(z,y) = pi;
        end
    end
end

%%
RIS = phased.URA("Element",ant, ...
    "Size", [N_ele N_ele], ...
    "ElementSpacing",[dx dy],"ArrayNormal","x");

% RIS transmitter part
taperedRISt = clone(RIS);
steeringVector_t = phased.SteeringVector("SensorArray",taperedRISt,'IncludeElementResponse',true);
taperedRIStq = clone(RIS);
steeringVectorq_t = phased.SteeringVector("SensorArray",taperedRIStq,'IncludeElementResponse',true);


%RIS receiver part
%taperedRISr = clone(RIS);
%taperedRISt1 = clone(RIS);
% startTaper_r = taperedRISr.Taper;
% steeringVector_r = phased.SteeringVector("SensorArray",taperedRISr);
%figure(), viewArray(RIS,'ShowTaper',true);
%% 
%tap = ones(N_ele, N_ele);
%taperedRISt.Taper = exp(1i*Phi);
taperedRISt.Taper = exp(1i*Phi); %contiuous phase matrix 
taperedRIStq.Taper = exp(1i*PhiTx); %quantized phase matrix
%taperedRISt1.Taper = exp(1i*Phi1);
%figure(), viewArray(taperedRISt,'ShowTaper',true);

%taperedRISt1.Taper = ones.*exp(1i*Phi1);

%Radiation pattern with continuous phase value
figure(), pattern(taperedRISt,f,-180:1:180,-90:1:90); %3D representation
figure(), pattern(taperedRISt,f,0,-90:1:90); %azimuth cut ->shows elevation angle
figure(), pattern(taperedRISt,f,-180:1:180,0); %elevation cut ->shows azimuth angle
[BW, ang] = beamwidth(taperedRISt,f,'dbDown',3); % -3 dB beamwidth 
maxG0Array = max(max(pattern(taperedRISt,f,'Type','directivity')));
%Radiation pattern with quantized phase value
% figure(), pattern(taperedRIStq,f,-180:1:180,-90:1:90);
% figure(), pattern(taperedRIStq,f,0,-90:1:90);
% figure(), pattern(taperedRIStq,f,-180:1:180,0);
[BWq, angq] = beamwidth(taperedRIStq,f,'dbDown',3);

antenna = phased.IsotropicAntennaElement("FrequencyRange",[f-B f+B]); %isotropic antenna for transmitter and receiver

%Definition of txsite and rxsite
tx_power_real = 1;

tx = txsite("cartesian", ...
    "Antenna",antenna, ...
    "TransmitterPower", tx_power_real,...
    "AntennaPosition",TX_loc, ...
    'TransmitterFrequency',f);

RIStx = txsite("cartesian", ...
    "Antenna",taperedRISt, ...
    "AntennaPosition",RIS_loc, ...
    'TransmitterFrequency',f);

rx = rxsite("cartesian", ...
    "Antenna",antenna, ...
    "AntennaPosition",RX_loc, ...
    "AntennaAngle",[0;90],"ReceiverSensitivity",-150);

RISrx = rxsite("cartesian", ...
    "Antenna",RIS, ...
    "AntennaPosition",RIS_loc, ...
    "ReceiverSensitivity",-150);

%% Evaluation points around Rx
epsilon = 0.005; %[m] distance between evaluation points
N_points = 7;
%x_near = linspace(RX_loc(1)-N_points/2*epsilon,RX_loc(1) +N_points/2*epsilon,N_points);
y_near = linspace(RX_loc(2)-N_points/2*epsilon,RX_loc(2) +N_points/2*epsilon,N_points);
z_near = linspace(RX_loc(3)-N_points/2*epsilon,RX_loc(3) +N_points/2*epsilon,N_points);
for k=1:N_points
    RX_loc_yVariations = [RX_loc(1) y_near(k) RX_loc(3)];
    RX_loc_zVariations = [RX_loc(1) RX_loc(2) z_near(k)];
    rx_yVariations{k} = rxsite("cartesian", ...
        "Antenna",antenna, ...
        "AntennaPosition",RX_loc_yVariations', ...
        "AntennaAngle",[0;90],"ReceiverSensitivity",-150);
    rx_zVariations{k} = rxsite("cartesian", ...
        "Antenna",antenna, ...
        "AntennaPosition",RX_loc_zVariations', ...
        "AntennaAngle",[0;90],"ReceiverSensitivity",-150);
end



% for kk=1:N_points
%     for kkk=1:N_points
%         RX_loc_near = [RX_loc(1) y_near(kk) z_near(kkk)];
%         rx_near{kk,kkk} = rxsite("cartesian", ...
%             "Antenna",antenna, ...
%             "AntennaPosition",RX_loc_near', ...
%             "AntennaAngle",[0;90],"ReceiverSensitivity",-150);
%     end
% end

%% Visualizzazione stanza
filename = "Ingenious Juttuli-Amberis (1).stl"; % Nome file formato stl
stldata = stlread(filename);

fontsize = 13;
equalaspect = true;
% scomposizione in lista di connettività e punti
M = stldata.ConnectivityList;
x0 = stldata.Points(:,1);
y0 = stldata.Points(:,2);
z0 = stldata.Points(:,3);

% scalatura valori
scale = 3e-2;

xmin = min(x0);
ymin =min(y0);
zmin = min(z0);

x = x0 - xmin;
y = y0 - ymin;
z = z0- zmin;

Stanza = triangulation(M,x.*scale, y.*scale, z.*scale);
xURC = 6.42318; xLLC = 0.033;
yURC = 6.45; yLLC = 0.03;
h_max = 3.0003;
dim = [xURC-xLLC yURC-yLLC h_max];

[~, axes1] = get_figure_axes(fontsize, equalaspect);
hold(axes1, 'on'); grid(axes1, 'on'); box(axes1, 'on');
trisurf(Stanza,'FaceAlpha', 0.3, 'EdgeAlpha', 0.2, 'Parent', axes1);
view([0 90]);

xlabel(axes1, 'X [m]');
ylabel(axes1, 'Y [m]');
zlabel(axes1, 'Z [m]');
hold(axes1, 'off');

%% Propagation model with 0 (pm0) and 1 reflection (pm1)
[permittivity,conductivity] = buildingMaterialPermittivity("brick",f);

pm1 = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","sbr", ...
    "AngularSeparation","low", ...
    "MaxNumReflections",1, ...
    "SurfaceMaterial","custom", ...
    "BuildingsMaterial","custom", ...
    "SurfaceMaterialConductivity",conductivity, ...
    "SurfaceMaterialPermittivity",permittivity, ...
    "BuildingsMaterialConductivity",conductivity, ...
    "BuildingsMaterialPermittivity",permittivity);

pm0 = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","sbr", ...
    "AngularSeparation","low", ...
    "MaxNumReflections",0);
% siteviewer("SceneModel",Stanza);
% show(tx,'ShowAntennaHeight', false);
% show(rx,'ShowAntennaHeight', false);
% show(RIStx,'ShowAntennaHeight', false);
% show(RISrx,'ShowAntennaHeight', false);
% pattern(tx,f); pattern(RIStx,f); 

%% Calcolo potenza ricevuta tra Tx e Rx senza RIS
tx_power_dBm_real=10*log10(tx_power_real)+30; %[W->dBm]
q0 = tx_power_dBm_real + 20*log10(lambda/(4*pi));% + maxG0;

Signal_Strength_dBm_real = sigstrength(rx,tx, ...
    'PropagationModel',pm1, ...
    'Map',Stanza,'Type','power'); %dBm
Signal_Strength_real = 10^((Signal_Strength_dBm_real-30)/10)*10^3; %mW

Tx_Rx_Distance = distance(tx, rx,'Map',Stanza); %LOS
            
% Calculates the pathloss exponent
pathloss_exponent = (q0 - Signal_Strength_real)/(20*log10(Tx_Rx_Distance));

rays = raytrace(tx,rx,pm1,"Map",Stanza);
%siteviewer("SceneModel",Stanza);
%raytrace(tx,rx,pm1);
%% Calcolo potenza ricevuta considerando la RIS (calcolo giusto perché considera solo una volta il FSPL)
%Tx_RIS_Distance = distance(tx, RISrx,'Map',Stanza);
%RISRx_Distance = distance(RIStx, rx,'Map',Stanza);
tx_power_Corr_dBm = tx_power_dBm_real + mag2db(d_T_RIS*d_RIS_R/(d_T_RIS+d_RIS_R))+mag2db(4*pi/lambda);

tx.TransmitterPower = 10^((tx_power_Corr_dBm-30)/10);
Signal_Strength_TxRIS_dBm_Corr = sigstrength(RISrx,tx, ...
    'PropagationModel',pm0, ...
    'Map',Stanza);
Signal_Strength_TxRIS_Corr = 10^((Signal_Strength_TxRIS_dBm_Corr-30)/10);
RIStx.TransmitterPower = Signal_Strength_TxRIS_Corr;
Signal_Strength_RISRx_dBm_Corr = sigstrength(rx,RIStx, ...
    'PropagationModel',pm0, ...
    'Map',Stanza);
Signal_Strength_RISRx_Corr = 10^((Signal_Strength_RISRx_dBm_Corr-30)/10)*10^3;

raysRIStx = raytrace(tx,RISrx,pm0,"Map",Stanza);
raysRISrx = raytrace(RIStx,rx,pm0,"Map",Stanza);
%raytrace(tx,RISrx,pm0);
%raytrace(RIStx,rx,pm0);
%% display risultati
disp(['Transmitted Power                    = ' , num2str(tx_power_dBm_real), ' dBm']);
disp(['Received Signal Strength Without RIS = ' , num2str(Signal_Strength_dBm_real), ' dBm']);
disp(['Received Signal Strength With RIS    = ' , num2str(Signal_Strength_RISRx_dBm_Corr), ' dBm']);

%% compute receiver power near the Rx

for k=1:N_points
    Signal_Strength_RISRx_dBm_Corr_yVariations(k) = sigstrength(rx_yVariations{k},RIStx, ...
        'PropagationModel',pm0, ...
        'Map',Stanza);
    Signal_Strength_RISRx_Corr_yVariations(k) = 10^((Signal_Strength_RISRx_dBm_Corr_yVariations(k)-30)/10)*10^3;
    Signal_Strength_RISRx_dBm_Corr_zVariations(k) = sigstrength(rx_zVariations{k},RIStx, ...
        'PropagationModel',pm0, ...
        'Map',Stanza);
    Signal_Strength_RISRx_Corr_zVariations(k) = 10^((Signal_Strength_RISRx_dBm_Corr_zVariations(k)-30)/10)*10^3;
end


% for kk=1:N_points
%     for kkk=1:N_points
%         Signal_Strength_RISRx_dBm_Corr_near(kk,kkk) = sigstrength(rx_near{kk,kkk},RIStx, ...
%             'PropagationModel',pm0, ...
%             'Map',Stanza);
%         Signal_Strength_RISRx_Corr_near(kk,kkk) = 10^((Signal_Strength_RISRx_dBm_Corr_near(kk,kkk)-30)/10)*10^3;
%     end
% end

%Compute the error beam pointing in y direction
py = polyfit(y_near,Signal_Strength_RISRx_dBm_Corr_yVariations,7);
y_near1 = linspace(RX_loc(2)-N_points/2*epsilon,RX_loc(2) +N_points/2*epsilon,1e5);
Signal_Strength_RISRx_dBm_Corr_yVariations1 = polyval(py,y_near1);

pz = polyfit(z_near,Signal_Strength_RISRx_dBm_Corr_zVariations,7);
z_near1 = linspace(RX_loc(3)-N_points/2*epsilon,RX_loc(3) +N_points/2*epsilon,1e5);
Signal_Strength_RISRx_dBm_Corr_zVariations1 = polyval(pz,z_near1);


[m1_y index_m1_y] = max(Signal_Strength_RISRx_dBm_Corr_yVariations1);
Deltay = y_near1(index_m1_y)-y_near(round(N_points/2));
Phi_Error = atand(Deltay/RX_loc(1));

[m1_z index_m1_z] = max(Signal_Strength_RISRx_dBm_Corr_zVariations1);
Deltaz = z_near1(index_m1_z)-z_near(round(N_points/2));
Theta_Error = atand(Deltaz/RX_loc(1));
%%

clear figure(1), figure(1); hold on
plot(y_near,Signal_Strength_RISRx_dBm_Corr_yVariations);
plot(y_near1,Signal_Strength_RISRx_dBm_Corr_yVariations1);
plot(RX_loc(2),Signal_Strength_RISRx_dBm_Corr,'x'); 
plot(linspace(y_near(1),y_near(end),1e5),Signal_Strength_RISRx_dBm_Corr*ones(1,1e5),'r--','LineWidth',1); 
 hold off
grid on, axis tight; xlabel('y direction [m]'), ylabel('Receiver Power [dBm]');
title('Received power vs $y$ direction','Interpreter','latex');
legend('Received Signal Strength','Power at Rx','Threeshold');


clear figure(2), figure(2); hold on
plot(z_near,Signal_Strength_RISRx_dBm_Corr_zVariations);
plot(z_near1,Signal_Strength_RISRx_dBm_Corr_zVariations1);
plot(RX_loc(3),Signal_Strength_RISRx_dBm_Corr,'x'); 
plot(linspace(z_near(1),z_near(end),1e5),Signal_Strength_RISRx_dBm_Corr*ones(1,1e5),'r--','LineWidth',1); 
 hold off
grid on, axis tight; xlabel('z direction [m]'), ylabel('Receiver Power [dBm]');
title('Received power vs $z$ direction','Interpreter','latex');
legend('Received Signal Strength','Power at Rx','Threeshold');

disp('Error of beam pointing');
disp(['Error in azimuth angle   = ' , num2str(Phi_Error), ' °']);
disp(['Error in elevation angle = ' , num2str(Theta_Error), ' °']);
% clear figure(1), figure(1); hold on
% counter = 1;
% 
% for kkk=1:N_points
%     plot(y_near,Signal_Strength_RISRx_dBm_Corr_near(:,kkk));
%     legends{counter} = sprintf('z = %.2f m',z_near(kkk));
%     counter =counter+1;
% end
% 
% plot(RX_loc(2),Signal_Strength_RISRx_dBm_Corr,'x'); legends{counter} =sprintf('Power at Rx');counter = counter+1;
% plot(linspace(RX_loc(2)-N_points/2*epsilon,RX_loc(2) +N_points/2*epsilon,1e5),Signal_Strength_RISRx_dBm_Corr*ones(1,1e5),'r--','LineWidth',1); legends{counter} =sprintf('threeshold');
% legend(legends); hold off
% grid on, axis tight; xlabel('y direction [m]'), ylabel('Receiver Power [dBm]');
% title('Received power vs $y$ direction','Interpreter','latex');
% 
% clear figure(2), figure(2); hold on
% counter = 1;
% 
% for kk=1:N_points
%     plot(z_near,Signal_Strength_RISRx_dBm_Corr_near(kk,:));
%     legends{counter} = sprintf('y = %.2f m',y_near(kk));
%     counter =counter+1;
% end
% 
% plot(RX_loc(3),Signal_Strength_RISRx_dBm_Corr,'x'); legends{counter} =sprintf('Power at Rx');counter = counter+1;
% plot(linspace(RX_loc(3)-N_points/2*epsilon,RX_loc(3) +N_points/2*epsilon,1e5),Signal_Strength_RISRx_dBm_Corr*ones(1,1e5),'r--','LineWidth',1); legends{counter} =sprintf('threeshold');
% legend(legends); hold off
% grid on, axis tight; xlabel('z direction [m]'), ylabel('Receiver Power [dBm]');
% title('Received power vs $z$ direction','Interpreter','latex');
%% |Tx-RIS Channel| 
TX_loc = TX_loc'; RX_loc = RX_loc'; RIS_loc = RIS_loc';
LOS = rays{1,1}.LineOfSight;

% i parametri sono calcolati in corrispondenza di f = 24.2 GHz
% n = Path Loss Exponent
% sigma = Shadow Fading Term (dB)
% b = Path Loss Parameter

% NLOS Parameters
n_NLOS=3.19;          % Path Loss Exponent (Indoor Office NLOS)
sigma_NLOS=abs(randn*8.29);      % Shadow Fading Term (dB) (Indoor Office NLOS)
b_NLOS=0.06;          % Path Loss Parameter (Indoor Office NLOS)
f0=24.2e9;            % Path Loss Parameter (GHz) (Indoor Office NLOS)

% LOS Parameters
n_LOS=1.73;           % Path Loss Exponent (Indoor Office LOS)
sigma_LOS=abs(randn*3.02);       % Shadow Fading Term (dB) (Indoor Office LOS)
b_LOS=0;              % Path Loss Parameter (Indoor Office NLOS)

% Generalized Element Radiation Pattern 
q=0.285;
%Gain=pi;
Gain = 10^(maxG0/10);

% Poisson distribution for cluster generation
lambda_p = 1.8; %Value referred to 28 GHz

%random phase variation
Phi_var = exp(1i*rand*2*pi);
%% STEP 1
% Compute probability LOS between T and RIS

% LOS Probability is Relatively Low for Indoors if d_T_RIS > 20
if RIS_loc(3)<TX_loc(3)   % for ground level RIS
    % InH LOS Probability
    if d_T_RIS<= 1.2
        p_LOS=1;
    elseif 1.2<d_T_RIS && d_T_RIS<6.5
        p_LOS=exp(-(d_T_RIS-1.2)/4.7);
    else
        p_LOS=0.32*exp(-(d_T_RIS-6.5)/32.6);
    end

    I_LOS=randsrc(1,1,[1,0;p_LOS 1-p_LOS]);

elseif RIS_loc(3)>=TX_loc(3) % for an RIS mounted at a high place (100% LOS)
    I_LOS=1;
end

I_LOS = 1; %Impose value equal to 1 for hypothesis

%% Calculate Tx Departure and RIS arrival angles to calculate array
% response vectors

% RIS arrival angles for LOS component
if I_LOS==1
    I_phi=sign(TX_loc(2)-RIS_loc(2));
    phi_T_RIS_LOS = I_phi* atand ( abs( RIS_loc(2)-TX_loc(2)) / abs(RIS_loc(1)-TX_loc(1)) ); %Tra RIS e Tx
   
    I_theta=sign(TX_loc(3)-RIS_loc(3));
    theta_T_RIS_LOS=I_theta * asind ( abs (RIS_loc(3)-TX_loc(3) ) / d_T_RIS );

    % Tx departure angles for LOS component
    I_phi_Tx=sign(TX_loc(2)-RIS_loc(2));
    phi_Tx_LOS = I_phi_Tx* atand ( abs( RIS_loc(2)-TX_loc(2)) / abs(RIS_loc(1)-TX_loc(1)) ); %tra Tx e RIS
    
    I_theta_Tx=sign(RIS_loc(3)-TX_loc(3));
    theta_Tx_LOS=I_theta_Tx * asind ( abs (RIS_loc(3)-TX_loc(3) )/ d_T_RIS );

    % Array Response Calculation (LOS)
    array_RIS_LOS=zeros(1,N_ele^2);
    
    counter2 = 1;
    counter3 = 1;
    for y=0:N_ele-1
        for z=0:N_ele-1
            %array_RIS_LOS(counter3)=exp(1i*k0*dx*(y*sind(theta_T_RIS_LOS) + z*sind(phi_T_RIS_LOS)*cosd(theta_T_RIS_LOS) )) ;
            array_RIS_LOS(counter3)=exp(1i*k0*dx*(y*cosd(90-theta_T_RIS_LOS) + z*sind(-phi_T_RIS_LOS)*sind(90-theta_T_RIS_LOS) )) ;
            counter3=counter3+1;
        end
    end

    array_Tx_LOS = zeros(1,numTx);
    
    for y=0:sqrt(numTx)-1
        for z=0:sqrt(numTx)-1
            array_Tx_LOS(counter2)=exp(1i*k0*dx*(z*sind(-phi_Tx_LOS)*sind(90-theta_Tx_LOS) + y*cosd(90-theta_Tx_LOS))) ;
            counter2=counter2+1;
        end
    end

    % Link Attentuation (LOS) 
    L_dB_LOS_wp=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_RIS)- sigma_LOS;
    L_LOS_wp=10^(L_dB_LOS_wp/10);

    %h_LOS_wp=sqrt(L_LOS_wp)*transpose(array_RIS_LOS)*array_Tx_LOS*exp(1i*rand*2*pi)*sqrt(Gain*(cosd(theta_T_RIS_LOS))^(2*q));
    h_LOS_wp=sqrt(L_LOS_wp)*transpose(array_RIS_LOS)*array_Tx_LOS*Phi_var*sqrt(Gain*(cosd(theta_T_RIS_LOS))^(2*q));
    %h_LOS_wp=sqrt(L_LOS_wp)*transpose(array_RIS_LOS)*array_Tx_LOS*sqrt(Gain*(cosd(theta_T_RIS_LOS))^(2*q));
else
    %h_LOS_wp = 0;
    h_LOS_wp = 0;
end

%% STEP 2
% Generate Clusters/Sub-rays, Azimuth/Elevation Departure Angles and Cluster Distances

for generate=1:100  % To ensure that at least one scatterer exist

    % Number of Clusters
    C=max([1,poissrnd(lambda_p)]);  % Poisson distributed
    %C = 0;

    % Number of Sub-rays per Cluster
    S=randi(30,1,C); % Uniformly distributed
    %S = 0;

    % Azimuth/Elevation Departure Angles
    phi_Tx=[ ];
    theta_Tx=[ ];
    phi_av=zeros(1,C);
    theta_av=zeros(1,C);
    for counter=1:C
        phi_av(counter)  = rand*180-90;     % mean azimuth departure angles of sub-rays (Laplacian distributed)
        theta_av(counter)= rand*90-45;      % mean elevation departure angles of sub-rays (Laplacian distributed)

        % cluster angles: First S(1) belongs to Cluster 1, Next S(2) belongs to Cluster 2....
        phi_Tx   = [phi_Tx,log(rand(1,S(counter))./rand(1,S(counter)))*sqrt(25/2) + phi_av(counter)];
        theta_Tx = [theta_Tx,log(rand(1,S(counter))./rand(1,S(counter)))*sqrt(25/2) + theta_av(counter)];
    end
    % Cluster Distances
    a_c=1+rand(1,C)*(d_T_RIS-1);    % Cluster distances uniform [1,d_T_RIS] 

    % Correction on Cluster Locations for Indoors
    % Room dimensions (Indoor Hotspot)
    %dim=[75,50,3.5];                  % x-y dimensions recommended by 5G Channel Model, height is assumed as 3.5 m
    
    Coordinates=zeros(C,3);          % for Clusters
    Coordinates2=zeros(sum(S),3);    % for Scatterers
    for counter=1:C
        loop = 1;
        Coordinates(counter,:)=[TX_loc(1) - a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
            TX_loc(2) - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
            TX_loc(3) + a_c(counter)*sind(theta_av(counter))] ;
        while Coordinates(counter,3)>dim(3) || Coordinates(counter,3)<0 ||  Coordinates(counter,2)>dim(2) ||  Coordinates(counter,2)<0  ||  Coordinates(counter,1)>dim(1) ||  Coordinates(counter,1)<0
            a_c(counter)=    0.8*a_c(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
            % Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
            % in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
            Coordinates(counter,:)=[TX_loc(1) - a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
                TX_loc(2) - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
                TX_loc(3) + a_c(counter)*sind(theta_av(counter))] ;
        end
    end
   
    a_c_rep=[];
    for counter3=1:C
        a_c_rep=[a_c_rep,repmat(a_c(counter3),1,S(counter3))];
    end
    for counter2=1:sum(S)
        Coordinates2(counter2,:)=[TX_loc(1) - a_c_rep(counter2)*cosd(theta_Tx(counter2))*cosd(phi_Tx(counter2)),...
            TX_loc(2) - a_c_rep(counter2)*cosd(theta_Tx(counter2))*sind(phi_Tx(counter2)),...
            TX_loc(3) + a_c_rep(counter2)*sind(theta_Tx(counter2))] ;
    end

    % Correction on Scatters
    % You may ignore the scatterers outside the walls for Indoors scenario
    ignore=[];

    for counter2=1:sum(S)
        if Coordinates(counter,3)>dim(3) || Coordinates(counter,3)<0 ||  Coordinates(counter,2)>dim(2) ||  Coordinates(counter,2)<0  ||  Coordinates(counter,1)>dim(1) ||  Coordinates(counter,1)<0
            ignore=[ignore,counter2];   % contains the indices of scatterers that can be ignored
        end
    end

    % updated indices
    indices=setdiff(1:sum(S),ignore);    % the set of active scatterer indices
    M_new=length(indices);               % number of IOs inside the room or above ground

    % Neccessary Loop to have at least one scatterer
    if M_new>0 % if M_new==0 --> all scatters are outside
        break  % break generate=1:100 % if M_new >0 we are OK (at least one scatter)
    end
end  

%% STEP 3 
% Calculate Arrival Angles for the RIS and the Link Distances
phi_cs_RIS=zeros(1,sum(S));
theta_cs_RIS=zeros(1,sum(S));
phi_Tx_cs=zeros(1,sum(S));
theta_Tx_cs=zeros(1,sum(S));
b_cs=zeros(1,sum(S));
d_cs=zeros(1,sum(S));

for counter2=indices

    b_cs(counter2)=norm(RIS_loc-Coordinates2(counter2,:));   % Distance between Scatterer and RIS
    d_cs(counter2)=a_c_rep(counter2)+b_cs(counter2);         % Total distance Tx-Scatterer-RIS

    I_phi=sign(Coordinates2(counter2,2)-RIS_loc(2));
    phi_cs_RIS(counter2)  = I_phi* atand ( abs( RIS_loc(2)-Coordinates2(counter2,2)) / abs(RIS_loc(1)-Coordinates2(counter2,1)) );
   
    I_theta=sign(Coordinates2(counter2,3)-RIS_loc(3));
    theta_cs_RIS(counter2)=I_theta * asind ( abs (RIS_loc(3)-Coordinates2(counter2,3) )/ b_cs(counter2) );

    I_phi_Tx_cs=sign(TX_loc(2)-Coordinates2(counter2,2));
    phi_Tx_cs(counter2) = I_phi_Tx_cs* atand ( abs( Coordinates2(counter2,2)-TX_loc(2)) / abs(Coordinates2(counter2,1)-TX_loc(1)) );

    I_theta_Tx_cs=sign(Coordinates2(counter2,3)-TX_loc(3));
    theta_Tx_cs(counter2)=I_theta_Tx_cs * asind ( abs (Coordinates2(counter2,3)-TX_loc(3) )/ a_c_rep(counter2) );
end

%% STEP 4
% Array Response Calculation
array_cs_RIS=zeros(sum(S),N_ele^2);

for counter2=indices
    counter3=1;
    for x=0:N_ele-1
        for y=0:N_ele-1
            %array_cs_RIS(counter2,counter3)=exp(1i*k0*dx*(x*sind(theta_cs_RIS(counter2)) + y*sind(phi_cs_RIS(counter2))*cosd(theta_cs_RIS(counter2)) )) ;
            array_cs_RIS(counter2,counter3)=exp(1i*k0*dx*(x*cosd(90-theta_cs_RIS(counter2)) + y*sind(-phi_cs_RIS(counter2))*sind(90-theta_cs_RIS(counter2)) )) ;
            counter3=counter3+1;
        end
    end
end

array_Tx_cs=zeros(sum(S),numTx);

for counter2 = indices
    counter3=1;
    for x=0:sqrt(numTx)-1
        for y = 0:sqrt(numTx)-1
            array_Tx_cs(counter2,counter3)=exp(1i*k0*dx*(x*sind(phi_Tx_cs(counter2))*cosd(theta_Tx_cs(counter2)) + y*sind(theta_Tx_cs(counter2)) )) ;
            counter3=counter3+1;
        end
    end
end
    
%% STEP 5 
% Calculate Link Attenuation and Generate Tx-RIS Channel (h) using 5G Channel Model

h_NLOS_wp=zeros(N_ele.^2,numTx);

beta=zeros(1,sum(S)); % to be reused for shared clusters 
shadow=beta;          % to be reused for shared clusters
for counter2=indices
    X_sigma=sigma_NLOS;

    %Lcs_dB_wp=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs(counter2))- X_sigma;
    Lcs_dB_wp=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs(counter2))- X_sigma;
    Lcs_wp=10^(Lcs_dB_wp/10);

    beta(counter2)=((randn+1i*randn)./sqrt(2));  % common complex gain for shared clusters
    shadow(counter2)=X_sigma;                    % commun shadow factor for shared clusters
    
    %h_NLOS_wp = h_NLOS_wp + beta(counter2)*sqrt(Gain*(cosd(theta_cs_RIS(counter2)))^(2*q))*sqrt(Lcs_wp)*transpose(array_cs_RIS(counter2,:))*array_Tx_cs(counter2,:);  % consider all scatters
    h_NLOS_wp = h_NLOS_wp + beta(counter2)*sqrt(Gain*(cosd(theta_cs_RIS(counter2)))^(2*q))*sqrt(Lcs_wp)*transpose(array_cs_RIS(counter2,:))*array_Tx_cs(counter2,:);  % consider all scatters
end

h_NLOS_wp=h_NLOS_wp.*sqrt(1/M_new);  % normalization
h_wp=h_NLOS_wp+h_LOS_wp; % include the LOS component (if any) h_LOS=0 when there is no LOS

%% STEPS 6-7 
% Generation of g (RIS-Rx Channel) - GENERATE A LOS CHANNEL

% Calculate Departure Angles Considering RIS and Rx Coordinates
%d_RIS_R=norm(RIS_loc-Rx_loc);

% Elevation Departure Angle
I_theta=sign(RX_loc(3) - RIS_loc(3));
theta_Rx_RIS=I_theta * asind( abs(RX_loc(3)-RIS_loc(3))/d_RIS_R ); % AoD of RIS

% Azimuth Departure Angle
I_phi=sign(RX_loc(2) - RIS_loc(2));
phi_Rx_RIS=I_phi * atand( abs(RX_loc(2)-RIS_loc(2))/ abs(RX_loc(1)-RIS_loc(1)) );

% AoA angles of Rx for g_LOS channel in an Indoor
phi_av_Rx  = rand*180-90;     % mean azimuth
theta_av_Rx= rand*180-90;      % mean elevation

phi_Rx   = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + phi_av_Rx];
theta_Rx = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + theta_av_Rx];

% Recalculate Array Response for Two Angles (in Rx direction)
%array_2=zeros(1,N_ele^2);
array_21=zeros(1,N_ele^2);

counter3=1;
for x=0:N_ele-1
    for y=0:N_ele-1
        %array_2(counter3)=exp(1i*k0*dx*(x*sind(theta_Rx_RIS) + y*sind(phi_Rx_RIS)*cosd(theta_Rx_RIS) )) ;
        array_2(counter3)=exp(1i*k0*dx*(x*cosd(90-theta_Rx_RIS) + y*sind(-phi_Rx_RIS)*sind(90-theta_Rx_RIS) )) ;
        counter3=counter3+1;
    end
end

array_Rx = zeros(1,numRx);

counter3=1;
for x=0:sqrt(numRx)-1
    for y=0:sqrt(numRx)-1
        array_Rx(counter3)=exp(1i*k0*dx*(x*sind(phi_Rx)*cosd(theta_Rx)+y*sind(theta_Rx)  )) ;
        counter3=counter3+1;
    end
end

% % LOS Link Attenuation

L_dB_LOS_2_wp=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_RIS_R)- sigma_LOS ;
L_LOS_2_wp=10^(L_dB_LOS_2_wp/10);

% Generate g (Pure LOS)
%g_wp = sqrt(Gain*(cosd(theta_Rx_RIS))^(2*q))*sqrt(L_LOS_2_wp)*transpose(array_2)*array_Rx*exp(1i*rand*2*pi);
g_wp = sqrt(Gain)*sqrt(L_LOS_2_wp)*transpose(array_2)*array_Rx*Phi_var;
%g_wp = sqrt(Gain)*transpose(array_2)*array_Rx*Phi_var;
%% STEP 8
% Generation of h_SISO

d_cs_tilde=zeros(1,sum(S));
h_SISO_NLOS=0;
h_SISO_NLOS_wp=0;
d_T_R = norm(TX_loc-RX_loc);

for counter2=indices

    % due to shared clusters d_cs_tilde ~ d_cs
    d_cs_tilde(counter2) = a_c_rep(counter2) + norm(Coordinates2(counter2,:)- RX_loc);

    I_phi_Tx_cs_SISO=sign(TX_loc(2)-Coordinates2(counter2,2));
    phi_Tx_cs_SISO(counter2) = I_phi_Tx_cs_SISO* atand ( abs( Coordinates2(counter2,2)-TX_loc(2)) / abs(Coordinates2(counter2,1)-TX_loc(1)) );
    
    I_theta_Tx_cs_SISO=sign(Coordinates2(counter2,3)-TX_loc(3));
    theta_Tx_cs_SISO(counter2)=I_theta_Tx_cs_SISO * asind ( abs (Coordinates2(counter2,3)-TX_loc(3) )/ a_c_rep(counter2) );

    % AoA for Rx in an Indoor
    phi_av_SISO(counter2)  = rand*180-90;     % mean azimuth
    theta_av_SISO(counter2)= rand*180-90;      % mean elevation

    phi_cs_Rx_SISO(counter2)   = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + phi_av_SISO(counter2)];
    theta_cs_Rx_SISO(counter2) = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + theta_av_SISO(counter2)];

    counter3=1;
    for x=0:sqrt(numRx)-1
        for y=0:sqrt(numRx)-1
            array_Rx_cs_SISO(counter2,counter3)=exp(1i*k0*dx*(x*sind(phi_cs_Rx_SISO(counter2))*cosd(theta_cs_Rx_SISO(counter2))+y*sind(theta_cs_Rx_SISO(counter2)) )) ;
            counter3=counter3+1;
        end
    end

    counter3=1;
    for x=0:sqrt(numTx)-1
        for y=0:sqrt(numTx)-1
            array_Tx_cs_SISO(counter2,counter3)=exp(1i*k0*dx*(x*sind(phi_Tx_cs_SISO(counter2))*cosd(theta_Tx_cs_SISO(counter2))+y*sind(theta_Tx_cs_SISO(counter2))  )) ;
            counter3=counter3+1;
        end
    end

    Lcs_dB_SISO_wp=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs_tilde(counter2))- shadow(counter2);
    Lcs_SISO_wp=10^(Lcs_dB_SISO_wp/10);

    % We consider the same complex path gain (small scale fading) with an excess phase (similar to array response)
    eta=k0* ( norm(Coordinates2(counter2,:)- RIS_loc) -  norm(Coordinates2(counter2,:)- RX_loc));

    h_SISO_NLOS_wp = h_SISO_NLOS_wp + beta(counter2)*exp(1i*eta)*sqrt(Lcs_SISO_wp)*transpose(array_Rx_cs_SISO(counter2,:))*array_Tx_cs_SISO(counter2,:);
    
end

h_SISO_NLOS_wp=h_SISO_NLOS_wp.*sqrt(1/M_new);


if RIS_loc(3) >= TX_loc(3)

    % % Include LOS component (Version 1.4)
    %     % InH LOS Probability
    if d_T_R<= 1.2
        p_LOS_3=1;
    elseif 1.2<d_T_R && d_T_R<6.5
        p_LOS_3=exp(-(d_T_R-1.2)/4.7);
    else
        p_LOS_3=0.32*exp(-(d_T_R-6.5)/32.6);
    end

    I_LOS_3=randsrc(1,1,[1,0;p_LOS_3 1-p_LOS_3]);

    % Do not recalculate, if T-RIS has LOS, we might have LOS for h_SISO as well (for ground level RIS)
    % If we have LOS for Tx-RIS, we have LOS for Tx-Rx
elseif RIS_loc(3) < TX_loc(3)  % RIS in the ground level
    I_LOS_3=I_LOS;
end
I_LOS_3 = 1;
if I_LOS_3==1

    L_SISO_LOS_dB_wp=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_R)- sigma_LOS;
    L_SISO_LOS_wp=10^(L_SISO_LOS_dB_wp/10);

    I_phi_Tx_SISO=sign(TX_loc(2)-RX_loc(2));
    phi_Tx_SISO = I_phi_Tx_SISO* atand ( abs( TX_loc(2)-RX_loc(2)) / abs(TX_loc(1)-RX_loc(1)) );

    I_theta_Tx_SISO=sign(RX_loc(3)-TX_loc(3));
    theta_Tx_SISO= I_theta_Tx_SISO* atand ( abs( RX_loc(3)-TX_loc(3)) / abs(d_T_R) );

    % AoA of Rx for Tx-Rx channel in an Indoor
    phi_av_SISO_LOS = rand*180-90;     % mean azimuth
    theta_av_SISO_LOS= rand*180-90;      % mean elevation

    phi_Rx_SISO  = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + phi_av_SISO_LOS];
    theta_Rx_SISO= [log(rand(1,1)./rand(1,1))*sqrt(25/2) + theta_av_SISO_LOS];

    counter3=1;
    for x=0:sqrt(numTx)-1
        for y=0:sqrt(numTx)-1
            array_Tx_SISO(counter3)=exp(1i*k0*dx*(x*sind(phi_Tx_SISO)*cosd(theta_Tx_SISO)+y*sind(theta_Tx_SISO)  )) ;
            counter3=counter3+1;
        end
    end
    
    counter3=1;
    for x=0:sqrt(numRx)-1
        for y=0:sqrt(numRx)-1
            array_Rx_SISO(counter3)=exp(1i*k0*dx*(x*sind(phi_Rx_SISO)*cosd(theta_Rx_SISO) +y*sind(theta_Rx_SISO))) ;
            counter3=counter3+1;
        end
    end
   
    h_SISO_LOS_wp= sqrt(L_SISO_LOS_wp)*exp(1i*rand*2*pi)*transpose(array_Rx_SISO)*array_Tx_SISO; 
end

%% consideriamo nullo il LOS
h_SISO_wp=h_SISO_NLOS_wp + h_SISO_LOS_wp; 
h_SISO_wp = 0;
%% total channel

%   g     = Vector of LOS channel coefficients between the RIS and the Rx
%  THETA  = Matrix of RIS element responses
%   h     = Vector of channel coefficients for the Tx-RIS link composed of M scatters
% h_SISO  = Characterizes the direct link channel between Tx and Rx
%   x     = Transmitted signal

PHI = Phi(:)';

THETA = zeros(N_ele^2,N_ele^2);
%THETAA = zeros(N_ele^2,N_ele^2);

for k = 1:N_ele^2
    THETA(k,k) = exp(1i.*PHI(k));
    %THETAA(k,k) = exp(1i.*PHI(k));
end

%canale considerando il valore del guadagno, ma senza potenza in
%trasmissione
%Tot_Channel_wp = g_wp' * THETA * h_wp + h_SISO_wp;
%Tot_Channel_wpp = g_wp' * THETAA * h_wp + h_SISO_wp;

Tot_Channel_wp = g_wp' * THETA * h_wp + h_SISO_wp; %consider LOS and NLOS component
%Tot_Channel_wpp1 = g_wp1' * THETAA * h_wp1 + h_SISO_wp;
Tot_Channel_LOS = g_wp' * THETA *h_LOS_wp; %consider just LOS component
%% channel model between Tx-Rx
d_cs_tilde_imag=zeros(1,sum(S));
h_SISO_NLOS_imag=0;
h_SISO_NLOS_wp_imag=0;
Rxx_loc_imag = [-RX_loc(1), RX_loc(2), RX_loc(3)];
d_T_R_imag = norm(TX_loc-Rxx_loc_imag);

for counter2=indices

    % due to shared clusters d_cs_tilde ~ d_cs
    d_cs_tilde_imag(counter2) = a_c_rep(counter2) + norm(Coordinates2(counter2,:)- Rxx_loc_imag);

    Lcs_dB_SISO_wp_imag=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs_tilde_imag(counter2))- shadow(counter2) - (rays{1}.PathLoss-mag2db(4*pi*d_T_R_imag/lambda));
    Lcs_SISO_wp_imag=10^(Lcs_dB_SISO_wp_imag/10);

    % We consider the same complex path gain (small scale fading) with an excess phase (similar to array response)
    eta_imag=k0* ( norm(Coordinates2(counter2,:)- RIS_loc) -  norm(Coordinates2(counter2,:)- Rxx_loc_imag));

    h_SISO_NLOS_wp_imag = h_SISO_NLOS_wp_imag + beta(counter2)*exp(1i*eta)*sqrt(Lcs_SISO_wp_imag)*transpose(array_Rx_cs_SISO(counter2,:))*array_Tx_cs_SISO(counter2,:);
    
end

h_SISO_NLOS_wp_imag=h_SISO_NLOS_wp_imag.*sqrt(1/M_new);



L_SISO_LOS_dB_wp_imag=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_R_imag)- sigma_LOS - (rays{1}.PathLoss-mag2db(4*pi*d_T_R_imag/lambda));
L_SISO_LOS_wp_imag=10^(L_SISO_LOS_dB_wp_imag/10);

I_phi_Tx_SISO_imag=sign(TX_loc(2)-Rxx_loc_imag(2));
phi_Tx_SISO_imag = I_phi_Tx_SISO_imag* atand ( abs( TX_loc(2)-Rxx_loc_imag(2)) / abs(TX_loc(1)-Rxx_loc_imag(1)) );

I_theta_Tx_SISO_imag=sign(Rxx_loc_imag(3)-TX_loc(3));
theta_Tx_SISO_imag= I_theta_Tx_SISO_imag* atand ( abs( Rxx_loc_imag(3)-TX_loc(3)) / abs(d_T_R_imag) );

% AoA of Rx for Tx-Rx channel in an Indoor
phi_av_SISO_LOS = rand*180-90;     % mean azimuth
theta_av_SISO_LOS= rand*180-90;      % mean elevation

h_SISO_LOS_wp_imag= sqrt(L_SISO_LOS_wp_imag)*Phi_var*transpose(array_Rx_SISO)*array_Tx_SISO;

h_SISO_wp_imag=h_SISO_NLOS_wp_imag + h_SISO_LOS_wp_imag; 

%% Modulation random bit-stream 64QAM/OFDM

scs = 30;
carrier = nrCarrierConfig('SubcarrierSpacing',scs);
p = 1;
pdsch = nrPDSCHConfig('NumLayers',p,'Modulation','64QAM');
[ind,info] = nrPDSCHIndices(carrier,pdsch);
numDataBits = info.G;
cws = randi([0 1],numDataBits,1);
% save("streamBit.mat","cws");
% cws = load("streamBit.mat");
% cws = cws.cws;
sym = nrPDSCH(carrier,pdsch,cws,'OutputDataType','single');
txGrid = nrResourceGrid(carrier,p);
txGrid(ind) = sym;
initialNSlot = carrier.NSlot;
cpl = 'extended';
[txWaveform,Info] = nrOFDMModulate(txGrid,scs,initialNSlot,'CyclicPrefix',cpl);

%% Rx signal after circular convolution between TxSignal and impulse response of channel
waveform_wp = conv2(txWaveform,Tot_Channel_wp,'same'); 
waveform_wp_LOS = conv2(txWaveform,Tot_Channel_LOS,'same');
waveform_wp_imag = conv2(txWaveform,h_SISO_wp_imag,'same'); 

%% Tx Bit-stream  
nrb = carrier.NSizeGrid;
bitTx = qamdemod(txGrid,64,'OutputType','bit');

%% Demodulation OFDM of RxSignal

grid_wp = nrOFDMDemodulate(waveform_wp,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_wp_LOS = nrOFDMDemodulate(waveform_wp_LOS,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_wp_imag = nrOFDMDemodulate(waveform_wp_imag,nrb,scs,initialNSlot,'CyclicPrefix',cpl);

%% Demodulation QAM of RxSignal 
bitRx_wp = qamdemod(grid_wp,64,'OutputType','bit');
bitRx_wp_LOS = qamdemod(grid_wp_LOS,64,'OutputType','bit');

bitRx_wp_imag = qamdemod(grid_wp_imag,64,'OutputType','bit');

%% difference between Tx and Rx bit
%Without TxPower, only Gain

[number_wp,ratio_wp] = biterr(bitTx,bitRx_wp);
BitErrorRate_wp = ratio_wp*100; %percentage of error bits
[number_wp_LOS,ratio_wp_LOS] = biterr(bitTx,bitRx_wp_LOS);
BitErrorRate_wp_LOS = ratio_wp_LOS*100; %percentage of error bits

[number_wp_imag,ratio_wp_imag] = biterr(bitTx,bitRx_wp_imag);
BitErrorRate_wp_imag = ratio_wp_imag*100; %percentage of error bits

disp('CASE 3: Channel model with clusters, scatters and random phenomena')
disp(['-> Bit Error Rate with RIS/conv  = ', num2str(BitErrorRate_wp), ' %. ']);
disp(['-> Bit Error Rate with RIS LOS  = ', num2str(BitErrorRate_wp_LOS), ' %. ']);

disp(['-> Bit Error Rate without RIS = ' , num2str(BitErrorRate_wp_imag), ' %. ']);
%% Integration of Tw Power in the model
%Per calcolare la potenza in ricezione espressa in dBm conoscendo la potenza in trasmissione
% espressa in dBm e la funzione di trasferimento del canale, è necessario considerare il guadagno
% o l'attenuazione introdotti dal canale.
% La potenza in ricezione può essere calcolata utilizzando la seguente formula:
%P_rx(dBm) = P_tx(dBm) + G_channel(dB)
%in which
%P_rx(dBm) = received power in dBm
%P_tx(dBm) = transmitted power in dBm
%G_channel(dB) = gain or attenuation of the channel in dB
%G_channel(dB) = 10log10(|H(f)|^2)
%with |H(f)| transfer function of the channel at frequency f

%TOT_CHANNEL = fft(Tot_Channel_wp1); %Transfer function of channel
Pw_dBm_H = 10*log10(abs(Tot_Channel_wp).^2);
Pw_dBm_H_LOS = 10*log10(abs(Tot_Channel_LOS).^2);
%Pw_dBm_H = 10*log10(abs(TOT_CHANNEL).^2); %Power associated to the channel (pag. 193 BOOK -> Teoria dei segnali)
%Pw_dBm_rx = tx_power_dBm_real + Pw_dBm_H;
%tx_power_Corr_dBm1 = tx_power_dBm_real + mag2db(d_T_RIS*d_RIS_R/(d_T_RIS+d_RIS_R))+mag2db(4*pi/lambda);
Pw_dBm_rx = tx_power_Corr_dBm + Pw_dBm_H;
Pw_dBm_rx_LOS = tx_power_Corr_dBm + Pw_dBm_H_LOS;
% TOT_CHANNEL_imag = fft(h_SISO_wp_imag); %Transfer function of channel
Pw_dBm_H_imag = 10*log10(abs(h_SISO_wp_imag).^2);
%Pw_dBm_H_imag = 10*log10(abs(TOT_CHANNEL_imag).^2); %Power associated to the channel (pag. 193 BOOK -> Teoria dei segnali)
Pw_dBm_rx_imag = tx_power_dBm_real + Pw_dBm_H_imag;

disp('Received power considering transmitted power');
disp(['Without RIS = ' , num2str(Pw_dBm_rx_imag), 'dBm']);
disp(['With RIS = ' , num2str(Pw_dBm_rx), 'dBm']);
disp(['With RIS LOS= ' , num2str(Pw_dBm_rx_LOS), 'dBm']);
%%
function [figureh, axesh] = get_figure_axes (fontsize, equalflag)
figureh = figure('Color', [1 1 1]);
axesh = axes('Parent',figureh);
    if equalflag
        set(axesh, 'DataAspectRatio', [1 1 1], 'FontSize', fontsize, 'FontWeight', 'bold',...
            'LineWidth', 2);
    else
        set(axesh, 'FontSize', 13, 'FontWeight', 'bold','LineWidth', 2);
    end

end
