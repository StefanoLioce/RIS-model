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
N_ele = 8; %sqrt of number of elemets - even number

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
% Ei = 1; %incident electric field
% Er = zeros(181,361); %initial value of reflected electric field
% Gamma = 1; %reflection amplitude
% 
% A = 1; %illuminating amplitude
% a = pi/6; %illuminating phase
% PHI = linspace(-180,180,361);
% THETA = linspace(-90,90,181);
% clear z, clear y;
% Phase_element = zeros(N_ele, N_ele);
% for yy=1:N_ele
%     for zz=1:N_ele
%         Tx_loc(zz,yy) = {[TX_loc(1); -(yy-1)*dx+TX_loc(2); (zz-1)*dy+TX_loc(3)]};
%         Rx_loc(zz,yy) = {[RX_loc(1); -(yy-1)*dx+RX_loc(2); (zz-1)*dy+RX_loc(3)]};
% 
%         phi_Txx(zz,yy) = sign((-(yy-1)*dx+TX_loc(2))-RIS_loc(2))*atand(abs(RIS_loc(2)-(-(yy-1)*dx+TX_loc(2)))/abs(RIS_loc(1)-TX_loc(1)));
%         theta_Txx(zz,yy) =  sign(((zz-1)*dy+TX_loc(3))-RIS_loc(3))*asind(abs(RIS_loc(3)-((zz-1)*dy+TX_loc(3)))/norm([TX_loc(1); -(yy-1)*dx+TX_loc(2); (zz-1)*dy+TX_loc(3)]-RIS_loc ));
% 
%         phi_Rxx(zz,yy) = sign((-(yy-1)*dx+RX_loc(2))-RIS_loc(2))*atand(abs((-(yy-1)*dx+RX_loc(2))-RIS_loc(2))/abs(RX_loc(1)-RIS_loc(1)));
%         theta_Rxx(zz,yy) =  sign(((zz-1)*dy+RX_loc(3))-RIS_loc(3))*asind(abs(((zz-1)*dy+RX_loc(3))-RIS_loc(3))/norm(RIS_loc-[RX_loc(1); -(yy-1)*dx+RX_loc(2); (zz-1)*dy+RX_loc(3)]));
%     end
% end
%%

Phase_element = zeros(N_ele,N_ele); %initialization of continuous phase matrix values
Phase_elementq = zeros(N_ele,N_ele); %initialization of 1 bit phase matrix values
Phase_elementqq = zeros(N_ele,N_ele); %initialization of2 bits phase matrix values
Azimuth = linspace(-90,90,181)';
Elevation = linspace(-90,90,181)';
Ei = zeros(length(Azimuth),length(Elevation)); %Array factor without phase 
Ei1 = zeros(length(Azimuth),length(Elevation)); %Array factor for RX location with continuous phase values
Ei1q = zeros(length(Azimuth),length(Elevation)); %Array factor for RX location with 1 bit phase values
Ei1qq = zeros(length(Azimuth),length(Elevation)); %Array factor for RX location with 2 bits phase values
pp = 0;
for az = 1:length(Azimuth)
    for el = 1:length(Elevation)
        for z = 1:N_ele
            for y = 1:N_ele
                clear pp;
                %Phi = k0*((x-1)*dx*sind(Elevation(el))*sind(Azimuth(az))+(y-1)*dx*sind(Elevation(el))*cosd(Azimuth(az)));
                pp = k0*((z-1)*dx*sind(90-Elevation(el))*sind(-Azimuth(az))+(y-1)*dx*cosd(90-Elevation(el)));
                while pp>=2*pi
                    pp=pp-2*pi;
                end
                while pp<=-2*pi
                    pp=pp+2*pi;
                end
                while pp <0
                    pp = pp +2*pi;
                end

                zn = (z-1)*dx; ym = (y-1)*dx;
                Phase_element(z,y) = k0*(zn*sind(90-theta_RX)*sind(-phi_RX)+ym*cosd(90-theta_RX));
                while Phase_element(z,y)>=2*pi
                    Phase_element(z,y)=Phase_element(z,y)-2*pi;
                end
                while Phase_element(z,y)<=-2*pi
                    Phase_element(z,y)=Phase_element(z,y)+2*pi;
                end
                if Phase_element(z,y) < 0
                    Phase_element(z,y) = Phase_element(z,y) + 2*pi;
                end

                if Phase_element(z,y)>pi/2 && Phase_element(z,y)<3*pi/2
                    Phase_elementq(z,y) = pi;
                else
                    Phase_elementq(z,y) = 0;
                end

                if Phase_element(z,y) >= pi/4 && Phase_element(z,y) <= 3/4*pi
                    Phase_elementqq(z,y) = pi/2;
                elseif Phase_element(z,y) > 5/4*pi && Phase_element(z,y)<=3*pi/2
                    Phase_elementqq(z,y) = 3*pi/2;
                elseif Phase_element(z,y)>3/4*pi && Phase_element(z,y) <=5/4*pi
                    Phase_elementqq(z,y) = pi;
                else
                    Phase_elementqq(z,y) = 0;
                end


                
                Ei(az,el) = Ei(az,el) + (exp(-1i*pp));
                Ei1(az,el) = Ei1(az,el) + (exp(-1i*pp).*exp(1i*Phase_element(z,y)));
                Ei1q(az,el) = Ei1q(az,el) + (exp(-1i*pp).*exp(1i*Phase_elementq(z,y)));
                Ei1qq(az,el) = Ei1qq(az,el) + (exp(-1i*pp).*exp(1i*Phase_elementqq(z,y)));
            end
        end

        %Ei(az,el) = Ei(az,el)*cos(Elevation(el));
       % Ei1(az,el) = Ei1(az,el)*cos(Elevation(el));
     
        
    end
end
E = (abs(10*log10((Ei))));

Ei1 = Ei1./max(max(Ei));
Ei1q = Ei1q./max(max(Ei));
Ei1qq = Ei1qq./max(max(Ei));
Ei = Ei./max(max(Ei));

%plot(Elevation,Ei(1,:))
figure(); 
subplot(2,1,1), plot(Azimuth,(Ei(:,ceil(length(Elevation)/2))),'LineWidth',3), grid on, axis tight, title('Elevation Cut, el = 0°'); xlabel('Azimuth angle'), ylabel('Normalized Electric Field');
subplot(2,1,2), plot(Elevation,(Ei(ceil(length(Azimuth)/2),:)),'LineWidth',3), grid on, axis tight, title('Azimuth Cut, az = 0°'); xlabel('Elevation angle'), ylabel('Normalized Electric Field');
%hold on, plot(Elevation,Ei1(91,:));

figure(); 
subplot(2,1,1);
plot(Azimuth,Ei1(:,ceil(length(Elevation)/2)+round(theta_RX)),'LineWidth',3);
hold on;
plot(Azimuth,Ei1q(:,ceil(length(Elevation)/2)+round(theta_RX)),'LineWidth',3);
plot(Azimuth,Ei1qq(:,ceil(length(Elevation)/2)+round(theta_RX)),'LineWidth',3);
hold off, grid on, axis tight, title('Elevation Cut, el = 0°');
xlabel('Azimuth angle'), ylabel('Normalized Electric Field');
legend('Continuous phase values','1 bit quantized phase values','2 bit quantized phase values');

subplot(2,1,2);
plot(Elevation,Ei1(ceil(length(Azimuth)/2)+round(phi_RX),:),'LineWidth',3); 
hold on;
plot(Elevation,Ei1q(ceil(length(Azimuth)/2)+round(phi_RX),:),'LineWidth',3);
plot(Elevation,Ei1qq(ceil(length(Azimuth)/2)+round(phi_RX),:),'LineWidth',3);
hold off, grid on, axis tight, title('Azimuth Cut, az = 0°'); 
xlabel('Elevation angle'), ylabel('Normalized Electric Field');
legend('Continuous phase values','1 bit quantized phase values','2 bit quantized phase values');
%%
disp('Received position : ');
disp(['Azimuth angle     : -> ',num2str(phi_RX),' °']);
disp(['Elevation angle   : -> ',num2str(theta_RX),' °']);
for z = 1:N_ele
    for y = 1:N_ele

        zn = (z-1)*dy; ym = (y-1)*dy;
      
        Phi(z,y) = k0*(ym*sind(90-theta_RX)*sind(-phi_RX)+zn*cosd(90-theta_RX));
        %Phii(z,y) = -k0*(ym*sind(90-theta_RX)*sind(-phi_RX)+zn*cosd(90-theta_RX)) + k0*norm([RIS_loc(1) ; RIS_loc(2)+ym ; RIS_loc(3)-zn]  -TX_loc);
        Phii(z,y) = k0*(ym*sind(90-theta_TX)*sind(-phi_TX)+zn*cosd(90-theta_TX));%+ k0*(ym*sind(90-theta_TX)*sind(-phi_TX)+zn*cosd(90-theta_TX));
         
        while Phi(z,y)>=2*pi
            Phi(z,y)=Phi(z,y)-2*pi;
        end
        while Phi(z,y)<=-2*pi
            Phi(z,y)=Phi(z,y)+2*pi;
        end
        if Phi(z,y) < 0
            Phi(z,y) = Phi(z,y) + 2*pi;
        end

        while Phii(z,y)>=2*pi
            Phii(z,y)=Phii(z,y)-2*pi;
        end
        while Phii(z,y)<=-2*pi
            Phii(z,y)=Phii(z,y)+2*pi;
        end
        if Phii(z,y) < 0
            Phii(z,y) = Phii(z,y) + 2*pi;
        end

        if Phi(z,y)>pi/2 && Phi(z,y)<3*pi/2
            Phiq(z,y) = pi;
        else
            Phiq(z,y) = 0;
        end

        if Phi(z,y) >= pi/4 && Phi(z,y) <= 3/4*pi
            Phiqq(z,y) = pi/2;
        elseif Phi(z,y) > 5/4*pi && Phi(z,y)<=3*pi/2
            Phiqq(z,y) = 3*pi/2;
        elseif Phi(z,y)>3/4*pi && Phi(z,y) <=5/4*pi
            Phiqq(z,y) = pi;
        else
            Phiqq(z,y) = 0;
        end


        if Phii(z,y)>pi/2 && Phii(z,y)<=3*pi/2
            Phiiq(z,y) = pi;
        else
            Phiiq(z,y) = 0;
        end
        

         if Phii(z,y) > pi/4 && Phii(z,y) <= 3/4*pi
            Phiiqq(z,y) = pi/2;
        elseif Phii(z,y) > 5/4*pi && Phii(z,y)<=3*pi/2
            Phiiqq(z,y) = 3*pi/2;
        elseif Phii(z,y)>3/4*pi && Phii(z,y) <=5/4*pi
            Phiiqq(z,y) = pi;
        else
            Phiiqq(z,y) = 0;
        end

        % if Phii(z,y)>=-pi/2 && Phii(z,y)<=pi/2
        %     Phiiq(z,y) = 0;
        % else
        %     Phiiq(z,y) = pi;
        % end
        
    end
end

%%
RIS = phased.URA("Element",ant, ...
    "Size", [N_ele N_ele], ...
    "ElementSpacing",[dx dy],"ArrayNormal","x");

% RIS transmitter part
taperedRISt = clone(RIS); % RIStx with continuous phase
steeringVector_t = phased.SteeringVector("SensorArray",taperedRISt,'IncludeElementResponse',true);

taperedRIStt = clone(RIS); % RISrx with continuous phase
steeringVector_tt = phased.SteeringVector("SensorArray",taperedRIStt,'IncludeElementResponse',true);

taperedRIStq = clone(RIS); %RIStx 1 bit quantized phase
taperedRIStqq = clone(RIS); %RIStx 2 bit quantized phase
steeringVectorq_t = phased.SteeringVector("SensorArray",taperedRIStq,'IncludeElementResponse',true);
steeringVectorqq_t = phased.SteeringVector("SensorArray",taperedRIStqq,'IncludeElementResponse',true);

taperedRISttq = clone(RIS); %RISrx 1 bit quantized phase
taperedRISttqq = clone(RIS); %RISrx 2 bit quantized phasenrap(
steeringVectorq_tt = phased.SteeringVector("SensorArray",taperedRISttq,'IncludeElementResponse',true);
steeringVectorqq_tt = phased.SteeringVector("SensorArray",taperedRISttqq,'IncludeElementResponse',true);


%RIS receiver part
%taperedRISr = clone(RIS);
%taperedRISt1 = clone(RIS);
% startTaper_r = taperedRISr.Taper;
% steeringVector_r = phased.SteeringVector("SensorArray",taperedRISr);
%figure(), viewArray(RIS,'ShowTaper',true);
%% 
%tap = ones(N_ele, N_ele);
%taperedRISt.Taper = exp(1i*Phi);
taperedRISt.Taper = exp(1i*Phi); %RIStx continuous phase matrix 
taperedRIStt.Taper = exp(1i*Phii); %RISrx continuous phase

taperedRIStq.Taper = exp(1i*Phiq); %RIStx 1 bit quantized phase matrix
taperedRIStqq.Taper = exp(1i*Phiqq); %RIStx 2 bit quantized phase matrix

taperedRISttq.Taper = exp(1i*Phiiq); %RISrx 1 bit quantized phase matrix
taperedRISttqq.Taper = exp(1i*Phiiqq); %RISrx 2 bit quantized phase matrix

%taperedRISt1.Taper = exp(1i*Phi1);
%figure(), viewArray(taperedRISt,'ShowTaper',true);

%taperedRISt1.Taper = ones.*exp(1i*Phi1);

%Radiation pattern with continuous phase value
%figure(), pattern(taperedRISt,f,-180:1:180,-90:1:90); %3D representation
% figure(), pattern(taperedRISt,f,0,-90:1:90); %azimuth cut ->shows elevation angle
% figure(), pattern(taperedRISt,f,-180:1:180,0); %elevation cut ->shows azimuth angle
[BW, ang] = beamwidth(RIS,f,'dbDown',3); % -3 dB beamwidth 
maxG0Array = max(max(pattern(taperedRISt,f,'Type','directivity')));
maxG0Arrayq = max(max(pattern(taperedRIStq,f,'Type','directivity')));
maxG0Arrayqq = max(max(pattern(taperedRIStqq,f,'Type','directivity')));
Gq = (maxG0Arrayq/(N_ele^2))/(maxG0Array/(N_ele^2));
Gqq = (maxG0Arrayqq/(N_ele^2))/(maxG0Array/(N_ele^2));
%Radiation pattern with quantized phase value
%figure(), pattern(taperedRIStq,f,-180:1:180,-90:1:90);
% figure(), pattern(taperedRIStq,f,0,-90:1:90);
% figure(), pattern(taperedRIStq,f,-180:1:180,0);
[BWq, angq] = beamwidth(taperedRIStq,f,'dbDown',3);

%figure(), pattern(taperedRIStqq,f,-180:1:180,-90:1:90);

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
RIStxq = txsite("cartesian", ...
    "Antenna",taperedRIStq, ...
    "AntennaPosition",RIS_loc, ...
    'TransmitterFrequency',f);
RIStxqq = txsite("cartesian", ...
    "Antenna",taperedRIStqq, ...
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
% epsilon = 0.01; %[m] distance between evaluation points
% N_points = 7;
% %x_near = linspace(RX_loc(1)-N_points/2*epsilon,RX_loc(1) +N_points/2*epsilon,N_points);
% y_near = linspace(RX_loc(2)-N_points/2*epsilon,RX_loc(2) +N_points/2*epsilon,N_points);
% z_near = linspace(RX_loc(3)-N_points/2*epsilon,RX_loc(3) +N_points/2*epsilon,N_points);
% for k=1:N_points
%     RX_loc_yVariations = [RX_loc(1) y_near(k) RX_loc(3)];
%     RX_loc_zVariations = [RX_loc(1) RX_loc(2) z_near(k)];
%     rx_yVariations{k} = rxsite("cartesian", ...
%         "Antenna",antenna, ...
%         "AntennaPosition",RX_loc_yVariations', ...
%         "AntennaAngle",[0;90],"ReceiverSensitivity",-150);
%     rx_zVariations{k} = rxsite("cartesian", ...
%         "Antenna",antenna, ...
%         "AntennaPosition",RX_loc_zVariations', ...
%         "AntennaAngle",[0;90],"ReceiverSensitivity",-150);
% end



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
x1_obs = 3.02; x2_obs = 6.42;
y1_obs = 2.94; y2_obs = 3.18;
z1_obs = 0; z2_obs = h_max;
%% Propagation model with 0 (pm0) and 1 reflection (pm1)
[permittivity,conductivity] = buildingMaterialPermittivity("brick",f);
epsilon0 = 8.854187817*10^-12;
mi0 = 4*pi*10^-7;
Z1 = sqrt(mi0/epsilon0);
Z2 = sqrt(mi0/epsilon0*permittivity);
Gamma_Wall = (Z2-Z1)/(Z2+Z1);
Attenuation_Wall = -20*log10(abs(Gamma_Wall));
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
% siteviewer("SceneModel",Stanza);
% show(tx,'ShowAntennaHeight', false);
% show(rx,'ShowAntennaHeight', false);
% show(RIStx,'ShowAntennaHeight', false);
% show(RISrx,'ShowAntennaHeight', false); pattern(RIStxq,f);
% siteviewer("SceneModel",Stanza);
% show(tx,'ShowAntennaHeight', false);
% show(rx,'ShowAntennaHeight', false);
% show(RIStx,'ShowAntennaHeight', false);
% show(RISrx,'ShowAntennaHeight', false); pattern(RIStxqq,f);

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
%% quantized

RIStxq.TransmitterPower = Signal_Strength_TxRIS_Corr;
Signal_Strength_RISRx_dBm_Corrq = sigstrength(rx,RIStxq, ...
    'PropagationModel',pm0, ...
    'Map',Stanza);
Signal_Strength_RISRx_Corrq = 10^((Signal_Strength_RISRx_dBm_Corrq-30)/10)*10^3;

RIStxqq.TransmitterPower = Signal_Strength_TxRIS_Corr;
Signal_Strength_RISRx_dBm_Corrqq = sigstrength(rx,RIStxqq, ...
    'PropagationModel',pm0, ...
    'Map',Stanza);
Signal_Strength_RISRx_Corrqq = 10^((Signal_Strength_RISRx_dBm_Corrqq-30)/10)*10^3;
%raytrace(tx,RISrx,pm0);
%raytrace(RIStx,rx,pm0);
%% display risultati
disp(['Transmitted Power                                    = ' , num2str(tx_power_dBm_real), ' dBm']);
disp(['Received Signal Strength Without RIS                 = ' , num2str(Signal_Strength_dBm_real), ' dBm']);
disp(['Received Signal Strength With RIS                    = ' , num2str(Signal_Strength_RISRx_dBm_Corr), ' dBm']);
disp(['Received Signal Strength With 1 bit quantized RIS    = ' , num2str(Signal_Strength_RISRx_dBm_Corrq), ' dBm']);
disp(['Received Signal Strength With 2 bit quantized RIS    = ' , num2str(Signal_Strength_RISRx_dBm_Corrqq), ' dBm']);
%% compute receiver power near the Rx

% for k=1:N_points
%     Signal_Strength_RISRx_dBm_Corr_yVariations(k) = sigstrength(rx_yVariations{k},RIStx, ...
%         'PropagationModel',pm0, ...
%         'Map',Stanza);
%     Signal_Strength_RISRx_Corr_yVariations(k) = 10^((Signal_Strength_RISRx_dBm_Corr_yVariations(k)-30)/10)*10^3;
%     Signal_Strength_RISRx_dBm_Corr_zVariations(k) = sigstrength(rx_zVariations{k},RIStx, ...
%         'PropagationModel',pm0, ...
%         'Map',Stanza);
%     Signal_Strength_RISRx_Corr_zVariations(k) = 10^((Signal_Strength_RISRx_dBm_Corr_zVariations(k)-30)/10)*10^3;
% end

% for k=1:N_points
%     Signal_Strength_RISRx_dBm_Corr_yVariationsq(k) = sigstrength(rx_yVariations{k},RIStxq, ...
%         'PropagationModel',pm0, ...
%         'Map',Stanza);
%     Signal_Strength_RISRx_Corr_yVariationsq(k) = 10^((Signal_Strength_RISRx_dBm_Corr_yVariationsq(k)-30)/10)*10^3;
%     Signal_Strength_RISRx_dBm_Corr_zVariationsq(k) = sigstrength(rx_zVariations{k},RIStxq, ...
%         'PropagationModel',pm0, ...
%         'Map',Stanza);
%     Signal_Strength_RISRx_Corr_zVariationsq(k) = 10^((Signal_Strength_RISRx_dBm_Corr_zVariationsq(k)-30)/10)*10^3;
% end

% for k=1:N_points
%     Signal_Strength_RISRx_dBm_Corr_yVariationsq(k) = sigstrength(rx_yVariations{k},RIStxq, ...
%         'PropagationModel',pm0, ...
%         'Map',Stanza);
%     Signal_Strength_RISRx_Corr_yVariationsq(k) = 10^((Signal_Strength_RISRx_dBm_Corr_yVariationsq(k)-30)/10)*10^3;
%     Signal_Strength_RISRx_dBm_Corr_zVariationsq(k) = sigstrength(rx_zVariations{k},RIStxq, ...
%         'PropagationModel',pm0, ...
%         'Map',Stanza);
%     Signal_Strength_RISRx_Corr_zVariationsq(k) = 10^((Signal_Strength_RISRx_dBm_Corr_zVariationsq(k)-30)/10)*10^3;
% end

% for kk=1:N_points
%     for kkk=1:N_points
%         Signal_Strength_RISRx_dBm_Corr_near(kk,kkk) = sigstrength(rx_near{kk,kkk},RIStx, ...
%             'PropagationModel',pm0, ...
%             'Map',Stanza);
%         Signal_Strength_RISRx_Corr_near(kk,kkk) = 10^((Signal_Strength_RISRx_dBm_Corr_near(kk,kkk)-30)/10)*10^3;
%     end
% end

%Compute the error beam pointing in y direction
% py = polyfit(y_near,Signal_Strength_RISRx_dBm_Corr_yVariations,7);
% y_near1 = linspace(RX_loc(2)-N_points/2*epsilon,RX_loc(2) +N_points/2*epsilon,1e5);
% Signal_Strength_RISRx_dBm_Corr_yVariations1 = polyval(py,y_near1);
% pz = polyfit(z_near,Signal_Strength_RISRx_dBm_Corr_zVariations,7);
% z_near1 = linspace(RX_loc(3)-N_points/2*epsilon,RX_loc(3) +N_points/2*epsilon,1e5);
% Signal_Strength_RISRx_dBm_Corr_zVariations1 = polyval(pz,z_near1);
% [m1_y index_m1_y] = max(Signal_Strength_RISRx_dBm_Corr_yVariations1);
% Deltay = y_near1(index_m1_y)-y_near(round(N_points/2));
% y_Rx_opt = RX_loc(2) + Deltay;
% if phi_RX <0
%     Phi_Error = phi_RX + sign(y_Rx_opt - RIS_loc(2))*atand((y_Rx_opt - RIS_loc(2))/(RX_loc(1)-RIS_loc(1)));
% else
%     Phi_Error = phi_RX - sign(y_Rx_opt - RIS_loc(2))*atand((y_Rx_opt - RIS_loc(2))/(RX_loc(1)-RIS_loc(1)));
% end
% [m1_z index_m1_z] = max(Signal_Strength_RISRx_dBm_Corr_zVariations1);
% Deltaz = z_near1(index_m1_z)-z_near(round(N_points/2));
% z_Rx_opt = RX_loc(3) + Deltaz;
% if theta_RX <0
%     Theta_Error = theta_RX + sign(z_Rx_opt-RIS_loc(3))*asind((z_Rx_opt-RIS_loc(3))/norm(RIS_loc - [RX_loc(1); RX_loc(2); z_Rx_opt]));
% else
%     Theta_Error = theta_RX - sign(z_Rx_opt-RIS_loc(3))*asind((z_Rx_opt-RIS_loc(3))/norm(RIS_loc - [RX_loc(1); RX_loc(2); z_Rx_opt]));
% end



% pyq = polyfit(y_near,Signal_Strength_RISRx_dBm_Corr_yVariationsq,7);
% y_near1q = linspace(RX_loc(2)-N_points/2*epsilon,RX_loc(2) +N_points/2*epsilon,1e5);
% Signal_Strength_RISRx_dBm_Corr_yVariations1q = polyval(pyq,y_near1q);
% pzq = polyfit(z_near,Signal_Strength_RISRx_dBm_Corr_zVariationsq,7);
% z_near1q = linspace(RX_loc(3)-N_points/2*epsilon,RX_loc(3) +N_points/2*epsilon,1e5);
% Signal_Strength_RISRx_dBm_Corr_zVariations1q = polyval(pzq,z_near1q);
% [m1_yq index_m1_yq] = max(Signal_Strength_RISRx_dBm_Corr_yVariations1q);
% Deltayq = y_near1q(index_m1_yq)-y_near(round(N_points/2));
% y_Rx_optq = RX_loc(2) + Deltayq;
% if phi_RX <0
%     Phi_Errorq = phi_RX + sign(y_Rx_optq - RIS_loc(2))*atand((y_Rx_optq - RIS_loc(2))/(RX_loc(1)-RIS_loc(1)));
% else
%     Phi_Errorq = phi_RX - sign(y_Rx_optq - RIS_loc(2))*atand((y_Rx_optq - RIS_loc(2))/(RX_loc(1)-RIS_loc(1)));
% end
% [m1_zq index_m1_zq] = max(Signal_Strength_RISRx_dBm_Corr_zVariations1q);
% Deltazq = z_near1q(index_m1_zq)-z_near(round(N_points/2));
% z_Rx_optq = RX_loc(3) + Deltazq;
% if theta_RX <0
%     Theta_Errorq = theta_RX + sign(z_Rx_optq-RIS_loc(3))*asind((z_Rx_optq-RIS_loc(3))/norm(RIS_loc - [RX_loc(1); RX_loc(2); z_Rx_optq]));
% else
%     Theta_Errorq = theta_RX - sign(z_Rx_optq-RIS_loc(3))*asind((z_Rx_optq-RIS_loc(3))/norm(RIS_loc - [RX_loc(1); RX_loc(2); z_Rx_optq]));
% end



%Compute the error beam pointing in y direction
% pyqq = polyfit(y_near,Signal_Strength_RISRx_dBm_Corr_yVariationsqq,7);
% y_near1qq = linspace(RX_loc(2)-N_points/2*epsilon,RX_loc(2) +N_points/2*epsilon,1e5);
% Signal_Strength_RISRx_dBm_Corr_yVariations1qq = polyval(pyqq,y_near1qq);
% 
% pzqq = polyfit(z_near,Signal_Strength_RISRx_dBm_Corr_zVariationsqq,7);
% z_near1qq = linspace(RX_loc(3)-N_points/2*epsilon,RX_loc(3) +N_points/2*epsilon,1e5);
% Signal_Strength_RISRx_dBm_Corr_zVariations1qq = polyval(pzqq,z_near1qq);
% 
% 
% [m1_yqq index_m1_yqq] = max(Signal_Strength_RISRx_dBm_Corr_yVariations1qq);
% Deltayqq = y_near1qq(index_m1_yqq)-y_near(round(N_points/2));
% y_Rx_optqq = RX_loc(2) + Deltayqq;
% if phi_RX <0
%     Phi_Errorqq = phi_RX + sign(y_Rx_optqq - RIS_loc(2))*atand((y_Rx_optqq - RIS_loc(2))/(RX_loc(1)-RIS_loc(1)));
% else
%     Phi_Errorqq = phi_RX - sign(y_Rx_optqq - RIS_loc(2))*atand((y_Rx_optqq - RIS_loc(2))/(RX_loc(1)-RIS_loc(1)));
% end
% [m1_zqq index_m1_zqq] = max(Signal_Strength_RISRx_dBm_Corr_zVariations1qq);
% Deltazqq = z_near1qq(index_m1_zqq)-z_near(round(N_points/2));
% z_Rx_optqq = RX_loc(3) + Deltazqq;
% if theta_RX <0
%     Theta_Errorqq = theta_RX + sign(z_Rx_optqq-RIS_loc(3))*asind((z_Rx_optqq-RIS_loc(3))/norm(RIS_loc - [RX_loc(1); RX_loc(2); z_Rx_optqq]));
% else
%     Theta_Errorqq = theta_RX - sign(z_Rx_optqq-RIS_loc(3))*asind((z_Rx_optqq-RIS_loc(3))/norm(RIS_loc - [RX_loc(1); RX_loc(2); z_Rx_optqq]));
% end

    %%

% figure(), hold on,
% plot(y_near,Signal_Strength_RISRx_dBm_Corr_yVariations,'LineWidth',3);
% plot(y_near1,Signal_Strength_RISRx_dBm_Corr_yVariations1,'LineWidth',3);
% 
% plot(y_near,Signal_Strength_RISRx_dBm_Corr_yVariationsq,'LineWidth',3);
% plot(y_near1q,Signal_Strength_RISRx_dBm_Corr_yVariations1q,'LineWidth',3);
% 
% plot(y_near,Signal_Strength_RISRx_dBm_Corr_yVariationsqq,'LineWidth',3);
% plot(y_near1qq,Signal_Strength_RISRx_dBm_Corr_yVariations1qq,'LineWidth',3);
% 
% plot(RX_loc(2),Signal_Strength_RISRx_dBm_Corrqq,'x','LineWidth',3); 
% plot(linspace(y_near(1),y_near(end),1e5),Signal_Strength_RISRx_dBm_Corrqq*ones(1,1e5),'r--','LineWidth',3); 
%  hold off
% grid on, axis tight; xlabel('y direction [m]'), ylabel('Receiver Power [dBm]');
% title('Received power vs $y$ direction','Interpreter','latex');
% legend('Received Signal Strength with continuous phase','Interpolated Receiver Signal Strength with continuous phase', ...
%     'Received Signal Strength with 1 bit quantized phase','Interpolated Receiver Signal Strength with 1 bit quantized phase', ...
%     'Received Signal Strength with 2 bit quantized phase','Interpolated Receiver Signal Strength with 2 bit quantizeds phase', ...
%     'Power at Rx','Threshold');
% 

% figure(), hold on
% plot(y_near,Signal_Strength_RISRx_dBm_Corr_yVariations,'LineWidth',3);
% plot(y_near1,Signal_Strength_RISRx_dBm_Corr_yVariations1,'LineWidth',3);
% plot(RX_loc(2),Signal_Strength_RISRx_dBm_Corr,'x','LineWidth',3); 
% plot(linspace(y_near(1),y_near(end),1e5),Signal_Strength_RISRx_dBm_Corr*ones(1,1e5),'r--','LineWidth',3); 
%  hold off
% grid on, axis tight; xlabel('y direction [m]'), ylabel('Receiver Power [dBm]');
% title('Received power vs $y$ direction and continuous phase values','Interpreter','latex');
% legend('Received Signal Strength','Interpolated Receiver Signal Strength','Power at Rx','Threshold');
% 
% 
% 
% 
% figure(), hold on
% plot(z_near,Signal_Strength_RISRx_dBm_Corr_zVariations,'LineWidth',3);
% plot(z_near1,Signal_Strength_RISRx_dBm_Corr_zVariations1,'LineWidth',3);
% plot(RX_loc(3),Signal_Strength_RISRx_dBm_Corr,'x','LineWidth',3); 
% plot(linspace(z_near(1),z_near(end),1e5),Signal_Strength_RISRx_dBm_Corr*ones(1,1e5),'r--','LineWidth',3); 
%  hold off
% grid on, axis tight; xlabel('z direction [m]'), ylabel('Receiver Power [dBm]');
% title('Received power vs $z$ direction and continuous phase values','Interpreter','latex');
% legend('Received Signal Strength','Interpolated Receiver Signal Strength','Power at Rx','Threshold');
% 
% disp('Error of beam pointing');
% disp(['Error in azimuth angle   = ' , num2str(Phi_Error), ' °']);
% disp(['Error in elevation angle = ' , num2str(Theta_Error), ' °']);



% figure(), hold on
% plot(z_near,Signal_Strength_RISRx_dBm_Corr_zVariationsq,'LineWidth',3);
% plot(z_near1q,Signal_Strength_RISRx_dBm_Corr_zVariations1q,'LineWidth',3);
% plot(RX_loc(3),Signal_Strength_RISRx_dBm_Corrq,'x','LineWidth',3); 
% plot(linspace(z_near(1),z_near(end),1e5),Signal_Strength_RISRx_dBm_Corrq*ones(1,1e5),'r--','LineWidth',3); 
%  hold off
% grid on, axis tight; xlabel('z direction [m]'), ylabel('Receiver Power [dBm]');
% title('Received power vs $z$ direction and 1 bit quantzed phase values','Interpreter','latex');
% legend('Received Signal Strength','Interpolated Receiver Signal Strength','Power at Rx','Threshold');
% 
% 
% figure(), hold on
% plot(y_near,Signal_Strength_RISRx_dBm_Corr_yVariationsq,'LineWidth',3);
% plot(y_near1q,Signal_Strength_RISRx_dBm_Corr_yVariations1q,'LineWidth',3);
% plot(RX_loc(2),Signal_Strength_RISRx_dBm_Corrq,'x','LineWidth',3); 
% plot(linspace(y_near(1),y_near(end),1e5),Signal_Strength_RISRx_dBm_Corrq*ones(1,1e5),'r--','LineWidth',3); 
%  hold off
% grid on, axis tight; xlabel('y direction [m]'), ylabel('Receiver Power [dBm]');
% title('Received power vs $y$ direction and 1 bit quantzed phase values','Interpreter','latex');
% legend('Received Signal Strength','Interpolated Receiver Signal Strength','Power at Rx','Threshold');




% disp('Error of beam pointing');
% disp(['Error in azimuth angle   = ' , num2str(Phi_Error), ' °']);
% disp(['Error in elevation angle = ' , num2str(Theta_Error), ' °']);

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

%% Valutiamo il ricevitore in differenti posizioni, su una circonferenza con centro la RIS
% Bandwidth at -3 dB of RIS is 5.33°. 90°/(5.33/2) = 34
% NN = ceil(90/BW);
% th = linspace(pi/2,pi,NN);
% distance = 1.5;
% xunit =distance * sin(th) + RIS_loc(1);
% yunit = distance * cos(th) + RIS_loc(2);
% % figure(), plot(xunit, yunit); grid on, axis tight;
% % clear rx; clear p; clear BWRx; clear angRx;
% for p =1:NN
% %     d_ris_r(p)=norm(RIS_loc-[xunit(p);yunit(p);RIS_loc(3)]); %distance between receiver and RIS
% %     theta_rx(p) =  I_theta * asind( abs(RIS_loc(3)-RIS_loc(3))/d_ris_r(p) ); %elevation angle between normal of RIS and receiver
% %     i_phi=sign(yunit(p) - RIS_loc(2));
% %     phi_rx(p)=i_phi * atand( abs(yunit(p)-RIS_loc(2))/ abs(xunit(p)-RIS_loc(1)) );%azimuth angle between normal of RIS and receiver
% %     clear zn, clear ym;
% %     phi = zeros(N_ele,N_ele);
% %     for z = 1:N_ele
% %         for y = 1:N_ele
% % 
% %             zn = (z-1)*dy; ym = (y-1)*dy;
% %             phi(z,y) = k0*(ym*sind(90-theta_rx(p))*sind(-phi_rx(p))+zn*cosd(90-theta_rx(p)));
% % 
% %             while phi(z,y)>=2*pi
% %                 phi(z,y)=phi(z,y)-2*pi;
% %             end
% %             while phi(z,y)<=-2*pi
% %                 phi(z,y)=phi(z,y)+2*pi;
% %             end
% %             if phi(z,y) < 0
% %                 phi(z,y) = phi(z,y) + 2*pi;
% %             end
% %             if phi(z,y)>pi/2 && phi(z,y)<=3*pi/2
% %                 phiq(z,y) = pi;
% %             else
% %                 phiq(z,y) = 0;
% %             end
% % 
% %             if phi(z,y) > pi/4 && phi(z,y) <= 3/4*pi
% %                 phiqq(z,y) = pi/2;
% %             elseif phi(z,y) > 5/4*pi && phi(z,y)<=3*pi/2
% %                 phiqq(z,y) = 3*pi/2;
% %             elseif phi(z,y)>3/4*pi && phi(z,y) <=5/4*pi
% %                 phiqq(z,y) = pi;
% %             else
% %                 phiqq(z,y) = 0;
% %             end
% % 
% %         end
% %     end
% %     clear taperedRIStqq.Taper;
% % 
% %     taperedRISt.Taper = exp(1i*phi);
% %     taperedRIStq.Taper = exp(1i*phiq);
% %     taperedRIStqq.Taper = exp(1i*phiqq);
% % 
%      rx(p) = rxsite("cartesian", ...
%          "Antenna",antenna, ...
%          "AntennaPosition",[xunit(p);yunit(p);RIS_loc(3)], ...
%          "AntennaAngle",[0;90],"ReceiverSensitivity",-150);
% %     RIStx.Antenna = taperedRISt;
% %     RIStxq.Antenna = taperedRIStq;
% %     RIStxqq.Antenna = taperedRIStqq;
% % 
% %     Signal_Strength_RISRx_dBm_corr(p) = sigstrength(rx(p),RIStx, ...
% %         'PropagationModel',pm0, ...
% %         'Map',Stanza);
% %     Signal_Strength_RISRx_corr(p) = 10^((Signal_Strength_RISRx_dBm_corr(p)-30)/10)*10^3;
% % 
% %     Signal_Strength_RISRx_dBm_corrq(p) = sigstrength(rx(p),RIStxq, ...
% %         'PropagationModel',pm0, ...
% %         'Map',Stanza);
% %     Signal_Strength_RISRx_corrq(p) = 10^((Signal_Strength_RISRx_dBm_corrq(p)-30)/10)*10^3;
% % 
% %     Signal_Strength_RISRx_dBm_corrqq(p) = sigstrength(rx(p),RIStxqq, ...
% %         'PropagationModel',pm0, ...
% %         'Map',Stanza);
% %     Signal_Strength_RISRx_corrqq(p) = 10^((Signal_Strength_RISRx_dBm_corrqq(p)-30)/10)*10^3;
% % 
% % 
% % 
% %     clear BWrx, clear angrx;
% %     clear BWrxq, clear angrxq;
% %     clear BWrxqq, clear angrxqq;
% %     [BWrx, angrx] = beamwidth(taperedRISt,f,'dbDown',3); % -3 dB beamwidth 
% %     BWRx(p,:) = BWrx; angRx(p,:) = angrx; 
% % 
% %     [BWrxq, angrxq] = beamwidth(taperedRIStq,f,'dbDown',3); % -3 dB beamwidth 
% %     BWRxq(p,:) = BWrxq; angRxq(p,:) = angrxq; 
% % 
% %     [BWrxqq, angrxqq] = beamwidth(taperedRIStqq,f,'dbDown',3); % -3 dB beamwidth 
% %     BWRxqq(p,:) = BWrxqq; angRxqq(p,:) = angrxqq;
% %     if p == 1 || p == NN
% %         figure(), pattern(taperedRIStqq,f,-180:1:180,-90:1:90);
% %     end
% end
% % 
% %%
% figure(); plot(xunit, yunit,'LineWidth',3), 
% hold on ,yyaxis right,plot(xunit, Signal_Strength_RISRx_dBm_corr,'LineWidth',3),
% plot(yunit, Signal_Strength_RISRx_dBm_corr,'LineWidth',3), hold off;
% grid on, axis tight, legend('Movement path of receiver','Receiver power in $x$ direction','Receiver power in $y$ direction','interpreter','latex');
% xlabel('x direction [m]'); yyaxis left, ylabel('y direction [m]'); yyaxis right, ylabel('Received Signal Strength [dBm]');
% title('Received power along circular path');


% siteviewer("SceneModel",Stanza);
% show(tx,'ShowAntennaHeight', false);
% show(rx,'ShowAntennaHeight', true);
% show(RIStx,'ShowAntennaHeight', true);
% %show(RISrx,'ShowAntennaHeight', true);
% pattern(tx,f); pattern(RIStx,f); 
% siteviewer("SceneModel",Stanza);
% show(tx,'ShowAntennaHeight', false);
% show(rx,'ShowAntennaHeight', false);
% show(RIStx,'ShowAntennaHeight', false);
% show(RISrx,'ShowAntennaHeight', false); pattern(RIStxq,f);
% siteviewer("SceneModel",Stanza);
% show(tx,'ShowAntennaHeight', false);
% show(rx,'ShowAntennaHeight', false);
% show(RIStx,'ShowAntennaHeight', false);
% show(RISrx,'ShowAntennaHeight', false); pattern(RIStxqq,f);


%%
% figure(), plot3(xunit,yunit,ones(1,NN)*min(Signal_Strength_RISRx_dBm_corrq),'LineWidth',3), hold on;
% plot3(xunit, yunit,Signal_Strength_RISRx_dBm_corr,'LineWidth',3);
% plot3(xunit, yunit,Signal_Strength_RISRx_dBm_corrq,'LineWidth',3);
% plot3(xunit, yunit,Signal_Strength_RISRx_dBm_corrqq,'LineWidth',3);
% 
% hold off
% grid on, axis tight, xlabel('x direction [m]'), ylabel('y direction [m]'); zlabel('Received Signal Strength [dBm]');
% legend('Circular path','Received Signal Strength with continuous phase','Received Signal Strength with 1 bit quantized phase','Received Signal Strength with 2 bits quantized phase');
% title('Received Signale Strength along circular path in $xy$ plane','interpreter','latex');
%% |Tx-RIS Channel| 
TX_loc = TX_loc'; RX_loc = RX_loc'; RIS_loc = RIS_loc';
LOS = rays{1,1}.LineOfSight;

% i parametri sono calcolati in corrispondenza di f = 24.2 GHz
% n = Path Loss Exponent
% sigma = Shadow Fading Term (dB)
% b = Path Loss Parameter
%%
% NLOS Parameters
n_NLOS=3.19;          % Path Loss Exponent (Indoor Office NLOS)
sigma_NLOS=0;%abs(randn*8.29);      % Shadow Fading Term (dB) (Indoor Office NLOS)
b_NLOS=0.06;          % Path Loss Parameter (Indoor Office NLOS)
f0=24.2e9;            % Path Loss Parameter (GHz) (Indoor Office NLOS)

% LOS Parameters
n_LOS=1.73;           % Path Loss Exponent (Indoor Office LOS)
sigma_LOS=0;%abs(randn*3.02);       % Shadow Fading Term (dB) (Indoor Office LOS)
b_LOS=0;              % Path Loss Parameter (Indoor Office NLOS)

% Generalized Element Radiation Pattern 
q=0.285;
Gain=pi;

% Poisson distribution for cluster generation
lambda_p = 1.8; %Value referred to 28 GHz

%random phase variation
Phi_var = 1;%exp(1i*rand*2*pi);

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
    % array_RIS_LOS=zeros(1,N_ele^2);
    counter2 = 1;
    % counter3 = 1;

    % for z=0:N_ele-1
    %     for y=0:N_ele-1
    %         clear pp, clear pp_q, clear pp_qq;
    %         %pp = k0*((y-1)*dx*sind(90-Elevation(ceil(length(Elevation)/2)+theta_Tx_LOS))*sind(-Azimuth(ceil(length(Azimuth)/2)+phi_Tx_LOS))+(z-1)*dx*cosd(90-Elevation(ceil(length(Elevation)/2)+theta_Tx_LOS)));
    % 
    %         pp = k0*dx*(z*cosd(90-theta_T_RIS_LOS) + y*sind(-phi_T_RIS_LOS)*sind(90-theta_T_RIS_LOS) );
    %         while pp>=2*pi
    %             pp=pp-2*pi;
    %         end
    %         while pp<=-2*pi
    %             pp=pp+2*pi;
    %         end
    %         if pp < 0
    %             pp = pp + 2*pi;
    %         end
    % 
    %         if pp>pi/2 && pp<=3*pi/2
    %             pp_q = pi;
    %         else
    %             pp_q = 0;
    %         end
    % 
    %         if pp > pi/4 && pp <= 3/4*pi
    %             pp_qq = pi/2;
    %         elseif pp > 5/4*pi && pp<=3*pi/2
    %             pp_qq = 3*pi/2;
    %         elseif pp>3/4*pi && pp <=5/4*pi
    %             pp_qq = pi;
    %         else
    %             pp_qq = 0;
    %         end
    % 
    % 
    % 
    %         %array_RIS_LOS(counter3)=exp(1i*k0*dx*(y*sind(theta_T_RIS_LOS) + z*sind(phi_T_RIS_LOS)*cosd(theta_T_RIS_LOS) )) ;
    %         array_RIS_LOS(counter3)=exp(1i*pp) ;
    %         array_RIS_LOS_q(counter3)=exp(1i*pp_q) ;
    %         array_RIS_LOS_qq(counter3)=exp(1i*pp_qq) ;
    % 
    % 
    %         %array_RIS_LOS(counter3)=exp(1i*k0*dx*(z*cosd(90-theta_T_RIS_LOS) + y*sind(-phi_T_RIS_LOS)*sind(90-theta_T_RIS_LOS) )) ;
    %         counter3=counter3+1;
    %     end
    % end

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
   

    array_RIS_LOS = ones(1,N_ele^2);
    h_LOS_wp=sqrt(L_LOS_wp)*transpose(array_RIS_LOS)*array_Tx_LOS*Phi_var.*Ei(ceil(length(Azimuth)/2)+round(phi_T_RIS_LOS),ceil(length(Elevation)/2)+round(theta_T_RIS_LOS));
    %h_LOS_wp=sqrt(L_LOS_wp)*transpose(array_RIS_LOS)*array_Tx_LOS*sqrt(Gain*(cosd(theta_T_RIS_LOS))^(2*q));
else
    h_LOS_wp = 0;
end

%% STEP 2
% Generate Clusters/Sub-rays, Azimuth/Elevation Departure Angles and Cluster Distances

% for generate=1:1000  % To ensure that at least one scatterer exist
% 
%     % Number of Clusters
%     C=max([1,poissrnd(lambda_p)]);  % Poisson distributed
%     save("C.mat","C");
% 
%     % Number of Sub-rays per Cluster
%     S=randi(30,1,C); % Uniformly distributed
%     save("S.mat","S");
% 
%     % Azimuth/Elevation Departure Angles
%     phi_Tx=[ ];
%     theta_Tx=[ ];
%     phi_av=zeros(1,C);
%     theta_av=zeros(1,C);
%     for counter=1:C
%         phi_av(counter)  = rand*180-90;     % mean azimuth departure angles of sub-rays (Laplacian distributed)
%         theta_av(counter)= rand*90-45;      % mean elevation departure angles of sub-rays (Laplacian distributed)
% 
%         % cluster angles: First S(1) belongs to Cluster 1, Next S(2) belongs to Cluster 2....
%         phi_Tx   = [phi_Tx,log(rand(1,S(counter))./rand(1,S(counter)))*sqrt(25/2) + phi_av(counter)];
%         theta_Tx = [theta_Tx,log(rand(1,S(counter))./rand(1,S(counter)))*sqrt(25/2) + theta_av(counter)];
%     end
%     save("phi_av.mat","phi_av");
%     save("theta_av.mat","theta_av");
%     save("phi_Tx.mat","phi_Tx");
%     save("theta_Tx.mat","theta_Tx");
%     % Cluster Distances
%     a_c=1+rand(1,C)*(d_T_RIS-1);    % Cluster distances uniform [1,d_T_RIS] 
%     save("a_c.mat","a_c");
%     % Correction on Cluster Locations for Indoors
%     % Room dimensions (Indoor Hotspot)
%     %dim=[75,50,3.5];                  % x-y dimensions recommended by 5G Channel Model, height is assumed as 3.5 m
% 
%     Coordinates=zeros(C,3);          % for Clusters
%     Coordinates2=zeros(sum(S),3);    % for Scatterers
%     for counter=1:C
%         loop = 1;
%         Coordinates(counter,:)=[TX_loc(1) - a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
%             TX_loc(2) - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
%             TX_loc(3) + a_c(counter)*sind(theta_av(counter))] ;
%         while Coordinates(counter,3)>dim(3) || Coordinates(counter,3)<0 ||  Coordinates(counter,2)>dim(2) ||  Coordinates(counter,2)<0  ||  Coordinates(counter,1)>dim(1) ||  Coordinates(counter,1)<0
%             a_c(counter)=    0.95*a_c(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
%             % Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
%             % in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
%             Coordinates(counter,:)=[TX_loc(1) - a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
%                 TX_loc(2) - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
%                 TX_loc(3) + a_c(counter)*sind(theta_av(counter))] ;
%         end
%         while Coordinates(counter,3)>=z1_obs && Coordinates(counter,3)<=z2_obs &&  Coordinates(counter,2)>=y1_obs &&  Coordinates(counter,2)<=y2_obs  &&  Coordinates(counter,1)>=x1_obs &&  Coordinates(counter,1)<=x2_obs
%             a_c(counter)=    0.95*a_c(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
%             % Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
%             % in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
%             Coordinates(counter,:)=[TX_loc(1) - a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
%                 TX_loc(2) - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
%                 TX_loc(3) + a_c(counter)*sind(theta_av(counter))] ;
%             if RX_loc(2)<y1_obs
%                 while Coordinates(counter,1)>= TX_loc(1) || Coordinates(counter,2)< RIS_loc(2)
%                     a_c(counter)=    0.95*a_c(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
%                     % Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
%                     % in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
%                     Coordinates(counter,:)=[TX_loc(1) - a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
%                         TX_loc(2) - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
%                         Coordinates(counter,3)] ;
%                 end
%             else
%                 while Coordinates(counter,1)>= TX_loc(1) && Coordinates(counter,2)> RIS_loc(2)
%                     a_c(counter)=    0.95*a_c(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
%                     % Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
%                     % in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
%                     Coordinates(counter,:)=[TX_loc(1) - a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
%                         TX_loc(2) - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
%                         Coordinates(counter,3)] ;
%                 end
%             end
% 
%         end
%         if RX_loc(2)<y1_obs
%                 while Coordinates(counter,1)>= TX_loc(1) || Coordinates(counter,2)< RIS_loc(2)
%                     a_c(counter)=    0.95*a_c(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
%                     % Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
%                     % in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
%                     Coordinates(counter,:)=[TX_loc(1) - a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
%                         TX_loc(2) - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
%                         Coordinates(counter,3)] ;
%                 end
%             else
%                 while Coordinates(counter,1)>= TX_loc(1) || Coordinates(counter,2)> RIS_loc(2)
%                     a_c(counter)=    0.95*a_c(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
%                     % Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
%                     % in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
%                     Coordinates(counter,:)=[TX_loc(1) - a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
%                        TX_loc(2) - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
%                         Coordinates(counter,3)] ;
%                 end
%         end
%     end
% 
%     a_c_rep=[];
%     for counter3=1:C
%         a_c_rep=[a_c_rep,repmat(a_c(counter3),1,S(counter3))];
%     end
%     for counter2=1:sum(S)
%         Coordinates2(counter2,:)=[TX_loc(1) - a_c_rep(counter2)*cosd(theta_Tx(counter2))*cosd(phi_Tx(counter2)),...
%             TX_loc(2) - a_c_rep(counter2)*cosd(theta_Tx(counter2))*sind(phi_Tx(counter2)),...
%             TX_loc(3) + a_c_rep(counter2)*sind(theta_Tx(counter2))] ;
%     end
%     save("a_c_rep.mat","a_c_rep");
%     % Correction on Scatters
%     % You may ignore the scatterers outside the walls for Indoors scenario
%     ignore=[];
% 
%     for counter2=1:sum(S)
%         if Coordinates2(counter,3)>dim(3) || Coordinates2(counter,3)<0 ||  Coordinates2(counter,2)>dim(2) ||  Coordinates2(counter,2)<0  ||  Coordinates2(counter,1)>dim(1) ||  Coordinates2(counter,1)<0
%             ignore=[ignore,counter2];   % contains the indices of scatterers that can be ignored
%         end
%         if Coordinates2(counter,3)>=z1_obs && Coordinates2(counter,3)<=z2_obs &&  Coordinates2(counter,2)>=y1_obs &&  Coordinates2(counter,2)<=y2_obs  &&  Coordinates2(counter,1)>=x1_obs &&  Coordinates2(counter,1)<=x2_obs
%             ignore=[ignore,counter2];   % contains the indices of scatterers that can be ignored
%         end
%         if Coordinates2(counter,1)>=TX_loc(1)
%             ignore = [ignore, counter2];
%         end
%     end
% 
%     % updated indices
%     indices=setdiff(1:sum(S),ignore);    % the set of active scatterer indices
%     M_new=length(indices);               % number of IOs inside the room or above ground
% 
%     % Neccessary Loop to have at least one scatterer
%     if M_new>0 % if M_new==0 --> all scatters are outside
%         break  % break generate=1:100 % if M_new >0 we are OK (at least one scatter)
%     end
% end  
% save("M_new.mat","M_new");
% save("indices.mat","indices");
% save("Coordinates.mat","Coordinates");
% save("Coordinates2.mat","Coordinates2");
% save("ignore.mat","ignore");
%%
C = load("C.mat"); C = C.C;
S = load("S.mat"); S = S.S;
phi_Tx = load("phi_Tx.mat"); phi_Tx = phi_Tx.phi_Tx;
theta_Tx = load("theta_Tx.mat"); theta_Tx = theta_Tx.theta_Tx;
a_c = load("a_c.mat"); a_c = a_c.a_c;
M_new = load("M_new.mat"); M_new = M_new.M_new;
indices = load("indices.mat"); indices = indices.indices;
Coordinates = load("Coordinates.mat"); Coordinates = Coordinates.Coordinates;
Coordinates2 = load("Coordinates2"); Coordinates2 = Coordinates2.Coordinates2;
ignore = load("ignore.mat"); ignore = ignore.ignore;
a_c_rep = load("a_c_rep.mat"); a_c_rep = a_c_rep.a_c_rep;
phi_av = load("phi_av.mat"); phi_av = phi_av.phi_av;
theta_av = load("theta_av.mat"); theta_av = theta_av.theta_av;
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
%array_cs_RIS=zeros(sum(S),N_ele^2);

% for counter2=indices
%     counter3=1;
%     for z=0:N_ele-1
%         for y=0:N_ele-1
%             %array_cs_RIS(counter2,counter3)=exp(1i*k0*dx*(x*sind(theta_cs_RIS(counter2)) + y*sind(phi_cs_RIS(counter2))*cosd(theta_cs_RIS(counter2)) )) ;
%             clear pp, clear pp_q, clear pp_qq;
%             pp = k0*dx*(z*cosd(90-theta_cs_RIS(counter2)) + y*sind(-phi_cs_RIS(counter2))*sind(90-theta_cs_RIS(counter2)) );
%             %%pp = k0*((y-1)*dx*sind(90-Elevation(ceil(length(Elevation)/2)+theta_cs_RIS(counter2)))*sind(-Azimuth(ceil(length(Azimuth)/2)+phi_cs_RIS))+(z-1)*dx*cosd(90-Elevation(ceil(length(Elevation)/2)+theta_cs_RIS)));
% 
%             while pp>=2*pi
%                 pp=pp-2*pi;
%             end
%             while pp<=-2*pi
%                 pp=pp+2*pi;
%             end
%             if pp < 0
%                 pp = pp + 2*pi;
%             end
% 
%             if pp>pi/2 && pp<=3*pi/2
%                 pp_q = pi;
%             else
%                 pp_q = 0;
%             end
% 
%             if pp > pi/4 && pp <= 3/4*pi
%                 pp_qq = pi/2;
%             elseif pp > 5/4*pi && pp<=3*pi/2
%                 pp_qq = 3*pi/2;
%             elseif pp>3/4*pi && pp <=5/4*pi
%                 pp_qq = pi;
%             else
%                 pp_qq = 0;
%             end
% 
% 
%             array_cs_RIS(counter2,counter3)=exp(-1i*pp) ;
%             array_cs_RIS_q(counter2,counter3)=exp(1i*pp_q) ;
%             array_cs_RIS_qq(counter2,counter3)=exp(1i*pp_qq) ;
% 
%             %array_cs_RIS(counter2,counter3)=exp(1i*k0*dx*(y*cosd(90-theta_cs_RIS(counter2)) + x*sind(-phi_cs_RIS(counter2))*sind(90-theta_cs_RIS(counter2)) )) ;
%             counter3=counter3+1;
%         end
%     end
% end

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
array_cs_RIS = ones(length(indices),N_ele^2);
for counter2=indices
    X_sigma=sigma_NLOS;

    Lcs_dB_wp=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs(counter2))- X_sigma;
    Lcs_wp=10^(Lcs_dB_wp/10);

    beta(counter2)=((randn+1i*randn)./sqrt(2));  % common complex gain for shared clusters
    shadow(counter2)=X_sigma;                    % commun shadow factor for shared clusters
    h_NLOS_wp = h_NLOS_wp + Ei(ceil(length(Azimuth)/2)+round(phi_cs_RIS(counter2)),ceil(length(Elevation)/2)+round(theta_cs_RIS(counter2))).*beta(counter2)*sqrt(Lcs_wp)*transpose(array_cs_RIS(counter2,:))*array_Tx_cs(counter2,:);  % consider all scatters
end

h_NLOS_wp=h_NLOS_wp.*sqrt(1/M_new);  % normalization

h_Cluster=(h_NLOS_wp+h_LOS_wp); % Tx-RIS with cluster
h_wCluster = h_LOS_wp;
%% STEPS 6-7 
% Generation of g (RIS-Rx Channel) - GENERATE A LOS CHANNEL

I_theta=sign(RX_loc(3) - RIS_loc(3));
theta_Rx_RIS=I_theta * asind( abs(RX_loc(3)-RIS_loc(3))/d_RIS_R ); % AoD of RIS
%theta_Rx_RIS=I_theta * asind( abs(RIS_loc(3)-RX_loc(3))/d_RIS_R );

% Azimuth Departure Angle
I_phi=sign(RX_loc(2) - RIS_loc(2));
phi_Rx_RIS=I_phi * atand( abs(RX_loc(2)-RIS_loc(2))/ abs(RX_loc(1)-RIS_loc(1)) );

% AoA angles of Rx for g_LOS channel in an Indoor
phi_av_Rx  = rand*180-90;     % mean azimuth
theta_av_Rx= rand*180-90;      % mean elevation

phi_Rx   = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + phi_av_Rx];
theta_Rx = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + theta_av_Rx];

% Recalculate Array Response for Two Angles (in Rx direction)
% array_2=zeros(1,N_ele^2);
% 
% counter3=1;
% for x=0:N_ele-1
%     for y=0:N_ele-1
%         %array_2(counter3)=exp(1i*k0*dx*(x*sind(theta_Rx_RIS) + y*sind(phi_Rx_RIS)*cosd(theta_Rx_RIS) )) ;
% 
%         % array_2(counter3)=exp(1i*k0*dx*(x*cosd(90-theta_Rx_RIS) + y*sind(-phi_Rx_RIS)*sind(90-theta_Rx_RIS) )) ;
%         clear pp, clear pp_q, clear pp_qq;
%         pp = k0*dx*(x*cosd(90-theta_Rx_RIS) + y*sind(-phi_Rx_RIS)*sind(90-theta_Rx_RIS) );
%         while pp>=2*pi
%             pp=pp-2*pi;
%         end
%         while pp<=-2*pi
%             pp=pp+2*pi;
%         end
%         if pp < 0
%             pp = pp + 2*pi;
%         end
% 
%         if pp>pi/2 && pp<=3*pi/2
%             pp_q = pi;
%         else
%             pp_q = 0;
%         end
% 
%         if pp > pi/4 && pp <= 3/4*pi
%             pp_qq = pi/2;
%         elseif pp > 5/4*pi && pp<=3*pi/2
%             pp_qq = 3*pi/2;
%         elseif pp>3/4*pi && pp <=5/4*pi
%             pp_qq = pi;
%         else
%             pp_qq = 0;
%         end
% 
% 
% 
%         array_2(x+1,y+1)=exp(1i*pp) ;
% 
%         %array_2(counter3)=exp(1i*k0*dx*(y*cosd(90-theta_Rx_RIS) + x*sind(-phi_Rx_RIS)*sind(90-theta_Rx_RIS) )) ;
%         % for p = 1:N_ele
%         %     array_2rx(p,counter3)=exp(1i*k0*dx*(x*cosd(90-theta_rx_RIS(p)) + y*sind(-phi_rx_RIS(p))*sind(90-theta_rx_RIS(p)) )) ;
%         %     %array_2rx(p,counter3)=exp(1i*k0*dx*(y*cosd(90-theta_rx_RIS(p)) + x*sind(-phi_rx_RIS(p))*sind(90-theta_rx_RIS(p)) )) ;
%         % end
%         counter3=counter3+1;
%     end
% end

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
array_2 = ones(1,N_ele^2);
g_wp = Ei1(ceil(length(Azimuth)/2)+round(phi_RX),ceil(length(Elevation)/2)+round(theta_RX)).*sqrt(L_LOS_2_wp)*transpose(array_2)*array_Rx*Phi_var;
g_wp_q = Ei1q(ceil(length(Azimuth)/2)+round(phi_RX),ceil(length(Elevation)/2)+round(theta_RX)).*sqrt(L_LOS_2_wp)*transpose(array_2)*array_Rx*Phi_var;
g_wp_qq = Ei1qq(ceil(length(Azimuth)/2)+round(phi_RX),ceil(length(Elevation)/2)+round(theta_RX)).*sqrt(L_LOS_2_wp)*transpose(array_2)*array_Rx*Phi_var;

% Calculate Departure Angles Considering RIS and Rx Coordinates
%d_RIS_R=norm(RIS_loc-Rx_loc);

% Elevation Departure Angle

% Phase_element20 = zeros(N_ele,N_ele); %initialization of continuous phase matrix values
% Phase_elementq20 = zeros(N_ele,N_ele); %initialization of 1 bit phase matrix values
% Phase_elementqq20 = zeros(N_ele,N_ele); %initialization of2 bits phase matrix values
% Ei20 = zeros(length(Azimuth),length(Elevation)); %Array factor without phase
% Ei120 = zeros(length(Azimuth),length(Elevation)); %Array factor for RX location with continuous phase values
% Ei1q20 = zeros(length(Azimuth),length(Elevation)); %Array factor for RX location with 1 bit phase values
% Ei1qq20 = zeros(length(Azimuth),length(Elevation)); %Array factor for RX location with 2 bits phase values
% pp = 0;
% 
% for p=1:NN
%     i_theta = sign(RIS_loc(3) - RIS_loc(3));
%     theta_rx_RIS(p)=i_theta * asind( abs(RIS_loc(3)-RIS_loc(3))/distance );
%     i_phi=sign(yunit(p) - RIS_loc(2));
%     phi_rx_RIS(p) = i_phi * atand( abs(yunit(p)-RIS_loc(2))/ abs(xunit(p)-RIS_loc(1)) );
%     for az = 1:length(Azimuth)
%         for el = 1:length(Elevation)
%             for z = 1:N_ele
%                 for y = 1:N_ele
%                     clear pp;
%                     %Phi = k0*((x-1)*dx*sind(Elevation(el))*sind(Azimuth(az))+(y-1)*dx*sind(Elevation(el))*cosd(Azimuth(az)));
%                     pp = k0*((z-1)*dx*sind(90-Elevation(el))*sind(-Azimuth(az))+(y-1)*dx*cosd(90-Elevation(el)));
%                     while pp>=2*pi
%                         pp=pp-2*pi;
%                     end
%                     while pp<=-2*pi
%                         pp=pp+2*pi;
%                     end
%                     while pp <0
%                         pp = pp +2*pi;
%                     end
% 
%                     zn = (z-1)*dx; ym = (y-1)*dx;
%                     Phase_element20(z,y) = k0*(zn*sind(90-theta_rx_RIS(p))*sind(-phi_rx_RIS(p))+ym*cosd(90-theta_rx_RIS(p)));
%                     while Phase_element20(z,y)>=2*pi
%                         Phase_element20(z,y)=Phase_element20(z,y)-2*pi;
%                     end
%                     while Phase_element20(z,y)<=-2*pi
%                         Phase_element20(z,y)=Phase_element20(z,y)+2*pi;
%                     end
%                     if Phase_element20(z,y) < 0
%                         Phase_element20(z,y) = Phase_element20(z,y) + 2*pi;
%                     end
% 
%                     if Phase_element20(z,y)>pi/2 && Phase_element20(z,y)<3*pi/2
%                         Phase_elementq20(z,y) = pi;
%                     else
%                         Phase_elementq20(z,y) = 0;
%                     end
% 
%                     if Phase_element20(z,y) >= pi/4 && Phase_element20(z,y) <= 3/4*pi
%                         Phase_elementqq20(z,y) = pi/2;
%                     elseif Phase_element20(z,y) > 5/4*pi && Phase_element20(z,y)<=3*pi/2
%                         Phase_elementqq20(z,y) = 3*pi/2;
%                     elseif Phase_element20(z,y)>3/4*pi && Phase_element20(z,y) <=5/4*pi
%                         Phase_elementqq20(z,y) = pi;
%                     else
%                         Phase_elementqq20(z,y) = 0;
%                     end
% 
% 
% 
%                     Ei20(az,el) = Ei20(az,el) + (exp(-1i*pp));
%                     Ei120(az,el) = Ei120(az,el) + (exp(-1i*pp).*exp(1i*Phase_element20(z,y)));
%                     Ei1q20(az,el) = Ei1q20(az,el) + (exp(-1i*pp).*exp(1i*Phase_elementq20(z,y)));
%                     Ei1qq20(az,el) = Ei1qq20(az,el) + (exp(-1i*pp).*exp(1i*Phase_elementqq20(z,y)));
%                 end
%             end
%         end
%     end
%     Ei120 = Ei120./max(max(Ei20));
%     Ei1q20 = Ei1q20./max(max(Ei20));
%     Ei1qq20 = Ei1qq20./max(max(Ei20));
% 
%     L_dB_LOS_2_wp20(p)=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(distance)- sigma_LOS ;
%     L_LOS_2_wp20(p)=10^(L_dB_LOS_2_wp20(p)/10);
% 
% 
% 
%     % Generate g (Pure LOS)
%     array_2 = ones(1,N_ele^2);
%     g_wp20(:,p) = Ei120(ceil(length(Azimuth)/2)+round(phi_rx_RIS(p)),ceil(length(Elevation)/2)+round(theta_rx_RIS(p))).*sqrt(L_LOS_2_wp20(p)).*transpose(array_2)*array_Rx*Phi_var;
%     g_wp_q20(:,p) = Ei1q20(ceil(length(Azimuth)/2)+round(phi_rx_RIS(p)),ceil(length(Elevation)/2)+round(theta_rx_RIS(p))).*sqrt(L_LOS_2_wp20(p)).*transpose(array_2)*array_Rx*Phi_var;
%     g_wp_qq20(:,p) = Ei1qq20(ceil(length(Azimuth)/2)+round(phi_rx_RIS(p)),ceil(length(Elevation)/2)+round(theta_rx_RIS(p))).*sqrt(L_LOS_2_wp20(p)).*transpose(array_2)*array_Rx*Phi_var;
% 
%     Phase_element20 = zeros(N_ele,N_ele); %initialization of continuous phase matrix values
%     Phase_elementq20 = zeros(N_ele,N_ele); %initialization of 1 bit phase matrix values
%     Phase_elementqq20 = zeros(N_ele,N_ele); %initialization of2 bits phase matrix values
%     Ei20 = zeros(length(Azimuth),length(Elevation)); %Array factor without phase
%     Ei120 = zeros(length(Azimuth),length(Elevation)); %Array factor for RX location with continuous phase values
%     Ei1q20 = zeros(length(Azimuth),length(Elevation)); %Array factor for RX location with 1 bit phase values
%     Ei1qq20 = zeros(length(Azimuth),length(Elevation)); %Array factor for RX location with 2 bits phase values
% 
% end



%% STEP 8
% Generation of h_SISO
RX_loc(1) = -RX_loc(1);
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

    Lcs_dB_SISO_wp=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs_tilde(counter2))- shadow(counter2) - Attenuation_Wall;
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

    L_SISO_LOS_dB_wp=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_R)- sigma_LOS - Attenuation_Wall;
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
   
    % h_SISO_LOS_wp= sqrt(L_SISO_LOS_wp)*exp(1i*rand*2*pi)*transpose(array_Rx_SISO)*array_Tx_SISO; 
    h_SISO_LOS_wp= sqrt(L_SISO_LOS_wp)*transpose(array_Rx_SISO)*array_Tx_SISO; 
end
h_SISO_wp=h_SISO_NLOS_wp + h_SISO_LOS_wp; 

% %%
% Phi = Phi';
% Phiq = Phiq';
% Phiqq = Phiqq';
% Phii = Phii';
% Phiiq = Phiiq';
% Phiiqq = Phiiqq';
% %% total channel
% 
% %   g     = Vector of LOS channel coefficients between the RIS and the Rx
% %  THETA  = Matrix of RIS element responses
% %   h     = Vector of channel coefficients for the Tx-RIS link composed of M scatters
% % h_SISO  = Characterizes the direct link channel between Tx and Rx
% %   x     = Transmitted signal
% 
% PHI = Phi(:); %PHI = PHI(:)';
% PHIq = Phiq(:); %PHIq = PHIq(:)';
% PHIqq = Phiqq(:); %PHIqq = PHIqq(:)';
% PHII = Phii(:); %PHI = PHI(:)';
% PHIIq = Phiiq(:); %PHIq = PHIq(:)';
% PHIIqq = Phiiqq(:); %PHIqq = PHIqq(:)';
% % 
% THETA = zeros(N_ele^2,N_ele^2);
% THETAq = zeros(N_ele^2,N_ele^2);
% THETAqq = zeros(N_ele^2,N_ele^2);
% THETAA = zeros(N_ele^2,N_ele^2);
% THETAAq = zeros(N_ele^2,N_ele^2);
% THETAAqq = zeros(N_ele^2,N_ele^2);
% %THETAA = zeros(N_ele^2,N_ele^2);
% % 
% for k = 1:N_ele^2
%     THETA(k,k) = exp(1i*PHI(k));
%     THETAq(k,k) = exp(1i*PHIq(k));
%     THETAqq(k,k) = exp(1i*PHIqq(k));
% end
% 
% %canale considerando il valore del guadagno, ma senza potenza in
% %trasmissione
% Tot_Channel_wCluster = transpose(g_wp) * THETA * h_wCluster; % Tx-RIS-Rx without Cluster
% Tot_Channel_wClusterq = transpose(g_wp) * THETAq * h_wCluster; % Tx-RIS-Rx without Cluster and 1 bit quantized phase
% Tot_Channel_wClusterqq = transpose(g_wp) * THETAqq * h_wCluster; % Tx-RIS-Rx without Cluster and 2 bit quantized phase
%%


Tot_Channel_wCluster = sqrt(Gain)*(transpose(g_wp)*h_wCluster);
Tot_Channel_wCluster_q = sqrt(Gain)*(transpose(g_wp_q)*h_wCluster);
Tot_Channel_wCluster_qq = sqrt(Gain)*(transpose(g_wp_qq)*h_wCluster);
Tt_Channel_wCluster = 20*log10(Tot_Channel_wCluster); Tt_Channel_wCluster = abs(Tt_Channel_wCluster);
Tt_Channel_wCluster_q = 20*log10(Tot_Channel_wCluster_q); Tt_Channel_wCluster_q = abs(Tt_Channel_wCluster_q);
Tt_Channel_wCluster_qq = 20*log10(Tot_Channel_wCluster_qq); Tt_Channel_wCluster_qq = abs(Tt_Channel_wCluster_qq);
Tt_Channel_imagwCluster = abs(20*log10(h_SISO_LOS_wp));
Tt_Channel_imagCluster = abs(20*log10(h_SISO_wp));


Tot_Channel_Cluster = sqrt(Gain)*transpose(g_wp)*h_Cluster;
Tot_Channel_Cluster_q = sqrt(Gain)*transpose(g_wp_q)*h_Cluster;
Tot_Channel_Cluster_qq = sqrt(Gain)*transpose(g_wp_qq)*h_Cluster;
Tt_Channel_Cluster = 20*log10(Tot_Channel_Cluster); Tt_Channel_Cluster = abs(Tt_Channel_Cluster);
Tt_Channel_Cluster_q = 20*log10(Tot_Channel_Cluster_q); Tt_Channel_Cluster_q = abs(Tt_Channel_Cluster_q);
Tt_Channel_Cluster_qq = 20*log10(Tot_Channel_Cluster_qq); Tt_Channel_Cluster_qq = abs(Tt_Channel_Cluster_qq);

tx_power_dBm_real =  linspace(0,30,7);
tx_power_Corr_dBm = tx_power_dBm_real + mag2db(d_T_RIS*d_RIS_R/(d_T_RIS+d_RIS_R))+mag2db(4*pi/lambda);

PrxwC = tx_power_Corr_dBm - Tt_Channel_wCluster;
PrxwC_q = tx_power_Corr_dBm - Tt_Channel_wCluster_q;
PrxwC_qq = tx_power_Corr_dBm - Tt_Channel_wCluster_qq;
Prx = tx_power_Corr_dBm - Tt_Channel_Cluster;
Prx_q = tx_power_Corr_dBm - Tt_Channel_Cluster_q;
Prx_qq = tx_power_Corr_dBm - Tt_Channel_Cluster_qq;


Prx_imagwCluster = tx_power_dBm_real - Tt_Channel_imagwCluster;
Prx_imagCluster = tx_power_dBm_real - Tt_Channel_imagCluster;


disp('Received power ');
disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster = ', num2str(Prx), ' dBm']);
disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster 1 BIT Quantization = ', num2str(Prx_q), ' dBm']);
disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster 2 BIT Quantization = ', num2str(Prx_qq), ' dBm']);
disp(['-> Tx-RIS-Rx Received Signal Strength without clusters = ', num2str(PrxwC), ' dBm']);
disp(['-> Tx-RIS-Rx Received Signal Strength without clusters 1 BIT Quantization = ', num2str(PrxwC_q), ' dBm']);
disp(['-> Tx-RIS-Rx Received Signal Strength without clusters 2 BIT Quantization = ', num2str(PrxwC_qq), ' dBm']);
disp(['-> Tx-Rx Received Signal Strength with Cluster = ', num2str(Prx_imagCluster), ' dBm']);
disp(['-> Tx-Rx Received Signal Strength without Cluster = ', num2str(Prx_imagwCluster), ' dBm']);
RSSI = [Prx;Prx_q; Prx_qq; PrxwC; PrxwC_q; PrxwC_qq; Prx_imagCluster; Prx_imagwCluster];
save("RSSI10.mat","RSSI");
%%
Pn = -100;
tx_power = 10.^((tx_power_dBm_real-30)/10);
tx_power_Corr = 10.^((tx_power_Corr_dBm-30)/10);
Tx_noise = 10^((Pn-30)/10);
SNR_C = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_Cluster)^2))-Pn;
SNR_Cq = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_Cluster_q)^2))-Pn;
SNR_Cqq= tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_Cluster_qq)^2))-Pn;

SNR_wC = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_wCluster)^2))-Pn;
SNR_wCq = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_wCluster_q)^2))-Pn;
SNR_wCqq =tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_wCluster_qq)^2))-Pn;

SNR_imagwC =tx_power_dBm_real+ (10*log10(abs(h_SISO_LOS_wp)^2))-Pn;
SNR_imagC = tx_power_dBm_real+(10*log10(abs(h_SISO_wp)^2))-Pn;
disp('Signal Noise Ratio with Pn = -100 dBm ');
disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster = ', num2str(SNR_C), ' dB']);
disp(['-> Tx-RIS-Rx Received Signal Strength with Clusterand 1 BIT Quantization = ', num2str(SNR_Cq), ' dB']);
disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster 2 BIT Quantization = ', num2str(SNR_Cqq), ' dB']);

disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster = ', num2str(SNR_wC), ' dB']);
disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster 1 BIT Quantization = ', num2str(SNR_wCq), ' dB']);
disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster 2 BIT Quantization = ', num2str(SNR_wCqq), ' dB']);

disp(['-> Tx-Rx Received Signal Strength with Cluster = ', num2str(SNR_imagC), ' dB']);
disp(['-> Tx-Rx Received Signal Strength without Cluster = ', num2str(SNR_imagwC), ' dB']);

R_imagwC = (((log2(1+SNR_imagwC))));
R_imagC = (((log2(1+SNR_imagC))));
R_C = (((log2(1+SNR_C))));
R_Cq = (((log2(1+SNR_Cq))));
R_Cqq = (((log2(1+SNR_Cqq))));
R_wC = (((log2(1+SNR_wC))));
R_wCq = (((log2(1+SNR_wCq))));
R_wCqq = ((log2(1+SNR_wCqq)));
% for xx = 1:7
%  R_imagwC(xx) = floor((mean(log2(1+(tx_power_dBm_real(xx)+ (10*log10(abs(h_SISO_LOS_wp)^2))-Pn)))));
%  R_imagC(xx) = floor((mean(log2(1+(tx_power_dBm_real(xx)+ (10*log10(abs(h_SISO_wp)^2))-Pn)))));
%  R_C(xx) = floor(abs((mean(log2(1+(tx_power_Corr_dBm(xx)+(20*log10(abs(mean(transpose(g_wp)) *  mean(h_Cluster))^2))-Pn)))))); 
%  R_Cq(xx) = floor(abs((mean(log2(1+(tx_power_Corr_dBm(xx)+(10*log10(abs(mean(transpose(g_wp_q)) * mean(h_Cluster))^2))-Pn))))));
%  R_Cqq(xx) = floor(abs((mean(log2(1+(tx_power_Corr_dBm(xx)+(10*log10(abs(mean(transpose(g_wp_qq)) * mean(h_Cluster))^2))-Pn))))));
%  R_wC(xx) = floor(abs((mean(log2(1+(tx_power_Corr_dBm(xx)+(10*log10(abs(mean(transpose(g_wp)) * mean(h_wCluster))^2))-Pn))))));
%  R_wCq(xx) = floor(abs((mean(log2(1+(tx_power_Corr_dBm(xx)+(10*log10(abs(mean(transpose(g_wp_q)) * mean(h_wCluster))^2))-Pn))))));
%  R_wCqq(xx) = floor(abs((mean(log2(1+(tx_power_Corr_dBm(xx)+(10*log10(abs(mean(transpose(g_wp_qq)) * mean(h_wCluster))^2))-Pn))))));
% end
disp('Ergodic Achievable rate');
disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster = ', num2str(R_C), ' bits/s/Hz']);
disp(['-> Tx-RIS-Rx Received Signal Strength with Clusterand 1 BIT Quantization = ', num2str(R_Cq), ' bits/s/Hz']);
disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster 2 BIT Quantization = ', num2str(R_Cqq), ' bits/s/Hz']);

disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster = ', num2str(R_wC), ' bits/s/Hz']);
disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster 1 BIT Quantization = ', num2str(R_wCq), ' bits/s/Hz']);
disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster 2 BIT Quantization = ', num2str(R_wCqq), ' bits/s/Hz']);

disp(['-> Tx-Rx Received Signal Strength with Cluster = ', num2str(R_imagC), ' bits/s/Hz']);
disp(['-> Tx-Rx Received Signal Strength without Cluster = ', num2str(R_imagwC), ' bits/s/Hz']);
EnhancedBitRate = [R_imagwC;R_wC; R_wCq; R_wCqq; R_imagC; R_C; R_Cq; R_Cqq];
save("EnhanceBIT10.mat","EnhancedBitRate");
figure();
plot(tx_power_dBm_real,R_imagwC,'LineWidth',3); hold on;
plot(tx_power_dBm_real,R_wC,'LineWidth',3); hold off;
grid on, axis tight, title('Enhanced Bit Rate vs Power');
xlabel('Transmitted power [dBm]'), ylabel('Enhanced Bit Rate [bits/s/Hz]');
legend('Without RIS','With RIS');
%%
EBR400 = load("EnhanceBIT20.mat"); EBR400 = EBR400.EnhancedBitRate;
EBR196 = load("EnhanceBIT14.mat"); EBR196 = EBR196.EnhancedBitRate;
EBR100 = load("EnhanceBIT10.mat"); EBR100 = EBR100.EnhancedBitRate;
tx_power_dBm_real = linspace(0,30,7);
%RSSI = [Prx;Prx_q; Prx_qq; PrxwC; PrxwC_q; PrxwC_qq; Prx_imagCluster; Prx_imagwCluster];

RSSI400 = load("RSSI20.mat"); RSSI400 = RSSI400.RSSI;
RSSI196 = load("RSSI14.mat"); RSSI196 = RSSI196.RSSI;
RSSI100 = load("RSSI10.mat"); RSSI100 = RSSI100.RSSI;
%RSSI Tx-RIS-Rx for 400/196/100 antenna element, without scatterers and
%continuous phase values
figure();
plot(tx_power_dBm_real,RSSI400(8,:),'LineWidth',3); 
hold on;
plot(tx_power_dBm_real,RSSI400(4,:),'LineWidth',3); 
plot(tx_power_dBm_real,RSSI196(4,:),'LineWidth',3);
plot(tx_power_dBm_real,RSSI100(4,:),'LineWidth',3);
grid on, axis tight;
xlabel('Transmitted Power [dBm]');
ylabel('Received Power [dBm]');
title('RSSI vs Transmitted power without scatterers and continuous phase values');
legend('Without RIS','400 unit elements','196 unit elements','100 unit elements');
hold off;
%RSSI Tx-RIS-Rx for 400/196/100 antenna element, with scatterers and
%continuous phase values
figure();
plot(tx_power_dBm_real,RSSI400(8,:),'LineWidth',3); 
hold on;
plot(tx_power_dBm_real,RSSI400(1,:),'LineWidth',3); 
plot(tx_power_dBm_real,RSSI196(1,:),'LineWidth',3);
plot(tx_power_dBm_real,RSSI100(1,:),'LineWidth',3);
grid on, axis tight;
xlabel('Transmitted Power [dBm]');
ylabel('Received Power [dBm]');
title('RSSI vs Transmitted power with scatterers and continuous phase values');
legend('Without RIS','400 unit elements','196 unit elements','100 unit elements');
hold off;

%ABR Tx-RIS-Rx for 400/196/100 antenna element, without scatterers and
%continuous phase values
figure();
plot(tx_power_dBm_real,EBR400(1,:),'LineWidth',3); 
hold on;
plot(tx_power_dBm_real,EBR400(2,:),'LineWidth',3); 
plot(tx_power_dBm_real,EBR400(5,:),'LineWidth',3);
plot(tx_power_dBm_real,EBR400(6,:),'LineWidth',3);
grid on, axis tight;
xlabel('Transmitted Power [dBm]');
ylabel('Achievable Bit Rate [bits/s/Hz]');
title('Achievable Bit Rate vs Transmitted power with/without scatterers');
legend('Without RIS and without scatterers','With RIS and without scatterers','Without RIS and with scatterers','With RIS and with scatterers');
hold off;
% ABR Tx-RIS-Rx for 400 antenna element, without scatterers, continuous, 1
% bit and 2 bits quantized phase values
%continuous phase values
figure();
plot(tx_power_dBm_real,EBR400(1,:),'LineWidth',3); 
hold on;
plot(tx_power_dBm_real,EBR400(6,:),'LineWidth',3); 
plot(tx_power_dBm_real,EBR196(6,:),'LineWidth',3);
plot(tx_power_dBm_real,EBR100(6,:),'LineWidth',3);
grid on, axis tight;
xlabel('Transmitted Power [dBm]');
ylabel('Achievable Bit Rate [bits/s/Hz]');
title('Achievable Bit Rate vs Transmitted power with scatterers and continuous phase values');
legend('Without RIS','400 unit elements','196 unit elements','100 unit elements');
hold off;
%%

figure();
plot(tx_power_dBm_real,EBR400(5,:),'LineWidth',3); 
hold on;
plot(tx_power_dBm_real,EBR400(6,:),'LineWidth',3); 
plot(tx_power_dBm_real,EBR196(6,:),'LineWidth',3);
plot(tx_power_dBm_real,EBR100(6,:),'LineWidth',3);
grid on, axis tight;
xlabel('Transmitted Power [dBm]');
ylabel('Achievable Bit Rate [bits/s/Hz]');
title('Achievable Bit Rate vs Transmitted power with scatterers');
legend('Without RIS','400 unit elements','196 unit elements','100 unit elements');
hold off;

figure();
plot(tx_power_dBm_real,EBR400(1,:),'LineWidth',3); 
hold on;
plot(tx_power_dBm_real,EBR400(3,:),'LineWidth',3); 
plot(tx_power_dBm_real,EBR196(3,:),'LineWidth',3);
plot(tx_power_dBm_real,EBR100(3,:),'LineWidth',3);
grid on, axis tight;
xlabel('Transmitted Power [dBm]');
ylabel('Achievable Bit Rate [bits/s/Hz]');
title('Achievable Bit Rate vs Transmitted power without Cluster and 1 bit quantized phase values');
legend('Without RIS','400 unit elements','196 unit elements','100 unit elements');
hold off;

figure();
plot(tx_power_dBm_real,EBR400(1,:),'LineWidth',3); 
hold on;
plot(tx_power_dBm_real,EBR400(6,:),'LineWidth',3); 
plot(tx_power_dBm_real,EBR400(7,:),'LineWidth',3);
plot(tx_power_dBm_real,EBR400(8,:),'LineWidth',3);
grid on, axis tight;
xlabel('Transmitted Power [dBm]');
ylabel('Achievable Bit Rate [bits/s/Hz]');
title('Achievable Bit Rate vs Transmitted power with scatterers and 400 antenna elements');
legend('Without RIS','Continuous phase values','1 bit quantized phase values','2 bits quantized phase values');
hold off;

figure();
plot(tx_power_dBm_real,EBR400(5,:),'LineWidth',3); 
hold on;
plot(tx_power_dBm_real,EBR400(7,:),'LineWidth',3); 
plot(tx_power_dBm_real,EBR196(7,:),'LineWidth',3);
plot(tx_power_dBm_real,EBR100(7,:),'LineWidth',3);
grid on, axis tight;
xlabel('Transmitted Power [dBm]');
ylabel('Achievable Bit Rate [bits/s/Hz]');
title('Achievable Bit Rate vs Transmitted power without Cluster and 1 bit quantized phase values');
legend('Without RIS','400 unit elements','196 unit elements','100 unit elements');
hold off;

figure();
plot(tx_power_dBm_real,EBR400(5,:),'LineWidth',3); 
hold on;
plot(tx_power_dBm_real,EBR400(8,:),'LineWidth',3); 
plot(tx_power_dBm_real,EBR196(8,:),'LineWidth',3);
plot(tx_power_dBm_real,EBR100(8,:),'LineWidth',3);
grid on, axis tight;
xlabel('Transmitted Power [dBm]');
ylabel('Achievable Bit Rate [bits/s/Hz]');
title('Achievable Bit Rate vs Transmitted power without Cluster and 2 bits quantized phase values');
legend('Without RIS','400 unit elements','196 unit elements','100 unit elements');
hold off;
%%
% Tot_Channel_Cluster = transpose(g_wp) *  h_Cluster; % Tx-RIS-Rx with Cluster
% Tot_Channel_Clusterq = transpose(g_wp_q) * h_Cluster_q; % Tx-RIS-Rx with Cluster and 1 bit quantized phase
% Tot_Channel_Clusterqq = transpose(g_wp_qq) * h_Cluster_qq; % Tx-RIS-Rx with Cluster and 2 bit quantized phase
% Tot_Channel_wCluster = transpose(g_wp) * h_wCluster; % Tx-RIS-Rx without Cluster
% Tot_Channel_wClusterq = transpose(g_wp_q) * h_wCluster_q; % Tx-RIS-Rx without Cluster and 1 bit quantized phase
% Tot_Channel_wClusterqq = transpose(g_wp_qq) * h_wCluster_qq; % Tx-RIS-Rx without Cluster and 2 bit quantized phase
% h_imag_Cluster =  h_SISO_wp;
% h_imag_wCluster =  h_SISO_LOS_wp;
% %% Optimization problem
% f = @(x,THETA1)(-(transpose(g_wp)  * (exp(1i*x).*THETA1) * h_Cluster));
% THETA1 = diag(ones(N_ele^2,1));
% fun = @(x)f(x,THETA1);
% x = fminsearch(fun,THETA);
% %%
% Tot_Channel_Cluster = abs(transpose(g_wp)  * THETA * h_Cluster); % Tx-RIS-Rx with Cluster
% Tot_Channel_Clusterq = abs(transpose(g_wp_q)  * THETA * h_Cluster_q); % Tx-RIS-Rx with Cluster and 1 bit quantized phase
% Tot_Channel_Clusterqq = abs(transpose(g_wp_qq) * THETA * h_Cluster_qq); % Tx-RIS-Rx with Cluster and 2 bit quantized phase
% Tot_Channel_wCluster = abs(transpose(g_wp) * THETA * h_wCluster); % Tx-RIS-Rx without Cluster
% Tot_Channel_wClusterq = abs(transpose(g_wp_q) * THETA * h_wCluster_q); % Tx-RIS-Rx without Cluster and 1 bit quantized phase
% Tot_Channel_wClusterqq = abs(transpose(g_wp_qq) * THETA * h_wCluster_qq); % Tx-RIS-Rx without Cluster and 2 bit quantized phase
% h_imag_Cluster =  abs(h_SISO_wp);
% h_imag_wCluster =  abs(h_SISO_LOS_wp);

%%

% for p =1:NN
%     for z = 1:N_ele
%         for y = 1:N_ele
% 
%             zn = (z-1)*dy; ym = (y-1)*dy;
%             phi(z,y) = k0*(ym*sind(90-theta_rx(p))*sind(-phi_rx(p))+zn*cosd(90-theta_rx(p)));
% 
%             while phi(z,y)>=2*pi
%                 phi(z,y)=phi(z,y)-2*pi;
%             end
%             while phi(z,y)<=-2*pi
%                 phi(z,y)=phi(z,y)+2*pi;
%             end
% 
%             if phi(z,y)>pi/2 && phi(z,y)<=3*pi/2
%                 phiq(z,y) = pi;
%             else
%                 phiq(z,y) = 0;
%             end
% 
%             if phi(z,y) > pi/4 && phi(z,y) <= 3/4*pi
%                 phiqq(z,y) = pi/2;
%             elseif phi(z,y) > 5/4*pi && phi(z,y)<=3*pi/2
%                 phiqq(z,y) = 3*pi/2;
%             elseif phi(z,y)>3/4*pi && phi(z,y) <=5/4*pi
%                 phiqq(z,y) = pi;
%             else
%                 phiqq(z,y) = 0;
%             end
% 
%         end
%     end
%     clear pHI, clear pHIq, clear pHIqq;
%     pHI = phi(:)';
%     pHIq = phiq(:)';
%     pHIqq = phiqq(:)';
% 
%     clear tHETA, clear tHETAq, clear tHETAqq;
%     for k = 1:N_ele^2
%         tHETA(k,k) = exp(1i.*pHI(k));
%         tHETAq(k,k) = exp(1i.*pHIq(k));
%         tHETAqq(k,k) = exp(1i.*pHIqq(k));
%     end
%     Tot_Channel_wprx(p,:) = g_wprx(:,p)' * tHETA * h_wCluster; %Tx-RIS-Rx without Cluster for 20 points
%     Tot_Channel_wprxq(p,:) = g_wprx(:,p)' * tHETAq * h_wCluster; %Tx-RIS-Rx without Cluster for 20 points and quantized phase
%     Tot_Channel_wprxqq(p,:) = g_wprx(:,p)' * tHETAqq * h_wCluster; %Tx-RIS-Rx without Cluster for 20 points and quantized phase
% end


% %% channel model between Tx-Rx
% d_cs_tilde_imag=zeros(1,sum(S));
% h_SISO_NLOS_imag=0;
% h_SISO_NLOS_wp_imag=0;
% Rxx_loc_imag = [-RX_loc(1), RX_loc(2), RX_loc(3)];
% d_T_R_imag = norm(TX_loc-Rxx_loc_imag);
% 
% for counter2=indices
% 
%     % due to shared clusters d_cs_tilde ~ d_cs
%     d_cs_tilde_imag(counter2) = a_c_rep(counter2) + norm(Coordinates2(counter2,:)- Rxx_loc_imag);
% 
%     Lcs_dB_SISO_wp_imag=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs_tilde_imag(counter2))- shadow(counter2) - (rays{1}.PathLoss-mag2db(4*pi*d_T_R_imag/lambda));
%     Lcs_SISO_wp_imag=10^(Lcs_dB_SISO_wp_imag/10);
% 
%     % We consider the same complex path gain (small scale fading) with an excess phase (similar to array response)
%     eta_imag=k0* ( norm(Coordinates2(counter2,:)- RIS_loc) -  norm(Coordinates2(counter2,:)- Rxx_loc_imag));
% 
%     % h_SISO_NLOS_wp_imag = h_SISO_NLOS_wp_imag + beta(counter2)*exp(1i*eta)*sqrt(Lcs_SISO_wp_imag)*transpose(array_Rx_cs_SISO(counter2,:))*array_Tx_cs_SISO(counter2,:);
%     h_SISO_NLOS_wp_imag = h_SISO_NLOS_wp_imag + beta(counter2)*sqrt(Lcs_SISO_wp_imag)*transpose(array_Rx_cs_SISO(counter2,:))*array_Tx_cs_SISO(counter2,:);
% 
% end
% 
% h_SISO_NLOS_wp_imag=h_SISO_NLOS_wp_imag.*sqrt(1/M_new);
% 
% 
% L_SISO_LOS_dB_wp_imag=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_R_imag)- sigma_LOS - (rays{1}.PathLoss-mag2db(4*pi*d_T_R_imag/lambda));
% L_SISO_LOS_wp_imag=10^(L_SISO_LOS_dB_wp_imag/10);
% 
% I_phi_Tx_SISO_imag=sign(TX_loc(2)-Rxx_loc_imag(2));
% phi_Tx_SISO_imag = I_phi_Tx_SISO_imag* atand ( abs( TX_loc(2)-Rxx_loc_imag(2)) / abs(TX_loc(1)-Rxx_loc_imag(1)) );
% 
% I_theta_Tx_SISO_imag=sign(Rxx_loc_imag(3)-TX_loc(3));
% theta_Tx_SISO_imag= I_theta_Tx_SISO_imag* atand ( abs( Rxx_loc_imag(3)-TX_loc(3)) / abs(d_T_R_imag) );
% 
% % AoA of Rx for Tx-Rx channel in an Indoor
% phi_av_SISO_LOS = rand*180-90;     % mean azimuth
% theta_av_SISO_LOS= rand*180-90;      % mean elevation
% 
% h_SISO_LOS_wp_imag= sqrt(L_SISO_LOS_wp_imag)*Phi_var*transpose(array_Rx_SISO)*array_Tx_SISO;
% 
% h_SISO_wp_imagCluster=h_SISO_NLOS_wp_imag + h_SISO_LOS_wp_imag; %Tx-Rx without RIS
% h_SISO_wp_imagwCluster= h_SISO_LOS_wp_imag; %Tx-Rx without RIS

%%
figure;
RX_loc(1) = -RX_loc(1);
plot3(TX_loc(1),TX_loc(2),TX_loc(3),'ko','LineWidth',3);grid on;hold on;
plot3(RX_loc(1),RX_loc(2),RX_loc(3),'kx','LineWidth',3);
% for pp = 1: NN
%     plot3(xunit(pp),yunit(pp),RIS_loc(3),'x','Color','r');
%     %text(xunit(pp),yunit(pp),RIS_loc(3), 'Rx','FontSize',14);
% end

%plot3(ones(1,N_points)*RX_loc(1),y_near,z_near,'x','Color','r');
%plot3(xunit,yunit,ones(1,NN)*RIS_loc(3),'LineWidth',3,'Color','b');
%text(1.3,1,'Circular Path of Rx','Color','blue','FontSize',14);
plot3(RIS_loc(1),RIS_loc(2),RIS_loc(3),'ksquare','LineWidth',3)
xlabel('x direction [m]');ylabel('y direction [m]');zlabel('z direction [m]');
text(TX_loc(1),TX_loc(2),TX_loc(3),'  Tx','FontSize',18);
text(RX_loc(1),RX_loc(2),RX_loc(3),'  Rx','FontSize',18);
text(RIS_loc(1),RIS_loc(2),RIS_loc(3),'  RIS','FontSize',18);
axis([xLLC xURC yLLC yURC 0 h_max]);

xx1 = x1_obs*ones(1000); yy1 = y1_obs*ones(1000); zz1 = linspace(z1_obs,z2_obs,1000);
xx2 = x1_obs*ones(1000); yy2 = y2_obs*ones(1000); zz2 = linspace(z1_obs,z2_obs,1000);
xx3 = x2_obs*ones(1000); yy3 = y1_obs*ones(1000); zz3 = linspace(z1_obs,z2_obs,1000);
xx4 = x2_obs*ones(1000); yy4 = y2_obs*ones(1000); zz4 = linspace(z1_obs,z2_obs,1000);
%xx2 = 6.42*ones(1000); yy2 = linspace(2.94,3.18,1000); zz2 = linspace(0,h_max,1000);
plot3(xx1,yy1,zz1,'r','LineWidth',3); plot3(xx2,yy2,zz2,'r','LineWidth',3); plot3(xx3,yy3,zz3,'r','LineWidth',3); plot3(xx4,yy4,zz4,'r','LineWidth',3);
xx11 = linspace(x1_obs,x2_obs,1000); yy11 = y1_obs*ones(1000); zz11 = z1_obs*ones(1000);
xx22 = linspace(x1_obs,x2_obs,1000); yy22 = y2_obs*ones(1000); zz22 = z1_obs*ones(1000);
xx33 = linspace(x1_obs,x2_obs,1000); yy33 = y1_obs*ones(1000); zz33 = z2_obs*ones(1000);
xx44 = linspace(x1_obs,x2_obs,1000); yy44 = y2_obs*ones(1000); zz44 = z2_obs*ones(1000);
plot3(xx11,yy11,zz11,'r','LineWidth',3); plot3(xx22,yy22,zz22,'r','LineWidth',3); plot3(xx33,yy33,zz33,'r','LineWidth',3); plot3(xx44,yy44,zz44,'r','LineWidth',3);
xx111 = x1_obs*ones(1000); yy111 = linspace(y1_obs,y2_obs,1000); zz111 = z1_obs*ones(1000);
xx222 = x1_obs*ones(1000); yy222 = linspace(y1_obs,y2_obs,1000); zz222 = z2_obs*ones(1000);
xx333 = x2_obs*ones(1000); yy333 = linspace(y1_obs,y2_obs,1000); zz333 = z1_obs*ones(1000);
xx444 = x2_obs*ones(1000); yy444 = linspace(y1_obs,y2_obs,1000); zz444 = z2_obs*ones(1000);
plot3(xx111,yy111,zz111,'r','LineWidth',3); plot3(xx222,yy222,zz222,'r','LineWidth',3); plot3(xx333,yy333,zz333,'r','LineWidth',3); plot3(xx444,yy444,zz444,'r','LineWidth',3);
for counter=1:C
    plot3(Coordinates(counter,1),Coordinates(counter,2),Coordinates(counter,3),'mdiamond','Color','b','LineWidth',2);
end
for counter2=1:sum(S)
    plot3(Coordinates2(counter2,1),Coordinates2(counter2,2),Coordinates2(counter2,3),'m.','Color','g','LineWidth',10);
end
hold off;
title('Simulated environment with scatterers');




%% Modulation random bit-stream 64QAM/OFDM

scs = 30;
carrier = nrCarrierConfig('SubcarrierSpacing',scs);
p = 1;
pdsch = nrPDSCHConfig('NumLayers',p,'Modulation','64QAM');
[ind,info] = nrPDSCHIndices(carrier,pdsch);
numDataBits = info.G;
cws = randi([0 1],numDataBits,1);
%save("streamBit.mat","cws");
%cws = load("streamBit.mat");
%cws = cws.cws;
sym = nrPDSCH(carrier,pdsch,cws,'OutputDataType','single');
txGrid = nrResourceGrid(carrier,p);
txGrid(ind) = sym;
initialNSlot = carrier.NSlot;
cpl = 'extended';
[txWaveform,Info] = nrOFDMModulate(txGrid,scs,initialNSlot,'CyclicPrefix',cpl);
%% Tx Bit-stream  
nrb = carrier.NSizeGrid;
gridTx = nrOFDMDemodulate(txWaveform,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
bitTx = qamdemod(gridTx,64,'OutputType','bit');
%% Rx signal after circular convolution between TxSignal and impulse response of channel
waveform_Cluster = sqrt(tx_power_Corr).*(txWaveform*(Tot_Channel_Cluster)); %channel Tx-RIS-Rx with clusters
waveform_Clusterq = sqrt(tx_power_Corr).*(txWaveform*(Tot_Channel_Cluster_q)); %channel Tx-RIS-Rx with clusters and 1 bit quantized phase
waveform_Clusterqq = sqrt(tx_power_Corr).*(txWaveform*Tot_Channel_Cluster_qq); %channel Tx-RIS-Rx with clusters and 2 bit quantized phase
waveform_wCluster = sqrt(tx_power_Corr).*(txWaveform*Tot_Channel_wCluster); % channel Tx-RIS-Rx without Clusters
waveform_wClusterq = sqrt(tx_power_Corr).*(txWaveform*Tot_Channel_wCluster_q); % channel Tx-RIS-Rx without Clusters and 1 bit quantized phase
waveform_wClusterqq = sqrt(tx_power_Corr).*(txWaveform*Tot_Channel_wCluster_qq); % channel Tx-RIS-Rx without Clusters and 2 bit quantized phase
waveform_imagCluster = txWaveform*h_SISO_wp; %Tx-Rx with Cluster
waveform_imagwCluster = txWaveform*h_SISO_LOS_wp; %Tx-Rx without Cluster
% for p=1:NN
%     clear waveform_wprx, clear grid_wprx, clear bitRx_wprx, clear number_wprx, clear ratio_wprx;
%     clear waveform_wprxq, clear grid_wprxq, clear bitRx_wprxq, clear number_wprxq, clear ratio_wprxq;
%     clear waveform_wprxqq, clear grid_wprxqq, clear bitRx_wprxqq, clear number_wprxqq, clear ratio_wprxqq;
% 
%     waveform_wprx = conv2(txWaveform,Tot_Channel_wprx(p),'same'); %channel Tx-RIS-Rx without clusters 20 points
%     grid_wprx = nrOFDMDemodulate(waveform_wprx,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
%     bitRx_wprx = qamdemod(grid_wprx,64,'OutputType','bit');
%     [number_wprx,ratio_wprx] = biterr(bitTx,bitRx_wprx);
%     BitErrorRate_wprx(p) = ratio_wprx*100; %percentage of error bits
% 
%     waveform_wprxq = conv2(txWaveform,Tot_Channel_wprxq(p),'same'); %channel Tx-RIS-Rx without clusters 20 points and quantized phase
%     grid_wprxq = nrOFDMDemodulate(waveform_wprxq,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
%     bitRx_wprxq = qamdemod(grid_wprxq,64,'OutputType','bit');
%     [number_wprxq,ratio_wprxq] = biterr(bitTx,bitRx_wprxq);
%     BitErrorRate_wprxq(p) = ratio_wprxq*100; %percentage of error bits
% 
%     waveform_wprxqq = conv2(txWaveform,Tot_Channel_wprxqq(p),'same'); %channel Tx-RIS-Rx without clusters 20 points and quantized phase
%     grid_wprxqq = nrOFDMDemodulate(waveform_wprxqq,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
%     bitRx_wprxqq = qamdemod(grid_wprxqq,64,'OutputType','bit');
%     [number_wprxqq,ratio_wprxqq] = biterr(bitTx,bitRx_wprxqq);
%     BitErrorRate_wprxqq(p) = ratio_wprxqq*100; %percentage of error bits
% end

%% Demodulation OFDM of RxSignal

grid_Cluster = nrOFDMDemodulate(waveform_Cluster,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_Clusterq = nrOFDMDemodulate(waveform_Clusterq,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_Clusterqq = nrOFDMDemodulate(waveform_Clusterqq,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_wCluster = nrOFDMDemodulate(waveform_wCluster,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_wClusterq = nrOFDMDemodulate(waveform_wClusterq,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_wClusterqq = nrOFDMDemodulate(waveform_wClusterqq,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_imagCluster = nrOFDMDemodulate(waveform_imagCluster,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_imagwCluster = nrOFDMDemodulate(waveform_imagwCluster,nrb,scs,initialNSlot,'CyclicPrefix',cpl);

%% Demodulation QAM of RxSignal 
bitRx_Cluster = qamdemod(grid_Cluster,64,'OutputType','bit');
bitRx_Clusterq = qamdemod(grid_Clusterq,64,'OutputType','bit');
bitRx_Clusterqq = qamdemod(grid_Clusterqq,64,'OutputType','bit');
bitRx_wCluster = qamdemod(grid_wCluster,64,'OutputType','bit');
bitRx_wClusterq = qamdemod(grid_wClusterq,64,'OutputType','bit');
bitRx_wClusterqq = qamdemod(grid_wClusterqq,64,'OutputType','bit');
bitRx_imagCluster = qamdemod(grid_imagCluster,64,'OutputType','bit');
bitRx_imagwCluster = qamdemod(grid_imagwCluster,64,'OutputType','bit');
%% difference between Tx and Rx bit
[number_Cluster,ratio_Cluster] = biterr(bitTx,bitRx_Cluster);
BitErrorRate_Cluster = ratio_Cluster*100; %percentage of error bits
[number_Clusterq,ratio_Clusterq] = biterr(bitTx,bitRx_Clusterq);
BitErrorRate_Clusterq = ratio_Clusterq*100; %percentage of error bits
[number_Clusterqq,ratio_Clusterqq] = biterr(bitTx,bitRx_Clusterqq);
BitErrorRate_Clusterqq = ratio_Clusterqq*100; %percentage of error bits
[number_wCluster,ratio_wCluster] = biterr(bitTx,bitRx_wCluster);
BitErrorRate_wCluster = ratio_wCluster*100; %percentage of error bits
[number_wClusterq,ratio_wClusterq] = biterr(bitTx,bitRx_wClusterq);
BitErrorRate_wClusterq = ratio_wClusterq*100; %percentage of error bits
[number_wClusterqq,ratio_wClusterqq] = biterr(bitTx,bitRx_wClusterqq);
BitErrorRate_wClusterqq = ratio_wClusterqq*100; %percentage of error bits
[number_imagCluster,ratio_imagCluster] = biterr(bitTx,bitRx_imagCluster);
BitErrorRate_imagCluster = ratio_imagCluster*100; %percentage of error bits
[number_imagwCluster,ratio_imagwCluster] = biterr(bitTx,bitRx_imagwCluster);
BitErrorRate_imagwCluster = ratio_imagwCluster*100; %percentage of error bits

disp('Bit Error Rate in different cases')
disp(['-> Channel Tx-RIS-Rx with clusters                           = ', num2str(BitErrorRate_Cluster), ' %. ']);
disp(['-> Channel Tx-RIS-Rx with clusters and 1 bit quantized phase       = ', num2str(BitErrorRate_Clusterq), ' %. ']);
disp(['-> Channel Tx-RIS-Rx with clusters and 2 bit quantized phase       = ', num2str(BitErrorRate_Clusterqq), ' %. ']);
disp(['-> Channel Tx-RIS-Rx without clusters                           = ', num2str(BitErrorRate_wCluster), ' %. ']);
disp(['-> Channel Tx-RIS-Rx without clusters and 1 bit quantized phase       = ', num2str(BitErrorRate_wClusterq), ' %. ']);
disp(['-> Channel Tx-RIS-Rx without clusters and 2 bit quantized phase       = ', num2str(BitErrorRate_wClusterqq), ' %. ']);
disp(['-> Channel Tx-Rx with clusters                               = ' , num2str(BitErrorRate_imagCluster), ' %. ']);
disp(['-> Channel Tx-Rx without clusters                               = ' , num2str(BitErrorRate_imagwCluster), ' %. ']);
%%
% Pn = -100; %dB
% 
% SignalRx_C = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_Cluster)^2));
% SignalRx_Cq = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_Clusterq)^2));
% SignalRx_Cqq= tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_Clusterqq)^2));
% 
% SignalRx_wC = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_wCluster)^2));
% SignalRx_wCq = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_wClusterq)^2));
% SignalRx_wCqq =tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_wClusterqq)^2));
% 
% SignalRx_imagwC =tx_power_dBm_real+ (10*log10(abs(h_imag_wCluster)^2));
% SignalRx_imagC = tx_power_dBm_real+(10*log10(abs(h_imag_Cluster)^2));
% disp('Received power ');
% disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster = ', num2str(SignalRx_C), ' dBm']);
% disp(['-> Tx-RIS-Rx Received Signal Strength with Clusterand 1 BIT Quantization = ', num2str(SignalRx_Cq), ' dBm']);
% disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster 2 BIT Quantization = ', num2str(SignalRx_Cqq), ' dBm']);
% 
% disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster = ', num2str(SignalRx_wC), ' dBm']);
% disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster 1 BIT Quantization = ', num2str(SignalRx_wCq), ' dBm']);
% disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster 2 BIT Quantization = ', num2str(SignalRx_wCqq), ' dBm']);
% 
% disp(['-> Tx-Rx Received Signal Strength with Cluster = ', num2str(SignalRx_imagC), ' dBm']);
% disp(['-> Tx-Rx Received Signal Strength without Cluster = ', num2str(SignalRx_imagwC), ' dBm']);
%%
Pn = -100;
SNR_C = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_Cluster)^2))-Pn;
SNR_Cq = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_Cluster_q)^2))-Pn;
SNR_Cqq= tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_Cluster_qq)^2))-Pn;

SNR_wC = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_wCluster)^2))-Pn;
SNR_wCq = tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_wCluster_q)^2))-Pn;
SNR_wCqq =tx_power_Corr_dBm+(10*log10(abs(Tot_Channel_wCluster_qq)^2))-Pn;

SNR_imagwC =tx_power_dBm_real+ (10*log10(abs(h_SISO_LOS_wp)^2))-Pn;
SNR_imagC = tx_power_dBm_real+(10*log10(abs(h_SISO_wp)^2))-Pn;
disp('Signal Noise Ratio with Pn = -100 dBm ');
disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster = ', num2str(SNR_C), ' dB']);
disp(['-> Tx-RIS-Rx Received Signal Strength with Clusterand 1 BIT Quantization = ', num2str(SNR_Cq), ' dB']);
disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster 2 BIT Quantization = ', num2str(SNR_Cqq), ' dB']);

disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster = ', num2str(SNR_wC), ' dB']);
disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster 1 BIT Quantization = ', num2str(SNR_wCq), ' dB']);
disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster 2 BIT Quantization = ', num2str(SNR_wCqq), ' dB']);

disp(['-> Tx-Rx Received Signal Strength with Cluster = ', num2str(SNR_imagC), ' dB']);
disp(['-> Tx-Rx Received Signal Strength without Cluster = ', num2str(SNR_imagwC), ' dB']);

R_imagwC = floor((mean(log2(1+SNR_imagwC))));
R_imagC = floor((mean(log2(1+SNR_imagC))));
R_C = floor((mean(log2(1+SNR_C))));
R_Cq = floor((mean(log2(1+SNR_Cq))));
R_Cqq = floor((mean(log2(1+SNR_Cqq))));
R_wC = floor((mean(log2(1+SNR_wC))));
R_wCq = floor((mean(log2(1+SNR_wCq))));
R_wCqq = floor(mean(log2(1+SNR_wCqq)));

% R_imagwC = floor((mean(log2(1+(tx_power_dBm_real+ (10*log10(abs(h_SISO_LOS_wp)^2))-Pn)))));
% R_imagC = floor((mean(log2(1+(tx_power_dBm_real+ (10*log10(abs(h_SISO_wp)^2))-Pn)))));
% R_C = floor(abs((mean(log2(1+(tx_power_Corr_dBm+(20*log10(abs(mean(transpose(g_wp)) *  mean(h_Cluster))^2))-Pn))))));
% R_Cq = floor(abs((mean(log2(1+(tx_power_Corr_dBm+(10*log10(abs(mean(transpose(g_wp_q)) * mean(h_Cluster))^2))-Pn))))));
% R_Cqq = floor(abs((mean(log2(1+(tx_power_Corr_dBm+(10*log10(abs(mean(transpose(g_wp_qq)) * mean(h_Cluster))^2))-Pn))))));
% R_wC = floor(abs((mean(log2(1+(tx_power_Corr_dBm+(10*log10(abs(mean(transpose(g_wp)) * mean(h_wCluster))^2))-Pn))))));
% R_wCq = floor(abs((mean(log2(1+(tx_power_Corr_dBm+(10*log10(abs(mean(transpose(g_wp_q)) * mean(h_wCluster))^2))-Pn))))));
% R_wCqq = floor(abs((mean(log2(1+(tx_power_Corr_dBm+(10*log10(abs(mean(transpose(g_wp_qq)) * mean(h_wCluster))^2))-Pn))))));

disp('Ergodic Achievable rate');
disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster = ', num2str(R_C), ' bits/s/Hz']);
disp(['-> Tx-RIS-Rx Received Signal Strength with Clusterand 1 BIT Quantization = ', num2str(R_Cq), ' bits/s/Hz']);
disp(['-> Tx-RIS-Rx Received Signal Strength with Cluster 2 BIT Quantization = ', num2str(R_Cqq), ' bits/s/Hz']);

disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster = ', num2str(R_wC), ' bits/s/Hz']);
disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster 1 BIT Quantization = ', num2str(R_wCq), ' bits/s/Hz']);
disp(['-> Tx-RIS-Rx Received Signal Strength without Cluster 2 BIT Quantization = ', num2str(R_wCqq), ' bits/s/Hz']);

disp(['-> Tx-Rx Received Signal Strength with Cluster = ', num2str(R_imagC), ' bits/s/Hz']);
disp(['-> Tx-Rx Received Signal Strength without Cluster = ', num2str(R_imagwC), ' bits/s/Hz']);

%% 
%creare una tabella in cui inserisco le 20 posizioni di test e mettiamo in
%una colonna la banda, nella seconda l'angolo min, nella terza l'angolo
%max, nella quarta il valore effettivo dell'angolo di ricezione. 
varType = ["double","double","double","double"];
varNames = ["Iteration Number", "Bit Error Rate", "Bit Error Rate 1 BIT", "Bit Error Rate 2 BIT"];
T = table('Size',[20,4],'VariableTypes',varType,'VariableNames',varNames);
T{:,1} = linspace(1,20,20)';
T{:,2} = BitErrorRate_wprx'; T{:,3} = BitErrorRate_wprxq'; T{:,4} = BitErrorRate_wprxqq';
T
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
Pw_dBm_Cluster = 10*log10(abs(Tot_Channel_Cluster).^2);
Pw_dBm_Clusterq = 10*log10(abs(Tot_Channel_Clusterq).^2);
Pw_dBm_Clusterqq = 10*log10(abs(Tot_Channel_Clusterqq).^2);
Pw_dBm_wCluster = 10*log10(abs(Tot_Channel_wCluster).^2);
Pw_dBm_wClusterq = 10*log10(abs(Tot_Channel_wClusterq).^2);
Pw_dBm_wClusterqq = 10*log10(abs(Tot_Channel_wClusterqq).^2);
Pw_dBm_imagCluster = 10*log10(abs(h_SISO_wp_imagCluster).^2);
Pw_dBm_imagwCluster = 10*log10(abs(h_SISO_wp_imagwCluster).^2);

%Pw_dBm_H = 10*log10(abs(TOT_CHANNEL).^2); %Power associated to the channel (pag. 193 BOOK -> Teoria dei segnali)

Pw_dBm_rx_Cluster = tx_power_Corr_dBm + Pw_dBm_Cluster;
Pw_dBm_rx_Clusterq = tx_power_Corr_dBm + Pw_dBm_Clusterq;
Pw_dBm_rx_Clusterqq = tx_power_Corr_dBm + Pw_dBm_Clusterqq;
Pw_dBm_rx_wCluster = tx_power_Corr_dBm + Pw_dBm_wCluster;
Pw_dBm_rx_wClusterq = tx_power_Corr_dBm + Pw_dBm_wClusterq;
Pw_dBm_rx_wClusterqq = tx_power_Corr_dBm + Pw_dBm_wClusterqq;
Pw_dBm_rx_imagCluster = tx_power_Corr_dBm + Pw_dBm_imagCluster;
Pw_dBm_rx_imagwCluster = tx_power_Corr_dBm + Pw_dBm_imagwCluster;

disp('Received power considering transmitted power');
disp(['Cluster and continuous phase = ' , num2str(Pw_dBm_rx_Cluster),' dBm']);
disp(['Cluster and 1 bit quantized phase = ' , num2str(Pw_dBm_rx_Clusterq),' dBm']);
disp(['Cluster and 2 bit quantized phase = ' , num2str(Pw_dBm_rx_Clusterqq),' dBm']);
disp(['Without cluster and continuous phase = ' , num2str(Pw_dBm_rx_wCluster),' dBm']);
disp(['Without cluster and 1 bit quantized phase = ' , num2str(Pw_dBm_rx_wClusterq),' dBm']);
disp(['Without cluster and 2 bit quantized phase = ' , num2str(Pw_dBm_rx_wClusterqq),' dBm']);
disp(['Tx-Rx with cluster = ', num2str(Pw_dBm_rx_imagCluster), ' dBm']);
disp(['Tx-Rx without cluster = ', num2str(Pw_dBm_rx_imagwCluster), ' dBm']);

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

