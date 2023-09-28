%%
clear all
close all
clc
%%
numTx = 1;
numRx = 1;
%% Visualizzazione stanza
filename = "Ingenious Juttuli-Amberis (1).stl"; % Nome file formato stl
stldata = stlread(filename);
%siteviewer("SceneModel","Ingenious Juttuli-Amberis (1).stl");

fontsize = 13;
equalaspect = true;
% scomposizione in lista di connettività e punti
M = stldata.ConnectivityList;
x0 = stldata.Points(:,1);
y0 = stldata.Points(:,2);
z0 = stldata.Points(:,3);

% scalatura valori
scale = 3e-2;
RIS_px = [0.0332; 3; 1.2]./scale;

xmin = min(x0);
ymin =min(y0);
zmin = min(z0);

x = x0 - xmin-RIS_px(1);
y = y0 - ymin-RIS_px(2);
z = z0- zmin-RIS_px(3);

Stanza = triangulation(M,x.*scale, y.*scale, z.*scale);
[~, axes1] = get_figure_axes(fontsize, equalaspect);
hold(axes1, 'on');
grid(axes1, 'on');
trisurf(Stanza,'FaceAlpha', 0.3, 'EdgeAlpha', 0.2, 'Parent', axes1);
xlabel(axes1, 'X [m]');
ylabel(axes1, 'Y [m]');
zlabel(axes1, 'Z [m]');
hold(axes1, 'off');
view([-45 30]);
axis equal

% Definizione estremi della stanza
xLLC = -2e-5; yLLC = -2.97;
xURC = 6.39221 ; yURC = 3.41855;
% oggetti x y dx dy
rects = [2.99  -0.06  3.4 0.24];

h_max = 3.5;
[~, axes1] = get_figure_axes(fontsize, equalaspect);
hold(axes1, 'on'); grid(axes1, 'on'); box(axes1, 'on');
trisurf(Stanza,'FaceAlpha', 0.3, 'EdgeAlpha', 0.2, 'Parent', axes1);
view([0 90]);

% Perimetrazione delle zone da escludere 
for ii = 1 :size(rects,1)
    rectangle(axes1, 'Position',rects(ii,:), 'EdgeColor','r', 'LineWidth',2);
end

xlabel(axes1, 'X [m]');
ylabel(axes1, 'Y [m]');
zlabel(axes1, 'Z [m]');
hold(axes1, 'off');

%% Dell'area rimasta se ne valuta la griglia --> tutte le posizioni calpestabili

res = .10;          %[m] risoluzione griglia 40 cm
antennaHeight = 0;  %[m] altezza dispositivo 1 m dal suolo

% definizione griglia
x_surf = (xLLC : res : xURC );
y_surf = (yLLC : res : yURC );
[X, Y] = meshgrid(x_surf, y_surf);
Z = antennaHeight.*ones(size(X));
%matrice con coordintate in basso a sinista e in alto a destra di ogni rettangolo
rectsabs = [rects(:,1), rects(:,2), rects(:,1) + rects(:,3), rects(:,2) + rects(:,4)];
mesh_nan_indexing = zeros(size(X,1), size(X, 2), size(rectsabs,1));
for ii = 1 : size(rectsabs,1)
    mesh_nan_indexing(:,:, ii) = X>=rectsabs(ii,1)&X<=rectsabs(ii,3)&Y>=rectsabs(ii,2)&Y<=rectsabs(ii,4);
end
mesh_nan_indices = any(mesh_nan_indexing,3);
%imagesc(mesh_nan_indices); colorbar;
Z(mesh_nan_indices) = NaN;
[~, axes1] = get_figure_axes(fontsize,equalaspect);
hold(axes1, 'on'); grid(axes1, 'on'); box(axes1, 'on');
trisurf(Stanza,'FaceAlpha', 0.3, 'EdgeAlpha', 0.2, 'Parent', axes1);
surf(X,Y,Z, 'Parent',axes1);
xlabel(axes1, 'X [m]');
ylabel(axes1, 'Y [m]');
zlabel(axes1, 'Z [m]');
hold(axes1, 'off');
view([0 90]);

%% Definizione dell'elemento base della RIS

c = physconst('LightSpeed');
f = 28e9; % [GHz]
Frequency = 28e9; %[GHz], valore modificabile di frequenza 
B= .5e9;
freq = f-B/2:.1e9:f+B/2;
lambda = c/f;
k0 = 2*pi/lambda; %numero dx'onda
Zopt = 50; %impedenza ottima dell'antenna

eps_r = 3.3;
h = 0.5e-3; % [m]

%Patch width

W = c/(2*f)*sqrt(2/(eps_r+1));

%Effective dielectric constant
eps_reff_r = (eps_r+1)/2 + (((eps_r-1)/2)*(1+12*h/W)^(-0.5));

%Extension length
DL = 0.412*h*((eps_reff_r+0.3)*(W/h+0.264))/((eps_reff_r-0.258)*(W/h+0.8));

%Patch length
L = c/(2*f*sqrt(eps_reff_r))-2*DL;

%Effective length of the patch
Leff = L +2*DL;

%Ground length
LG = 6*h+L;

%Ground width
WG = 6*h+W;

%Resonant frequency of dominant TM010 mode
fr = c/(2*L*sqrt(eps_r));
%Resonant frequency of dominant TM010 mode to include edge effects
frc = c/(2*Leff*sqrt(eps_r));
%Fringe factor
q = frc/fr;

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

for i = 1:length(y)
    Rin(i) = Rin0*(cos(pi/L*y(i))^2);
end

% trovare il valore più prossimo a 50 ohm

Diff = abs(Rin - ones(1,length(y)).*Zopt);
[m1 index_m1] = min(Diff(1:floor(length(Diff)/2)));
[m2 index_m2] = min(Diff(floor(length(Diff)/2):end));
Z50_1 = Rin(index_m1);
y50_1 = y(index_m1);
Z50_2 = Rin(index_m2 + round(length(Diff)/2)-1);
y50_2 = y(index_m2 + round(length(Diff)/2)-1);
Z50 = [Z50_1 Z50_2];
y50 = [y50_1 y50_2];

% Plot of the find points

% figure(); plot(y/L,Rin/Rin0);
% grid on; axis tight; 
% xlabel('y/L');
% ylabel('Rin/Rin0');
% title('Normalized input resistance');
% hold on

% plot(y/L,(ones(1,length(y))).*Zopt/Rin0);
% plot(y50/L, Z50/Rin0,'gx');
% plot([y50_1/L y50_1/L] ,[1 0],'r--','LineWidth',2)
% text(y50_1/L,0.73,'\leftarrow First Value at 50 Ohm /','Interpreter','tex',...
%     'Color','r','FontSize',10);
% plot([y50_2/L y50_2/L] ,[1 0],'r--','LineWidth',2)
% text(.49,0.73,'Second Value at 50 Ohm \leftarrow ','Interpreter','tex',...
%     'Color','r','FontSize',10);
% legend ('Normalized input resistance(y)' , 'Line at 50 Ohm' ,  'Point at 50 Ohm');
% hold off

%
% disp(['Value of impedance close to 50 Ohm is ', num2str(Z50_1), ' Ohm']);
% disp(['Value of distance from y0 = 0 is y = ' num2str(y50_1), ' meters']);

% Imponiamo il feed point in corrispondenza del valore ottimale dell'impedenza diingresso

%theoretical model
ant.FeedOffset = [y50(1)- L/2 , 0];

p = PatternPlotOptions;
p.Transparency = 0.5;
%figure(); pattern(ant,f,'patternOptions',p); %Radiation pattern

%ant.TiltAxis='y';
ant.Tilt=0;
G0 = pattern(ant,f,'Type','realizedgain');
maxG0 = max(G0(:));
disp(['Max antenna TX/RX Gain: ' num2str(maxG0) '[dBi]']);
ant.Tilt = 90; ant.TiltAxis = [0 1 0];
%figure(); pattern(ant,f,'patternOptions',p);

%% interdistanza e numero elementin RIS
dx = lambda/2;
dy = lambda/2;
N_ele = 8; %potenza del 2

%% Posizione Tx e Rx conoscendo angoli di incidenza, riflessione e distanza
theta_Rx = -1.7326; %elevation receiver angle considering the center in the RIS (posso considerare anche il primo elemento in alto a sinistra)
theta_Tx = -1.1467; %elevation transmitter angle considering the center in the RIS

phi_Tx = -25.7444;
phi_Rx = 14.5407;
d_Tx = 4.9969;
d_Rx = 3.3074;

xTx = d_Tx*cosd(phi_Tx);
yTx = d_Tx*sind(phi_Tx)*cosd(theta_Tx);
zTx = d_Tx*sind(theta_Tx)*sind(phi_Tx);

xRx = d_Rx*cosd(phi_Rx);
yRx = d_Rx*sind(phi_Rx)*cosd(theta_Rx);
zRx = d_Rx*sind(theta_Rx)*sind(phi_Rx);

% xRx = d_Rx*sind(theta_Rx)*cosd(phi_Rx);
% yRx = d_Rx*sind(theta_Rx)*sind(phi_Rx);
% zRx = d_Rx*cosd(theta_Rx);

% [xTx,yTx,zTx] = sph2cart(phi_Tx*pi/180,elevation_Tx*pi/180,d_Tx);
Tx_loc = [xTx, yTx, zTx];
% [xRx,yRx,zRx] = sph2cart(phi_Rx*pi/180,elevation_Rx*pi/180,d_Rx);
Rx_loc = [xRx, yRx, zRx];

%% Posizione Tx e Rx conoscendo le reali posizioni
xRx = 3.2; yRx= 0.83 ; zRx= -0.1;
xTx = 4.5; yTx= -2.17; zTx= -0.1;
TX_loc = [xTx, yTx, zTx];
RX_loc = [xRx, yRx, zRx];

xRIS = 0; yRIS = 0; zRIS =0; RIS_xyz = [xRIS yRIS zRIS];
d_T_RIS = norm(TX_loc-RIS_xyz); 
d_RIS_R=norm(RIS_xyz-RX_loc);
% rTX = sqrt(xTx^2+yTx^2+zTx^2);
% phi_TX = atand((yTx-yRIS)/(xTx-xRIS)); %angolo di partenza dal Tx
% comp_theta_TX = acosd((zTx-zRIS)/(sqrt(xTx^2+yTx^2+zTx^2)));
% theta_TX = 90-comp_theta_TX;

% rRX = sqrt(xRx^2+yRx^2+zRx^2);
% phi_RX = atand(yRx/xRx);
% comp_theta_RX = acosd(zRx/(sqrt(xRx^2+yRx^2+zRx^2)));
% theta_RX = 90-comp_theta_RX;

I_phi=sign(yTx-yRIS);
phi_TX = I_phi* atand ( abs( yRIS-yTx) / abs(xRIS-xTx) ); %Tra RIS e Tx
   
I_theta=sign(zTx-zRIS);
theta_TX=I_theta * asind ( abs (zRIS-zTx ) / d_T_RIS );

I_theta=sign(zRx - zRIS);
theta_RX=I_theta * asind( abs(zRx-zRIS)/d_RIS_R ); % AoD of RIS

% Azimuth Departure Angle
I_phi=sign(yRx - yRIS);
phi_RX=I_phi * atand( abs(yRx-yRIS)/ abs(xRx-xRIS) );

%% Posizione relativa del ricevitore rispetto i vari elementi della RIS
%come origine è stato scelto il primo elemento in alto a sinistra

% verifica della posizione del Tx e Rx se è già occupato da un oggetto
if (xTx >=rects(1,1) && xTx<=rects(1,1)+rects(1,3)) && (yTx >=rects(1,2) && yTx<=rects(1,2)+rects(1,4))
    disp('ERRORE: LA POSIZIONE DEL TRASMETTITORE COINCIDE CON LA POSIZIONE DELLA PARETE');
end
if (xRx >=rects(1,1) && xRx<=rects(1,1)+rects(1,3)) && (yRx >=rects(1,2) && yRx<=rects(1,2)+rects(1,4))
    disp('ERRORE: LA POSIZIONE DEL RICEVITORE COINCIDE CON LA POSIZIONE DELLA PARETE');
end

for yy=1:N_ele
    for zz=1:N_ele
        Tx_loc(zz,yy) = {[xTx, -(yy-1)*dx+yTx, -(zz-1)*dy+zTx]};
        Rx_loc(zz,yy) = {[xRx, -(yy-1)*dx+yRx, -(zz-1)*dy+zRx]};
        
        rTx(zz,yy) = sqrt(xTx^2+((yy-1)*dx+yTx)^2+((zz-1)*dy+zTx)^2);
        phi_Tx(zz,yy) = sign(((yy-1)*dx+yTx)-yRIS)*atand(abs(yRIS-((yy-1)*dx+yTx))/abs(xRIS-xTx));
       
        theta_Tx(zz,yy) =  sign(((zz-1)*dy+zTx)-zRIS)*asind(abs(zRIS-((zz-1)*dy+zTx)/norm([xTx, -(yy-1)*dx+yTx, -(zz-1)*dy+zTx]-RIS_xyz )));

        rRx(zz,yy) = sqrt(xRx^2+((yy-1)*dx+yRx)^2+((zz-1)*dy+zRx)^2);
        phi_Rx(zz,yy) = sign(((yy-1)*dx+yRx)-yRIS)*atand(abs(((yy-1)*dx+yRx)-yRIS)/abs(xRx-xRIS));
        theta_Rx(zz,yy) = sign(((zz-1)*dy+zRx)-zRIS)*asind(abs(((zz-1)*dy+zRx)-zRIS)/norm(RIS_xyz-[xRx, -(yy-1)*dx+yRx, -(zz-1)*dy+zRx]));
    end
end

%% Calcolo campo elettrico dell'onda riflessa dalla RIS al ricevitore
%nelle condizioni di far field posso considerare gli angoli di incidenza e
%riflessione uguali. Se supponiamo che theta sia 0 o 90° il campo è nullo
%(NON CORRETTO)
Ei = 1; %campo elettrico incidente
Er = 0;
% Er1 = 0; Er2 = 0;
Gamma_dB = -3*rand(N_ele,N_ele); %dB
Gamma = exp(Gamma_dB/10);
clear i;
% theta_RX_ = linspace(-90,90,181);
% phi_RX_ = linspace(-180,180,361);
% Er1_ = zeros(length(theta_RX_),length(phi_RX_));
% Er2_ = zeros(length(theta_RX_),length(phi_RX_));
% for n = 1:N_ele
%     for m = 1:N_ele
%         PHI1(m,n) = -k0*sind(theta_Rx(m,n))*((m-1)*dx*cosd(phi_Rx(m,n))+(n-1)*dy*sind(phi_Rx(m,n))) + k0*sind(theta_Tx(m,n))*((m-1)*dx*cosd(phi_Tx(m,n))+(n-1)*dy*sind(phi_Tx(m,n)));
%         PHI2(m,n) = -k0*sind(theta_Rx(m,n))*((m-1)*dx*cosd(phi_Rx(m,n))+(n-1)*dy*sind(phi_Rx(m,n)))+k0*sqrt(xTx.^2 + (-(m-1)*dx+yTx).^2 +(-(n-1)*dx+zTx).^2);
%         while PHI1(m,n)>=2*pi
%             PHI1(m,n)=PHI1(m,n)-2*pi;
%         end
% 
%         PHI1_comp(m,n) = cos(PHI1(m,n))+i*sin(PHI1(m,n));
% 
%         if PHI1(m,n)> -pi/2 && PHI1(m,n)< pi/2
%             PHI1_q(m,n) = 0;
%         else
%             PHI1_q(m,n) = 180;
%         end
% 
%         while PHI2(m,n)>=2*pi
%             PHI2(m,n)=PHI2(m,n)-2*pi;
%         end
% 
%         PHI2_comp(m,n) = cos(PHI2(m,n))+i*sin(PHI2(m,n));
% 
%         if PHI2(m,n)> -pi/2 && PHI2(m,n)< pi/2
%             PHI2_q(m,n) = 0;
%         else
%             PHI2_q(m,n) = 180;
%         end
% 
%         Er1 = Er1+ Gamma(m,n)*(cosd(PHI1_q(m,n))+i*sind(PHI1_q(m,n)))*Ei*exp((-i*k0*sind(theta_Rx(m,n)))*((m-1)*dx*cosd(phi_Rx(m,n))+(n-1)*dy*sind(phi_Rx(m,n))));
%         Er2 = Er2+ Gamma(m,n)*(cosd(PHI2_q(m,n))+i*sind(PHI2_q(m,n)))*Ei*exp((-i*k0*sind(theta_Rx(m,n)))*((m-1)*dx*cosd(phi_Rx(m,n))+(n-1)*dy*sind(phi_Rx(m,n))));
%     end
% end
% 
% Er1 = Er1*cosd(theta_Rx(m,n))*cosd(theta_Tx(m,n));
% Er2 = Er2*cosd(theta_Rx(m,n))*cosd(theta_Tx(m,n));
% 
% for k = 1:length(theta_RX_)
%     for kk = 1:length(phi_RX_)
%         for n = 1:N_ele
%             for m = 1:N_ele
%                 Er1_(k,kk) = Er1_(k,kk)+ Gamma(m,n)*(cosd(PHI1_q(m,n))+i*sind(PHI1_q(m,n)))*Ei*exp((-i*k0*sind(theta_RX_(k)))*((m-1)*dx*cosd(phi_RX_(kk))+(n-1)*dy*sind(phi_RX_(kk))));
%                 Er2_(k,kk) = Er2_(k,kk)+ Gamma(m,n)*(cosd(PHI2_q(m,n))+i*sind(PHI2_q(m,n)))*Ei*exp((-i*k0*sind(theta_RX_(k)))*((m-1)*dx*cosd(phi_RX_(kk))+(n-1)*dy*sind(phi_RX_(kk))));
%             end
%         end
% 
%         Er1_(k,kk) = Er1_(k,kk)*cosd(theta_RX_(k)-phi_RX)*cosd(theta_TX);
%         Er2_(k,kk) = Er2_(k,kk)*cosd(theta_RX_(k)-phi_RX)*cosd(theta_TX);
%     end
% end




for n = 1:N_ele
    for m = 1:N_ele
        %PHI(m,n) = -k0*sind(theta_Rx(m,n))*((m-1)*dx*cosd(phi_Rx(m,n))+(n-1)*dy*sind(phi_Rx(m,n))) + k0*sind(theta_Tx(m,n))*((m-1)*dx*cosd(phi_Tx(m,n))+(n-1)*dy*sind(phi_Tx(m,n)));
        PHI(m,n) = -k0*sind(theta_Rx(m,n))*((m-1)*dx*cosd(phi_Rx(m,n))+(n-1)*dy*sind(phi_Rx(m,n)))+k0*sqrt(xTx.^2 + (-(m-1)*dx+yTx).^2 +(-(n-1)*dx+zTx).^2);
        while PHI(m,n)>=2*pi
            PHI(m,n)=PHI(m,n)-2*pi;
        end

        PHI_comp(m,n) = cos(PHI(m,n))+i*sin(PHI(m,n));

        if PHI(m,n)> -pi/2 && PHI(m,n)< pi/2
            PHI_q(m,n) = 0;
        else
            PHI_q(m,n) = 180;
        end
        Er = Er+ Gamma(m,n)*(cosd(PHI_q(m,n))+i*sind(PHI_q(m,n)))*Ei*exp((-i*k0*sind(theta_Rx(m,n)))*((m-1)*dx*cosd(phi_Rx(m,n))+(n-1)*dy*sind(phi_Rx(m,n))));
    end
end

Er = Er*cosd(theta_Rx(m,n))*cosd(theta_Tx(m,n));

%% RIS Phased URA (si può usare solo nelle condizioni di far field,
% perchè non posso settare un valore di fase differente per i vari elementi)
RIS = phased.URA("Element",ant, ...
    "Size", [N_ele N_ele], ...
    "ElementSpacing",[dx dy],"ArrayNormal","x");

%Scan the antenna beam by applying a taper for a range of angles.
% For each angle, update the radiation pattern in Site Viewer.
% This approach of scanning the beam produces different patterns than physically
% rotating the antenna, as could be achieved by setting AntennaAngle of the
% transmitter site. This step is used to validate the orientation of the antenna's
% main beam.

% Get the starting array taper
taperedRISt = clone(RIS);
taperedRISt_quant = clone(RIS);
startTaper_t = taperedRISt.Taper;
nbar = 2; 
sll = -30;
steeringVector_t = phased.SteeringVector("SensorArray",taperedRISt,'IncludeElementResponse',true);

taperedRISr = clone(RIS);
startTaper_r = taperedRISr.Taper;
steeringVector_r = phased.SteeringVector("SensorArray",taperedRISr);
%% 
% The azimuth angle must be between –180° and 180°, and the elevation angle must be between –90° and 90°
% Define angles over which to perform sweep
azsweep = -180:15:180;
elsweep = -90:15:90;

% Set up tapering window and steering vector
N = N_ele*N_ele;
sltaper = taylorwin(N,nbar,sll)';


% Sweep the angles and show the antenna pattern for each
for az = azsweep
    for el=elsweep
        sv = steeringVector(f,[az; el]);
        taperedRIS.Taper = sltaper.*sv';
        
        %pattern(taperedRIS,f,azsweep,elsweep,'PropagationSpeed',c,'Type','powerdb', ...
          %  'CoordinateSystem','rectangular','Weights',sv);
    end
end
figure(), viewArray(taperedRIS,'Title','Tapered RIS','ShowTaper',true);
figure(), viewArray(RIS,'Title','NON Tapered RIS','ShowTaper',true);
%% 
%twinz = taylorwin(N_ele,nbar,sll);

% along the y axis
%twiny = taylorwin(N_ele,nbar,sll);

% Get the total taper values by multiplying the vectors of both dimensions
%tap = twinz*twiny';
tap = ones(N_ele, N_ele);
taperedRISt.Taper = ones.*PHI_comp';
taperedRISt_quant.Taper = ones.*exp(-1i*PHI');
% Apply the taper
%taperedRISt.Taper = tap.*(reshape((steeringVector_t(f,[-phi_RX;theta_RX])),[N_ele,N_ele]));
%taperedRISt.Taper = taperedRISt.Taper.*tap.*(reshape((steeringVector_t(f,theta_RX)),[N_ele,N_ele]));

% taperedRISr.Taper = tap.*(reshape((steeringVector_r(f,-phi_TX)),[N_ele,N_ele]));
% taperedRISr.Taper = taperedRISr.Taper.*tap.*(reshape((steeringVector_r(f,theta_TX)),[N_ele,N_ele]));

% figure(), viewArray(taperedRISt,'Title','Tapered RIS Tx','ShowTaper',true);
% figure(), viewArray(taperedRISr,'Title','Tapered RIS Rx','ShowTaper',true);
% figure(), viewArray(RIS,'Title','NON Tapered RIS','ShowTaper',true);
figure(), pattern(taperedRISt,f);
figure(), pattern(taperedRISt_quant,f);
figure(), pattern(RIS,f);
DirectivityRIS = max(max(pattern(taperedRISt,f)));
antenna = phased.IsotropicAntennaElement("FrequencyRange",[f-B f+B]);
%% RIS rectangular Array
%Con questo approccio non riesco a determinare il patter, OUT OF MEMORY

ant.Tilt = 0; ant.TiltAxis = [1 0 0];
RIS_ant = rectangularArray("Element",ant,"Size", [N_ele N_ele], ...
    "Tilt", 90,"TiltAxis", [0 1 0],...
    "RowSpacing",dy,"ColumnSpacing",dx);
RIS_ant_r = clone(RIS_ant);
RIS_ant_t = clone(RIS_ant_r);
RIS_ant_t.PhaseShift = (PHI_q(:))';
show(RIS_ant_r);
%% Modello della RIS considerando una RIS infinita

phi_obs = 0;
theta_obs = 20;
elevation_obs = 90 - theta_obs;
radius_obs = 100*lambda;
[x_obs,y_obs,z_obs] = sph2cart(phi_obs*pi/180,elevation_obs*pi/180,radius_obs);
obs_loc = [x_obs, y_obs, z_obs];
Nia = 50;
ifa=rectangularArray('Element',ant,'Size',[2*Nia+1 2*Nia+1],...
    'ColumnSpacing',ant.GroundPlaneLength,'RowSpacing',...
    ant.GroundPlaneWidth); %infinite RIS

irs = infiniteArray(Element=ant);
irs.ScanAzimuth = phi_obs;
irs.ScanElevation = elevation_obs;
numSummationTerms(irs,20);
figure(), show(irs); title("Infinite IRS");

AFdB=arrayFactor(ifa,f,irs.ScanAzimuth, irs.ScanElevation); % in dB
AF=10^(0.1*AFdB);

%% definizione delle antenne mediante ARRAYCONFIG
%Non posso specificare l'elemento radiante. Viene considerando un dipolo.
Tx_ant = arrayConfig("Size",[1 1]);
Rx_ant = arrayConfig("Size",[1 1]);
RIS_ant = arrayConfig("Size", [N_ele N_ele],"ElementSpacing", [dx dy]);
helperViewArray(Tx_ant); helperViewArray(Rx_ant); 
helperViewArray(RIS_ant); 

%% Modello con ConformalArray NON CONSIDERARLO
% The conformalArray class creates an antenna array using any element from the antenna
% or array library. You can also specify an array of any arbitrary geometry, such as a
% circular array, a nonplanar array, an array with nonuniform geometry,
% or a conformal array of arrays.Conformal arrays are used in:
% Direction-finding systems that use circular arrays or stacked circular arrays
% Aircraft systems due to surface irregularities or mechanical stress
RIS_v2 = rectangularArray("Element",ant,"Size", [N_ele N_ele], ...
    "RowSpacing",dy,"ColumnSpacing",dx, ...
    "Tilt",90);

Tx_RIS = conformalArray("Element",({ant RIS_v2}), ...
    "ElementPosition",[(Tx_loc{1})'; ([0 0 0])'], ...
    "AmplitudeTaper", [1 0 0 0 0 0]);

RIS_Rx = conformalArray("Element",({RIS_v2 ant}), ...
    "ElementPosition",[(Rx_loc{1})'; (rx.AntennaPosition)'], ...
    "AmplitudeTaper", [1 0]);

layout(Tx_RIS); layout(RIS_Rx);
show(Tx_RIS); show(RIS_Rx);
%figure(); pattern(RIS_v2,f,'patternOptions',p);

%Posso studiare il problema in due fasi. Prima simulando la propagazione
%dal Tx alla RIS e poi dalla RIS al Rx. 

%Tx-RIS Array Coupling Model
S_Tx = sparameters(Tx_RIS,f);
%RIS-Rx Array Coupling Model
S_Rx = sparameters(RIS_Rx,f);

%% phased Radiator and phased collector
%Usando queste due funzioni possiamo convertire una sequenza in un segnale
%elettromagnetico e possiamo catturare quest'ultimo per riconvertirlo in
%una sequenza. 

%COMBINERADIATEDSIGNALS enables the coherent summation of the radiated signals from
% all elements of an array to produce plane waves.
% Set this property to false to obtain individual radiated signal for
% each radiating element.
ULA = phased.ULA("Element",ant, ...
    "ElementSpacing",dx, ...
    "NumElements",N_ele);

%% !! per modificare la fase si usa il taper, ma non capisco come vengono settati!!
nbar=2; %number (nbar) of nearly constant-level sidelobes adjacent to the mainlobe.
sll = -20; %maximum sidelobe level in dB relative to the mainlobe peak.

% along the z axis
twinz = taylorwin(N_ele,nbar,sll);

% along the y axis
%twiny = taylorwin(N_ele,nbar,sll);

% Get the total taper values by multiplying the vectors of both dimensions
tap = twinz%*twiny.';

% Apply the taper
taperedULA = clone(ULA);
taperedULA.Taper = tap;

viewArray(taperedULA,'Title','Tapered URA','ShowTaper',true);

%%
URA = phased.ReplicatedSubarray("Subarray",ULA, ...
    "GridSpacing",dy, ...
    "SubarraySteering","Phase", ...
    "GridSize",[N_ele 1]);
viewArray(URA); pattern(URA,f);

%Definizione RIS trasmittente
RIS_tx = phased.Radiator("Sensor",URA, ...
    "OperatingFrequency",f, ...
    "CombineRadiatedSignals",true);

%Definizione trasmettitore
%antenna = phased.IsotropicAntennaElement("FrequencyRange",[f-B f+B]);
Tx = phased.Radiator("Sensor",antenna, ...
    "OperatingFrequency",f);

%Definizione ricevitore
Rx = phased.Radiator("Sensor",antenna, ...
    "OperatingFrequency",f);

%Definizione RIS ricevente
RIS_rx = phased.Collector("Sensor",URA, ...
    "OperatingFrequency",f, ...
    "SensorGainMeasure","dB");

% Definite le antenne, posso captare il segnale x ricevuto da una data
% direzione rappresentata dagli angoli
%NON SPECIFICO LA DISTANZA
y = RIS_rx(x,[Phi_Tx theta_TX]); 

%% definizione dei siti del trasmettitore e ricevitore

tx_power_real = 1;

tx = txsite("cartesian", ...
    "Antenna",antenna, ...
    "TransmitterPower", tx_power_real,...
    "AntennaPosition",(Tx_loc{1})', ...
    'TransmitterFrequency',f);

rx = rxsite("cartesian", ...
    "Antenna",antenna, ...
    "AntennaPosition",(Rx_loc{1})', ...
    "AntennaAngle",[0;90],"ReceiverSensitivity",-150);

%% Definizione del sito della RIS costituito da Tx e Rx sovrapposti.
% sv = steeringVector(f,[phi_Tx(1,1); theta_Tx(1,1)]);
% RIS.Taper = sltaper.*sv';

%taperedRISt
RIStx = txsite("cartesian", ...
    "Antenna",taperedRISt, ...
    "AntennaPosition",[0; floor(N_ele/2)*dx; -floor(N_ele/2)*dy], ...
    'TransmitterFrequency',f);

% la potenza trasmessa della RIS deve essere considerata quella ricevuta
% dal trasmettitore e opportunamente attenuata dalla RIS (1-3 dB)

%taperedRISr
RISrx = rxsite("cartesian", ...
    "Antenna",RIS, ...
    "AntennaPosition",[0; floor(N_ele/2)*dx; -floor(N_ele/2)*dy], ...
    "ReceiverSensitivity",-150);
%%
siteviewer("SceneModel",Stanza);
show(tx,'ShowAntennaHeight', false);
show(rx,'ShowAntennaHeight', false);
show(RIStx,'ShowAntennaHeight', false);
show(RISrx,'ShowAntennaHeight', false);
pattern(tx,f); pattern(RIStx,f); 

figure(); pattern(RIS,f,[-180:180],0, ...
    'PropagationSpeed',c, ...
    'CoordinateSystem','polar', ...
    'Type','powerdb','Normalize',true);
figure(); pattern(taperedRISt,f,[-180:180],0, ...
    'PropagationSpeed',c, ...
    'CoordinateSystem','polar', ...
    'Type','powerdb','Normalize',true);
figure(); pattern(taperedRISr,f,0,[-90:90], ...
    'PropagationSpeed',c, ...
    'CoordinateSystem','polar', ...
    'Type','powerdb','Normalize',true);
%% modello di propagazione (considerando 0 e 1 riflessione)

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

%% 
rays = raytrace(tx,rx,pm1);%"Map",siteviewer("SceneModel",Stanza));
rays = rays{1,1};
%Show the obstruction to the line-of-sight path.
%siteviewer("SceneModel",Stanza);
los(tx,rx);

raysRIS = raytrace(tx,RISrx,pm0,"Map",siteviewer("SceneModel",Stanza));

% plot(rays,'Colormap',jet,'ColorLimits',[50, 95]);
% tx_power_dBm-rays.PathLoss

% The comm.RayTracingChannel System object filters a signal through a multipath fading
% channel that is defined by propagation rays.
rtChan = comm.RayTracingChannel(rays,tx,RISrx);
rtChan.SampleRate = 300e6;
rtChan.ReceiverVirtualVelocity = [0.1; 0.1; 0];

%Use the showProfile object function to visualize the power delay profile (PDP),
% angle of departure (AoD) and angle of arrival (AoA) of the rays in the channel.
% In the visualization, the PDP has taken into account the transmit and receive
% array pattern gains in addition to the path loss for each ray.
% showProfile(rtChan);
% %QUELLO CHE VEDO NON è CORRETTO
% 
% numTx = rtChanInfo.NumTransmitElements;
% numRx = rtChanInfo.NumReceiveElements;

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
%siteviewer("SceneModel",Stanza);
%raytrace(tx,rx,pm1);
rays = raytrace(tx,rx,pm1,"Map",Stanza);
            
%% Calcolo potenza ricevuta considerando la RIS (calcolo errato perché mi considera due volte il FSPL)
Signal_Strength_TxRIS_dBm_err = sigstrength(RISrx,tx, ...
    'PropagationModel',pm0, ...
    'Map',Stanza);
%Attenuazione introdotta dalla RIS
%Reflection loss within bandwidth < 3 dB
Att_RISdB = -3*rand(1);
Att_RIS = db2mag(Att_RISdB);

Signal_Strength_TxRIS_err = 10^((Signal_Strength_TxRIS_dBm_err-30)/10);

Tx_RIS_Distance = distance(tx, RISrx,'Map',Stanza);
% los(tx,RISrx);
% los(RIStx,rx);
RISRx_Distance = distance(RIStx, rx,'Map',Stanza);

RIStx.TransmitterPower = Signal_Strength_TxRIS_err;
RIStx.SystemLoss = 3*rand(1);

Signal_Strength_RISRx_dBm_err = sigstrength(rx,RIStx, ...
    'PropagationModel',pm0, ...
    'Map',Stanza);
Signal_Strength_RISRx_err = 10^((Signal_Strength_RISRx_dBm_err-30)/10)*10^3;
%% Calcolo potenza ricevuta considerando la RIS (calcolo giusto perché considera solo una volta il FSPL)
Tx_RIS_Distance = distance(tx, RISrx,'Map',Stanza);
RISRx_Distance = distance(RIStx, rx,'Map',Stanza);
tx_power_Corr_dBm = tx_power_dBm_real + mag2db(Tx_RIS_Distance*RISRx_Distance/(Tx_RIS_Distance+RISRx_Distance))+mag2db(f)-147.55;
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
%raytrace(tx,RISrx,pm0);
%raytrace(RIStx,rx,pm0);
raysRISrx = raytrace(tx,RISrx,pm0,"Map",Stanza);
raysRIStx = raytrace(RIStx,rx,pm0,"Map",Stanza);

%% display risultati
disp(['Transmitted Power                    = ' , num2str(tx_power_dBm_real), ' dBm']);
disp(['Received Signal Strength Without RIS = ' , num2str(Signal_Strength_dBm_real), ' dBm']);
disp(['Received Signal Strength With RIS    = ' , num2str(Signal_Strength_RISRx_dBm_Corr), ' dBm']);
%% Visualizzazione dei fasci incidenti e riflessi

siteviewer("SceneModel",Stanza);

% raytrace(tx,rx,pm1);
% raytrace(tx,RISrx,pm0);
% raytrace(RIStx,rx,pm0);

%% Modellazione probabilistica del canale

%% |Tx-RIS Channel| 

LOS = rays{1,1}.LineOfSight;

% i parametri sono calcolati in corrispondenza di f = 24.2 GHz
% n = Path Loss Exponent
% sigma = Shadow Fading Term (dB)
% b = Path Loss Parameter
xTx = TX_loc(1) + RIS_px(1)*scale; yTx = TX_loc(2) + RIS_px(2)*scale; zTx = TX_loc(3) + RIS_px(3)*scale;
xRx = RX_loc(1) + RIS_px(1)*scale; yRx = RX_loc(2) + RIS_px(2)*scale; zRx = RX_loc(3) + RIS_px(3)*scale;

Tx_loc = [xTx yTx zTx];
Rx_loc = [xRx yRx zRx];

xRIS = RIS_px(1)*scale; yRIS = RIS_px(2)*scale; zRIS = RIS_px(3)*scale;
RIS_xyz = [xRIS yRIS zRIS];
% NLOS
n_NLOS=3.19;          % Path Loss Exponent (Indoor Office NLOS)
sigma_NLOS=8.29;      % Shadow Fading Term (dB) (Indoor Office NLOS)
b_NLOS=0.06;          % Path Loss Parameter (Indoor Office NLOS)
f0=24.2e9;            % Path Loss Parameter (GHz) (Indoor Office NLOS)

% LOS
n_LOS=1.73;           % Path Loss Exponent (Indoor Office LOS)
sigma_LOS=3.02;       % Shadow Fading Term (dB) (Indoor Office LOS)
b_LOS=0;              % Path Loss Parameter (Indoor Office NLOS)

% Element Radiation Pattern Parameters
q=0.285;
%Gain=pi;
Gain = maxG0;

%% Generazione cluster con distribuzione Poissoniana

lambda_p = 1.8; %Valore riferito a 28 GHz

% Number of Clusters
%C=max([1,poissrnd(lambda_p)]);  % Poisson distributed

% Number of Sub-rays per Cluster
%S=randi([1 30],1,C); % Uniform distributed

%% STEP 1
% Calculate Tx-RIS LOS distance  and Generate LOS Component for h 

d_T_RIS = Tx_RIS_Distance;
% LOS Probability is Relatively Low for Indoors if d_T_RIS > 20
if zRIS<zTx   % for ground level RIS
    % InH LOS Probability
    if d_T_RIS<= 1.2
        p_LOS=1;
    elseif 1.2<d_T_RIS && d_T_RIS<6.5
        p_LOS=exp(-(d_T_RIS-1.2)/4.7);
    else
        p_LOS=0.32*exp(-(d_T_RIS-6.5)/32.6);
    end

    I_LOS=randsrc(1,1,[1,0;p_LOS 1-p_LOS]);

elseif zRIS>=zTx % for an RIS mounted at a high place (100% LOS)
    I_LOS=1;
end

%% Calculate Tx Departure and RIS arrival angles to calculate array
% response vectors

% RIS arrival angles for LOS component
if I_LOS==1
    I_phi=sign(yTx-yRIS);
    phi_T_RIS_LOS = I_phi* atand ( abs( yRIS-yTx) / abs(xRIS-xTx) ); %Tra RIS e Tx
   
    I_theta=sign(zTx-zRIS);
    theta_T_RIS_LOS=I_theta * asind ( abs (zRIS-zTx ) / d_T_RIS );

    % Tx departure angles for LOS component
    I_phi_Tx=sign(yTx-yRIS);
    phi_Tx_LOS = I_phi_Tx* atand ( abs( yRIS-yTx) / abs(xRIS-xTx) ); %tra Tx e RIS
    % phi_Tx_LOS = I_phi_Tx* atand ( abs( xRIS-xTx) / abs(yRIS-yTx) ); %tra Tx e RIS
    
    I_theta_Tx=sign(zRIS-zTx);
    theta_Tx_LOS=I_theta_Tx * asind ( abs (zRIS-zTx )/ d_T_RIS );

    % Array Response Calculation (LOS)
    array_RIS_LOS=zeros(1,N_ele^2);

    counter2 = 1;
    counter3 = 1;
    for y=0:N_ele-1
        for z=0:N_ele-1
            % array_RIS_LOS(counter3)=exp(1i*k0*dx*(x*sind(theta_T_RIS_LOS) + y*sind(phi_T_RIS_LOS)*cosd(theta_T_RIS_LOS) )) ;
            array_RIS_LOS(counter3)=exp(1i*k0*dx*(y*sind(theta_T_RIS_LOS) + z*sind(phi_T_RIS_LOS)*cosd(theta_T_RIS_LOS) )) ;
            counter3=counter3+1;
        end
    end

    array_Tx_LOS = zeros(1,numTx);

    for y=0:sqrt(numTx)-1
        for z=0:sqrt(numTx)-1
            % array_Tx_LOS(counter2)=exp(1i*k0*dx*(x*sind(phi_Tx_LOS)*cosd(theta_Tx_LOS) + y*sind(theta_Tx_LOS))) ;
            array_Tx_LOS(counter2)=exp(1i*k0*dx*(y*sind(phi_Tx_LOS)*cosd(theta_Tx_LOS) + z*sind(theta_Tx_LOS))) ;
            counter2=counter2+1;
        end
    end

    % Link Attentuation (LOS) - teniamo conto anche della potenza trasmessa
    L_dB_LOS=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_RIS)- randn*sigma_LOS + tx_power_dBm_real;
    L_LOS=10^(L_dB_LOS/10);

    L_dB_LOS_wp=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_RIS)- randn*sigma_LOS;
    L_LOS_wp=10^(L_dB_LOS_wp/10);

    h_LOS=sqrt(L_LOS)*transpose(array_RIS_LOS)*array_Tx_LOS*exp(1i*rand*2*pi)*sqrt(Gain*(cosd(theta_T_RIS_LOS))^(2*q));
    h_LOS_wp=sqrt(L_LOS_wp)*transpose(array_RIS_LOS)*array_Tx_LOS*exp(1i*rand*2*pi)*sqrt(Gain*(cosd(theta_T_RIS_LOS))^(2*q));
    h_LOS_wd=sqrt(L_LOS)*transpose(array_RIS_LOS)*array_Tx_LOS*exp(1i*rand*2*pi)*sqrt(DirectivityRIS);
    h_LOS_wpwd=sqrt(L_LOS_wp)*transpose(array_RIS_LOS)*array_Tx_LOS*exp(1i*rand*2*pi)*sqrt(DirectivityRIS);
   
else
    h_LOS=0;
    h_LOS_wp = 0;
    h_LOS_wd =0;
    h_LOS_wpwd =0;
end

%% STEP 2
% Generate Clusters/Sub-rays, Azimuth/Elevation Departure Angles and Cluster Distances

for generate=1:100  % To ensure that at least one scatterer exist

    % Number of Clusters
    C=max([1,poissrnd(lambda_p)]);  % Poisson distributed

    % Number of Sub-rays per Cluster
    S=randi(30,1,C); % Uniformly distributed

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
    dim = [xURC-xLLC yURC-yLLC h_max];

    Coordinates=zeros(C,3);          % for Clusters
    Coordinates2=zeros(sum(S),3);    % for Scatterers
    for counter=1:C
        loop = 1;
        Coordinates(counter,:)=[xTx + a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
            yTx - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
            zTx + a_c(counter)*sind(theta_av(counter))] ;
        % Coordinates(counter,:)=[yTx + a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
        %     -xTx - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
        %     zTx + a_c(counter)*sind(theta_av(counter))] ;
        %while Coordinates(counter,3)>dim(3) || Coordinates(counter,3)<0 ||  Coordinates(counter,2)>dim(2) ||  Coordinates(counter,2)<0  ||  Coordinates(counter,1)>dim(1) ||  Coordinates(counter,1)<0
        while Coordinates(counter,3)>dim(3) || Coordinates(counter,3)<0 ||  Coordinates(counter,2)>dim(2) ||  Coordinates(counter,2)<0  ||  Coordinates(counter,1)>dim(1) ||  Coordinates(counter,1)<0
            a_c(counter)=    0.8*a_c(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
            % Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
            % in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
            Coordinates(counter,:)=[xTx + a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
                yTx - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
                zTx + a_c(counter)*sind(theta_av(counter))] ;
             % Coordinates(counter,:)=[yTx + a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
             %    -xTx - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
             %    zTx + a_c(counter)*sind(theta_av(counter))] ;
        end
    end
   
    a_c_rep=[];
    for counter3=1:C
        a_c_rep=[a_c_rep,repmat(a_c(counter3),1,S(counter3))];
    end
    for counter2=1:sum(S)
        Coordinates2(counter2,:)=[xTx + a_c_rep(counter2)*cosd(theta_Tx(counter2))*cosd(phi_Tx(counter2)),...
            yTx - a_c_rep(counter2)*cosd(theta_Tx(counter2))*sind(phi_Tx(counter2)),...
            zTx + a_c_rep(counter2)*sind(theta_Tx(counter2))] ;
        % Coordinates2(counter2,:)=[yTx - a_c_rep(counter2)*cosd(theta_Tx(counter2))*cosd(phi_Tx(counter2)),...
        %     xTx - a_c_rep(counter2)*cosd(theta_Tx(counter2))*sind(phi_Tx(counter2)),...
        %     zTx + a_c_rep(counter2)*sind(theta_Tx(counter2))] ;
    end

    % Correction on Scatters
    % You may ignore the scatterers outside the walls for Indoors scenario
    ignore=[];

    for counter2=1:sum(S)
        %if Coordinates2(counter2,3)>dim(3) || Coordinates2(counter2,3)<0 ||  Coordinates2(counter2,2)>dim(2) ||  Coordinates2(counter2,2)<0  ||  Coordinates2(counter2,1)>dim(1) ||  Coordinates2(counter2,1)<0
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

    b_cs(counter2)=norm(RIS_xyz-Coordinates2(counter2,:));   % Distance between Scatterer and RIS
    d_cs(counter2)=a_c_rep(counter2)+b_cs(counter2);         % Total distance Tx-Scatterer-RIS

    I_phi=sign(Coordinates2(counter2,2)-yRIS);
    phi_cs_RIS(counter2)  = I_phi* atand ( abs( yRIS-Coordinates2(counter2,2)) / abs(xRIS-Coordinates2(counter2,1)) );
    % I_phi=sign(-Coordinates2(counter2,1)+xRIS);
    % phi_cs_RIS(counter2)  = I_phi* atand ( abs( -xRIS+Coordinates2(counter2,1)) / abs(yRIS-Coordinates2(counter2,2)) );

    I_theta=sign(Coordinates2(counter2,3)-zRIS);
    theta_cs_RIS(counter2)=I_theta * asind ( abs (zRIS-Coordinates2(counter2,3) )/ b_cs(counter2) );

    I_phi_Tx_cs=sign(yTx-Coordinates2(counter2,2));
    phi_Tx_cs(counter2) = I_phi_Tx_cs* atand ( abs( Coordinates2(counter2,2)-yTx) / abs(Coordinates2(counter2,1)-xTx) );
    %  I_phi_Tx_cs=sign(-xTx+Coordinates2(counter2,1));
    % phi_Tx_cs(counter2) = I_phi_Tx_cs* atand ( abs( Coordinates2(counter2,1)-xTx) / abs(Coordinates2(counter2,2)-yTx) );
    

    I_theta_Tx_cs=sign(Coordinates2(counter2,3)-zTx);
    theta_Tx_cs(counter2)=I_theta_Tx_cs * asind ( abs (Coordinates2(counter2,3)-zTx )/ a_c_rep(counter2) );
end

%% STEP 4
% Array Response Calculation
array_cs_RIS=zeros(sum(S),N_ele^2);
for counter2=indices
    counter3=1;
    for x=0:N_ele-1
        for y=0:N_ele-1
            % array_cs_RIS(counter2,counter3)=exp(1i*k0*dx*(x*sind(theta_cs_RIS(counter2)) + y*sind(phi_cs_RIS(counter2))*cosd(theta_cs_RIS(counter2)) )) ;
            array_cs_RIS(counter2,counter3)=exp(1i*k0*dx*(x*sind(theta_cs_RIS(counter2)) + y*sind(phi_cs_RIS(counter2))*cosd(theta_cs_RIS(counter2)) )) ;
            counter3=counter3+1;
        end
    end
end

array_Tx_cs=zeros(sum(S),numTx);

for counter2 = indices
    counter3=1;
    for x=0:sqrt(numTx)-1
        for y = 0:sqrt(numTx)-1
            % array_Tx_cs(counter2,counter3)=exp(1i*k0*dx*(x*sind(phi_Tx_cs(counter2))*cosd(theta_Tx_cs(counter2)) + y*sind(theta_Tx_cs(counter2)) )) ;
            array_Tx_cs(counter2,counter3)=exp(1i*k0*dx*(x*sind(phi_Tx_cs(counter2))*cosd(theta_Tx_cs(counter2)) + y*sind(theta_Tx_cs(counter2)) )) ;
            counter3=counter3+1;
        end
    end
end
    
%% STEP 5 
% Calculate Link Attenuation and Generate Tx-RIS Channel (h) using 5G Channel Model

h_NLOS=zeros(N_ele.^2,numTx);
h_NLOS_wd=zeros(N_ele.^2,numTx);
h_NLOS_wp=zeros(N_ele.^2,numTx);
h_NLOS_wpwd=zeros(N_ele.^2,numTx);

beta=zeros(1,sum(S)); % to be reused for shared clusters 
shadow=beta;          % to be reused for shared clusters
for counter2=indices
    X_sigma=randn*sigma_NLOS;
    Lcs_dB=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs(counter2))- X_sigma + tx_power_dBm_real;%-Att_RISdB;
    Lcs=10^(Lcs_dB/10);

    Lcs_dB_wp=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs(counter2))- X_sigma;%-Att_RISdB;
    Lcs_wp=10^(Lcs_dB_wp/10);

    beta(counter2)=((randn+1i*randn)./sqrt(2));  % common complex gain for shared clusters
    shadow(counter2)=X_sigma;                    % commun shadow factor for shared clusters
    
    h_NLOS = h_NLOS + beta(counter2)*sqrt(Gain*(cosd(theta_cs_RIS(counter2)))^(2*q))*sqrt(Lcs)*transpose(array_cs_RIS(counter2,:))*array_Tx_cs(counter2,:);  % consider all scatters
    h_NLOS_wp = h_NLOS_wp + beta(counter2)*sqrt(Gain*(cosd(theta_cs_RIS(counter2)))^(2*q))*sqrt(Lcs_wp)*transpose(array_cs_RIS(counter2,:))*array_Tx_cs(counter2,:);  % consider all scatters
    h_NLOS_wd = h_NLOS_wd + beta(counter2)*sqrt(DirectivityRIS*(cosd(theta_cs_RIS(counter2)))^(2*q))*sqrt(Lcs)*transpose(array_cs_RIS(counter2,:))*array_Tx_cs(counter2,:);
    h_NLOS_wpwd = h_NLOS_wpwd + beta(counter2)*sqrt(DirectivityRIS*(cosd(theta_cs_RIS(counter2)))^(2*q))*sqrt(Lcs_wp)*transpose(array_cs_RIS(counter2,:))*array_Tx_cs(counter2,:);

end
h_NLOS=h_NLOS.*sqrt(1/M_new);  % normalization
h=h_NLOS+h_LOS; % include the LOS component (if any) h_LOS=0 when there is no LOS

h_NLOS_wd=h_NLOS_wd.*sqrt(1/M_new);  % normalization
h_wd=h_NLOS_wd+h_LOS_wd; % include the LOS component (if any) h_LOS=0 when there is no LOS


h_NLOS_wp=h_NLOS_wp.*sqrt(1/M_new);  % normalization
h_wp=h_NLOS_wp+h_LOS_wp; % include the LOS component (if any) h_LOS=0 when there is no LOS

h_NLOS_wpwd=h_NLOS_wpwd.*sqrt(1/M_new);  % normalization
h_wpwd=h_NLOS_wpwd+h_LOS_wpwd;

%% STEPS 6-7 
% Generation of g (RIS-Rx Channel) - GENERATE A LOS CHANNEL

% Calculate Departure Angles Considering RIS and Rx Coordinates
d_RIS_R=norm(RIS_xyz-Rx_loc);

% Elevation Departure Angle
I_theta=sign(zRx - zRIS);
theta_Rx_RIS=I_theta * asind( abs(zRx-zRIS)/d_RIS_R ); % AoD of RIS

% Azimuth Departure Angle
I_phi=sign(yRx - yRIS);
phi_Rx_RIS=I_phi * atand( abs(yRx-yRIS)/ abs(xRx-xRIS) );
% I_phi=sign(-xRx + xRIS);
% phi_Rx_RIS=I_phi * atand( abs(-xRx+RIS)/ abs(yRx-yRIS) );

% AoA angles of Rx for g_LOS channel in an Indoor
phi_av_Rx  = rand*180-90;     % mean azimuth
theta_av_Rx= rand*180-90;      % mean elevation

phi_Rx   = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + phi_av_Rx];
theta_Rx = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + theta_av_Rx];

% Recalculate Array Response for Two Angles (in Rx direction)
array_2=zeros(1,N_ele^2);
counter3=1;
for x=0:N_ele-1
    for y=0:N_ele-1
        % array_2(counter3)=exp(1i*k0*dx*(x*sind(theta_Rx_RIS) + y*sind(phi_Rx_RIS)*cosd(theta_Rx_RIS) )) ;
        array_2(counter3)=exp(1i*k0*dx*(x*sind(theta_Rx_RIS) + y*sind(phi_Rx_RIS)*cosd(theta_Rx_RIS) )) ;
        counter3=counter3+1;
    end
end

array_Rx = zeros(1,numRx);

counter3=1;
for x=0:sqrt(numRx)-1
    for y=0:sqrt(numRx)-1
        % array_Rx(counter3)=exp(1i*k0*dx*(x*sind(phi_Rx)*cosd(theta_Rx)+y*sind(theta_Rx)  )) ;
        array_Rx(counter3)=exp(1i*k0*dx*(x*sind(phi_Rx)*cosd(theta_Rx)+y*sind(theta_Rx)  )) ;
        counter3=counter3+1;
    end
end

% % LOS Link Attenuation
L_dB_LOS_2=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_RIS_R)- randn*sigma_LOS + tx_power_dBm_real;
L_LOS_2=10^(L_dB_LOS_2/10);

L_dB_LOS_2_wp=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_RIS_R)- randn*sigma_LOS ;
L_LOS_2_wp=10^(L_dB_LOS_2/10);

% Generate g (Pure LOS)
g=sqrt(Gain*(cosd(theta_Rx_RIS))^(2*q))*sqrt(L_LOS_2)*transpose(array_2)*array_Rx*exp(1i*rand*2*pi);
g_wd=sqrt(DirectivityRIS)*sqrt(L_LOS_2)*transpose(array_2)*array_Rx*exp(1i*rand*2*pi);
g_wp = sqrt(Gain*(cosd(theta_Rx_RIS))^(2*q))*sqrt(L_LOS_2_wp)*transpose(array_2)*array_Rx*exp(1i*rand*2*pi);
g_wpwd = sqrt(DirectivityRIS)*sqrt(L_LOS_2_wp)*transpose(array_2)*array_Rx*exp(1i*rand*2*pi);

%% STEP 8
% Generation of h_SISO

d_T_R=norm(Tx_loc-Rx_loc);

d_cs_tilde=zeros(1,sum(S));
h_SISO_NLOS=0;
h_SISO_NLOS_wp=0;

for counter2=indices

    % due to shared clusters d_cs_tilde ~ d_cs
    d_cs_tilde(counter2) = a_c_rep(counter2) + norm(Coordinates2(counter2,:)- Rx_loc);


    I_phi_Tx_cs_SISO=sign(yTx-Coordinates2(counter2,2));
    phi_Tx_cs_SISO(counter2) = I_phi_Tx_cs_SISO* atand ( abs( Coordinates2(counter2,2)-yTx) / abs(Coordinates2(counter2,1)-xTx) );
    % I_phi_Tx_cs_SISO=sign(-xTx+Coordinates2(counter2,1));
    % phi_Tx_cs_SISO(counter2) = I_phi_Tx_cs_SISO* atand ( abs( -Coordinates2(counter2,1)+xTx) / abs(Coordinates2(counter2,2)-yTx) );

    I_theta_Tx_cs_SISO=sign(Coordinates2(counter2,3)-zTx);
    theta_Tx_cs_SISO(counter2)=I_theta_Tx_cs_SISO * asind ( abs (Coordinates2(counter2,3)-zTx )/ a_c_rep(counter2) );

    % AoA for Rx in an Indoor
    phi_av_SISO(counter2)  = rand*180-90;     % mean azimuth
    theta_av_SISO(counter2)= rand*180-90;      % mean elevation

    phi_cs_Rx_SISO(counter2)   = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + phi_av_SISO(counter2)];
    theta_cs_Rx_SISO(counter2) = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + theta_av_SISO(counter2)];

    counter3=1;
    for x=0:sqrt(numRx)-1
        for y=0:sqrt(numRx)-1
            % array_Rx_cs_SISO(counter2,counter3)=exp(1i*k0*dx*(x*sind(phi_cs_Rx_SISO(counter2))*cosd(theta_cs_Rx_SISO(counter2))+y*sind(theta_cs_Rx_SISO(counter2)) )) ;
            array_Rx_cs_SISO(counter2,counter3)=exp(1i*k0*dx*(x*sind(phi_cs_Rx_SISO(counter2))*cosd(theta_cs_Rx_SISO(counter2))+y*sind(theta_cs_Rx_SISO(counter2)) )) ;
            counter3=counter3+1;
        end
    end

    counter3=1;
    for x=0:sqrt(numTx)-1
        for y=0:sqrt(numTx)-1
            % array_Tx_cs_SISO(counter2,counter3)=exp(1i*k0*dx*(x*sind(phi_Tx_cs_SISO(counter2))*cosd(theta_Tx_cs_SISO(counter2))+y*sind(theta_Tx_cs_SISO(counter2))  )) ;
            array_Tx_cs_SISO(counter2,counter3)=exp(1i*k0*dx*(x*sind(phi_Tx_cs_SISO(counter2))*cosd(theta_Tx_cs_SISO(counter2))+y*sind(theta_Tx_cs_SISO(counter2))  )) ;
            counter3=counter3+1;
        end
    end
    
    Lcs_dB_SISO=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs_tilde(counter2))- shadow(counter2) + tx_power_dBm_real;
    Lcs_SISO=10^(Lcs_dB_SISO/10);

    Lcs_dB_SISO_wp=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs_tilde(counter2))- shadow(counter2);
    Lcs_SISO_wp=10^(Lcs_dB_SISO_wp/10);

    % We consider the same complex path gain (small scale fading) with an excess phase (similar to array response)
    eta=k0* ( norm(Coordinates2(counter2,:)- RIS_xyz) -  norm(Coordinates2(counter2,:)- Rx_loc));

    h_SISO_NLOS = h_SISO_NLOS + beta(counter2)*exp(1i*eta)*sqrt(Lcs_SISO)*transpose(array_Rx_cs_SISO(counter2,:))*array_Tx_cs_SISO(counter2,:);
    h_SISO_NLOS_wp = h_SISO_NLOS_wp + beta(counter2)*exp(1i*eta)*sqrt(Lcs_SISO_wp)*transpose(array_Rx_cs_SISO(counter2,:))*array_Tx_cs_SISO(counter2,:);
end

h_SISO_NLOS=h_SISO_NLOS.*sqrt(1/M_new);  % normalization
h_SISO_NLOS_wp=h_SISO_NLOS_wp.*sqrt(1/M_new);

if zRIS >= zTx

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
elseif zRIS < zTx  % RIS in the ground level
    I_LOS_3=I_LOS;
end

if I_LOS_3==1
    L_SISO_LOS_dB=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_R)- randn*sigma_LOS + tx_power_dBm_real;
    L_SISO_LOS=10^(L_SISO_LOS_dB/10);

    L_SISO_LOS_dB_wp=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_R)- randn*sigma_LOS;
    L_SISO_LOS_wp=10^(L_SISO_LOS_dB_wp/10);

    I_phi_Tx_SISO=sign(yTx-yRx);
    phi_Tx_SISO = I_phi_Tx_SISO* atand ( abs( yTx-yRx) / abs(xTx-xRx) );
    % I_phi_Tx_SISO=sign(-xTx+xRx);
    % phi_Tx_SISO = I_phi_Tx_SISO* atand ( abs( -xTx+xRx) / abs(yTx-yRx) );

    I_theta_Tx_SISO=sign(zRx-zTx);
    theta_Tx_SISO= I_theta_Tx_SISO* atand ( abs( zRx-zTx) / abs(d_T_R) );


    % AoA of Rx for Tx-Rx channel in an Indoor
    phi_av_SISO_LOS = rand*180-90;     % mean azimuth
    theta_av_SISO_LOS= rand*180-90;      % mean elevation

    phi_Rx_SISO  = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + phi_av_SISO_LOS];
    theta_Rx_SISO= [log(rand(1,1)./rand(1,1))*sqrt(25/2) + theta_av_SISO_LOS];


    counter3=1;
    for x=0:sqrt(numTx)-1
        for y=0:sqrt(numTx)-1
            % array_Tx_SISO(counter3)=exp(1i*k0*dx*(x*sind(phi_Tx_SISO)*cosd(theta_Tx_SISO)+y*sind(theta_Tx_SISO)  )) ;
            array_Tx_SISO(counter3)=exp(1i*k0*dx*(x*sind(phi_Tx_SISO)*cosd(theta_Tx_SISO)+y*sind(theta_Tx_SISO)  )) ;
            counter3=counter3+1;
        end
    end
    
    counter3=1;
    for x=0:sqrt(numRx)-1
        for y=0:sqrt(numRx)-1
            % array_Rx_SISO(counter3)=exp(1i*k0*dx*(x*sind(phi_Rx_SISO)*cosd(theta_Rx_SISO) +y*sind(theta_Rx_SISO))) ;
            array_Rx_SISO(counter3)=exp(1i*k0*dx*(x*sind(phi_Rx_SISO)*cosd(theta_Rx_SISO) +y*sind(theta_Rx_SISO))) ;
            counter3=counter3+1;
        end
    end
    
    h_SISO_LOS= sqrt(L_SISO_LOS)*exp(1i*rand*2*pi)*transpose(array_Rx_SISO)*array_Tx_SISO;
    h_SISO_LOS_wp= sqrt(L_SISO_LOS_wp)*exp(1i*rand*2*pi)*transpose(array_Rx_SISO)*array_Tx_SISO;
else
    h_SISO_LOS=0;
    h_SISO_LOS_wp=0;
end

%% consideriamo nullo il LOS
h_SISO_LOS = 0; h_SISO_LOS_wp = 0;
h_SISO=h_SISO_NLOS + h_SISO_LOS; 
h_SISO_wp=h_SISO_NLOS_wp + h_SISO_LOS_wp; 

%% total channel

%   g     = Vector of LOS channel coefficients between the RIS and the Rx
%  THETA  = Matrix of RIS element responses
%   h     = Vector of channel coefficients for the Tx-RIS link composed of M scatters
% h_SISO  = Characterizes the direct link channel between Tx and Rx
%   x     = Transmitted signal
% PHI_vet =[];
PHI1 = PHI';
% 
% for k=1:N_ele
%     PHI_vet = [PHI_vet PHI1(k,:)];
% end
% 
% THETA = zeros(N_ele^2,N_ele^2);
% 
% for k = 1:N_ele^2
%     THETA(k,k) = exp(i*PHI_vet(N_ele^2-k+1));
% end
% % PHI = rad2deg(PHI_comp(:));
% % for k = 1:N_ele^2
% %     THETA(k,k) = exp(1i*PHI(k));
% % end
% %THETA = Att_RIS.*THETA;
% %Tot_Channel = reshape(g,N_ele,N_ele)*THETA*reshape(h,N_ele,N_ele) + h_SISO;
% Tot_Channel = g' * THETA * h + h_SISO;
% Tot_Channel_wd = g_wd' * THETA * h_wd + h_SISO;
% Tot_Channel_wp = g_wp' * THETA * h_wp + h_SISO_wp;
% Tot_Channel_w1p = g_wp' * THETA * h + h_SISO;

PHI_vet_trasp =[];

for k=1:N_ele
    PHI_vet_trasp = [PHI_vet_trasp PHI1(:,k)'];
end

THETA_trasp = zeros(N_ele^2,N_ele^2);

for k = 1:N_ele^2
    THETA_trasp(k,k) = exp(i*PHI_vet_trasp(N_ele^2-k+1));
end

%PHI_trasp = rad2deg(PHI_comp(:))';
% THETA_trasp = zeros(N_ele^2,N_ele^2); 
% for k = 1:N_ele^2
%     THETA_trasp(k,k) = exp(1i*PHI_trasp(k));
% end

Tot_Channel_trasp = g' * THETA_trasp * h + h_SISO;
Tot_Channel_wd_trasp = g_wd' * THETA_trasp * h_wd + h_SISO;
Tot_Channel_wp_trasp = g_wp' * THETA_trasp * h_wp + h_SISO_wp;
Tot_Channel_w1p_trasp = g_wp' * THETA_trasp * h + h_SISO;
%% tot channel version 2

PHI_vet2 =[];

for k=1:N_ele
    PHI_vet2 = [PHI_vet2 PHI1(k,:)];
end

THETA2 = zeros(N_ele^2,N_ele^2);

for k = 1:N_ele^2
    THETA2(k,k) = exp(i*PHI_vet2(k));
end

Tot_Channel2 = g' * THETA2 * h + h_SISO;
Tot_Channel_wd2 = g_wd' * THETA2 * h_wd + h_SISO;
Tot_Channel_wp2 = g_wp' * THETA2 * h_wp + h_SISO_wp;
Tot_Channel_w1p2 = g_wp' * THETA2 * h + h_SISO;

PHI_vet_trasp2 =[];

for k=1:N_ele
    PHI_vet_trasp2 = [PHI_vet_trasp2 PHI1(:,k)'];
end

THETA_trasp2 = zeros(N_ele^2,N_ele^2);

for k = 1:N_ele^2
    THETA_trasp2(k,k) = exp(i*PHI_vet_trasp2(k));
end

Tot_Channel_trasp2 = g' * THETA_trasp2 * h + h_SISO;
Tot_Channel_wd_trasp2 = g_wd' * THETA_trasp2 * h_wd + h_SISO;
Tot_Channel_wp_trasp2 = g_wp' * THETA_trasp2 * h_wp + h_SISO_wp;
Tot_Channel_w1p_trasp2 = g_wp' * THETA_trasp2 * h + h_SISO;

%% Ergodic achievable rate of a communication system
Noise_power_dBm = -100; %dBm
p = abs(Tot_Channel_wp_trasp)^2 * db2mag(tx_power_dBm_real) / db2mag(Noise_power_dBm); %Signal Nooise Ratio
p_dB = mag2db(p);
R = mean(log2(1+p)); %[bits/s/Hz]

%% Modello canale senza RIS
%AGGIUNGERE ATTENUAZIONE DEL MURO

%d_T_R_WRIS=rays{1,1}.PropagationDistance;
Rx_WRIS = Rx_loc{1}; Rx_WRIS(1) = -Rx_WRIS(1);
d_T_R_WRIS=norm(Tx_loc{1}-Rx_WRIS);

d_cs_tilde_WRIS=zeros(1,sum(S));
h_SISO_NLOS_WRIS=0;

for counter2=indices

    % due to shared clusters d_cs_tilde ~ d_cs
    d_cs_tilde_WRIS(counter2) = a_c_rep(counter2) + norm(Coordinates2(counter2,:)- Rx_WRIS);


    I_phi_Tx_cs_SISO_WRIS=sign(yTx-Coordinates2(counter2,2));
    phi_Tx_cs_SISO_WRIS(counter2) = I_phi_Tx_cs_SISO_WRIS* atand ( abs( Coordinates2(counter2,2)-yTx) / abs(Coordinates2(counter2,1)-xTx) );

    I_theta_Tx_cs_SISO_WRIS=sign(Coordinates2(counter2,3)-zTx);
    theta_Tx_cs_SISO_WRIS(counter2)=I_theta_Tx_cs_SISO_WRIS * asind ( abs (Coordinates2(counter2,3)-zTx )/ a_c_rep(counter2) );

    % AoA for Rx in an Indoor
    phi_av_SISO_WRIS(counter2)  = rand*180-90;     % mean azimuth
    theta_av_SISO_WRIS(counter2)= rand*180-90;      % mean elevation

    phi_cs_Rx_SISO_WRIS(counter2)   = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + phi_av_SISO_WRIS(counter2)];
    theta_cs_Rx_SISO_WRIS(counter2) = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + theta_av_SISO_WRIS(counter2)];

    counter3=1;
    for x=0:sqrt(numRx)-1
        for y=0:sqrt(numRx)-1
            array_Rx_cs_SISO_WRIS(counter2,counter3)=exp(1i*k0*dx*(x*sind(phi_cs_Rx_SISO_WRIS(counter2))*cosd(theta_cs_Rx_SISO_WRIS(counter2))+y*sind(theta_cs_Rx_SISO_WRIS(counter2)) )) ;
            counter3=counter3+1;
        end
    end

    counter3=1;
    for x=0:sqrt(numTx)-1
        for y=0:sqrt(numTx)-1
            array_Tx_cs_SISO_WRIS(counter2,counter3)=exp(1i*k0*dx*(x*sind(phi_Tx_cs_SISO_WRIS(counter2))*cosd(theta_Tx_cs_SISO_WRIS(counter2))+y*sind(theta_Tx_cs_SISO_WRIS(counter2))  )) ;
            counter3=counter3+1;
        end
    end

    Lcs_dB_SISO_WRIS=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs_tilde_WRIS(counter2))- shadow(counter2) + tx_power_dBm_real;
    Lcs_SISO_WRIS=10^(Lcs_dB_SISO_WRIS/10);

    % We consider the same complex path gain (small scale fading) with an excess phase (similar to array response)
    eta_WRIS=k0* ( norm(Coordinates2(counter2,:)- RIS_xyz) -  norm(Coordinates2(counter2,:)- Rx_loc{1}));

    h_SISO_NLOS_WRIS = h_SISO_NLOS_WRIS + beta(counter2)*exp(1i*eta_WRIS)*sqrt(Lcs_SISO_WRIS)*transpose(array_Rx_cs_SISO_WRIS(counter2,:))*array_Tx_cs_SISO_WRIS(counter2,:);
end

h_SISO_NLOS_WRIS=h_SISO_NLOS_WRIS.*sqrt(1/M_new);  % normalization

if I_LOS_3==1
    L_SISO_LOS_dB_WRIS=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_R_WRIS)- randn*sigma_LOS + tx_power_dBm_real;
    L_SISO_LOS_WRIS=10^(L_SISO_LOS_dB_WRIS/10);

    I_phi_Tx_SISO_WRIS=sign(yTx-Rx_WRIS(2));
    phi_Tx_SISO_WRIS = I_phi_Tx_SISO_WRIS* atand ( abs( yTx-Rx_WRIS(2)) / abs(xTx-Rx_WRIS(1)) );

    I_theta_Tx_SISO_WRIS=sign(zRx-zTx);
    theta_Tx_SISO_WRIS= I_theta_Tx_SISO_WRIS* atand ( abs( zRx-zTx) / abs(d_T_R_WRIS) );


    % AoA of Rx for Tx-Rx channel in an Indoor
    phi_av_SISO_LOS_WRIS = rand*180-90;     % mean azimuth
    theta_av_SISO_LOS_WRIS= rand*180-90;      % mean elevation

    phi_Rx_SISO_WRIS  = [log(rand(1,1)./rand(1,1))*sqrt(25/2) + phi_av_SISO_LOS_WRIS];
    theta_Rx_SISO_WRIS= [log(rand(1,1)./rand(1,1))*sqrt(25/2) + theta_av_SISO_LOS_WRIS];


    counter3=1;
    for x=0:sqrt(numTx)-1
        for y=0:sqrt(numTx)-1
            array_Tx_SISO_WRIS(counter3)=exp(1i*k0*dx*(x*sind(phi_Tx_SISO_WRIS)*cosd(theta_Tx_SISO_WRIS)+y*sind(theta_Tx_SISO_WRIS)  )) ;
            counter3=counter3+1;
        end
    end

    counter3=1;
    for x=0:sqrt(numRx)-1
        for y=0:sqrt(numRx)-1
            array_Rx_SISO_WRIS(counter3)=exp(1i*k0*dx*(x*sind(phi_Rx_SISO_WRIS)*cosd(theta_Rx_SISO_WRIS) +y*sind(theta_Rx_SISO_WRIS))) ;
            counter3=counter3+1;
        end
    end



    h_SISO_LOS_WRIS= sqrt(L_SISO_LOS_WRIS)*exp(1i*rand*2*pi)*transpose(array_Rx_SISO_WRIS)*array_Tx_SISO_WRIS;
else
    h_SISO_LOS_WRIS=0;
end

h_SISO_WRIS=h_SISO_NLOS_WRIS + h_SISO_LOS_WRIS; 

%%
% Valutiamo se la controparte nella comunicazione riceve o meno il segnale
% Otteniamo dei collegamenti verdi se il dispositivo è coperto dal
% diagramma di radiazione, altrimenti otteniamo un collegamento rosso. 
sc = [0 0.3 0];
link(rx,tx,pm1,"SuccessColor",sc);
link(RISrx,tx,pm0,"SuccessColor",sc);
link(rx,RIStx,pm0,"SuccessColor",sc);

%% Modello con planeWaveExcitation
%The planeWaveExcitation object creates an environment in which a plane wave
% excites an antenna or array. Plane wave excitation is a scattering solution
% that solves the receiver antenna problem.

for mm = 1:numel(theta_inc)

    % Construct a planeWaveExcitation object with the IRS as element
    dir=[sind(theta_inc(mm))*cosd(phi_inc);...
        sind(theta_inc(mm))*sind(phi_inc);...
        -cosd(theta_inc(mm))];
    pol=[cosd(theta_inc(mm))*cosd(phi_inc);...
        cosd(theta_inc(mm))*sind(phi_inc);...
        sind(theta_inc(mm))];
    pw=planeWaveExcitation('Element',irs);
    pw.Direction=dir;
    pw.Polarization=pol;

    % Compute incident field upon the IRS
    Ein = pol;
    Einc = dot(Ein,pol/norm(pol));

    % Compute the outward scattered field from the IRS
    [Eo,~] = EHfields(pw,f,rx.AntennaPosition);
    Eo = Eo*AF;
    Eobs = dot(Eo,pol/norm(pol));

    % Compute the reflection coefficient of the IRS
    MagReflection(1,mm) = abs(Eobs/Einc);
    PhaseReflection(1,mm) = (angle(Eobs/Einc))*180/pi;
end
%% Trasmissione e ricezione del segnale modulato 64QAM/OFDM

scs = 30;
carrier = nrCarrierConfig('SubcarrierSpacing',scs);
p = 1;
pdsch = nrPDSCHConfig('NumLayers',p,'Modulation','64QAM');
[ind,info] = nrPDSCHIndices(carrier,pdsch);
numDataBits = info.G;
cws = randi([0 1],numDataBits,1);
sym = nrPDSCH(carrier,pdsch,cws,'OutputDataType','single');
txGrid = nrResourceGrid(carrier,p);
txGrid(ind) = sym;
initialNSlot = carrier.NSlot;
cpl = 'extended';
[txWaveform,Info] = nrOFDMModulate(txGrid,scs,initialNSlot,'CyclicPrefix',cpl);


% waveform = conv2(txWaveform,Tot_Channel,'same');
% waveform_wd = conv2(txWaveform,Tot_Channel_wd,'same'); %circular convolution
% waveform_wp = conv2(txWaveform,Tot_Channel_wp,'same'); %circular convolution
% waveform_w1p = conv2(txWaveform,Tot_Channel_w1p,'same'); %circular convolution
waveform_trasp = conv2(txWaveform,Tot_Channel_trasp,'same'); %circular convolution
waveform_wd_trasp = conv2(txWaveform,Tot_Channel_wd_trasp,'same'); %circular convolution
waveform_wp_trasp = conv2(txWaveform,Tot_Channel_wp_trasp,'same'); %circular convolution
waveform_w1p_trasp = conv2(txWaveform,Tot_Channel_w1p_trasp,'same'); %circular convolution
% waveform2 = conv2(txWaveform,Tot_Channel2,'same'); %circular convolution
% waveform_wd2 = conv2(txWaveform,Tot_Channel_wd2,'same'); %circular convolution
% waveform_wp2 = conv2(txWaveform,Tot_Channel_wp2,'same'); %circular convolution
% waveform_w1p2 = conv2(txWaveform,Tot_Channel_w1p2,'same'); %circular convolution
% waveform_trasp2 = conv2(txWaveform,Tot_Channel_trasp2,'same'); %circular convolution
% waveform_wd_trasp2 = conv2(txWaveform,Tot_Channel_wd_trasp2,'same'); %circular convolution
% waveform_wp_trasp2 = conv2(txWaveform,Tot_Channel_wp_trasp2,'same'); %circular convolution
% waveform_w1p_trasp2 = conv2(txWaveform,Tot_Channel_w1p_trasp2,'same'); %circular convolution

nrb = carrier.NSizeGrid;
bitTx = qamdemod(txGrid,64,'OutputType','bit');

% grid = nrOFDMDemodulate(waveform,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
% grid_wd = nrOFDMDemodulate(waveform_wd,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
% grid_wp = nrOFDMDemodulate(waveform_wp,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
% grid_w1p = nrOFDMDemodulate(waveform_w1p,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_trasp = nrOFDMDemodulate(waveform_trasp,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_wd_trasp = nrOFDMDemodulate(waveform_wd_trasp,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_wp_trasp = nrOFDMDemodulate(waveform_wp_trasp,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
grid_w1p_trasp = nrOFDMDemodulate(waveform_w1p_trasp,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
% grid2 = nrOFDMDemodulate(waveform2,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
% grid_wd2 = nrOFDMDemodulate(waveform_wd2,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
% grid_wp2 = nrOFDMDemodulate(waveform_wp2,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
% grid_w1p2 = nrOFDMDemodulate(waveform_w1p2,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
% grid_trasp2 = nrOFDMDemodulate(waveform_trasp2,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
% grid_wd_trasp2 = nrOFDMDemodulate(waveform_wd_trasp2,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
% grid_wp_trasp2 = nrOFDMDemodulate(waveform_wp_trasp2,nrb,scs,initialNSlot,'CyclicPrefix',cpl);
% grid_w1p_trasp2 = nrOFDMDemodulate(waveform_w1p_trasp2,nrb,scs,initialNSlot,'CyclicPrefix',cpl);

% bitRx = qamdemod(grid,64,'OutputType','bit');
% bitRx_wd = qamdemod(grid_wd,64,'OutputType','bit');
% bitRx_wp = qamdemod(grid_wp,64,'OutputType','bit');
% bitRx_w1p = qamdemod(grid_w1p,64,'OutputType','bit');
bitRx_trasp = qamdemod(grid_trasp,64,'OutputType','bit');
bitRx_wd_trasp = qamdemod(grid_wd_trasp,64,'OutputType','bit');
bitRx_wp_trasp = qamdemod(grid_wp_trasp,64,'OutputType','bit');
bitRx_w1p_trasp = qamdemod(grid_w1p_trasp,64,'OutputType','bit');
% bitRx2 = qamdemod(grid2,64,'OutputType','bit');
% bitRx_wd2 = qamdemod(grid_wd2,64,'OutputType','bit');
% bitRx_wp2 = qamdemod(grid_wp2,64,'OutputType','bit');
% bitRx_w1p2 = qamdemod(grid_w1p2,64,'OutputType','bit');
% bitRx_trasp2 = qamdemod(grid_trasp2,64,'OutputType','bit');
% bitRx_wd_trasp2 = qamdemod(grid_wd_trasp2,64,'OutputType','bit');
% bitRx_wp_trasp2 = qamdemod(grid_wp_trasp2,64,'OutputType','bit');
% bitRx_w1p_trasp2 = qamdemod(grid_w1p_trasp2,64,'OutputType','bit');

% difference between Tx and Rx bit
% [number,ratio] = biterr(bitTx,bitRx);
% BitErrorRate = ratio*100; %percentage of error bits
[number_trasp,ratio_trasp] = biterr(bitTx,bitRx_trasp);
BitErrorRate_trasp = ratio_trasp*100; %percentage of error bits

disp(['Modulation OFDM with 64 QAM with RIS ']);
% disp(['Bit Error Rate Considering TxPower between Tx-RIS and RIS-Rx = ' , num2str(BitErrorRate), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_trasp), ' %']);
disp(['Bit Error Rate Considering TxPower between Tx-RIS and RIS-Rx = ' , num2str(BitErrorRate_trasp), ' %']);

% [number_wd,ratio_wd] = biterr(bitTx,bitRx_wd);
% BitErrorRate_wd = ratio_wd*100; %percentage of error bits
[number_wd_trasp,ratio_wd_trasp] = biterr(bitTx,bitRx_wd_trasp);
BitErrorRate_wd_trasp = ratio_wd_trasp*100; %percentage of error bits
% disp(['Bit Error Rate Considering Directivity and TxPower between Tx-RIS and RIS-Rx = ' , num2str(BitErrorRate_wd), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_wd_trasp), ' %']);
disp(['Bit Error Rate Considering Directivity and TxPower between Tx-RIS and RIS-Rx = ' , num2str(BitErrorRate_wd_trasp), ' %']);

% [number_wp,ratio_wp] = biterr(bitTx,bitRx_wp);
% BitErrorRate_wp = ratio_wp*100; %percentage of error bits
[number_wp_trasp,ratio_wp_trasp] = biterr(bitTx,bitRx_wp_trasp);
BitErrorRate_wp_trasp = ratio_wp_trasp*100; %percentage of error bits
% disp(['Bit Error Rate Without TxPower and considering the Gain = ' , num2str(BitErrorRate_wp), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_wp_trasp), ' %']);
disp(['Bit Error Rate Without TxPower and considering the Gain = ' , num2str(BitErrorRate_wp_trasp), ' %']);

% [number_w1p,ratio_w1p] = biterr(bitTx,bitRx_w1p);
% BitErrorRate_w1p = ratio_w1p*100; %percentage of error bits
[number_w1p_trasp,ratio_w1p_trasp] = biterr(bitTx,bitRx_w1p_trasp);
BitErrorRate_w1p_trasp = ratio_w1p_trasp*100; %percentage of error bits
% disp(['Bit Error Rate Considering Gain and TxPower between Tx-RIS = ' , num2str(BitErrorRate_w1p), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_w1p_trasp), ' %']);
disp(['Bit Error Rate Considering Gain and TxPower between Tx-RIS = ' , num2str(BitErrorRate_w1p_trasp), ' %']);

% 
% % show results for version 2
% [number2,ratio2] = biterr(bitTx,bitRx2);
% BitErrorRate2 = ratio2*100; %percentage of error bits
% [number_trasp2,ratio_trasp2] = biterr(bitTx,bitRx_trasp2);
% BitErrorRate_trasp2 = ratio_trasp2*100; %percentage of error bits
% 
% disp(['Bit Error Rate Considering TxPower in g and h = ' , num2str(BitErrorRate2), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_trasp2), ' %']);
% 
% [number_wd2,ratio_wd2] = biterr(bitTx,bitRx_wd2);
% BitErrorRate_wd2 = ratio_wd2*100; %percentage of error bits
% [number_wd_trasp2,ratio_wd_trasp2] = biterr(bitTx,bitRx_wd_trasp2);
% BitErrorRate_wd_trasp2 = ratio_wd_trasp2*100; %percentage of error bits
% 
% disp(['Bit Error Rate Considering Directivity = ' , num2str(BitErrorRate_wd2), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_wd_trasp2), ' %']);
% 
% [number_wp2,ratio_wp2] = biterr(bitTx,bitRx_wp2);
% BitErrorRate_wp2 = ratio_wp2*100; %percentage of error bits
% [number_wp_trasp2,ratio_wp_trasp2] = biterr(bitTx,bitRx_wp_trasp2);
% BitErrorRate_wp_trasp2 = ratio_wp_trasp2*100; %percentage of error bits
% 
% disp(['Bit Error Rate Without TxPower = ' , num2str(BitErrorRate_wp2), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_wp_trasp2), ' %']);
% 
% [number_w1p2,ratio_w1p2] = biterr(bitTx,bitRx_w1p2);
% BitErrorRate_w1p2 = ratio_w1p2*100; %percentage of error bits
% [number_w1p_trasp2,ratio_w1p_trasp2] = biterr(bitTx,bitRx_w1p_trasp2);
% BitErrorRate_w1p_trasp2 = ratio_w1p_trasp2*100; %percentage of error bits
% 
% disp(['Bit Error Rate Considering TxPower only in h = ' , num2str(BitErrorRate_w1p2), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_w1p_trasp2), ' %']);


%% Generazione segnale 

% Parameters for QAM modulation per subcarrier
bps = 4; %number of bits per simbol
M = 2^bps; % modulation order QAM

% Create OFDM modulator and demodulator objects 
fftLen = 64; %defyne the number of subcarriers
%fftLen = 256;
cpLen = fftLen/4; %cyclic prefix -> to make an OFDM signal intensive
                  % to time dispersion on the radio channel
nSym = 1; %number of symbol per RE (Resource Element)
%numGuardBandCarriers = [2; 2];
numGuardBandCarriers = [5; 4];
%pilotCarrierIdx = [20:15:120, 155:10:240]';
pilotCarrierIdx = [10:5:60]';
 if mod(fftLen,2) == 0
    NullIdx = (fftLen/2)+1; %NFFT is even.
else
    NullIdx = (fftLen+1)/2; %NFFT is odd.
end

numDataCarriers = ...
    fftLen - sum(numGuardBandCarriers) - length(pilotCarrierIdx) - 1;
%nullIdx = [1:6, 20:15:120, NullIdx, 155:10:240, 256-4:256];
nullIdx = [1:5, 10:5:30, NullIdx, 35:5:60, 64-3:64];
nt    = 1;   % Number of transmit antennas

ofdmMod = comm.OFDMModulator( ...
    "FFTLength",fftLen, ....
    "NumGuardBandCarriers",numGuardBandCarriers, ...
    "InsertDCNull",true, ...
    "PilotInputPort",true, ...
    "PilotCarrierIndices",pilotCarrierIdx, ...
    "CyclicPrefixLength",cpLen, ...
    "NumSymbols",nSym, ...
    "NumTransmitAntennas",nt);

cd = comm.ConstellationDiagram( ...    
    "ReferenceConstellation", qammod((0:M-1)', M), ...
    "XLimits", [-sqrt(M) sqrt(M)], ...
    "YLimits", [-sqrt(M) sqrt(M)]);
%show normalized constellation map
%scatterplot(qammod((0:M-1)', M)); 

clear dataIn;
dataIn = randi([0 1],1e3,1); %generated bit
Nsym = ceil(length(dataIn)/bps); %numero simboli di partenza
dataIn = [dataIn; zeros(Nsym*bps-length(dataIn),1)]; %pedding di zero per avere un numero intero di simboli

N_OFDMSym = ceil(Nsym/numDataCarriers); %numero simboli OFDM, cioè numero delle colonne
                                        
dataIn = [dataIn; zeros((N_OFDMSym*numDataCarriers - Nsym)*bps,1)]; %pedding per avere un numero 
                                                %intero di simboli OFDM

qamSig = qammod(dataIn,M,'InputType','bit');
%cd(qamSig);

qamSig = reshape(qamSig,[numDataCarriers,ceil(length(qamSig)/numDataCarriers)]);
clear inSig;
for j = 1:size(qamSig,2)
    inSig(:,j) = ofdmmod(qamSig(:,j),fftLen, ...
        cpLen,nullIdx');
end

RxSignal = conv2(inSig,Tot_Channel,'same'); %circular convolution
RxSignal_wd = conv2(inSig,Tot_Channel_wd,'same'); %circular convolution
RxSignal_wp = conv2(inSig,Tot_Channel_wp,'same'); %circular convolution
RxSignal_w1p = conv2(inSig,Tot_Channel_w1p,'same'); %circular convolution
RxSignal_trasp = conv2(inSig,Tot_Channel_trasp,'same'); %circular convolution
RxSignal_wd_trasp = conv2(inSig,Tot_Channel_wd_trasp,'same'); %circular convolution
RxSignal_wp_trasp = conv2(inSig,Tot_Channel_wp_trasp,'same'); %circular convolution
RxSignal_w1p_trasp = conv2(inSig,Tot_Channel_w1p_trasp,'same'); %circular convolution

clear outSig;
for j = 1:size(RxSignal,2)
    outSig(:,j) = ofdmdemod(RxSignal(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig = outSig(:);
dataOut = qamdemod(outSig,M,'OutputType','bit');
%cd(outSig);

clear outSig_wd;
for j = 1:size(RxSignal_wd,2)
    outSig_wd(:,j) = ofdmdemod(RxSignal_wd(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_wd = outSig_wd(:);
dataOut_wd = qamdemod(outSig_wd,M,'OutputType','bit');

clear outSig_wp;
for j = 1:size(RxSignal_wp,2)
    outSig_wp(:,j) = ofdmdemod(RxSignal_wp(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_wp = outSig_wp(:);
dataOut_wp = qamdemod(outSig_wp,M,'OutputType','bit');

clear outSig_w1p;
for j = 1:size(RxSignal_w1p,2)
    outSig_w1p(:,j) = ofdmdemod(RxSignal_w1p(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_w1p = outSig_w1p(:);
dataOut_w1p = qamdemod(outSig_w1p,M,'OutputType','bit');

%trasposto
clear outSig_trasp;
for j = 1:size(RxSignal_trasp,2)
    outSig_trasp(:,j) = ofdmdemod(RxSignal_trasp(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_trasp = outSig_trasp(:);
dataOut_trasp = qamdemod(outSig_trasp,M,'OutputType','bit');

clear outSig_wd_trasp;
for j = 1:size(RxSignal_wd_trasp,2)
    outSig_wd_trasp(:,j) = ofdmdemod(RxSignal_wd_trasp(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_wd_trasp = outSig_wd_trasp(:);
dataOut_wd_trasp = qamdemod(outSig_wd_trasp,M,'OutputType','bit');

clear outSig_wp_trasp;
for j = 1:size(RxSignal_wp_trasp,2)
    outSig_wp_trasp(:,j) = ofdmdemod(RxSignal_wp_trasp(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_wp_trasp = outSig_wp_trasp(:);
dataOut_wp_trasp = qamdemod(outSig_wp_trasp,M,'OutputType','bit');

clear outSig_w1p_trasp;
for j = 1:size(RxSignal_w1p_trasp,2)
    outSig_w1p_trasp(:,j) = ofdmdemod(RxSignal_w1p_trasp(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_w1p_trasp = outSig_w1p_trasp(:);
dataOut_w1p_trasp = qamdemod(outSig_w1p_trasp,M,'OutputType','bit');

%% received signal version 2
RxSignal2 = conv2(inSig,Tot_Channel2,'same'); %circular convolution
RxSignal_wd2 = conv2(inSig,Tot_Channel_wd2,'same'); %circular convolution
RxSignal_wp2 = conv2(inSig,Tot_Channel_wp2,'same'); %circular convolution
RxSignal_w1p2 = conv2(inSig,Tot_Channel_w1p2,'same'); %circular convolution
RxSignal_trasp2 = conv2(inSig,Tot_Channel_trasp2,'same'); %circular convolution
RxSignal_wd_trasp2 = conv2(inSig,Tot_Channel_wd_trasp2,'same'); %circular convolution
RxSignal_wp_trasp2 = conv2(inSig,Tot_Channel_wp_trasp2,'same'); %circular convolution
RxSignal_w1p_trasp2 = conv2(inSig,Tot_Channel_w1p_trasp2,'same'); %circular convolution

clear outSig2;
for j = 1:size(RxSignal2,2)
    outSig2(:,j) = ofdmdemod(RxSignal2(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig2 = outSig2(:);
dataOut2 = qamdemod(outSig2,M,'OutputType','bit');
%cd(outSig);

clear outSig_wd2;
for j = 1:size(RxSignal_wd2,2)
    outSig_wd2(:,j) = ofdmdemod(RxSignal_wd2(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_wd2 = outSig_wd2(:);
dataOut_wd2 = qamdemod(outSig_wd2,M,'OutputType','bit');

clear outSig_wp2;
for j = 1:size(RxSignal_wp2,2)
    outSig_wp2(:,j) = ofdmdemod(RxSignal_wp2(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_wp2 = outSig_wp2(:);
dataOut_wp2 = qamdemod(outSig_wp2,M,'OutputType','bit');

clear outSig_w1p2;
for j = 1:size(RxSignal_w1p2,2)
    outSig_w1p2(:,j) = ofdmdemod(RxSignal_w1p2(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_w1p2 = outSig_w1p2(:);
dataOut_w1p2 = qamdemod(outSig_w1p2,M,'OutputType','bit');

%trasposto
clear outSig_trasp2;
for j = 1:size(RxSignal_trasp2,2)
    outSig_trasp2(:,j) = ofdmdemod(RxSignal_trasp2(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_trasp2 = outSig_trasp2(:);
dataOut_trasp2 = qamdemod(outSig_trasp2,M,'OutputType','bit');

clear outSig_wd_trasp2;
for j = 1:size(RxSignal_wd_trasp2,2)
    outSig_wd_trasp2(:,j) = ofdmdemod(RxSignal_wd_trasp2(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_wd_trasp2 = outSig_wd_trasp2(:);
dataOut_wd_trasp2 = qamdemod(outSig_wd_trasp2,M,'OutputType','bit');

clear outSig_wp_trasp2;
for j = 1:size(RxSignal_wp_trasp2,2)
    outSig_wp_trasp2(:,j) = ofdmdemod(RxSignal_wp_trasp2(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_wp_trasp2 = outSig_wp_trasp2(:);
dataOut_wp_trasp2 = qamdemod(outSig_wp_trasp2,M,'OutputType','bit');

clear outSig_w1p_trasp2;
for j = 1:size(RxSignal_w1p_trasp2,2)
    outSig_w1p_trasp2(:,j) = ofdmdemod(RxSignal_w1p_trasp2(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig_w1p_trasp2 = outSig_w1p_trasp2(:);
dataOut_w1p_trasp2 = qamdemod(outSig_w1p_trasp2,M,'OutputType','bit');
%% SINR
SINR = tx_power_dBm_real - mag2db(abs(mean(mean(RxSignal_wp_trasp))));
%% 
RxSignal2 = conv2(inSig,h_SISO_WRIS,'same'); %circular convolution

for j = 1:size(RxSignal2,2)
    outSig2(:,j) = ofdmdemod(RxSignal2(:,j),fftLen, ...
        cpLen,0,nullIdx');
end

outSig2 = outSig2(:);
dataOut2 = qamdemod(outSig2,M,'OutputType','bit');

%% difference between Tx and Rx bit
[number,ratio] = biterr(dataIn,dataOut);
BitErrorRate = ratio*100; %percentage of error bits
[number_trasp,ratio_trasp] = biterr(dataIn,dataOut_trasp);
BitErrorRate_trasp = ratio_trasp*100; %percentage of error bits

disp(['Modulation OFDM with ', num2str(M),' QAM with RIS ']);
disp(['Bit Error Rate Considering TxPower in g and h = ' , num2str(BitErrorRate), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_trasp), ' %']);

[number_wd,ratio_wd] = biterr(dataIn,dataOut_wd);
BitErrorRate_wd = ratio_wd*100; %percentage of error bits
[number_wd_trasp,ratio_wd_trasp] = biterr(dataIn,dataOut_wd_trasp);
BitErrorRate_wd_trasp = ratio_wd_trasp*100; %percentage of error bits

disp(['Bit Error Rate Considering Directivity = ' , num2str(BitErrorRate_wd), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_wd_trasp), ' %']);

[number_wp,ratio_wp] = biterr(dataIn,dataOut_wp);
BitErrorRate_wp = ratio_wp*100; %percentage of error bits
[number_wp_trasp,ratio_wp_trasp] = biterr(dataIn,dataOut_wp_trasp);
BitErrorRate_wp_trasp = ratio_wp_trasp*100; %percentage of error bits

disp(['Bit Error Rate Without TxPower = ' , num2str(BitErrorRate_wp), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_wp_trasp), ' %']);

[number_w1p,ratio_w1p] = biterr(dataIn,dataOut_w1p);
BitErrorRate_w1p = ratio_w1p*100; %percentage of error bits
[number_w1p_trasp,ratio_w1p_trasp] = biterr(dataIn,dataOut_w1p_trasp);
BitErrorRate_w1p_trasp = ratio_w1p_trasp*100; %percentage of error bits

disp(['Bit Error Rate Considering TxPower only in h = ' , num2str(BitErrorRate_w1p), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_w1p_trasp), ' %']);

%% show results for version 2
[number2,ratio2] = biterr(dataIn,dataOut2);
BitErrorRate2 = ratio2*100; %percentage of error bits
[number_trasp2,ratio_trasp2] = biterr(dataIn,dataOut_trasp2);
BitErrorRate_trasp2 = ratio_trasp2*100; %percentage of error bits

disp(['Modulation OFDM with ', num2str(M),' QAM with RIS ']);
disp(['Bit Error Rate Considering TxPower in g and h = ' , num2str(BitErrorRate2), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_trasp2), ' %']);

[number_wd2,ratio_wd2] = biterr(dataIn,dataOut_wd2);
BitErrorRate_wd2 = ratio_wd2*100; %percentage of error bits
[number_wd_trasp2,ratio_wd_trasp2] = biterr(dataIn,dataOut_wd_trasp2);
BitErrorRate_wd_trasp2 = ratio_wd_trasp2*100; %percentage of error bits

disp(['Bit Error Rate Considering Directivity = ' , num2str(BitErrorRate_wd2), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_wd_trasp2), ' %']);

[number_wp2,ratio_wp2] = biterr(dataIn,dataOut_wp2);
BitErrorRate_wp2 = ratio_wp2*100; %percentage of error bits
[number_wp_trasp2,ratio_wp_trasp2] = biterr(dataIn,dataOut_wp_trasp2);
BitErrorRate_wp_trasp2 = ratio_wp_trasp2*100; %percentage of error bits

disp(['Bit Error Rate Without TxPower = ' , num2str(BitErrorRate_wp2), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_wp_trasp2), ' %']);

[number_w1p2,ratio_w1p2] = biterr(dataIn,dataOut_w1p2);
BitErrorRate_w1p2 = ratio_w1p2*100; %percentage of error bits
[number_w1p_trasp2,ratio_w1p_trasp2] = biterr(dataIn,dataOut_w1p_trasp2);
BitErrorRate_w1p_trasp2 = ratio_w1p_trasp2*100; %percentage of error bits

disp(['Bit Error Rate Considering TxPower only in h = ' , num2str(BitErrorRate_w1p2), ' %. If we consider the PHI transposed we obtain a BER = ' , num2str(BitErrorRate_w1p_trasp2), ' %']);

%% 
[number2,ratio2] = biterr(dataIn,dataOut2);
BitErrorRate2 = ratio2*100; %percentage of error bits
disp(['Modulation OFDM with ', num2str(M),' QAM without RIS']);
disp(['Bit Error Rate = ' , num2str(BitErrorRate2), ' %']);

%% Equalization during signal recovery
pathDelays = [0 4e-3 8e-3 12e-3];
sRate = 1000;
sampIdx = round(pathDelays/(1/sRate)) + 1;
hImp = complex(zeros(fftLen,1));
hImp(sampIdx) = mean(Tot_Channel,1);
hall = fftshift(fft(hImp));
dataIdx = setdiff((1:fftLen)',nullIdx');
h_datasubcarr = hall(dataIdx);
eqSig = outSig ./ h_datasubcarr;
cdScope(eqSig);
%% %% Generazione segnale
% % Parameters for QAM modulation per subcarrier
% bitsPerSymbol = 6; %number of bits per simbol
% modOrder = 2^bitsPerSymbol; % modulation order QAM
% % Create OFDM modulator and demodulator objects 
% fftLen = 256; %defyne the number of subcarriers
% cpLen = fftLen/4; %cyclic prefix -> to make an OFDM signal intensive
%                   % to time dispersion on the radio channel
% nSym = 1; %number of symbol per RE
% numGuardBandCarriers = [8; 7];
% pilotCarrierIdx = [20:15:120, 155:10:240]';
% numDataCarriers = ...
%     fftLen - sum(numGuardBandCarriers) - length(pilotCarrierIdx) - 1;
% ofdmMod = comm.OFDMModulator( ...
%     "FFTLength",fftLen, ....
%     "NumGuardBandCarriers",numGuardBandCarriers, ...
%     "InsertDCNull",true, ...
%     "PilotInputPort",true, ...
%     "PilotCarrierIndices",pilotCarrierIdx, ...
%     "CyclicPrefixLength",cpLen, ...
%     "NumSymbols",nSym, ...
%     "NumTransmitAntennas",numTx);
% ofdmDemod = comm.OFDMDemodulator(ofdmMod);
% ofdmDemod.NumReceiveAntennas = numRx;
% showResourceMapping(ofdmMod);
% cd = comm.ConstellationDiagram( ...    
%     "ReferenceConstellation", qammod((0:modOrder-1)', modOrder), ...
%     "XLimits", [-sqrt(modOrder) sqrt(modOrder)], ...
%     "YLimits", [-sqrt(modOrder) sqrt(modOrder)]);
% %show normalized constellation map
% scatterplot(qammod((0:modOrder-1)', modOrder)); 
% if mod(fftLen,2) == 0
%     NullIdx = (fftLen/2)+1; %NFFT is even.
% else
%     NullIdx = (fftLen+1)/2; %NFFT is odd.
% end
% numDataCarriers = ...
%     fftLen - sum(numGuardBandCarriers) - length(pilotCarrierIdx) - 1;
% nullIdx = [1:8, 20:15:120, NullIdx, 155:10:240, 256-6:256];
% data = randi([0 modOrder-1],numDataCarriers,nSym,numTx); %generazione simboli
% txSig = qammod(data,modOrder,'UnitAveragePower',true); %conversione bit in simboli QAM 
% % xp = randi([0 1],20*bitsPerSymbol,1);
% % pilotIn = qammod(xp,modOrder,'InputType','bit');
% % x = ofdmMod(txSig,pilotIn); %segnale modulato trasmesso
% outSig = ofdmmod(txSig,fftLen,cpLen,nullIdx'); %il numero di colonne indicano il numero di simboli OFDM trasmessi
% %Pass the signal through a noisy channel
% %rxSig = awgn(txSig,25);
% % txSig è una matrice di numeri complessi
% %la trasformo in un segnale moltiplicando la parte reale per il coseno e la
% %parte immaginaria per -seno
% T_bit = 1/Frequency; %bit duration
% T_sym = bitsPerSymbol*T_bit; %symbol duration
% f_sym = 1/T_sym;
% t = 0;
% TX_SIG = 0;
% for k = 1:length(txSig)
%     TX_SIG =[TX_SIG; real(txSig(k))*cos(2*pi*Frequency*((k-1)*T_sym:1e-12:k*T_sym)')+imag(txSig(k))*sin(2*pi*Frequency*((k-1)*T_sym:1e-12:k*T_sym)')];
%     t = [t,((k-1)*T_sym:1e-12:k*T_sym)];
% end
% figure(); plot(t,TX_SIG);
% %
% avgPower = mean(abs(rxSig).^2);
% sizeTx = size(txSig);
% scatterplot(reshape(txSig,1,sizeTx(1)*sizeTx(2))); 
% title('64-QAM Tx modulation');
% cd(txSig);
% % for i=1:length(sym)
% %     Re = real(sym(i));
% %     Im = imag(sym(i));
% % end
% %la parte reale va moltiplicata per cos(2pift)
% %la parte immaginaria va sfasata di pi/2 (sfasamento introdotto dalla
% %presenza dello j) e poi moltiplicato per sin(2pift)
% %filterCoeffs = rcosdesign(0.35, 4, bitsPerCarrier);
% %y = filter(filterCoeffs, 1, upsample(sym,bitsPerCarrier));
% %Create an error rate calculation object to compute bit error rate (BER).
% errRate = comm.ErrorRate;
% %Convert Eb/No value to anSNR.
% codeRate = 3/4;
% EbNo = 30;              % in dB
% SNR = convertSNR(EbNo,"ebno", ...
%   "BitsPerSymbol",bitsPerSymbol, ...
%   "CodingRate",codeRate);
% SNRLin = 10^(SNR/10);      % Linear
% %
% y = conv2(outSig,Tot_Channel,'same'); %segnale modulato ricevuto
% % TXSIG = fft(txSig);
% % TOT_CHANNEL = fft(Tot_Channel);
% % Y = TXSIG.* TOT_CHANNEL;
% % y = ifft(Y);
% %dataOut= ofdmdemod(y,fftLen,cpLen,cpLen,nullIdx');
% dataOut= ofdmdemod(y,fftLen,cpLen,cpLen,nullIdx');
% yd = qamdemod(dataOut,modOrder);
% %yp = qamdemod(pilotOut,modOrder);

%% funzione per il plot

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

function [mappe, xmap, ymap] = create_var_map(Z, rx_test, test_data)

xmap=       NaN(size(Z,1),size(Z,2));
ymap=       NaN(size(Z,1),size(Z,2));
mappe = struct;

%maps polulation
for jj = 1 : length(rx_test)
    for ii=1:size(test_data.n,1)
        mappe(jj).nmap(test_data.indices(ii,1),test_data.indices(ii,2),:)=test_data.n(ii,jj);
        mappe(jj).strengthmap(test_data.indices(ii,1),test_data.indices(ii,2),:)=test_data.SStrength(ii,jj);
        xmap(test_data.indices(ii,1),test_data.indices(ii,2),:)=test_data.Location(ii,1);
        ymap(test_data.indices(ii,1),test_data.indices(ii,2),:)=test_data.Location(ii,2);
    end
end
end

function helperViewArray(array)
%HELPERVIEWARRAY Visualize array geometry
%   helperViewArray(ARRAY) visualizes the array geometry of the arrayConfig
%   object, ARRAY.

%   Copyright 2020 The MathWorks, Inc. 

% Get element positions
[elPos, elSpacing] = getElementPosition(array);
arraySize = array.Size;

% Plot elements in y-z plane
figure('Position', [360 360 320 320])
scatter(elPos(2,:), elPos(3,:), 100, 'filled');
hold on; axis equal; grid off;

% Label elements
if arraySize(2) > 1 
    textXOffset = elSpacing(2)/10;
else
    textXOffset = elSpacing(1)/10;
end

if arraySize(1) > 1
    textYOffset = elSpacing(1)/10;
else
    textYOffset = elSpacing(2)/10;
end

for i = 1:size(elPos, 2)
    text(elPos(2,i) + textXOffset, elPos(3,i) - textYOffset, num2str(i));
end

% Remove ticks
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.XTick = [];
ax.XTickLabel = [];
ax.YTick = [];
ax.YTickLabel = [];

title('Array Geometry');
xlabel('y'); ylabel('z');
if arraySize(2) == 1
    xlim([-elSpacing(1), elSpacing(1)]);
else
    xlim([-elSpacing(2)*arraySize(2)/2, elSpacing(2)*arraySize(2)/2]);
end

if arraySize(1) == 1
    ylim([-elSpacing(2), elSpacing(2)]);
else
    ylim([-elSpacing(1)*arraySize(1)/2, elSpacing(1)*arraySize(1)/2]);
end

end

function [elPos, elSpacing] = getElementPosition(cfgArray)

elSpacing = repmat(cfgArray.ElementSpacing, ...
    [1, 2/length(cfgArray.ElementSpacing)]);

% Calculate element positions along rows and columns
numCol = cfgArray.Size(2);
numRow = cfgArray.Size(1);    
rowPos = (0:numCol-1)*elSpacing(2);
rowPos = rowPos - rowPos(end)/2;
colPos = (0:numRow-1)*elSpacing(1);
if numCol > 1
    colPos =  colPos(end)/2 - colPos;
else
    colPos =  colPos - colPos(end)/2;
end

% Formulate the position grid on the plane where the array panel lies
expRowPos = kron(rowPos, ones(1, numRow));
expColPos = repmat(colPos, 1, numCol);

% Formulate [x;y;z] positions
numEl = prod(cfgArray.Size);

elPos = [zeros(1, numEl); expRowPos; expColPos];

end




