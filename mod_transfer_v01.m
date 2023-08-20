clc;
clear;
colordef black; format long;
close all;

%% Parameter Setting

lambda = 432.0; %546.1; % Unit : nm
sk16_schott = 1.62286; %%“N”-glass as an environmentally friendly alternative to conventional lead and arsenic-containing glass types.
n_bk7_schott = 1.51872; % SCHOTT N-BK7
surface_num = 6;
distance = [0.01, 0.02, 0.035, 0.045, 0.061, 0.072, 0.3];   % Unit : mm
material = [1, n_bk7_schott, 1, sk16_schott, 1, n_bk7_schott, 1]; % Unit : mm
y_radius = [0.3, -0.3, -0.245, 0.245, 0.3, -0.3]; % Unit : mm
aperture = 0.075;   % Unit : mm


ang_x = 0;
ang_y = 0;
cross_diameter_num = 201;

%------------------------ Application Switch ------------------------%

Use_Paraxial_Solve = 1;         % 0 = OFF, 1 = ON

View_Lens = 1;                  % 0 = OFF, 1 = ON
    viewplane = 1;              % 1 = XZ, 2 = YZ, 3 = 3D
    display_line = 21;
Spot_Diagram = 1;               % 0 = OFF, 1 = ON
Transmission_Plane = 1;         % 0 = OFF, 1 = ON
Line_Spread_Function = 1;       % 0 = OFF, 1 = ON
MTF = 1;                        % 0 = OFF, 1 = ON
Point_Spread_Function = [1; 1]; % 0 = OFF, 1 = ON; perspective:[yz, xy]

%% Source Setting
lambda = lambda*1e-6;   % nm -> mm
[s_x, s_y, s_z, L, M, N] = light_source_setting(aperture,distance,cross_diameter_num,ang_x,ang_y);

%% Calculate Paraxial Focal Length
[BFL, EFL] = paraxial_focal_length(surface_num,distance,material,y_radius);
disp(['BFL = ',num2str(BFL),', EFL = ',num2str(EFL)])
save BFL BFL;
if Use_Paraxial_Solve == 1
    distance(end+1) = BFL-distance(end);
    material(end+1) = 1;
    y_radius(end+1) = inf;
end

%% Lens data
Lens.lambda = lambda*1e-6;
Lens.surface_num = surface_num;
Lens.distance = distance;
Lens.material = material;
Lens.y_radius = y_radius;
Lens.aperture = aperture;

%% Ray Tracing
curvature = 1./y_radius;
s_x_all = cell(1,numel(distance)); s_y_all = cell(1,numel(distance)); s_z_all = cell(1,numel(distance));
delta = zeros(size(s_x,1),size(s_x,2));
Opti_Path_Diff = zeros(size(s_x,1),size(s_x,2));

for i = 1:numel(distance)
    if i == numel(distance)
        z0 = ones(size(z0,1),size(z0,2))*sum(distance);
        x0 = s_x+(L./N).*(z0-s_z);
        y0 = s_y+(M./N).*(z0-s_z);

        x = [s_x;x0];    y = [s_y;y0];    z = [s_z;z0];
        s_x_all{i} = x;  s_y_all{i} = y;  s_z_all{i} = z;
    else
        z0 = s_z+distance(i)-delta;
        x0 = s_x+(L./N).*(z0-s_z);
        y0 = s_y+(M./N).*(z0-s_z);

        B = N-curvature(i).*(L.*x0+M.*y0);
        C = curvature(i).*(x0.^2+y0.^2);
        delta = C./(B+sqrt(B.^2-curvature(i).*C));

        if i <= surface_num
            Opti_Path_Diff = Opti_Path_Diff+abs(delta);
        end

        x1 = x0+L.*delta; y1 = y0+M.*delta; z1 = z0+N.*delta;
        x = [s_x;x1];    y = [s_y;y1];    z = [s_z;z1];
        s_x_all{i} = x;  s_y_all{i} = y;  s_z_all{i} = z;

        CosInc = sqrt(B.^2-curvature(i).*C);
        nTrans_CosTrans = sqrt((material(i+1).^2)-((material(i).^2).*(1-CosInc.^2)));
        k = curvature(i).*(nTrans_CosTrans-material(i).*CosInc);

        L_Trans = (material(i).*L-k.*x1)./material(i+1); L = L_Trans;
        M_Trans = (material(i).*M-k.*y1)./material(i+1); M = M_Trans;
        N_Trans = sqrt(1-(L_Trans.^2+M_Trans.^2));       N = N_Trans;

        s_x = x1; s_y = y1; s_z = z1;
    end
end

Data = data_reshape(s_x_all,s_y_all,s_z_all,cross_diameter_num);

%% View Lens
if View_Lens == 1
    display_tools.view_lens(Lens, Data, display_line, viewplane)
end

%% Spot Diagram
if Spot_Diagram == 1
    display_tools.spot_diagram(Data)
end

%% Optical Path Difference (OPD) at Transmission Plane
trans_plane_data = trans_plane_position_and_optical_path(surface_num, distance, material, Data, L, M, N);

if Transmission_Plane == 1
    display_tools.transmission_plane(trans_plane_data)
end

%% Paraxial Solve
if Use_Paraxial_Solve == 1
    focal_plane_position = sum(distance(end-1:end));
else
    focal_plane_position = distance(end);
end
diffra_limit = diffraction_limit(lambda,aperture,BFL,EFL,focal_plane_position);

%% Line Spread Function
if Line_Spread_Function == 1
    LSF_data = line_spread_function(lambda, aperture, trans_plane_data, focal_plane_position);
    display_tools.line_spread_function(LSF_data,diffra_limit)
end

%% MTF
if MTF == 1
    if Line_Spread_Function == 0
        LSF_data = line_spread_function(lambda, aperture, trans_plane_data, focal_plane_position);
    end
    display_tools.MTF(Lens,LSF_data,diffra_limit,ang_y)
end

%% Point Spread Function
if sum(Point_Spread_Function) > 0
    PSF_data = point_spread_function(lambda, aperture, trans_plane_data, focal_plane_position, Point_Spread_Function);
    display_tools.point_spread_function(Point_Spread_Function,PSF_data,trans_plane_data,focal_plane_position)
end

 % % %   data reshape %%%%
function Data = data_reshape(s_x_all,s_y_all,s_z_all,cross_diameter_num)
for i = 1:numel(s_x_all)
    Data.X_1{i} = reshape(s_x_all{i}(1:cross_diameter_num,:),1,cross_diameter_num^2);
    Data.X_1{i}(2,:) = reshape(s_x_all{i}(cross_diameter_num+1:end,:),1,cross_diameter_num^2);
    Data.Y_1{i} = reshape(s_y_all{i}(1:cross_diameter_num,:),1,cross_diameter_num^2);
    Data.Y_1{i}(2,:) = reshape(s_y_all{i}(cross_diameter_num+1:end,:),1,cross_diameter_num^2);
    Data.Z_1{i} = reshape(s_z_all{i}(1:cross_diameter_num,:),1,cross_diameter_num^2);
    Data.Z_1{i}(2,:) = reshape(s_z_all{i}(cross_diameter_num+1:end,:),1,cross_diameter_num^2);

    Data.X_2{1,i} = s_x_all{i}(1:cross_diameter_num,:);
    Data.X_2{2,i} = s_x_all{i}(cross_diameter_num+1:end,:);
    Data.Y_2{1,i} = s_y_all{i}(1:cross_diameter_num,:);
    Data.Y_2{2,i} = s_y_all{i}(cross_diameter_num+1:end,:);
    Data.Z_2{1,i} = s_z_all{i}(1:cross_diameter_num,:);
    Data.Z_2{2,i} = s_z_all{i}(cross_diameter_num+1:end,:);
end
end

%%light source setting %%%%
function [s_x, s_y, s_z, L, M, N] = light_source_setting(aperture,distance,cross_diameter_num,ang_x,ang_y)
if mod(cross_diameter_num,2)==0
    cross_diameter_num = cross_diameter_num+1;
end

x = linspace(-aperture/2,aperture/2,cross_diameter_num);
y = linspace(-aperture/2,aperture/2,cross_diameter_num);

[s_x,s_y] = meshgrid(x,y);
r = sqrt(s_x.^2+s_y.^2);

s_x(r>aperture/2) = nan; s_y(r>aperture/2) = nan;
% s_x = reshape(s_x,1,numel(s_x)); s_y = reshape(s_y,1,numel(s_y));
% s_x(isnan(s_x))=[]; s_y(isnan(s_y))=[];
s_z = zeros(size(s_x,1),size(s_x,2));
L = ones(size(s_x,1),size(s_x,2))*sind(ang_x); M = ones(size(s_x,1),size(s_x,2))*sind(ang_y); N = sqrt(1-((L.^2)+(M.^2)));

if ang_x ~= 0 || ang_y ~= 0
    s_x = s_x+(L./N).*(-distance(1));
    s_y = s_y+(M./N).*(-distance(1));
end
end


% % % Diffraction limit%%%%%%%%%%%%%
function diffra_limit = diffraction_limit(lambda,aperture,BFL,EFL,focal_plane_position)

z = focal_plane_position+(EFL-BFL);
size = round(2*aperture^2/(lambda*z))*2*10;
if mod(size,2)==0
    size = size+1;
end

y = linspace(-aperture/2,aperture/2,size);
f_number = z/aperture;
cutoff_freq = 1/(lambda*f_number);

theata = atan(aperture/2/z);
theata_all = linspace(-theata,theata,size);
I = (sin((aperture/2)*2*pi*sin(theata_all)./lambda)./((aperture/2)*2*pi*sin(theata_all)./lambda)).^2;
I(isnan(I)) = 1;
diffra_limit.I = I;

LSF = I;
OTF = fftshift(fft(LSF));
MTF = abs(OTF);
MTF = MTF./max(MTF);


a = linspace(-size/2,size/2,size);
f = a*(1/aperture);

diffra_limit.MTF = MTF;
diffra_limit.f = f;
diffra_limit.y = y;
diffra_limit.cutoff_freq = cutoff_freq;


end

% % % LSF %%%
function LSF_data = line_spread_function(lambda, aperture, trans_plane_data, focal_plane_position)
%%
k0 = 2*pi/lambda;
center_index = (size(trans_plane_data.y,2)+1)/2;
trans_plane_data.OP(isnan(trans_plane_data.OP)) = 0;
trans_plane_data.y(isnan(trans_plane_data.y)) = 0;

amp = trans_plane_data.OP(:,center_index); amp(amp~=0) = 1;
phase = k0*trans_plane_data.OP(:,center_index);

phase_Mask = amp.*exp(1i*phase);

%%
y = trans_plane_data.y(:,center_index);
% aperture_radius = max(trans_plane_data.y(:,center_index));
aperture_radius = aperture/2;
Monitor_y = linspace(-aperture_radius,aperture_radius,length(y));
% Monitor_y = linspace(-0.2,0.2,length(y));


%%
Distance= focal_plane_position;

Monitor = zeros(length(Monitor_y),1);
% const = 1/(1i*lambda);
parfor i = 1:length(Monitor_y)
%     R = ((Monitor_y(i)-y).^2)/(2*Distance);
%     Monitor(i) = sum((phase_Mask.*exp(1i*k0*R)+Distance),'all');

    R = sqrt(Distance^2+(Monitor_y(i)-y).^2);
    Monitor(i) = sum((phase_Mask.*exp(1i*k0*R)),'all');
%     Monitor(i) = sum((phase_Mask.*exp(1i*k0*R).*Distance./R.^2),'all');
end

LSF_data.Monitor_y = Monitor_y;
LSF_data.Intensity = abs(Monitor);
LSF_data.Intensity_normalize = LSF_data.Intensity/sum(amp,'all');
LSF_data.Power = abs(Monitor).^2;
LSF_data.Power_normalize = LSF_data.Power/sum(amp,'all')^2;
LSF_data.Strehl_ratio = max(LSF_data.Power/sum(amp,'all').^2);
end

% % % paraxial focal length  %%%%%

function [BFL, EFL] = paraxial_focal_length(surface_num,distance,material,sur_radius)
M = eye(2,2);

for i = 1:surface_num
    M_traslation = [1, distance(i); 0, 1];
    M_refraction = [1, 0; (material(i)-material(i+1))/(material(i+1)*sur_radius(i)), material(i)/material(i+1)];
    M = M_refraction * M_traslation * M;
end

A = M(1,1); B = M(1,2); C = M(2,1); D = M(2,2);
q = -A/C;
s = (1-A)/C;
f_s = q-s;

BFL = q; EFL = f_s;
end

% %%  % psf %%%%

function PSF_data = point_spread_function(lambda, aperture, trans_plane_data, propagate_distance, perspective_switch)

k0 = (2*pi/lambda);
amp = trans_plane_data.OP;
amp(~isnan(amp))=1; amp(isnan(amp))=0;

trans_plane_data.x(isnan(trans_plane_data.x)) = 0;
trans_plane_data.y(isnan(trans_plane_data.y)) = 0;
trans_plane_data.OP(isnan(trans_plane_data.OP)) = 0;

if perspective_switch(1) == 1
% %    y_bound = [min(min(trans_plate_data.y)),max(max(trans_plate_data.y))];
    y_bound = [-aperture/2,aperture/2];

    Monitor_Boundary = [y_bound(1), y_bound(end), size(trans_plane_data.y,1); ...
        (propagate_distance+propagate_distance*0.3)/200, propagate_distance+propagate_distance*0.3, (propagate_distance+propagate_distance*0.3)/200];   %[Y_Min, Y_Max, Y_Step; Z_Min, Z_Max, Z_Increment]
    Monitor1_y = linspace(Monitor_Boundary(1,1),Monitor_Boundary(1,2),Monitor_Boundary(1,3));
    Monitor1_z = Monitor_Boundary(2,1):Monitor_Boundary(2,3):Monitor_Boundary(2,2);


    phase = k0*trans_plane_data.OP;
    phase_mask = amp.*exp(1i*phase);

    Monitor = zeros(length(Monitor1_y),length(Monitor1_z));
    tic
    parfor i = 1:length(Monitor1_z)
        Monitor_xy = zeros(length(phase_mask),1);
        for ii = 1:length(Monitor1_y)
            R = sqrt(Monitor1_z(i).^2+(trans_plane_data.y-Monitor1_y(ii)).^2+(trans_plane_data.x).^2);
            Monitor_xy(ii) = sum((phase_mask.*exp(1i*k0*R)),'all');
        end
        Monitor(:,i) = Monitor_xy;
    end
    toc

    [~, index_z] = min(abs(Monitor1_z-trans_plane_data.dz-propagate_distance));

    PSF_data.Monitor1_y = Monitor1_y;
    PSF_data.Monitor_z = Monitor1_z;
    PSF_data.Intensity = abs(Monitor);
    PSF_data.Intensity_normalize = PSF_data.Intensity/max(max(PSF_data.Intensity));
    PSF_data.Power = abs(Monitor).^2;
    PSF_data.Power_normalize = PSF_data.Power/max(max(PSF_data.Power));
    PSF_data.Strehl_ratio = max(max(PSF_data.Power(:,index_z)))/sum(amp,'all').^2;
end

% % % XY plane
if perspective_switch(2) == 1
    %-------------------------------------------------------
    x_bound = [-aperture/2,aperture/2];
    y_bound = [-aperture/2,aperture/2];

    Monitor2_x = linspace(x_bound(1),x_bound(end),size(trans_plane_data.x,1));
    Monitor2_y = linspace(y_bound(1),y_bound(end),size(trans_plane_data.y,1));
    %-------------------------------------------------------
    phase = k0*trans_plane_data.OP;
    phase_mask = amp.*exp(1i*phase);

    Monitor_xy = zeros(length(Monitor2_x),length(Monitor2_y));
    tic
    for i = 1:length(Monitor2_x)
        parfor ii = 1:length(Monitor2_y)
            R = sqrt(propagate_distance^2+(trans_plane_data.y-Monitor2_y(ii)).^2+(trans_plane_data.x-Monitor2_x(i)).^2);
            Monitor_xy(i,ii) = sum((phase_mask.*exp(1i*k0*R)),'all');
        end
    end
    toc

    PSF_data.Monitor2_x = Monitor2_x;
    PSF_data.Monitor2_y = Monitor2_y;
    PSF_data.Intensity_xy = abs(Monitor_xy);
    PSF_data.Intensity_normalize_xy = PSF_data.Intensity_xy/max(max(PSF_data.Intensity_xy));
    PSF_data.Power_xy = abs(Monitor_xy).^2;
    PSF_data.Power_normalize_xy = PSF_data.Power_xy/max(max(PSF_data.Power_xy));
    PSF_data.Strehl_ratio_xy = max(max(PSF_data.Power_xy))/sum(amp,'all').^2;
end
end

% % %
% % transmission plane position
function [trans_plane_data] = trans_plane_position_and_optical_path(surface_num, distance, material, Data, L, M, N)
OP = zeros(size(Data.X_2{1,1},1),size(Data.X_2{1,1},2));
for i = 1:surface_num
    OP = OP+sqrt((Data.X_2{2,i}-Data.X_2{1,i}).^2+(Data.Y_2{2,i}-Data.Y_2{1,i}).^2 ...
             +(Data.Z_2{2,i}-Data.Z_2{1,i}).^2)*material(i);
end

x0 = Data.X_2{2,surface_num}; y0 = Data.Y_2{2,surface_num}; z0 = Data.Z_2{2,surface_num};
dz = max(max(z0))-z0;
x1 = x0+(L./N).*dz;
y1 = y0+(M./N).*dz;
z1 = max(max(z0));
OP = OP+sqrt((x1-x0).^2+(y1-y0).^2+(z1-z0).^2)*material(surface_num+1);

trans_plane_data.x = x1;
trans_plane_data.y = y1;
trans_plane_data.z = z1;
trans_plane_data.dz = z1-sum(distance(1:surface_num));
trans_plane_data.OP = OP;
end
