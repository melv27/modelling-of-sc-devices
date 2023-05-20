%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving 1D Poisson + Drift Diffusion semiconductor eqns using
%                    Gummel algorithm
%
%                 
%               
%
%     The code as is will calculate and plot a JV curve
%     as well as carrier densities, current densities, and electric field
%     distributions of a pn junction.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%%----- SETUP SIMULATION -----------------------------------------------------

% set temperature
T0       = 300;     % [K]

%secs1d_silicon_material_properties;
device.material = silicon_material_properties(T0);

% physical constants and parameters
secs1d_physical_constants;

device.doping.NA = 5E17; % [m^3]
device.doping.ND = 1E17; % [m^3]

% set device geometry
device.geometry.length         = 50e-6; % [m]
device.geometry.p_layer_length = 25e-6;  % [m]

% set device mesh

% uniform mesh
device.mesh.Nelements = 1000;
device.mesh.x = linspace (0, device.geometry.length, device.mesh.Nelements+1)';
% number of mesh cells
device.mesh.sinodes= [1:length(device.mesh.x)];
% position of junction
device.mesh.xm = device.geometry.p_layer_length; % device.geometry.length /2;

% set control parameters for simulation flow
% tolerances for convergence checks
itercontrol.tol = 1e-4;
itercontrol.maxit =  5000;
itercontrol.ptol = 1e-15;
itercontrol.pmaxit = 1000;

% if TRUE this parameter enable the use of previously calculate quantities as 
% initial guess for self-consistent solution of p,n, and psi
% if FALSE the thermal equilibrium is chosen as initial guess
% we will stay, for the time being with the thermal equilibrium as initial guess
itercontrol.prev_guess = false;

% set external voltage here!
voltage_step = 0;      %[V] 
voltage_start = 0;     %[V]
voltage_end = 0; %[V]

number_voltages = floor((voltage_end-voltage_start)/voltage_step)+1;


% this array stores all voltages for which the currents will be calculated
%voltage_ramp = linspace(voltage_start,voltage_end, number_voltages);
voltage_ramp = [0, 0.7];

find_zero_voltage = find(voltage_ramp == 0);
if (find_zero_voltage > 0)
    voltage_ramp(find(voltage_ramp == 0)) = voltage_step/4;
end;

current_ramp = zeros(1,length(voltage_ramp));
Efield_data = nan(device.mesh.Nelements,2); % dummy matrix to save electric field data in for plotting
potential_data = nan(device.mesh.Nelements+1,2); % same as Efield_data but for the potential


% as the current might be very small for small voltages, it is useful to calculate
% the current for large voltage first and proceed with decreasing voltages
%    you are welcome to test, whether the simulated current changes if one loops 
%    through voltages in increasing order

for bias_voltage_num = length(voltage_ramp):-1:1
    
    % ask for tighter convergence for small voltages
    if (voltage_ramp(bias_voltage_num) < 0.2)
        itercontrol.tol = 1e-7;
    else
        itercontrol.tol = 1e-4;
    end;

    itercontrol.prev_guess = false;

    disp(['Calculate j at voltage V=',num2str(voltage_ramp(bias_voltage_num)),' V']);

    [current_ramp(bias_voltage_num),profile, it, res] = current4voltage(voltage_ramp(bias_voltage_num),T0,device,itercontrol);
    
    if length(voltage_ramp) == 2
        Efield_data(:,bias_voltage_num) = profile.efield;
        potential_data(:,bias_voltage_num) = profile.psi;
    end;

end

figure(1)
scale = 1;
plot2micron = 1E6; % scale from meter to micrometer
set(1,'Position', [13 500 435 320]);
title({['Energy and potential profiles'],['@ T= ',num2str(T0)]}); 
hold on;
plot(device.mesh.x*plot2micron, profile.psi, 'LineWidth',2,...
     'Color', [0.5*scale 0.2*scale 0]); 
plot(device.mesh.x*plot2micron, profile.EFn,  'LineWidth',2,'Color',...
      [0 0 scale]);
plot(device.mesh.x*plot2micron, profile.EFp, 'LineWidth',2,'Color',...
      [scale 0 0]);
plot(device.mesh.x*plot2micron, profile.Ec, 'LineWidth',1,'Color',...
      [0 0 scale]);
plot(device.mesh.x*plot2micron, profile.Ev, 'LineWidth',1,'Color',...
      [scale 0 0]);
xlabel('position / {\mu m}');
ylabel('potential or energy'); 
legend(['{\Psi} /V'],['{E_{F,n}} /eV'],['{E_{F,p}} /eV'],...
       ['{E_{C}} /eV'],['{E_{V}} /eV'])
my_legend = legend;
my_legend.Location = 'northeastoutside';    
axis tight;
%saveas(figure(1),['plots','/pn_jn_NA',num2str(device.doping.NA),'exc1_inf.png']);
hold off;

figure(2)
scale = 1;
plot2micron = 1E6; % scale from meter to micrometer
set(2,'Position', [13 100 435 320]);
%hold on;
semilogy(device.mesh.x*plot2micron, profile.n, 'LineWidth',2,...
     'Color', [0 0 1]); 
title({['Charge density profiles'],['@ T= ',num2str(T0)]}); 
xlabel('position / {\mu m}');
ylabel('density / {m^{-3}} '); 
legend('n');
my_legend = legend;
my_legend.Location = 'northeastoutside';    
axis tight;
%saveas(figure(2),['plots','/n_profiles_NA',num2str(device.doping.NA),'exc1_inf.png']);
%hold off;

figure(3)
scale = 1;
plot2micron = 1E6; % scale from meter to micrometer
set(3,'Position', [500 100 435 320]);
%hold on
semilogy(device.mesh.x*plot2micron, profile.p,  'LineWidth',2,'Color',...
      [1 0 0]);
title({['Charge density profiles'],['@ T= ',num2str(T0)]}); 
xlabel('position / {\mu m}');
ylabel('density / {m^{-3}} '); 
legend('p');
my_legend = legend;
my_legend.Location = 'northeastoutside';    
axis tight;
%saveas(figure(3),['plots','/p_profiles_NA',num2str(device.doping.NA),'exc1_inf.png']);
%hold off

figure(4)
set(4,'Position', [1000 500 435 320]);
title({['Current density vs voltage'],['@ T= ',num2str(T0),' K']}); 
hold on;
plot(voltage_ramp, abs(current_ramp))
legend ('Jtot' );
xlabel('voltage  / V');
ylabel('current density  / A{m^{-2}} '); 
axis tight;
%saveas(figure(4),['plots','/current_density_NA',num2str(device.doping.NA),'exc1_inf.png']);
hold off;

figure(5)
set(5,'Position', [500 500 435 320]);
title({['Current density vs voltage'],['@ T= ',num2str(T0),' K']}); 
hold on;
plot(voltage_ramp, abs(current_ramp),'LineWidth',3);
scatter(voltage_ramp, abs(current_ramp));
legend ('Jtot' );
xlabel('voltage  / V');
ylabel('current density  / A{m^{-2}} '); 
%this plots the current on a logarithmic scale
set(gca,'yscale','log');
%saveas(figure(5),['plots','/current_density_log_NA',num2str(device.doping.NA),'exc1_inf.png']);
hold off;

% for the sake of comparing two different voltages, all plots except for
% J-V curves are only plotted if exactly 2 voltage values are used in the
% voltage_ramp vector

if length(voltage_ramp) == 2
    for i = length(voltage_ramp):-1:1
        label = ['V= ', num2str(voltage_ramp(i)), ' V'];
        figure(6)
        set(6,'Position', [1000 500 435 320]);
        title({['Electric field vs position'],['@ T= ',num2str(T0),' K']}); 
        hold on;
        plot(device.mesh.x(1:1000)*plot2micron, Efield_data(:,i),  'LineWidth',2,'DisplayName', label)
        %legend (['V= ', num2str(voltage_ramp(i)), ' V']);
        xlabel('position / {\mu m}');
        ylabel('electric field / V/m '); 
        axis tight;
        %saveas(figure(6),['plots','/Efield.png']);
        hold off;
        legend show
        
        figure(7)
        set(7,'Position', [1000 500 435 320]);
        title({['Electrostatic potential vs position'],['@ T= ',num2str(T0),' K']}); 
        hold on;
        plot(device.mesh.x*plot2micron, potential_data(:,i),  'LineWidth',2,'DisplayName', label)
        %legend (['V= ', num2str(voltage_ramp(i)), ' V']);
        xlabel('position / {\mu m}');
        ylabel('electrostatic potential / V '); 
        axis tight;
        %saveas(figure(7),['plots','/potential.png']);
        hold off;
        legend show
    end
end

%close all;
%Efield_data(:,1)
