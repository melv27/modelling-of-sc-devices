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

%clear all;
close all;

%%----- SETUP SIMULATION -----------------------------------------------------

% set temperature
T0       = 300;     % [K]

%secs1d_silicon_material_properties;
device.material = silicon_material_properties(T0);

% physical constants and parameters
secs1d_physical_constants;
%NA_array = [5e23,1e23,1e22,1e21,1e20]
device.doping.NA = 5e23; % [m^3]
device.doping.ND = 1e23; % [m^3]

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
voltage_step = 0.01;      %[V] 
voltage_start = -0.2;     %[V]
voltage_end = 1.2; %[V]

number_voltages = floor((voltage_end-voltage_start)/voltage_step)+1;


% this array stores all voltages for which the currents will be calculated
voltage_ramp = linspace(voltage_start,voltage_end, number_voltages);
%voltage_ramp = [0, 0.5];

find_zero_voltage = find(voltage_ramp == 0);
if (find_zero_voltage > 0)
    voltage_ramp(find(voltage_ramp == 0)) = voltage_step/4;
end;

current_ramp = zeros(1,length(voltage_ramp));
Efield_data = nan(device.mesh.Nelements,2); % dummy matrix to save electric field data in for plotting
potential_data = nan(device.mesh.Nelements+1,2); % same as Efield_data but for the potential
n_data = nan(device.mesh.Nelements+1,2); % same for n charge carrier density
p_data = nan(device.mesh.Nelements+1,2); % same for p charge carrier density
EFn_data = nan(device.mesh.Nelements+1,2); % fermi level in n-doped region
EFp_data = nan(device.mesh.Nelements+1,2); % fermi level in p-doped region
Ec_data = nan(device.mesh.Nelements+1,2); % conduction band
Ev_data = nan(device.mesh.Nelements+1,2); % valence band


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
        n_data(:,bias_voltage_num) = profile.n;
        p_data(:,bias_voltage_num) = profile.p;
        EFn_data(:,bias_voltage_num) = profile.EFn;
        EFp_data(:,bias_voltage_num) = profile.EFp;
        Ec_data(:,bias_voltage_num) = profile.Ec;
        Ev_data(:,bias_voltage_num) = profile.Ev;
    end;
end;



% for the sake of comparing two different voltages, all plots except for
% J-V curves are only plotted if exactly 2 voltage values are used in the
% voltage_ramp vector

if length(voltage_ramp) == 2
    for i = length(voltage_ramp):-1:1
        label1 = ['E for V= ', num2str(voltage_ramp(i)), ' V'];
        label2 = ['\Psi for V= ', num2str(voltage_ramp(i)), ' V'];
        
        if i == 1
            fig1 = 1; fig2 = 2; fig3 = 3;
        else
            fig1 = 4; fig2 = 5; fig3 = 6;
        end

        figure(fig1)
        scale = 1;
        plot2micron = 1E6; % scale from meter to micrometer
        set(fig1,'Position', [13 100 435 320]);
        hold on;
        semilogy(device.mesh.x*plot2micron, n_data(:,i), 'LineWidth',2,...
             'Color', [0 0 1]); % 'DisplayName', 'n'
        title({['Charge density profiles'],['@ T= ',num2str(T0), ' K', ' and V=', ...
            num2str(voltage_ramp(i)), ' V']},'FontSize', 14); 
        xlabel('position / {\mu m}','FontSize', 12);
        ylabel('density / {m^{-3}} ','FontSize', 12); 
        legend(['n'],'FontSize', 12);
        my_legend = legend;
        my_legend.Location = 'northeastoutside';
        axis tight
        xlim([24.5 25.5]);
        saveas(figure(fig1),['plots','/n_V',num2str(voltage_ramp(i)),...
            '_NA',num2str(device.doping.NA),'exc2.png']);
        hold off;
        %legend.location = 'northeastoutside';
        %legend show

        figure(fig2)
        scale = 1;
        plot2micron = 1E6; % scale from meter to micrometer
        set(fig2,'Position', [13 100 435 320]);
        hold on;
        semilogy(device.mesh.x*plot2micron, p_data(:,i), 'LineWidth',2,...
             'Color', [1 0 0]);
        title({['Charge density profiles'],['@ T= ',num2str(T0), ' K', ' and V=', ...
            num2str(voltage_ramp(i)), ' V']},'FontSize', 14); 
        xlabel('position / {\mu m}','FontSize', 12);
        ylabel('density / {m^{-3}} ','FontSize', 12); 
        legend(['p'],'FontSize', 12);
        my_legend = legend;
        my_legend.Location = 'northeastoutside';
        axis tight
        xlim([24.5 25.5]);
        saveas(figure(fig2),['plots','/p_V',num2str(voltage_ramp(i)),...
            '_NA',num2str(device.doping.NA),'exc2.png']);
        hold off;
        %legend.location = 'northeastoutside';
        %legend show
        
        figure(fig3)
        set(fig3,'Position', [13 500 435 320]);
        title({['Energy level diagram'],['@ T= ',num2str(T0), ' K', ' and V=', ...
            num2str(voltage_ramp(i)), ' V']},'FontSize', 14); 
        hold on;
        plot(device.mesh.x*plot2micron, EFn_data(:,i),  'LineWidth',2,'Color',...
              [0 0 scale]);
        plot(device.mesh.x*plot2micron, EFp_data(:,i), 'LineWidth',2,'Color',...
              [scale 0 0]);
        plot(device.mesh.x*plot2micron, Ec_data(:,i), 'LineWidth',1,'Color',...
              [0 0 scale]);
        plot(device.mesh.x*plot2micron, Ev_data(:,i), 'LineWidth',1,'Color',...
              [scale 0 0]);
        xlabel('position / {\mu m}','FontSize', 12);
        ylabel('potential or energy','FontSize', 12); 
        legend(['{E_{F,n}} /eV'],['{E_{F,p}} /eV'],...
               ['{E_{C}} /eV'],['{E_{V}} /eV'],'FontSize', 12)
        my_legend = legend;
        my_legend.Location = 'northeastoutside';    
        axis tight;
        xlim([24.5 25.5]);
        saveas(figure(fig3),['plots','/pn_jn_V',...
            num2str(voltage_ramp(i)),'_NA',num2str(device.doping.NA),'exc2.png']);
        hold off;

        figure(7)
        set(7,'Position', [1000 500 435 320]);
        title({['Electric field vs position'],['@ T= ',num2str(T0),' K']},'FontSize', 14); 
        hold on;
        plot(device.mesh.x(1:1000)*plot2micron, Efield_data(:,i),  'LineWidth',2,'DisplayName', label1)
        %legend (['V= ', num2str(voltage_ramp(i)), ' V']);
        xlabel('position / {\mu m}','FontSize', 12);
        ylabel('electric field / V/m ','FontSize', 12); 
        axis tight;
        xlim([24.5 25.5]);
        legend
        hold off;


        figure(8)
        set(8,'Position', [1000 500 435 320]);
        title({['Electrostatic potential vs position'],['@ T= ',num2str(T0),' K']},'FontSize', 14); 
        hold on;
        plot(device.mesh.x*plot2micron, potential_data(:,i),  'LineWidth',2,'DisplayName', label2)
        %legend (['V= ', num2str(voltage_ramp(i)), ' V']);
        xlabel('position / {\mu m}','FontSize', 12);
        ylabel('electrostatic potential / V ','FontSize', 12); 
        axis tight;
        xlim([24.5 25.5]);
        legend
        hold off;

    end
    saveas(figure(7),['plots','/Efield_V',num2str(voltage_ramp(2)),'_NA',...
        num2str(device.doping.NA),'exc2.png']);
    saveas(figure(8),['plots','/potential_V',num2str(voltage_ramp(2)),'_NA',...
        num2str(device.doping.NA),'exc2.png']);
end
 
if length(voltage_ramp) ~= 2
    figure(9)
    set(9,'Position', [1000 500 435 320]);
    title({['Current density vs voltage'],['@ T= ',num2str(T0),' K']},'FontSize', 14); 
    hold on;
    plot(voltage_ramp, abs(current_ramp))
    legend ('Jtot','FontSize', 12);
    xlabel('voltage  / V','FontSize', 12);
    ylabel('current density  / A{m^{-2}} ','FontSize', 12); 
    axis tight;
    saveas(figure(9),['plots','/current_density_NA',num2str(device.doping.NA),'exc2.png']);
    hold off;
    
    figure(10)
    set(10,'Position', [500 500 435 320]);
    title({['Current density vs voltage'],['@ T= ',num2str(T0),' K']},'FontSize', 14); 
    hold on;
    plot(voltage_ramp, abs(current_ramp),'LineWidth',3);
    scatter(voltage_ramp, abs(current_ramp));
    legend ('Jtot','FontSize', 12);
    xlabel('voltage  / V','FontSize', 12);
    ylabel('current density  / A{m^{-2}} ','FontSize', 12); 
    % this plots the current on a logarithmic scale
    set(gca,'yscale','log');
    saveas(figure(10),['plots','/current_density_log_NA',num2str(device.doping.NA),'exc2.png']);
    hold off;

end
current_ramp_inf = current_ramp;
close all;
