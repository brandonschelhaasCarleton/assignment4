set(0, 'DefaultFigureWindowStyle', 'Docked')

clc
clear
close all

% ELEC 4700 - Assignment 4
% Brandon Schelhaas
% 101036851
% Programmed on MATLAB 2020a

global G C F;

doAll = 1;

% Assignment 1 Parameter Definitions
regionLength = 200; % nm
regionWidth = 100; % nm
boxes{1}.x = [80 120];
boxes{1}.y = [60 100];
boxes{1}.sigma = 0.01;
boxes{2}.x = [80 120];
boxes{2}.y = [0 40];
boxes{2}.sigma = 0.01;

% Set ny and nx, in nanometers, due to regionLength and regionWidth
ny = regionWidth;
nx = regionLength;
G = sparse(ny*ny, nx*nx);
Bv = zeros(nx*ny, 1); % voltage boundary

% ---------------------------
%        Part 1 and 2
% ---------------------------
doPart1and2 = 0;
if doPart1and2 | doAll
    % Part A - plot density, observe at 0.8V

    vin = 0.1:0.1:10;
    current_matrix = zeros(ny, nx);
    % i1 = zeros(ny, boxes{1}.x(1));
    % i2 = zeros(boxes{1}.y(2)-boxes{1}.y(1), boxes{1}.x(2)-boxes{1}.x(1));
    % i3 = zeros(boxes{1}.y(1)-boxes{2}.y(2), boxes{1}.x(2)-boxes{1}.x(1));
    % i4 = zeros(boxes{2}.y(2)-boxes{2}.y(1), boxes{2}.x(2)-boxes{2}.x(1));
    % i5 = zeros(ny, nx-boxes{1}.x(2));

    % Create vectors for mapped current and resistances
    i1_vec = zeros(1, length(vin));
    i2_vec = zeros(1, length(vin));
    i3_vec = zeros(1, length(vin));
    i4_vec = zeros(1, length(vin));
    i5_vec = zeros(1, length(vin));
    r1_vec = zeros(1, length(vin));
    r2_vec = zeros(1, length(vin));
    r3_vec = zeros(1, length(vin));
    r4_vec = zeros(1, length(vin));
    r5_vec = zeros(1, length(vin));
    i_vec = zeros(1, length(vin));
    r_vec = zeros(1, length(vin));

    % Run DC Sweep on Assignment 2 solution
    for i = 1:length(vin)
        % Apply the input voltage
        volt_app.x = vin(i);
        
        % Run Assignment 2's solution to find current in region
        assign2
        
        % Find the magnitude of current for each mesh
        current_matrix(:,:) = sqrt(current.x(:,:).*current.x(:,:) + current.y(:,:).*current.y(:,:));
        
        % Map the region's current into the separated regions
        i1 = current_matrix(1:ny, 1:boxes{1}.x(1));
        i2 = current_matrix(boxes{1}.y(1)+1:ny, boxes{1}.x(1)+1:boxes{1}.x(2));
        i3 = current_matrix(boxes{2}.y(2)+1:boxes{1}.y(1), boxes{1}.x(1)+1:boxes{1}.x(2));
        i4 = current_matrix(1:boxes{2}.y(2), boxes{1}.x(1)+1:boxes{1}.x(2));
        i5 = current_matrix(1:ny, boxes{1}.x(2)+1:nx);
        
        % Average the separated regions' currents
        i1_vec(i) = mean(i1, 'all');
        i2_vec(i) = mean(i2, 'all');
        i3_vec(i) = mean(i3, 'all');
        i4_vec(i) = mean(i4, 'all');
        i5_vec(i) = mean(i5, 'all');
        
        % Find Region 1's resistance
        v_high = volt_2d(:,1);
        v_low = volt_2d(:, boxes{1}.x(1));
        r1_vec(i) = (mean(v_high, 'all') - mean(v_low, 'all'))/(i1_vec(i));
        
        % Find Region 2's resistance
        v_high = volt_2d(boxes{1}.y(1)+1:ny, boxes{1}.x(1)+1);
        v_low = volt_2d(boxes{1}.y(1)+1:ny, boxes{1}.x(2));
        r2_vec(i) = (mean(v_high, 'all') - mean(v_low, 'all'))/(i2_vec(i));
        
        % Find Region 3's resistance
        v_high = volt_2d(boxes{2}.y(2)+1:boxes{1}.y(1), boxes{1}.x(1)+1);
        v_low = volt_2d(boxes{2}.y(2)+1:boxes{1}.y(1), boxes{1}.x(2));
        r3_vec(i) = (mean(v_high, 'all') - mean(v_low, 'all'))/(i3_vec(i));
        
        % Find Region 4's resistance
        v_high = volt_2d(1:boxes{2}.y(2), boxes{1}.x(1)+1);
        v_low = volt_2d(1:boxes{2}.y(2), boxes{1}.x(2));
        r4_vec(i) = (mean(v_high, 'all') - mean(v_low, 'all'))/(i4_vec(i));
        
        % Find Region 5's resistance
        v_high = volt_2d(:,boxes{1}.x(2)+1);
        v_low = volt_2d(:, nx);
        r5_vec(i) = (mean(v_high, 'all') - mean(v_low, 'all'))/(i5_vec(i));
        
        % Combine the series and parallel resistances
        r_vec(i) = r1_vec(i) + (1 / ( (1/r2_vec(i)) + (1/r3_vec(i)) + (1/r4_vec(i)) ) ) + r5_vec(i);
        
        % Calculate the current into the region
        i_vec(i) = vin(i)/r_vec(i);
        
        fprintf('%1d / %2d \n', i, length(vin));
    end
    
    % Plot current vs. input voltage
    figure
    plot(vin, i_vec)
    title('Current vs. Voltage');
    xlabel('Input Voltage [V]'); ylabel('Current [A]');
    
    % Fit a resistance value to the measured current
    calculatedRes = mean(r_vec, 'all');
    linearFit = polyfit(vin, i_vec, 1);
    fittedRes = 1/linearFit(1);
end

% ---------------------------
%           Part 3
% ---------------------------

% KCL Equations:
% (1) (v1-v2)/R1 + sC1(v1-v2) + Iin = 0
% (2) (v2-v1)/R1 + sC1(v2-v1) + v2/r2 + IL = 0
% (3) v3/r3 - IL = 0
% (4) (v4-v5)/r4 + Ix = 0
% (5) (v5-v4)/r4 + v5/ro = 0
% (6) v2-v3 = L (dIL/dt) --> -sL + v2-v3 = 0
% (7) v4 - (alpha/r3)*v3 = 0    (create a VCVS)
% (8) v1 = Vin

% Resources that helped create stamps
% (1) http://jefftrull.github.io/eda/eigen/c++/mna/2013/05/22/MNA-Eigen.html
% (2) https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119078388.app2

doPart3b = 1;
if doPart3b | doAll
    % Define number of nodes in circuit
    numNodes = 5;

    % Reset G, C, F matrices
    G = zeros(numNodes, numNodes);
    C = zeros(numNodes, numNodes);
    F = zeros(numNodes, 1);

    % Define component values
    R1 = 1;
    R2 = 2;
    if doPart1and2 | doAll
        R3 = fittedRes;
    else
        R3 = 10;
%         R3 = 186.9;
    end
    R4 = 0.1;
    Ro = 1000;
    C1 = 0.25;
    L1 = 0.2;
    alpha = 100;

    part3
    
    G
    C
end

% ---------------------------
%           Part 4
% ---------------------------

doPart4 = 1;
if doPart4 | doAll
    % Define number of nodes in circuit
    numNodes = 5;

    % Define component values
    R1 = 1;
    R2 = 2;
    if doPart1and2 | doAll
        R3 = fittedRes;
    else
        R3 = 10;
%         R3 = 186.9;
    end
    R4 = 0.1;
    Ro = 1000;
    C1 = 0.25;
    L1 = 0.2;
    alpha = 100;
    
    % Setup simulation parameters and sources
    simTime = 1; % 1s
    numSteps = 1000;
    dt = simTime/numSteps;
    time = linspace(0, 1, numSteps);
    freq = 1/0.03;
%     freq = 8;
%     freq = 100;
    
    % Create step, sine, and gauss sources
    v_step = zeros(1, numSteps);
    v_sin = zeros(1, numSteps);
    v_gauss = zeros(1, numSteps);
    
    step_transition = 0.03; % 0.03s
    v_step = double(time >= step_transition);
    
    v_sin = sin(2*pi*freq*time);

    gauss_dist = makedist('Normal', 'mu', 0.06, 'sigma', 0.03);
    gauss_pulse = pdf(gauss_dist, time); % makes a max amplitude of ~13.3
    v_pulse = (gauss_pulse.*sin(pi*time))/max(gauss_pulse.*sin(pi*time));

    plotSources = 0;
    if plotSources
        % Plot sources
        figure
        sgtitle('Plot of Input Voltage Sources');
        subplot(3,1,1)
        plot(time, v_step);
        title('Input Source 1: Step Voltage');
        xlabel('time [s]'); ylabel('Voltage [V]');
        xlim([0 simTime]); ylim([0 max(v_step)*1.1]);
        subplot(3,1,2)
        plot(time, v_sin);
        title('Input Source 2: Sine Wave');
        xlabel('time [s]'); ylabel('Voltage [V]');
        xlim([0 simTime]); ylim([min(v_sin)*1.1 max(v_sin)*1.1]);
        subplot(3,1,3)
        plot(time, v_pulse)
        title('Input Source 3: Gaussian Pulse');
        xlabel('time [s]'); ylabel('Voltage [V]');
        xlim([0 simTime]); ylim([min(v_pulse)*-1.1 max(v_pulse)*1.1]);
    end

    for k = 1:3 %3 sources
        if k == 1
            v_in = v_step;
        elseif k == 2
            v_in = v_sin;
        else
            v_in = v_pulse;
        end
        
        part4
    end
    
end

% ---------------------------
%           Part 5
% ---------------------------

doPart5 = 1;
if doPart5 | doAll
    % Define number of nodes in circuit
    numNodes = 5;

    % Define component values
    R1 = 1;
    R2 = 2;
    if doPart1and2 | doAll
        R3 = fittedRes;
    else
        R3 = 10;
%         R3 = 186.9;
    end
    R4 = 0.1;
    Ro = 1000;
    C1 = 0.25;
    L1 = 0.2;
    alpha = 100;
    
    % Setup simulation parameters and sources
    simTime = 1; % 1s
    numSteps = 1000;
    dt = simTime/numSteps;
    time = linspace(0, 1, numSteps);
    freq = 1/0.03;
    
    gauss_dist = makedist('Normal', 'mu', 0.06, 'sigma', 0.03);
    gauss_pulse = pdf(gauss_dist, time); % makes a max amplitude of ~13.3
    v_pulse = (gauss_pulse.*sin(pi*time))/max(gauss_pulse.*sin(pi*time));
    
    v_in = v_pulse;  
    part5
end

