% Normal Simulation

vout = zeros(1, numSteps);
Cn = 0.00001;
In = 0.001;
for i = 1:numSteps
    G = zeros(numNodes, numNodes);
    C = zeros(numNodes, numNodes);
    F = zeros(numNodes, 1);
    
    % Setup G, C, and F matrices through stamp functions
    setResistor(1, 2, R1); % Set R1 -- V(1)
    setResistor(2, 0, R2); % Set R2 -- V(2)
    setResistor(3, 0, R3); % Set R3 -- V(3)
    setResistor(4, 5, R4); % Set R4 -- V(4)
    setResistor(5, 0, Ro); % Set R5 -- V(5)
    setCapacitor(1, 2, C1); % Set C1
    setCapacitor(0, 3, Cn); % Set Cn
    setInductor(2, 3, L1); % Set L1 -- V(6) or IL in V
    setVCVS(3, 0, 4, 0, alpha/R3); % Set VCVS (V = alpha*I3 = (alpha/R3)*V3 -- V(7) or Ix in V (current through VCVS)

    % Everytime voltage source is setup, G matrix is changed, so circuit needs to be redone (put above)
    setVoltageSource(1, 0, v_in(i)); % V(8) or Iin in V
    setCurrentSource(3, 0, In);
    
    % Solve for unknowns
    if i == 1
        V_new = (G + (C./dt)) \ F; % no old voltage, so C/dt * 0 = 0
    else
        V_new = (G + (C./dt)) \ (F + (C./dt)*V_old);
    end
    
    % Take output voltage from unknowns
    vout(i) = V_new(5); % V(5) = vout
    
    % Set old to new voltage for next timestep
    V_old = V_new;
end

% print C (report requirement)
C

% set frequency x-axis
fshift = (-length(v_in)/2:length(v_in)/2-1);

% plot time and freq domain
figure
sgtitle('Noise Current Source: Constant Magnitude');
subplot(1, 2, 1)
plot(time, v_in, time, vout);
title('Pulse Response');
xlabel('time [s]'); ylabel('Voltage [V]');
legend('Vin', 'Vout');
subplot(1, 2, 2)
plot(fshift, mag2db(abs(fftshift(fft(v_in)))));
hold on;
plot(fshift, mag2db(abs(fftshift(fft(vout)))));
hold on;
title('Frequency Spectrum');
xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
legend('Vin', 'Vout');

% Add random noise

% Generate random currents and plot
rand_currents = In * randn(1, length(v_in));
figure
histogram(rand_currents)
title('Random Noise Current Distribution');
xlabel('In [A]'); ylabel('Number');

for i = 1:numSteps
    G = zeros(numNodes, numNodes);
    C = zeros(numNodes, numNodes);
    F = zeros(numNodes, 1);
    
    % Setup G, C, and F matrices through stamp functions
    setResistor(1, 2, R1); % Set R1 -- V(1)
    setResistor(2, 0, R2); % Set R2 -- V(2)
    setResistor(3, 0, R3); % Set R3 -- V(3)
    setResistor(4, 5, R4); % Set R4 -- V(4)
    setResistor(5, 0, Ro); % Set R5 -- V(5)
    setCapacitor(1, 2, C1); % Set C1
    setCapacitor(0, 3, Cn); % Set Cn
    setInductor(2, 3, L1); % Set L1 -- V(6) or IL in V
    setVCVS(3, 0, 4, 0, alpha/R3); % Set VCVS (V = alpha*I3 = (alpha/R3)*V3 -- V(7) or Ix in V (current through VCVS)

    % Everytime voltage source is setup, G matrix is changed, so circuit needs to be redone (put above)
    setVoltageSource(1, 0, v_in(i)); % V(8) or Iin in V
    setCurrentSource(3, 0, rand_currents(i));
    
    if i == 1
        V_new = (G + (C./dt)) \ F; % no old voltage, so C/dt * 0 = 0
    else
        V_new = (G + (C./dt)) \ (F + (C./dt)*V_old);
    end
    vout(i) = V_new(5); % V(5) = vout
    
    V_old = V_new;
end

fshift = (-length(v_in)/2:length(v_in)/2-1);

figure
sgtitle('Noise Current Source: Normally-Distributed');
subplot(1, 2, 1)
plot(time, v_in, time, vout);
title('Pulse Response');
xlabel('time [s]'); ylabel('Voltage [V]');
legend('Vin', 'Vout');
subplot(1, 2, 2)
plot(fshift, mag2db(abs(fftshift(fft(v_in)))));
hold on;
plot(fshift, mag2db(abs(fftshift(fft(vout)))));
hold on;
title('Frequency Spectrum');
xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
legend('Vin', 'Vout');

% Vary Cn

Cn = [0.00001 0.1 0.01 0.001];
figure
for k = 1:length(Cn)
    for i = 1:numSteps
        G = zeros(numNodes, numNodes);
        C = zeros(numNodes, numNodes);
        F = zeros(numNodes, 1);

        % Setup G, C, and F matrices through stamp functions
        setResistor(1, 2, R1); % Set R1 -- V(1)
        setResistor(2, 0, R2); % Set R2 -- V(2)
        setResistor(3, 0, R3); % Set R3 -- V(3)
        setResistor(4, 5, R4); % Set R4 -- V(4)
        setResistor(5, 0, Ro); % Set R5 -- V(5)
        setCapacitor(1, 2, C1); % Set C1
        setCapacitor(0, 3, Cn(k)); % Set Cn
        setInductor(2, 3, L1); % Set L1 -- V(6) or IL in V
        setVCVS(3, 0, 4, 0, alpha/R3); % Set VCVS (V = alpha*I3 = (alpha/R3)*V3 -- V(7) or Ix in V (current through VCVS)

        % Everytime voltage source is setup, G matrix is changed, so circuit needs to be redone (put above)
        setVoltageSource(1, 0, v_in(i)); % V(8) or Iin in V
        setCurrentSource(3, 0, rand_currents(i));

        if i == 1
            V_new = (G + (C./dt)) \ F; % no old voltage, so C/dt * 0 = 0
        else
            V_new = (G + (C./dt)) \ (F + (C./dt)*V_old);
        end
        vout(i) = V_new(5); % V(5) = vout

        V_old = V_new;
    end
    
    sgtitle('Varying Cn Investigation');
    if k == 1
        subplot(4, 2, 1)
        plot(time, v_in, time, vout);
        title('Pulse Response (Default)');
        xlabel('time [s]'); ylabel('Voltage [V]');
        legend('Vin', 'Vout');
        subplot(4, 2, 2)
        plot(fshift, mag2db(abs(fftshift(fft(v_in)))));
        hold on;
        plot(fshift, mag2db(abs(fftshift(fft(vout)))));
        hold on;
        title(['Frequency Spectrum (Default) | Cn = ', num2str(Cn(1))]);
        xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
        legend('Vin', 'Vout');
    elseif k == 2
        subplot(4, 2, 3)
        plot(time, v_in, time, vout);
        title('Pulse Response');
        xlabel('time [s]'); ylabel('Voltage [V]');
        legend('Vin', 'Vout');
        subplot(4, 2, 4)
        plot(fshift, mag2db(abs(fftshift(fft(v_in)))));
        hold on;
        plot(fshift, mag2db(abs(fftshift(fft(vout)))));
        hold on;
        title(['Frequency Spectrum | Cn = ', num2str(Cn(2))]);
        xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
        legend('Vin', 'Vout');
    elseif k == 3
        subplot(4, 2, 5)
        plot(time, v_in, time, vout);
        title('Pulse Response');
        xlabel('time [s]'); ylabel('Voltage [V]');
        legend('Vin', 'Vout');
        subplot(4, 2, 6)
        plot(fshift, mag2db(abs(fftshift(fft(v_in)))));
        hold on;
        plot(fshift, mag2db(abs(fftshift(fft(vout)))));
        hold on;
        title(['Frequency Spectrum | Cn = ', num2str(Cn(3))]);
        xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
        legend('Vin', 'Vout');
    else
        subplot(4, 2, 7)
        plot(time, v_in, time, vout);
        title('Pulse Response');
        xlabel('time [s]'); ylabel('Voltage [V]');
        legend('Vin', 'Vout');
        subplot(4, 2, 8)
        plot(fshift, mag2db(abs(fftshift(fft(v_in)))));
        hold on;
        plot(fshift, mag2db(abs(fftshift(fft(vout)))));
        hold on;
        title(['Frequency Spectrum | Cn = ', num2str(Cn(4))]);
        xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
        legend('Vin', 'Vout');
    end
end

% Vary timestep
numSteps = [1000 100 10000];
dt = simTime./numSteps;
figure
Cn = 0.00001;
for k = 1:length(numSteps)
    vout = zeros(1, numSteps(k));
    
    time = linspace(0, 1, numSteps(k));
    gauss_pulse = pdf(gauss_dist, time); % makes a max amplitude of ~13.3
    v_pulse = (gauss_pulse.*sin(pi*time))/max(gauss_pulse.*sin(pi*time));
    v_in = v_pulse;
    
    rand_currents = In * randn(1, length(v_in));
    
    for i = 1:numSteps(k)
               
        G = zeros(numNodes, numNodes);
        C = zeros(numNodes, numNodes);
        F = zeros(numNodes, 1);

        % Setup G, C, and F matrices through stamp functions
        setResistor(1, 2, R1); % Set R1 -- V(1)
        setResistor(2, 0, R2); % Set R2 -- V(2)
        setResistor(3, 0, R3); % Set R3 -- V(3)
        setResistor(4, 5, R4); % Set R4 -- V(4)
        setResistor(5, 0, Ro); % Set R5 -- V(5)
        setCapacitor(1, 2, C1); % Set C1
        setCapacitor(0, 3, Cn); % Set Cn
        setInductor(2, 3, L1); % Set L1 -- V(6) or IL in V
        setVCVS(3, 0, 4, 0, alpha/R3); % Set VCVS (V = alpha*I3 = (alpha/R3)*V3 -- V(7) or Ix in V (current through VCVS)

        % Everytime voltage source is setup, G matrix is changed, so circuit needs to be redone (put above)
        setVoltageSource(1, 0, v_in(i)); % V(8) or Iin in V
        setCurrentSource(3, 0, rand_currents(i));

        if i == 1
            V_new = (G + (C./dt(k))) \ F; % no old voltage, so C/dt * 0 = 0
        else
            V_new = (G + (C./dt(k))) \ (F + (C./dt(k))*V_old);
        end
        vout(i) = V_new(5); % V(5) = vout

        V_old = V_new;
    end
    
    fshift = (-length(v_in)/2:length(v_in)/2-1);
    
    sgtitle('Changing Timestep Investigation');    
    if k == 1
        subplot(3, 2, 1)
        plot(time, v_in, time, vout);
        title('Pulse Response (Default)');
        xlabel('time [s]'); ylabel('Voltage [V]');
        legend('Vin', 'Vout');
        subplot(3, 2, 2)
        plot(fshift, mag2db(abs(fftshift(fft(v_in)))));
        hold on;
        plot(fshift, mag2db(abs(fftshift(fft(vout)))));
        hold on;
        title(['Frequency Spectrum (Default) | dt = ', num2str(dt(1))]);
        xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
        legend('Vin', 'Vout');
    elseif k == 2
        subplot(3, 2, 3)
        plot(time, v_in, time, vout);
        title('Pulse Response');
        xlabel('time [s]'); ylabel('Voltage [V]');
        legend('Vin', 'Vout');
        subplot(3, 2, 4)
        plot(fshift, mag2db(abs(fftshift(fft(v_in)))));
        hold on;
        plot(fshift, mag2db(abs(fftshift(fft(vout)))));
        hold on;
        title(['Frequency Spectrum | dt = ', num2str(dt(2))]);
        xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
        legend('Vin', 'Vout');
    else
        subplot(3, 2, 5)
        plot(time, v_in, time, vout);
        title('Pulse Response');
        xlabel('time [s]'); ylabel('Voltage [V]');
        legend('Vin', 'Vout');
        subplot(3, 2, 6)
        plot(fshift, mag2db(abs(fftshift(fft(v_in)))));
        hold on;
        plot(fshift, mag2db(abs(fftshift(fft(vout)))));
        hold on;
        title(['Frequency Spectrum | dt = ', num2str(dt(3))]);
        xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
        legend('Vin', 'Vout');
    end
end
