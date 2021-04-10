% Transient simulation (time domain)

vout = zeros(1, numSteps);
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
    setInductor(2, 3, L1); % Set L1 -- V(6) or IL in V
    setVCVS(3, 0, 4, 0, alpha/R3); % Set VCVS (V = alpha*I3 = (alpha/R3)*V3 -- V(7) or Ix in V (current through VCVS)

    % Everytime voltage source is setup, G matrix is changed, so circuit needs to be redone (put above)
    setVoltageSource(1, 0, v_in(i)); % V(8) or Iin in V
    
    % Solve for unknowns
    if i == 1
        V_new = (G + (C./dt)) \ F; % no old voltage, so C/dt * 0 = 0
    else
        V_new = (G + (C./dt)) \ (F + (C./dt)*V_old);
    end
    
    % Take vout from unknowns
    vout(i) = V_new(5); % V(5) = vout
    
    % Set v_old to v_new for next timestep
    V_old = V_new;
end

% Set the frequency x-axis
fshift = (-length(v_in)/2:length(v_in)/2-1);

% Plot time domain and frequency domain responses
figure
subplot(1, 2, 1)
plot(time, v_in, time, vout);
if v_in == v_step
    title('Step Response');
elseif v_in == v_sin
    title('Sine Response');
else
    title('Pulse Response');
end
xlabel('time [s]'); ylabel('Voltage [V]');
xlim([min(time) max(time)]);
legend('Vin', 'Vout');
subplot(1, 2, 2)
plot(fshift, mag2db(abs(fftshift(fft(v_in)))));
hold on;
plot(fshift, mag2db(abs(fftshift(fft(vout)))));
hold on;
title('Frequency Spectrum');
xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
legend('Vin', 'Vout');

