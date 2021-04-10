% 3bi) DC Sweep
vin = -10:0.1:10;
% vin = vin';
v3 = zeros(length(vin), 1);
vout = zeros(length(vin), 1);
for i = 1:length(vin)
    % Need to reset G, C, F matrices every time because of stamping
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
    setVoltageSource(1, 0, vin(i)); % V(8) or Iin in V

    % Solve for unknowns
    V = G\F;

    % V contains [v1 v2 v3 v4 vout IL Ix Iin]
    vout(i) = V(5); % V(5) = vout
    v3(i) = V(3); % V(3) = v3
end

% Plot Vout and V3 vs. Vin
figure
plot(vin, vout, vin, v3);
title('DC Simulation : Vout vs. Vin');
xlabel('Vin [V]'); ylabel('Vout [V]');
legend('Vout', 'V3')

% 3bii) AC Sweep
freqs = 0:0.1:50; % 1 to 50 Hz
omegas = (2*pi).*freqs;
vin = 1; % Set vin to 1 for frequency domain
vout = zeros(length(freqs),1);
gain = zeros(length(freqs),1);

% Changed vin, so reset the circuit and matrices
G = zeros(numNodes, numNodes);
C = zeros(numNodes, numNodes);
F = zeros(numNodes, 1);
setResistor(1, 2, R1);
setResistor(2, 0, R2);
setResistor(3, 0, R3); 
setResistor(4, 5, R4); 
setResistor(5, 0, Ro);
setCapacitor(1, 2, C1);
setInductor(2, 3, L1);
setVCVS(3, 0, 4, 0, alpha/R3);
setVoltageSource(1, 0, vin);

% Actual AC Sweep
for i = 1:length(freqs)
    V = (G + j*omegas(i).*C)\F;
    vout(i) = V(5);
end
gain = 20.*log(vout./vin);

% Plot AC Sweep
figure
plot(omegas, gain);
title('Gain(w) = Vo/V1');
xlabel('w [rads/s]'); ylabel('Gain [dB]');
xlim([0 max(omegas)]);

% 3biii) Random perturbations on C
randomNums = 0.05*randn(500, 1); % random dist with std = 0.05
randomNums = randomNums + 1; % shift the distribution so there's no negatives
gain = zeros(length(randomNums), 1);
omega = pi;

% Solve for random perturbations
for i = 1:length(randomNums)
    V = (G + j*omega.*C.*randomNums(i))\F; % add random variation to C
    gain(i) = 20.*log(abs(V(5))); %find magnitude of complex number
end

% Plot the histogram of perturbations
figure
subplot(1,2,1)
histogram(C1*randomNums);
title('Capacitance Values with Perturbations');
xlabel('Capacitance [F]'); ylabel('Number');
subplot(1,2,2)
histogram(gain);
title('Gain(w) = Vo/V1 With C Perturbations');
xlabel('Gain [dB]'); ylabel('Number');