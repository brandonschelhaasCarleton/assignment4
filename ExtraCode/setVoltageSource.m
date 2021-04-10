function setVoltageSource(node1, node2, val)
% Node1 = input node
% Node2 = output node
% val = Vin value;

global G C F;

% Have to add new index for Iin
newIndex = size(C,1) + 1;
G(newIndex, newIndex) = 0;
C(newIndex, newIndex) = 0;
F(newIndex) = 0;

% Create stamp for voltage source in G matrix
% Help from resources in main (+/- for voltage nodes)
if node1 ~= 0 
    G(newIndex, node1) = G(newIndex, node1) + 1;
    G(node1, newIndex) = G(node1, newIndex) + 1;
end
if node2 ~= 0
    G(newIndex, node2) = G(newIndex, node2) - 1;
    G(node2, newIndex) = G(node2, newIndex) - 1;
end

% Vin is a source ==> put value in the newIndex
F(newIndex) = val;
end

