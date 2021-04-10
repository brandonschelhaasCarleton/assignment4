function setInductor(node1, node2, val)
% Node1 = input node
% Node2 = output node
% val = inductor value;

global G C F;

% Have to add new index to accommodate inductor current
newIndex = size(C,1) + 1;
G(newIndex, newIndex) = 0;
C(newIndex, newIndex) = 0;
F(newIndex) = 0;

% Create stamp for inductor
% Help from resources in main (+/- for voltage nodes)
if node1 ~= 0
    G(newIndex, node1) = G(newIndex, node1) + 1;
    G(node1, newIndex) = G(node1, newIndex) + 1;
end
if node2 ~= 0
    G(newIndex, node2) = G(newIndex, node2) - 1;
    G(node2, newIndex) = G(node2, newIndex) - 1;
end

% -sL => -L in C
C(newIndex, newIndex) = C(newIndex, newIndex) - val;
end

