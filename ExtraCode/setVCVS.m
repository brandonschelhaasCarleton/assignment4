function setVCVS(node1_in, node2_in, node1_out, node2_out, val)
% Node1_in = input node for controlling branch
% Node2_in = output node for controlling branch
% Node1_out = input node for dependent branch
% Node2_out = output node for dependent branch
% val = coefficient value;

global G C F;

% Add new index for Ix => the current flowing into VCVS
newIndex = size(C,1) + 1;
G(newIndex, newIndex) = 0;
C(newIndex, newIndex) = 0;
F(newIndex) = 0;

% Create stamp for resistor in C matrix
% Help from resources in main
if node1_in ~= 0
    G(newIndex, node1_in) = G(newIndex, node1_in) - val;
end
if node2_in ~= 0
    G(newIndex, node2_in) = G(newIndex, node2_in) + val;
end
if node1_out ~= 0
    G(newIndex, node1_out) = 1;
    G(node1_out, newIndex) = 1;
end
if node2_out ~= 0
    G(newIndex, node2_out) = -1;
    G(node2_out, newIndex) = -1;
end
end

