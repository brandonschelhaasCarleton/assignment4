function setCapacitor(node1, node2, val)
% Node1 = input node
% Node2 = output node
% val = capacitor value;

global C;

% Create stamp for capacitor in C matrix
if node1 ~= 0 && node2 == 0 % node 2 ground
    C(node1, node1) = C(node1, node1) + val;
end
if node1 == 0 && node2 ~= 0 % node 1 ground
    C(node2, node2) = C(node2, node2) + val;
end
if node1 ~= 0 && node2 ~= 0 % KCL
    C(node1, node1) = C(node1, node1) + val;
    C(node2, node2) = C(node2, node2) + val;
    C(node1, node2) = C(node1, node2) - val;
    C(node2, node1) = C(node2, node1) - val;
end
end

