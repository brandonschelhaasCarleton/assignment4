function setResistor(node1, node2, val)
% Node1 = input node
% Node2 = output node
% val = resistor value;

global G;
    
% Transform resistor into conductance for G matrix
value = 1/val;

% Create stamp for resistor in G matrix
if node1 ~= 0 && node2 == 0 % node 2 ground
    G(node1, node1) = G(node1, node1) + value;
end
if node1 == 0 && node2 ~= 0 % node 1 ground
    G(node2, node2) = G(node2, node2) + value;
end
if node1 ~= 0 && node2 ~= 0 % from general KCL
    G(node1, node1) = G(node1, node1) + value;
    G(node2, node2) = G(node2, node2) + value;
    G(node1, node2) = G(node1, node2) - value;
    G(node2, node1) = G(node2, node1) - value;
end
end

