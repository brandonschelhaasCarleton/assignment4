function setCurrentSource(node1, node2, val)
% Node1 = input node
% Node2 = output node
% val = Vin value;

% node1 --> source --> node2
% (current flow out of node1 and into node2)

global F;

% Create stamp for current source in F matrix
% Help from resources in main
if node1 ~= 0 
    F(node1) = F(node1) - val;
end
if node2 ~= 0
    F(node2) = F(node2) + val;
end

end


