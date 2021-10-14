function Lvf = Lie_Derivative(Q,V,f)

% V = L_{V}{f} : Lie derivative of f along V
% f: function (scalar)
% V: vector
% Q: Coordinates which we differentiate with respect to (for example Q =
% [x;y;z] for Cartesian coordinates)

Lvf = 0;
for i = 1:length(V)
    J = diff(f,Q(i));
    Lvf = Lvf+V(i)*J;
end

Lvf = simplify(Lvf);

end