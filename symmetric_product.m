function sym = symmetric_product(Q,f,g1,g2)

% sym = [g1,[f,g2]] symmetric product
% Q: Coordinates which we differentiate with respect to (for example Q =
% [x;y;z] for Cartesian coordinates)

V = Lie_Bracket(Q,f,g2);
sym = Lie_Bracket(Q,g1,V);

end