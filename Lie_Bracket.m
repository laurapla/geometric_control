function V = Lie_Bracket(Q,f,g)

% V = [f,g] Lie Bracket
% Q: Coordinates which we differentiate with respect to (for example Q =
% [x;y;z] for Cartesian coordinates)

for i = 1:length(f)
    for j = 1:length(f)
        J_f(i,j) = diff(f(i),Q(j));
        J_g(i,j) = diff(g(i),Q(j));
    end
end

V = simplify(J_g*f-J_f*g);

end