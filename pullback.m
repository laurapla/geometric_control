function F = pullback(Q,f,G,U,variable,order)

% F = pullback of the vector field f along the flow phi_t^G of the
% time-varying vector field G
% Q: Coordinates which we differentiate with respect to (for example Q =
% [x;y;z] for Cartesian coordinates) / states of the system
% f = drift vector
% G = matrix with the control vectors in columns
% U = vector with the inputs in columns
% variable = integration variable
% order = order of the approximation

F = f;
n = size(G,2); % number of control vectors and inputs

syms s;
U = subs(U,variable,s);

for i = 1:order
    for j = 1:n
        if i==1
            F = F+Lie_Bracket(Q,G(:,j),f)*int(U(:,j),s,0,variable);
        elseif i==2
            sum2 = 0;
            for l = 1:n
                syms s2;
                Uj = subs(U(:,j),s,s2);
                sum2 = sum2+Lie_Bracket(Q,G(:,j),Lie_Bracket(Q,G(:,l),f))*int(Uj*int(U(:,l),s,0,s2),s2,0,variable);
            end
            F = simplify(F+sum2);
        end
    end
end

end