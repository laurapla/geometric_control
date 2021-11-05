function F = pullback(Q,f,G,var,N)

% F = pullback of the vector field f along the flow phi_t^G of the
% time-varying vector field G
% Q: Coordinates which we differentiate with respect to (for example Q =
% [x;y;z] for Cartesian coordinates) / states of the system
% f = drift vector
% G = matrix with the control vectors in columns
% var = integration variable
% N = order of the approximation

F = f; % Initialize the pullback with the first term

s = sym('s',[N 1]);
s = [var; s]; % Integration variables

% First term of the Gpullback
V1 = Lie_Bracket(Q,subs(G,var,s(2)),f);
Vint1 = int(V1,s(2),0,s(1));
F = F+simplify(Vint1); % add the term to the pullback

% Second to N terms of the Gpullback
i = 2;
while i<=N
    
    % Compute the Lie Bracket
    V2 = Lie_Bracket(Q,subs(G,var,s(i+1)),V1);
    Vint2 = V2;
    
    % Integrate the Lie Bracket
    for j = i:-1:1
        Vint2 = int(Vint2,s(j+1),0,s(j));
    end
    
    F = F+simplify(Vint2);
    V1 = V2;
    i = i+1;
    
end

F = simplify(F);

end