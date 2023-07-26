function expansion = Taylor_expansion(L,Q,Qeq,dq,order)

% Multivariable Taylor series expansion
% expansion: final expansion of the function L
% L: function to be expanded
% Q: variables of expansion
% Qeq: values of the variables around which the expansion is performed
% dq: Q-Qeq
% order: order of the Taylor series

N = length(Q);
expansion = L;

% Zeroth term
for i = 1:N
    expansion = subs(expansion,Q(i),Qeq(i));
end

% Terms with derivatives up to the order
j = 1;
while j<=order
    
    if j==1
        der1 = L;
    end
    
    % Differentiate with respect to the states
    der2 = 0;
    for i = 1:N
        derivative = diff(der1,Q(i));
        der2 = der2+derivative*dq(i)/factorial(j)*factorial(j-1);
    end
    
    % Substitute by the equilibrium values
    integr = der2;
    for i = 1:N
        integr = subs(integr,Q(i),Qeq(i));
    end
    
    expansion = expansion+integr;
    der1 = der2;
    j = j+1;
    
end

% expansion = simplify(expansion);

end