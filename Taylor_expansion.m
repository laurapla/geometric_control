function exp_avg = Taylor_expansion(L,Q,Qeq,dq,order)

% Multivariable Taylor series expansion
% expansion: final expansion of the function L
% L: function to be expanded
% Q: variables of expansion
% Qeq: values of the variables around which the expansion is performed
% dq: Q-Qeq
% order: order of the Taylor series

N = length(Q);
expansion = L;
syms w t;
syms epsil real;

w = w/epsil;

% Zeroth term
for i = 1:N
    if i~=6 && i~=7
        expansion = simplify(subs(expansion,Q(i),Qeq(i)));
    end
end

exp_avg = simplify(expand(int(expansion,t,0,2*pi/w)*w/(2*pi)));

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
        if i~=6 && i~=7
            integr = simplify(subs(integr,Q(i),Qeq(i)));
        end
    end
    
    int_avg = simplify(expand(int(integr,t,0,2*pi/w)*w/(2*pi)));
    exp_avg = exp_avg+int_avg;
    der1 = simplify(der2);
    j = j+1;
    
end

exp_avg = simplify(exp_avg);

end