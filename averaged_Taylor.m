function dPsi = averaged_Taylor(Psi,Q,A,phi)

N = size(Q,1); % Number of states

dPsi = 0;
for i = 1:N
    der1 = diff(Psi,Q(i));
    for q = 1:N
        der2 = diff(der1,Q(q));
        dPsi = dPsi+der2*A(i)*A(q)*cos(phi(i)-phi(q))/2;
    end
end

dPsi = simplify(dPsi);

end