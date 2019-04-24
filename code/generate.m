function [Vout] = generate(G, C, b, Vin, In, Nrpt, delta_t, OutputNode)

Vold = zeros(length(b),1);
Vout = zeros(length(b),1);

for i=1:Nrpt
    % Overwrite b matrix with noisy Vin
    b(7) = Vin(i);
    
    % Overwrite b matrix with noisy In
    b(3) = -In(i);
    b(4) = In(i);
    
    A = C/delta_t + G;
    B = b + C*Vold/delta_t;
    x = A\B;
    
    Vout(i) = x(OutputNode);
    
    Vold = x;
end

end

