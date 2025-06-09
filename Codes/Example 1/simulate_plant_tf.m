function y = simulate_plant_tf(u,N,Numerator,Denominator)

u_prev = zeros(2,1);
y_prev = zeros(2,1);

u_new = zeros(N,1);
y_new = zeros(N,1);

for i = 1:N
    
    u_new(i) = u(i);
        
    u_prev = [u_new(i);u_prev(1:length(u_prev)-1)];
    y_prev = [y_new(i);y_prev(1:length(y_prev)-1)];
    
    y_new(i) = -Denominator*y_prev + Numerator*u_prev;
    
end

% y = y_new;
y = (180/pi)*y_new + 5;

end