function y = simulate_plant_ss(u,N,ss_info)
    
    Ap = ss_info{1};
    Bp = ss_info{2};
    Cp = ss_info{3};
    Dp = ss_info{4};
    
    x(:,1) = [0 0]';
    
    y = zeros(N,1);
    
    for i = 1:N
        
        x(:,i+1) = Ap*x(:,i) + Bp*u(i);
        
        y(i) = Cp*x(:,i) + Dp*u(i);
        
    end

end