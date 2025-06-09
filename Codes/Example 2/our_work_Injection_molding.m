function out = our_work_Injection_molding()

%% Weighting functions
s = tf('s');
Wd_tf = 0.1*(0.01/5)*((s+5)/(s+0.01));
Wr_tf = 20*(0.01/5)*((s+5)/(s+0.01));
Wu_tf = 0.5*(10/0.5)*((s+0.01)/(s+5));

freq = 200;   % (default: 200)
T = 1/freq;

Wd_dis = c2d(Wd_tf,T,'tustin');
Wr_dis = c2d(Wr_tf,T,'tustin');
Wu_dis = c2d(Wu_tf,T,'tustin');

Wd_ss = ss(Wd_dis);
Awd = Wd_ss.A;
Bwd = Wd_ss.B;
Cwd = Wd_ss.C;
Dwd = Wd_ss.D;

Wr_ss = ss(Wr_dis);
Awr = Wr_ss.A;
Bwr = Wr_ss.B;
Cwr = Wr_ss.C;
Dwr = Wr_ss.D;

Wu_ss = ss(Wu_dis);
Awu = Wu_ss.A;
Bwu = Wu_ss.B;
Cwu = Wu_ss.C;
Dwu = Wu_ss.D;

%% Plant (Injection molding process)
Ap = [ 1.582 -0.592
         1      0   ];

Bp_u = [ 1
         0 ];
      
Bp_dist = [    0
            -54.26 ];      

Cp = [1.69 1.419];

Dp = 0;

sys_dis_ss = ss(Ap,Bp_u,Cp,Dp);

ss_info = {sys_dis_ss.A sys_dis_ss.B sys_dis_ss.C sys_dis_ss.D};

sys_dis_tf = tf(sys_dis_ss);
Numerator = sys_dis_tf.Numerator{1};
Numerator = Numerator(2:end);
Denominator = sys_dis_tf.Denominator{1};
Denominator = Denominator(2:end);

%% Past I/O data preparation from plant model
N_data = 2000;
u_id = idinput(N_data); % PRBS input
% u_id = ones(N_data,1);
% y_id = simulate_plant_tf(u_id,N_data,Numerator,Denominator);   % plant output
y_id = simulate_plant_ss(u_id,N_data,ss_info);   % plant output

alpha = 0;
snr_dB = 50;  % (default: 50)
[integrated_noisy_y_id,scaled_integrated_noise,noise] = add_noise_func(y_id,alpha,snr_dB);

% figure(2);
% subplot(3,1,1);
% plot(y_id,'r--');
% hold on;
% plot(integrated_noisy_y_id,'b');
% grid on;
% ylabel('\bfy');
% legend('\bfHealthy output','\bfNoisy output','Orientation','Horizontal')
% 
% subplot(3,1,2);
% plot(scaled_integrated_noise,'b--');
% grid on;
% ylabel('\bf\xi');
% legend('\bfScaled integrated noise');
% 
% subplot(3,1,3);
% plot(u_id,'b');
% ylim([-1.5,1.5]);
% grid on;
% ylabel('\bfu');
% xlabel('\bfSampling time');
% legend('\bfPRBS input');

%% Simulation parameters
ts = 0;
dt = T;
t_f = 50;  % (default: 50)
nSteps = floor(t_f/dt);  % steps of simulation loop
time = (ts:nSteps-1)*dt;    % time vector

p = 10; % order of past data (default: 10)
f = p;  % order of future data (default: 10)

gamma_min = 1.550871;   % (default: 0.301217)
gamma = (1.1)*gamma_min;    % upper bound of cost function   

%% Lw and Lu calculation
[Lw,Lu] = Lw_Lu_calculation(u_id,integrated_noisy_y_id,p,f,'v');

%% Toeplitz matrices (markov parameters)
Hd = Dwd*eye(f);
Hr = Dwr*eye(f);
Hu = Dwu*eye(f);

for j1 = 1:f
    for i1 = 1:f-1
        if(i1+j1)>f
            break;
        end
        Hd(i1+j1,j1) = Cwd*(Awd^(i1-1))*Bwd;
        Hr(i1+j1,j1) = Cwr*(Awr^(i1-1))*Bwr;
        Hu(i1+j1,j1) = Cwu*(Awu^(i1-1))*Bwu;
    end
end

%% Extended observability matrices
GAMMA_d = zeros(f,1);
GAMMA_r = zeros(f,1);
GAMMA_u = zeros(f,1);
for i2 = 1:f
    GAMMA_d(i2) = Cwd*(Awd^(i2-1));
    GAMMA_r(i2) = Cwr*(Awr^(i2-1));
    GAMMA_u(i2) = Cwu*(Awu^(i2-1));
end

%% Constant matrices
F_l = ones(f,1);
GAMMA_L = tril(ones(f,f));

K1 = [1 0];
K2 = [0 1];
Kref =  repmat(K1,f,1);
Kd =  repmat(K2,f,1);
Ku = eye(f,f);
Hw = Hr*Kref - Hr*Hd*Kd;     
   
OMEGA_cell = {   Hw        -Hr*GAMMA_L*Lu*Ku    -Hr*GAMMA_L*Lw     -Hr*GAMMA_d      GAMMA_r     zeros(f,1)
               zeros(f,2)           Hu            zeros(f,2*p)       zeros(f,1)     zeros(f,1)    GAMMA_u  };

OMEGA_cell_tr = cellfun(@transpose,OMEGA_cell,'UniformOutput',false);
OMEGA_cell_tr = OMEGA_cell_tr';

[m,n] = size(OMEGA_cell);
M_cell = cell(n,n);
for i = 1:n
    for j = 1:n        
        M_cell{i,j} = OMEGA_cell_tr{i,1}*OMEGA_cell{1,j} + OMEGA_cell_tr{i,2}*OMEGA_cell{2,j};        
    end
end

MM = M_cell{1,1} - M_cell{1,2}*(M_cell{2,2}^-1)*M_cell{2,1};
gamma_min = sqrt(abs(max(eig(MM))));

fprintf('*********************************************************\n');
fprintf('*\t\t\t\t\t Our work:\t\t\t\t\t\t\t*\n');
fprintf('*\t\t\t\t gamma_min = %f\t\t\t\t\t*\n',gamma_min);
fprintf('*\t\t\tIf you have changed parameters,\t\t\t\t*\n');
fprintf('*\t\tcopy above value to line 113 and run again!\t\t*\n');
fprintf('*********************************************************\n\n');

M_cell{1,1} = M_cell{1,1} - (gamma^2)*eye(size(M_cell{1,1}));

sai_cell = { -Hr*F_l
            zeros(f,1)};

N_cell = cell(n,1);        
for i = 1:n
    N_cell{i,1} = OMEGA_cell_tr{i,1}*sai_cell{1,1} + OMEGA_cell_tr{i,2}*sai_cell{2,1};
end

%% Controller constants
M_11 = M_cell{1,1};
M_12 = M_cell{1,2};
M_13 = M_cell{1,3};
M_14 = M_cell{1,4};
M_15 = M_cell{1,5};
M_16 = M_cell{1,6};
M_21 = M_cell{2,1};
M_22 = M_cell{2,2};
M_23 = M_cell{2,3};
M_24 = M_cell{2,4};
M_25 = M_cell{2,5};
M_26 = M_cell{2,6};

N_11 = N_cell{1,1};
N_21 = N_cell{2,1};

coeff = -((M_22 - M_21*(M_11^-1)*M_12)^-1)*[-M_21*(M_11^-1) eye(f)];

pr_coeff = coeff*[M_13 M_14 M_15 M_16
                  M_23 M_24 M_25 M_26];
              
add_coeff = coeff*[N_11
                   N_21];                          
               
%% Initial values
up_past = zeros(2,1);
yp_past = zeros(2,1);
% up_past = u_id(end-2+1:end,1);
% yp_past = u_id(end-2+1:end,1);

u_past_prev = zeros(p,1);
% u_past_prev = u_id(end-p+1:end,1);
u_past_new = u_past_prev;

y_past_prev = zeros(p,1);
% y_past_prev = y_id(end-p+1:end,1);
y_past_new = y_past_prev;

wp_prev = [u_past_prev;y_past_prev];

uf_prev = zeros(f,1);
uf_new = uf_prev;

y_t = 0;

uf = zeros(nSteps,1);
yf = zeros(nSteps+1,1);
 
yp = zeros(nSteps+1,1);

x_wd = zeros(nSteps+1,1);
x_wr = zeros(nSteps+1,1);
x_wu = zeros(nSteps+1,1);

xp = zeros(2,nSteps+2);

delta_uf_prev = zeros(f,1);
delta_uf_new = delta_uf_prev;

setPoint = zeros(nSteps+1,1);

% dist = zeros(nSteps+1,1);

% noise = randn(nSteps+1,1);
noise = filter([1 -alpha],[1 -1],noise);

dist_mat = load('dist.mat');
dist = dist_mat.ans(2,:);

%% Simulation
for k = 1:nSteps

%    setPoint(k) = 0.6 + 0.2*(sin(2*pi*k/500)+sin(2*pi*k/1000)+sin(2*pi*k/1500));
   setPoint(k) = (50*(-1)^(round(k/2000)))*(time(k)<=20) + (50*sin(k*pi/1000)+20*cos(k*pi/500))*(time(k)>20 & time(k)<=30) + (50*(-1)^(round(k/2000)))*(time(k)>30)+70;
%    setPoint(k) = (100*(-1)^(round(k/5000)));
   
%    dist(k+1) = 0.001*(time(k)>=35 && time(k)<=45);
%    dist(k+1) = 0.001*(time(k)>=30 && time(k)<=40);
    
   wp_new = [u_past_new;y_past_new];

   delta_wp = wp_new - wp_prev;
    
   uf_prev = uf_new;
   
   delta_uf_prev = delta_uf_new;
   
   delta_uf_new = 1*(pr_coeff*[delta_wp' x_wd(k) x_wr(k) x_wu(k)]' + add_coeff*y_t);
   
   uf_new = uf_prev + delta_uf_new;
%    uf_new = 1;
   uf(k+1) = uf_new(1);
   
   y_hat_f = F_l*y_t + GAMMA_L*Lw*delta_wp + GAMMA_L*Lu*delta_uf_new;
   yf(k+1) = y_hat_f(1);
   
   % Plant simulation (for validating subspace predictor)
%    yp(k+1) = -Denominator*yp_past + Numerator*up_past;

    Ap_z = ss_info{1};
    Bp_z = ss_info{2};
    Cp_z = ss_info{3};
    Dp_z = ss_info{4};
    
    xp(:,k+2) = Ap_z*xp(:,k+1) + Bp_z*(uf(k+1)+dist(k+1));
    
    yp(k+1) = Cp_z*xp(:,k+1) + Dp_z*uf(k+1);
   
%    noise_scale = snr_scale(yp(k+1),noise(k+1),snr_dB);
%    yp(k+1) = yp(k+1) + 0*noise_scale.*noise(k+1);
   
   up_past = [uf(k+1);up_past(1:length(up_past)-1)];
   yp_past = [yp(k+1);yp_past(1:length(yp_past)-1)];
      
   y_t = yp(k+1);
   
   % Updating states of weighting functions
   x_wd(k+1) = Awd*x_wd(k) + Bwd*(dist(k));
   x_wr(k+1) = Awr*x_wr(k) + Bwr*(setPoint(k) - yp(k));
   x_wu(k+1) = Awu*x_wu(k) + Bwu*delta_uf_prev(1);
   
   u_past_prev = u_past_new;
   y_past_prev = y_past_new;
   wp_prev = [u_past_prev;y_past_prev]; 
   
   u_past_new = [uf(k+1);u_past_new(1:length(u_past_new)-1)];
   y_past_new = [yp(k+1);y_past_new(1:length(y_past_new)-1)];
    
end

setPoint = setPoint(1:end-1);
dist = dist(1:end-1);
uf = uf(1:end-1);
yf = yf(1:end-1);
yp = yp(1:end-1);

out.time = time;
out.uf = uf;
out.yp = yp;
out.setPoint = setPoint;

end