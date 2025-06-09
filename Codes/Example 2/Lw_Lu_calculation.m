function [Lw,Lu] = Lw_Lu_calculation(u_id,y_id,p,f,type)

if nargin<5
    type = 'n'; 
end
if size(u_id,1) > size(u_id,2)
    u_id=u_id'; 
end
if size(y_id,1) > size(y_id,2)
    y_id=y_id';
end

m = size(u_id,1);
l = size(y_id,1);
Nsysid = length(y_id);

if strcmp(type,'v')
    u_id = filter([1 -1],1,u_id);
    y_id = filter([1 -1],1,y_id);
end

%% Generate linear predictors
Up = hankel(u_id(1:p*m),u_id(p*m:(Nsysid-f)*m));
Up = Up(:,1:m:end);
Uf = hankel(u_id(p*m+1:(p+f)*m),u_id((p+f)*m:Nsysid*m));
Uf = Uf(:,1:m:end);
Yp = hankel(y_id(1:p*l),y_id(p*l:(Nsysid-f)*l));
Yp = Yp(:,1:l:end);
Yf = hankel(y_id(p*l+1:(p+f)*l),y_id((p+f)*l:Nsysid*l));
Yf = Yf(:,1:l:end);
Wp = [Up; Yp];
[R,Q] = RQ_Gram([Wp;Uf;Yf]); %#ok

%% Generate linear predictors, contd.
L = R(p*(m+l)+f*m+1:end,1:p*(m+l)+f*m)*pinv(R(1:p*(m+l)+f*m,1:p*(m+l)+f*m));
Lw = L(:,1:p*(m+l));
Lu = L(:,p*(m+l)+1:end);
