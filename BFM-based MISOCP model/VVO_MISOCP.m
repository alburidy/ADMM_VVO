% Programmed by A. Alburidy and L. Fan
% arburidy@gmail.com
% If you find this code useful for your research, please cite our paper at:
% https://github.com/alburidy/ADMM-VVO-Optimization
%==========================
%#ok<*CLCLS,*NOPTS,*FNDSB>
clear; close all; clear classes; clc;
cvx_solver_settings -clearall
%--------------------------
% Loading case study data
fetching_system_data; 
%--------------------------
nbr = size(mpc.branch,1);
nn=size(mpc.bus,1);
r = mpc.branch(:,3);
x = mpc.branch(:,4);
z = r+j*x;
Pd = mpc.bus(2:end,3)/mpc.baseMVA;  % active Load power demands
Qd = mpc.bus(2:end,4)/mpc.baseMVA;  % reactive Load power demands
fr = mpc.branch(:,1);   % Vector of sending end nodes
to = mpc.branch(:,2);   % Vector of receiving end nodes
%---------------------------
% Matrices to recover angles later
C_matrix=zeros(nn,nbr);
for n=1:nbr
    for m=1:nn
        if m==mpc.branch(n,1)
            C_matrix(m,n)=1;
        elseif m==mpc.branch(n,2)
            C_matrix(m,n)=-1;
        end
    end
end
B_matrix=C_matrix(2:end,:);
%---------------------------
%Voltage Limits
Vmin=0.9^2;
Vmax=1.1^2;
%---------------------------
T=1; % Horizons


% Voltage Regulator:
tmin_VR=0.9;
tmax_VR=1.1;
K_VR = 32;
Dt_VR=(tmax_VR-tmin_VR)/K_VR;
Positions_VR = 0:K_VR;
for k=1:length(Positions_VR)
   W(k) = (tmin_VR+(k-1)*Dt_VR)^2; %squared ratio
end


tic
cvx_begin quiet
cvx_solver mosek

variables V(nn,T) 
variable L(nn,T)
variable P(nn,T)
expression Pline(nn, T)
variable Q(nn, T)
expression Qline(nn, T)
expression Ploss(nn)

variable Qcb(length(Nc), T)
variable C(length(Nc), T) integer

variables VR(length(VR_bus), T) 
variable o(K_VR+1,length(VR_bus), T) binary
variable y(K_VR+1,length(VR_bus), T)
variables ratio_VR(length(VR_bus),T)


Ploss=0;
for t= 1: T
Ploss = Ploss + sum(r.*L(to, t));
end



minimize  Ploss;

subject to


for i= 1: T


for k= 1: length(Nc)
    C(k,i)<=Kc(k);
    C(k,i)>=0;
    Qcb(k,i)==Qc(k)*C(k,i)/Kc(k);
end


%Voltage Limits
    for k= 1: nn
           Vmin <= V(k, i) <= Vmax
    end
    for k= 1: length(VR_bus)
           Vmin <= VR(k,i) <= Vmax
    end
 

% OPF
Pline(:,i) = zeros(nbr+1, 1);
Qline(:,i) = zeros(nbr+1, 1);
for k = nbr: -1: 1
    Pline(to(k), i) = Pd(k) + r(k).*L(to(k),i);

    if ismember(to(k),Nc) == 1 
        Qline(to(k), i) = Qd(k) + x(k).*L(to(k),i)- Qcb(find(Nc==to(k)),i);
    else
        Qline(to(k), i) = Qd(k) + x(k).*L(to(k),i);
    end
end


for k = nbr: -1: 1
    Pline(fr(k),i) = Pline(fr(k),i) + Pline(to(k),i); 
    Qline(fr(k),i) = Qline(fr(k),i) + Qline(to(k),i);
end

P(:,i) == Pline(:,i);
Q(:,i) == Qline(:,i); 





%Voltage
for k=1:nbr
    if ismember(k,VR_bus) == 1 && ismember(to(k),VR_bus_to) == 1
        V(to(k), i) == VR(find(VR_bus==k),i)-2*(r(k)*P(to(k),i)+x(k)*Q(to(k),i)) + (r(k)^2+x(k)^2)*L(to(k),i); 
    else
        V(to(k), i) == V(fr(k), i)-2*(r(k)*P(to(k),i)+x(k)*Q(to(k),i)) + (r(k)^2+x(k)^2)*L(to(k),i); 
    end
end

%Current
for k = 1: nbr
    if ismember(k,VR_bus) == 1 && ismember(to(k),VR_bus_to) == 1
       L(to(k),i) + VR(find(VR_bus==k),i) >= norm([2*P(to(k),i); 2*Q(to(k),i); L(to(k),i)-VR(find(VR_bus==k),i)]);  
    else
       L(to(k),i) + V(fr(k),i) >= norm([2*P(to(k),i); 2*Q(to(k),i); L(to(k),i)-V(fr(k),i)]); 
    end
end


V(1,i)==1;

%VR
for k = 1: length(VR_bus)
    VR(k,i) == sum(y(:,k,i).*W');
    V(VR_bus(k),i).*ones(K_VR+1,1)-(ones(K_VR+1,1)-o(:,k,i)).*Vmax <= y(:,k,i);
    V(VR_bus(k),i).*ones(K_VR+1,1)-(ones(K_VR+1,1)-o(:,k,i)).*Vmin >= y(:,k,i);
    o(:,k,i).*Vmin <= y(:,k,i);
    o(:,k,i).*Vmax >= y(:,k,i);
    sum(o(:,k,i)) == 1;
    ratio_VR(k,i) == sum(o(:,k,i).*W')
end

end

cvx_end
toc

% Recoving voltage, current and angles
U_tap = (sqrt(V(VR_bus))./sqrt(V(VR_bus_to))-0.89375)./0.00625;
Sline=Pline+j*Qline;
I_line=zeros(nbr,1);
v_ph=zeros(nn,1);
v_ph(1,1)=1*exp(j*0);

clear j
for k = 1: nbr
    if ismember(k,VR_bus) == 1 && ismember(to(k),VR_bus_to) == 1
        beta(k,1) = angle(VR(find(VR_bus==k))-((r(k)-j*x(k)))*(Pline(k+1)+j*Qline(k+1)));
        I_line(k,1)=(1/VR(find(VR_bus==k)))*conj(Sline(to(k)))*v_ph(fr(k));
        v_ph(to(k))=v_ph(fr(k))-z(k)*I_line(k);
    else
        beta(k,1) = angle(V(fr(k))-((r(k)-j*x(k)))*(Pline(k+1)+j*Qline(k+1)));
        I_line(k,1)=(1/V(fr(k)))*conj(Sline(to(k)))*v_ph(fr(k));
        v_ph(to(k))=v_ph(fr(k))-z(k)*I_line(k);
    end
end



theta_MISOCP= [beta'*inv(B_matrix)]';


V_mag = sqrt(V);
I_meg = sqrt(L);
% %==================================
% Run MATPOWER for comparieson
mpc.bus(Nc',6)=(C.*Cstp)*mpc.baseMVA; % Feeding the resulted CB settings to matpower
mpc.branch(br_oltc,9)=(0.89375+U_tap.*0.00625); % Feeding the resulted VRs and LTCs settings to matpower
% Run MATPOWER
mtpwr = runpf(mpc,mpoption('verbose',0,'out.all',0));
if mtpwr.success==1
% Computing matpower's current
[Y, Yi, Yj] = makeYbus(mpc);
matpower_V_mag=mtpwr.bus(:,8);
matpower_theta=deg2rad(mtpwr.bus(:,9));
% Building nodal voltage phasors
V_phsor=matpower_V_mag.*exp(1i*matpower_theta);
matpower_I_mag=abs(Yi*V_phsor);

Pline_matpower=mtpwr.branch(:,14);
Qline_matpower=mtpwr.branch(:,15);

theta_MISOCP=[0;theta_MISOCP];

Diff_volt=abs((V_mag-matpower_V_mag)./matpower_V_mag).*100;
Diff_current=abs((I_meg(2:end)-matpower_I_mag)./matpower_I_mag).*100;
Diff_theta=abs((theta_MISOCP-matpower_theta)./matpower_theta).*100;

if max(matpower_V_mag)<=1.1 && min(matpower_V_mag)>=0.9
figure(3)
stairs(V_mag,'LineWidth',1);
xlim([1 nn])
grid on
hold on 
stairs(matpower_V_mag,'LineWidth',1);
hold off
legend('MISOCP','Matpower');
title('Nodal Voltage Magnitudes');
xlabel('Nodes');
ylabel('Nodal Voltage Magnitude (pu)');


figure(4)
stairs(theta_MISOCP,'LineWidth',1);
hold on 
xlim([1 nn])
grid on
stairs(matpower_theta,'LineWidth',1);
hold off
legend('MISOCP','Matpower');
title('Nodal Voltage Angles');



figure(5)
stairs(Diff_volt,'LineWidth',1);
xlim([1 nn])
grid on
title('Voltage Magnitude Deviation Ratio');
xlabel('Nodes');
ylabel('Deviation Ratio (%)');


figure(6)
stairs(Diff_theta,'LineWidth',1);
xlim([1 nn]);
grid on
title('Angle values Deviation Ratio');
xlabel('Nodes');
ylabel('Deviation Ratio (%)');


figure(7)
bar([Pline_matpower Pline(2:end)*mpc.baseMVA])
legend('Matpower','MISOCP')
title('Active power flow');


figure(8)
bar([Qline_matpower Qline(2:end)*mpc.baseMVA])
legend('Matpower','MISOCP')
title('Reactive Power Flow');
xlabel('Branchs');
ylabel('Reactive Power Flow (kVARs)');


figure(9)
bar([I_meg(2:end) matpower_I_mag])
grid on
legend('MISOCP','Matpower');
title('Branches Current flow');

figure(10)
polarplot(angle(I_line),abs(I_line),'LineWidth',1)
hold on
polarplot(angle(Yi*V_phsor),matpower_I_mag,'LineWidth',1)
hold off
legend('MISOCP','Matpower')
title('Branch Current Flow (pu)')

Total_power_loss_relaxed_BFM = cvx_optval*1e3*mpc.baseMVA
Total_power_loss_MATPOWER = sum(real(get_losses(mtpwr)))*1e3
Voltage_diff = max(abs(V_mag-matpower_V_mag))
else
    disp('These settings violate on or more of the constraints!');
end
else
    disp('MATPOWER Numerically Failed');
end