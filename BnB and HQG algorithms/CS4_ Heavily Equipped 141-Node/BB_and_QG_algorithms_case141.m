% Programmed by A. Alburidy and L. Fan
% arburidy@gmail.com
% If you find this code useful for your research, please cite our paper at:
% https://github.com/alburidy/ADMM-VVO-Optimization
%==========================
yalmip('clear');
clear; clc; close all;
%--------------------------
CaseName ='case141';
mpc = loadcase(CaseName);    % Enter any matpower radial case
%--------------------------
SD.oltc_l=[1 2;43 44;42 54;23 24];% Enter LTCs locations
SD.cb_l=[8;30;40;48;60;65;70;80;85;90;100;112];       % Enter SCBs locations
%--------------------------
Qc_max=repmat(1200,12,1);  % Enter Qc max. capacity in kVARs
SD.Cstp = 0.03;              % Enter Qc incremental step change value.
SD.Cinr = repmat(4,12,1);         % Enter max. integer number for bsh
%--------------------------
% Enter resulted ADMM total power losses (kWs) "for comparison"
P_loss_ADMM =0.378125712724398*1e3;
%###################################################
% You do NOT need to change anything after this line
%---------------------------------------------------
SD.tpn=size(SD.oltc_l,1);
SD.cbn=size(SD.cb_l,1);
SD.Qc_max=Qc_max/(1e3*mpc.baseMVA);
nn=size(mpc.bus,1);
nbr=size(mpc.branch,1);
Pd=(mpc.bus(2:end,3)/mpc.baseMVA);  % active Load power demands
Qd=(mpc.bus(2:end,4)/mpc.baseMVA);  % reactive Load power demands
Smax=mpc.branch(:,6)./mpc.baseMVA;  % load lines maximum loading level
Smax(Smax==0) = 999;
[Y, Yi, Yj] = makeYbus(mpc);
Gi=real(Yi); Bi=imag(Yi);
%--------------------------
[~, SD.br_oltc]=ismember(SD.oltc_l,mpc.branch(:,1:2),'rows');
SD.id_br_no_t=setdiff([1:size(mpc.branch,1)]',SD.br_oltc,'stable');
br_no_t=setdiff(mpc.branch(:,1:2),SD.oltc_l,'rows','stable');
SD.fr_no_t=br_no_t(:,1);
SD.to_no_t=br_no_t(:,2);
SD.id_no_t1=sub2ind(size(Gi),SD.id_br_no_t, SD.to_no_t);
SD.id_t1=sub2ind(size(Gi),SD.br_oltc,SD.oltc_l(:,2));
SD.id_no_t2=sub2ind(size(Gi),SD.id_br_no_t, SD.fr_no_t);
SD.id_t2=sub2ind(size(Gi),SD.br_oltc,SD.oltc_l(:,1));
%--------------------------
ij=zeros(nbr,nn);
ji=zeros(nbr,nn);
for n=1:nbr %branches
    for m=1:nn %buses
        if m==mpc.branch(n,1)
            ij(n,m)=-1;
        elseif m==mpc.branch(n,2)
            ji(n,m)=-1;
        end
    end
end
SD.ij=ij(:,2:end);
SD.ji=ji(:,2:end);
%####################################################################
% Algorithm starts here
%####################################################################
v     = sdpvar(nn,1);
theta = sdpvar(nn,1);
u     = intvar(SD.tpn,1); %integer variable
u_c   = intvar(SD.cbn,1); %integer variable

assign(v,1);
assign(theta,0);
assign(u,17);
assign(u_c,0);

Pij(SD.id_br_no_t,1)=...
Gi(SD.id_no_t2).*v(SD.fr_no_t).^2+v(SD.fr_no_t).*v(SD.to_no_t).*(Gi(SD.id_no_t1).*cos(theta(SD.fr_no_t)-theta(SD.to_no_t))+...
Bi(SD.id_no_t1).*sin(theta(SD.fr_no_t)-theta(SD.to_no_t)));

Pij(SD.br_oltc,1)=...
Gi(SD.id_t2).*(v(SD.oltc_l(:,1))./(0.89375+u*0.00625)).^2+(v(SD.oltc_l(:,1))./(0.89375+u*0.00625)).*v(SD.oltc_l(:,2)).*(Gi(SD.id_t1).*cos(theta(SD.oltc_l(:,1))-theta(SD.oltc_l(:,2)))+...
Bi(SD.id_t1).*sin(theta(SD.oltc_l(:,1))-theta(SD.oltc_l(:,2))));

Pji(SD.id_br_no_t,1)=...
Gi(SD.id_no_t2).*v(SD.to_no_t).^2+v(SD.fr_no_t).*v(SD.to_no_t).*(Gi(SD.id_no_t1).*cos(theta(SD.to_no_t)-theta(SD.fr_no_t))+...
Bi(SD.id_no_t1).*sin(theta(SD.to_no_t)-theta(SD.fr_no_t)));

Pji(SD.br_oltc,1)=...
Gi(SD.id_t2).*v(SD.oltc_l(:,2)).^2+(v(SD.oltc_l(:,1))./(0.89375+u*0.00625)).*v(SD.oltc_l(:,2)).*(Gi(SD.id_t1).*cos(theta(SD.oltc_l(:,2))-theta(SD.oltc_l(:,1)))+...
Bi(SD.id_t1).*sin(theta(SD.oltc_l(:,2))-theta(SD.oltc_l(:,1))));

Qij(SD.id_br_no_t,1)=...
-Bi(SD.id_no_t2).*v(SD.fr_no_t).^2+v(SD.fr_no_t).*v(SD.to_no_t).*(Gi(SD.id_no_t1).*sin(theta(SD.fr_no_t)-theta(SD.to_no_t))-...
Bi(SD.id_no_t1).*cos(theta(SD.fr_no_t)-theta(SD.to_no_t)));

Qij(SD.br_oltc,1)=...
-Bi(SD.id_t2).*(v(SD.oltc_l(:,1))./(0.89375+u*0.00625)).^2+(v(SD.oltc_l(:,1))./(0.89375+u*0.00625)).*v(SD.oltc_l(:,2)).*(Gi(SD.id_t1).*sin(theta(SD.oltc_l(:,1))-theta(SD.oltc_l(:,2)))-...
Bi(SD.id_t1).*cos(theta(SD.oltc_l(:,1))-theta(SD.oltc_l(:,2))));

Qji(SD.id_br_no_t,1)=...
-Bi(SD.id_no_t2).*v(SD.to_no_t).^2+v(SD.fr_no_t).*v(SD.to_no_t).*(Gi(SD.id_no_t1).*sin(theta(SD.to_no_t)-theta(SD.fr_no_t))-...
Bi(SD.id_no_t1).*cos(theta(SD.to_no_t)-theta(SD.fr_no_t)));

Qji(SD.br_oltc,1)=...
-Bi(SD.id_t2).*v(SD.oltc_l(:,2)).^2+(v(SD.oltc_l(:,1))./(0.89375+u*0.00625)).*v(SD.oltc_l(:,2)).*(Gi(SD.id_t1).*sin(theta(SD.oltc_l(:,2))-theta(SD.oltc_l(:,1)))-...
Bi(SD.id_t1).*cos(theta(SD.oltc_l(:,2))-theta(SD.oltc_l(:,1))));

Ploss(SD.id_br_no_t,1)=...
-Gi(SD.id_no_t1).*(v(SD.fr_no_t).^2+v(SD.to_no_t).^2-2.*v(SD.fr_no_t).*v(SD.to_no_t).*cos(theta(SD.fr_no_t)-theta(SD.to_no_t)));

Ploss(SD.br_oltc,1)=...
-Gi(SD.id_t1).*((v(SD.oltc_l(:,1))./(0.89375+u*0.00625)).^2+v(SD.oltc_l(:,2)).^2-2.*(v(SD.oltc_l(:,1))./(0.89375+u*0.00625)).*v(SD.oltc_l(:,2)).*cos(theta(SD.oltc_l(:,1))-theta(SD.oltc_l(:,2))));
%--------------------------------------   
Pi=SD.ij'*Pij+SD.ji'*Pji;
Qi=SD.ij'*Qij+SD.ji'*Qji;

% Include CBs
Qi(SD.cb_l-1,1)=Qi(SD.cb_l-1,1)+[(u_c.*SD.Cstp).*v(SD.cb_l).^2];
%--------------------------------------

Constraints =[v(1)==1;
              theta(1)==0;
              0.9<=v<=1.1;
              -pi<=theta<=pi;
              %------------
              Pi-Pd==0;
              Qi-Qd==0;
              %------------
              (Pij.^2+Qij.^2)<=Smax.^2;
              %------------
              1<=u<=33;
              0<=u_c<=SD.Cinr];

Objective = sum(Ploss);

total_time=tic;
Sol = optimize(Constraints,Objective,sdpsettings('usex0',1,'solver','knitro','knitro.optionsfile','knitro_options.txt'));
timer=toc(total_time);

flag = Sol.problem;
if flag == 0
    Objective=value(Objective);
    P_loss_MINLP=value(sum(Ploss))*mpc.baseMVA*1e3;
    v=value(v);
    theta=value(theta);
    u=value(u);
    u_c=value(u_c);
    
    %####################################################################
    % Prints the results
    %####################################################################
    disp('Obj. Value | max(V) | min(V) | u ');
    disp('---------------------------------------------------|');
    fprintf(['%8.5f %10.4f %8.4f |',...
        repmat(' %2.1f ',1, SD.tpn+SD.cbn),'\n\n\n']',Objective,max(v),min(v),[u;u_c]');
    disp('---------------------------------------------------|');
    fprintf('Total execution time: %d minutes and %f seconds \n\n\n', floor(timer/60), rem(timer,60));
    
    %Run matpower power flow in the system original status(for comparison)
    matpwr_original=runpf(mpc,mpoption('verbose',0,'out.all',0));
    % Calculating total power losses for both cases:
    P_loss_original = sum(real(get_losses(matpwr_original)))*1e3;
    
    disp('---------------------------------------------------|');
    disp('P_loss ADMM | P_loss MINLP | P_loss Original Status');
    fprintf('%10.4f %14.4f %16.4f \n',P_loss_ADMM,P_loss_MINLP,P_loss_original);
else
    Sol.info
end