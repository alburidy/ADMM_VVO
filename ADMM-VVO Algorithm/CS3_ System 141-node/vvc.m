% Programmed by A. Alburidy and L. Fan
% arburidy@gmail.com
% If you find this code useful for your research, please cite our paper at:
% https://github.com/alburidy/ADMM-VVO-Optimization
%==========================
function [obj_x,P_loss,vm,Qc,v,theta,flag_x,Time_x]=vvc(SD,itr,nn,rho,u,u_c,lambda,Gi,Bi,Pd,Qd,Smax,beta)

v=sdpvar(nn,1);
theta=sdpvar(nn,1);
vm=sdpvar(SD.tpn,1);
Qc=sdpvar(SD.cbn,1);

if itr<1
    assign(v,1);
    assign(vm,1);
    assign(theta,0);
    assign(Qc,0);
end

Pij(SD.id_br_no_t,1)=...
Gi(SD.id_no_t2).*v(SD.fr_no_t).^2+v(SD.fr_no_t).*v(SD.to_no_t).*(Gi(SD.id_no_t1).*cos(theta(SD.fr_no_t)-theta(SD.to_no_t))+...
Bi(SD.id_no_t1).*sin(theta(SD.fr_no_t)-theta(SD.to_no_t)));

Pij(SD.br_oltc,1)=...
Gi(SD.id_t2).*vm.^2+vm.*v(SD.oltc_l(:,2)).*(Gi(SD.id_t1).*cos(theta(SD.oltc_l(:,1))-theta(SD.oltc_l(:,2)))+...
Bi(SD.id_t1).*sin(theta(SD.oltc_l(:,1))-theta(SD.oltc_l(:,2))));

Pji(SD.id_br_no_t,1)=...
Gi(SD.id_no_t2).*v(SD.to_no_t).^2+v(SD.fr_no_t).*v(SD.to_no_t).*(Gi(SD.id_no_t1).*cos(theta(SD.to_no_t)-theta(SD.fr_no_t))+...
Bi(SD.id_no_t1).*sin(theta(SD.to_no_t)-theta(SD.fr_no_t)));

Pji(SD.br_oltc,1)=...
Gi(SD.id_t2).*v(SD.oltc_l(:,2)).^2+vm.*v(SD.oltc_l(:,2)).*(Gi(SD.id_t1).*cos(theta(SD.oltc_l(:,2))-theta(SD.oltc_l(:,1)))+...
Bi(SD.id_t1).*sin(theta(SD.oltc_l(:,2))-theta(SD.oltc_l(:,1))));

Qij(SD.id_br_no_t,1)=...
-Bi(SD.id_no_t2).*v(SD.fr_no_t).^2+v(SD.fr_no_t).*v(SD.to_no_t).*(Gi(SD.id_no_t1).*sin(theta(SD.fr_no_t)-theta(SD.to_no_t))-...
Bi(SD.id_no_t1).*cos(theta(SD.fr_no_t)-theta(SD.to_no_t)));

Qij(SD.br_oltc,1)=...
-Bi(SD.id_t2).*vm.^2+vm.*v(SD.oltc_l(:,2)).*(Gi(SD.id_t1).*sin(theta(SD.oltc_l(:,1))-theta(SD.oltc_l(:,2)))-...
Bi(SD.id_t1).*cos(theta(SD.oltc_l(:,1))-theta(SD.oltc_l(:,2))));

Qji(SD.id_br_no_t,1)=...
-Bi(SD.id_no_t2).*v(SD.to_no_t).^2+v(SD.fr_no_t).*v(SD.to_no_t).*(Gi(SD.id_no_t1).*sin(theta(SD.to_no_t)-theta(SD.fr_no_t))-...
Bi(SD.id_no_t1).*cos(theta(SD.to_no_t)-theta(SD.fr_no_t)));

Qji(SD.br_oltc,1)=...
-Bi(SD.id_t2).*v(SD.oltc_l(:,2)).^2+vm.*v(SD.oltc_l(:,2)).*(Gi(SD.id_t1).*sin(theta(SD.oltc_l(:,2))-theta(SD.oltc_l(:,1)))-...
Bi(SD.id_t1).*cos(theta(SD.oltc_l(:,2))-theta(SD.oltc_l(:,1))));


Ploss(SD.id_br_no_t,1)=...
-Gi(SD.id_no_t1).*(v(SD.fr_no_t).^2+v(SD.to_no_t).^2-2.*v(SD.fr_no_t).*v(SD.to_no_t).*cos(theta(SD.fr_no_t)-theta(SD.to_no_t)));

Ploss(SD.br_oltc,1)=...
-Gi(SD.id_t1).*(vm.^2+v(SD.oltc_l(:,2)).^2-2.*vm.*v(SD.oltc_l(:,2)).*cos(theta(SD.oltc_l(:,1))-theta(SD.oltc_l(:,2))));
   
%--------------------------------------   
Pi=SD.ij'*Pij+SD.ji'*Pji;
Qi=SD.ij'*Qij+SD.ji'*Qji;

% Include CBs
Qi(SD.cb_l-1,1)=Qi(SD.cb_l-1,1)+Qc; % i.e. to install CB on node 8, put 7
%--------------------------------------

Constraints_x =[v(1)==1;
                theta(1)==0;
                0.9<=v<=1.1;
                -pi<=theta<=pi;
                (0.9/1.1)-0.1<=vm<=(1.1/0.9)+0.1;
                %------------
                Pi-Pd==0;
                Qi-Qd==0;
                %------------
                (Pij.^2+Qij.^2)<=Smax.^2;
                %------------
                -1<=Qc<=SD.Qc_max.*(1.1)^2+1;
                ];

Objective_x = sum(Ploss)+(rho/2)*norm([vm-(v(SD.oltc_l(:,1))./(0.89375+u*0.00625));(Qc-((u_c*SD.Cstp).*v(SD.cb_l).^2))*beta]+lambda,2)^2;

Sol_x = optimize(Constraints_x,Objective_x,sdpsettings('solver','ipopt','verbose',0));

flag_x=Sol_x.problem;
if flag_x == 0
    obj_x=value(Objective_x);
    P_loss=value(sum(Ploss));
    v=value(v);
    vm=value(vm);
    theta=value(theta);
    Qc=value(Qc);
    Time_x=Sol_x.solvertime;
elseif flag_x == 1 || 4
    disp('let us relax Vm and Qs constraints even more');
    Constraints_x_relaxed =[v(1)==1;
                theta(1)==0;
                0.9<=v<=1.1;
                -pi<=theta<=pi;
                (0.9/1.1)-2<=vm<=(1.1/0.9)+2;
                %------------
                Pi-Pd==0;
                Qi-Qd==0;
                %------------
                (Pij.^2+Qij.^2)<=Smax.^2;
                %------------
                -5<=Qc<=SD.Qc_max.*(1.1)^2+5];
   Sol_x = optimize(Constraints_x_relaxed,Objective_x,sdpsettings('solver','ipopt','verbose',0));

   flag_x=Sol_x.problem;
   if flag_x == 0
    obj_x=value(Objective_x);
    P_loss=value(sum(Ploss));
    v=value(v);
    vm=value(vm);
    theta=value(theta);
    Qc=value(Qc);
   end
else
    Sol_x.info
end
end