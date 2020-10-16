% Programmed by A. Alburidy and L. Fan
% arburidy@gmail.com
% If you find this code useful for your research, please cite our paper at:
% https://github.com/alburidy/ADMM-VVO-Optimization
%==========================
function [Obj_z,u,u_c,flag_z,Time_z]=oltc_optimal_value(SD,itr,v,vm,rho,lambda,Qc,beta)

u= sdpvar(SD.tpn,1);
u_c=sdpvar(SD.cbn,1);

if itr<1
    assign(u,17);
    assign(u_c,0);
end

Objective_z = (rho/2)*norm([vm-(v(SD.oltc_l(:,1))./(0.89375+u*0.00625));(Qc-((u_c*SD.Cstp).*v(SD.cb_l).^2))*beta]+lambda,2)^2;

Constraints_z = [1<=u<=33;
                 0<=u_c<=SD.Cinr];

Sol_z = optimize(Constraints_z,Objective_z,sdpsettings('solver','ipopt','verbose',0));

flag_z=Sol_z.problem;

if flag_z == 0
    Obj_z=value(Objective_z);
    u=value(u);
    u_c=value(u_c);
    Time_z=Sol_z.solvertime;
else
    Sol_z.info
end

end