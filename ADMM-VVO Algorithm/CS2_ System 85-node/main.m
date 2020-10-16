% Programmed by A. Alburidy and L. Fan
% arburidy@gmail.com
% If you find this code useful for your research, please cite our paper at:
% https://github.com/alburidy/ADMM-VVO-Optimization
%==========================
yalmip('clear');
clear; clc; close all;
%--------------------------
mpc = loadcase('case85');   % Enter any matpower radial case
%--------------------------
SD.oltc_l=[1 2;26 27;57 58];% Enter LTCs locations
SD.cb_l=[12;33;68];         % Enter SCBs locations
%--------------------------
Qc_max = [400;800;600]; % Enter Qc max. capacity in kVARs
SD.Cstp = 0.02;         % Enter Qc incremental step change value.
SD.Cinr = [2;4;3];      % Enter max. integer number for bsh
%--------------------------
SD.tpn=size(SD.oltc_l,1);
SD.cbn=size(SD.cb_l,1);
SD.Qc_max=Qc_max/(1e3*mpc.baseMVA);
%--------------------------
% Setting rho and lambda initial values:
lambda = [repmat(1e-1,SD.tpn,1);repmat(1e-2,SD.cbn,1)];
rho = 15;
beta=5e-1; % SCB Scaling Factor
%--------------------------
% Setting rho updating criteria:
M=10;
T=2;
Kf=100;
%###################################################
% You do NOT need to change anything after this line
%---------------------------------------------------
% Model parameters
epison_1=1e-6;
epison_2=1e-6;
itr_max=250;
test_1=1;
test_2=1;
itr=0;
u=repmat(17,SD.tpn,1);
u_c=zeros(SD.cbn,1);
%--------------------------
nn=size(mpc.bus,1);
nbr=size(mpc.branch,1);
Pd=(mpc.bus(2:end,3)/mpc.baseMVA);  % active Load power demands
Qd=(mpc.bus(2:end,4)/mpc.baseMVA);  % reactive Load power demands
Smax=mpc.branch(:,6)./mpc.baseMVA;  % load lines maximum loading level
Smax(Smax==0) = 2.5;
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
%--------------------------------------
% Ploting the convergence progress while iterating:
    figure(1)
    colors = hsv(SD.tpn+SD.cbn);
    % Naming legends of the convergence plot
    for k=1:SD.tpn
        legends{k,1}=sprintf('LTC# %d',k);
    end
    i1=0;
    for k=SD.tpn+1:SD.tpn+SD.cbn
        i1=i1+1;
        legends{k,1}=sprintf('CB on Node %d',SD.cb_l(i1));
    end
%--------------------------------------
disp('iteration# | obj. Value | max(V) | min(V) | rho Value | primal residual | dual residual | u ');
disp('--------------------------------------------------------------------------------------------');
total_time=tic;
while (test_1>epison_1 || test_2>epison_2) && itr<=itr_max
    % The primal update step:
    [obj_x,P_loss,vm,Qc,v,theta,flag_x,Time_x]=vvc(SD,itr,nn,rho,u,u_c,lambda,Gi,Bi,Pd,Qd,Smax,beta);
    % The slack update step:
    [Obj_z,u_n,u_c_n,flag_z,Time_z]=oltc_optimal_value(SD,itr,v,vm,rho,lambda,Qc,beta);
    
    u_n=round(value(u_n),0);
    u_c_n=round(value(u_c_n),0);
    
    % if either one of the problems did NOT converge STOP!
    if (flag_x ~= 0 || flag_z ~= 0)
        fprintf('flag_x:%1.0f flag_z:%1.0f\n',flag_x,flag_z);
        break;
    end
        
    % Lagrangia multiplier update:
    lambda_u = lambda+value([vm-(v(SD.oltc_l(:,1))./(0.89375+u_n*0.00625));(Qc-((u_c_n*SD.Cstp).*v(SD.cb_l).^2))*beta]);
    
    % Termination checks:
    % 1- Check the primal residual:
    test_1 = value(norm(lambda_u-lambda));
    % 2- Check the dual residual:
    test_2 = value(rho*norm([u_n;u_c_n]-[u;u_c]));
    
    % iteration counter
    itr=itr+1;
    
    % Print current iteration results
    fprintf(['%6.0f %12.5f %12.4f %8.4f %10.3f %15f %13.0f      |',...
        repmat(' %2.0f ',1, SD.tpn+SD.cbn),'\n']',itr,P_loss*mpc.baseMVA*1e3,max(v),min(v),rho,test_1,test_2,[u;u_c]');
    
    % Record history for graph and analysis
    history.obj(itr,1)=obj_x;
    history.test_1(itr,1)=test_1;
    history.test_2(itr,1)=test_2;
    history.lambda(itr,:)=lambda';
    history.convergence(itr,:)=value([vm-(v(SD.oltc_l(:,1))./(0.89375+u*0.00625));(Qc-((u_c*SD.Cstp).*v(SD.cb_l).^2))*beta]);
    history.rho(itr,1)=rho;
    
    Executing_time(itr,1) = Time_x+Time_z;
    
    % Assign Updated values to variales:
    lambda = value(lambda_u);
    u=u_n;
    u_c=u_c_n;
    
    % Updating penalty parameter:
    if itr >=4 && mod(Kf,itr)==0
    if test_1>M*test_2
        rho=rho*T;
    elseif test_2>M*test_1
        rho=rho/T;
    end
    end

    % Ploting the convergence progress
    PLOT = plot(1:itr,history.convergence(1:itr,:),'LineWidth', 1.5);
    set(PLOT, {'color'}, num2cell(colors, 2));
    grid on;
    legend(legends);
    set(legend,'position',[0.6,0.3,0.1,0.1]);
    xlim([1 itr+1]);
    title('Convergence Progress'); xlabel('iteration #');
    hold on
end
timer=toc(total_time);
xlim([1 itr]);
hold off
disp('--------------------------------------------------------------------------------------------');
if  (itr<=itr_max && flag_x==0 && flag_z==0)
    
    fprintf('The whole model execution time: %d minutes and %f seconds\n\n', floor(timer/60), rem(timer,60));
    fprintf('Total Solver execution time: %d minutes and %f seconds\n', floor(sum(Executing_time)/60), rem(sum(Executing_time),60));
    
    mpc.branch(SD.br_oltc,9)=(0.89375+u.*0.00625);
    mpc.bus(SD.cb_l,6)=(u_c.*SD.Cstp)*mpc.baseMVA;
    
    matpwr=runpf(mpc,mpoption('verbose',0,'out.all',0));
    sum(real(get_losses(matpwr)))-P_loss*mpc.baseMVA
    max(abs(matpwr.bus(:,8)-v))
    max(abs(deg2rad(matpwr.bus(:,9))-theta))
    plot_convergance;
elseif itr>itr_max
    disp('Maximum iteration steps is reached!');
end
