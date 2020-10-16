% Programmed by A. Alburidy and L. Fan
% arburidy@gmail.com
% If you find this code useful for your research, please cite our paper at:
% https://github.com/alburidy/ADMM-VVO-Optimization
%==========================
% ask the user to choose a case study
picked_sys = input(['Which case study you want to test?',...
       '\n Type 1 for the 85-Node test.',...
       '\n Type 2 for the 141-Node test',...
       '\n Type 3 for the heavily equipped 141-Node test. \n']);


if picked_sys == 1
    mpc = loadcase('case85');   % load system data using MATPOWER's function 'loadcase'
    Nc=[12,33,68];                       % SCBs locations
    Qc=[400,800,600]./(1e3*mpc.baseMVA); % Qc max. capacity in kVARs
    Kc=[2,4,3];                          % max. integer number for bsh
    Cstp=0.02;                           % Qc incremental step change value.
    VR_bus = [1;26;57]; 	% Enter LTCs locations 'from bus'
    VR_bus_to = [2,27,58];  % Enter LTCs locations 'to bus'
    br_oltc = [1;26;57];    % Endexing LTCs locations
    
elseif picked_sys == 2
    
    mpc = loadcase('case141');   % load system data using MATPOWER's function 'loadcase'
    Nc=[8,48,90,112];            % SCBs locations
    Qc=[1200,1200,600,600]./(1e3*mpc.baseMVA);  % Qc max. capacity in kVARs
    Kc=[4,4,2,2];                               % max. integer number for bsh
    Cstp=0.03;                                  % Qc incremental step change value.
    VR_bus = [1;42;43];     % Enter LTCs locations 'from bus'
    VR_bus_to = [2,54,44];  % Enter LTCs locations 'to bus'
    br_oltc = [1;43;53];    % Endexing LTCs locations

else
    mpc = loadcase('case141');   % load system data using MATPOWER's function 'loadcase'
    Nc=[8;30;40;48;60;65;70;80;85;90;100;112];  % SCBs locations
    Qc=repmat(1200,12,1)./(1e3*mpc.baseMVA);    % Qc max. capacity in kVARs
    Kc=repmat(4,12,1);                          % max. integer number for bsh
    Cstp=0.03;                                  % Qc incremental step change value.
    VR_bus = [1;43;42;23];          % Enter LTCs locations 'from bus'
    VR_bus_to = [2;44;54;24];       % Enter LTCs locations 'to bus'
    br_oltc = [1;43;53;23];         % Endexing LTCs locations
end
