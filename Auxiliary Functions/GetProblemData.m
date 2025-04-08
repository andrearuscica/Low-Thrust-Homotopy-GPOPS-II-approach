function [problem_data] = GetProblemData(trajectory)
GTOC4Asteroids = readtable('GTOC4.txt');

% Convert Asteroid ID's to a string for manipulation
Asteroid_Names = string(GTOC4Asteroids{:,1});
Asteroid_Names = erase(Asteroid_Names,"'");  % Remove ' for string comparison
pars_struct = load("1.GTOC_Data/pars.mat");
pars = pars_struct.pars;

if trajectory ==  3  %'Moscow'
    Flag = 3;
    arc_num = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';
    '18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';
    '35';'36';'37';'38';'39';'40';'41';'42'; '43';'44';'45';'46';'47';'48'};
    asteroids_num = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';
    '18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';
    '35';'36';'37';'38';'39';'40';'41';'42'; '43';'44';'45';'46';'47';'48';'49'};

elseif trajectory ==  1  %'Johnson'
    Flag = 1;
    arc_num = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';
    '18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';
    '35';'36';'37';'38';'39';'40';'41';'42'; '43';'44';'45'};
    asteroids_num = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';
    '18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';
    '35';'36';'37';'38';'39';'40';'41';'42'; '43';'44';'45';'46'};
elseif trajectory ==  4  %'Moscow Original'
    Flag = 4;
    arc_num = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';
    '18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';
    '35';'36';'37';'38';'39';'40';'41';'42'; '43';'44'};
    asteroids_num = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';
    '18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';
    '35';'36';'37';'38';'39';'40';'41';'42'; '43';'44';'45'};
elseif trajectory == 2   %'Barbee'
    Flag = 2;
    arc_num = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';
    '18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';
    '35';'36';'37';'38';'39'};
    asteroids_num = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';
    '18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';
    '35';'36';'37';'38';'39';'40'};
end

[solution] = Data_GTOC4_Solution(Flag, Asteroid_Names);

% Loading txt file with Asteroid Data
GTOC4Asteroids = readtable("1.GTOC_Data/GTOC4.txt");
Asteroids_OE = GTOC4Asteroids(solution.Asteroid_Seq_ID,3:8);

%Asteroid Flyby Dates
Ref_Date = 54800; %[MJD]
Flyby_Dates = zeros(length(solution.Asteroid_Seq_ID),1);
for i = [1:length(solution.Asteroid_Seq_ID)]
    Flyby_Dates(i) = solution.Flyby_Dates(1,i);
end

ToFs = zeros(length(solution.Asteroid_Seq_ID)-1,1);
for i = [1:length(solution.Asteroid_Seq_ID)-1]
    ToFs(i) = (Flyby_Dates(i+1) - Flyby_Dates(i))*86400;
end

Asteroids_OEs = zeros(length(solution.Asteroid_Seq_ID),6);
for i = [1:length(solution.Asteroid_Seq_ID)]
    Asteroids_OEs(i,:) = table2array(Asteroids_OE(i,:));
end

%asteroids cartesian states
CartStates = zeros(length(solution.Asteroid_Seq_ID),6);
for i = [1:length(solution.Asteroid_Seq_ID)]
    [r,v] = CartesianState(Asteroids_OEs(i,:), (Flyby_Dates(i)-Ref_Date)*86400, pars);
    CartStates(i,:) = [r,v];
end

EarthOE = [0.999988049532578, 1.671681163e-2, 0.8854353e-3, 175.4064769, 287.61577546, 257.60683707];
[rE,vE] = CartesianState(EarthOE,(solution.Launch_Date-Ref_Date+800)*86400, pars);
dvLaunch = [-0.205100, -1.436312, 0.072829];
vLaunch = [vE + dvLaunch];
ToF_launch = (Flyby_Dates(1)-solution.Launch_Date)*86400;
[a,p,e,ERROR,VI_E,VF_E,TPAR,THETA] = lambertMR(rE,CartStates(1,1:3),((Flyby_Dates(1)-solution.Launch_Date)*86400),pars.mu_sun,0,0,0,1);


arc_ri = zeros(length(solution.Asteroid_Seq_ID)-1,3);
arc_rf = zeros(length(solution.Asteroid_Seq_ID)-1,3);
arc_VI = zeros(length(solution.Asteroid_Seq_ID)-1,3);
arc_VF = zeros(length(solution.Asteroid_Seq_ID)-1,3);
masses = zeros(length(solution.Asteroid_Seq_ID)-1,1);
MEE_end = zeros(length(solution.Asteroid_Seq_ID)-1,6);
for i = [1:length(solution.Asteroid_Seq_ID)-1]
[a,p,e,ERROR,VI,VF,TPAR,THETA] = lambertMR(CartStates(i,1:3), CartStates(i+1,1:3), ToFs(i), pars.mu_sun, 0, 0, 0, 1);
arc_ri(i,1:3) = CartStates(i,1:3);
arc_rf(i,1:3) = CartStates(i+1,1:3);
arc_VI(i,1:3) = VI;
arc_VF(i,1:3) = VF;
masses(i,1) = solution.SC_Mass(i);
[aOE0,eOE0,i0,RAAN0,w0,MA0,nu0] = orbitalElements(arc_rf(i,:),arc_VF(i,:),pars);
arc_OE = [aOE0,eOE0,i0,RAAN0,w0,nu0];
arc_MEE_end = orbital2equinoctial(arc_OE);
MEE_end(i,:) = arc_MEE_end;
end


problem_data.arcs = table(arc_num,arc_ri,arc_rf,arc_VI,arc_VF,ToFs,masses);
problem_data.MEE = table(arc_num,MEE_end);
problem_data.Asteroids_CSV = CartStates;
problem_data.Asteroids_OE = table(asteroids_num,Asteroids_OEs,Flyby_Dates);
problem_data.EarthArc = VF_E; % To be changed taking into account Earth injection with vinf and the little thrusting
problem_data.Launch_state_ToF = [rE, vLaunch, ToF_launch];

end