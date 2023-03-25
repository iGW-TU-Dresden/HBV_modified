clear all
close all

load input_parameters.mat %Input time series including warm-up
load Date.mat

%Input separation
    P=InFlow{1};
    ETpot=InFlow{2};   
    T=InFlow{2};   

%Catchment area
A=155380; %m³

%Soil parameters
parBETA = Param(1);
parFC = Param(2);

parLP = Param(3);

%Hillslope & Groundwater parameters (groundwater directly passed to the
%stream)
parK0 = 10^Param(4);
parK1 = 10^Param(5);

parUZL = Param(6);

%Riparian zone parameters (feed by hillslope discharge)
parK0_RZ = 1; %Excess water directly passed to the stream
parK1_RZ = 10^Param(7);

parUZL_RZ = Param(8);

%Rain % Snow parameters
parPCORR = 1;%Param(11); %can ssumed to be 1
parTT = 0;%Param(12); %can assumed to be 0 (doesn't matter)
parCFMAX = 5;%Param(13); %can assumed to be 5 (doesn't matter)
parSFCF = 1;%Param(14); %can assumed to be 1 (doesn't matter)
parCFR = .01;%Param(15); %can assumed to be 0,01 (doesn't matter)
parCWH = .1;%Param(16); %can assumed to be 0.1 (doesn't matter)

% State variables
SNOWPACK = nan(size(P)); SNOWPACK(1,:) = 0.0001;
MELTWATER = nan(size(P)); MELTWATER(1,:) = 0.0001;
SM = nan(size(P)); SM(1,:) = 0.0001;
SUZ = nan(size(P)); SUZ(1,:) = 0.0001;
SLZ = nan(size(P)); SLZ(1,:) = 0.0001;
ETact = nan(size(P)); ETact(1,:) = 0.0001;
Qsim =zeros(size(P)); Qsim(1,:) = 0.0001;
Q0 = zeros(size(P)); 
Q1 = zeros(size(P)); 

% State variables Riparian Zone
SUZ_RZ = nan(size(P)); SUZ_RZ(1,:) = 0.0001;
SLZ_RZ = nan(size(P)); SLZ_RZ(1,:) = 0.0001;
Q0_RZ = zeros(size(P));
Q1_RZ = zeros(size(P)); 

for ii = length(P(1,:)):-1:1
    P(:,ii) = parPCORR(ii).*P(:,ii);
    SNOW(:,ii) = P(:,ii);
    RAIN(:,ii) = P(:,ii);
    RAIN(T(:,ii)<parTT(ii),ii) = 0;
    SNOW(T(:,ii)>=parTT(ii),ii) = 0;
    SNOW(:,ii) = SNOW(:,ii).*parSFCF(ii);
end

for t = 2:length(P)
    
    % Snow
    SNOWPACK(t,:) = SNOWPACK(t-1,:)+SNOW(t,:);
    melt = parCFMAX .* (T(t,:)-parTT);
    melt(melt<0) = 0;
    melt = min(melt,SNOWPACK(t,:));
    MELTWATER(t,:) = MELTWATER(t-1,:)+melt;
    SNOWPACK(t,:) = SNOWPACK(t,:)-melt;
    refreezing = parCFR .* parCFMAX .* (parTT-T(t,:));
    refreezing(refreezing<0) = 0;
    refreezing = min(refreezing,MELTWATER(t,:));
    SNOWPACK(t,:) = SNOWPACK(t,:)+refreezing;
    MELTWATER(t,:) = MELTWATER(t,:)-refreezing;
    tosoil = MELTWATER(t,:) - (parCWH .* SNOWPACK(t,:));
    tosoil(tosoil<0) = 0;
    MELTWATER(t,:) = MELTWATER(t,:)-tosoil;
    
    % Soil and evaporation
    soil_wetness = (SM(t-1,:)./parFC).^parBETA;
    soil_wetness(soil_wetness<0) = 0;
    soil_wetness(soil_wetness>1) = 1;
    recharge = (RAIN(t,:)+tosoil) .* soil_wetness;
    SM(t,:) = SM(t-1,:)+RAIN(t,:)+tosoil-recharge;
    excess = SM(t,:)-parFC;
    excess(excess<0) = 0;
    SM(t,:) = SM(t,:)-excess;
    evapfactor = SM(t,:) ./ (parLP .* parFC);
    evapfactor(evapfactor<0) = 0;
    evapfactor(evapfactor>1) = 1;
    ETact(t,:) = ETpot(t,:).*evapfactor;
    ETact(t,:) = min(SM(t,:), ETact(t,:));
    SM(t,:) = SM(t,:)-ETact(t,:);
    
    % Hillslope & Groundwater boxes
    SUZ(t,:) = SUZ(t-1,:)+recharge+excess;
    Q0(t,:) = parK0 .* max(SUZ(t,:)-parUZL, 0.0);
    SUZ(t,:) = SUZ(t,:)-Q0(t,:);
    Q1(t,:) = parK1.*SUZ(t,:);
    SUZ(t,:) = SUZ(t,:)-Q1(t,:);
    
    % Riparian Zone boxes
    SUZ_RZ(t,:) = SUZ_RZ(t-1,:)+Q0(t,:); %Inflow from hillslope discharge
    Q0_RZ(t,:) = parK0_RZ .* max(SUZ_RZ(t,:)-parUZL_RZ, 0.0);
    SUZ_RZ(t,:) = SUZ_RZ(t,:)-Q0_RZ(t,:);
    Q1_RZ(t,:) = parK1_RZ.*SUZ_RZ(t,:);
    SUZ_RZ(t,:) = SUZ_RZ(t,:)-Q1_RZ(t,:);
    
    %Total discharge
    Qsim(t,:) = Q1(t,:) + Q0_RZ(t,:) + Q1_RZ(t,:); %Groundwater + Hillslope(bypassing the RZ) + Riparian Zone
end

simQ=Qsim/1000/86400*A;%mm/d->m³/s  Simulated catchment discharge

plot(Datesim,simQ)
datetick('x',10)
