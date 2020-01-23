%%
% Author: Tobias Gerken, tug15@psu.edu
% Version: 0.1 (2020_01_22)

% Script for running on ACI Cluster to compute variogram and bootstrap
% samples

function funOut = CalculateBootstrap(cs, nboot)

%% Personal header to adjust paths 
if strcmp(getenv('computername'),'DESKTOP-45CVB98')
    addpath(genpath('C:\Users\tobia\OneDrive - The Pennsylvania State University\Projects\ACT-America\Code\GitCode'));
    basepath = 'C:\Users\tobia\';
elseif strcmp(getenv('computername'),'DESKTOP-114H9OU')
    addpath(genpath('C:\Users\tobia\OneDrive - The Pennsylvania State University\Projects\ACT-America\Code\GitCode'));
    basepath = 'C:\Users\tobia\';
elseif strcmp(getenv('computername'),'DESKTOP-A2GKIRA')
    addpath(genpath('D:\OneDrive - The Pennsylvania State University\Projects\ACT-America\Code\GitCode'));
    basepath = 'D:\';
elseif strcmp(getenv('computername'),  'E2-MET-WKDT013')
    addpath(genpath('C:\Users\tug15\OneDrive - The Pennsylvania State University\Projects\ACT-America\Code\GitCode'));
     basepath = 'C:\Users\tug15\';  
else
    basepath = '.\';
end
%set(0,'DefaultLegendAutoUpdate','off') % from 2017 onwards legends autoupdate, add to startup.m
%set(0, 'DefaultFigureVisible', 'on')
%% Global Settings 

% use some kind of geometric scaling for variogram bins 
n =40;
r = 1.2;
a = 1;
s = a*r.^(0:n-1);
s = [0 unique(floor(s))];
s = s(s<500);
edges = s ;
maxd = edges(end);

clear s a n r

%if ~exist('Cases')
%    error('need to pass Cases argument')
%end
%disp(Cases)
%disp(nboot)

if cs ==1
Cases = {'Regions'};
elseif cs ==2
Cases = {'Airmasses'};
elseif cs ==3
    Cases = {'Seasons'};
end

Levels = {'ABL','LFT','HFT'} ;   
Regions = {'NEMA','SC','MWe'} ;
Seasons = {'WIN','SPR','SUM', 'FAL'} ;
Airmasses = {'Fair', 'Cold', 'Warm'};

Nme.Regions = Regions; 
Nme.Seasons = Seasons;
Nme.Airmasses = Airmasses;

% Bootstrapping parameters
%nboot = 1000; % number of bootstraps (allow some additionals for non-convergence)
BlockLengthMeth = 'Optimal' ;% calculates blocklength based on autocorr, default is simple which is n^(1/3)

% Set models

models = {'exponential', ...
          'gaussian'};
modeltype = {'unbounded',...
             'unbounded'};

modFun = {@(b,h)b(3)+b(2)*(1-exp(-h./b(1))), ... % exponential model
          @(b,h)b(3)+b(2)*(1-exp(-(h.^2)/(b(1)^2))), ...% gaussian model
          };
% select which models to do 
DoModels = {'exponential'};

% 
modifier = [1 1 1 ; 0.5 1 1; 2 1 1; 1 0.5 1; 1 2 1 ;0.5 0.5 1; 2 2 1; 0.5 2 1; 2 0.5 1] ;
            
%% Load data and 

% Load Data 
DataDir = ['./InDir/'] ;
OutDir = ['./OutDir/'] ;
File = 'AllLegs.mat' ;

load([DataDir File]);

CO2.Obs = [AllLegs.CO2_Obs];
CO2.WRF = [AllLegs.CO2_WRF]; % for now because I am dumb, remove after rerun of AllLegs 
CO2.CT = [AllLegs.CO2_CT]; % for now because I am dumb, remove after rerun of AllLegs 

CO2.Res_WRF = [CO2.Obs-CO2.WRF];
CO2.Res_CT = [CO2.Obs-CO2.CT];

%% Set flags for data selection
t = AllLegs.time;
StartDay = floor(min(t));
EndDay = floor(max(t));

Flag.Level =  nan(size(CO2.Obs)) ;
Flag.Level(AllLegs.z_AGL >= 4000) = 3 ;
Flag.Level(AllLegs.z_AGL < 4000 & AllLegs.z_AGL >=1500) = 2 ;
Flag.Level(AllLegs.z_AGL <1500) = 1 ;

Flag.Region =  nan(size(CO2.Obs)) ;
%Flag.Region(AllLegs.Lat >=35.77 &  AllLegs.Lon >=-87.5  )=  1 ;
Flag.Region(AllLegs.Lat >=33.75 &  AllLegs.Lon >=-87.5  )=  1 ;
Flag.Region(AllLegs.Lat < 37.00 &  AllLegs.Lon < -84.39 )=  2 ;
Flag.Region(AllLegs.Lat >=37.00 &  AllLegs.Lon < -87.5  )=  3 ;

Flag.Airmass = nan(size(CO2.Obs)) ;
Flag.Airmass(AllLegs.AMF == 0) = 1 ;
Flag.Airmass(AllLegs.AMF == 1) = 2 ;
Flag.Airmass(AllLegs.AMF == 2) = 3 ;

% Seasons: Winter : Fall = 1:4
m= month(AllLegs.time);
Flag.Season =  nan(size(CO2.Obs)) ;
Flag.Season(ismember(m,[12 1 2]))=1;
Flag.Season(ismember(m,[3:5]))=2;
Flag.Season(ismember(m,[6:8]))=3;
Flag.Season(ismember(m,[9:11]))=4;
% Regions
Flag.Region =  nan(size(CO2.Obs)) ;
Flag.Region(AllLegs.Lat >=33.75 &  AllLegs.Lon >=-87.5  )=  1 ;
Flag.Region(AllLegs.Lat < 37.00 &  AllLegs.Lon < -84.39 )=  2 ;
Flag.Region(AllLegs.Lat >=37.00 &  AllLegs.Lon < -87.5  )=  3 ;
% Airmass Flags 
Flag.Airmass = nan(size(CO2.Obs)) ;
Flag.Airmass(AllLegs.AMF == 0) = 1 ;
Flag.Airmass(AllLegs.AMF == 1) = 2 ;
Flag.Airmass(AllLegs.AMF == 2) = 3 ;


%% Assemble and save distance pairs over all days 


% loop over days
Lat =  AllLegs.Lat;
Lon =  AllLegs.Lon;
z =  AllLegs.z_AGL;

for mod = {'CT', 'WRF'} 

    if strcmp(mod{1}, 'WRF')
        c  = CO2.Res_WRF;
    elseif strcmp(mod{1}, 'CT')
        c  = CO2.Res_WRF;
    end

% Select data by season and level 
    for ii = 1:length(Seasons)
        for jj = 1:length(Levels)
        
            iCase = (Flag.Season == ii);
            iLev  = (Flag.Level == jj);
        
            dummy(1,:) = Lat(iLev & iCase);
            dummy(2,:) = Lon(iLev & iCase);
            dummy(3,:) = t(iLev & iCase);
            dummy(4,:) = c(iLev & iCase);
        
            dummy=dummy(:,isfinite(dummy(4,:)));
        
            D.(mod{1}).(Seasons{ii}).(Levels{jj})=dummy;
            clear dummy
        end
    end
    % Do the Same for Region
    for ii = 1:length(Regions)
        for jj = 1:length(Levels)
        
            iCase = (Flag.Region == ii);
            iLev  = (Flag.Level == jj);
        
            dummy(1,:) = Lat(iLev & iCase);
            dummy(2,:) = Lon(iLev & iCase);
            dummy(3,:) = t(iLev & iCase);
            dummy(4,:) = c(iLev & iCase);
        
            dummy=dummy(:,isfinite(dummy(4,:)));
        
            D.(mod{1}).(Regions{ii}).(Levels{jj})=dummy;
            clear dummy
        end
    end
    % Do the Same for Airmass
    for ii = 1:length(Airmasses)
        for jj = 1:length(Levels)
        
            iCase = (Flag.Airmass == ii);
            iLev  = (Flag.Level == jj);
        
            dummy(1,:) = Lat(iLev & iCase);
            dummy(2,:) = Lon(iLev & iCase);
            dummy(3,:) = t(iLev & iCase);
            dummy(4,:) = c(iLev & iCase);
        
            dummy=dummy(:,isfinite(dummy(4,:)));
        
            D.(mod{1}).(Airmasses{ii}).(Levels{jj})=dummy;
            clear dummy
        end
    end
end

% get rid of some variables
clear Flag CO2 Lat Lon

%% Now start the actual calculation loop
for c = Cases
    for nm = 1:length(Nme.(c{1}))     
        for mod = {'CT', 'WRF'}
            for Level = Levels
                disp([Nme.(c{1}){nm} ' ' mod{1} '_' Level{1} ' Nboot: ' num2str(nboot)]  )
                tic
                % select the appropriate data 
                 Data = D.(mod{1}).(Nme.(c{1}){nm}).(Level{1}) ; 
            
                    n = length(Data);
                if strcmp(BlockLengthMeth, 'Optimal') 
                    b = ceil(opt_block_length(Data(4,:)'));
                    b = b(2);
                else
                    b = round(length(Data)^(1/3));
                end
        
                nb = ceil(n/b);
        
        
                % 1) Do the best guess for level and season
                %       - Calculate empirical varigram over all data 
                %       - calculate a best fit using a specified model 
                %       - Plot results 
        
                % % do preparations 
                DataCirc = [Data  Data(:,1:b(1))];
                N_tot = nan(nboot+1, length(edges)-1);
                Gamma = nan(nboot+1, length(edges)-1);
                G_sigma = nan(nboot+1, length(edges)-1);
                for m = DoModels
                    Gamma_pred.(m{1}) =  nan(nboot+1, length(edges)-1);
                end
                % draw the block bootstrap samples 
                boots  = randi(n,[nb nboot]);
        
                % Now loop over samples with the first being the original data 
        
                for sample =  1:nboot+1 
                    disp(sample)
                    if sample == 1
                        DataSample = Data; % The first sample is the original data to get a best guess
                    else % otherwise use the bootstrap samples 
                        iii = nan(b,nb);
                        for i = 1:nb
                            iii(:,i) = boots(i,sample-1):boots(i,sample-1)+b-1;
                        end
                        iii=iii(1:n);
                        DataSample= DataCirc(:,iii);
                    end
            
                    % Allocate arrays for speed
                    N = nan(60, length(edges)-1);
                    V = nan(60, length(edges)-1);
                    Var = nan(60, length(edges)-1);
        
                    % Calculate distances for each day 
                    ct =1;
                    for day = StartDay:EndDay  
                        datestr(day);
                        dTh = 0.2 ; % Allow for some overflow into next day, since flights were conducted from morning to evening UTC
                        iDy = (DataSample(3,:)>=day+ dTh & DataSample(3,:) <day + 1 + dTh) ;
    
                        if sum(iDy)==0
                            continue
                        end
            
                        % calculate distances 
                        lt = DataSample(1,iDy); 
                        lo = DataSample(2,iDy); 
                        y  = DataSample(4,iDy);
                        x(:,2) = distance(lt,lo(1), lt(1), lo(1), referenceEllipsoid('wgs84'))/1000;
                        x(:,1) = distance(lt(1),lo, lt(1), lo(1), referenceEllipsoid('wgs84'))/1000;
    
                        d = variogram_mod(x,y', 'maxdist', maxd, 'plot',false, 'edges',edges);
                        Distance = d.distance ;
                        N(ct,:) = d.num;
                        V(ct,:) = d.val;
                        Var(ct,:) = d.variance;
    
                        ct = ct+1;
                        clear x % clear to avoid issues 
                    end % day 

                    % Now assemble overall statistics 
                    N_tot(sample,:) = nansum(N(:,:));
                    Gamma(sample,:) = nansum(V(:,:).*N(:,:))./N_tot(sample,:);

                    G_sigma(sample,:) = sqrt(...
                                        1/(sum(N_tot(sample,:))-1) .* ...
                                        ( ... 
                                        nansum( (N(:,:)-1).*V(:,:) + ...
                                        N(:,:).*V(:,:).^2 ) - ...
                                        nansum(N(:,:)).*Gamma(sample,:).^2 ...
                                        ));
                    end % Sample 
                
                    % extract best guess variogram
                    %Gamma_org.(Season{1}).(Level{1}) =Gamma(1,:);
                    %N_tot_org.(Season{1}).(Level{1}) =N_tot(1,:);
                    %G_sigma_org.(Season{1}).(Level{1}) =G_sigma(1,:);
          
                    % Now fit the best model and calculate the bootstrapped CI
                    for m = DoModels
                        Index = find(contains(models,m{1}));
            
                        for sample = 1:size(Gamma,1)
                            i_maxd=find(isfinite(Gamma(sample,:)),1,'last');
                            start = [Distance(i_maxd)*2/3 , max(Gamma(sample,:)) 0];
                            
                            if strcmp(modeltype{Index}, 'unbounded')
                                start(1) = start(1)/3;
                            end
                            repeat = true;
                            pp=1;
                            while repeat
                                try
                                    model= fitnlm(Distance,Gamma(sample,:),modFun{Index},start.*modifier(pp,:),'Weight',1./G_sigma(sample,:));
                                    repeat = false;
                                catch
                                    pp = pp+1;
                                end
                            end
                            
                            % save parameters  
                            Coefs.(m{1})(sample,:) =  model.Coefficients.Estimate;
                            Gamma_pred.(m{1})(sample,:) =  predict(model,Distance,'Simultaneous',true);
                            if sample == 1
                                model_data.(m{1}) = model;
                            end
                
                        end
                        CI.(m{1}) = quantile(Gamma_pred.(m{1})(2:end,:),[0.05, 0.95]);
                    end
                    toc
                    timeElapsed = tic;
                    save([OutDir 'VarOutput_' mod{1} '_' Nme.(c{1}){nm} '_' Level{1} '.mat'], 'CI', 'Gamma', 'Gamma_pred', 'Coefs', 'model_data', 'G_sigma', 'N_tot', 'b','boots', 'timeElapsed')
                    
            end
        end
     end
end

funOut = 1;
end