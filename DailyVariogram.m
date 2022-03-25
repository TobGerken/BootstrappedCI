% Author: Tobias Gerken, tug15@psu.edu
% Version: 0.1 (2019_03_19)
% Rev: Now adapted for new WRF output files
%

clear all; close all
warning('off','all')
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
set(0,'DefaultLegendAutoUpdate','off') % from 2017 onwards legends autoupdate, add to startup.m
set(0, 'DefaultFigureVisible', 'on')

%% set varigram parameters 

n =40;
r = 1.2;
a = 1;
s = a*r.^(0:n-1);
s = [0 unique(floor(s))];
s = s(s<500);
s(end)=500;
edges = s ;
maxd = edges(end);

clear s a n r


%% load data 

DataDir = [basepath 'OneDrive - The Pennsylvania State University\Projects\ACT-America\Data\Processing\ExtractedData\'] ;
PlotDir = [basepath 'OneDrive - The Pennsylvania State University\Projects\ACT-America\Code\ACT_ModelDataPaper\Plots\'] ;
Dir = [basepath 'OneDrive - The Pennsylvania State University\Projects\ACT-America\Code\ACT_ModelDataPaper\Code\'] ;


File = 'AllLegs.mat' ;
load([DataDir File])
% do this for WRF only now
CO2.Obs = [AllLegs.CO2_Obs];
CO2.WRF = [AllLegs.CO2_WRF];
CO2.CT = [AllLegs.CO2_CT];

CO2.Res_WRF = [CO2.WRF-CO2.Obs];
CO2.Res_CT  = [CO2.CT-CO2.Obs];

% Only retain finite residuals
t = AllLegs.time;
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

Cases = {'All', 'Season', 'Region', 'Airmass', 'Level'} ;
Level = {'ABL','LFT','HFT'} ;   
Region = {'NEMA','SC','MWe'} ;
Season = {'WIN','SPR','SUM', 'FAL'} ;
Airmass = {'Fair', 'Cold', 'Warm'};

Nme.Region = Region; 
Nme.Level = Level;
Nme.Season = Season;
Nme.Airmass = Airmass;

Label.Region = {'N.East Mid-Atl.', 'South Central', 'Mid-West'};
Label.Season = {'DJF', 'MAM', 'JJA', 'SON'};
Label.Airmass = {'Fair', 'Cold Sec', 'Warm Sec'};


StartDay = floor(min(t));
EndDay = floor(max(t));

%%
Var = struct();
set(0, 'DefaultFigureVisible', 'off')
for model = {'CT','WRF'}
    
    % preallocate vectors 
    
    Var.distance = 0.5*(edges(1:end-1)+edges(2:end)); 
    
    for lg = {'ABL','LFT','HFT'}  
        Var.(model{1}).(lg{1}).gamma=nan([200, numel( Var.distance)])  ;
        Var.(model{1}).(lg{1}).num=nan([200, numel( Var.distance)])  ;
        Var.(model{1}).(lg{1}).gamma_var=nan([200, numel( Var.distance)]) ;
    end   
    
    
Lat =  AllLegs.Lat;
Lon =  AllLegs.Lon;
c = CO2.(['Res_' model{1}]);
z =  AllLegs.z_AGL;

% remove absolute outliers defined as outside [0.0001 0.9999] quantiles
c(c<quantile(c,.0005)) = NaN;
c(c>quantile(c,.9995)) = NaN;

ct =1;

for day = StartDay:EndDay
    datestr(day)
    dTh = 0.2 ; % Allow for some overflow into next day, since flights were conducted from morning to evening UTC
    iDy = (t>=day+ dTh & t<day + 1 + dTh) ;
    if ~sum(iDy)==0  % only do days with data 
       
    % Check whether this is true 
    lt = Lat(iDy);
    lo = Lon(iDy);
    zD = z(iDy);
    y  = c(iDy);
     
    [ dx, dy ] = latlon2xy( lt, lo, nanmean(lt), nanmean(lo), true );
    clear x
    x(:,1) = dx;
    x(:,2) = dy;
    
    iFnt = isfinite(dx) & isfinite(dy) & isfinite(y) & isfinite(zD) ;

    Var.Date(ct) = day;
    Var.Region(ct) = mode(Flag.Region(iDy));
    Var.Season(ct) = mode(Flag.Season(iDy));
   
    
    for lg = {'ABL','LFT','HFT'} 

        iSel = iFnt;
        if strcmp(lg{1},'ABL')
            iSel = (iFnt & Flag.Level(iDy) ==  1) ;
        elseif strcmp(lg{1},'LFT')
            iSel = (iFnt & Flag.Level(iDy) ==  2) ;
        elseif strcmp(lg{1},'HFT')
            iSel = (iFnt & Flag.Level(iDy) ==  3) ;
        end
        if sum(iSel)>0
            %d = variogram_mod(x(iSel,:),y(iSel)', 'maxdist', maxd, 'plot',false, 'edges',edges);
            %Var.(lg{1}).val(ct,:) = d.val;
            %Var.(lg{1}).num(ct,:) = d.num;
            
            d = variogram_mod(x(iSel,:),y(iSel)', 'maxdist', maxd, 'plot',false, 'edges',edges);

            Var.(model{1}).(lg{1}).num(ct,:) = d.num;
            Var.(model{1}).(lg{1}).gamma(ct,:) = d.val;
            Var.(model{1}).(lg{1}).gamma_var(ct,:) = d.variance;
        end
    end
    ct = ct+1;
    end
end
end


for model = {'CT', 'WRF'}
for lg = {'ABL','LFT','HFT'}  
    Var.(model{1}).(lg{1}).gamma =    Var.(model{1}).(lg{1}).gamma(1:ct-1,:)  ;
    Var.(model{1}).(lg{1}).gamma_var= Var.(model{1}).(lg{1}).gamma_var(1:ct-1,:) ;
    Var.(model{1}).(lg{1}).num      = Var.(model{1}).(lg{1}).num(1:ct-1,:) ;

end   
end
save('Var')

%% Now do overall variogram 

col.WRF = [0.75 0 0];
col.CT  = [0 0 0.75];
col2.CT = [0    0.4470    0.7410];
col2.WRF= [0.8500    0.3250    0.0980];

figure('units','inches','position',[0 0 4 6])
ha = tight_subplot(3,1,[0.02 0.02], [0.05 0.05], [0.1 0.1])

ct = 1    
for level = {'ABL','LFT','HFT'}
    axes(ha(ct))
    hold on; box on;
    text(0.01, 0.9, ['\bf' level{1}], 'units','normalized')
    text(0.8, 0.9, ['\bf CT' ], 'color', col.CT',  'units','normalized')
    text(0.8, 0.8, ['\bf WRF' ], 'color', col.WRF',  'units','normalized')
    for model = {'WRF', 'CT'}
        
    N   =  Var.(model{1}).(level{1}).num;
    G   =  Var.(model{1}).(level{1}).gamma;
    V   =  Var.(model{1}).(level{1}).gamma_var;
    
    N_tot = nansum(N(:,:));
    Gamma = nansum(G(:,:).*N(:,:))./N_tot;

    G_sigma   = sqrt(...
                     1./(N_tot-1) .* ...
                    ( ... 
                     nansum( (N(:,:)-1).*V(:,:) + ...
                     N(:,:).*G(:,:).^2 ) - ...
                     N_tot.*Gamma.^2 ...
                     ));
                 

    line(Var.distance, nanmedian(G),'color',col.(model{1}), 'linestyle','--')
    errorbar(Var.distance, Gamma, G_sigma, 'o','color',col.(model{1}),'markerfacecolor',col.(model{1}))
    end
    yl = get(gca,'ylim')
    set(gca,'ylim',[0 yl(2)], 'xlim', [0 500])
ct = ct+1
end
set(ha(1:2), 'xticklabel','')
axes(ha(2))
ylabel('\bf \gamma(x) (ppm^2)')
axes(ha(3))
xlabel('\bf Distance x (km)')
print('Var_All.png','-dpng','-r300')

%% Plot all individually

    N   =  Var.(model{1}).ABL.num;
    G   =  Var.(model{1}).ABL.gamma;
    V   =  Var.(model{1}).ABL.gamma_var;
    
        N_tot = nansum(N(:,:));
    Gamma = nansum(G(:,:).*N(:,:))./N_tot;

    G_sigma   = sqrt(...
                     1./(N_tot-1) .* ...
                    ( ... 
                     nansum( (N(:,:)-1).*V(:,:) + ...
                     N(:,:).*G(:,:).^2 ) - ...
                     N_tot.*Gamma.^2 ...
                     ));
                 
    figure
    subplot(2,1,1)
    title([model{1} ' - ' level{1}])
    hold on; box on
    plot(Var.distance, G,'color','k')
    plot(Var.distance, Gamma,'color','r', 'linewidth',2)
    plot(Var.distance, nanmedian(G),'color','b', 'linewidth',2)
    ylabel('\gamma (ppm^2)')
    xlabel('Distance (km)')
    set(gca,'ylim',[0 50])
    subplot(2,1,2)
    hold on; box on
    plot(Var.distance, G,'color','k')
    plot(Var.distance, Gamma,'color','r', 'linewidth',2)
    plot(Var.distance, nanmedian(G),'color','b', 'linewidth',2)
    ylabel('\gamma (ppm^2)')
    xlabel('Distance (km)')  
    set(gca','yscale','log')
    print('VarExample.png','-dpng','-r300')
    

%% Now do this for season and region

for cs = {'Region', 'Season'}
 ny=length(Nme.(cs{1}))   
figure('units','inches','position',[0 0 10 6])
ha = tight_subplot(3,ny,[0.02 0.04], [0.05 0.05], [0.05 0.05])
   
ct = 1    
for level = {'ABL','LFT','HFT'}
    for run = 1:length(Nme.(cs{1}))
        axes(ha(ct))
         hold on; box on;
         i_run = (Var.(cs{1}) == run)
        if strcmp(level,'ABL')
            title([Nme.(cs{1}){run} ' n = ' num2str(sum(i_run))])
        end
         
    text(0.01, 0.9, ['\bf' level{1}], 'units','normalized')
    if ct == 4
        text(0.8, 0.9, ['\bf CT' ], 'color', col.CT',  'units','normalized')
        text(0.8, 0.8, ['\bf WRF' ], 'color', col.WRF',  'units','normalized')
    end
    for model = {'WRF', 'CT'}
        
    N   =  Var.(model{1}).(level{1}).num(i_run,:);
    G   =  Var.(model{1}).(level{1}).gamma(i_run,:);
    V   =  Var.(model{1}).(level{1}).gamma_var(i_run,:);
    
    N_tot = nansum(N(:,:));
    Gamma = nansum(G(:,:).*N(:,:))./N_tot;

    G_sigma   = sqrt(...
                     1./(N_tot-1) .* ...
                    ( ... 
                     nansum( (N(:,:)-1).*V(:,:) + ...
                     N(:,:).*G(:,:).^2 ) - ...
                     N_tot.*Gamma.^2 ...
                     ));
                 

    line(Var.distance, nanmedian(G),'color',col.(model{1}), 'linestyle','--')
    errorbar(Var.distance, Gamma, G_sigma, 'o','color',col.(model{1}),'markerfacecolor',col.(model{1}))
    end
    yl = get(gca,'ylim')
    set(gca,'ylim',[0 yl(2)], 'xlim', [0 500])
ct = ct+1
end
end
set(ha(1:ct-ny-1), 'xticklabel','')
axes(ha(ny+1))
ylabel('\bf \gamma(x) (ppm^2)')
axes(ha(end-1))
%xlabel('\bf Distance x (km)')
print(['Var_' cs{1} '.png'],'-dpng','-r300')
end


%% Fit exponential variogram
% Now do this for season and region
modifier = [1 1 1 ; 0.5 1 1; 2 1 1; 1 0.5 1; 1 2 1 ;0.5 0.5 1; 2 2 1; 0.5 2 1; 2 0.5 1] ;
modFun = @(b,h)b(3)+b(2)*(1-exp(-h./b(1))) ;

for cs = {'Region', 'Season'}
 ny=length(Nme.(cs{1}))   
figure('units','inches','position',[0 0 10 6])
ha = tight_subplot(3,ny,[0.02 0.04], [0.05 0.05], [0.05 0.05])

ct = 1    
for level = {'ABL','LFT','HFT'}
    for run = 1:length(Nme.(cs{1}))
        axes(ha(ct))
         hold on; box on;
         i_run = (Var.(cs{1}) == run);
        if strcmp(level,'ABL')
            title([Nme.(cs{1}){run} ' n = ' num2str(sum(i_run))])
        end
         
    text(0.01, 0.9, ['\bf' level{1}], 'units','normalized')
    if ct == 4
        text(0.8, 0.9, ['\bf CT' ], 'color', col.CT',  'units','normalized')
        text(0.8, 0.8, ['\bf WRF' ], 'color', col.WRF',  'units','normalized')
    end
    for model = {'WRF', 'CT'}
        
    N   =  Var.(model{1}).(level{1}).num(i_run,:);
    G   =  Var.(model{1}).(level{1}).gamma(i_run,:);
    V   =  Var.(model{1}).(level{1}).gamma_var(i_run,:);
    
    N_tot = nansum(N(:,:));
    Gamma = nansum(G(:,:).*N(:,:))./N_tot;

    G_sigma   = sqrt(...
                     1./(N_tot-1) .* ...
                    ( ... 
                     nansum( (N(:,:)-1).*V(:,:) + ...
                     N(:,:).*G(:,:).^2 ) - ...
                     N_tot.*Gamma.^2 ...
                     ));
                 
    % fit exponential model
     start = [Var.distance(end)*2/3 , max(Gamma) 0];
           
     
                            %if strcmp(modeltype{Index}, 'unbounded')
                                start(1) = start(1)/3;
                            %end
                            repeat = true;
                            pp=1;
                            while repeat
                                try
                                    varMod= fitnlm(Var.distance,Gamma,modFun,start.*modifier(pp,:),'Weight',1./G_sigma);
                                    repeat = false;
                                catch
                                    pp = pp+1;
                                    if pp> length(modifier)
                                        error('model did not converge')
                                    end
                                end
                            end
    
    Coefs =  varMod.Coefficients.Estimate;
    Gamma_pred =  predict(varMod,Var.distance','Simultaneous',true);
                 
    line(Var.distance, nanmedian(G),'color',col.(model{1}), 'linestyle','--')
    line(Var.distance, Gamma, 'linestyle','none', 'marker', 'o','color',col.(model{1}),'markerfacecolor',col.(model{1}))
    
    line(Var.distance, Gamma_pred,'color',col.(model{1}), 'linestyle','-', 'linewidth',2)
    if strcmp(model{1},'WRF')
        text(0.02, 0.8, [num2str(round(Coefs(1)*3)) ' ' num2str(Coefs(2),'%5.2f') ' ' num2str(Coefs(3),'%5.2f')], ...
            'units','normalized','color',col.WRF)
    else
        text(0.02, 0.7, [num2str(round(Coefs(1)*3)) ' ' num2str(Coefs(2),'%5.2f') ' ' num2str(Coefs(3),'%5.2f')], ...
            'units','normalized','color',col.CT)
    end
    end
    yl = get(gca,'ylim')
    set(gca,'ylim',[0 yl(2)], 'xlim', [0 500])
ct = ct+1
end
end
set(ha(1:ct-ny-1), 'xticklabel','')
axes(ha(ny+1))
ylabel('\bf \gamma(x) (ppm^2)')
axes(ha(end-1))
%xlabel('\bf Distance x (km)')
print(['Var_' cs{1} '_model.png'],'-dpng','-r300')
end
    %% Daily variogram 

%% Redo figure from 