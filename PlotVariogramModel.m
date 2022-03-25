% Make figure and calculate statistics for the variogram

clear all; close all;
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

%%

n =40;
r = 1.2;
a = 1;
s = a*r.^(0:n-1);
s = [0 unique(floor(s))];
s = s(s<500);
edges = s ;
maxd = edges(end);
distance = (edges(1:end-1)+edges(2:end))/2;

ParentDir = [basepath 'OneDrive - The Pennsylvania State University\Projects\ACT-America\Code\ACT_ModelDataPaper\VarigramBootstrap\']
cd(ParentDir)
DataDir =  ['./OutDir/'] ;

Cases = {'Regions', 'Seasons', 'Airmasses'};
Models = {'CT', 'WRF'};

Regions = {'NEMA','SC','MWe'} ;
Seasons = {'WIN','SPR','SUM', 'FAL'} ;
Airmasses = {'Fair', 'Cold', 'Warm'};

Nme.Regions = Regions; 
Nme.Seasons = Seasons;
Nme.Airmasses = Airmasses;

Levels = {'ABL','LFT','HFT'} ;  
col.WRF = [0.75 0 0];
col.CT  = [0 0 0.75];
col2.CT = [0    0.4470    0.7410];
col2.WRF= [0.8500    0.3250    0.0980];

% for cse = Cases
%     fig = figure('units','inches','position',[0 0 9 9]);
%     ha = tight_subplot(3, length(Nme.(cse{1})),[0.02 0.02],[0.1 0.1], [0.1 0.1]);
%     ct =0 ;   
%    
%     for item = Nme.(cse{1}) 
%         for level = Levels
%             % load data and plot
%             ct = ct+1;
%             axes(ha(ct))
%             hold on;
%             for model = Models    
%                 try
%                     D = load([DataDir 'VarOutput_' model{1} '_' Nme.(cse{1}){1} '_' level{1} '.mat']);
%                 catch
%                     continue
%                 end
%                     % do plot
%                     errorbar(distance, D.Gamma, D.G_sigma, 'color',col.(model{1}), 'o', 'markerfacecolor', col.(model{1}))
%                     
%                     line(distance,D.Gamma_pred.exponential(1,:),'color',col.(model{1}), 'linewidth',2)
%                     clear x y 
%             end
%         end
%     end
% end


%% Fit exponential variogram

plot_eb = false
% Now do this for season and region
%modifier = [1 1 1 ; 0.5 1 1; 2 1 1; 1 0.5 1; 1 2 1 ;0.5 0.5 1; 2 2 1; 0.5 2 1; 2 0.5 1] ;
%modFun = @(b,h)b(3)+b(2)*(1-exp(-h./b(1))) ;
for cis = 1:2
for cs = Cases % loop over cases
    % set figure details
    ny=length(Nme.(cs{1}))  ;  
    figure('units','inches','position',[0 0 10 6])
    ha = tight_subplot(3,ny,[0.02 0.04], [0.05 0.05], [0.05 0.05]);

    ct = 1    ;
    for level = {'ABL','LFT','HFT'} % now loop categories and models
        for run = 1:length(Nme.(cs{1}))    
            axes(ha(ct))
            hold on; box on;
            % add titles and legends
            if ct<=ny
                title(Nme.(cs{1}){run})
            end
            text(0.01, 0.9, ['\bf' level{1}], 'units','normalized')
            if ct == ny
                text(0.8, 0.9, ['\bf CT' ], 'color', col.CT',  'units','normalized')
                text(0.8, 0.8, ['\bf WRF' ], 'color', col.WRF',  'units','normalized')
            end
            
            for model = {'WRF', 'CT'}
                % now do actual plotting 
                try
                    D = load([DataDir 'VarOutput_' model{1} '_' Nme.(cs{1}){run} '_' level{1} '.mat']);
                    disp([DataDir 'VarOutput_' model{1} '_' Nme.(cs{1}){run} '_' level{1} '.mat'])
                catch
                    continue
                end
                 
    %           N   =  Var.(model{1}).(level{1}).num(i_run,:);
                %G   =  Var.(model{1}).(level{1}).gamma(i_run,:);
                %V   =  Var.(model{1}).(level{1}).gamma_var(i_run,:);
    
                N_tot = D.N_tot(1,:) ;
                Gamma = D.Gamma(1,:) ;
                G_sigma = D.G_sigma(1,:) ; 
                
                Gamma_pred = D.Gamma_pred.exponential(1,:);
                
                Gamma_BS = D.Gamma_pred.exponential(2:end,:);
                
                % recalculate bootstrapp following MIT 
                % <<<
                deltastar =bsxfun(@minus, Gamma_BS', Gamma_pred');
                d=quantile(deltastar,[0.1,0.9],2);
                CI_MIT_l = Gamma_pred + d(:,1)';
                CI_MIT_u = Gamma_pred + d(:,2)';
                % >>>
                
                if cis ==1
                    CI_l = [D.CI.exponential(1,:) D.CI.exponential(1,end)];
                    CI_u = [D.CI.exponential(2,:) D.CI.exponential(2,end)];
                elseif cis ==2
                    CI_l = [CI_MIT_l CI_MIT_l(end)];
                    CI_u = [CI_MIT_u CI_MIT_u(end)];
                    Valid.(model{1}).(Nme.(cs{1}){run}).(level{1}) =  Gamma_BS(:,30) >= D.CI.exponential(1,30)/2 & Gamma_BS(:,30) <= D.CI.exponential(2,30)*2
                end
                clear D
                x = [ distance 500];  
                patch([x fliplr(x)], [CI_l fliplr(CI_u)], [x x]*0, ...
                    'facecolor', col2.(model{1}), 'facealpha',0.5)
                
                line(distance, Gamma, 'linestyle','none', 'marker', 'o','color',col.(model{1}),'markerfacecolor',col.(model{1}))
                if plot_eb
                    errorbar(distance, Gamma, G_sigma, 'o','color',col.(model{1}),'markerfacecolor',col.(model{1}))
                end
                line(distance, Gamma_pred,'color',col.(model{1}), 'linestyle','-', 'linewidth',2)
                %if strcmp(model{1},'WRF')
                %    text(0.02, 0.8, [num2str(round(Coefs(1)*3)) ' ' num2str(Coefs(2),'%5.2f') ' ' num2str(Coefs(3),'%5.2f')], ...
                %    'units','normalized','color',col.WRF)
                %else
                %    text(0.02, 0.7, [num2str(round(Coefs(1)*3)) ' ' num2str(Coefs(2),'%5.2f') ' ' num2str(Coefs(3),'%5.2f')], ...
                %    'units','normalized','color',col.CT)
                %end
            end
            yl = get(gca,'ylim');
            set(gca,'ylim',[0 yl(2)], 'xlim', [0 500])
            ct = ct+1;
        end
    end
    
set(ha(1:ct-ny-1), 'xticklabel','')
axes(ha(ny+1))
ylabel('\bf \gamma(x) (ppm^2)')
axes(ha(end-1))
%xlabel('\bf Distance x (km)')
if plot_eb
    print(['Var_' cs{1} 'CI_model_' num2str(cis) '.png'],'-dpng','-r300')
else
   print(['Var_' cs{1} 'CI_model_' num2str(cis) '_noeb.png'],'-dpng','-r300')
end
end
end

%% Histogram of nugget
for cs = Cases % loop over cases
    % set figure details
    ny=length(Nme.(cs{1}))  ;  
    figure('units','inches','position',[0 0 10 6])
    ha = tight_subplot(3,ny,[0.05 0.04], [0.05 0.05], [0.05 0.05]);

    ct = 1    ;
    for level = {'ABL','LFT','HFT'} % now loop categories and models
        for run = 1:length(Nme.(cs{1}))    
            axes(ha(ct))
            hold on; box on;
            % add titles and legends
            if ct<=ny
                title(Nme.(cs{1}){run})
            end
            text(0.01, 0.9, ['\bf' level{1}], 'units','normalized')
            if ct == ny
                text(0.8, 0.9, ['\bf CT' ], 'color', col.CT',  'units','normalized')
                text(0.8, 0.8, ['\bf WRF' ], 'color', col.WRF',  'units','normalized')
            end
            
            for model = {'WRF', 'CT'}
                % now do actual plotting 
                try
                    D = load([DataDir 'VarOutput_' model{1} '_' Nme.(cs{1}){run} '_' level{1} '.mat']);
                    disp([DataDir 'VarOutput_' model{1} '_' Nme.(cs{1}){run} '_' level{1} '.mat'])
                catch
                    continue
                end
                 
    %           N   =  Var.(model{1}).(level{1}).num(i_run,:);
                %G   =  Var.(model{1}).(level{1}).gamma(i_run,:);
                %V   =  Var.(model{1}).(level{1}).gamma_var(i_run,:);
    
                N_tot = D.N_tot(1,:) ;
                Gamma = D.Gamma(1,:) ;
                G_sigma = D.G_sigma(1,:) ; 
                
                Gamma_pred = D.Gamma_pred.exponential(1,:);
                
                Gamma_BS = D.Gamma_pred.exponential(2:end,:);
                
                % >>>
                Nugget = D.Coefs.exponential(2:1001,3) ;
                Nugget_Data = D.Coefs.exponential(1,3) ;

                Nugget = Nugget(Nugget> quantile(Nugget, 0.01) & Nugget < quantile(Nugget, 0.99));
                
                hold on
                [counts, edges] = histcounts(Nugget,  25, 'normalization','probability');
                counts = counts *.98;
                histogram('BinEdges', edges, 'BinCounts', counts,   'facecolor', col2.(model{1}), 'facealpha',0.5)
                line(Nugget_Data.*[1 1], [0 0.3], 'color', col.(model{1}), 'linestyle','--')
            end
            %yl = get(gca,'ylim');
            set(gca,'ylim',[0 0.3])
            ct = ct+1;
        end
    end
    
%set(ha(1:ct-ny-1), 'xticklabel','')
axes(ha(ny+1))
ylabel('\bf Probability')
axes(ha(end-1))
xlabel('\bf Nugget (ppm^2)')
print(['HistNugget_' cs{1} '.png'],'-dpng','-r300')
end

close all
%% Plot x-y plots of ranges vs sill
for cs = Cases % loop over cases
    % set figure details
    ny=length(Nme.(cs{1}))  ;  
    figure('units','inches','position',[0 0 10 6])
    ha = tight_subplot(3,ny,[0.05 0.04], [0.1 0.05], [0.08 0.02]);

    ct = 1    ;
    for level = {'ABL','LFT','HFT'} % now loop categories and models
        for run = 1:length(Nme.(cs{1}))    
            axes(ha(ct))
            hold on; box on;
            % add titles and legends
            if ct<=ny
                title(Nme.(cs{1}){run})
            end
            text(0.01, 0.9, ['\bf' level{1}], 'units','normalized')
            if ct == ny
                text(0.8, 0.9, ['\bf CT' ], 'color', col.CT',  'units','normalized')
                text(0.8, 0.8, ['\bf WRF' ], 'color', col.WRF',  'units','normalized')
            end
            
            for model = {'WRF', 'CT'}
                % now do actual plotting 
                try
                    D = load([DataDir 'VarOutput_' model{1} '_' Nme.(cs{1}){run} '_' level{1} '.mat']);
                    disp([DataDir 'VarOutput_' model{1} '_' Nme.(cs{1}){run} '_' level{1} '.mat'])
                catch
                    continue
                end
                 
    %           N   =  Var.(model{1}).(level{1}).num(i_run,:);
                %G   =  Var.(model{1}).(level{1}).gamma(i_run,:);
                %V   =  Var.(model{1}).(level{1}).gamma_var(i_run,:);
    
                N_tot = D.N_tot(1,:) ;
                Gamma = D.Gamma(1,:) ;
                G_sigma = D.G_sigma(1,:) ; 
                
                Gamma_pred = D.Gamma_pred.exponential(1,:);
                
                Gamma_BS = D.Gamma_pred.exponential(2:end,:);
                
                % >>>
                Sill = D.Coefs.exponential(2:1001,2) ;
                Range = D.Coefs.exponential(2:1001,1)*3 ;

                
                Sill_Data.(model{1}) = D.Coefs.exponential(1,2) ;
                Range_Data.(model{1}) = D.Coefs.exponential(1,1)*3 ;
                Nugget_Data.(model{1}) = D.Coefs.exponential(1,3)*3 ;
                
                Sills.(model{1}).(Nme.(cs{1}))(ct) = D.Coefs.exponential(1,2) ;
                Ranges.(model{1}).(Nme.(cs{1}))(ct) = D.Coefs.exponential(1,1)*3 ;
                Nuggets.(model{1}).(Nme.(cs{1}))(ct) = D.Coefs.exponential(1,3) ;
                
                hold on
                
                scatter(Range,Sill,'.','markeredgecolor',col2.(model{1}));

                
                y_dummy = quantile(Sill,[0.25 0.75]);
                x_dummy = quantile(Range,[0.25 0.75]);
                if strcmp(model{1}, 'WRF')
                    xl(1) = x_dummy(1)/6;
                    xl(2) = x_dummy(2)*6;
                    yl(1) = y_dummy(1)/6;
                    yl(2) = y_dummy(2)*6;
                else
                    xl(1) = min(xl(1), x_dummy(1)/6);
                    xl(2) = max(xl(2), x_dummy(2)*6);
                    yl(1) = min(yl(1), y_dummy(1)/6);
                    yl(2) = max(yl(2), y_dummy(2)*6);
                end
            end
            for model = {'WRF', 'CT'} 
                scatter(Range_Data.(model{1}), Sill_Data.(model{1}),'o','filled', 'markeredgecolor',col.(model{1}), 'markerfacecolor',col.(model{1}))
            end
            %yl = get(gca,'ylim');
            set(gca,'ylim',yl, 'xlim',xl)
            ct = ct+1;
        end
    end
    
%set(ha(1:ct-ny-1), 'xticklabel','')
axes(ha(ny+1))
ylabel('\bf Sill (ppm^2)')
axes(ha(end-1))
xlabel('\bf Range x (km)')
print(['Scatter_RangeSill' cs{1} '.png'],'-dpng','-r300')
end

