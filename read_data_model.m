clear all

fc1 = brewermap(5,'Set1');
fc2 = brighten(fc1,.8);
M_1 = fc1(3,:);
M_2 = fc1(1,:);
M_3 = fc1(2,:);
M_4 = fc1(4,:);
M_5 = fc1(5,:);
m_1 = fc2(1,:);
m_2 = fc2(2,:);
m_3 = fc2(3,:);

fold  = '/Volumes/D_Visioni2/GeoMIP/'; %Change with directory where you downloaded data
n_exp = {'G6sulfur','G6solar','ssp585','ssp245'}; %Experiments name
n_mod = {'CNRM','IPSL','CESM2','UKESM1','CanESM5','NESM3','MPI-ESM1-2-LR'}; %Available models

%%% Set up of available ensemble members %%%
n_ens{1,1}=[1,2,3];n_ens{1,2}=[1];n_ens{1,3}=n_ens{1,1};n_ens{1,4}=n_ens{1,1};
n_ens{2,1}=[1];n_ens{2,2}=n_ens{2,1};n_ens{2,3}=n_ens{2,1};n_ens{2,4}=n_ens{2,1};
n_ens{3,1}=[1,2];n_ens{3,2}=[1,2];n_ens{3,3}=[1];n_ens{3,4}=[1];
n_ens{4,1}=[1,4,8];n_ens{4,2}=[1,4,8];n_ens{4,3}=[1,4,8];n_ens{4,4}=[1,4,8];
n_ens{5,1}=[1,2,3];n_ens{5,2}=n_ens{5,1};n_ens{5,3}=n_ens{5,1};n_ens{5,4}=n_ens{5,1};
n_ens{6,1}=[1,2];n_ens{6,2}=n_ens{6,1};n_ens{6,3}=n_ens{6,1};n_ens{6,4}=n_ens{6,1};
n_ens{7,1}=[2,3];n_ens{7,2}=n_ens{7,1};n_ens{7,3}=n_ens{7,1};n_ens{7,4}=n_ens{7,1};

%if the beginning year is different than 2015 (some GeoMIP models, not all)
y_start = [0,5,0,5,0,0,0; 0,5,0,5,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0];

%variable names
var1 = {'rsdt','rsut','rlut','rsutcs','rlutcs','ts','od550aer','pr'};
%variable explanation
var_l = {'Down SW TOA','Up SW TOA','Up LW TOA','Up SW CS TOA','Up LW CS TOA','TS','AOD','precipitation'};
%in case a change of unit is necessary
corr = [1,1,1,1,1,1,1,86400];

%select variable to analyse
iv = 8;

cc = 0;

figure(1)
fprintf('Model\tens size\t \n')
set(gcf, 'Position',  [200, 200, 1000, 1000])
for im=1:4
    fprintf([n_mod{im} '\t(' ')\t'])
    subplot(4,2,im)
    
    for ie=1:4

        if iv==7 & im==2; var1{iv} = 'od550aerso'; end         
        if iv==7 & im~=2; var1{iv} = 'od550aer'; end
        
        fprintf([n_exp{ie} ' '])
        nne = n_ens{im,ie};
        for i=1:length(nne)
            list = dir([fold n_exp{ie} '/' n_mod{im} '/' var1{iv} '_*' n_mod{im} ...
                '*' n_exp{ie} '_r' num2str(nne(i)) '*.nc']);
            if isempty(list)
                x = 0;
                if i==1, fprintf('(N) '); end
            else
                if i==1, fprintf('(Y) '); end
                for j=1:length(list)
                    x2 = ncread([list(j).folder '/' list(j).name],var1{iv})*corr(iv);              
                    if length(list)>1 && j>1
                        x = cat(3,x,x2);
                    else
                        x = x2;
                    end
                end
            end
            if (im==3 & i == 1 & ie == 1 & isempty(list)==0)
                xsize = size(x);
                xlngth = NaN(xsize(1),xsize(2),xsize(3)+5*12);
                xlngth(:,:,5*12+1:size(xlngth,3)) = x;
                clear x
                x = xlngth;
            end
            if (im==3 & i == 1 & ie == 2 & isempty(list)==0)
                xsize = size(x);
                xlngth = NaN(xsize(1),xsize(2),xsize(3)+4*12);
                xlngth(:,:,4*12+1:size(xlngth,3)) = x;
                clear x
                x = xlngth;
            end
            if length(nne)>1
                x_e(:,:,:,i) = x;
            end
        end
        if length(nne)>1
            if i==1, x_e = zeros([size(x),length(nne)]); end
            x = x_e;
        end
        if isempty(list), else
            lat = ncread([list(1).folder '/' list(1).name],'lat');
            lon = ncread([list(1).folder '/' list(1).name],'lon');
            time = ncread([list(1).folder '/' list(1).name],'time');
        end
        
        %%
        
        if isempty(list)
        else
            cc = cc+1;
            n_exp2{cc} = n_exp{ie}; %skip over missing experiments, if there are
            ww = cos(lat/180*pi);
            x_zon = squeeze(nanmean(x,1)); %average over longitude
            siz = size(x);
            if length(nne)>1
                x = squeeze(nanmean(x,4)); %do ensemble mean for saving maps
            end
            x_map_20{ie,im} = squeeze(nanmean(x(:,:,((siz(3)/12-20)*12-1):siz(3)),3)); %save map data for last 2 decades
            LAT{im} = lat;
            LON{im} = lon;
            
            xs = size(x_zon);
            
            x_zon = reshape(x_zon,xs(1),12,xs(2)/12,length(nne)); %divide by years
            
            x_y = squeeze(nanmean(x_zon,2)); % get annual mean
            
            if iv==7 % for AOD, take mean of last 10 years for lat distr
                x_ye = squeeze(nanmean(x_y,3));
                x_ym = squeeze(nanmean(x_ye(:,...
                    (size(x_ye,2)-10):size(x_ye,2)),2));
                x_lat{im,ie} = x_ym;
            end
            
            x_mon = squeeze(nanmean(x_zon(:,:,(xs(2)/12-19):xs(2)/12,:),3)); % get seasonal mean, last 20 y
            
            x_glob = squeeze(nansum(x_y.*ww)/sum(ww)); %do global average

            x_glob_mon{im,ie} = squeeze(nansum(x_mon.*ww)/sum(ww)); %do global average
            
            x_glob_mon_lat{im,ie} = squeeze(nanmean(x_mon,3));
            
  
                
            if length(nne)>1 %do ensemble average and plot single runs
                if ie<=2 & xs(2)/12>81
                    x_glob(1:xs(2)/12-81,:) = NaN;
                end  
                x_glob_m = nanmean(x_glob,2);
                for imi=1:length(nne)
                    plot((1:(xs(2)/12))+y_start(ie,im)+2014,x_glob(:,imi),'Linewidth',1,'Color',fc2(ie,:));
                    hold on
                end
            else
                if ie<=2 & xs(2)/12>81
                    x_glob(1:xs(2)/12-81) = NaN;
                end                  
                x_glob_m = x_glob;
            end
            
            x_glob_m_save(im,ie,:) = x_glob_m((6-y_start(ie,im)):length(x_glob_m));
            x_dec = reshape(x_glob_m((7-y_start(ie,im)):length(x_glob_m)),10,8);
            x_decade(ie,im,:) = nanmean(x_dec,1);
            
            p(cc) = plot((1:(xs(2)/12))+y_start(ie,im)+2014,x_glob_m,'Linewidth',2,'Color',fc1(ie,:));
            hold on
            
            %set(gca,'XTick',[2020:10:2100],'YTick',[0:360],'FontSize',18,'FontWeight','Bold')
        end
        
        clear x x_e
        
    end
    axis([2015 2100 -Inf Inf])
    title([var_l{iv} ' - ' n_mod{im}],'FontSize',18,'FontWeight','Bold')
    set(gca,'Linewidth',2,'FontSize',16,'FontWeight','Bold')
    %    legend(p(:),n_exp2{:},'Location','best')
    
    clear p n_exp2 xt
    cc = 0;
    fprintf(['\n'])
end

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['figures/singlemodels_' var1{iv} '.eps'])
save([var1{iv} '_timeseries.mat'],'x_glob_m_save')
if iv == 7
    save([var1{iv} '_latdistr.mat'],'x_lat')
end
%%
figure(2)
set(gcf, 'Position',  [200, 200, 900, 700])
for ie=1:4
    for im=1:4
        if mean(x_glob_m_save(im,ie,:),3)==0, x_glob_m_save(im,ie,:)=NaN; end
    end
    hold on
    box on
    std_glob_m = squeeze(nanstd(x_glob_m_save(:,ie,:),0,1))/sqrt(6);
    avg_glob_m = squeeze(nanmean(x_glob_m_save(:,ie,:),1));
    Y_std = [avg_glob_m'-std_glob_m',fliplr(avg_glob_m'+std_glob_m')];
    X = [2020:1:2100,fliplr(2020:1:2100)];
    fill(X,Y_std,fc1(ie,:),'LineStyle','none','facealpha',.2)
    p(ie) = plot(2020:1:2100,squeeze(nanmean(x_glob_m_save(:,ie,:),1)),...
        'Linewidth',2,'Color',fc1(ie,:));
end

for ie=1:4
    plot(2020:1:2100,squeeze(nanmean(x_glob_m_save(:,ie,:),1)),...
        'Linewidth',2,'Color',fc1(ie,:));
    hold on
end
legend(p(:),n_exp{:},'Location','best')
set(gca,'Linewidth',2,'FontSize',16,'FontWeight','Bold')
set(gca,'XTick',[2020:10:2100],'YTick',[0:400],'FontSize',18,'FontWeight','Bold')
axis([2020 2100 -Inf Inf])
title([var_l{iv} ' - model average'],'FontSize',18,'FontWeight','Bold')
set(gcf,'renderer','painters')
print(gcf,'-depsc2',['figures/allmodels_' var1{iv} '.eps'])