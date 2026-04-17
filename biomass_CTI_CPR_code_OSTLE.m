% biomass_CTI.m
% Written by Clare Ostle 2026

% Analysis of North Atlantic copepod assemblage biomass, abundance, and
% Community Temperature Index (CTI) based on Beaugrand (2002) groupings.
%
% Assemblages: warm-temperate, pseudo-oceanic temperate, cold-temperate, subarctic
% Period: 1958-2022, from the Continuous Plankton Recorder (CPR) Survey
%

clearvars

addpath /Users/Clare/Documents/MATLAB/CPR/data/functions
addpath /Users/Clare/Documents/MATLAB/File_Exchange/colormaps
addpath /Users/Clare/Documents/MATLAB/File_Exchange/regress
addpath /Users/Clare/Documents/MATLAB/functions/autocorrelation/

%% Parameters
start_year = 1958;
end_year   = 2022;

max_lon = 20;  min_lon = -50;
max_lat = 70;  min_lat = 40;

lat = -89.5:1:89.5;
lon = -179.5:1:179.5;

ln = (lon <= max_lon) & (lon >= min_lon);
ln_NAtl = lon(ln);
lt = (lat <= max_lat) & (lat >= min_lat);
lt_NAtl = lat(lt);

%% Load CPR database
load /Users/Clare/Documents/MATLAB/CPR/data/database_extracts/NorthAtlantic_CPR_1958_2025_20251128.mat
Database_extract = taxa_TT_pivot;

%% Define copepod assemblage IDs (Beaugrand 2002)
all_NAtlantic_copepod_IDs = [53; 419; 51; 45; 9; 47; 59; 60; 58; 412; ...
    411; 10588; 408; 380; 313; 10603; 10604; 304; 48; 50; 10596; 6; ...
    61; 41; 370; 57; 5; 55; 407; 321; 10605; 40; 44; 56; 42];

warm_temperate_IDs         = [53; 419; 51; 45; 9; 47; 59; 60; 58; 412; 411; 10588; 408; 380; 313];
southern_shelf_edge_IDs    = [10603; 10604; 304; 48];
pseudo_oceanic_temperate_IDs = [50; 10596; 6; 61; 41];
cold_temperate_IDs         = [370; 57; 5; 55];
subarctic_IDs              = [407; 321; 10605; 40];
arctic_IDs                 = [44; 56; 42];

% Extract copepod subset for export
ind = ismember(taxa_names.accepted_id, all_NAtlantic_copepod_IDs);
Copepod_Taxa_names = taxa_names(ind,:);
latlon = logical(ones(3,1));
ind_counts = [latlon; ind];
Copepod_counts = taxa_TT_pivot(:, ind_counts);

%% Compare STIs across assemblages
load '/Users/Clare/Documents/MATLAB/CPR/data/monthly_grids/NAtlantic_STIs.mat';

assemblage_id_sets = {warm_temperate_IDs, southern_shelf_edge_IDs, ...
    pseudo_oceanic_temperate_IDs, cold_temperate_IDs, subarctic_IDs, arctic_IDs};
assemblage_sti_names = {'Warm','Southern','Pseudo','Cold','Subarctic','Arctic'};

mean_STIs = nan(1, 6);
for a = 1:6
    idx = ismember(taxa_names.accepted_id, assemblage_id_sets{a});
    mean_STIs(a) = mean(STImedian(idx));
end
mean_warm_STI      = mean_STIs(1);
mean_southern_STI  = mean_STIs(2);
mean_pseudo_STI    = mean_STIs(3);
mean_cold_STI      = mean_STIs(4);
mean_subarctic_STI = mean_STIs(5);
mean_Arctic_STI    = mean_STIs(6);

%% Data coverage mask (require >=75% valid years per cell)
load /Users/Clare/Documents/MATLAB/CPR/ALLNAtlantic/biomass/cpr_annual_nans_from_obana_sus_2022.mat

n_years = (end_year - start_year) + 1;
max_missing_years = floor(0.25 * n_years);
cpr_not_enough = cpr_annual_nans > max_missing_years;

%% Mapping parameters
max_lon = 10;  min_lon = -45;
max_lat = 65;  min_lat = 43;

lon_ind = lon <= max_lon & lon >= min_lon;
lat_ind = lat <= max_lat & lat >= min_lat;

% Assemblage file/name lookup used throughout
assemblage_files = {'monthly_cpr_Warm-temperate copepods 2022_gridded.mat', ...
                    'monthly_cpr_Pseudo-oceanic temperate copepods 2022_gridded.mat', ...
                    'monthly_cpr_Cold-temperate copepods 2022_gridded.mat', ...
                    'monthly_cpr_Subarctic copepods 2022_gridded.mat'};
assemblage_labels = {'Warm-temperate', 'Pseudo-oceanic temperate', ...
                     'Cold-temperate', 'Subarctic copepods'};

%% Plot mean assemblage abundance (Figure 55)
figure(55);
for group = 1:4
    load(['/Users/Clare/Documents/MATLAB/CPR/data/monthly_grids/' assemblage_files{group}])

    var2plot = squeeze(nanmean(monthly_cpr(:, lat_ind, lon_ind), 1));
    subplot(2,2,group);
    h = worldmap([min_lat max_lat], [min_lon max_lon]);
    pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, var2plot, 'FaceColor','flat');
    p = findobj(h,'type','patch');
    set(p,'FaceColor',[0.8 0.8 0.8]);
    geoshow('landareas.shp', 'FaceColor', 'k')

    hcb = colorbar('westoutside');
    hcb.Label.String = 'Counts per 3m^3 (log_{10}(x+1))';
    hcb.Label.FontSize = 14;
    caxis([0 1]);
    drawnow;
    hcb.Position(1) = hcb.Position(1) - 0.05;

    name_clean = strtrim(regexprep(assemblage_labels{group}, '\s*\d+\s*', ' '));
    title(name_clean)
end

%% Sampling effort
figure
effort_reg = cpr_sampleno(:, lat_ind, lon_ind);
effort_total = squeeze(nansum(effort_reg, 1));

ind_no_samples = effort_total == 0;
cpr_not_enough(ind_no_samples) = 1;

subplot(2,1,1)
h = worldmap([min_lat max_lat], [min_lon max_lon]);
pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, log10(effort_total + 1), 'FaceColor','flat');
p = findobj(h, 'type','patch');
set(p, 'FaceColor',[0.8 0.8 0.8]);
geoshow('landareas.shp', 'FaceColor','k');
gridm('off');
cbh = colorbar('westoutside');
cbh.Label.String = 'log_{10}(sample count + 1)';
cbh.Label.FontSize = 11;
drawnow;
cbh.Position(1) = cbh.Position(1) - 0.05;
title('CPR sampling effort', 'FontSize',13, 'FontWeight','bold');

effort_monthly = squeeze(nansum(nansum(effort_reg, 2), 3));
effort_ann = sum(reshape(effort_monthly, 12, []), 1)';
date_v_eff = datevec(time);
uniq_yrs = unique(date_v_eff(:,1));
time_eff = double(uniq_yrs(1:numel(effort_ann)));

subplot(2,1,2)
bar(time_eff, effort_ann, 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none', 'FaceAlpha',0.8);
ylabel('Annual sample count', 'FontSize',11, 'FontWeight','bold');
xlabel('Year', 'FontSize',11);
title('CPR sampling effort', 'FontSize',13, 'FontWeight','bold');
set(gca, 'FontSize',10);
xlim([1958 2025]); ylim([0 max(effort_ann)*1.1]); box on

%% Objective-mapped decadal abundance maps
figure;
ff = 1;

for group = 1:4
    [cpr_monthly_obj, ~, cpr_obj_error] = obj_map_CPR(assemblage_files{group}, ...
        max_lon, min_lon, max_lat, min_lat);

    load(['/Users/Clare/Documents/MATLAB/CPR/data/monthly_grids/' assemblage_files{group}])
    date_v = datevec(time);
    start = 1958;

    for yy = 1:5
        year_ind = date_v(:,1) >= start & date_v(:,1) < start + 13;
        var2plot = squeeze(nanmean(cpr_monthly_obj(year_ind,:,:), 1));
        var2plot(cpr_not_enough) = NaN;
        mask = ~isnan(var2plot);
        Cplot = var2plot;
        Cplot(~mask) = 0;

        subplot(5, 4, ff);
        h = worldmap([min_lat max_lat], [min_lon max_lon]);
        hs = pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, Cplot);
        mlabel off; plabel off
        set(hs, 'EdgeColor','none', 'FaceAlpha','texturemap', ...
            'AlphaData',double(mask), 'AlphaDataMapping','none');
        p = findobj(h,'type','patch');
        set(p,'FaceColor',[0.8 0.8 0.8]);
        geoshow('landareas.shp', 'FaceColor', [.5 .5 .5])

        name_clean = strtrim(regexprep(assemblage_labels{group}, '\s*\d+\s*', ' '));
        if ff <= 4
            text(-1800000, 8300000, name_clean, 'FontWeight','Bold', 'FontSize',14);
        end

        if start == 2010
            title([num2str(start) ' - ' num2str(date_v(end,1))], 'fontsize',14);
        else
            title([num2str(start) ' - ' num2str(start+12)], 'fontsize',14);
        end
        caxis([0 1]);
        start = start + 13;
        ff = ff + 4;
    end

    % Reset ff for next group
    if group < 4
        ff = group + 1;
    end
end

hcb = colorbar;
set(hcb, 'Position', [.93 .2 0.03 0.6]);
hcb.Label.String = 'Counts per 3m^3 (log_{10}(x+1))';
hcb.Label.FontSize = 14;
caxis([0 1]);

%% Abundance trend maps
start_year = 1958;
figure;

for group = 1:4
    load(['/Users/Clare/Documents/MATLAB/CPR/data/monthly_grids/' assemblage_files{group}])

    [cpr_annual_obana_trend, ~, ~, ~, ~, ~] = obj_map_TREND_CPR( ...
        assemblage_files{group}, max_lon, min_lon, max_lat, min_lat, start_year, end_year);

    subplot(2,2,group)
    h = worldmap([min_lat max_lat], [min_lon max_lon]);
    pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, cpr_annual_obana_trend);
    hold on
    colormap(gca, b2r(min(cpr_annual_obana_trend(:)), max(cpr_annual_obana_trend(:))));
    cbh = colorbar('Location','westoutside');
    p = findobj(h,'type','patch');
    set(p,'FaceColor',[0.8 0.8 0.8]);
    geoshow('landareas.shp', 'FaceColor', 'k')
    gridm('off');

    name_clean = regexprep(assemblage_labels{group}, '\d[0-9_]+\d', '');
    title(['Abundance of ' name_clean])
    set(get(cbh,'ylabel'), 'string', 'log_{10}(organisms per 3m^3 +1).yr^-^1', 'fontsize',14);
    cbh.Ruler.Exponent = 0;
end

%% Load pre-computed biomass outputs
load /Users/Clare/Documents/MATLAB/CPR/ALLNAtlantic/biomass/biomass_gridded_copepod_assemblages.mat

%% Biomass trend maps (5-degree grid vs objective mapping, Figure 100)
figure(100);

assemblage_outputs_list = {z_warm_temperate, z_temperate_output, z_cold_output, z_subarctic_output};
assemblage_labels_short = {'Warm-temperate','Pseudo-oceanic temperate','Cold-temperate','Subarctic'};

for group = 1:4
    z = assemblage_outputs_list{group};

    % 5-degree grid trends
    latCenters = z.monthly_cpr_5deg.latCenters;
    lonCenters = z.monthly_cpr_5deg.lonCenters;
    monthTimes = z.monthly_cpr_5deg.monthTimes;
    X = z.monthly_cpr_5deg.monthlyMean;

    yrs = year(monthTimes);
    uniqY = unique(yrs);
    nLat = size(X,1); nLon = size(X,2);
    annualMean = nan(nLat, nLon, numel(uniqY));
    for yi = 1:numel(uniqY)
        annualMean(:,:,yi) = nanmean(X(:,:,yrs == uniqY(yi)), 3);
    end

    cpr_annual_obana_trend = nan(nLat, nLon);
    for iLat = 1:nLat
        for iLon = 1:nLon
            y = squeeze(annualMean(iLat,iLon,:));
            x = double(uniqY);
            ok = ~isnan(y);
            if sum(ok) >= 5
                p = polyfit(x(ok), double(y(ok)), 1);
                cpr_annual_obana_trend(iLat,iLon) = p(1);
            end
        end
    end

    % 5-degree map
    subplot(2,4,group)
    h = worldmap([min_lat max_lat], [min_lon max_lon]);
    pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, cpr_annual_obana_trend);
    hold on
    colormap(gca, b2r(min(cpr_annual_obana_trend(:)), max(cpr_annual_obana_trend(:))));
    cbh = colorbar;
    p = findobj(h,'type','patch');
    set(p,'FaceColor',[0.8 0.8 0.8]);
    geoshow('landareas.shp', 'FaceColor', [.5 .5 .5])
    gridm('off');
    title(sprintf('%s Carbon Mass %d-%d 5deg', assemblage_labels_short{group}, start_year, end_year), 'FontSize',14);
    set(gca,'fontsize',14);

    % Objective-mapped (masked) map
    var2plot = z.cpr_annual_obana_trend;
    var2plot(cpr_not_enough) = NaN;
    mask = ~isnan(var2plot);
    Cplot = var2plot;
    Cplot(~mask) = 0;

    subplot(2,4,group+4)
    h = worldmap([min_lat max_lat], [min_lon max_lon]);
    hs = pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, Cplot);
    set(hs, 'EdgeColor','none', 'FaceAlpha','texturemap', ...
        'AlphaData',double(mask), 'AlphaDataMapping','none');
    hold on
    colormap(gca, b2r(min(z.cpr_annual_obana_trend(:)), max(z.cpr_annual_obana_trend(:))));
    cbh = colorbar;
    p = findobj(h,'type','patch');
    set(p,'FaceColor',[0.8 0.8 0.8]);
    geoshow('landareas.shp', 'FaceColor', [.5 .5 .5])
    gridm('off');
    title(sprintf('%s Carbon Mass %d-%d ObjMap', assemblage_labels_short{group}, start_year, end_year), 'FontSize',14);
    set(gca,'fontsize',14);
end

%% Biomass trend maps per assemblage (Figure 200, masked + centred colourbar)
figure(200);

for group = 1:4
    z = assemblage_outputs_list{group};

    var2plot = z.cpr_annual_obana_trend;
    var2plot(cpr_not_enough) = NaN;
    mask = ~isnan(var2plot);
    Cplot = var2plot;
    Cplot(~mask) = 0;

    subplot(2,2,group)
    h = worldmap([min_lat max_lat], [min_lon max_lon]);
    hs = pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, Cplot);
    set(hs, 'EdgeColor','none', 'FaceAlpha','texturemap', ...
        'AlphaData',double(mask), 'AlphaDataMapping','none');
    hold on

    clim_val = max(abs(Cplot(mask)));
    caxis(gca, [-clim_val clim_val]);
    colormap(gca, b2r(-clim_val, clim_val));
    cbh = colorbar('westoutside');
    p = findobj(h,'type','patch');
    set(p,'FaceColor',[0.8 0.8 0.8]);
    geoshow('landareas.shp', 'FaceColor', 'k')
    gridm('off');
    title(sprintf('%s Carbon Mass %d-%d', assemblage_labels_short{group}, start_year, end_year), 'FontSize',14);
    set(get(cbh,'ylabel'), 'string', 'log_{10}(\mugC per 3m^3+1).yr^-^1', 'fontsize',14);
    cbh.Ruler.Exponent = 0;
    set(gca,'fontsize',14);
    drawnow;
    cbh.Position(1) = cbh.Position(1) - 0.05;
end

%% CTI gridding
load /Users/Clare/Documents/MATLAB/CPR/ALLNAtlantic/biomass/cti_outputs.mat

cti_out_ALL_NAtl_copepods = cti_grid_trend(taxa_TT_pivot, taxa_names, ...
    all_NAtlantic_copepod_IDs, [], max_lon, min_lon, max_lat, min_lat, ...
    start_year, end_year, 'North Atlantic Copepods', 1);

%% Combined abundance, biomass, CTI trend maps
[cpr_annual_obana_trend, ~, ~, ~, ~, ~] = obj_map_TREND_CPR( ...
    'monthly_cpr_NAtlantic_copepods 2022_gridded.mat', ...
    max_lon, min_lon, max_lat, min_lat, start_year, end_year);
cpr_annual_obana_trend(cpr_not_enough) = NaN;

figure
% Abundance
subplot(3,1,1)
h = worldmap([min_lat max_lat], [min_lon max_lon]);
pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, cpr_annual_obana_trend);
hold on
colormap(gca, b2r(min(cpr_annual_obana_trend(:)), max(cpr_annual_obana_trend(:))));
cbh = colorbar;
p = findobj(h,'type','patch');
set(p,'FaceColor',[0.8 0.8 0.8]);
geoshow('landareas.shp', 'FaceColor', 'k'); gridm('off');
title('Copepod Abundance');
set(get(cbh,'ylabel'), 'string', 'log_{10}(organisms per 3m^3 +1).yr^-^1', 'fontsize',14);
cbh.Ruler.Exponent = 0;

% Biomass
z_copepod_biomass_output.cpr_annual_obana_trend(cpr_not_enough) = NaN;

subplot(3,1,2)
h = worldmap([min_lat max_lat], [min_lon max_lon]);
pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, z_copepod_biomass_output.cpr_annual_obana_trend);
hold on
colormap(gca, b2r(min(z_copepod_biomass_output.cpr_annual_obana_trend(:)), ...
    max(z_copepod_biomass_output.cpr_annual_obana_trend(:))));
cbh = colorbar;
p = findobj(h,'type','patch');
set(p,'FaceColor',[0.8 0.8 0.8]);
geoshow('landareas.shp', 'FaceColor', 'k'); gridm('off');
title('Copepod Biomass');
set(get(cbh,'ylabel'), 'string', 'log_{10}(\mugC per 3m^3+1).yr^-^1', 'fontsize',14);
cbh.Ruler.Exponent = 0;

% CTI
cti_out_ALL_NAtl_copepods.cpr_annual_trend(cpr_not_enough) = NaN;

subplot(3,1,3)
h = worldmap([min_lat max_lat], [min_lon max_lon]);
pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, cti_out_ALL_NAtl_copepods.cpr_annual_trend);
hold on
colormap(gca, b2r(min(cti_out_ALL_NAtl_copepods.cpr_annual_trend(:)), ...
    max(cti_out_ALL_NAtl_copepods.cpr_annual_trend(:))));
cbh = colorbar;
p = findobj(h,'type','patch');
set(p,'FaceColor',[0.8 0.8 0.8]);
geoshow('landareas.shp', 'FaceColor', 'k'); gridm('off');
title('Copepod CTI');
set(get(cbh,'ylabel'), 'string', '^{\circ}C.yr^{-1}', 'fontsize',14);
cbh.Ruler.Exponent = 0;

%% -----------------------------------------------------------------------
%  STATISTICAL ANALYSIS
%  Trend tests with autocorrelation correction (mod_glover_AC)
%  GLS AR(1) regression for CTI~SST
%% -----------------------------------------------------------------------

% Biomass trend
ln = (lon <= max_lon) & (lon >= min_lon);
lt = (lat <= max_lat) & (lat >= min_lat);

monthly_bio = z_copepod_biomass_output.monthly_cpr(:, lt, ln);
var_bio     = squeeze(nanmean(nanmean(monthly_bio, 2), 3));
var_bio_ann = squeeze(nanmean(reshape(var_bio, 12, []), 1))';
time_bio    = uniqY;

stats_biomass = calc_trend_AC(double(time_bio), var_bio_ann);

fprintf('\n=== BIOMASS TREND ===\n');
fprintf('Sen slope:  %.4f log10(ugC/3m3) yr-1\n', stats_biomass.slope_Sen);
fprintf('Per decade: %.3f\n', stats_biomass.slope_Sen*10);
fprintf('Pct change: %.1f%%\n', stats_biomass.pct_change);
fprintf('AC-corrected p: %.4f  (Neff=%.1f)\n', stats_biomass.p_AC, stats_biomass.N_eff);
fprintf('Mann-Kendall z: %.3f, p = %.4f\n', stats_biomass.MK_z, stats_biomass.MK_p);

% CTI trend
monthly_cti_reg = cti_out_ALL_NAtl_copepods.monthly_cti.monthly_cpr(:, lt, ln);
var_cti         = squeeze(nanmean(nanmean(monthly_cti_reg, 2), 3));
var_cti_ann     = squeeze(nanmean(reshape(var_cti, 12, []), 1))';
time_cti        = time_bio;

stats_cti = calc_trend_AC(double(time_cti), var_cti_ann);

fprintf('\n=== CTI TREND ===\n');
fprintf('Sen slope:  %.4f degC yr-1\n', stats_cti.slope_Sen);
fprintf('Per decade: %.3f\n', stats_cti.slope_Sen*10);
fprintf('AC-corrected p: %.4f  (Neff=%.1f)\n', stats_cti.p_AC, stats_cti.N_eff);
fprintf('Mann-Kendall z: %.3f, p = %.4f\n', stats_cti.MK_z, stats_cti.MK_p);

% Assemblage carbon mass trends
load /Users/Clare/Documents/MATLAB/CPR/ALLNAtlantic/biomass/biomass_gridded_copepod_assemblages.mat

assemblage_names = {'Warm-temperate','Pseudo-oceanic temperate','Cold-temperate','Subarctic'};
assemblage_outputs = {z_warm_temperate, z_temperate_output, z_cold_output, z_subarctic_output};

fprintf('\n=== ASSEMBLAGE CARBON MASS TRENDS ===\n');
fprintf('%-30s  %8s  %8s  %6s  %6s\n', 'Assemblage','Sen/yr','Sen/dec','AC_p','MK_p');

stats_assemblage = struct();
for g = 1:4
    monthly_int = assemblage_outputs{g}.monthly_cpr(:, lt, ln);
    vv     = squeeze(nanmean(nanmean(monthly_int, 2), 3));
    vv_ann = squeeze(nanmean(reshape(vv, 12, []), 1))';

    st = calc_trend_AC(double(uniqY), vv_ann);
    stats_assemblage(g).name  = assemblage_names{g};
    stats_assemblage(g).stats = st;

    fprintf('%-30s  %+8.4f  %+8.3f  %6.4f  %6.4f\n', ...
        assemblage_names{g}, st.slope_Sen, st.slope_Sen*10, st.p_AC, st.MK_p);
end

%% CTI ~ SST relationship (GLS with AR1)
load('/Users/Clare/Documents/MATLAB/DATA/ICOADS/month_mean_sst_1deg.mat');

lon_ind_sst = (lon_sst >= min_lon) & (lon_sst <= max_lon);
lat_ind_sst = (lat_sst >= min_lat) & (lat_sst <= max_lat);
sst_reg     = sst(:, lat_ind_sst, lon_ind_sst);
sst_ts      = squeeze(nanmean(nanmean(sst_reg, 2), 3));

time_sst_dv = datevec(time_sst);
year_sst    = time_sst_dv(:,1);
period_ind  = (year_sst >= start_year) & (year_sst <= end_year);
sst_ts      = sst_ts(period_ind);
year_sst    = year_sst(period_ind);

uniq_years = unique(year_sst);
n_yrs      = numel(uniq_years);
sst_ann    = nan(n_yrs, 1);
for yi = 1:n_yrs
    sst_ann(yi) = nanmean(sst_ts(year_sst == uniq_years(yi)));
end
time_sst = uniq_years;

fprintf('SST: %d years (%d-%d)\n', n_yrs, time_sst(1), time_sst(end));

% Align CTI and SST on common years
[~, ia, ib]  = intersect(time_cti, time_sst);
cti_aligned  = var_cti_ann(ia);
sst_aligned  = sst_ann(ib);

% OLS reference
X_ols = [ones(size(sst_aligned)), sst_aligned];
[b_ols, ~, ~, ~, stats_ols] = regress(cti_aligned, X_ols);
fprintf('\nCTI ~ SST (OLS): slope=%.4f, R2=%.3f, p=%.4f\n', b_ols(2), stats_ols(1), stats_ols(3));

% AR1 from OLS residuals
resid_ols = cti_aligned - X_ols * b_ols;
corr_matrix = corrcoef(resid_ols(1:end-1), resid_ols(2:end));
rho = corr_matrix(1,2);
fprintf('AR1 rho = %.3f\n', rho);

% Prais-Winsten GLS
n_t = numel(cti_aligned);
W = eye(n_t);
for i = 2:n_t
    W(i,i-1) = -rho;
end
W(1,1) = sqrt(1 - rho^2);

cti_gls = W * cti_aligned;
X_gls   = W * X_ols;
[b_gls, ~, ~, ~, stats_gls] = regress(cti_gls, X_gls);

fprintf('\nCTI ~ SST (GLS AR1): slope=%.4f, R2=%.3f, p=%.4f\n', b_gls(2), stats_gls(1), stats_gls(3));

% Tracking ratio
stats_sst      = calc_trend_AC(double(time_sst), sst_ann);
tracking_ratio = stats_cti.slope_Sen / stats_sst.slope_Sen;
fprintf('\nTracking ratio (CTI/SST): %.3f\n', tracking_ratio);

% Cross-correlation
max_lag = 10;
[xc, lags] = xcorr(zscore(cti_aligned), zscore(sst_aligned), max_lag, 'coeff');
[~, peak_idx] = max(xc);
best_lag = lags(peak_idx);
fprintf('CTI-SST peak correlation at lag = %d years\n', best_lag);

%% Biomass + CTI time series with trend lines
load '/Users/Clare/Documents/MATLAB/CPR/data/monthly_grids/monthly_cpr_NAtlantic_copepods 2022_gridded.mat'

C = z_copepod_biomass_output.cpr_annual_obana_trend;
C(cpr_not_enough) = NaN;
mask = ~isnan(C);
Cplot = C;
Cplot(~mask) = 0;

figure;

% Biomass map
subplot(2,2,1)
h = worldmap([min_lat max_lat], [min_lon max_lon]);
hs = pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, Cplot);
set(hs, 'EdgeColor','none', 'FaceAlpha','texturemap', ...
    'AlphaData',double(mask), 'AlphaDataMapping','none');
hold on
clim_val = max(abs(Cplot(mask)));
caxis(gca, [-clim_val clim_val]);
colormap(gca, b2r(-clim_val, clim_val));
cbh = colorbar('Location','westoutside');
p = findobj(h,'type','patch');
set(p,'FaceColor',[0.8 0.8 0.8]);
geoshow('landareas.shp', 'FaceColor', 'k'); gridm('off');
title('Copepod Biomass');
set(get(cbh,'ylabel'), 'string', 'log_{10}(\mugC per 3m^3+1).yr^-^1', 'fontsize',14);
cbh.Ruler.Exponent = 0;
drawnow;
cbh.Position(1) = cbh.Position(1) - 0.05;

% Biomass time series with AC-corrected trend
ln = (lon <= max_lon) & (lon >= min_lon);
lt = (lat <= max_lat) & (lat >= min_lat);
monthly_cpr_NAtl_int = z_copepod_biomass_output.monthly_cpr(:, lt, ln);
var2plot = squeeze(nanmean(monthly_cpr_NAtl_int, 2));
var2plot = squeeze(nanmean(var2plot, 2));
var2plot = reshape(var2plot, 12, []);
var2plot = squeeze(nanmean(var2plot, 1))';
var2plot_raw = var2plot;
var2plot = zscore(var2plot_raw);
time = reshape(time, 12, []);
time = time(1,:)';

subplot(2,2,2)
H2 = bar(time, var2plot);
datetick('x','yyyy');
ylim([-3 3])
set(H2,'EdgeColor', 'none');
drawnow; pause(0.05);
H2.Face.ColorType = 'truecoloralpha';
H2.Face.ColorData(4) = 255 * 0.3;

msk = ~isnan(var2plot(:));
coef = polyfit(time(~msk), var2plot(~msk), 1);
h = refline(coef(1), coef(2));
h.LineWidth = 2;

% AC-corrected significance
s = regstats(time(~msk), var2plot_raw(~msk), 'linear', 'all');
[my,~,~,smy] = lsqfity(time, var2plot_raw);
P_NEW = mod_glover_AC(time, var2plot_raw);
P_CRIT = tinv(0.975, P_NEW);
P_NEW_SE_SLOPE = smy / sqrt(P_NEW);
P_NEW_TB_SLOPE = abs(my) / P_NEW_SE_SLOPE;
INDEX_SIG = P_NEW_TB_SLOPE < P_CRIT;

if INDEX_SIG
    disp('Copepod Biomass Trend is not significant after accounting for autocorrelation');
else
    disp('Copepod Biomass Trend is significant after accounting for autocorrelation');
end
set(gca,'fontsize',12);
ylabel('zscore')
title('Copepod Biomass','fontsize',14);

% CTI map
C = cti_out_ALL_NAtl_copepods.cpr_annual_trend;
C(cpr_not_enough) = NaN;
mask = ~isnan(C);
Cplot = C;
Cplot(~mask) = 0;

subplot(2,2,3)
h = worldmap([min_lat max_lat], [min_lon max_lon]);
hs = pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, Cplot);
set(hs, 'EdgeColor','none', 'FaceAlpha','texturemap', ...
    'AlphaData',double(mask), 'AlphaDataMapping','none');
hold on
clim_val = max(abs(Cplot(mask)));
caxis(gca, [-clim_val clim_val]);
colormap(gca, b2r(-clim_val, clim_val));
cbh = colorbar('Location','westoutside');
p = findobj(h,'type','patch');
set(p,'FaceColor',[0.8 0.8 0.8]);
geoshow('landareas.shp', 'FaceColor', 'k'); gridm('off');
title('Copepod CTI');
set(get(cbh,'ylabel'), 'string', '^{\circ}C.yr^{-1}', 'fontsize',14);
cbh.Ruler.Exponent = 0;
drawnow;
cbh.Position(1) = cbh.Position(1) - 0.05;

% CTI time series
monthly_cpr_NAtl_int = cti_out_ALL_NAtl_copepods.monthly_cti.monthly_cpr(:, lt, ln);
var2plot = squeeze(nanmean(monthly_cpr_NAtl_int, 2));
var2plot = squeeze(nanmean(var2plot, 2));
var2plot = reshape(var2plot, 12, []);
var2plot = squeeze(nanmean(var2plot, 1))';
var2plot_raw = var2plot;
var2plot = zscore(var2plot_raw);

subplot(2,2,4)
H2 = bar(time, var2plot);
datetick('x','yyyy');
ylim([-3 3])
set(H2,'EdgeColor', 'none');
drawnow; pause(0.05);
H2.Face.ColorType = 'truecoloralpha';
H2.Face.ColorData(4) = 255 * 0.3;

msk = ~isnan(var2plot(:));
coef = polyfit(time(~msk), var2plot(~msk), 1);
h = refline(coef(1), coef(2));
h.LineWidth = 2;

[my,~,~,smy] = lsqfity(time, var2plot_raw);
P_NEW = mod_glover_AC(time, var2plot_raw);
P_CRIT = tinv(0.975, P_NEW);
P_NEW_SE_SLOPE = smy / sqrt(P_NEW);
P_NEW_TB_SLOPE = abs(my) / P_NEW_SE_SLOPE;
INDEX_SIG = P_NEW_TB_SLOPE < P_CRIT;

if INDEX_SIG
    disp('Copepod CTI Trend is not significant after accounting for autocorrelation');
else
    disp('Copepod CTI Trend is significant after accounting for autocorrelation');
end
set(gca,'fontsize',12);
ylabel('zscore')
title('Copepod CTI','fontsize',14);

%% Summary table export
Variable = {'Biomass (log10 ugC)'; 'CTI (degC)'; 'SST (degC)'; ...
    assemblage_names{1}; assemblage_names{2}; ...
    assemblage_names{3}; assemblage_names{4}; ...
    'Tracking ratio (CTI/SST)'};

all_stats_arr = [stats_biomass; stats_cti; stats_sst; ...
    stats_assemblage(1).stats; stats_assemblage(2).stats; ...
    stats_assemblage(3).stats; stats_assemblage(4).stats];

Sen_slope      = [arrayfun(@(s) s.slope_Sen, all_stats_arr); tracking_ratio];
Sen_per_decade = [Sen_slope(1:7) * 10; NaN];
Pct_change     = [arrayfun(@(s) s.pct_change, all_stats_arr); NaN];
N_eff          = [arrayfun(@(s) s.N_eff, all_stats_arr); NaN];
AC_p           = [arrayfun(@(s) s.p_AC, all_stats_arr); NaN];
MK_p           = [arrayfun(@(s) s.MK_p, all_stats_arr); NaN];
Significant    = [AC_p(1:7) < 0.05; NaN];

T = table(Variable, Sen_slope, Sen_per_decade, Pct_change, N_eff, AC_p, MK_p, Significant);
outpath = '/Users/Clare/Documents/MATLAB/CPR/ALLNAtlantic/biomass/figures/summary_stats.csv';
% writetable(T, outpath);

%% -----------------------------------------------------------------------
%  CLIMATE INDICES
%% -----------------------------------------------------------------------

% NAO
load('/Users/Clare/Documents/MATLAB/DATA/NAO/NAO_SPG_indices.mat');

% SPG Index (PC2, ORAS5-derived)
T_spg         = readtable('/Users/Clare/Documents/MATLAB/DATA/SPG/SPG_index_1958_2025_v2.xlsx');
time_spg_long = double(T_spg.Year(:));
spg_long      = double(T_spg.PC2(:));
spg_long      = spg_long ./ max(abs(spg_long));

valid_spg     = ~isnan(time_spg_long) & ~isnan(spg_long);
time_spg_long = time_spg_long(valid_spg);
spg_long      = spg_long(valid_spg);

% SSS Index
T_sss    = readtable('/Users/Clare/Documents/MATLAB/DATA/SSS/SSS_index_1958_2025_v3.xlsx');
time_sss = double(T_sss.Year(:));
sss_long = double(T_sss.MeanSSS(:));
sss_long = zscore(sss_long);

valid_sss = ~isnan(time_sss) & ~isnan(sss_long);
time_sss  = time_sss(valid_sss);
sss_long  = sss_long(valid_sss);

%% 5-panel driver/response overview
figure('Position', [100 100 900 1100]);

subplot(5,1,1)
bar(double(time_nao), nao_ann, 'FaceColor',[0.5 0.5 0.5], 'EdgeColor','none');
ylabel('NAO index'); title('Winter NAO (Hurrell)');
set(gca,'FontSize',12); yline(0,'k--'); xlim([1958 2025]);

subplot(5,1,2)
bar(time_spg_long, spg_long, 'FaceColor',[0.2 0.4 0.8], 'EdgeColor','none');
ylabel('SPG PC2'); title('Subpolar Gyre Index (PC2, 1958-2025)');
set(gca,'FontSize',12); yline(0,'k--'); xlim([1958 2025]);

subplot(5,1,3)
bar(time_sss, sss_long, 'FaceColor',[0.2 0.7 0.6], 'EdgeColor','none');
ylabel('z-Score'); title('Sea Surface Salinity');
set(gca,'FontSize',12); yline(0,'k--'); xlim([1958 2025]);

subplot(5,1,4)
bar(time_bio, zscore(var_bio_ann), 'FaceColor',[0.8 0.3 0.3], 'EdgeColor','none');
ylabel('z-score'); title('Copepod Biomass');
set(gca,'FontSize',12); yline(0,'k--'); xlim([1958 2025]);

subplot(5,1,5)
bar(time_cti, zscore(var_cti_ann), 'FaceColor',[0.9 0.6 0.1], 'EdgeColor','none');
ylabel('z-score'); title('CTI');
set(gca,'FontSize',12); yline(0,'k--'); xlim([1958 2025]);
xlabel('Year');

sgtitle('North Atlantic drivers and copepod responses 1958-2025', ...
    'FontSize',13, 'FontWeight','bold');

%% CTI~SST residual analysis
[~, ia, ib] = intersect(time_cti, time_sst);
cti_aligned  = var_cti_ann(ia);   cti_aligned = cti_aligned(:);
sst_aligned  = sst_ann(ib);       sst_aligned = sst_aligned(:);
t_plot       = double(time_cti(ia));  t_plot = t_plot(:);

cti_predicted = b_gls(1) + b_gls(2) .* sst_aligned;
cti_residuals = cti_aligned - cti_predicted;
resid_sd      = std(cti_residuals, 'omitnan');
resid_mean    = mean(cti_residuals, 'omitnan');

large_pos = t_plot(cti_residuals > resid_mean + resid_sd);
large_neg = t_plot(cti_residuals < resid_mean - resid_sd);
fprintf('Years CTI > SST prediction by >1 SD: '); fprintf('%d ', large_pos); fprintf('\n');
fprintf('Years CTI < SST prediction by >1 SD: '); fprintf('%d ', large_neg); fprintf('\n');

%% CTI-SST deviation figure
figure('Position', [100 100 1000 500]);

subplot(2,1,1)
yyaxis left
plot(t_plot, zscore(cti_aligned), '-o', 'Color',[0.9 0.5 0.1], 'LineWidth',1.5, 'MarkerSize',3);
hold on
plot(t_plot, movmean(zscore(cti_aligned), 5, 'omitnan'), '-', 'Color',[0.9 0.5 0.1], 'LineWidth',2.5);
ylabel('CTI (z-score)', 'Color',[0.9 0.5 0.1]);
ax = gca; ax.YColor = [0.9 0.5 0.1];
yyaxis right
plot(t_plot, zscore(sst_aligned), '-o', 'Color',[0.2 0.4 0.8], 'LineWidth',1.5, 'MarkerSize',3);
plot(t_plot, movmean(zscore(sst_aligned), 5, 'omitnan'), '-', 'Color',[0.2 0.4 0.8], 'LineWidth',2.5);
ylabel('SST (z-score)', 'Color',[0.2 0.4 0.8]);
ax.YColor = [0.2 0.4 0.8];
xregion(1986, 1995, 'FaceColor',[0.8 0.8 0.8], 'FaceAlpha',0.3, 'EdgeColor','none');
xregion(2018, 2022, 'FaceColor',[1.0 0.8 0.6], 'FaceAlpha',0.3, 'EdgeColor','none');
xlabel('Year'); title('A  CTI and SST - standardised');
legend({'CTI','CTI 5-yr mean','SST','SST 5-yr mean','NAO/SPG shift ~1990','Heatwave ~2020'}, ...
    'Location','northwest','FontSize',9);
set(gca,'FontSize',11); xlim([1955 2025]);

subplot(2,1,2)
bar(t_plot, cti_residuals, 'FaceColor',[0.5 0.5 0.5], 'EdgeColor','none');
hold on
yline(resid_mean + resid_sd, '--r', 'LineWidth',1.5);
yline(resid_mean - resid_sd, '--b', 'LineWidth',1.5);
yline(0, 'k', 'LineWidth',0.5);
xregion(1986, 1995, 'FaceColor',[0.8 0.8 0.8], 'FaceAlpha',0.3, 'EdgeColor','none');
xregion(2018, 2022, 'FaceColor',[1.0 0.8 0.6], 'FaceAlpha',0.3, 'EdgeColor','none');
for yy = 1:numel(t_plot)
    if abs(cti_residuals(yy)) > abs(resid_mean) + 1.5*resid_sd
        text(t_plot(yy), cti_residuals(yy) + sign(cti_residuals(yy))*0.02, ...
            num2str(t_plot(yy)), 'FontSize',8, 'HorizontalAlignment','center');
    end
end
xlabel('Year'); ylabel('CTI residual (deg C)');
title('B  CTI residuals from SST regression');
legend('Location','northwest','FontSize',9);
set(gca,'FontSize',11); xlim([1955 2025]);
sgtitle('CTI-SST deviations: periods where biological change diverges from physical forcing', ...
    'FontSize',13, 'FontWeight','bold');

%% Residual correlations with climate drivers (AC-corrected)
% NAO
[~, ia_nao, ib_nao] = intersect(double(time_nao), t_plot);
nao_res     = nao_ann(ia_nao);
resid_nao   = cti_residuals(ib_nao);
r_nao_resid = corr(nao_res(:), resid_nao(:), 'rows','complete');

% SPG PC2
[~, ia_spg, ib_spg] = intersect(time_spg_long, t_plot);
spg_res     = spg_long(ia_spg);
resid_spg   = cti_residuals(ib_spg);
r_spg_resid = corr(spg_res(:), resid_spg(:), 'rows','complete');

% SSS
[~, ia_sss_r, ib_sss_r] = intersect(time_sss, t_plot);
sss_res       = sss_long(ia_sss_r);
resid_sss_r   = cti_residuals(ib_sss_r);
r_sss_resid   = corr(sss_res(:), resid_sss_r(:), 'rows','complete');

fprintf('\n=== RESIDUAL CORRELATIONS ===\n');
fprintf('NAO vs CTI residuals: r = %.3f\n', r_nao_resid);
fprintf('SPG vs CTI residuals: r = %.3f\n', r_spg_resid);
fprintf('SSS vs CTI residuals: r = %.3f\n', r_sss_resid);

% AC-corrected p-values for SPG and SSS
t_spg_com = time_spg_long(ia_spg);
[my_spg, ~, ~, smy_spg] = lsqfity(spg_res(:), resid_spg(:));
P_spg = mod_glover_AC(t_spg_com(:), resid_spg(:));
t_stat_spg = abs(my_spg) / (smy_spg / sqrt(P_spg));
p_spg_AC = 2 * (1 - tcdf(t_stat_spg, P_spg));

t_sss_com = time_sss(ia_sss_r);
[my_sss, ~, ~, smy_sss] = lsqfity(sss_res(:), resid_sss_r(:));
P_sss = mod_glover_AC(t_sss_com(:), resid_sss_r(:));
t_stat_sss = abs(my_sss) / (smy_sss / sqrt(P_sss));
p_sss_AC = 2 * (1 - tcdf(t_stat_sss, P_sss));

fprintf('SPG AC-corrected p=%.4f (Neff=%.1f)\n', p_spg_AC, P_spg);
fprintf('SSS AC-corrected p=%.4f (Neff=%.1f)\n', p_sss_AC, P_sss);

%% SPG + CTI residuals comparison
load /Users/Clare/Documents/MATLAB/DATA/SPG/SPG_Nicolas_1958_2020.mat

figure('Position',[100 100 1000 600]);

subplot(2,1,1)
bar(time_spg_long, spg_long, 'FaceColor',[0.2 0.4 0.8],'EdgeColor','none');
hold on; yline(0,'k--');
yyaxis right
plot(time_hatun, hatun_index, 'k-', 'LineWidth',2);
xregion(1993, 2016, 'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.3,'EdgeColor','none');
yyaxis left
ylabel('SPG PC2 (v2)','FontSize',12);
yyaxis right
ylabel('Hatun index','FontSize',12);
title('A  Subpolar Gyre Index PC2  1958-2025','FontSize',13,'FontWeight','bold');
xlim([1955 2026]); set(gca,'FontSize',11);

subplot(2,1,2)
yyaxis left
bar(t_plot, cti_residuals, 'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.6);
ylabel('CTI residual (deg C)','FontSize',12); hold on; yline(0,'k');
yyaxis right
plot(time_spg_long, spg_long, 'b-', 'LineWidth',2);
ylabel('SPG PC2','FontSize',12,'Color','b');
ax = gca; ax.YColor = 'b';
xregion(1993, 2016, 'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.3,'EdgeColor','none');
xlabel('Year','FontSize',12);
title('B  CTI residuals vs SPG PC2','FontSize',13,'FontWeight','bold');
legend({'CTI residual','SPG PC2','Satellite era'},'Location','northwest','FontSize',10);
xlim([1955 2026]); set(gca,'FontSize',11);

sgtitle('SPG influence on copepod community thermalization', ...
    'FontSize',13, 'FontWeight','bold');

%% Tropicalization / deborealization decomposition (Chust et al. 2024)
assemblage_STI         = [mean_warm_STI, mean_pseudo_STI, mean_cold_STI, mean_subarctic_STI];
assemblage_names_short = {'Warm-temperate','Pseudo-oceanic','Cold-temperate','Subarctic'};
CTI_mean               = mean(var_cti_ann, 'omitnan');
thermal_bias           = assemblage_STI - CTI_mean;
abundance_change       = [stats_assemblage(1).stats.slope_Sen, ...
                          stats_assemblage(2).stats.slope_Sen, ...
                          stats_assemblage(3).stats.slope_Sen, ...
                          stats_assemblage(4).stats.slope_Sen];

process = cell(1,4);
for aa = 1:4
    if     thermal_bias(aa) > 0 && abundance_change(aa) > 0; process{aa} = 'Tropicalization';
    elseif thermal_bias(aa) < 0 && abundance_change(aa) < 0; process{aa} = 'Deborealization';
    elseif thermal_bias(aa) < 0 && abundance_change(aa) > 0; process{aa} = 'Borealization';
    elseif thermal_bias(aa) > 0 && abundance_change(aa) < 0; process{aa} = 'Detropicalization';
    end
    fprintf('%s: STI=%.2f, bias=%.2f, slope=%.4f -> %s\n', ...
        assemblage_names_short{aa}, assemblage_STI(aa), thermal_bias(aa), ...
        abundance_change(aa), process{aa});
end

intensity        = thermal_bias .* abundance_change;
trop_intensity   = sum(intensity(strcmp(process,'Tropicalization')));
deb_intensity    = sum(intensity(strcmp(process,'Deborealization')));
total_intensity  = sum(abs(intensity));

fprintf('\nTropicalization:  %.4f (%.1f%%)\n', trop_intensity, 100*trop_intensity/total_intensity);
fprintf('Deborealization:  %.4f (%.1f%%)\n', abs(deb_intensity), 100*abs(deb_intensity)/total_intensity);

% Pie chart
process_values = [abs(trop_intensity), abs(deb_intensity)];
process_pct    = 100 * process_values / sum(process_values);

figure;
p = pie(process_pct);
p(1).FaceColor = [0.8 0.2 0.2];  p(3).FaceColor = [1.0 0.6 0.6];
p(2).String = sprintf('Tropicalization\n%.1f%%', process_pct(1));
p(4).String = sprintf('Deborealization\n%.1f%%', process_pct(2));
p(2).FontSize = 13; p(4).FontSize = 13;
title({'Underlying processes driving CTI change','North Atlantic copepods 1958-2022'}, ...
    'FontSize',13, 'FontWeight','bold');
annotation('textbox', [0.15 0.02 0.7 0.08], ...
    'String','cf. Chust et al. (2024): Tropicalization 54%, Deborealization 18% across European seas', ...
    'FontSize',9, 'EdgeColor',[0.7 0.7 0.7], 'HorizontalAlignment','center', 'BackgroundColor','w');

%% -----------------------------------------------------------------------
%  FIGURE 4: Biological response to ocean warming
%  A: CTI+SST time series  B: CTI~SST scatter  C: Stacked drivers  D: Tracking ratio
%% -----------------------------------------------------------------------

load('/Users/Clare/Documents/MATLAB/DATA/AMO/monthly_AMO.mat');
years_amo = monthly_AMO_unsmoothed(:,1);
amo_ann   = nanmean(monthly_AMO_unsmoothed(:,2:13), 2);
time_amo  = years_amo;

figure('Position', [100 100 1400 900]);

% Panel A: CTI and SST time series
ax_a = subplot(2,2,1);
[~, ia, ib] = intersect(time_cti, time_sst);
cti_plot  = var_cti_ann(ia);
sst_plot  = sst_ann(ib);
t_plot_ts = time_cti(ia);
ylim_cti  = [min(cti_plot)-0.1, max(cti_plot)+0.1];
ylim_sst  = [min(sst_plot)-0.1, max(sst_plot)+0.1];



yyaxis left
h1 = plot(t_plot_ts, cti_plot, '-o', 'Color',[0.9 0.5 0.1], ...
    'LineWidth',1.5, 'MarkerSize',3, 'MarkerFaceColor',[0.9 0.5 0.1]);
h2 = plot(t_plot_ts, movmean(cti_plot, 5, 'omitnan'), '-', ...
    'Color',[0.9 0.5 0.1], 'LineWidth',3);
ylabel('CTI (\circC)', 'FontSize',11, 'Color',[0.9 0.5 0.1]);
ylim(ylim_cti); ax = gca; ax.YColor = [0.9 0.5 0.1];

yyaxis right
h3 = plot(t_plot_ts, sst_plot, '-o', 'Color',[0.2 0.4 0.8], ...
    'LineWidth',1.5, 'MarkerSize',3, 'MarkerFaceColor',[0.2 0.4 0.8]);
h4 = plot(t_plot_ts, movmean(sst_plot, 5, 'omitnan'), '-', ...
    'Color',[0.2 0.4 0.8], 'LineWidth',3);
ylabel('SST (\circC)', 'FontSize',11, 'Color',[0.2 0.4 0.8]);
ylim(ylim_sst); ax.YColor = [0.2 0.4 0.8];

xlabel('Year','FontSize',12); set(gca,'FontSize',10);
legend([h1,h2,h3,h4], {'CTI','CTI (5-yr mean)','SST','SST (5-yr mean)'}, ...
    'Location','northwest','FontSize',8);
xlim([1958 2023]); box on

% Panel B: CTI ~ SST scatter + GLS fit
ax_b = subplot(2,2,2);
scatter(sst_plot, cti_plot, 70, t_plot_ts, 'filled', ...
    'MarkerEdgeAlpha',0.6, 'MarkerFaceAlpha',0.7, 'MarkerEdgeColor',[0.4 0.4 0.4]);
colormap(gca, cool(max(t_plot_ts)-min(t_plot_ts)+1));
cb = colorbar; cb.Label.String = 'Year'; cb.FontSize = 10;
hold on
xl = [min(sst_plot)-0.1, max(sst_plot)+0.1];
h_gls = plot(xl, b_gls(1) + b_gls(2).*xl, 'k-', 'LineWidth',2.5, 'DisplayName','GLS fit (AR1)');
str = sprintf('Slope = %.3f deg C CTI / deg C SST\nR2 = %.2f\nTracking ratio = %.3f', ...
    b_gls(2), stats_gls(1), tracking_ratio);
text(min(sst_plot)+0.05, max(cti_plot)-0.08, str, 'FontSize',9, ...
    'BackgroundColor','w', 'EdgeColor',[0.7 0.7 0.7], 'Margin',4);
xlabel('SST (\circC)','FontSize',11); ylabel('CTI (\circC)','FontSize',11);
legend(h_gls, 'Location','southeast','FontSize',9);
set(gca,'FontSize',10); box on

% Panel C: 5 stacked sub-panels (AMO, SPG, SSS, CTI, Biomass)
ax_c_main = subplot(2,2,3);
pos_c = get(ax_c_main, 'Position');
delete(ax_c_main);

n_stacks    = 5;
spacing     = 0.007;
height_each = (pos_c(4) - (n_stacks-1)*spacing) / n_stacks;

% C1: AMO
k = 1;
ax_c1 = axes('Position', [pos_c(1), pos_c(2)+(n_stacks-k)*(height_each+spacing), pos_c(3), height_each]);
[~, ia_amo_c, ~] = intersect(time_amo, (1960:2022)');
amo_c = amo_ann(ia_amo_c); t_amo_c = time_amo(ia_amo_c);
amo_ylim = [min(amo_c)-0.05, max(amo_c)+0.05];

bar(t_amo_c, amo_c, 'FaceColor',[0.7 0.3 0.3],'EdgeColor','none','FaceAlpha',0.8);
yline(0,'k--','LineWidth',0.8);
ylabel('AMO','FontSize',9,'FontWeight','bold');
set(ax_c1,'FontSize',9,'XTickLabel',[]); xlim([1958 2023]); ylim(amo_ylim); box on

% C2: SPG PC2
k = 2;
ax_c2 = axes('Position', [pos_c(1), pos_c(2)+(n_stacks-k)*(height_each+spacing), pos_c(3), height_each]);
[~, ia_spg_c, ~] = intersect(time_spg_long, (1958:2022)');
spg_c = spg_long(ia_spg_c); t_spg_c = time_spg_long(ia_spg_c);
spg_ylim = [min(spg_c)-0.5, max(spg_c)+0.5];

bar(t_spg_c, spg_c, 'FaceColor',[0.2 0.4 0.8],'EdgeColor','none','FaceAlpha',0.8);
yline(0,'k--','LineWidth',0.8);
ylabel('SPG PC2','FontSize',9,'FontWeight','bold');
set(ax_c2,'FontSize',9,'XTickLabel',[]); xlim([1958 2023]); ylim(spg_ylim); box on

% C3: SSS
k = 3;
ax_c3 = axes('Position', [pos_c(1), pos_c(2)+(n_stacks-k)*(height_each+spacing), pos_c(3), height_each]);
[~, ia_sss_c, ~] = intersect(time_sss, (1958:2022)');
sss_c = sss_long(ia_sss_c); t_sss_c = time_sss(ia_sss_c);
sss_ylim = [min(sss_c)-0.5, max(sss_c)+0.5];

bar(t_sss_c, sss_c, 'FaceColor',[0.2 0.7 0.6],'EdgeColor','none','FaceAlpha',0.8);
yline(0,'k--','LineWidth',0.8);
ylabel('SSS (z)','FontSize',9,'FontWeight','bold');
set(ax_c3,'FontSize',9,'XTickLabel',[]); xlim([1958 2023]); ylim(sss_ylim); box on

% C4: CTI z-score
k = 4;
ax_c4 = axes('Position', [pos_c(1), pos_c(2)+(n_stacks-k)*(height_each+spacing), pos_c(3), height_each]);
cti_z = zscore(var_cti_ann);
cti_ylim = [min(cti_z)-0.1, max(cti_z)+0.1];

bar(double(time_cti), cti_z, 'FaceColor',[0.9 0.5 0.1],'EdgeColor','none','FaceAlpha',0.8);
yline(0,'k--','LineWidth',0.8);
ylabel('CTI (z)','FontSize',9,'FontWeight','bold');
set(ax_c4,'FontSize',9,'XTickLabel',[]); xlim([1958 2023]); ylim(cti_ylim); box on

% C5: Biomass z-score
k = 5;
ax_c5 = axes('Position', [pos_c(1), pos_c(2)+(n_stacks-k)*(height_each+spacing), pos_c(3), height_each]);
bio_z = zscore(var_bio_ann);
bio_ylim = [min(bio_z)-0.1, max(bio_z)+0.1];

bar(double(time_bio), bio_z, 'FaceColor',[0.8 0.3 0.3],'EdgeColor','none','FaceAlpha',0.8);
yline(0,'k--','LineWidth',0.8);
ylabel('Biomass (z)','FontSize',9,'FontWeight','bold');
xlabel('Year','FontSize',9);
set(ax_c5,'FontSize',9); xlim([1958 2023]); ylim(bio_ylim); box on

% Panel D: Tracking ratio context
ax_d = subplot(2,2,4);
lit_labels = {'Birds (Devictor 2012)', 'Plants (Devictor 2012)', ...
              sprintf('N. Atl. Copepods\n(this study)')};
lit_values = [0.57, 0.24, tracking_ratio];
lit_cols   = {[0.6 0.8 0.4], [0.3 0.7 0.3], [0.9 0.5 0.1]};

hold on
for bb = 1:3
    barh(bb, lit_values(bb), 0.5, 'FaceColor',lit_cols{bb}, ...
        'EdgeColor',[0.3 0.3 0.3], 'LineWidth',1.2);
end
xline(1.0, '--k', 'LineWidth',2);
for bb = 1:3
    text(lit_values(bb)+0.03, bb, sprintf('%.2f', lit_values(bb)), ...
        'FontSize',10, 'VerticalAlignment','middle', 'FontWeight','bold');
end
yticks(1:3); yticklabels(lit_labels);
xlabel('Tracking ratio (CTI / SST trend)', 'FontSize',11);
xlim([0 1.1]); set(gca,'FontSize',10); box on

fprintf('Figure 4 saved\n');

%% -----------------------------------------------------------------------
%  FIGURE: Biomass & CTI maps + dual-axis time series with SST
%% -----------------------------------------------------------------------

col_biomass = [0.6 0.2 0.5];
col_cti     = [0.9 0.5 0.1];
col_sst     = [0.8 0.1 0.1];

% Biomass map
C = z_copepod_biomass_output.cpr_annual_obana_trend;
C(cpr_not_enough) = NaN;
mask_map = ~isnan(C);
Cplot = C;
Cplot(~mask_map) = 0;

figure('Position', [100 100 1200 800]);

subplot(2,2,1)
h = worldmap([min_lat max_lat], [min_lon max_lon]);
hs = pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, Cplot);
set(hs, 'EdgeColor','none', 'FaceAlpha','texturemap', ...
    'AlphaData',double(mask_map), 'AlphaDataMapping','none');
hold on
clim_val = max(abs(Cplot(mask_map)));
caxis(gca, [-clim_val clim_val]);
colormap(gca, b2r(-clim_val, clim_val));
cbh = colorbar('Location','westoutside');
p = findobj(h,'type','patch');
set(p,'FaceColor',[0.8 0.8 0.8]);
geoshow('landareas.shp', 'FaceColor','k'); gridm('off');
title('Copepod Biomass');
set(get(cbh,'ylabel'), 'string', 'log_{10}(\mugC per 3m^3+1).yr^{-1}', 'fontsize',12);
cbh.Ruler.Exponent = 0;
drawnow;
cbh.Position(1) = cbh.Position(1) - 0.05;

% Biomass + SST dual-axis
ln = (lon <= max_lon) & (lon >= min_lon);
lt = (lat <= max_lat) & (lat >= min_lat);
monthly_bio = z_copepod_biomass_output.monthly_cpr(:, lt, ln);
var_bio     = squeeze(nanmean(nanmean(monthly_bio, 2), 3));
var_bio_ann = squeeze(nanmean(reshape(var_bio, 12, []), 1))';

[~, ia_bio, ib_bio] = intersect(uniqY, time_sst);
bio_plot   = var_bio_ann(ia_bio);
sst_plot   = sst_ann(ib_bio);
t_plot_bio = double(uniqY(ia_bio));

subplot(2,2,2)
hold on
yl_bio = [min(bio_plot)-0.05, max(bio_plot)+0.05];
yyaxis left
h1 = plot(t_plot_bio, bio_plot, '-o', 'Color',col_biomass, 'LineWidth',1.5, 'MarkerSize',3, 'MarkerFaceColor',col_biomass);
h2 = plot(t_plot_bio, movmean(bio_plot, 5, 'omitnan'), '-', 'Color',col_biomass, 'LineWidth',3);
ylabel('Biomass (log_{10} \mugC m^{-3})', 'FontSize',11, 'Color',col_biomass);
ylim(yl_bio);
ax = gca; ax.YColor = col_biomass;
yyaxis right
h3 = plot(t_plot_bio, sst_plot, '-o', 'Color',col_sst, 'LineWidth',1.5, 'MarkerSize',3, 'MarkerFaceColor',col_sst);
h4 = plot(t_plot_bio, movmean(sst_plot, 5, 'omitnan'), '-', 'Color',col_sst, 'LineWidth',3);
ylabel('SST (\circC)', 'FontSize',11, 'Color',col_sst);
ax.YColor = col_sst;
xlabel('Year','FontSize',12);
title('Copepod Biomass and SST', 'FontSize',14);
xlim([1958 2023]); box on

% CTI map
C = cti_out_ALL_NAtl_copepods.cpr_annual_trend;
C(cpr_not_enough) = NaN;
mask_map = ~isnan(C);
Cplot = C;
Cplot(~mask_map) = 0;

subplot(2,2,3)
h = worldmap([min_lat max_lat], [min_lon max_lon]);
hs = pcolorm(min_lat:1:max_lat, min_lon:1:max_lon, Cplot);
set(hs, 'EdgeColor','none', 'FaceAlpha','texturemap', ...
    'AlphaData',double(mask_map), 'AlphaDataMapping','none');
hold on
clim_val = max(abs(Cplot(mask_map)));
caxis(gca, [-clim_val clim_val]);
colormap(gca, b2r(-clim_val, clim_val));
cbh = colorbar('Location','westoutside');
p = findobj(h,'type','patch');
set(p,'FaceColor',[0.8 0.8 0.8]);
geoshow('landareas.shp', 'FaceColor','k'); gridm('off');
title('Copepod CTI');
set(get(cbh,'ylabel'), 'string', '^{\circ}C.yr^{-1}', 'fontsize',12);
cbh.Ruler.Exponent = 0;
drawnow;
cbh.Position(1) = cbh.Position(1) - 0.05;

% CTI + SST dual-axis
monthly_cti_reg = cti_out_ALL_NAtl_copepods.monthly_cti.monthly_cpr(:, lt, ln);
var_cti         = squeeze(nanmean(nanmean(monthly_cti_reg, 2), 3));
var_cti_ann     = squeeze(nanmean(reshape(var_cti, 12, []), 1))';

[~, ia_cti, ib_cti] = intersect(uniqY, time_sst);
cti_plot     = var_cti_ann(ia_cti);
sst_plot_cti = sst_ann(ib_cti);
t_plot_cti   = double(uniqY(ia_cti));

subplot(2,2,4)
hold on
yl_cti = [min(cti_plot)-0.1, max(cti_plot)+0.1];
yyaxis left
h1 = plot(t_plot_cti, cti_plot, '-o', 'Color',col_cti, 'LineWidth',1.5, 'MarkerSize',3, 'MarkerFaceColor',col_cti);
h2 = plot(t_plot_cti, movmean(cti_plot, 5, 'omitnan'), '-', 'Color',col_cti, 'LineWidth',3);
ylabel('CTI (\circC)', 'FontSize',11, 'Color',col_cti);
ylim(yl_cti);
ax = gca; ax.YColor = col_cti;
yyaxis right
h3 = plot(t_plot_cti, sst_plot_cti, '-o', 'Color',col_sst, 'LineWidth',1.5, 'MarkerSize',3, 'MarkerFaceColor',col_sst);
h4 = plot(t_plot_cti, movmean(sst_plot_cti, 5, 'omitnan'), '-', 'Color',col_sst, 'LineWidth',3);
ylabel('SST (\circC)', 'FontSize',11, 'Color',col_sst);
ax.YColor = col_sst;
xlabel('Year','FontSize',12);
title('Copepod CTI and SST', 'FontSize',14);
set(gca,'FontSize',11);
xlim([1958 2023]); box on

%% -----------------------------------------------------------------------
%  FIGURE: 8-panel stacked overview
%  CTI | Biomass | C. finmarchicus | C. helgolandicus | AMO | NAO | SPG | SSS
%% -----------------------------------------------------------------------

col_cti     = [0.9 0.5 0.1];
col_biomass = [0.6 0.2 0.5];
col_cfin    = [0.2 0.4 0.7];
col_chel    = [0.8 0.4 0.3];
col_amo     = [0.7 0.3 0.3];
col_nao     = [0.5 0.5 0.5];
col_spg     = [0.2 0.4 0.8];
col_sss     = [0.2 0.7 0.6];



% C. finmarchicus
load('/Users/Clare/Documents/MATLAB/CPR/data/monthly_grids/monthly_cpr_Calanus finmarchicus 2022_gridded.mat');
ln_sp = (lon <= max_lon) & (lon >= min_lon);
lt_sp = (lat <= max_lat) & (lat >= min_lat);
cfin_monthly = monthly_cpr(:, lt_sp, ln_sp);
cfin_ts  = squeeze(nanmean(nanmean(cfin_monthly, 2), 3));
cfin_ann = squeeze(nanmean(reshape(cfin_ts, 12, []), 1))';
cfin_z   = zscore(cfin_ann);
date_v_sp = datevec(time);
time_sp   = unique(date_v_sp(:,1));
time_cfin = double(time_sp(1:numel(cfin_ann)));

% C. helgolandicus
load('/Users/Clare/Documents/MATLAB/CPR/data/monthly_grids/monthly_cpr_Calanus helgolandicus 2022_gridded.mat');
chel_monthly = monthly_cpr(:, lt_sp, ln_sp);
chel_ts  = squeeze(nanmean(nanmean(chel_monthly, 2), 3));
chel_ann = squeeze(nanmean(reshape(chel_ts, 12, []), 1))';
chel_z   = zscore(chel_ann);
time_chel = time_cfin;

% Figure
n_panels   = 8;
fig_height = 160 * n_panels;
figure('Position', [100 50 1000 fig_height]);

left       = 0.10;
width      = 0.85;
gap        = 0.008;
top_margin = 0.04;
bot_margin = 0.06;
usable     = 1 - top_margin - bot_margin;
h_each     = (usable - (n_panels-1)*gap) / n_panels;

panel_data = { ...
    double(time_cti), zscore(var_cti_ann), col_cti,     'CTI (z)';
    double(time_bio), zscore(var_bio_ann), col_biomass,  'Biomass (z)';
    time_cfin,        cfin_z,              col_cfin,     '{\itC. fin.} (z)';
    time_chel,        chel_z,              col_chel,     '{\itC. hel.} (z)';
    double(time_amo), amo_ann,             col_amo,      'AMO';
    double(time_nao), nao_ann,             col_nao,      'NAO';
    time_spg_long,    spg_long,            col_spg,      'SPG';
    time_sss,         sss_long,            col_sss,      'SSS (z)'};

ax_arr = gobjects(n_panels, 1);
for k = 1:n_panels
    ax_arr(k) = axes('Position', [left, 1-top_margin-k*h_each-(k-1)*gap, width, h_each]);

    t_k = panel_data{k,1};
    v_k = panel_data{k,2};
    idx_k = (t_k >= xl(1)) & (t_k <= xl(2));

    yl = [min(v_k(idx_k))-0.3, max(v_k(idx_k))+0.3];
    bar(t_k(idx_k), v_k(idx_k), 'FaceColor',panel_data{k,3}, 'EdgeColor','none', 'FaceAlpha',0.8);
    hold on; yline(0,'k--','LineWidth',0.8);
    ylabel(panel_data{k,4}, 'FontSize',11, 'FontWeight','bold');
    set(gca, 'FontSize',10, 'TickLength',[0 0]);
    xlim(xl); ylim(yl); box on

    if k < n_panels
        set(gca, 'XTickLabel',[]);
    else
        xlabel('Year','FontSize',11);
    end
end
linkaxes(ax_arr, 'x');
fprintf('8-panel figure created\n');

%% Console summary
fprintf('\n========================================================================\n');
fprintf('FIGURE 4 SUMMARY STATISTICS\n');
fprintf('========================================================================\n');
fprintf('CTI~SST GLS slope: %.4f degC CTI per degC SST\n', b_gls(2));
fprintf('R2: %.3f,  p: %.4f\n', stats_gls(1), stats_gls(3));
fprintf('Tracking ratio: %.3f (~%.0f%% of SST warming rate)\n', tracking_ratio, 100*tracking_ratio);
fprintf('SPG vs CTI residuals: r=%.3f, AC p=%.4f\n', r_spg_resid, p_spg_AC);
fprintf('SSS vs CTI residuals: r=%.3f, AC p=%.4f\n', r_sss_resid, p_sss_AC);
fprintf('========================================================================\n\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function stats = calc_trend_AC(time_yr, y)
% Autocorrelation-corrected trend test with Mann-Kendall and Theil-Sen
    ok = ~isnan(y) & ~isnan(time_yr);
    t  = time_yr(ok);
    yy = y(ok);
    n  = numel(t);

    [my, by, ~, smy] = lsqfity(t, yy);
    stats.slope_OLS = my;
    stats.intercept = by;

    pairs = nchoosek(1:n, 2);
    slopes = (yy(pairs(:,2)) - yy(pairs(:,1))) ./ (t(pairs(:,2)) - t(pairs(:,1)));
    stats.slope_Sen = median(slopes);

    P_NEW    = mod_glover_AC(t, yy);
    P_CRIT   = tinv(0.975, P_NEW);
    SE_slope = smy / sqrt(P_NEW);
    t_stat   = abs(my) / SE_slope;
    stats.p_AC   = 2 * (1 - tcdf(t_stat, P_NEW));
    stats.sig_AC = t_stat >= P_CRIT;
    stats.N_eff  = P_NEW;
    stats.N_raw  = n;

    S = 0;
    for ii = 1:n-1
        S = S + sum(sign(yy(ii+1:end) - yy(ii)));
    end
    varS = (n*(n-1)*(2*n+5)) / 18;
    if     S > 0; z_mk = (S-1) / sqrt(varS);
    elseif S < 0; z_mk = (S+1) / sqrt(varS);
    else;         z_mk = 0;
    end
    stats.MK_S = S;
    stats.MK_z = z_mk;
    stats.MK_p = 2 * (1 - normcdf(abs(z_mk)));

    early = mean(yy(t <= t(1) + 10));
    late  = mean(yy(t >= t(end) - 10));
    stats.pct_change = 100 * (late - early) / abs(early);
end

