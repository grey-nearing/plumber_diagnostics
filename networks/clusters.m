%function [Mc,Sc,Mi,Si] = clusters(TD,MI);

Ai = zeros(Nx,Ny)/0;
Si = zeros(Nx,Ny,Ns)/0;
Mi = zeros(Nx,Ny,Nm)/0;
Sc = zeros(Nx,Ny)/0;
Mc = zeros(Nx,Ny)/0;

%% model groups
%clear Im;
%Im{1} = [1,2];
%Im{2} = [3];
%Im{3} = [4];
%Im{4} = [5,6];
%Im{5} = [7,8];
%Im{6} = [9];
%Im{7} = [10,11,12];
%Im{8} = [13,14];
%Ngroups = length(Im);

for x = 1:Nx
    for y = 1:Ny
        
% --- all stats ----------------------------------------------------------
        
        % all cluster center
        td = squeeze(TD(x,y,:,:));
        mi = squeeze(MI(y,:,:));
        
        % remove any bad site
        Is = find(any(~isnan([td,mi]')));
        J = find(all(isnan([td,mi]'))); td(J,:) = []; mi(J,:) = [];
        assert(isempty(J))
        
        % remove any bad model
        Im = find(all(~isnan([td;mi])));
        J = find(any(isnan([td;mi]))); td(:,J) = []; mi(:,J) = [];
        
        [Nsites,Nmodels] = size(td);
        [x,y,Nsites,Nmodels]
        
        % find cluster centers
        Xa = mean(td(:));
        Ya = mean(mi(:));
        
        % all cluster mean distance
        Ai(x,y) = mean( ((td(:)-Xa).^2 + (mi(:)-Ya).^2).^(1/2) );
        
% --- model stats --------------------------------------------------------
        
        % init storage
        Xm = zeros(Nmodels,1)/0;
        Ym = zeros(Nmodels,1)/0;
        
        % loop through models
        for m = 1:Nmodels
            
            %   % remove missing models from group
            %   Idex = Im{m};
            %   Ir = [];
            %   for ii = 1:length(Idex)
            %    for jj = 1:length(J)
            %     if Idex(ii) == J(jj); Ir = [Ir,ii]; end;
            %    end
            %   end
            %   Idex(Ir) = [];
            %   if isempty(Idex)
            %    Mi(x,y,m) = 0/0;
            %    Ym(m) = 0/0;
            %    Xm(m) = 0/0;
            %    continue
            %   end
            
            % pull model group
            tdg = td(:,m); tdg = tdg(:);
            mig = mi(:,m); mig = mig(:);
            
            % find cluster centers
            Xm(m) = mean(tdg);
            Ym(m) = mean(mig);
            
            % mean distance of each model from their center
            Mi(x,y,m) = mean( ((tdg-Xm(m)).^2 + (mig-Ym(m)).^2).^(1/2) );
        end
        
        % distance of cluster centers to mean
        Mca(x,y) = mean( ((Xm-Xa).^2 + (Ym-Ya).^2).^(1/2) );
        
        % mean distance between centers
        Mc(x,y) = 0;
        count = 0;
        for m = 1:Nmodels
            for mm = 1:m-1
                Mc(x,y) = Mc(x,y) + ((Xm(m)-Xm(mm))^2 + (Ym(m)-Ym(mm))^2)^(1/2);
                count = count + 1;
            end
        end
        Mc(x,y) = Mc(x,y)/count;
        
% --- site stats ---------------------------------------------------------
        
        % init storage
        Xs = zeros(Nsites,1)/0;
        Ys = zeros(Nsites,1)/0;
        
        % loop through sites
        for s = 1:Nsites
            
            % find cluster centers
            Xs(s) = mean(td(s,:));
            Ys(s) = mean(mi(s,:));
            
            % mean distance of each model from their center
            Si(x,y,Is(s)) = mean( ((td(s,:)-Xs(s)).^2 + (mi(s,:)-Ys(s)).^2).^(1/2) );
        end
        
        % distance of cluster centers to mean
        Sca(x,y) = mean( ((Xs-Xa).^2 + (Ys-Ya).^2).^(1/2) );
        
        % mean distance between centers
        Sc(x,y) = 0;
        count = 0;
        for s = 1:Nsites
            for ss = 1:s-1
                Sc(x,y) = Sc(x,y) + ((Xs(s)-Xs(ss))^2 + (Ys(s)-Ys(ss))^2)^(1/2);
                count = count + 1;
            end
        end
        Sc(x,y) = Sc(x,y)/count;
        
    end % y
end % x

nanmean((repmat(Ai,[1,1,Nm])-Mi)./repmat(Ai,[1,1,Nm]),3)
nanmean((repmat(Ai,[1,1,Ns])-Si)./repmat(Ai,[1,1,Ns]),3)

(nanmean(Mi,3)-nanmean(Si,3))./nanmean(Si,3)
(nanmean(Mi,3)-nanmean(Si,3))./Ai

(Mc-Sc)./Ai

Mca-Sca

%
figure(4);close(4);figure(4);
set(gcf,'color','w');

%bar((nanmean(Mi,3)-nanmean(Si,3))./nanmean(Si,3))
subplot(2,1,1)
a = nanmean(Mi,3)./Ai;
bar(a(:,7:9))
grid on;
legend(varNames{7:9},'location','sw');
set(gca,'xticklabel',varNames)
title('mean dist to model centers','fontsize',20);
set(gca,'ylim',[0,1]);
ylabel('Fracitonal Distance','fontsize',18);

subplot(2,1,2)
a = nanmean(Si,3)./Ai;
bar(a(:,7:9))
grid on;
set(gca,'xticklabel',varNames)
title('mean dist to site centers','fontsize',20);
set(gca,'ylim',[0,1]);
xlabel('Sending Variable','fontsize',18);
ylabel('Fracitonal Distance','fontsize',18);

% save sites figure
figure(4);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/Figure 5 - Clusters.png');





