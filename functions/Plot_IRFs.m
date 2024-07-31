DarkGrey  = [219 219 219]./255;
LightGrey = [235 235 235]./255;

LightBlue = [102,102,255]./255;
DarkBlue  = [051 051 255]./255;

LightRed = [255,153,153]./255;
DarkRed  = [255 0 0]./255;


hmaxtoplot = 19;
bands      = [16,50,84];

joint      = 1;         %=0 --> point-wise estimates; =1 --> joint estimates
ik         = 0;


[~,nVars,~,~] = size(Draws_IRFs);


p68 = prctile(Draws_IRFs(:,:,1:hmaxtoplot+1,:),bands,4);

if ~isempty(Draws_IRFs_narrative)
    p68_2 = prctile(squeeze(Draws_IRFs_narrative(:,:,1:hmaxtoplot+1,:)),bands,4);
end

IRFs_to_plot                              = Draws_BE;
% IRFs_to_plot(cumulateWhich,:,:)           = cumsum(IRFs_to_plot(cumulateWhich,:,:),3);
IRFs_to_plot_narrative                    = Draws_BE_narrative;
% IRFs_to_plot_narrative(cumulateWhich,:,:) = cumsum(IRFs_to_plot_narrative(cumulateWhich,:,:),3);

Bands_to_plot                      = Draws_Bands;
% Bands_to_plot(cumulateWhich,:,:,:) = cumsum(Bands_to_plot(cumulateWhich,:,:,:),3);

Bands_to_plot_narrative                      = Draws_Bands_narrative;
% Bands_to_plot_narrative(cumulateWhich,:,:,:) = cumsum(Bands_to_plot_narrative(cumulateWhich,:,:,:),3);

f  = figure('Units','normalized','Position',[.025 .1 .8 .825],'PaperType','A4');
tl = tiledlayout(nVars,nshocks,'TileSpacing','compact');
for i = 1:nVars
    for j = 1:nshocks
        nexttile; hold on
        xband = [0:hmaxtoplot hmaxtoplot:-1:0];            
        if joint == 1  % simulatanous confidence bands               
            if ik == 1  %           based on Inoue/Kilian
                p01   = plot(0:hmaxtoplot,squeeze(Draws_JCS(i,j,1:hmaxtoplot+1,1:1000)),'Color',LightRed,'Linewidth',1); hold on
%                 yband = [max(squeeze(Draws_JCS(ii,jj,1:hmaxtoplot+1,:)),[],2)' min(squeeze(Draws_JCS(ii,jj,hmaxtoplot+1:-1:1,:)),[],2)'];                        
            else  %           based on Montiel-Olea/Plagborg-Möller       
                yband = [squeeze(Bands_to_plot(i,j,1:hmaxtoplot+1,2))' squeeze(Bands_to_plot(i,j,hmaxtoplot+1:-1:1,1))']; 
                p11   = fill(xband,yband,LightRed,'FaceAlpha',.3,'EdgeColor',LightRed,'EdgeAlpha',.8); hold on
            end
        else        % pointwise bands
            p14  = plot(0:hmaxtoplot,squeeze(p68(i,j,:,[1 3])),'k--','Linewidth',1.5);
        end
        if  numel(unique([sign(min(yband)) sign(max(yband,[],'all'))])) > 1, zl = 1; end

        if ~isempty(Draws_IRFs_narrative)
            if joint == 1           % simulatanous confidence bands               
                if ik == 1  %           based on Inoue/Kilian
                    p01   = plot(0:hmaxtoplot,squeeze(Draws_JCS_narrative(i,j,1:hmaxtoplot+1,:)),'Color',LightBlue,'Linewidth',1);   
%                     yband = [max(squeeze(Draws_JCS_narrative(i,j,1:hmaxtoplot+1,:)),[],2)' min(squeeze(Draws_JCS_narrative(i,j,hmaxtoplot+1:-1:1,:)),[],2)'];        
                else        %            based on Montiel-Olea/Plagborg-Möller   
                    yband = [squeeze(Bands_to_plot_narrative(i,j,1:hmaxtoplot+1,2))' squeeze(Bands_to_plot_narrative(i,j,hmaxtoplot+1:-1:1,1))']; 
                    p21   = fill(xband,yband,LightBlue,'FaceAlpha',.3,'EdgeColor',LightBlue,'EdgeAlpha',.8); hold on    
                end
            else  % pointwise bands               
                p24  = plot(0:hmaxtoplot,squeeze(p68_2(i,j,:,[1 3])),'b--','Linewidth',1.5);
            end
        end

        % pointwise means
        p13  = plot(0:hmaxtoplot,squeeze(p68(i,j,:,2)),'Color',DarkRed,'Linewidth',1.5,'Linestyle','--'); 
        p23  = plot(0:hmaxtoplot,squeeze(p68_2(i,j,:,2)),'Color',DarkBlue,'Linewidth',1.5,'Linestyle','--'); 
        
        
        p12 = plot(0:hmaxtoplot,squeeze(IRFs_to_plot(i,j,1:hmaxtoplot+1)),'Color',DarkRed,'Linewidth',3,'LineStyle','-'); hold on             % Bayes estimates
        p22 = plot(0:hmaxtoplot,squeeze(IRFs_to_plot_narrative(i,j,1:hmaxtoplot+1)),'Color',DarkBlue,'Linewidth',3,'LineStyle','-.'); hold on            % Bayes estimates
        
        if  numel(unique([sign(min(yband)) sign(max(yband,[],'all'))])) > 1, zl = 1; else zl = 0; end

        if zl == 1
            yline(0,'k','LineWidth',2)
        end
        set(gca,'FontName','Times','FontSize',14,'XTick',(0:round(hmaxtoplot/4.33):hmaxtoplot)');
        
        axis tight; box off; grid on;
        if i == 1, title(shockNames(j),'FontName','Times','FontSize',14); end
        if j == 1, ylabel(varNames(i),'FontName','Times','FontSize',14,'FontWeight','normal'); end

        set(gca,'XTickLabel',num2str((0:round(hmaxtoplot/4.33):hmaxtoplot)'))
    end
end
nexttile(nshocks*nVars-nVars+1); xlabel('Months','FontName','Times','FontSize',14);
nexttile(nshocks*nVars-nVars+2); xlabel('Months','FontName','Times','FontSize',14);
nexttile(nshocks*nVars-nVars+3); xlabel('Months','FontName','Times','FontSize',14);
nexttile(nshocks*nVars-nVars+4); xlabel('Months','FontName','Times','FontSize',14);
% l = legend([p12 p22 p13 p23],'Bayes estimate (signs)','Bayes estimate (signs & narratives)','Posterior median (signs)','Posterior median (signs & narratives)');
l = legend([p12 p22 p13 p23],'Bayes estimate (signs)','Bayes estimate (narratives)','Posterior median (signs)','Posterior median (narratives)');
l.Box = 'off'; l.FontSize = 14;  l.NumColumns = 4;
rect = [0.1, 0.001, .8 .05];
set(l, 'Position', rect)

if save_figs == 1
    exportgraphics(f,'Output/Revision/IRFs_dlIP_tlag.pdf'); %close
end


