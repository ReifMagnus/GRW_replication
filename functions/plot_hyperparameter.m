%% this script plots the optimized hyperparameters according to Giannone et al, Restat

f  = figure('Units','normalized','Position',[.1 .1 .3 .5]);
tl = tiledlayout(3,2,'Padding','tight','TileSpacing','tight');
nexttile(1:2);
histogram(P(:,1),'FaceColor','k','FaceAlpha',.5); hold on; grid on; box off;
title('$\lambda$','Interpreter','latex','FontSize',16)
nexttile;
histogram(P(:,2),'FaceColor','k','FaceAlpha',.5); hold on; grid on; box off;
title('$\chi_1$','Interpreter','latex','FontSize',16)
nexttile;
histogram(P(:,3),'FaceColor','k','FaceAlpha',.5); hold on; grid on; box off;
title('$\chi_2$','Interpreter','latex','FontSize',16)
nexttile;
histogram(P(:,4),'FaceColor','k','FaceAlpha',.5); hold on; grid on; box off;
title('$\chi_3$','Interpreter','latex','FontSize',16)
nexttile;
histogram(P(:,5),'FaceColor','k','FaceAlpha',.5); hold on; grid on; box off;
title('$\rho$','Interpreter','latex','FontSize',16)
if save_figs == 1
    exportgraphics(f,'Output/Revision/GLPvalues_dlIP_tlag.pdf'); %close
end
