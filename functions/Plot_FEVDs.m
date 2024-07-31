%% Forecast Error Variance Decompositions With Bands
bands      = [16,50,84];

numSavedDraws     = size(A0_save,3);
numSavedNarrative = size(A0_narrative,3);

hmax = 100;

% draws_FEVDs           = nan(n,n,hmaxtoplot+1,numSavedDraws); % variable,shock,horizon,draw
% draws_FEVDs_narrative = nan(n,n,hmaxtoplot+1,numSavedNarrative); % variable,shock,horizon,draw

% for draw = 1:numSavedDraws
%     IRFs                    = getIRFs(Beta_save(:,:,draw),A0_save(:,:,draw),exog,n,p,hmax);
%     draws_FEVDs(:,:,:,draw) = varianceDecompositionOfVAR(IRFs,hmaxtoplot);
% end
% 
% if ~isempty(Beta_narrative)
%     for draw = 1:numSavedNarrative
%         IRFs                              = getIRFs(Beta_narrative(:,:,draw),A0_narrative(:,:,draw),exog,n,p,hmax);
%         draws_FEVDs_narrative(:,:,:,draw) = varianceDecompositionOfVAR(IRFs,hmaxtoplot);
%     end
% end
%-----------------------------------------------------------------------------------------------
% based on joint posterior

aux      = getIRFs(Beta_save(:,:,I),A0_save(:,:,I),exog,n,p,hmax);   
FEVD_BE  = varianceDecompositionOfVAR(aux,hmax).*100;         % bayes estimator

for i = 1:size(Draws_JCS,4)
    aux             = getIRFs(Beta_save(:,:,S(i)),A0_save(:,:,S(i)),exog,n,p,hmax);   
    FEVD_set(:,:,:,i) = varianceDecompositionOfVAR(aux,hmax).*100;         
end

FEVD_p16 = min(FEVD_set,[],4);
FEVD_p84 = max(FEVD_set,[],4);


for i = 1:size(Draws_JCS_narrative,4)
    aux                         = getIRFs(Beta_narrative(:,:,Sn(i)),A0_narrative(:,:,Sn(i)),exog,n,p,hmax);   
    FEVD_set_narrative(:,:,:,i) = varianceDecompositionOfVAR(aux,hmax).*100;        
end

aux    = getIRFs(Beta_narrative(:,:,In),A0_narrative(:,:,In),exog,n,p,hmax);
FEVD_n = varianceDecompositionOfVAR(aux,hmax).*100;               % bayes estimator
        
FEVD_n_p16 = min(FEVD_set_narrative,[],4);
FEVD_n_p84 = max(FEVD_set_narrative,[],4);

%% Produce FEVD table for Latex
%%
fprintf('\\begin{table}[t]\n')
fprintf('\\centering\n')
fprintf('\\small\n')
fprintf('\\caption{Contribution of structural shocks to FEVD (in \\%%)}\\label{table:fevd}\n')
fprintf('\\begin{tabularx}{\\textwidth}{mnnnn}\n')
fprintf('\\hline\\hline\\noalign{\\smallskip}\n')
fprintf('& \\textbf{Flow supply\\newline shock} & \\textbf{Flow \\newline demand shock}  & \\textbf{Storage demand\\newline shock}  & \\textbf{Preference\\newline shock}  \\\\\n')
fprintf('\\cline{2-5} \\noalign{\\smallskip}\\noalign{\\smallskip}\n')
fprintf('Gas net import  &  %.1f &  %.1f  &  %.1f  &  %.1f  \\tabularnewline\n',FEVD_n(1,1,100),FEVD_n(1,2,100),FEVD_n(1,3,100),FEVD_n(1,4,100))
fprintf('growth &  [%.1f, %.1f] &  [%.1f, %.1f]  &  [%.1f, %.1f]  &  [%.1f, %.1f]  \\tabularnewline\n',FEVD_n_p16(1,1,100),FEVD_n_p84(1,1,100),FEVD_n_p16(1,2,100),FEVD_n_p84(1,2,100),FEVD_n_p16(1,3,100),FEVD_n_p84(1,3,100),FEVD_n_p16(1,4,100),FEVD_n_p84(1,4,100))
fprintf('Industrial production &  %.1f  &  %.1f  &  %.1f  &  %.1f \\tabularnewline\n',FEVD_n(2,1,100),FEVD_n(2,2,100),FEVD_n(2,3,100),FEVD_n(2,4,100))
fprintf('growth &  [%.1f, %.1f] &  [%.1f, %.1f]  &  [%.1f, %.1f]  &  [%.1f, %.1f]  \\tabularnewline\n',FEVD_n_p16(2,1,100),FEVD_n_p84(2,1,100),FEVD_n_p16(2,2,100),FEVD_n_p84(2,2,100),FEVD_n_p16(2,3,100),FEVD_n_p84(2,3,100),FEVD_n_p16(2,4,100),FEVD_n_p84(2,4,100))
fprintf('Real gas price  &  %.1f &  %.1f  &  %.1f  &  %.1f \\tabularnewline\n',FEVD_n(3,1,100),FEVD_n(3,2,100),FEVD_n(3,3,100),FEVD_n(3,4,100))
fprintf('growth &  [%.1f, %.1f] &  [%.1f, %.1f]  &  [%.1f, %.1f]  &  [%.1f, %.1f]  \\tabularnewline\n',FEVD_n_p16(2,1,100),FEVD_n_p84(2,1,100),FEVD_n_p16(2,2,100),FEVD_n_p84(3,2,100),FEVD_n_p16(3,3,100),FEVD_n_p84(3,3,100),FEVD_n_p16(3,4,100),FEVD_n_p84(3,4,100))
fprintf('\\multirow{2}{*}{Gas inventories} &  %.1f &  %.1f  &  %.1f  &  %.1f \\tabularnewline\n',FEVD_n(4,1,100),FEVD_n(4,2,100),FEVD_n(4,3,100),FEVD_n(4,4,100))
fprintf(' &  [%.1f, %.1f] &  [%.1f, %.1f]  &  [%.1f, %.1f]  &  [%.1f, %.1f]  \\tabularnewline\n',FEVD_n_p16(4,1,100),FEVD_n_p84(4,1,100),FEVD_n_p16(4,2,100),FEVD_n_p84(4,2,100),FEVD_n_p16(4,3,100),FEVD_n_p84(4,3,100),FEVD_n_p16(4,4,100),FEVD_n_p84(4,4,100))
fprintf('\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n')
fprintf('\\end{tabularx}\n')
fprintf('\\parbox{\\textwidth}{\\footnotesize \\textbf{Notes:} Variance decomposition based on Bayes estimate of impulse responses in Figure \\ref{figure:IRFs} for models satisfying both conventional and narrative sign restrictions. 68\\%% error bands in brackets. Unconditional variances are approximated by setting the forecast horizon to $h=100$ months.}\n')
fprintf('\\end{table}\n')
