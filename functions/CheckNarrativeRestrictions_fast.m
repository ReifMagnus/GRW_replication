function [check] = CheckNarrativeRestrictions_fast(IRFs,n,shocks,narrative_signs,narrative_contributions,startperiod,endperiod,shock,variable,findnarrativecontrib)
%% CheckNarrativeRestrictions: Checks whether narrative restrictions are satisfied

narrative_restrictions       = narrative_signs(narrative_signs~=0);
check_narrative_sign         = (sum(narrative_restrictions == sign(shocks(narrative_signs~=0))) == length(narrative_signs(narrative_signs~=0)));
check_narrative_contrib_temp = zeros(length(findnarrativecontrib),1);

for restriction = 1:length(findnarrativecontrib)

    Tstart                  = startperiod(restriction);
    Tend                    = endperiod(restriction);
    [~,contributions_restr] = getHDs_fast(IRFs,n,shocks(Tstart:Tend,:),variable(restriction));

    switch narrative_contributions(findnarrativecontrib(restriction))
        case 1 % Soft Restriction (The shock is the most important contributor)
            check_narrative_contrib_temp(restriction) =  (abs(contributions_restr(shock(restriction))) == max(abs(contributions_restr)));
        case -1 % Soft Restriction (The shock is the least important contributor)
            check_narrative_contrib_temp(restriction) =  (abs(contributions_restr(shock(restriction))) == min(abs(contributions_restr)));
        case 2 % Soft Restriction (The shock is the overwhelming contributor)
            check_narrative_contrib_temp(restriction) =  (abs(contributions_restr(shock(restriction))) > sum(abs(contributions_restr(1:end ~= shock(restriction)))));
        case -2 % Soft Restriction (The shock is not the overwhelming contributor)
            check_narrative_contrib_temp(restriction) =  (abs(contributions_restr(shock(restriction))) < sum(abs(contributions_restr(1:end ~= shock(restriction)))));
        case 3 % cumulative impact is negative 
            check_narrative_contrib_temp(restriction) =  contributions_restr(shock(restriction)) < 0;
        case 4 % Soft Restriction (The shock is not the least important contributor)
            check_narrative_contrib_temp(restriction) =  (abs(contributions_restr(shock(restriction))) ~= min(abs(contributions_restr)));            
    end
end

check_narrative_contrib = (sum(check_narrative_contrib_temp) == length(findnarrativecontrib));
check = check_narrative_sign*check_narrative_contrib;

end

