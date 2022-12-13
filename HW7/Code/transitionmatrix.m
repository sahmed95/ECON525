function F = transitionmatrix(choice_prob,B, states)
   F = zeros(8,8);
   for i =1:length(choice_prob)
        for j=1:length(choice_prob)  
            m_1 = states(i, 3); 
            m_2 = states(j,3); 
            a1_old = states(i,1);
            a2_old = states(i,2);
            a1_new = states(j,1); 
            a2_new = states(j,2); 
            if m_1 == 1
                k = 1;
            else 
                k = 2;
            end
            if m_2 == 1
                l = 1; 
            else 
                l = 2; 
            end 
            ind1 = find(states(:,1) == a2_old); 
            state1 = states(ind1,:); 
            ind2 = find(state1(:,2) == a1_old); 
            state2 = state1(ind2,:);
            ind3 = find(state2(:,3) == m_1);
            n = ind1(ind2(ind3));
            if a1_new == 0
                prob_1 = 1-choice_prob(i); 
            else 
                prob_1 = choice_prob(i); 
            end
            if a2_new == 0
                prob_2 = 1-choice_prob(n); 
            else 
                prob_2 = choice_prob(n);
            end
            F(i,j) = B(k,l)*prob_1*prob_2; 
        end 
    end 
end