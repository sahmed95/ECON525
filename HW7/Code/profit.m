function pi = profit(lambda, delta, phi, choice_prob,states, agent)
    pi = zeros(8,1);
    if agent == 1
        b = 2; 
    else 
        b = 1; 
    end
    for i=1:length(pi)
        m = states(i,3); 
        a1 = states(i,agent);
        a2 = states(i,b);
        ind1 = find(states(:,agent) == a2); 
        state1 = states(ind1,:); 
        ind2 = find(state1(:,b) == a1); 
        state2 = state1(ind2,:);
        ind3 = find(state2(:,3) == m);
        n = ind1(ind2(ind3));
        if agent == 2 
            n = i; 
        end 
        prob = choice_prob(n); 
        val1 = lambda*m-(1-a1)*phi; 
        val2 = val1-delta; 
        pi(i) = prob*val2+(1-prob)*val1; 
    end 
end
    