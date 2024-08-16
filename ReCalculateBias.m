for trial_i = 1:150   
 % calculate bias
    if results(trial_i,6) == 1  
        results(trial_i,21) = CalAngleDiff90(abs(results(trial_i,4)-results(trial_i,8)));       
        results(trial_i,22) = CalAngleDiff90(abs(results(trial_i,5)-results(trial_i,8)));
        results(trial_i,23) = CalAngleDiff90(abs(results(trial_i,5)-results(trial_i,9)));
        results(trial_i,24) = CalAngleDiff90(abs(results(trial_i,4)-results(trial_i,9)));
    else
        results(trial_i,21) = CalAngleDiff90(abs(results(trial_i,4)-results(trial_i,9)));       
        results(trial_i,22) = CalAngleDiff90(abs(results(trial_i,5)-results(trial_i,9)));
        results(trial_i,23) = CalAngleDiff90(abs(results(trial_i,5)-results(trial_i,8)));
        results(trial_i,24) = CalAngleDiff90(abs(results(trial_i,4)-results(trial_i,8)));
    end

    results(trial_i,25) = CalAngleDiff90(abs(results(trial_i,4)-results(trial_i,5)));
    results(trial_i,30) = CalAngleDiff90(abs(results(trial_i,8)-results(trial_i,9)));
    results(trial_i,31) = results(trial_i,30)-results(trial_i,25);

    if (results(trial_i,21) + results(trial_i,22)) - results(trial_i,25) > 0.000001
        results(trial_i,16) = results(trial_i,21);
    else
        results(trial_i,16) = -results(trial_i,21);
    end

    if (results(trial_i,23) + results(trial_i,24)) - results(trial_i,25) > 0.000001
        results(trial_i,17) = results(trial_i,23);
    else
        results(trial_i,17) = -results(trial_i,23);
    end

    % remember wrong order
    if (results(trial_i,22) < results(trial_i,21)) && (results(trial_i,24) < results(trial_i,23))
        results(trial_i,26) = 1;
    end

    % remember wrong order
    if results(trial_i,8)>=results(trial_i,4) && results(trial_i,8)>=results(trial_i,5) && results(trial_i,9)>=results(trial_i,5) && results(trial_i,9)>=results(trial_i,4)
        results(trial_i,26) = 1;
    end

    % remember wrong order
    if results(trial_i,8)<=results(trial_i,4) && results(trial_i,8)<=results(trial_i,5) && results(trial_i,9)<=results(trial_i,5) && results(trial_i,9)<=results(trial_i,4)
        results(trial_i,26) = 1;
    end

    if results(trial_i,6) == 1 
        results(trial_i,27) = results(trial_i,16);
        results(trial_i,28) = results(trial_i,17);
    else
        results(trial_i,27) = results(trial_i,17);
        results(trial_i,28) = results(trial_i,16);
    end

end