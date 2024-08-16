% inputori ~ [0 180] 
function outputori = CalAngleDiff90(inputori)
    if inputori > 90
        outputori = 180-inputori;
    else
        outputori = inputori;
    end
end