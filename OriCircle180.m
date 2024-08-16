% inputori ~ [-180 360] 

function outputori = OriCircle180(inputori)
    if inputori > 180
        outputori = inputori - 180;
    elseif inputori < 0
        outputori = inputori + 180;
    else
        outputori = inputori;
    end
end