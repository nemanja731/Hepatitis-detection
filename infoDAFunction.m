function InfoDA = infoDAFunction(InfoD, X, Y, type)
    numSamples = length(X);
    numSteps = 100;
    switch type
        case 'D  '
        case 'AGE'
            up = 78;
            low = 10;
            step = (up - low)/numSteps;
            X = floor(X/step)*step;
        case 'BIL'
            up = 4.9;
            low = 0.3;
            step = (up - low)/numSteps;
            X = floor(X/step)*step;
        case 'ALK'
            up = 187;
            low = 26;
            step = (up - low)/numSteps;
            X = floor(X/step)*step;
        case 'SGO'
            up = 331;
            low = 14;
            step = (up - low)/numSteps;
            X = floor(X/step)*step;
        case 'ALB'
            up = 5.1;
            low = 2.1;
            step = (up - low)/numSteps;
            X = floor(X/step)*step;
        case 'PRO'
            up = 100;
            low = 20;
            step = (up - low)/numSteps;
            X = floor(X/step)*step;
    end
    values = unique(X);
    InfoDA = 0;
    for v = 1:length(values)
        Z = Y(X == values(v));
        Z(Z == 2) = 0;
        InfoDA = InfoDA + sum(X == values(v))/numSamples*entropy(Z);
    end
    InfoDA = InfoD - InfoDA;
end