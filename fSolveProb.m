function difference=fSolveProb(females,males)

%% Parameters

incomeFemale=[1 2 3];
incomeMale=[2 3 4];

nGridFemale=length(incomeFemale);
nGridMale=length(incomeMale);

piFemale=[1/3 1/3 1/3];
piMale=[1/3 1/3 1/3];

bbeta=0.9;

mmu=1/2;

pf=females/sum(females);
pm=males/sum(males);

%% Critical Values

nGridFemale=length(incomeFemale);
nGridMale=length(incomeMale);

valueMarriageFemaleSecond=zeros(nGridFemale,nGridMale);
valueMarriageMaleSecond=zeros(nGridFemale,nGridMale);

for i=1:nGridFemale
    for j=1:nGridMale
        valueMarriageFemaleSecond(i,j)=fValueMarriageSecondFemale(incomeFemale(i),incomeMale(j),.5);
         valueMarriageMaleSecond(i,j)=fValueMarriageSecondMale(incomeFemale(i),incomeMale(j),.5);
    end
end

% Theta for Men and Mu for Women, when he has a lower reservation value

criticalMuFemale=zeros(nGridFemale,nGridMale);
tthetaCheckMale=zeros(nGridFemale,nGridMale);

for i=1:nGridFemale
    for j=1:nGridMale
        fCritical1=@(values) fCriticalValueCase1(values(1),values(2),i,j,incomeFemale,incomeMale,bbeta,pm,pf,...
            valueMarriageFemaleSecond,valueMarriageMaleSecond);
        v=fsolve(fCritical1,[mmu,0]);
        criticalMuFemale(i,j)=v(1);
        tthetaCheckMale(i,j)=v(2);
    end
end

% Theta for Women and Mu for men , when she has a lower reservation value

criticalMuMale=zeros(nGridFemale,nGridMale);
tthetaCheckFemale=zeros(nGridFemale,nGridMale);
 
for i=1:nGridFemale
    for j=1:nGridMale
        fCritical1=@(values) fCriticalValueCase2(values(1),values(2),i,j,incomeFemale,incomeMale,bbeta,pm,pf,...
            valueMarriageFemaleSecond,valueMarriageMaleSecond);
        v=fsolve(fCritical1,[mmu,0]);
        criticalMuMale(i,j)=v(1);
        tthetaCheckFemale(i,j)=v(2);
    end
end


%% Probability of a match

tthetaFemaleFirst=zeros(nGridFemale,nGridMale);
tthetaMaleFirst=zeros(nGridFemale,nGridMale);

for i=1:nGridFemale
    for j=1:nGridMale
        tthetaFemaleFirst(i,j)=log(incomeFemale(i)/(mmu*(incomeFemale(i)+incomeMale(j))))+...
    bbeta*(sum(valueMarriageFemaleSecond(i,:).*pm)-fValueMarriageSecondFemale(incomeFemale(i),incomeMale(j),...
    mmu));
        tthetaMaleFirst(i,j)=log(incomeMale(j)/((1-mmu)*(incomeFemale(i)+incomeMale(j))))+...
    bbeta*(sum(valueMarriageMaleSecond(:,j).*pf')-fValueMarriageSecondMale(incomeFemale(i),incomeMale(j),...
    mmu));
    end
end

probMatch=zeros(nGridFemale,nGridMale);

for i=1:nGridFemale
    for j=1:nGridMale
        if tthetaFemaleFirst(i,j)>=tthetaMaleFirst(i,j)
                probMatch(i,j)=exp(-tthetaCheckMale(i,j))/(1+exp(-tthetaCheckMale(i,j)));
        else
                probMatch(i,j)=exp(-tthetaCheckFemale(i,j))/(1+exp(-tthetaCheckFemale(i,j)));
        end
    end
end

singleFemale=zeros(nGridFemale,1);
singleMale=zeros(nGridMale,1);

for i=1:nGridFemale
    singleFemale(i)=(1- sum(probMatch(i,:).*piMale))*piFemale(i);
end

for j=1:nGridMale
    singleMale(j)=(1- sum(probMatch(:,j).*piFemale'))*piMale(j);
end

difference=sum(abs((singleFemale/sum(singleFemale))'-pf))+sum(abs((singleMale/sum(singleMale))'-pm));

end









