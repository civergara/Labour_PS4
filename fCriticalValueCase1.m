function F=fCriticalValueCase1(mmu,ttheta,i,j,incomeFemale,incomeMale,bbeta,pm,pf,valueMarriageFemaleSecond,...
    valueMarriageMaleSecond)

F=[log(mmu)+bbeta*fValueMarriageSecondFemale(incomeFemale(i),incomeFemale(j),mmu)-log(incomeFemale(i)/...
    (incomeFemale(i)+incomeMale(j)))+ttheta-bbeta*sum(valueMarriageFemaleSecond(i,:).*pm);...
    log(incomeMale(j)/((1-mmu)*(incomeFemale(i)+incomeMale(j))))+...
    bbeta*(sum(valueMarriageMaleSecond(:,j).*pf')-fValueMarriageSecondMale(incomeFemale(i),incomeMale(j),...
    mmu))]; 

end
