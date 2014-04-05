function valueMale=fValueMarriageSecondMale(yf,ym,mmu)

thresholdFemaleSecond= max(log(yf/(mmu*(yf+ym))),-.5);
thresholdMaleSecond= max(log(ym/((1-mmu)*(yf+ym))),-.5);

int1=@(ttheta) ttheta.*exp(-ttheta)./((1+exp(-ttheta)).^2);
int2=@(ttheta) log(ym+yf*(1-exp(-ttheta))).*exp(-ttheta)./(1+exp(-ttheta)).^2;

if thresholdFemaleSecond>=thresholdMaleSecond
        valueMale=log((1-mmu)*(yf+ym))*exp(-thresholdFemaleSecond)/(1+exp(-thresholdFemaleSecond))+...
        quadgk(int1,0,Inf)+...
        quadgk(int2,0,thresholdFemaleSecond)+log(ym)/2;
else
        valueMale=log((1-mmu)*(yf+ym))*exp(-thresholdMaleSecond)/(1+exp(-thresholdMaleSecond))+...
        quadgk(int1,thresholdMaleSecond,Inf)+log(ym)/(1+exp(-thresholdMaleSecond));
end
  
end

