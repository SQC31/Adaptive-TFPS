function [l1,l2,l3]=myfind(A)
[l10,l20,l30]=find(~A);
[l11,l21,l31]=find(A);
if size(A,1)==1
    l1=[l10,l11]'; l2=[l20,l21]'; l3=[1-l30,l31]';
else
    l1=[l10;l11]; l2=[l20;l21]; l3=[1-l30;l31];
end

end