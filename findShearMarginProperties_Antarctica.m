function [Amat,gsizemat,tmpmat] = findShearMarginProperties_Antarctica(SRmat,theta,Tmat,SMBmat,p,D,n,height_ratio)

Amat = zeros(size(SRmat));
gsizemat = zeros(size(SRmat));
tmpmat = zeros(size(SRmat));
for i=1:size(SRmat,1)
    for j=1:size(SRmat,2)
        elat = SRmat(i,j);
        H = Tmat(i,j);
        smb = SMBmat(i,j);
        [Amat(i,j),gsizemat(i,j),tmpmat(i,j)] = findPropertiesofShearMargin(elat,theta,H,p,D,n,smb,height_ratio);
    end
    fprintf('Iteration %d of %d done \n',i,size(SRmat,1))
end

end