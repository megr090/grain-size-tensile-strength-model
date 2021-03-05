function [viscositymat,Amat,zeta,reduc,gsizemat,tmpmat,z,H] = findShearMarginProperties_Antarctica_DProf(SRmat,festmat,Tmat,SMBmat,p,D,n,i,j)

f = festmat(i,j);
elat = SRmat(i,j);
H = Tmat(i,j);
smb = SMBmat(i,j);
if H<1
    z = [0:1];
elseif H<100
    z = [0:1:H]; % m
else
    z = [0:1:H]; % m
end
viscositymat = zeros(length(z));
Amat = zeros(length(z));
gsizemat = zeros(length(z));
tmpmat = zeros(length(z));

[viscositymat,Amat,zeta,reduc,gsizemat,tmpmat] = findPropertiesofShearMargin_DProf(elat,f,H,p,D,n,smb,z);

end