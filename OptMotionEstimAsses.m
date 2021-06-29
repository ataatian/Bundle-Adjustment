
clear
clc

load AliGt myXgt myYgt myXsongt myYsongt
load Kmatrix K
W= 60;
load Points3Dreffile Points3Dref
Points3Dref= Points3Dref(:,1:W);
ScaleEstimate= 0.5199;

myXg= [myXgt(:,1) myXgt(:,4)];
myYg= [myYgt(:,1) myYgt(:,4)];
NumPoints= size(myXg,1);
NoiseLevel= 2;
ANGstatGT= [   -0.0772
                0.3656
               -0.0344];
TstatGT= [    -0.4630
              -0.1830
               0.1499];
for sample= 1:100
    for j= 1:2
        ok=1;
        while(ok) 
            XNoise= randn(NumPoints,1)*NoiseLevel;
% %             figure(111)
% %             subplot(1,2,1)
% %             hist(XNoise,13);
            if (abs(mean(XNoise(:)))<0.01) & (abs(std(XNoise(:))-NoiseLevel)<0.02)
                ok=0;
            end
        end
        ok=1;
        while(ok) 
            YNoise= randn(NumPoints,1)*NoiseLevel;
% %             subplot(1,2,2)
% %             hist(YNoise,13);
            if (abs(mean(YNoise(:)))<0.01) & (abs(std(YNoise(:))-NoiseLevel)<0.02)
                ok=0;
            end
        end
        myX(:,j)= myXg(:,j)+XNoise;
        myY(:,j)= myYg(:,j)+YNoise;
    end
    [R_E,T_E,R_G,T_G,N]= CamMotionExNormLastSINGLEforStatistics(myX,myY,K);
    ANGstat(:,sample)= rodrigues(R_E);
    Tstat(:,sample)=   T_E/norm(T_E)*ScaleEstimate;
    pix1= [myX(:,1)';myY(:,1)';ones(1, NumPoints)];
    pix2= [myX(:,2)';myY(:,2)';ones(1, NumPoints)];
    pix1= inv(K)*pix1;
    pix2= inv(K)*pix2;
    pix1= [pix1(1,:)./pix1(3,:);pix1(2,:)./pix1(3,:);pix1(3,:)./pix1(3,:)];
    pix2= [pix2(1,:)./pix2(3,:);pix2(2,:)./pix2(3,:);pix2(3,:)./pix2(3,:)];
    [XP_E, lambda_E]= compute3DStructure(pix1(1:2,1:W), pix2(1:2,1:W),R_E, Tstat(:,sample));
    Points3Dstat(:,:,sample)= XP_E(:,:,1);
end
STDTstat= [std(Tstat')]'
STDANGstat= [std((ANGstat)')]'
ErrmeanTstat= [mean(Tstat')]'-TstatGT
ErrmeanANGstat= [mean(ANGstat')]'-ANGstatGT

for i= 1:W
    setPoints2(:,:)= Points3Dstat(:,i,:);
    STDPoints(:,i)= std(setPoints2');
    V= [mean(setPoints2')]';
    ErrmeanPoints(:,i)= V(1:3)  - Points3Dref(:,i);
end
STDPoints(1:3,:)

ErrmeanPoints(1:3,:)

MeanSTDall= [mean(STDPoints(1:3,:)')]'
MeanErrall= [mean(ErrmeanPoints(1:3,:)')]'




W= 60;
load Points3Dreffile Points3Dref
Points3Dref= Points3Dref(:,1:W);
figure
plot3(Points3Dref(1,:),Points3Dref(2,:),Points3Dref(3,:),'.b'),hold on
plot3(0,0,0,'+k');
    M= [Points3Dref(:,4) Points3Dref(:,22) Points3Dref(:,30)]';
    MM= [M [+1;+1;+1]];
    myRank= rank(MM);
    [U,S,V] = svd(MM);
    MYno1= V(:,4)/V(4,4);
    Theno1= MYno1(1:3);

    dist= 1/norm(Theno1);
    NormalVector= Theno1*dist^2;
    plot3([0 NormalVector(1)],[0 NormalVector(2)],[0 NormalVector(3)])
    
    
    
    
    