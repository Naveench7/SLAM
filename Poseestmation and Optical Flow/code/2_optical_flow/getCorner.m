function res = getCorner(id)
%% CHANGE THE NAME OF THE FUNCTION TO getCorner
    %% Input Parameter Description
    % id = List of all the AprilTag ids detected in the current image(data)
    
    %% Output Parameter Description
    % res = List of the coordinates of the 4 corners (or 4 corners and the
    % centre) of each detected AprilTag in the image in a systematic method
     Tag=[0, 12, 24, 36, 48, 60, 72, 84,  96;1, 13, 25, 37, 49, 61, 73, 85,  97;2, 14, 26, 38, 50, 62, 74, 86,  98;3, 15, 27, 39, 51, 63, 75, 87,  99;4, 16, 28, 40, 52, 64, 76, 88, 100;5, 17, 29, 41, 53, 65, 77, 89, 101;6, 18, 30, 42, 54, 66, 78, 90, 102;7, 19, 31, 43, 55, 67, 79, 91, 103;8, 20, 32, 44, 56, 68, 80, 92, 104;9, 21, 33, 45, 57, 69, 81, 93, 105;10, 22, 34, 46, 58, 70, 82, 94, 106;11, 23, 35, 47, 59, 71, 83, 95, 107];
    x=0;
    fptlx=zeros(1,12);
    for i=0:11
        if i==0
            fptlx(i+1)=x;
        else
            x=x+2*0.152;
            fptlx(i+1)=x;
        end
    end
    fptly=zeros(size(fptlx));   
    fpblx=0.152+fptlx;
    fpbly=fptly;
    fptrx=fptlx;
    fptry=0.152+fptly;
    fpbrx=fpblx;
    fpbry=0.152+fpbly;
    fptl=cat(1,fptlx,fptly);
    fptr=cat(1,fptrx,fptry);
    fpbl=cat(1,fpblx,fpbly);
    fpbr=cat(1,fpbrx,fpbry);
    sptl=cat(1,fptlx,2*0.152+fptly);
    sptr=cat(1,fptlx,3*0.152+fptly);
    spbl=cat(1,fpblx,2*0.152+fpbly);
    spbr=cat(1,fpblx,3*0.152+fpbly);
    tptl=cat(1,fptlx,4*0.152+fptly);
    tptr=cat(1,fptlx,5*0.152+fptly);
    tpbl=cat(1,fpblx,4*0.152+fpbly);
    tpbr=cat(1,fpblx,5*0.152+fpbly);
    frptl=cat(1,fptlx,5*0.152+fptly+0.178);
    frptr=cat(1,fptlx,6*0.152+fptly+0.178);
    frpbl=cat(1,fpblx,5*0.152+fpbly+0.178);
    frpbr=cat(1,fpblx,6*0.152+fpbly+0.178);
    fiptl=cat(1,fptlx,7*0.152+fptly+0.178);
    fiptr=cat(1,fptrx,8*0.152+fptly+0.178);
    fipbl=cat(1,fpblx,7*0.152+fpbly+0.178);
    fipbr=cat(1,fpblx,8*0.152+fpbly+0.178);
    siptl=cat(1,fptlx,9*0.152+fptly+0.178);
    siptr=cat(1,fptlx,10*0.152+fptly+0.178);
    sipbl=cat(1,fpblx,9*0.152+fpbly+0.178);
    sipbr=cat(1,fpblx,10*0.152+fpbly+0.178);
    septl=cat(1,fptlx,10*0.152+fptly+0.178+0.178);
    septr=cat(1,fptlx,11*0.152+fptly+0.178+0.178);
    sepbl=cat(1,fpblx,10*0.152+fpbly+0.178+0.178);
    sepbr=cat(1,fpblx,11*0.152+fpbly+0.178+0.178);
    eptl=cat(1,fptlx,12*0.152+fptly+0.178+0.178);
    eptr=cat(1,fptlx,13*0.152+fptly+0.178+0.178);
    epbl=cat(1,fpblx,12*0.152+fpbly+0.178+0.178);
    epbr=cat(1,fpblx,13*0.152+fpbly+0.178+0.178);
    nptl=cat(1,fptlx,14*0.152+fptly+0.178+0.178);
    nptr=cat(1,fptlx,15*0.152+fptly+0.178+0.178);
    npbl=cat(1,fpblx,14*0.152+fpbly+0.178+0.178);
    npbr=cat(1,fpblx,15*0.152+fpbly+0.178+0.178);
    pfourx=zeros(length(id),1);
    pfoury=zeros(length(id),1);
    pthreex=zeros(length(id),1);
    pthreey=zeros(length(id),1);
    ptwox=zeros(length(id),1);
    ptwoy=zeros(length(id),1);
    ponex=zeros(length(id),1);
    poney=zeros(length(id),1);
    for j=1:length(id)
        [r,c]=find(Tag==id(j));
        if c==1
            pfourx(j,1)=fptl(1,r);
            pfoury(j,1)=fptl(2,r);
            pthreex(j,1)=fptr(1,r);
            pthreey(j,1)=fptr(2,r);
            ptwox(j,1)=fpbr(1,r);
            ptwoy(j,1)=fpbr(2,r);
            ponex(j,1)=fpbl(1,r);
            poney(j,1)=fpbl(2,r);
        elseif c==2
            pfourx(j,1)=sptl(1,r);
            pfoury(j,1)=sptl(2,r);
            pthreex(j,1)=sptr(1,r);
            pthreey(j,1)=sptr(2,r);
            ptwox(j,1)=spbr(1,r);
            ptwoy(j,1)=spbr(2,r);
            ponex(j,1)=spbl(1,r);
            poney(j,1)=spbl(2,r);
        elseif c==3
            pfourx(j,1)=tptl(1,r);
            pfoury(j,1)=tptl(2,r);
            pthreex(j,1)=tptr(1,r);
            pthreey(j,1)=tptr(2,r);
            ptwox(j,1)=tpbr(1,r);
            ptwoy(j,1)=tpbr(2,r);
            ponex(j,1)=tpbl(1,r);
            poney(j,1)=tpbl(2,r);
        elseif c==4
            pfourx(j,1)=frptl(1,r);
            pfoury(j,1)=frptl(2,r);
            pthreex(j,1)=frptr(1,r);
            pthreey(j,1)=frptr(2,r);
            ptwox(j,1)=frpbr(1,r);
            ptwoy(j,1)=frpbr(2,r);
            ponex(j,1)=frpbl(1,r);
            poney(j,1)=frpbl(2,r);
        elseif c==5
            pfourx(j,1)=fiptl(1,r);
            pfoury(j,1)=fiptl(2,r);
            pthreex(j,1)=fiptr(1,r);
            pthreey(j,1)=fiptr(2,r);
            ptwox(j,1)=fipbr(1,r);
            ptwoy(j,1)=fipbr(2,r);
            ponex(j,1)=fipbl(1,r);
            poney(j,1)=fipbl(2,r);
        elseif c==6
            pfourx(j,1)=siptl(1,r);
            pfoury(j,1)=siptl(2,r);
            pthreex(j,1)=siptr(1,r);
            pthreey(j,1)=siptr(2,r);
            ptwox(j,1)=sipbr(1,r);
            ptwoy(j,1)=sipbr(2,r);
            ponex(j,1)=sipbl(1,r);
            poney(j,1)=sipbl(2,r);
        elseif c==7
            pfourx(j,1)=septl(1,r);
            pfoury(j,1)=septl(2,r);
            pthreex(j,1)=septr(1,r);
            pthreey(j,1)=septr(2,r);
            ptwox(j,1)=sepbr(1,r);
            ptwoy(j,1)=sepbr(2,r);
            ponex(j,1)=sepbl(1,r);
            poney(j,1)=sepbl(2,r);
        elseif c==8
            pfourx(j,1)=eptl(1,r);
            pfoury(j,1)=eptl(2,r);
            pthreex(j,1)=eptr(1,r);
            pthreey(j,1)=eptr(2,r);
            ptwox(j,1)=epbr(1,r);
            ptwoy(j,1)=epbr(2,r);
            ponex(j,1)=epbl(1,r);
            poney(j,1)=epbl(2,r);
        elseif c==9
            pfourx(j,1)=nptl(1,r);
            pfoury(j,1)=nptl(2,r);
            pthreex(j,1)=nptr(1,r);
            pthreey(j,1)=nptr(2,r);
            ptwox(j,1)=npbr(1,r);
            ptwoy(j,1)=npbr(2,r);
            ponex(j,1)=npbl(1,r);
            poney(j,1)=npbl(2,r);
        end
    end
    pzerox=pfourx+0.076;
    pzeroy=pfoury+0.076;
    res=[pzerox pzeroy;ponex poney;ptwox ptwoy;pthreex pthreey;pfourx pfoury];
end