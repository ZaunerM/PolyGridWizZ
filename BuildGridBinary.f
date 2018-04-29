c----------------------------------------------------------------------
      program Intake_Gridgen
c----------------------------------------------------------------------
c

      implicit none

      integer, parameter :: rk = 8!selected_real_kind(10,100)

      integer :: i,j,k,xb,zb,recnum,nxt,nzt,ntot,ifi

      integer :: xm, xp, ym, yp, zm, zp, npt3, npt1

      real(kind=rk) :: xpr,ypr,zpr, dummy

      integer, parameter :: nbX = 3,nbZ = 1
      integer :: bcxm(nbX,nbZ),bcxp(nbX,nbZ),
     &           bcym(nbX,nbZ),bcyp(nbX,nbZ),
     &           bczm(nbX,nbZ),bczp(nbX,nbZ)


      integer :: nxb(nbX,nbZ), nyb(nbX,nbZ),nzb(nbX,nbZ),
     &           xmintf(nbX,nbZ),xpintf(nbX,nbZ),
     &           zmintf(nbX,nbZ),zpintf(nbX,nbZ)

      integer, parameter :: halo=3 !nplus=355,itmax=25 !!nz=1

      integer :: ib

      integer :: nxpb(nbX,nbZ),nypb(nbX,nbZ),nzpb(nbX,nbZ)

      real(kind=rk) :: xlb(nbX,nbZ),ylb(nbX,nbZ),zlb(nbX,nbZ)

      real(kind=rk), allocatable :: xx(:,:), yy(:,:)
      real(kind=rk), allocatable :: xx1(:,:), yy1(:,:)
      real(kind=rk), allocatable :: xx2(:,:), yy2(:,:)
      real(kind=rk), allocatable :: xx3(:,:), yy3(:,:)
      real(kind=rk), allocatable :: tmp1(:,:), tmp2(:,:)
      
      character(75) :: filename, outdir, tail

      outdir ='./'
      
! Set flow info and boundary conditions
      xm=1;xp=2;ym=3;yp=4;zm=5;zp=6
      
      bcxm(1,1) = 1 ! Interface
      bcxm(2,1) = 1 ! Interface
      bcxm(3,1) = 1 ! Interface
c.
      bcxp(1,1) = 3
      bcxp(2,1) = 1 ! Interface
      bcxp(3,1) = 3
c.      
      bcym(1,1) = 1 ! Interface
      bcym(2,1) = 1
      bcym(3,1) = 1 ! Interface
c.      
      bcyp(1,1) = 2
      bcyp(2,1) = 2
      bcyp(3,1) = 2
      
cnn domain sizes = to 1
      xlb(:,:) = 1.0d0
      ylb(:,:) = 1.0d0
      zlb(:,:) = 1.0d0
c
      open(57,file='Airfoil_2D.bin',form='unformatted')

cnn write num of blocks
      write(57) nbX*nbZ

      do zb=1,nbZ
       do xb=1,nbX
        ib=nbX*(zb-1)+xb
        
        if(ib.eq.1)then
        open(11,file='Bl1.dat',form='formatted')
        elseif(ib.eq.2)then
        open(11,file='Bl2.dat',form='formatted')
        else
        open(11,file='Bl3.dat',form='formatted')
        end if
        nzpb(:,:)=1
        read(11,*) nxpb(xb,zb), nypb(xb,zb)

        print *,'Just read', nxpb(xb,zb)
        
        close(11)
        
        if(ib .eq. 2) then
          write(57) nxpb(xb,zb)-1,nypb(xb,zb)-1,nzpb(xb,zb)-1
        else
          write(57) nxpb(xb,zb)-1,nypb(xb,zb)-1,nzpb(xb,zb)-1
        end if

cnn write num of procs
        if(ib.eq.1)then
          write(57) 5,6,1 !2,1,1
        elseif(ib.eq.2)then
          write(57) 14,6,1 !4,1,1 
        else
          write(57) 5,6,1 !2,1,1 
        end if

        write(57) xlb(xb,zb),ylb(xb,zb),zlb(xb,zb)

        xpr=0.0d0; ypr=0.0d0; zpr=0.0d0
        write(57) xpr,ypr,zpr

        write(57) bcxm(xb,zb),bcxp(xb,zb),bcym(xb,zb),
     &            bcyp(xb,zb),bczm(xb,zb),bczp(xb,zb)
cnn write interfaces stuff

cnn INTERF ALL POINTS
cnn
!        write(57)0,0,0
cnn

       if(ib.eq.1) then
        write(57)2,2,3,1,2 
       elseif(ib.eq.2) then
        write(57)2,1,3,1,1 
       else
        write(57)2,2,1,2,2
       end if
!       if(ib.eq.1) then
!        write(57)1,2,1
!       elseif(ib.eq.2) then
!        write(57)2,1,3,1,1
!       else
!        write(57)1,2,2
!       end if
       npt1=8
       npt3=8
       if(ib.eq.1) then
        write(57)xm,npt1+1,nypb(xb,zb),2,1,1,1 
        write(57)ym,1,nxpb(xb,zb),1,1,1,1 
       elseif(ib.eq.2) then
        write(57)xm,2,nypb(xb,zb),npt1+1,1,1,1 
        write(57)xp,2,nypb(xb,zb),npt3+1,1,1,1 
       else
        write(57)xm,npt3+1,nypb(xb,zb),2,1,1,1
        write(57)ym,1,nxpb(xb,zb),1,1,1,1
       end if

cnn end looping blocks for initial info
       end do
      end do
cnn ###################################
c. Open Matlab grid 
        print *,'start loop through blocks and allocate' 
        open(47,file='grid_matlab.bin',access='direct',
     &       form='unformatted',recl=8)
        recnum = 0
      print *,'grid_matlab.bin opened'
c.
      do zb=1,nbZ
       do xb=1,nbX
        ib=nbX*(zb-1)+xb
        print *,'ib=',ib 
        if(ib.eq.1)then
        open(11,file='Bl1.dat',form='formatted')
        allocate(xx1(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
        allocate(yy1(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
        elseif(ib.eq.2)then
        print *,'allocate block 2'
        open(11,file='Bl2.dat',form='formatted')
        allocate(tmp1(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
        allocate(tmp2(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
        allocate(xx2(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
        allocate(yy2(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
        print *,'Allocated: Bl2.dat - xx2',1-halo,':',nxpb(xb,zb)+halo
        print *,'Allocated: Bl2.dat - xx2',1-halo,':',nypb(xb,zb)+halo
        print *,'Allocated: Bl2.dat - yy2',1-halo,':',nxpb(xb,zb)+halo
        print *,'Allocated: Bl2.dat - yy2',1-halo,':',nypb(xb,zb)+halo

        else
        open(11,file='Bl3.dat',form='formatted')
        allocate(xx3(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
        allocate(yy3(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
        end if
     
        print *,'Allocated: Bl3.dat - xx3',1-halo,':',nxpb(xb,zb)+halo 
        print *,'Allocated: Bl3.dat - xx3',1-halo,':',nypb(xb,zb)+halo
 
        nzpb(:,:)=1
        read(11,*) nxpb(xb,zb),nypb(xb,zb)
        
        if(ib.eq.1)then
          do i=1,nxpb(xb,zb)
          do j=1,nypb(xb,zb)
          read(11,*) xx1(i,j), dummy, yy1(i,j)
          end do
          end do
          print *,'xx1 read in'
        elseif(ib.eq.2)then
          do i=1,nxpb(xb,zb)
          do j=1,nypb(xb,zb)
          read(11,*) tmp1(i,j), dummy, tmp2(i,j)
          end do
          end do
          do i=1,nxpb(xb,zb)
          do j=1,nypb(xb,zb)
          xx2(i,j) = tmp1(i,j)
          yy2(i,j) = tmp2(i,j)
          end do
          end do
          print *,'xx2 read in'
!          nxpb(xb,zb) = nxpb(xb,zb)-2
        else
          do i=1,nxpb(xb,zb)
          do j=1,nypb(xb,zb)
          read(11,*) xx3(i,j), dummy, yy3(i,j)
          end do
          end do
          print *,'xx3 read in'
        end if
        
        close(11)
        print *,'start print N',ib,xb,zb
        print *,'N',ib,nxpb(xb,zb),nypb(xb,zb)

       end do
      end do
      
      print *,'xx & yy created'
cnn Add the halos
c. xm - Bl1
!      print *,'test TElow:', xx1(1,npt1),xx2(1,1)
!      print *,'test TElow:', yy1(1,npt1),yy2(1,1)
!
!      print *,'test TEup:' ,xx3(1,npt3),xx2(nxpb(2,1),1)
!      print *,'test TEup:' ,yy3(1,npt3),yy2(nxpb(2,1),1)
       print *,'add halos'
      do i=1-halo,0
        xx1(i,:) = xx1(1,:)
        yy1(i,:) = yy1(1,:)
      end do

      ifi = halo+1
      do i=1-halo,0
        ifi = ifi-1
        xx1(i,npt1:nypb(1,1)+halo) = xx2(ifi,1:nypb(2,1)+halo)
        yy1(i,npt1:nypb(1,1)+halo) = yy2(ifi,1:nypb(2,1)+halo)
      end do
c. xm - Bl3
      do i=1-halo,0
        xx3(i,:) = xx3(1,:)
        yy3(i,:) = yy3(1,:)
      end do

      ifi = halo+1
      do i=1-halo,0
        ifi = ifi-1
        xx3(i,npt3:nypb(3,1)+halo) = 
     &                   xx2(nxpb(2,1)-ifi+1,1:nypb(2,1)+halo)
        yy3(i,npt3:nypb(3,1)+halo) = 
     &                   yy2(nxpb(2,1)-ifi+1,1:nypb(2,1)+halo)
      end do

      print *,'BL3xm done'

c. xp - Bl1
      do i=nxpb(1,1)+1,nxpb(1,1)+halo
        xx1(i,:) = xx1(nxpb(1,1),:)
        yy1(i,:) = yy1(nxpb(1,1),:)
      end do
c. xp - Bl3
      do i=nxpb(3,1)+1,nxpb(3,1)+halo
        xx3(i,:) = xx3(nxpb(3,1),:)
        yy3(i,:) = yy3(nxpb(3,1),:)
      end do

      print *,'BL1&3xp done'

c. ym - Bl1
      ifi = halo+1 !+1+1 for sharp trailing edge case!!!
      do j=1-halo,0
        ifi = ifi-1
        xx1(:,j) = xx3(:,ifi)
        yy1(:,j) = yy3(:,ifi)
      end do
c. ym - Bl3
      ifi = halo+1 !+1+1 for sharp trailing edge case!!!
      do j=1-halo,0
        ifi = ifi-1
        xx3(:,j) = xx1(:,ifi)
        yy3(:,j) = yy1(:,ifi)
      end do

      print *,'BL1&3ym done'

c. yp - Bl1
      do j=nypb(1,1)+1,nypb(1,1)+halo
        xx1(:,j) = xx1(:,nypb(1,1))
        yy1(:,j) = yy1(:,nypb(1,1))
      end do
c. yp - Bl3
      do j=nypb(3,1)+1,nypb(3,1)+halo
        xx3(:,j) = xx3(:,nypb(3,1))
        yy3(:,j) = yy3(:,nypb(3,1))
      end do

      print *,'BL1&3yp done'

c. Block 2
c. xm 
      ifi = halo+1
      do i=1-halo,0
        ifi = ifi-1
        print *, nypb(2,1)+halo-1, '=', nypb(1,1)+halo-npt1
        xx2(i,1-halo:nypb(2,1)+halo) = xx1(ifi,npt1-halo:nypb(1,1)+halo)
        yy2(i,1-halo:nypb(2,1)+halo) = yy1(ifi,npt1-halo:nypb(1,1)+halo)
      end do
c. ym
      do j=1-halo,0
        xx2(:,j) = xx2(:,1)
        yy2(:,j) = yy2(:,1)
      end do
c. yp 
      do j=nypb(2,1)+1,nypb(2,1)+halo
        xx2(:,j) = xx2(:,nypb(2,1))
        yy2(:,j) = yy2(:,nypb(2,1))
      end do
c. xm 
!      ifi = halo+1
!      do i=1-halo,0
!        ifi = ifi-1
!        xx2(i,1:nypb(2,1)+halo) = xx1(ifi,npt1:nypb(1,1)+halo)
!        yy2(i,1:nypb(2,1)+halo) = yy1(ifi,npt1:nypb(1,1)+halo)
!      end do
c. xp 
      ifi = 0
      do i=nxpb(2,1)+1,nxpb(2,1)+halo
        ifi = ifi+1
        print *, nypb(2,1)+halo-1, '=', nypb(3,1)+halo-npt3
        xx2(i,1-halo:nypb(2,1)+halo) = xx3(ifi,npt3-halo:nypb(3,1)+halo)
        yy2(i,1-halo:nypb(2,1)+halo) = yy3(ifi,npt3-halo:nypb(3,1)+halo)
      end do
c. xm 
      ifi = halo+1
      do i=1-halo,0
        ifi = ifi-1
        print *, nypb(2,1)+halo-1, '=', nypb(1,1)+halo-npt1
        xx2(i,1-halo:nypb(2,1)+halo) = xx1(ifi,npt1-halo:nypb(1,1)+halo)
        yy2(i,1-halo:nypb(2,1)+halo) = yy1(ifi,npt1-halo:nypb(1,1)+halo)
      end do
c.
      print *,'test TElow:', xx1(1,npt1),xx2(1,1)
      print *,'test TElow:', yy1(1,npt1),yy2(1,1)

      print *,'test TEup:' ,xx3(1,npt3),xx2(nxpb(2,1),1)
      print *,'test TEup:' ,yy3(1,npt3),yy2(nxpb(2,1),1)
c. Check grid
      do zb=1,nbZ
       do xb=1,nbX
        ib=nbX*(zb-1)+xb
        allocate(xx(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
        allocate(yy(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
c.
        if(ib.eq.1)then
          open(11,file='Bl1h.dat',form='formatted')
          xx = xx1
          yy = yy1
        elseif(ib.eq.2)then
          open(11,file='Bl2h.dat',form='formatted')
          xx = xx2
          yy = yy2
        else
          open(11,file='Bl3h.dat',form='formatted')
          xx = xx3
          yy = yy3
        end if
c.
        write(11,*) nxpb(xb,zb)+2*halo,nypb(xb,zb)+2*halo
c.
        do i=1-halo,nxpb(xb,zb)+halo
         do j=1-halo,nypb(xb,zb)+halo
          write(11,*) xx(i,j),yy(i,j)
         end do
        end do
        close(11)
c.
        deallocate(xx,yy)
        end do
       end do
c.
      do zb=1,nbZ
       do xb=1,nbX
        ib=nbX*(zb-1)+xb
c.
        allocate(xx(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
        allocate(yy(1-halo:nxpb(xb,zb)+halo, 
     &              1-halo:nypb(xb,zb)+halo))
        if(ib.eq.1)then
          xx = xx1
          yy = yy1
        elseif(ib.eq.2)then
          xx = xx2
          yy = yy2
        else
          xx = xx3
          yy = yy3
        end if
        
        write(57) ((xx(i,j),i=1-halo,nxpb(xb,zb)+halo)
     &                     ,j=1-halo,nypb(xb,zb)+halo)
        write(57) ((yy(i,j),i=1-halo,nxpb(xb,zb)+halo)
     &                     ,j=1-halo,nypb(xb,zb)+halo)

c...
c. Matlab grid 
c.
        do k=1,nzpb(xb,zb)
         do j=1-halo,nypb(xb,zb)+halo
          do i=1-halo,nxpb(xb,zb)+halo
           recnum=recnum+1
           write(47,rec=recnum) xx(i,j)
          end do
         end do
        end do

        do k=1,nzpb(xb,zb)
         do j=1-halo,nypb(xb,zb)+halo
          do i=1-halo,nxpb(xb,zb)+halo
           recnum=recnum+1
           write(47,rec=recnum) yy(i,j)
          end do
         end do
        end do

        deallocate(xx,yy)

cnn end looping blocks
        end do
       end do
cnn ###############

      close(57)
      close(47)

c  
 
      stop
c
      end
