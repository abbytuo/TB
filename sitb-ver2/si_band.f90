!      module scalesys
!      implicit none
!      integer:: no=4     !num of atoms;num of orbitals per atom
!      integer:: ncell_x=14,ncell_y=14,ncell_z=14       !num of cells in each direction,even and the cells are centered at O
!      integer:: ncell,na,dimen                    !total num of cells;atoms
!      integer:: npcell,n_pcell=2,dimen_pcell                    !num of primitive cells; num of atoms per primitive cell; degree of freedom per primitive cell
!      double precision,pointer:: x(:,:)          ! coordinates of cells (the left front lowest atom)
!      double precision,pointer:: xsi(:,:)               ! coordinates of Si atoms

!      double precision bondlength,rmax                 !bondlength; radius limit of interactions
!      integer,parameter:: num_k=51      !num of kpoints
!       double precision,parameter:: rcuts=35.d0 !scale of system:radius coord limit to the left lowest atom of each cell
!      end module
       module band_dos
       use scalesys
       implicit none
!--------------band-----------------------
       integer:: n_k=51
       double precision,pointer:: k_rec(:,:)
       integer:: n_kplot=100
       double precision,pointer:: kplot(:),band(:,:) 

!--------------DOS-----------------------
       double precision dos_s,dos_p
       end module
       program tb_diamond
       use scalesys
       use diamond
       use band_dos
       implicit none
       double precision,pointer:: H(:,:,:),S(:,:,:)
       double precision phase,phase_factor(2)

       double complex,pointer:: zr(:,:)  !eigenvalue and eigenvectors of H_pcom
       double precision,pointer:: eigen(:)

       double precision:: H_unit(no,no,2),S_unit(no,no,2)
       double complex:: H_pcom(2*no,2*no),S_pcom(2*no,2*no)
       double complex:: S_pcomh(2*no,2*no)
       double complex,pointer:: H_w(:)
       DOUBLE precision:: H_temp(no),S_temp(no)
       DOUBLE precision:: vector(3),unitv(3)
       DOUBLE precision:: overlap_H,overlap_S

       DOUBLE precision:: es,ep
       DOUBLE precision:: twopai=datan(1.d0)*8.d0

       INTEGER ipiv(2*no)
       INTEGER:: lwork=2*no

       INTEGER i,ii,j,jj,k,l,m,n
       INTEGER i1,j1,ik,jk,kk,info
       INTEGER len_H
       DOUBLE precision rr,factor
       DOUBLE complex,pointer:: work(:)
       DOUBLE precision,pointer:: dwork(:)

       CHARACTER(len=20) filename
       open(11,file='para.in',status='old')
       open(12,file='kpoints.in',status='old')
       open(13,file='dos.in',status='old')

       open(101,file='band.dat')
       open(102,file='overlap.dat')
       open(103,file='output.dat')
!**********************************************************
       read(11,*)
       read(11,*)rmax
!
       close(11)
! 
       read(13,*)
       read(13,*) dos_s,dos_p
       close(13)
       call mkdiamond 
       if(npcell==0) npcell=na/n_pcell
!
       dimen=na*no
       print*,'dimen=',dimen
       allocate(H(dimen,dimen,2))
       allocate(S(dimen,dimen,2))
       do i=1,dimen
         do j=1,dimen
           do ik=1,2
              H(i,j,ik)=CMPLX(0.d0,0.d0)
              S(i,j,ik)=cmplx(0.d0,0.d0)
           enddo
         enddo
       enddo
!
       allocate(k_rec(n_k,3))
       allocate(kplot(n_kplot))
       allocate(band(dimen,n_kplot))
!       
       kplot(1)=0
       read(12,*)(k_rec(1,j),j=1,3)
       do i=2,n_k
          read(12,*)(k_rec(i,j),j=1,3)
          kplot(i)=(k_rec(i,1)-k_rec(i-1,1))**2+&
     & (k_rec(i,2)-k_rec(i-1,2))**2+(k_rec(i,3)-k_rec(i-1,3))**2 
          kplot(i)=dsqrt(kplot(i))+kplot(i-1)
       enddo
       close(12)
!       
       factor=bondlength*4.d0/dsqrt(3.d0)
!      do i=1,n_k
!        do j=1,3
!           k_rec(i,j)=k_rec(i,j)/factor
!        enddo
!      enddo
       dimen_pcell=no*n_pcell
       len_H=dimen_pcell*(dimen_pcell+1)/2
       allocate(H_w(len_H))
       allocate(work(2*dimen-1))
       allocate(dwork(3*dimen-2))
       allocate(eigen(dimen_pcell))
       allocate(zr(dimen_pcell,dimen_pcell))
!********************************************
       print*,'n_kpoints=',n_k
       print*,'npcell=',npcell
       print*,'rmax=',rmax
       do jj=1,dimen_pcell
          write(filename,"(I2)") jj
          open(jj,file='band'//trim(adjustl(filename)))
       enddo
       do k=1,n_k 
!
          call onsite(es,ep)
          do i=1,n_pcell
!            i1=npcell*(i-1)+1
             do j=i,n_pcell
!               j1=npcell*(j-1)
                do ik=1,no
                  do jk=1,no
                    do kk=1,2
                       H_unit(ik,jk,kk)=0.d0
                       S_unit(ik,jk,kk)=0.d0
                    enddo
                  enddo
                enddo
                if(i.ne.j) goto 1
                H_unit(1,1,1)=es
                H_unit(2,2,1)=ep
                H_unit(3,3,1)=ep
                H_unit(4,4,1)=ep
                S_unit(1,1,1)=1.0d0
                S_unit(2,2,1)=1.0d0
                S_unit(3,3,1)=1.0d0
                S_unit(4,4,1)=1.0d0
                write(103,*)'test: es ep= ',es,ep
 1              continue 
                do m=0,npcell-1
                   rr=0.d0
                   do jj=1,3
!                     vector(jj)=xsi(j1+m,jj)-xsi(i1,jj)
                      vector(jj)=xsi(j+n_pcell*m,jj)-xsi(i,jj)
                      rr=rr+vector(jj)*vector(jj)
                   enddo
                   rr=dsqrt(rr)
                   do jj=1,3
                      unitv(jj)=vector(jj)/rr
                   enddo
!cc              
                   if(rr.le.0.1d0) cycle
                   rr=rr*factor
                   if(rr.gt.rmax) cycle
!cc               
                   call integral(rr,H_temp,S_temp)
                   phase=vector(1)*k_rec(k,1)+vector(2)*k_rec(k,2)+&
     & vector(3)*k_rec(k,3)
!实部             
                   phase_factor(1)=dcos(phase*twopai)
!虚部             
                   phase_factor(2)=dsin(phase*twopai)
                   do ik=1,2
!first order interaction correction for H_onsite
                      H_unit(1,1,ik)=H_unit(1,1,ik)+&
     & H_temp(1)*phase_factor(ik)
                      S_unit(1,1,ik)=S_unit(1,1,ik)+&
     & S_temp(1)*phase_factor(ik)
                     
                      H_unit(1,2,ik)=H_unit(1,2,ik)+&
     & H_temp(2)*unitv(1)*phase_factor(ik)
                      S_unit(1,2,ik)=S_unit(1,2,ik)+&
     & S_temp(2)*unitv(1)*phase_factor(ik)
                     
                      H_unit(1,3,ik)=H_unit(1,3,ik)+&
     & H_temp(2)*unitv(2)*phase_factor(ik)
                      S_unit(1,3,ik)=S_unit(1,3,ik)+&
     & S_temp(2)*unitv(2)*phase_factor(ik)
                       
                      H_unit(1,4,ik)=H_unit(1,4,ik)+&
     & H_temp(2)*unitv(3)*phase_factor(ik)
                      S_unit(1,4,ik)=S_unit(1,4,ik)+&
     & S_temp(2)*unitv(3)*phase_factor(ik)
                     
                      overlap_H=unitv(1)*unitv(1)*H_temp(3)+&
     & (1.d0-unitv(1)*unitv(1))*H_temp(4)
                      overlap_S=unitv(1)*unitv(1)*S_temp(3)+&
     & (1.d0-unitv(1)*unitv(1))*S_temp(4)
                      H_unit(2,2,ik)=H_unit(2,2,ik)+&
     & overlap_H*phase_factor(ik)
                      S_unit(2,2,ik)=S_unit(2,2,ik)+&
     & overlap_S*phase_factor(ik)
                     
                      overlap_H=unitv(1)*unitv(2)*(H_temp(3)-H_temp(4))
                      overlap_S=unitv(1)*unitv(2)*(S_temp(3)-S_temp(4))
                      H_unit(2,3,ik)=H_unit(2,3,ik)+&
     & overlap_H*phase_factor(ik)
                      S_unit(2,3,ik)=S_unit(2,3,ik)+&
     & overlap_S*phase_factor(ik)
                     
                      overlap_H=unitv(1)*unitv(3)*(H_temp(3)-H_temp(4))
                      overlap_S=unitv(1)*unitv(3)*(S_temp(3)-S_temp(4))
                      H_unit(2,4,ik)=H_unit(2,4,ik)+&
     & overlap_H*phase_factor(ik)
                      S_unit(2,4,ik)=S_unit(2,4,ik)+&
     & overlap_S*phase_factor(ik)
                     
                      overlap_H=unitv(2)*unitv(2)*H_temp(3)+&
     & (1.d0-unitv(2)*unitv(2))*H_temp(4)
                      overlap_S=unitv(2)*unitv(2)*S_temp(3)+&
     & (1.d0-unitv(2)*unitv(2))*S_temp(4)
                      H_unit(3,3,ik)=H_unit(3,3,ik)+&
     & overlap_H*phase_factor(ik)
                      S_unit(3,3,ik)=S_unit(3,3,ik)+&
     & overlap_S*phase_factor(ik)
                     
                      overlap_H=unitv(2)*unitv(3)*(H_temp(3)-H_temp(4))
                      overlap_S=unitv(2)*unitv(3)*(S_temp(3)-S_temp(4))
                      H_unit(3,4,ik)=H_unit(3,4,ik)+&
     & overlap_H*phase_factor(ik)
                      S_unit(3,4,ik)=S_unit(3,4,ik)+&
     & overlap_S*phase_factor(ik)
                     
                      overlap_H=unitv(3)*unitv(3)*H_temp(3)+&
     & (1.d0-unitv(3)*unitv(3))*H_temp(4)
                      overlap_S=unitv(3)*unitv(3)*S_temp(3)+&
     & (1.d0-unitv(3)*unitv(3))*S_temp(4)
                      H_unit(4,4,ik)=H_unit(4,4,ik)+&
     & overlap_H*phase_factor(ik)
                      S_unit(4,4,ik)=S_unit(4,4,ik)+&
     & overlap_S*phase_factor(ik)
                   enddo
                enddo
                do ii=1,no-1 
                  do jj=ii+1,no
                    do kk=1,2
                       rr=1.d0
                       if(ii.eq.1)rr=-1.d0    !s-orbital
                       H_unit(jj,ii,kk)=H_unit(ii,jj,kk)*rr
                       S_unit(jj,ii,kk)=S_unit(ii,jj,kk)*rr
                    enddo
                  enddo
                enddo
!cc            
                do ii=1,no
                  do jj=1,no
                    do kk=1,2
                       ik=no*(i-1)+ii
                       jk=no*(j-1)+jj 
                       H(ik,jk,kk)=H_unit(ii,jj,kk)
                       S(ik,jk,kk)=S_unit(ii,jj,kk)
                    enddo
                  enddo
                enddo
             enddo
          enddo
          write(103,*)
!        
          do ii=1,dimen_pcell
            do jj=1,dimen_pcell
               H_pcom(ii,jj)=CMPLX(H(ii,jj,1),H(ii,jj,2))
               S_pcom(ii,jj)=CMPLX(S(ii,jj,1),S(ii,jj,2))
            enddo
          enddo
!	  
          do ii=1,dimen_pcell
            do jj=1,dimen_pcell
               S_pcomh(ii,jj)=S_pcom(ii,jj)
            enddo
          enddo
          n=dimen_pcell
          call zgetrf( n, n, S_pcomh, n, ipiv, info)
          call zgetri( n, S_pcomh, n, ipiv, work, lwork, info)
!	  
	  write(103,*) 'H matrix elements'
	  write(103,*) 'real:'
          do i=1,dimen_pcell
	    write(103,10) (real(H_pcom(i,j)),j=1,dimen_pcell)
          enddo
	  write(103,*) 'imag:'
          do i=1,dimen_pcell
	    write(103,10) (aimag(H_pcom(i,j)),j=1,dimen_pcell)
          enddo
	  write(103,*)

          n=dimen_pcell
          call zhpev('N','U',n,H_w,eigen,zr,n,work,dwork,INFO)

	  write(103,*) 'k-points:'
          write(103,*) k, k_rec(k,1),k_rec(k,2),k_rec(k,3)
	  write(103,*)
          do jj=1,dimen_pcell
             band(jj,k)=eigen(jj)
             write(jj,*) kplot(k), eigen(jj)
          enddo
	  write(103,*) 'band'
	  write(103,10) kplot(k),(eigen(jj),jj=1,dimen_pcell)
	  write(103,*)
	  write(101,10) kplot(k),(eigen(jj),jj=1,dimen_pcell)
 10       format(10f10.5)
       enddo
       do jj=1,dimen_pcell
          close(jj)
       enddo
! end loop for k points.
!      call Fermi(eigen,fermi,eban)
       deallocate(H)
       deallocate(S)
       deallocate(k_rec)
       deallocate(kplot)
       deallocate(band)
       deallocate(H_w)
       deallocate(work)
       deallocate(dwork)
       deallocate(eigen)
       deallocate(zr)
       close(101)
       close(102)
       close(103)
       stop
       end
