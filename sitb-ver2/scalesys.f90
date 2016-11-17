       module scalesys
       implicit none
       integer,parameter:: no=4     !num of atoms;num of orbitals per atom
       integer:: ncell_x=9,ncell_y=9,ncell_z=9       !num of cells in each direction,even and the cells are centered at O
       integer:: ncell,na,dimen                    !total num of conventional cells;atoms;total dimension
       integer:: npcell,dimen_pcell                    !num of primitive cells; num of atoms per primitive cell; degree of freedom per primitive cell
       double precision bondlength,rmax                 !bondlength; radius limit of interactions
       integer,parameter:: num_k=51      !num of kpoints
!       double precision,parameter:: rcuts=35.d0 !scale of system:radius coord limit to the left lowest atom of each cell
       end module

        module diamond
        use scalesys
        implicit none

        integer:: n_cell=8,n_pcell=2  ! num of atoms per conventional cell;primitive cell
        
        double precision,pointer:: xsi(:,:)               ! coordinates of Si atoms
          
        contains
        subroutine Si_crystal_primitive_cell
        implicit none
          integer i,j,k,m,n
	  double precision rr
          double precision,pointer:: x(:,:)          ! coordinates of cells (the left front lowest atom)
          double precision,pointer:: xb(:,:)      ! coordinates of basis axis corresponding to general(original)
          double precision,pointer:: xbpos(:,:),xpos(:,:)               ! coordinates of Si atoms
	  common /pos_/ xbpos,xpos

          allocate(x(ncell_x*ncell_y*ncell_z,3))
          allocate(xsi(ncell_x*ncell_y*ncell_z*n_pcell,3))
	  allocate(xb(3,3))

          open(10,file='lattice_primitive.in')
	  read(10,*) bondlength
	  read(10,*) n_pcell
	  do i=1,3
	     read(10,*) (xb(i,j),j=1,3)
	  enddo
	  do i=1,n_pcell
	     read(10,*) (xsi(i,j),j=1,3)
	  enddo
	  close(10)
          rr=bondlength*4.d0/dsqrt(3.d0)
!         do i=1,3
!            xb(i)=xb(i)*rr
!            yb(i)=yb(i)*rr
!            zb(i)=zb(i)*rr
!         enddo
!         do i=1,n_pcell
!           do j=1,3
!              xsi(i,j)=xsi(i,j)*rr
!           enddo
!         enddo
          m=1
          do i=-ncell_x/2,(ncell_x-1)-(ncell_x-1)/2
            do j=-ncell_y/2,(ncell_y-1)-(ncell_y-1)/2
              do k=-ncell_z/2,(ncell_z-1)-(ncell_z-1)/2
		  if(i==0.and.j==0.and.k==0) cycle
                  m=m+1
		  do n=1,3
                     x(m,n)=i*xb(1,n)+j*xb(2,n)+k*xb(3,n)
		  enddo
              enddo
            enddo
          enddo
          npcell=m

	  do i=1,3
            do j=1,3
	       xb(i,j)=xb(i,j)*ncell_x
            enddo
	  enddo
          do i=1,npcell
            do k=1,n_pcell
	      do j=1,3
                m=(i-1)*n_pcell+k
                xsi(m,j)=x(i,j)+xsi(k,j)
              enddo
            enddo
          enddo
          na=m
       print*,'npcell=',npcell
       print*,'n_pcell=',n_pcell
       print*,'na=',na
          allocate(xbpos(3,3))
          allocate(xpos(ncell_x*ncell_y*ncell_z*n_pcell,3))
	  do i=1,3
            do j=1,3
               xbpos(i,j)=xb(i,j)*rr
            enddo
	  enddo
	  do i=1,na
            do j=1,3
               xpos(i,j)=xsi(i,j)*rr
            enddo
	  enddo
          deallocate(x)
	  deallocate(xb)
        return
        end subroutine

        subroutine Si_crystal_conventional_cell
        implicit none
          double precision:: ax=5.431d0,ay=5.431d0,az=5.431d0 ! lattice const
          double precision,pointer:: x(:,:)          ! coordinates of cells (the left front lowest atom)
          double precision,pointer:: xb(:),yb(:),zb(:)      ! coordinates of basis axis corresponding to general(original)
          integer i,j,k,m,n
          allocate(x(ncell_x*ncell_y*ncell_z,3))
          allocate(xsi(ncell_x*ncell_y*ncell_z*n_cell,3))
	  allocate(xb(3))
	  allocate(yb(3))
	  allocate(zb(3))
          
          xb(1)=ax*ncell_x
          xb(2)=0.d0
          xb(3)=0.d0
          yb(1)=0.d0
          yb(2)=ay*ncell_y
          yb(3)=0.d0
          zb(1)=0.d0
          zb(2)=0.d0
          zb(3)=az*ncell_z
          x(1,1)=0.d0
          x(1,2)=0.d0
          x(1,3)=0.d0

          do i=1,ncell_x
            do j=1,ncell_y
              do k=1,ncell_z
                  m=k+(j-1)*ncell_z+(i-1)*ncell_y*ncell_z
                  x(m,1)=(i-1)*ax
                  x(m,2)=(j-1)*ay
                  x(m,3)=(k-1)*az
              enddo
            enddo
          enddo
          ncell=ncell_x*ncell_y*ncell_z

          xsi(1,1)=0.d0
          xsi(1,2)=0.d0
          xsi(1,3)=0.d0
          xsi(2,1)=0.d0
          xsi(2,2)=0.5d0
          xsi(2,3)=0.5d0
          xsi(3,1)=0.5d0
          xsi(3,2)=0.d0
          xsi(3,3)=0.5d0
          xsi(4,1)=0.5d0
          xsi(4,2)=0.5d0
          xsi(4,3)=0.d0

          do m=1,4
            do k=1,3  
              xsi(m+4,k)=xsi(m,k)+0.25d0
            enddo
          enddo
          do m=1,8
            do k=1,3  
              xsi(m,k)=xsi(m,k)*ax
            enddo
          enddo

          do i=1,ncell_x*ncell_y*ncell_z
            do k=1,n_cell
	      do j=1,3
                m=(i-1)*n_cell+k
                xsi(m,j)=x(i,j)+xsi(k,j)
              enddo
            enddo
          enddo
          na=ncell_x*ncell_y*ncell_z*n_cell
        return
        end subroutine
        
        end module


        subroutine mkdiamond
!       program main
        use scalesys
        use diamond
        implicit none
        double precision,pointer:: xbpos(:,:),xpos(:,:)               ! coordinates of Si atoms
	common /pos_/ xbpos,xpos
        integer i,j,k,m

          open(1,file='pos.vasp')
          write(1,*) 'Si'
          write(1,*) 1
  
!         call Si_crystal_conventional_cell
          call Si_crystal_primitive_cell
          write(1,*) (xbpos(1,j),j=1,3)
          write(1,*) (xbpos(2,j),j=1,3)
          write(1,*) (xbpos(3,j),j=1,3)
          write(1,*) 'Si'
          write(1,*) na
          write(1,*) 'C'

          do m=1,na
              write(1,*) xpos(m,1),xpos(m,2),xpos(m,3)
          enddo
          close(1)
          deallocate(xbpos)
	  deallocate(xpos)
        return
        end 
