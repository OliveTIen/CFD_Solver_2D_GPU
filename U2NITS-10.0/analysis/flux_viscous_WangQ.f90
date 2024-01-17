!!===========================================
!!  viscous flux computation
!!-------------------------------------------
subroutine flux_viscous
    use mainvar
    implicit none
    integer:: i,j,k,ig,NL,NR,n1,n2
    real*8 :: xy(2),dxyL(2),dxyR(2),vvh(MaxOrder), &
              tgd(2,Nvar),Flu_Vis(Nvar),Flu(Nvar),vvhx(MaxOrder,2),&
              Yl(Nsp-1),Yr(Nsp-1),u_l(Nvar),u_r(Nvar)
    real*8 :: dx,dy,sav,sav1,sav2,sav1n,sav2n,&
              rl,rr,ul,ur,vl,vr,pl,pr,tl,tr,uu,vv,pp,tem,&
              Cpl,gamal,Rcpcvl,Cpr,gamar,Rcpcvr,Cp,gama,Rcpcv,&
              dudx,dudy,dvdx,dvdy,dpdx,dpdy,dtdx,dtdy,dtdn,vmulc,akmu,vlac,div,&
              txx,txy,tyy,btx,bty,dxk
    real*8,parameter:: yita=0.5
    
    ! for each face boundary
    do i=1,  NF
      NL= Neighbor(1,i)
      NR= Neighbor(2,i)
      DX=     VECX(1,i)
      DY=     VECX(2,i)
      n1= FNode(1,i)
      n2= FNode(2,i)
      sav1= dx ! 法向量
      sav2= dy
      sav = max(sqrt(dx*dx+dy*dy),small)
      sav1n=sav1/sav
      sav2n=sav2/sav

      ! for each Gauss Integration Point
      do ig=1, NG
        !!------------------------------
        !! UL and UR
        !!------------------------------
        xy(:)= Coor(:,n1)+ wxx(ig)*(Coor(:,n2)-Coor(:,n1)) ! calculate Gauss Integration Point coordinate
      
        !!! UL
        dxyL(:)= xy(:)- CellXY(:,NL) ! dxyL: n direction
        call fbaseAll(dxyL,vvh,NL) ! vvh: non-dimensional n direction
        ! calculate Conservation variables at cell boundary
        do k=1,Nvar 
          u_l(k)=PA(k,NL)+ sum(gradC(:,k,NL)*vvh(:))
        enddo
      
        pl= u_l(1)
        ul= u_l(2)
        vl= u_l(3)
        tl= u_l(4)
        Yl(:) = u_l(5:Nvar)
        call ComputeGasParameter(Yl, Cpl, gamal, Rcpcvl) ! calculate Rcpcvl by Y, Cp, gamma?
        rl= pl/(Rcpcvl*tl)
        ! if error occurs, set boundary value with cell center value
        if( rl<0. .or. pl<0. .or. tl<0. .or. minval(Yl)<0. .or. maxval(Yl)>1.)then
          pl= PA(1,NL)
          ul= PA(2,NL)
          vl= PA(3,NL)
          tl= PA(4,NL)
          Yl(:)= PA(5:Nvar,NL)
          call ComputeGasParameter(Yl, Cpl, gamal, Rcpcvl)
          rl= pl/(Rcpcvl*tl)
        endif
      
        !!! UR
        !! inner boundary, 
        if( NR>0 )then
          dxyR(:)= xy(:)- (CellXY(:,NR)-SideOfs(:,i))
          call fbaseAll(dxyR,vvh,NR)

          do k=1,Nvar
            u_r(k)=PA(k,NR)+ sum(gradC(:,k,NR)*vvh(:))
          enddo
      
          pr= u_r(1)
          ur= u_r(2)
          vr= u_r(3)
          tr= u_r(4)
          Yr(:)= u_r(5:Nvar)
          call ComputeGasParameter(Yr, Cpr, gamar, Rcpcvr)
          rr= pr/(Rcpcvr*tr)

          if(rr<0. .or. pr<0. .or. tr<0. .or. minval(Yr)<0. .or. maxval(Yr)>1. )then
            pr= PA(1,NR)
            ur= PA(2,NR)
            vr= PA(3,NR)
            tr= PA(4,NR)
            Yr(:)= PA(5:Nvar,NR)
            call ComputeGasParameter(Yr, Cpr, gamar, Rcpcvr)
            rr= pr/(Rcpcvr*tr)
          endif

        !! wall boundary
        elseif( NR==WallID )then  
          rr= rl; pr= pl; tr= tl
          ul= 0.; vl= 0.
          ur= 0.; vr= 0.
          Yr= Yl; Cpr= Cpl; Rcpcvr= Rcpcvl; gamar= gamal

        !! other boundary
        else
          rr= rl; pr= pl; tr= tl
          ur= ul; vr= vl
          Yr= Yl; Cpr= Cpl; Rcpcvr= Rcpcvl; gamar= gamal
        endif
      
        !!------------------------------
        !! viscous flux
        !!------------------------------
        tem= 0.5*(tl+tr ) ! temperature
        uu=  0.5*(ul+ur )
        vv=  0.5*(vl+vr )
        pp=  0.5*(pl+pr )
        call ComputeGasParameter(0.5*(Yl+Yr),Cp,gama,Rcpcv)
        rr= pp/(Rcpcv*tem)  ! T = p/((gamma-1)*rho) -> rho = p/((gamma-1)*T)
      
        !! 对内部单元，tgd取左右平均值加yita*sav1n*(u_r(k)-u_l(k))/dxk
        !! 对边界单元，tgd取左值
        ! tgd - gradient matrix
        ! NL - left cell index
        call fbaseGrad_dxdy( dxyL,vvhx,NL )
        ! 上述代码等价于：
        !  dx=1./reflen(1,NL)
        !  dy=1./reflen(2,NL)
        !  
        !  vvhx(:,1)=(/ 1.,0./)*dx
        !  vvhx(:,2)=(/ 0.,1./)*dy

        do k =1, Nvar
          ! GradC(MaxOrder,NVar,NCell)
          ! gradC(1:MaxOrder,k,i)= matmul(AMLS(1:MaxOrder,1:KFace(i),i),C(1:KFace(i),k)), in file "reconstruction.f90"
          tgd(1,k)= sum( gradC(:,k,NL)*vvhx(:,1) )
          tgd(2,k)= sum( gradC(:,k,NL)*vvhx(:,2) )
        enddo
      
        !! NR
        if(NR>0)then
          dxk=min(vol(NL),vol(NR))/sav
          call fbaseGrad_dxdy( dxyR,vvhx,NR )
          !! gradient obtained by dGRP formula
          do k =1, Nvar
            tgd(1,k)= 0.5*(tgd(1,k)+ sum( gradC(:,k,NR)*vvhx(:,1) ) ) + &
                      yita*sav1n*(u_r(k)-u_l(k))/dxk
            tgd(2,k)= 0.5*(tgd(2,k)+ sum( gradC(:,k,NR)*vvhx(:,2) ) ) + &
                      yita*sav2n*(u_r(k)-u_l(k))/dxk
          enddo
        endif
        
        dpdx= tgd(1,1)
        dpdy= tgd(2,1)

        dudx= tgd(1,2)
        dudy= tgd(2,2)
        
        dvdx= tgd(1,3)
        dvdy= tgd(2,3)
        
        dtdx= tgd(1,4)
        dtdy= tgd(2,4)
        
        if( NR==WallID )then
          dtdn= dtdx*sav1n+ dtdy*sav2n
          dtdx= dtdx- dtdn*sav1n
          dtdy= dtdy- dtdn*sav2n
        endif
        
        if( IfConstantViscousity )then
          vmulc= umu_give !! 1.0
        else
          !C_T  = 110.4/tem_ref
          !abs(tem)**1.5 *(1.+C_T)/(tem+C_T)
          vmulc=1.458*abs(tem)**1.5/(tem+110.4)*1.0d-6
        endif
        
        vlac =-2./3.*vmulc
        div= dudx+dvdy
        txx= 2.*vmulc*dudx+vlac*div ! tau_xx = 2*mu*u_x + lambda*(u_x+v_y)
        tyy= 2.*vmulc*dvdy+vlac*div
        txy= vmulc*(dudy+dvdx)      ! tau_xy = mu*(u_y+v_x)
        
        akmu= cp*vmulc/prl  !! gamma*mu/p; vmulc: mu
        
        btx=akmu*dtdx+ uu*txx+ vv*txy
        bty=akmu*dtdy+ uu*txy+ vv*tyy
      
        Flu_Vis(1)= 0.
        Flu_Vis(2)=-(sav1*txx+ sav2*txy)
        Flu_Vis(3)=-(sav1*txy+ sav2*tyy)
        Flu_Vis(4)=-(sav1*btx+ sav2*bty)

        ! diffusion terms needed to be added
        Flu_Vis(4)= Flu_Vis(4) - rr*Cp*gas_ck*(sav1*dtdx + sav2*dtdy)
        Flu_Vis(5:Nvar)= -rr*gas_cd*(sav1*tgd(1,5:Nvar)+sav2*tgd(2,5:Nvar))  

        !!-------------------------------------------------

        Flu(:)= Flu_Vis(:)     !! /Re_in: for non-dimensionalize

        do k=1, Nvar
        hWA(k,NL)= hWA(k,NL) - Flu(k)*wgg(ig) ! numeral flux in element
        if( NR>0 .and. Fproperty(i)/=100)&
        hWA(k,NR)= hWA(k,NR) + Flu(k)*wgg(ig)
        enddo

      enddo
    
    enddo
      
      return
end subroutine

