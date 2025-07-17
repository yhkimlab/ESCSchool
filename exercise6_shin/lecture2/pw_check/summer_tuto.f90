!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM summer_tuto
  !-----------------------------------------------------------------------
  !
  ! 
  !
  USE io_files,  ONLY : prefix, tmp_dir
  USE mp_global, ONLY : mp_startup
  USE mp_pools,  ONLY : npool
  USE control_flags, ONLY : gamma_only
  USE environment,   ONLY : environment_start, environment_end
  USE wvfct,     ONLY : nbnd
  USE klist,     ONLY : nkstot, two_fermi_energies
  USE noncollin_module, ONLY : noncolin, i_cons
  USE lsda_mod,  ONLY : nspin
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE mp_images, ONLY : intra_image_comm
  USE mp,                ONLY : mp_bcast, mp_sum, mp_barrier
  USE parameters,ONLY : npk
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER (len=256) :: filband, filp, outdir
  LOGICAL :: lsigma(4), lsym, lp, no_overlap, plot_2d
  INTEGER :: spin_component, firstk, lastk
  INTEGER :: ios
  !
  NAMELIST / summ / outdir, prefix 
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'summ' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'

  ios = 0
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, summ, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
     !
  ENDIF
  !
 write(stdout,*) ios, tmp_dir
  !
  CALL mp_bcast( ios, ionode_id, intra_image_comm )
  IF (ios /= 0) CALL errore ('summ', 'reading summ namelist', abs(ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, intra_image_comm )
  CALL mp_bcast( prefix, ionode_id, intra_image_comm )

  !
  CALL read_file()
  !
  !
  CALL openfil_pp()
  !
  write(stdout,*)' check 1'
  CALL check_norm()

  CALL check_kin()

  CALL check_rho_Hatree()
  !
  CALL environment_end ( 'summ' )
  !
  CALL stop_pp
  STOP
END PROGRAM summer_tuto
!
!-----------------------------------------------------------------------
SUBROUTINE check_norm ()
  !-----------------------------------------------------------------------
  !
  !
  USE kinds,                ONLY : dp
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : at
  USE constants,            ONLY : rytoev
  USE gvect,                ONLY : g, ngm
  USE klist,                ONLY : xk, nks, nkstot, ngk, igk_k
  USE io_files,             ONLY : iunpun, nwordwfc, iunwfc
  USE wvfct,                ONLY : nbnd, et, npwx
  USE uspp,                 ONLY : nkb, vkb
  USE uspp_param,           ONLY : upf, nh, nhm
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions, ONLY : evc
  USE io_global,            ONLY : ionode, ionode_id, stdout
 USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,                ONLY : mp_bcast, mp_sum, mp_barrier
  USE mp,                   ONLY : mp_bcast
  USE mp_images,            ONLY : intra_image_comm
  USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                   deallocate_bec_type, becp
  USE uspp_init,            ONLY : init_us_2

  IMPLICIT NONE
  ! becp   : <psi|beta> at current  k-point
  INTEGER :: ibnd, jbnd, iter, i, ik, ig,jpw,  npw
  LOGICAL :: mask(nbnd,nbnd)
  REAL(DP)  ::rnfac(nbnd,nks)

  
  !
  DO ik = 1, nks
     !
     !   read eigenfunctions
     !
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
     !
     ! calculate becp = <psi|beta>, needed to compute spsi = S|psi>
     !
     npw = ngk(ik)

      rnfac(:,:) = 0.0d0
       do ibnd = 1 , nbnd
       do jbnd = 1 , nbnd
         do jpw = 1, npw
         rnfac(jbnd,ibnd) = rnfac(jbnd,ibnd) +dble(evc(jpw,ibnd)*CONJG(evc(jpw,jbnd)))
         enddo
         if(noncolin) then
         do jpw = npwx+1, npwx+npw
         rnfac(jbnd,ibnd) = rnfac(jbnd,ibnd) +dble(evc(jpw,ibnd)*CONJG(evc(jpw,jbnd)))
         enddo
         endif
      enddo
    enddo
         call mp_sum(rnfac,intra_bgrp_comm)

         DO ibnd = 1,nbnd
        DO jbnd = 1,nbnd
          write(stdout,*) "j,i" , ibnd,jbnd 
          write(stdout,*) "<j|i>",rnfac(jbnd,ibnd) 
        ENDDO
        ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE check_norm
!-----------------------------------------------------------------------
SUBROUTINE check_kin ()
  !-----------------------------------------------------------------------
  !
  !
  USE kinds,                ONLY : dp
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : at
  USE constants,            ONLY : rytoev
  USE gvect,                ONLY : g, ngm
  USE klist,                ONLY : xk, nks, nkstot, ngk, igk_k
  USE io_files,             ONLY : iunpun, nwordwfc, iunwfc
  USE wvfct,                ONLY : nbnd, et, npwx,g2kin,wg
  USE uspp,                 ONLY : nkb, vkb
  USE uspp_param,           ONLY : upf, nh, nhm
  USE noncollin_module,     ONLY : noncolin, npol
  USE mp_pools,             ONLY : inter_pool_comm
  USE wavefunctions, ONLY : evc
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
 USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,                ONLY : mp_bcast, mp_sum, mp_barrier
  USE mp_images,            ONLY : intra_image_comm
  USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                   deallocate_bec_type, becp
  USE uspp_init,            ONLY : init_us_2

  IMPLICIT NONE
  COMPLEX(DP) ::hpsi(npwx*npol,nbnd)
  ! becp   : <psi|beta> at current  k-point
  INTEGER :: ibnd, jbnd, iter, i, ik, ig,jpw,  npw,ipw
  REAL(DP)  :: Energy,h_tmp

 Energy =0.0d0
  !
  DO ik = 1, nks
     !
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
     !
     ! calculate becp = <psi|beta>, needed to compute spsi = S|psi>
     !
     npw = ngk(ik)

     call  g2_kin(ik)

     hpsi = DCMPLX(0.0d0,0.0d0)
     DO ibnd = 1, nbnd

        hpsi (1:npw, ibnd) = hpsi(1:npw, ibnd) + g2kin (1:npw) * evc (1:npw, ibnd)
        IF ( noncolin ) THEN
           hpsi (npwx+1:npwx+npw, ibnd) = hpsi(npwx+1:npwx+npw,ibnd) + g2kin (1:npw) * evc (npwx+1:npwx+npw, ibnd)
        END IF
     END DO


     do ibnd = 1, nbnd
        h_tmp = 0.0d0

          do ipw = 1, npw
            h_tmp = h_tmp + dble(conjg(evc(ipw,ibnd))*hpsi(ipw,ibnd))
          enddo
        
          if(noncolin) then
            do ipw = npwx+1, npwx+npw
             h_tmp = h_tmp + dble(conjg(evc(ipw,ibnd))*hpsi(ipw,ibnd))
            enddo
          endif
             call mp_sum(h_tmp,intra_bgrp_comm)
 
         Energy = Energy + h_tmp*wg(ibnd,ik) 
      enddo

    enddo

    call mp_sum(Energy,inter_pool_comm)

     write(stdout,*) "Kinetic energy",Energy 


  !
  RETURN
  !
END SUBROUTINE check_kin
!-----------------------------------------------------------------------
SUBROUTINE check_rho_Hatree ()
  !-----------------------------------------------------------------------
  !
  !
  USE kinds,                ONLY : dp
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : at
  USE constants,            ONLY : rytoev
  USE gvect,                ONLY : g, ngm, gstart,gg
  USE klist,                ONLY : xk, nks, nkstot, ngk, igk_k
  USE io_files,             ONLY : iunpun, nwordwfc, iunwfc
  USE wvfct,                ONLY : nbnd, et, npwx,g2kin,wg
  USE uspp,                 ONLY : nkb, vkb
  USE uspp_param,           ONLY : upf, nh, nhm
  USE noncollin_module,     ONLY : noncolin, npol
  USE mp_pools,             ONLY : inter_pool_comm
  USE wavefunctions, ONLY : evc
    USE constants,         ONLY : fpi, e2
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
 USE mp_bands,  ONLY : intra_bgrp_comm
  USE fft_interfaces,    ONLY : fft_interpolate
 USE mp,                ONLY : mp_bcast, mp_sum, mp_barrier
  USE mp_images,            ONLY : intra_image_comm
  USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                   deallocate_bec_type, becp
  USE uspp_init,            ONLY : init_us_2
    USE fft_base,         ONLY : dffts,dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  USE lsda_mod,  ONLY : nspin
  USE wavefunctions, ONLY : psic
    USE cell_base, ONLY : omega,tpiba2
   USE lsda_mod, ONLY : lsda, nspin, isk

  IMPLICIT NONE
  COMPLEX(DP) ::hpsi(npwx*npol,nbnd)
  ! becp   : <psi|beta> at current  k-point
  INTEGER :: ibnd, jbnd, iter, i, ik, ig,jpw,  npw,ipw, is
  REAL(DP)  :: Energy,h_tmp

  COMPLEX(DP) :: sdpsi(npwx*npol,nbnd,nks)
  REAL(DP)     :: sdrho_r(dfftp%nnr,nspin),charge,ehart
  COMPLEX(DP), Allocatable  ::cwfcs(:,:,:),cwfcs2(:,:,:)
  REAL(kind=DP) :: grids(dffts%nnr),gridd(dfftp%nnr)
  COMPLEX(DP)   :: sdrho_g(ngm,nspin)
  REAL(DP)              :: fac
  REAL(DP)              :: rgtot_re, rgtot_im



  allocate(cwfcs(dffts%nnr,nbnd,nks))
        grids = 0
        gridd = 0
        sdrho_r(:,:) = 0.d0
        sdrho_g(:,:) = 0.d0
        is = 1
        cwfcs = DCMPLX(0.0d0,0.0d0)

        do ik = 1,nks
         npw = ngk(ik)
         IF ( lsda ) is = isk(ik)
          do i = 1 , nbnd
           if(wg(i,ik) ==0.0d0)   go to 39
             CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
             psic(1:dffts%nnr)=0.d0
             psic(dffts%nl(igk_k(1:npw,ik)))  = evc(1:npw,i)
             CALL invfft ('Wave', psic, dffts)
             cwfcs(1:dffts%nnr,i,ik)= psic(1:dffts%nnr)
39           continue 
          enddo


          grids = 0
          gridd = 0
          do i = 1, nbnd
!            do i = 1 ,nbnd
             if(wg(i,ik) ==0.0d0)   go to 49
             grids(1:dffts%nnr)=grids(1:dffts%nnr) +DBLE(cwfcs(1:dffts%nnr,i,ik) & 
                      *CONJG(cwfcs(1:dffts%nnr,i,ik)))*wg(i,ik)
49           continue 
          enddo
             call fft_interpolate(dffts,grids,dfftp,gridd)     
             sdrho_r(1:dfftp%nnr,1) =sdrho_r(1:dfftp%nnr,1)+gridd(1:dfftp%nnr)
           if(nspin==2) then
             if(is==1) sdrho_r(1:dfftp%nnr,2) =sdrho_r(1:dfftp%nnr,2)+gridd(1:dfftp%nnr)
             if(is==2) sdrho_r(1:dfftp%nnr,2) =sdrho_r(1:dfftp%nnr,2)-gridd(1:dfftp%nnr)
           endif
         enddo
     

     call mp_sum(sdrho_r, inter_pool_comm) 

         sdrho_r(:,:) = sdrho_r(:,:) /omega
    
     deallocate(cwfcs)


     do is = 1, nspin
         psic(:)=(0.d0,0.d0)
         psic(1:dfftp%nnr)=dcmplx(sdrho_r(1:dfftp%nnr,is),0.d0)
         CALL fwfft ('Rho', psic, dfftp)
         sdrho_g(1:ngm,is)=psic(dfftp%nl(1:ngm))
     enddo 

  charge = 0.D0
  !
!!!! Take G=0 value
  IF ( gstart == 2 ) THEN
     !
     charge = omega*REAL( sdrho_g(1,1) )
     !
  ENDIF
  !
  CALL mp_sum( charge, intra_bgrp_comm )


     write(stdout,*) "Charge",charge 

      ehart = 0.0d0

       DO ig = gstart, ngm
           !
           fac = 1.D0 / gg(ig) 
           !
           rgtot_re = REAL(  sdrho_g(ig,1) )
           rgtot_im = AIMAG( sdrho_g(ig,1) )
           !
           ehart = ehart + ( rgtot_re**2 + rgtot_im**2 ) * fac
           !
       ENDDO
!$omp end parallel do
     fac = e2 * fpi / tpiba2
     ehart = ehart * fac
     ehart = ehart * 0.5D0 * omega

  CALL mp_sum( ehart, intra_bgrp_comm )
     write(stdout,*) "Hartree",ehart
  RETURN
  !
END SUBROUTINE check_rho_Hatree
