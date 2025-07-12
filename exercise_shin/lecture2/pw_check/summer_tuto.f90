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
  ! See files INPUT_BANDS.* in Doc/ directory for usage
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
  NAMELIST / bands / outdir, prefix, filband, filp, spin_component, lsigma,&
                       lsym, lp, filp, firstk, lastk, no_overlap, plot_2d
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
     READ (5, bands, iostat = ios)
     !
     lsigma(4)=.false.
     tmp_dir = trimcheck (outdir)
     !
  ENDIF
  !
  !
  CALL mp_bcast( ios, ionode_id, intra_image_comm )
  IF (ios /= 0) CALL errore ('bands', 'reading bands namelist', abs(ios) )
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
  CALL summ_band()
  !
  CALL environment_end ( 'summ' )
  !
  CALL stop_pp
  STOP
END PROGRAM summer_tuto
!
!-----------------------------------------------------------------------
SUBROUTINE summ_band (filband, spin_component, lsigma, no_overlap)
  !-----------------------------------------------------------------------
  !
  !    This routine writes the band energies on a file. The routine orders
  !    the eigenvalues using the overlap of the eigenvectors to give
  !    an estimate crossing and anticrossing of the bands. This simplified
  !    method works in many, but not in all the cases.
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
  USE mp,                   ONLY : mp_bcast
  USE mp_images,            ONLY : intra_image_comm
  USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                   deallocate_bec_type, becp
  USE uspp_init,            ONLY : init_us_2

  IMPLICIT NONE
  CHARACTER (len=*) :: filband
  INTEGER, INTENT(IN) :: spin_component
  LOGICAL, INTENT(IN) :: lsigma(4), no_overlap

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
END SUBROUTINE summ_band
