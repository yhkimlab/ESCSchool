# Makefile for postprocessing (PP)

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../Modules/ \
          $(MOD_FLAG)../../dft-d3/

PPOBJS = \
adduscore.o \
addusdens1d.o \
add_shift_cc.o \
add_shift_lc.o \
add_shift_us.o \
chdens_module.o \
chdens_bspline.o \
chdens_bspline_sphere.o \
compute_ppsi.o \
compute_sigma_avg.o \
cube.o \
dosg.o \
summer_tuto.o \
do_initial_state.o \
do_shift_ew.o \
d_matrix_nc.o \
d_matrix_so.o \
elf.o \
fermisurfer_common.o \
ggen1d.o \
local_dos.o \
local_dos_mag.o \
local_dos1d.o \
openfil_pp.o\
partialdos.o \
paw_postproc.o \
projections_mod.o \
projwfc_box.o \
punch_plot.o \
punch_rism.o \
read_proj.o \
smallgk.o \
solvdens.o \
stm.o \
stop_pp.o \
sum_band_kin.o \
sym_band.o \
vasp_read_chgcar_mod.o \
vasp_xml_module.o \
wannier_proj.o \
wannier_enrg.o \
work_function.o \
write_p_avg.o \
write_proj.o \
write_io_header.o \
write_hamiltonians.o \
globalmod.o \
fouriermod.o \
idwmod.o \
oscdft_pp_mod.o \
oscdft_et_mod.o \
xc_vdW_scale_mod.o \
fft_interpolation_mod.o  

PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a 
QEMODS = $(BASEMODS)

MODULES = $(PWOBJS) $(QEMODS)

TLDEPS= pwlibs

all : tldeps open_grid.x average.x bands.x summer_tuto.x band_interpolation.x dos.x epsilon.x initial_state.x fs.x \
      plan_avg.x plotband.x plotproj.x plotrho.x pmw.x pp.x projwfc.x \
      pawplot.x sumpdos.x pw2wannier90.x pw2critic.x pw2gw.x pw2gt.x \
      wannier_ham.x wannier_plot.x molecularpdos.x \
      pw2bgw.x wfck2r.x fermi_velocity.x fermi_proj.x ppacf.x pprism.x \
      d3hess.x \
      oscdft_pp.x oscdft_et.x


libpp.a : $(PPOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

pp_example.x: pp_example.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ pp_example.o libpp.a $(MODULES) $(QELIBS) 

open_grid.x : open_grid.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		open_grid.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

average.x : average.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		average.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

bands.x : bands.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		bands.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

summer_tuto.x : summer_tuto.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		summer_tuto.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

band_interpolation.x : band_interpolation.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		band_interpolation.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

d3hess.x : d3hess.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		d3hess.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

dos.x : dos.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		dos.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pawplot.x : pawplot.o libpp.a libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		pawplot.o libpp.a $(MODULES) $(QELIBS)

plan_avg.x : plan_avg.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		plan_avg.o libpp.a $(MODULES) \
		$(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

plotband.x : plotband.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		plotband.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

plotproj.x : plotproj.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		plotproj.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

plotrho.x : plotrho.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		plotrho.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pmw.x : poormanwannier.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		poormanwannier.o libpp.a $(MODULES) \
		$(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pp.x : postproc.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		postproc.o libpp.a $(MODULES) \
		$(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

projwfc.x : projwfc.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		projwfc.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw2wannier90.x : pw2wannier90.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		pw2wannier90.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw2critic.x : pw2critic.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		pw2critic.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

wfck2r.x : wfck2r.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		wfck2r.o libpp.a $(MODULES)  $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

initial_state.x : initial_state.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		initial_state.o libpp.a $(MODULES) \
		$(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw2gw.x : pw2gw.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		pw2gw.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

sumpdos.x : sumpdos.o 
	$(LD) $(LDFLAGS) -o $@ sumpdos.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

epsilon.x : epsilon.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ epsilon.o libpp.a $(MODULES) \
	$(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

wannier_ham.x : wannier_ham.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
	wannier_ham.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

wannier_plot.x : wannier_plot.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
	wannier_plot.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

molecularpdos.x : molecularpdos.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		molecularpdos.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw2bgw.x : pw2bgw.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		pw2bgw.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

#bgw2pw.x : bgw2pw.o libpp.a $(MODULES) 
#	$(LD) $(LDFLAGS) -o $@ \
#		bgw2pw.o libpp.a $(MODULES) $(QELIBS)
#	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

fs.x : fermisurface.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		fermisurface.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

fermi_velocity.x : fermi_velocity.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		fermi_velocity.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

fermi_proj.x : fermi_proj.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		fermi_proj.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

# Added by Yang Jiao
ppacf.x : ppacf.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		ppacf.o libpp.a $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw2gt.x : pw2gt.o libpp.a $(MODULES) 
	$(LD) $(LDFLAGS) -o $@ \
		pw2gt.o libpp.a $(MODULES) $(QELIBS)

pprism.x : postrism.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		postrism.o libpp.a $(MODULES) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

# ***OSCDFT BEGIN***
oscdft_pp.x : oscdft_pp.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		oscdft_pp.o libpp.a $(MODULES) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )
oscdft_et.x : oscdft_et.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		oscdft_et.o libpp.a $(MODULES) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )
# ***OSCDFT END***

clean :
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.L libpp.a

include make.depend
