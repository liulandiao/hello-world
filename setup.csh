echo "#########################"
echo "####RUN15PP200GEV LIB####"
echo "##G4 SIMULATION VERSION##"
echo "#########################"

unsetenv OFFLINE_MAIN
unsetenv ONLINE_MAIN
unsetenv LD_LIBRARY_PATH
unsetenv ROOTSYS
unsetenv G4_MAIN

source /opt/phenix/bin/phenix_setup.csh new
#source /opt/phenix/bin/phenix_setup.csh ana.500
#source /opt/phenix/bin/phenix_setup.csh pro.105

cd $G4_MAIN/bin
source geant4.csh
cd -
#source ${G4_MAIN}/bin/geant4.csh

setenv G4SYSTEM Geant4-10.0.2

#VMC
setenv VGM_INSTALL /gpfs/mnt/gpfs02/phenix/hhj/hhj1/snowball/local/vgm-4.01/vgm-build
setenv USE_VGM 1
setenv G4INSTALL ${G4_MAIN}

#HEPMC
setenv HEPMC "/opt/phenix/stow/HepMC"
setenv PYTHIA8 "/afs/rhic.bnl.gov/x8664_sl6/opt/phenix/stow/pythia8180"
setenv FFTW3 "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/fftw-3.3.4"
setenv LD_LIBRARY_PATH ${G4_MAIN}/lib:${LD_LIBRARY_PATH}
setenv LD_LIBRARY_PATH /gpfs/mnt/gpfs02/phenix/hhj/hhj1/snowball/local/geant4_vmc/lib/tgt_linux:${LD_LIBRARY_PATH}
setenv LD_LIBRARY_PATH /gpfs/mnt/gpfs02/phenix/hhj/hhj1/snowball/local/vgm-4.01/vgm-build/lib64:${LD_LIBRARY_PATH}

setenv MY_INSTALL "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/install" 
setenv LD_LIBRARY_PATH ${MY_INSTALL}/lib:${LD_LIBRARY_PATH}
setenv PATH ${MY_INSTALL}/bin:${PATH}

echo "##LD_LIBRARY_PATH##" 
echo ${LD_LIBRARY_PATH}
