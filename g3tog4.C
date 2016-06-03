void g3tog4(
        int nevents = 1000, 
        std::string inputname = "oscar.input"
        ){

    //gSystem->Load("libjprof.so");
    //prof *Pr = new prof;

    gROOT->LoadMacro("g4libs.C");
    g4libs();

    gSystem->Load("libfun4all.so");	// framework + reco modules
    gSystem->Load("libfvtx_subsysreco.so");
    gSystem->Load("libsimreco.so");
    gSystem->Load("libPHHepMCNode.so");
    //gSystem->Load("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/install/lib/libPHG3toG4.so");
    gSystem->Load("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/install/lib/libPHG3toG4.so");
    //gSystem->Load("libPHG3toG4.so");

    //gSystem->Load("../PHHepMCNode/install/lib/libPHHepMC.so");
    //gSystem->Load("libPHHepMC.so");
    //gSystem->Load("libPHHepMCNode.so");

    Fun4AllServer *se = Fun4AllServer::instance();
    se->Verbosity(0);

    SubsysReco *sync = new SyncSimreco();
    se->registerSubsystem(sync);

    PHG3toG4 *myG3toG4 = new PHG3toG4("PHG3toG4","Sim3D++.root","geom_run15_v1.root");
    //PHG3toG4 *myG3toG4 = new PHG3toG4("PHG3toG4","Sim3D++.root","geom_new_nooverlaps.root");
    myG3toG4->Verbosity(3);
    //myG3toG4->SetBField(false); //zero-field sim
    //myG3toG4->SetPhysicsList("QGSP_BERT");
    //myG3toG4->SetPhysicsList("FTFP_BERT");

    PHG3toG4ActivityTrigger *trig = new PHG3toG4ActivityTrigger();
    //trig->AddParticles("9999"); //Inclusive mode
    //trig->AddParticles("13,321,211,2212");
    //trig->AddParticles("13");
    trig->AddDetector("MUI");
    trig->AddDetectorMinHits("MUI",1);
    //trig->AddDetector("MUT");
    //trig->AddDetectorMinHits("MUT",1);
    //trig->SetDoAbsParticlesOnly();
    trig->PrintConfig();
    myG3toG4->RegisterTrigger(trig);


    //myG3toG4->UseParticleGun("mu+",5/*GeV*/);
    //myG3toG4->UseParticleGunForwardOnly();
    //myG3toG4->UseParticleGunRandomVtx(0,0,0,0,-10,10);
    myG3toG4->ReportEvery(1);

    se->registerSubsystem(myG3toG4);

    // DST output manager
    //Fun4AllDstInputManager *dst_input_mgr  = new Fun4AllDstInputManager("HepMC","DST");
    //se->registerInputManager(dst_input_mgr);
    //se->fileopen( dst_input_mgr->Name(), inputname.c_str() );
    //Fun4AllHepMCInputManager *in = new Fun4AllHepMCInputManager("DSTIN");
    //in->ReadOscar(1);
    //in->UseWriteableNode();
    //se->registerInputManager(in);
    //in->fileopen(inputname);

    Fun4AllOscarInputManager *in = new Fun4AllOscarInputManager("OSCARIN");
    in->Verbosity(1);
    se->registerInputManager(in);
    in->fileopen(inputname.c_str());

    se->run(nevents);

    se->End();
}

