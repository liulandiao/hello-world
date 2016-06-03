// $Id: Fun4FVTX_RecoDST_sim.C,v 1.7 2014/10/01 16:39:32 slash Exp $
/*!
  mutoo dst analysis loop for mutoo slowsim simulated dst, possibly embeded
  into either simulated or real data dst
  */
void Fun4FVTX_RecoDST_sim(
        int nEvents = 0,
        char *signalfile = "dst_sim.root",
        char *backgroundfile = "background.root", // 0,
        //char *backgroundfile = NULL, // 0,
        char *ana_file = "muon_ana_ntuples.root",
        char *dstfile = "dstembed_out.root",
        char *prdffile = "data_out.prdf",
        char *ndstfile = "ndstembed_out.root",
        const char *singlepdstout = "singlemuon_embed_pdst.root",
        const char *dimuonpdstout = "dimuon_embed_pdst.root",
        float rCutFactor = 1.0,
        int run_number = 434153 // Run15pAu200 
        //int run_number = 437107 // Run15pAl200
        //int run_number = 421716 // Run15pp200 
        //int run_number = 408921 // Run14AuAu200
        )
{

    // dump configuration
    cout << endl;
    cout << "Fun4Muons_RecoDST_sim" << endl;
    if( signalfile ) cout << "Fun4Muons_RecoDST_sim - signalfile : " << signalfile << endl;
    if( backgroundfile ) cout << "Fun4Muons_RecoDST_sim - backgroundfile : " << backgroundfile << endl;
    if( dstfile ) cout << "Fun4Muons_RecoDST_sim - dstfile : " << dstfile << endl;
    if( prdffile ) cout << "Fun4Muons_RecoDST_sim - prdffile : " << prdffile << endl;
    if( ndstfile ) cout << "Fun4Muons_RecoDST_sim - ndstfile : " << ndstfile << endl;

    bool is_pp = true;
    bool is_sim = false;
    if ( backgroundfile ){
        is_sim = false;
    }
    else {
        is_sim = true;
    }

    // use SVX software
    bool use_SVX = false;

    // use FVTX software
    bool use_FVTX = true;
    bool FVTX_mcreco = false;    // Take MC hits, optionally smear, make tracks and fit
    bool FVTX_realreco = true;   // Take FVTX clusters, make coords and tracks and fit

    // use RPC software
    bool use_RPC = false;

    // dst/ndst/pdst output
    bool write_dst = false;
    bool write_ndst = false;
    bool write_pdst = true;

    // MUTOO evaluator
    bool do_evaluation = true;
    bool fill_mutr_mc_info = true;

    // code profiler
    bool do_profile = false;

    // perfect detector geometry (no dead channels/HV; 100% muid efficiency)
    bool use_perfect_detector = false;

    /*
       when this flag is true, the mutr/muid hit maps merged after the embedding
       are packed again to a "merged" simulated PRDF
       */
    bool write_merged_prdf = false;

    // load libraries
    gSystem->Load("libfun4all");
    gSystem->Load("libmutoo_subsysreco");
    //gSystem->Load("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/13.run14/04.runQA/09.embedding_study/mutoo_subsysreco/install/lib/libmutoo_subsysreco");
    //gSystem->Load("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/12.hadron_cocktail/07.sngmuon/fvtxoo_modules_fix3/install/lib/libfvtxoo_modules");
    //gSystem->Load("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/12.hadron_cocktail/07.sngmuon/fvtx_subsysreco/install/lib/libfvtx_subsysreco");
    gSystem->Load("libfvtx_subsysreco");
    gSystem->Load("libfun4allfuncs_muons");
    gSystem->Load("libfun4allfuncs" );
    gSystem->Load("libcteval" );
    if( (ndstfile && write_ndst) || (write_pdst) ){
        gSystem->Load("libMWGOO");
    }

    if( do_evaluation )
    {
        // evaluation requires some extra reconstruction modules
        // from muon_subsysreco
        gSystem->Load("libmuon_subsysreco");
    }

    gSystem->Load("libMWG_interface");

    if( do_profile ) {
        gSystem->Load("libjprof");
        prof *Pr = new prof;
    }

    ///////////////////////////////////////////
    // recoConsts setup
    //////////////////////////////////////////
    recoConsts *rc = recoConsts::instance();
    rc->set_IntFlag("SVXACTIVE", 1);
    rc->set_IntFlag("RUNNUMBER", run_number);

    // Muon ext vtx verbosity
    // this allows to print which vertex is used and its value
    unsigned int use_vertex;
    if (is_sim){
        use_vertex = MuonUnpackSim::SIGNAL; // "SIGNAL", "BACKGROUND" or "RECO"
    }
    else{
        //use_vertex = MuonUnpackSim::RECO; // "SIGNAL", "BACKGROUND" or "RECO"
        use_vertex = MuonUnpackSim::SIGNAL; // "SIGNAL", "BACKGROUND" or "RECO"
    }

    TMutExtVtx::get().set_verbosity( MUTOO::ALOT );
    if (use_vertex == MuonUnpackSim::SIGNAL){
        TMutExtVtx::get().set_priorities( "SIM", 0 );
    }
    else{
        TMutExtVtx::get().set_priorities( "FVTX", 0 );
        TMutExtVtx::get().set_priorities( "SVX_PRECISE", 1 );
        TMutExtVtx::get().set_priorities( "BBC", 2 );
    }

    // Read the FVTX dead channel map and geometry from the database:
    // changes track finding parameters and vertex weight used in fvtx-mutr match
    TFvtxGlobalParCntrl::set_bool_par("is_pp",is_pp);   
    TFvtxGlobalParCntrl::set_bool_par("use_svx",use_SVX);   

    // Determine what DB parameters to read:
    TFvtxGlobalParCntrl::set_bool_par("is_sim",is_sim);
    TFvtxGlobalParCntrl::set_bool_par("deadmap_use_calibration_database",true);
    TFvtxGlobalParCntrl::set_bool_par("deadmap_use_production_map_for_sim",true);
    TFvtxGlobalParCntrl::set_bool_par("geom_use_calibration_database",false);
    if (is_sim){
        //TFvtxGlobalParCntrl::set_string_par("geom_root_file_path","/phenix/u/jinhuang/work/FVTX/miliped_work/millipede/");
        TFvtxGlobalParCntrl::set_string_par("geom_root_file_path","./");
        TFvtxGlobalParCntrl::set_string_par("geom_root_file_name","fvtx_geometry.root");
    }
    else{
        //TFvtxGlobalParCntrl::set_string_par("geom_root_file_path","/phenix/u/workarea/brooks/simulation/pisa2000/wrk/");
        //TFvtxGlobalParCntrl::set_string_par("geom_root_file_name", "fvtxgeom_v7.root");
        TFvtxGlobalParCntrl::set_string_par("geom_root_file_path","./");
        TFvtxGlobalParCntrl::set_string_par("geom_root_file_name","fvtx_geometry.root");
    }

    // database control initialization
    if( use_perfect_detector )
    {

        // tell muid tube efficiency must be set to 1
        rc->set_DoubleFlag("MUIOO_TUBE_EFF",1.0);

        // make perfect mutr detector
        TMutDatabaseCntrl::set_database_access("disable_HV",false);
        TMutDatabaseCntrl::set_database_access("attenuated_channels",false);
        TMutDatabaseCntrl::set_database_access("dead_channels",false);

    } else {

        /*
           to have muid inneficiencies turned on you need to have local copies of
           tube_eff_north_default.txt
           tube_eff_south_default.txt
           */

        ///////////////////////////////////////////
        // MUID setup
        //////////////////////////////////////////
        std::string muid_eff_south("muid_tube_eff_south_Run15pAu200.txt");
        std::string muid_eff_north("muid_tube_eff_north_Run15pAu200.txt");

        {
            std::cout << "using local two pack efficiciency files:" << muid_eff_south << ", "  << muid_eff_north << std::endl;
            TMuiHVMask::set_mode(TMuiHVMask::FROM_FILE);
            TMuiHVMask::set_filename_south(muid_eff_south.c_str());
            TMuiHVMask::set_filename_north(muid_eff_north.c_str());
        }

        /*
           to have mutr HV mask turned on, you need to have the following set to true
           and a local copy of mut.disabledAnodes.dat
           */

        TMutDatabaseCntrl::set_database_access("use_local_dead_HV_file",true); //mut.disabledAnodes.dat
        TString hv_file = Form("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/07.genericQA/01.mutr_database/Run15pAu200/mut.disabledAnodes_Run15pAu200_updated/mut.disabledAnodes.dat_run%d",run_number);
        TMutDatabaseCntrl::set_filename("use_local_dead_HV_file", hv_file.Data());
        TMutDatabaseCntrl::set_database_access("use_local_dead_wire_file",true); //mut.disabledWires.dat
        TMutDatabaseCntrl::set_database_access("use_local_attenuated_chnl_file",true); //AttenuatedChannels.txt
        TMutDatabaseCntrl::set_database_access("use_local_dead_channel_file",true); //DeadChannels.dat
        TMutDatabaseCntrl::set_database_access("use_local_landau_parameters_file",true);
        TMutDatabaseCntrl::set_filename("use_local_landau_parameters_file", "landau_parameters.txt");
        TMutDatabaseCntrl::set_database_access("use_local_st1_autocad_file",true); //0xx_x.dat
        TMutDatabaseCntrl::set_database_access("use_local_dcm_map_file",true); //mut.St1ChannelMap.North
        TMutDatabaseCntrl::set_database_access("use_local_octant_survey_file",true); //MutOctantPositions.dat
    }

    //TMutDatabaseCntrl::set_verbosity( TMutDatabaseCntrl::GLOBAL_ALIGN, TMutDatabaseCntrl::NONE );
    //TMutDatabaseCntrl::set_verbosity( TMutDatabaseCntrl::INTERNAL_ALIGN, TMutDatabaseCntrl::NONE );
    //TMutDatabaseCntrl::print();

    ///////////////////////////////////////////
    // Make the Server
    //////////////////////////////////////////
    Fun4AllServer *se = Fun4AllServer::instance();
    se->Verbosity(0);

    ///////////////////////////////////////////
    // Subsystems
    //////////////////////////////////////////

    // counter
    se->registerSubsystem( new MuonCounter() );

    //MuonTrigFilter *filter_bbcq = new MuonTrigFilter("BBC_CHARGE",MuonTrigFilter::BBC_CHARGE,MuonTrigFilter::ABORT_EVENT);
    //filter_bbcq->set_bbc_charge_window(0,150);
    //se->registerSubsystem( filter_bbcq );

    // this part should no be necessary unless you run on an 'old' simulated DST
    // (for which VtxSimreco did not run in Fun4Muons_PISA
    bool use_vtx_sim_reco = false;
    if( use_vtx_sim_reco )
    {
        gSystem->Load("libsimreco_base");
        VtxSimreco* vtx_reco = new VtxSimreco();
        vtx_reco->SmearZ( true );
        vtx_reco->UseXY( false );
        vtx_reco->OverwriteBBC( true );
        //vtx_reco->OverwriteBBC( false );
        vtx_reco->ZVertexSigma( 0.00 );
        se->registerSubsystem( vtx_reco );
    }

    // this runs the SVX reconstruction (stand-alone)
    if( use_SVX )
    {
        // VTX Reconstruction software:
        gSystem->Load("libsimreco");
        //gSystem->Load("libembedreco");
        gSystem->Load("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/embedreco/install/lib/libembedreco");
        gSystem->Load("libsvx");

        SvxParManager *svxpar = new SvxParManager();
        svxpar->Verbosity(0);
        //svxpar->set_BeamCenter(0.,0.);
        svxpar->set_OffsetVtxToCnt(0.0, 0.0, 0.0);
        svxpar->set_OffsetEastToWest(0.0, 0.0, 0.0);
        svxpar->set_ReadGeoParFromFile(1);
        svxpar->set_GeometryFileName("svxPISA.par");
        se->registerSubsystem(svxpar);

        if ( is_sim ){
            SvxSimulator *svxsim = new SvxSimulator();
            svxsim->Verbosity(0);
            se->registerSubsystem(svxsim);
        }

        if ( backgroundfile )
        {
            SvxEmbedSimhit *svxembed = new SvxEmbedSimhit();
            svxembed->set_StripixelNoise(0.0); // no noise
            svxembed->Verbosity(0);
            se->registerSubsystem(svxembed);
        }

        /*
           SvxDecode *svxdecode = new SvxDecode();
           svxdecode->Verbosity(0);
           svxdecode->includePixel(true);
           svxdecode->includeStripixel(true);
           svxdecode->setAdcOffset(24);
           svxdecode->setAdcCutoff(-24);
           se->registerSubsystem(svxdecode);
           */

        SvxApplyHotDead *svxhotdead = new SvxApplyHotDead();
        svxhotdead->Verbosity(0);
        se->registerSubsystem(svxhotdead);

        SvxReco *svxrec = new SvxReco();
        svxrec->Verbosity(0);
        //svxrec->Load_ThresholdFile("threshold.h");
        svxrec->set_ThisIsSimulation();
        //svxrec->set_UseStripThresholdDatbase(true);
        //svxrec->Load_ThresholdFile("threshold_ideal.h");
        svxrec->set_StripixelAdcSumThreshold(0);
        se->registerSubsystem(svxrec);

        // Skip for inclusive HF muon analysis
        SvxPriVertexSeedFinder *svxvtxseedfinder = new SvxPriVertexSeedFinder();
        svxvtxseedfinder->Verbosity(0);
        //se->registerSubsystem(svxvtxseedfinder);

        // Skip for inclusive HF muon analysis
        SvxStandAloneReco *svxstandalone = new SvxStandAloneReco();
        svxstandalone->Verbosity(0);
        svxstandalone->setVertexRecoFlag(2);
        svxstandalone->setPPFlag(true);
        //se->registerSubsystem(svxstandalone);

        // Skip for inclusive HF muon analysis
        SvxPrimVertexFinder *svxprimvtxfinder = new SvxPrimVertexFinder();
        svxprimvtxfinder->Verbosity(0);
        //se->registerSubsystem(svxprimvtxfinder);
    }


    // Mutr/Muid unpacking
    MuonUnpackSim* muon_unpack = new MuonUnpackSim();
    // Treat even MC backgrounds as "real" backgrounds because there isn't a utility
    // to merge MC information. In this case, only the "signal" MC information will
    // show up in the evaluation ntuples:
    if( backgroundfile ) muon_unpack->SetMode( MuonUnpackSim::MC_SIGNAL_REAL_BG );
    else muon_unpack->SetMode( MuonUnpackSim::MC_SIGNAL_NO_BG );

    muon_unpack->SetVtxToUse( use_vertex );

    //muon_unpack->set_flag( MuonUnpackSim::NO_CHARGE_SMEAR, true );
    se->registerSubsystem( muon_unpack );

    // RPC unpacking/reconstruction
    if( use_RPC )
    {
        gSystem->Load("librpc_subsysreco");
        //if( backgroundfile ) se->registerSubsystem( new RpcUnpackSim( "RPCUNPACKSIM", RpcUnpackSim::MC_SIGNAL_REAL_BG ) );
        if( backgroundfile ) se->registerSubsystem( new RpcUnpackSim( "RPCUNPACKSIM", RpcUnpackSim::MC_SIGNAL_MC_BG ) );
        else  se->registerSubsystem( new RpcUnpackSim( "RPCUNPACKSIM", RpcUnpackSim::MC_SIGNAL_NO_BG ) );
        se->registerSubsystem( new RpcReco() );
    }

    if( use_FVTX )
    {
        if( FVTX_mcreco ) {

            FvtxUnpackSim* fvtx_unpack = new FvtxUnpackSim();
            if( backgroundfile ) fvtx_unpack->SetMode( FvtxUnpackSim::MC_SIGNAL_MC_BG );
            else fvtx_unpack->SetMode( FvtxUnpackSim::MC_SIGNAL_NO_BG );
            se->registerSubsystem( fvtx_unpack );

            se->registerSubsystem( new FvtxRecoMC() );

        } else if( FVTX_realreco ) {

            FvtxUnpackSim* fvtx_unpack = new FvtxUnpackSim();
            // Treat even MC backgrounds as "real" backgrounds because there isn't a utility
            // to merge MC information. In this case, only the "signal" MC information will
            // show up in the evaluation ntuples:
            //if( backgroundfile ) fvtx_unpack->SetMode( FvtxUnpackSim::MC_SIGNAL_REAL_BG );
            //else fvtx_unpack->SetMode( FvtxUnpackSim::MC_SIGNAL_NO_BG );
            //No need to do embedding FVTX hit for inclusive HF muon analysis
            fvtx_unpack->SetMode( FvtxUnpackSim::MC_SIGNAL_NO_BG );
            se->registerSubsystem( fvtx_unpack );

            FvtxReco *fvtx_reco = new FvtxReco();
            //fvtx_reco->run_after_burner(false);            // Adds VTX hits to tracking
            //fvtx_reco->set_use_svx_cluster(true);          // Uses smeared PISA hits if false
            //fvtx_reco->set_do_mutr_matching(false);         // Match FVTX-MuTr tracks and refit or not
            //fvtx_reco->set_finder(3);                       // Use Hough track finder
            se->registerSubsystem( fvtx_reco );

            //TMutNode<mFvtxKalFitPar>::find_node(se->topNode(),"mFvtxKalFitPar")->set_do_evaluation(true);
            //TMutNode<mFvtxKalFitPar>::find_node(se->topNode(),"mFvtxKalFitPar")->set_vtx_phi_error(0.005);
            //TMutNode<mFvtxKalFitPar>::find_node(se->topNode(),"mFvtxKalFitPar")->set_vtx_r_error(0.015);
            //TMutNode<mFvtxKalFitPar>::find_node(se->topNode(),"mFvtxKalFitPar")->set_vtx_z_error(0.012);
            //TMutNode<mMutKalFitWithSiliRealPar>::find_node(se->topNode(),"mMutKalFitWithSiliRealPar")->set_do_evaluation(true);
            //TMutNode<mMutKalFitWithSiliRealPar>::find_node(se->topNode(),"mMutKalFitWithSiliRealPar")->set_vtx_phi_error(0.005);
            //TMutNode<mMutKalFitWithSiliRealPar>::find_node(se->topNode(),"mMutKalFitWithSiliRealPar")->set_vtx_r_error(0.015);
            //TMutNode<mMutKalFitWithSiliRealPar>::find_node(se->topNode(),"mMutKalFitWithSiliRealPar")->set_vtx_z_error(0.012);
            //TMutNode<mFvtxKalFitPar>::find_node(se->topNode(),"mFvtxKalFitPar")->set_vtx_phi_error(1.005);
            //TMutNode<mFvtxKalFitPar>::find_node(se->topNode(),"mFvtxKalFitPar")->set_vtx_r_error(1.015);
            //TMutNode<mFvtxKalFitPar>::find_node(se->topNode(),"mFvtxKalFitPar")->set_vtx_z_error(1.012);

            // Set the vertex sigma that is use in FVTX-MuTr track matching:
            if (is_pp)
            {
                //TMutNode<mMutKalFitWithSiliRealPar>::find_node(se->topNode(),"mMutKalFitWithSiliRealPar")->set_vtx_sigma_z(25.0);
                //TMutNode<mMutKalFitWithSiliRealPar>::find_node(se->topNode(),"mMutKalFitWithSiliRealPar")->set_vtx_sigma_xy(10.0);
                //TMutNode<mFvtxFindTrackPar>::find_node(se->topNode(),"mFvtxFindTrackPar")->set_rCutFactor(rCutFactor);
            }
            else
            {
                //TMutNode<mMutKalFitWithSiliRealPar>::find_node(se->topNode(),"mMutKalFitWithSiliRealPar")->set_vtx_sigma_z(1.0);
                //TMutNode<mMutKalFitWithSiliRealPar>::find_node(se->topNode(),"mMutKalFitWithSiliRealPar")->set_vtx_sigma_xy(0.5);
            }


            // primary vertex from FVTX
            //if ( use_vertex == MuonUnpackSim::RECO ){
            if ( 0 ){
                FvtxPrimVertex* fvtxprimvtx = new FvtxPrimVertex();
                //fvtxprimvtx->set_fvtx_Rres(1.0);  // changed to 1cm, recomended by Mike
                fvtxprimvtx->set_source(FvtxPrimVertex::Tracks);
                //fvtxprimvtx->set_fitter_active(true);
                fvtxprimvtx->set_clustering(FvtxPrimVertex::AllInOne); // AllInOne for single vertex, Cesars for multiple verteces
                //fvtxprimvtx->set_bbcz_window(100.0); // FVTX-BBC vertex matching window

                // Not clear this has an effect for simulated events:
                //fvtxprimvtx->set_mutr_use_fvtx( false );  // Do not use FVTX as primary vertex
                //fvtxprimvtx->set_mutr_use_fvtx( true );  // Do not use FVTX as primary vertex

                se->registerSubsystem(fvtxprimvtx);
            }
        }

        // Mutr/Muid reconstruction
        // Note:  MuonDev must be called before FvtxRecoWithMut

        MuiooReco* muioo_reco = new MuiooReco();
        // Larger cut values needed if you are simulating CuAu data sets:
        //muioo_reco->set_asymm_cut_par(100);
        muioo_reco->set_max_occupancy_per_arm(1000);
        se->registerSubsystem( muioo_reco );

        MuonDev* muon_dev = new MuonDev();
        //muon_dev->set_flag(MuonDev::SKIP_VERTEX_FIT, 1);
        se->registerSubsystem( muon_dev );

        // Perform FVTX-Mutr track matching and refit track:
        FvtxRecoWithMut *fvtx_reco_withmut = new FvtxRecoWithMut();
        se->registerSubsystem( fvtx_reco_withmut );
        }

        // official muon evaluation modules
        if( do_evaluation )
        {
            //MuonAnaTuples* ana_tuple = new MuonAnaTuples( "MUONANATUPLES" );
            //se->registerSubsystem( ana_tuple );
            //ana_tuple->set_event_filename( "event_vtx.txt");
            //ana_tuple->set_flags( MuonAnaTuples::VERTEX | MuonAnaTuples::CENTRALITY);// | MuonAnaTuples::REACTION_PLANE);

            //se->registerSubsystem( new MuonAnaTuples() );
            //se->registerSubsystem( new MuonEval() );

            /*
               cout << "Fun4Muons_RecoDST_sim - registering residual evaluation" << endl;
               const char* resolution_file = "MutooResEval.root";
               MutooResEval* res_eval = new MutooResEval( "MUTOORESEVAL", resolution_file );
               res_eval->set_refit( true );
               res_eval->set_magnets_on( true );
               se->registerSubsystem( res_eval );
               */
        }
        // Fill FVTX evaluation ntuples:
        if( do_evaluation )
        {
            //se->registerSubsystem( new MuonEval() );
            MuonEval* mueval = new MuonEval();
            mueval->set_flags(0);
            se->registerSubsystem(mueval);

            if ( use_FVTX ){
                //FvtxEval* fvtxeval = new FvtxEval("FvtxEval","fvtx_eval_embed_pisa.root");
                //se->registerSubsystem(fvtxeval);
                //
                FvtxMCEval* fvtxeval_mc = new FvtxMCEval("FvtxMCEval","fvtx_mc_eval_embed_pisa.root");
                se->registerSubsystem(fvtxeval_mc);
            }
        }

        /*
           global Reco must be run after the reconstruction
           but before the evalutations for all nodes/ntuples to be filled
           properly
           */
        se->registerSubsystem( new GlobalReco() );
        se->registerSubsystem( new GlobalReco_muons() );

        // prdf packer
        if( prdffile && write_merged_prdf )
            se->registerSubsystem( new MuonPackPRDF() );

        ///////////////////////////////////////////
        // Input manager
        ///////////////////////////////////////////
        rc->set_CharFlag("EMBED_MC_TOPNODE","SIGNAL");
        rc->set_CharFlag("EMBED_REAL_TOPNODE","BACKGROUND");
        Fun4AllInputManager *signal = new Fun4AllNoSyncDstInputManager( "SIGNAL_IM", "DST", "SIGNAL" );
        se->registerInputManager( signal );
        se->fileopen( signal->Name(), signalfile );

        // try load simulated background, if requested
        if( backgroundfile )
        {
            Fun4AllInputManager *background = new Fun4AllNoSyncDstInputManager( "BACKGROUND_IM", "DST", "BACKGROUND" );
            se->registerInputManager(background);
            se->fileopen( background->Name(),backgroundfile );
        }

        /*
           try load background also as a head dst, to get global detectors
           this should not interfere with mutoo nodes (even if duplicated) since they
           are overwritten by the PHMapManager
           */
        if( backgroundfile )
        {
            cout << "Fun4Muons_RecoDST_sim - reoppening BACKGROUND under DST node to get additionnal nodes copied.\n";
            Fun4AllInputManager *background = new Fun4AllNoSyncDstInputManager( "BACKGROUND_DST_IM", "DST" );
            se->registerInputManager(background);
            se->fileopen( background->Name(),backgroundfile );

        } else if( signalfile ) {
            cout << "Fun4Muons_RecoDST_sim - reoppening SIGNAL under DST node to get additionnal nodes copied.\n";
            Fun4AllInputManager *signal = new Fun4AllNoSyncDstInputManager( "SIGNAL_DST_IM", "DST" );
            se->registerInputManager( signal );
            se->fileopen( signal->Name(), signalfile );
        }

        ///////////////////////////////////////////
        // Output manager
        ///////////////////////////////////////////

        // dst output manager
        if( dstfile && write_dst )
        {

            Fun4AllDstOutputManager *dstManager  = new Fun4AllDstOutputManager("DSTOUT",dstfile);
            se->registerOutputManager(dstManager);

            // Header and vertex nodes
            dstManager->AddNode("RunHeader");
            dstManager->AddNode("EventHeader");
            dstManager->AddNode("VtxOut");
            dstManager->AddNode("BbcOut");
            dstManager->AddNode("BbcRaw");
            dstManager->AddNode("ZdcOut");
            dstManager->AddNode("ZdcRaw");
            dstManager->AddNode("TrigLvl1");

            dstManager->AddNode("TMuiPseudoLL1");

            // muioo nodes
            dstManager->AddNode("TMCPrimary");
            dstManager->AddNode("TMuiMCHitO");
            dstManager->AddNode("TMuiHitO");
            dstManager->AddNode("TMuiClusterO");
            dstManager->AddNode("TMui1DRoadO");
            dstManager->AddNode("TMuiRoadO");
            dstManager->AddNode("TMuiPseudoBLTO");
            dstManager->AddNode("TMuiPseudoLL1");

            // mutoo nodes
            dstManager->AddNode("TMutMCHit");
            dstManager->AddNode("TMutMCTrk");
            dstManager->AddNode("TMutHit");
            dstManager->AddNode("TMutClus");
            dstManager->AddNode("TMutCoord");
            dstManager->AddNode("TMutGapCoord");
            dstManager->AddNode("TMutStub");
            dstManager->AddNode("TMutTrk");
            dstManager->AddNode("TMutVtx");

            // rpc nodes
            if( use_RPC )
            {
                dstManager->AddNode("TRpcMCHit");
                dstManager->AddNode("TRpcHit");
                dstManager->AddNode("TRpcClus");
                dstManager->AddNode("TRpcCoord");
                dstManager->AddNode("TRpcTrk");
            }

            // Fvtx nodes
            if( use_FVTX )
            {
                dstManager->AddNode( "TFvtxMCHit" );
                dstManager->AddNode( "TFvtxMCTrk" );
            }

            // From EVA node
            dstManager->AddNode("header");
            dstManager->AddNode("fkin");
            dstManager->AddNode("primary");
            dstManager->AddNode("pythia");
        }

        // ndst output manager
        if( ndstfile && write_ndst )
        {

            Fun4AllDstOutputManager *ndstManager = new Fun4AllDstOutputManager( "MWGOUT",ndstfile );
            se->registerOutputManager(ndstManager);

            ndstManager->AddNode("RunHeader");
            ndstManager->AddNode("EventHeader");
            ndstManager->AddNode("header");
            ndstManager->AddNode("BbcOut");
            ndstManager->AddNode("PHGlobal");
            ndstManager->AddNode("PHGlobal_MUON");
            ndstManager->AddNode("PHMuoTracksOO");
            ndstManager->AddNode("TrigLvl1");
            ndstManager->AddNode("TMuiPseudoBLTO");
            ndstManager->AddNode("TMuiPseudoLL1");
            ndstManager->AddNode("VtxOut");

            ndstManager->AddNode("RpSumXYObject");
            ndstManager->AddNode("ReactionPlaneObject");

        }

        // picodst output manager
        if( write_pdst )
        {
            // global Reco
            se->registerSubsystem( new GlobalReco() );

            // MWG
            PHInclusiveNanoCuts *MWGcuts = new MWGInclusiveNanoCutsv2();
            se->registerSubsystem(new MWGFvtxReco(MWGcuts));

            //gSystem->Load("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/13.run14/04.runQA/09.embedding_study/picoDST_object/install/lib/libpicodst_object.so");
            //gSystem->Load("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/12.hadron_cocktail/07.sngmuon/picoDST_object/install/lib/libpicodst_object.so");
            gSystem->Load("libpicodst_object.so");
            mFillSingleMuonContainer* msngl = new mFillSingleMuonContainer();
            //msngl->set_bbcz_cut(20.0);
            //msngl->set_pz_cut(1.0);
            //msngl->set_lastgap_cut(2);
            //msngl->set_DG0_cut(30);
            //msngl->set_DDG0_cut(20);
            //msngl->set_chi2_cut(20);
            //msngl->set_nidhits_cut(1);
            msngl->set_is_sim(true);
            se->registerSubsystem(msngl);

            /*
               mFillDiMuonContainer* mdi = new mFillDiMuonContainer(false); // do not make mixed events
               se->registerSubsystem(mdi);
               mdi->set_bbcz_cut(100.0);
               mdi->set_mass_cut(0.5);
               mdi->set_is_sim(false);
               mdi->set_is_pp(is_pp);
               */

            /*
               Fun4AllOutputManager *outdimu = new Fun4AllDstOutputManager("Outdimu",dimuonpdstout);
               outdimu->AddNode("Sync");
               outdimu->AddNode("DiMuonContainer");
               outdimu->AddNode("PHPythiaHeader");
               outdimu->AddNode("PHPythia");
               outdimu->AddEventSelector("mFillDiMuonContainer");
               if(fill_mutr_mc_info) outdimu->AddNode("MCDiMuonContainer"); // Mutr MC information container
               se->registerOutputManager(outdimu);
               */

            if(fill_mutr_mc_info)
            {
                //********************* Mutr MC information Module ******************//
                //mFillMCSingleMuonContainer* msngl_mc = new mFillMCSingleMuonContainer();
                mFillMCSingleMuonFvtxContainer* msngl_mc = new mFillMCSingleMuonFvtxContainer();
                se->registerSubsystem(msngl_mc);

                //mFillMCDiMuonContainer* mdi_mc = new mFillMCDiMuonContainer();
                //se->registerSubsystem(mdi_mc);
            }

            Fun4AllOutputManager *outsmu = new Fun4AllDstOutputManager("Outsmu",singlepdstout);
            outsmu->AddNode("SingleMuonContainer");
            outsmu->AddNode("Sync");
            outsmu->AddNode("VtxOut");
            outsmu->AddNode("PHGlobal");
            outsmu->AddNode("RunHeader");
            outsmu->AddNode("EventHeader");
            //outsmu->AddNode("PHPythiaHeader");
            //outsmu->AddNode("PHPythia");
            outsmu->AddEventSelector("mFillSingleMuonContainer");
            //if (fill_mutr_mc_info) outsmu->AddNode("MCSingleMuonContainer"); // Mutr MC information container
            if (fill_mutr_mc_info) outsmu->AddNode("MCSingleMuonFvtxContainer"); // Mutr MC information container
            se->registerOutputManager(outsmu);

        }

        // prdf output
        if( prdffile && write_merged_prdf )
        {
            Fun4AllPrdfOutputManager *prdf_io  = new Fun4AllPrdfOutputManager("PRDFOUT",prdffile);
            prdf_io->InitPrdfNode( se->topNode() );
            se->registerOutputManager(prdf_io);
        }

        ///////////////////////////////////////////
        // Analyze the Data.
        //////////////////////////////////////////
        MuonUtil::dump_process_time();
        gSystem->ListLibraries();
        se->run(nEvents);
        se->End();
        MuonUtil::dump_process_time();
        cout << "Completed reconstruction." << endl;

    }
