#include <iostream>
#include <iomanip>

using namespace std;

double sqr(const double &x) {return x*x;}

void make_pisa_oscar_input(int npart=1000, int myseed=0){

    gRandom->SetSeed(myseed);

    //input file name, inlcuding Eta and PT distribution of each PID particle
    char *fname = "input_hist_Run8dAu200.root";
    char hname[300];

    //pid in PISA
    int PISApid[8] = {8, 9, 11, 12, 10, 16, 14, 15};
    //pid in PDG
    int PDGpid[8] = {211, -211, 321, -321, 130, 310, 2212, -2212};
    //Eta range
    int etamin[14] = {-240, -220, -200, -180, -160, -140, -120, 100, 120, 140, 160, 180, 200, 220};
    int etamax[14] = {-220, -200, -180, -160, -140, -120, -100, 120, 140, 160, 180, 200, 220, 240};

    //particle name
    const char* pid_label[8] = {"pi^{+}", "pi^{-}", "K^{+}", "K^{-}", "K_{L}^{0}", "K_{S}^{0}", "p", "pbar"};
    //particle mass
    double mass[8] = {0.139570, 0.139570, 0.493667, 0.493667, 0.497614, 0.497614, 0.938272, 0.938272};

    //read the file
    TFile *infile = new TFile(fname);

    //prepare histograms
    TH1D *PID = infile->Get("PIDmod");
    TH1D *ETA[8];
    TH1D *hPT[8][14];

    double Npion = 0, Nkaon = 0, Nproton = 0;

    /////////////////
    //LOADING HISTS//
    /////////////////
    for (int ipid=0; ipid<8; ipid++){
        sprintf(hname, "ETA_PID%d", PISApid[ipid]);
        cout << "LOADING : " << hname << endl;
        ETA[ipid] = (TH1D*)infile->Get(hname);
        for (int ieta=0; ieta<14; ieta++){
            sprintf(hname, "PT4_PID%d_ETA%d_%d", PISApid[ipid], etamin[ieta], etamax[ieta]);
            cout << "LOADING : " << hname << endl;
            hPT[ipid][ieta] = (TH1D*)infile->Get(hname);
        }
    }

    //BBC-Z distribution from Run15pAu200MB BBCLL1 +/- 30 cm triggered events 
    TF1 *fBBCZ = new TF1("fBBCZ","pol3",-40,40);
    fBBCZ->SetParameters(3.62e+09, 2.52e+07, -1.43e+06, -19172.7);

    double pid, pT, eta;
    int _pid, _eta;
    double _ptot, _theta, _vertex, _phi;
    double _xvtx, _yvtx, _zvtx;
    TVector3 vec;

    TRandom3 *rand = new TRandom3();
    rand->SetSeed(myseed);

    ofstream fout;
    fout.open("oscar.input");

    ///////////////////
    //MAKE SOME DELAY//
    ///////////////////
    double bignum = 1.0e7*rand->Rndm();
    cout << bignum << "\twaiting for loop start..." << endl; 

    int pidbin, etabin;

    for (double icount=0; icount<bignum; icount++){
    }

    //////////////////
    //SET SEED AGAIN//
    //////////////////
    //generate uniformly-distributed floating points between 0 and 1
    gRandom->SetSeed(rand->Rndm());
    //rand->SetSeed(0);

    cout << "starting particle generation..." << endl;

    fout << "# OSC1999A" << endl;
    fout << "# final_id_p_x" << endl;
    fout << "# SimName 1.0" << endl;
    fout << "#" << endl;

    for (int ii=0; ii<npart; ii++){
        ///////////
        //GET PID//
        ///////////
        //Return a random number distributed according the histogram bin contents
        pid = PID->GetRandom();
        _pid = (int)pid;
        for (int ipid=0; ipid<8; ipid++){
            if ( PISApid[ipid]==_pid ){
                pidbin = ipid;
                break;
            }
        }

        if ( _pid==8 ) Npion++;
        else if ( _pid==11 ) Nkaon++;
        else if ( _pid==14 ) Nproton++;

        ///////////
        //GET ETA//
        ///////////
        //Return a random number distributed according the histogram bin contents
        eta = ETA[pidbin]->GetRandom(); 
        for (int ieta=0; ieta<14; ieta++){
            if ( 100*eta>etamin[ieta] && 100*eta<etamax[ieta] ){
                etabin = ieta;
                break;
            }
        }

        //////////
        //GET PT//
        //////////
        //Return a random number distributed according the histogram bin contents
        pT = hPT[pidbin][etabin]->GetRandom();

        _theta = 2*TMath::ATan( TMath::Exp(-eta) )*180.0/TMath::Pi();
        _phi = rand->Rndm() * 360.0;
        _ptot = pT / TMath::Sin(TMath::Pi()*_theta/180.0);

        vec.SetPtEtaPhi(pT, eta, _phi);
        const double energy = sqrt(sqr(vec.X()) + sqr(vec.Y()) + sqr(vec.Z()) + sqr(mass[pidbin]));

        _xvtx = _yvtx = 0.0;
        _zvtx = fBBCZ->GetRandom(); 

        fout << "0  1" << endl;

        fout << "0" << " "
            << PDGpid[pidbin] << " "
            << "0" << " "
            << vec.X() << " "
            << vec.Y() << " "
            << vec.Z() << " "
            << energy << " "
            << mass[pidbin] << " "
            << _xvtx*1e13 << " "
            << _yvtx*1e13 << " "
            << _zvtx*1e13 << " "
            << "0" << endl;
        fout << "0  0" << endl;
    }

    cout << "PION : " << Npion;
    cout << "\tKAON : " << Nkaon;
    cout << "\tPROTON : " << Nproton << endl;
    }
