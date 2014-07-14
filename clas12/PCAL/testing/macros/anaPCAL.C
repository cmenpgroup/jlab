void anaPCAL(string input_file="in.root", Long64_t max=0, Long64_t min =0)
{
    if (!TClassTable::GetDict("PCAL")) {
        gSystem->Load("libclas12banks.so");
    }
    if (!TClassTable::GetDict("PCALHB")) {
        gSystem->Load("libclas12banks.so");
    }

    cout << "Using inout file: " << input_file << endl;
    TFile* f = TFile::Open(input_file.c_str());
    if (f == 0) {
        printf("Error: cannot open file\n");
        return;
    }

    TTree* tree = (TTree*) f->Get("CLASEVENT");
    
    TClonesArray* PCALstore = new TClonesArray("PCAL", 1);
    TClonesArray* PCALHBstore = new TClonesArray("PCALHB", 1);
    
    TBranch* PCALbranch = tree->GetBranch("PCAL");
    TBranch* PCALHBbranch = tree->GetBranch("PCALHB");

    PCALbranch->SetAutoDelete(kFALSE);
    PCALbranch->SetAddress(&PCALstore);
    
    PCALHBbranch->SetAutoDelete(kFALSE);
    PCALHBbranch->SetAddress(&PCALHBstore);

    Long64_t nentries = tree->GetEntries();

    Long64_t StopEvent;
    Long64_t StartEvent;
    
    if(max){
      StopEvent = max;
    }else{
      StopEvent = nentries;
    }
    
    if(min){
    StartEvent = min;
    }else{
      StartEvent = 0;
    }

    for (Long64_t ev = StartEvent; ev < StopEvent; ++ev) {
        
        PCALstore->Clear();
	PCALHBstore->Clear();
        
        tree->GetEntry(ev);
        
        Int_t PCALrows = PCALstore->GetEntriesFast();
        Int_t PCALHBrows = PCALHBstore->GetEntriesFast();
	
        cout << "Event " << ev << endl;
	
        cout <<"*******PCAL data*******:"<<endl;
        for (Long64_t row = 0; row < PCALrows; ++row) {
            PCAL* data = (PCAL*) PCALstore->At(row);
            cout << "***Row***: "<< row << endl;
            cout << "ETot: " << data->ETot << "  E: " << data->vx << endl;
            cout << "x_avg: " << data->x_avg << "  y_avg: " << data->y_avg << "  z_avg: " << data->z_avg << endl;
            cout << "lx_avg: " << data->lx_avg << "  ly_avg: " << data->ly_avg << "  lz_avg: " << data->lz_avg << endl;
            cout << "t_avg: " << data->t_avg << endl;
            cout << "pid: " << data->pid << endl;
            cout << "vx: " << data->vx << "  vy: " << data->vy << "  vz: " << data->vz << endl;
            
            //cout << "mpid: " << data->mpid << endl;
            //cout << "mvx " << data->mvx << endl;
            //cout << "mvy " << data->mvy << endl;
            //cout << "mvz " << data->mvz << endl;
            cout << "Sector: " << data->sector << endl;
            cout << "Layer: " << data->layer << endl;
            cout << "View: " << data->view << endl;
            cout << "Strip: " << data->strip << endl;
            cout << "ADC: " << data->ADC << endl;
            cout << "TDC: " << data->TDC << endl;
	    cout << endl;
        }
        
        cout <<"*******PCALHB data:*******"<<endl;
        for (Long64_t row = 0; row < PCALHBrows; ++row) {
            PCALHB* data2 = (PCALHB*) PCALHBstore->At(row);
            cout << "***Row***: " << row << endl;
            cout << "E__hit: " << data2->E__hit << endl;
            //cout << "dE_hit: " << data2->dE_hit << endl;
            cout << "t_hit: " << data2->t_hit << endl;
            //cout << "dt_hit: " << data2->dt_hit << endl;
            cout << "i_hit: " << data2->i_hit << "  j_hit: " << data2->j_hit << endl;
            cout << "di_hit: " << data2->di_hit << "  dj_hit: " << data2->dj_hit << endl;
            //cout << "x_hit: " << data2->x_hit << endl;
            //cout << "y_hit: " << data2->y_hit << endl;
            //cout << "z_hit: " << data2->z_hit << endl;
            //cout << "dx_hit: " << data2->dx_hit << endl;
            //cout << "dy_hit: " << data2->dy_hit << endl;
            //cout << "dz_hit: " << data2->dz_hit << endl;
            cout << "u2_hit: " << data2->u2_hit << "  v2_hit: " << data2->v2_hit << "  w2_hit: " << data2->w2_hit << endl;
            cout << "u3_hit: " << data2->u3_hit <<  "  v3_hit: " << data2->v3_hit << "  w3_hit: " << data2->w3_hit << endl;
            cout << "u4_hit: " << data2->u4_hit << "  v4_hit: " << data2->v4_hit << "  w4_hit: " << data2->w4_hit << endl;
            cout << "center_U: " << data2->center_U << "  center_V: " << data2->center_V << "  center_W: " << data2->center_W << endl;
            cout << "path_U: " << data2->path_U << "  path_V: " << data2->path_V << "  path_W: " << data2->path_W << endl;
            //cout << "CH21: " << data2->CH21 << endl;
            //cout << "CH22: " << data2->CH22 << endl;
            cout << "Sector: " << data2->Sector << endl;
            cout << "Layer: " << data2->Layer << endl;
            cout << "Nstrip_U: " << data2->Nstrip_U << "  Nstrip_V: " << data2->Nstrip_V << "  Nstrip_W: " << data2->Nstrip_W << endl;
            //cout << "MatchID1: " << data2->MatchID1 << endl;
            //cout << "MatchID2: " << data2->MatchID2 << endl;
            cout << "istat: " << data2->istat << endl;
            cout << endl;

        }

    }
}
