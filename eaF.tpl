 // this version estimates annual F
GLOBALS_SECTION
   #include <time.h>
DATA_SECTION
  init_int mode
  init_number quota
  init_number WeightMULTI
  init_number WeightDEPM
  init_number WeightM
  init_number WeightH
  init_number WeightDEV
  init_number p1
  init_number p2
  init_number p3
  init_number p4
  init_number p5
  init_number p6
  init_number p7
  init_number p8
  init_int fyear
  init_int dyear
  init_int lyear
  init_int pyear
  init_int fage
  init_int lage
  init_int LegSize
  init_number PropFem
  init_number SDcorrect
  init_number Mprior
  init_number MpriorSD
  init_number steepprior
  init_number steeppriorSD
  init_vector qgrow(fyear,lyear)
  init_vector age(fage,lage)
  init_vector SexMat(fage,lage)
  init_vector retain1(fage,lage)
  init_vector retain2(fage,lage)
  init_vector WtF(fage,lage)
  init_vector WtM(fage,lage)
  init_vector OCat(fyear,pyear)
  init_vector OCatSD(fyear,pyear)
  init_number Ndepm
  init_vector lowDEPM(dyear,lyear)
  init_vector DEPM(dyear,lyear)
  init_vector uppDEPM(dyear,lyear)
  init_vector cvDEPM(dyear,lyear)
  init_vector sdDEPM(dyear,lyear)
  init_vector lnDEPM(dyear,lyear)
  init_vector lnsdDEPM(dyear,lyear)
  init_vector NageSample(dyear,lyear)
  init_matrix SurAge(dyear,lyear,fage,lage)
  init_number checknum
  !! ad_comm::change_datafile_name("prop_age_composition.dat");
  init_matrix p_age(dyear,lyear,fage,lage)
  !! ad_comm::change_datafile_name("pred_prop_age_composition.dat");
  init_matrix hatp_age(dyear,lyear,fage,lage)
  !! ad_comm::change_datafile_name("value_of_phi.dat");
  init_number init_phi
    
  int a
  int i
  int j
  int seed
  int sim
   
PARAMETER_SECTION
  //MARK
  init_number Rstar(p1)
  init_number Q(p2)
  init_number Fo(p3)
  init_number M(p4)
  init_number steep(p5)
  init_number SelA50(p6)
  init_number SelSlope(p6)
  init_vector dev(fyear,lyear,p7)
  init_vector fut_dev(lyear+1,pyear,p7)
  init_number devSD(-5)
  init_number fut_devSD(-5)
  init_vector Fy(fyear,pyear,p8)
  sdreport_vector CA(fyear,pyear)
  matrix F(fyear,pyear,fage,lage)
  matrix Z(fyear,pyear,fage,lage)
  matrix S(fyear,pyear,fage,lage)
  number Ftemp
  vector vulner(fage,lage)
  vector retain(fage,lage)
  matrix N(fyear,pyear,fage,lage)
  vector VN(fage,lage)
  vector PR(fage,lage)

  matrix VulNumbers(fyear,pyear,fage,lage)
  matrix VulNos(fyear,pyear,fage,lage)
  matrix VulBio(fyear,pyear,fage,lage)
  matrix spawBio(fyear,pyear,fage,lage)
  matrix matBio(fyear,pyear,fage,lage)
  vector TspawBio(fyear,pyear)
  vector TvulBio(fyear,pyear)
  sdreport_vector SSB(fyear,pyear)
  sdreport_number VEstar
  number VMatBio
  matrix ECage(fyear,pyear,fage,lage)
  matrix C(fyear,pyear,fage,lage)
  vector ECat(fyear,pyear)
  sdreport_vector Rec(fyear,pyear)
  sdreport_vector Fy_eff(fyear,pyear)
  sdreport_vector propSSB(fyear,pyear)
  sdreport_number SRRa
  sdreport_number SRRb
  sdreport_number Estar
  vector sums(dyear,lyear)
  matrix PrVulAge(dyear,lyear,fage,lage)
  matrix ObsAge(dyear,lyear,fage,lage)
  matrix ObsPropAge(dyear,lyear,fage,lage)
  vector multi(dyear,lyear)
  vector N_age_comp(dyear,lyear)
  matrix PropAged(dyear,lyear,fage,lage)
  number crate
  number pen3
  number aa
  number bb
  number cc
  number dd
  number mm
  number nn
  number pp
  number rr
  number var
  number ppp
  number normdev
  number Npen
  vector phiy(dyear,lyear)
  number phi
  number temp1
  number temp2
  objective_function_value ff
 
PRELIMINARY_CALCS_SECTION
  if (checknum!=1234) { cout << "checknum " << checknum << endl;  exit(1); }; 
  
  sim=0;
  for (i=dyear; i<=lyear; i++) {
  if  (sum(SurAge(i))>0) { PropAged(i)=SurAge(i)/sum(SurAge(i));}
  else                   { PropAged(i)=0; }}
  for (i=dyear; i<=lyear; i++) { 
  if (sum(SurAge(i))>0) {N_age_comp(i)=sum(SurAge(i));}
  else                  {N_age_comp(i)=0.0; }}

  // calculate rge McAllaster and Ianell scalor phi
  phi=init_phi;
  phiy=0.0;
  //cout << "McAllister  phi="<< phi << endl;
  for (i=dyear; i<=lyear; i++) {
  temp1=0.0; temp2=0.0;
  for (j=fage; j<=lage; j++) { 
  temp1 +=(hatp_age(i,j)*(1.0-hatp_age(i,j))); 
  temp2 +=square(p_age(i,j)-hatp_age(i,j)); }
  if  (N_age_comp(i)>0) {phiy(i) +=temp1/temp2; }
  else                  {phiy(i)=0.0;}}
  phi= sum(phiy)/sum(N_age_comp);
  //cout << "McAllister  phi="<< phi << endl;
        
PROCEDURE_SECTION
  Npen=0.0;
  get_numbers_at_age();
  get_age_structure();
  estimate_effective_fmort();
  evaluate_the_objective_function();  

FUNCTION get_retention
  if (i<2001)  {for (j=fage; j<=lage; j++) {retain(j)= retain1(j); } }  
  if (i>=2001) {for (j=fage; j<=lage; j++) {retain(j)= retain2(j); } }   

FUNCTION get_numbers_at_age
 for (j=fage; j<=lage; j++) { vulner(j)=1.0/(1.0+exp(-1.0*SelSlope*(j-SelA50))); }
 
 // calculate relative fished numbers and the fished initial spawning biomass
 PR(fage)=1.0; Estar=0.0;
 for (j=fage;j<=lage-1;j++) {PR(j+1)=PR(j)*mfexp(-M-Fo*vulner(j)*retain1(j));}
 // plus group removed 23/06/15
 //PR(lage)=PR(lage)/(1-mfexp(-M-Fo*vulner(lage)*retain1(lage)));  //make plus group
 for (j=fage;j<=lage;j++)   {Estar+=PR(j)*WtF(j)*PropFem*SexMat(j); }

 // calcualte the initial unfished numbers and unfished SSB.   
 VN(fage)=Rstar;  VEstar=0.0;  
 for (j=fage;j<lage-1;j++) {VN(j+1)=VN(j)*mfexp(-M); }
 // plus group removed 23/06/15
 //VN(lage)=VN(lage)/(1.0-mfexp(-M)); //plus group
 for (j=fage;j<=lage;j++)   {VEstar+=VN(j)*WtF(j)*PropFem*SexMat(j); }

 // calculate the B-H stock recruitment parameters and first year recruitment
 SRRa=(VEstar/Rstar)*(1.0-((steep-0.2)/(0.8*steep)));
 SRRb=(steep-0.2)/(0.8*steep*Rstar);
 Rec(fyear)=(Estar-SRRa)/(SRRb*Estar);

 // set up the numbers of fish in the first year
 N(fyear,fage)= Rec(fyear);
 for (j=fage;j<lage;j++) {
   N(fyear,j+1)=N(fyear,j)*mfexp(-M-Fo*vulner(j)*retain1(j));
 }
 // plus group removed 23/06/15
 //N(fyear,lage)=N(fyear,lage)/(1.0-mfexp(-M-Fo*vulner(lage)*retain1(lage))); //make plus group

 // calculate the numbers of fish in each year and age 
 TvulBio=0.0;  SSB=0.0;  Npen=0.0;
 // FA=Fy;
  for (i=fyear; i<pyear; i++) {
    get_retention();
    estimated_catch();  
    // **************** calculate numbers of fish  ************************ ************
    for (j=fage;j<lage;j++) {
      N(i+1,j+1)=N(i,j)*mfexp(-M-Fy(i)*vulner(j)*retain(j));
    }
    // plus group removed 23/06/15
    //N(i+1,lage)=N(i+1,lage)/(1-mfexp(-M-Fy(i)*vulner(lage-1)*retain(lage))); //make lage a plus group
    for (j=fage;j<=lage;j++){
      SSB(i)+=PropFem*N(i,j)*SexMat(j)*WtF(j);
    }
    Rec(i)=SSB(i)/(SRRa+SRRb*SSB(i));
    //deviations - no bias adjustment for est years
    //bias adj on future years
    if (i<=lyear) {Rec(i)=Rec(i)*mfexp(dev(i))-0.5*devSD*devSD; }
    if (i>lyear)  {Rec(i)=Rec(i)*mfexp(fut_dev(i)-0.5*fut_devSD*fut_devSD); }
    N(i+1,fage)=Rec(i);
  } // end year

  for (i=pyear; i<=pyear; i++) {
    for (j=fage+1;j<=lage;j++) {
      SSB(i) +=N(i,j)*PropFem*SexMat(j)*WtF(j);
    }
    Rec(i)=SSB(i)/(SRRa+SRRb*SSB(i));
    //bias adj on future years
    Rec(i)=Rec(i)*mfexp(fut_dev(i)-0.5*fut_devSD*fut_devSD);
  }

  for (i=pyear; i<=pyear; i++) {
    get_retention();
    estimated_catch();
  }

  for (i=fyear; i<=pyear; i++) {
    propSSB(i)=SSB(i)/VEstar*100;
  }
 
FUNCTION estimated_catch
   CA(i)=0.0; 
   for (j=fage; j<=lage; j++) {
     F(i,j)=Fy(i)*vulner(j)*retain(j);
     Z(i,j)=F(i,j)+M;
     S(i,j)=1.0-exp(-Z(i,j));
     C(i,j)=N(i,j)*0.5*(WtF(j)+WtM(j))*S(i,j)*(F(i,j)/Z(i,j)); // catch  
     CA(i) +=C(i,j); 
   }
  
FUNCTION get_age_structure
 // Age structure from model;
  for (i=fyear; i<=pyear; i++) {
    for (j=fage;j<=lage;j++){
      VulNumbers(i,j)=N(i,j)*vulner(j); 
  }}
  for (i=dyear;i<=lyear;i++)  {
    //PrVulAge is matrix (year by age)
    PrVulAge(i)=VulNumbers(i)/sum(VulNumbers(i));  
    if (NageSample(i)>0) {
      ObsPropAge(i)=SurAge(i)/sum(SurAge(i));
    } else  ObsPropAge(i)=0.0;  //i.e. prevent dividing by zero for years with no sampling
    ObsAge(i)=SurAge(i);    
  }

  multi=0;
  for (i=dyear;i<=lyear;i++) { 
    for (j=fage; j<=lage; j++) {
      multi(i)+=-phi*SurAge(i,j)*log(PrVulAge(i,j));     
  }}

FUNCTION estimate_effective_fmort
  for (i=fyear; i<=pyear; i++) {
    Fy_eff(i)=0.0;
    if (i<2001)  {
      for (j=5;j<=11;j++){
        Fy_eff(i)+=Fy(i)*vulner(j)*retain1(j);
      }
    }else{
      for (j=5;j<=11;j++){
        Fy_eff(i)+=Fy(i)*vulner(j)*retain2(j);
      }
    }   
    Fy_eff(i)=Fy_eff(i)/(11.0-5.0+1.0);
  }
  
FUNCTION evaluate_the_objective_function
  //age comp likelihood
  aa=0;
  for (i=dyear;i<=lyear;i++) { 
    aa+=multi(i); // log likelihood
  }
  //depm likelihood
  // assume depm index (or biomass estimate) is lognormally distributed)
  cc=0;
  for (i=dyear;i<=lyear;i++) { 
    if (lnDEPM(i)>0) {
        cc+=0.5*square((log(Q*SSB(i))-lnDEPM(i))/lnsdDEPM(i)) + 0.5*log(square(lnsdDEPM(i))*2.0 * M_PI)+lnDEPM(i);
    }
  }
  //catch likelihood
  dd=0.0 ;
  for (i=fyear; i<=pyear; i++)   {
    dd +=0.5*square((OCat(i)-CA(i))/OCatSD(i)) + 0.5*log(square(OCatSD(i))*(2.0*M_PI));
  }
 
  //recruitment deviations 
  pp=0.0;   
  for (i=fyear; i<=lyear; i++)   {
    pp+=0.50*(square(dev(i)/devSD)); 
  }
  //future recruitment deviations
  ppp=0.0;
  for (i=lyear+1; i<=pyear; i++)   {
    ppp += 0.5*(square(fut_dev(i)/fut_devSD));
  }
  
  //constrain rec devs to sum to zero
  rr=0.0; 
  for (i=fyear; i<=lyear; i++)   {
    rr += dev(i);
  }
  rr=1000.0*rr*rr;

  //prior on natural Mortality
  mm=0.0;
  if (active(M)){
    mm=0.5*square((M-Mprior)/MpriorSD)-0.5*log(square(MpriorSD)*(2.0*M_PI));   
  }
  // prior on steepness
  nn=0.0;
  if (active(steep)){
    nn=0.5*square((steep-steepprior)/steeppriorSD) -0.5*log(square(steeppriorSD)*(2.0*M_PI));  
  }
 
  //random walk not used
  //var=Fysd*Fysd;
  //rr+=0.5*log(2.0*M_PI*var) + 0.5*square(Fy(fyear)-Fy0)/var;
  //for(i=fyear+1; i<=lyear; i++){
  //  rr +=0.5*log(2.0*M_PI*var) + 0.5*square(Fy(i)-Fy(i-1))/var;
  //}
  
  //objective function
  ff = 0.0;
  ff = WeightMULTI*aa+WeightDEPM*cc + dd; //Negative LL
  //add penalties
  ff += WeightDEV*pp + ppp + rr;
  ff += WeightM*mm + WeightH*nn + 1000.0*Npen;

  if (sd_phase()) {
   cout <<"OCat(2016)="<<OCat(2016)<< " phi "<<phi <<" aa.multi "<<aa<<" cc.DEPM "<<cc<< " dd.catch "<<dd<<endl;
   cout << "pp " << pp << " ppp " << ppp << " rr " << rr<< " ff.overalln "<<ff<< endl;
  }
   
  // ===== projection=======do 1,000,000 cycles and save every 400th (with first 500 discarded)
  
  if (mceval_phase()){
    sim=sim+1.0;
    if (sim>500){
      cout << steep << " " <<  M << " " <<  quota << " " << VEstar << " " << Rstar << " ";
      for (i=fyear;i<=pyear;i++) {cout <<  SSB(i) << " ";} 
      for (i=fyear;i<=pyear;i++) {cout <<  Fy(i) << " ";} 
      cout << endl;
    } 
  }


REPORT_SECTION

  report<< "WtF" << endl;
  report<< WtF << endl;
  report<< "WtM" << endl;
  report<< WtM << endl;
  report<< "vulner " << endl;
  report<< vulner << endl;
  report << "Fy";
  for (i=fyear; i<=pyear; i++)  { report << " " << i ; }
  report <<  endl;
  for (i=fyear; i<=pyear; i++)  { report << " " << Fy(i) ; }
  report <<  endl;
  report << "Observed Catch " << endl;
  for (i=fyear; i<=pyear; i++)  { report << " " << i ; } report <<  endl;
  report << " " << OCat << endl;
  report << "TvulBio"<< endl;
  for (i=fyear; i<=pyear; i++)  { report << " " << i ; }  report <<  endl;
  report <<  " " ;
  for (i=fyear; i<=pyear; i++)  { report << TvulBio(i) << " "; }
  report <<  endl;
  report << "SSB" << endl;
  for (i=fyear; i<=pyear; i++)  { report << " " << i ; }  report <<  endl;
  report <<  " " ;
  for (i=fyear; i<=pyear; i++)  { report << SSB(i) << " "; }
  report <<  endl;
  report << "observed then predicted ages " << endl; 
  for (j=fage; j<=lage; j++){  report << " " << age(j);}
  report << endl;
  for (i=dyear; i<=lyear; i++){
    for (j=fage; j<=lage; j++){  report << " " << ObsPropAge(i,j);}
    report << endl;
    for (j=fage; j<=lage; j++){  report << " " << PrVulAge(i,j);}
    report << endl;
  } 
  report << "CA" << endl;
  report << CA << endl;
  report << "DEPM" << endl;
  report << DEPM << endl;
  report << lowDEPM << endl;
  report << uppDEPM << endl;
  report << "sdDEPM" << endl;
  report << sdDEPM << endl;
  report << "SSage   SSdepm  SScatch PenM  PenH  recdev futrecdev  recconstrain" << endl;
  report << aa << " "  << cc << " " << dd << " " << mm << " " << nn << " " << pp << " "<< ppp << " " << rr << " " << endl;
  report << "Wt_age*SSaa Wt_depm*SSdepm   WeightM*PenM   WeightH*PenH  ff" << endl;
  report << WeightMULTI*aa << " "  << WeightDEPM*cc << " " << WeightM*mm << " " << WeightH*nn <<  " " << ff << endl;
  report << "Estar  Rstar " << Estar << " " << Rstar << endl;
  report << "SSB" << endl;
  report << SSB << endl;
  report << "Rec" << endl;
  report << Rec << endl;
  report << "SRRa SRRb " << SRRa << " " << SRRb << endl;
  report << "BandH_recruit" << endl;
  for (i=fyear; i<=pyear; i++){  report <<" "<< SSB(i)/(SRRa+SRRb*SSB(i));}
  report << endl;
  report << "ObsPropAge" << endl;
  report << ObsPropAge << endl;
  report << "PrVulAge" << endl;
  report <<PrVulAge << endl;
  report << "Retain1" << endl;
  report << retain1 << endl;
  report << "Retain2" << endl;
  report << retain2 << endl;
  report << "NageSample" << endl;
  report << NageSample << endl;
   
  ofstream out3("prop_age_composition.dat");
  out3 <<  ObsPropAge<<endl; 
  out3.close();
 
  ofstream out4("pred_prop_age_composition.dat");
  out4 << PrVulAge <<endl; 
 
  ofstream out5("value_of_phi.dat");
  out5 << phi << endl;
  out5.close();
  
  // RUNTIME_SECTION
  // convergence criteria 0.00001
TOP_OF_MAIN_SECTION
 arrmblsize=800000;
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(300);
                                           

