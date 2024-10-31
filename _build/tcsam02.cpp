#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
    #include <math.h>
    #include <time.h>
    #include <limits>
    #include <admodel.h>
    #include "TCSAM.hpp"
    #ifdef PRINT2B1
        #undef PRINT2B1
    #endif
    #define PRINT2B1(o) std::cout<<(o)<<std::endl; rpt::echo<<(o)<<std::endl;
    #ifdef PRINT2B2
        #undef PRINT2B2
    #endif
    #define PRINT2B2(t,o) std::cout<<(t)<<(o)<<std::endl; rpt::echo<<(t)<<(o)<<std::endl;
    adstring model  = tcsam::MODEL;
    adstring modVer = tcsam::VERSION; 
    time_t start,finish;
    //model objects
    ModelConfiguration*  ptrMC; //ptr to model configuration object
    ModelParametersInfo* ptrMPI;//ptr to model parameters info object
    ModelOptions*        ptrMOs;//ptr to model options object
    ModelDatasets*       ptrMDS;//ptr to model datasets object
    ModelDatasets*       ptrSimMDS;//ptr to simulated model datasets object
    OFLResults*          ptrOFLResults;//ptr to OFL results object for MCMC calculations
    //MSE objects
    MSE_OpModInfo* ptrOMI=0; //ptr to MSE OpModInfo object
    //population dynamics objects
    int runAlt = 0;
    PopDyInfo*    pPDI=0;//  population dynamics info
    CatchInfo*    pCDI=0;//  catch info
    PopProjector* pPPr=0;//  population projector
    //dimensions for R output
    adstring yDms;
    adstring xDms;
    adstring mDms;
    adstring sDms;
    adstring fDms;
    adstring vDms;
    adstring ypDms;
    adstring zbDms;
    adstring zpDms;
    adstring zcDms;
    //file streams and filenames
    long ctrMCMC = 0;    //counter for mcmc output
    std::ofstream mcmc;  //stream for mcmc output
    //filenames
    adstring fnMCMC = "tcsam02.MCMC.R";
    adstring fnConfigFile;//configuration file
    adstring fnPin;       //pin file
    //runtime flags (0=false)
    int jitter     = 0;//use jittering for initial parameter values
    int resample   = 0;//use resampling for initial parameter values
    int mcevalOn   = 0;//flag indicating model is being run in mceval phase
    int mseMode    = 0;//flag indicating model is being run in an MSE
    int mseOpModMode  = 0;//flag indicating model is being run in an MSE in operating model mode
    int mseEstModMode = 0;//flag indicating model is being run in an MSE in estimation model mode 
    int usePin     = 0;//flag to initialize parameter values using a pin file
    int doRetro    = 0;//flag to facilitate a retrospective model run
    int fitSimData = 0;//flag to fit model to simulated data calculated in the PRELIMINARY_CALCs section
    int doOFL      = 0;///<flag (0/1) to do OFL calculations
    int doTAC      = 0;///<calculate TAC using harvest control rule indicated by value of doTAC
    int doDynB0    = 0;//flag to run dynamic B0 calculations after final phase
    int yRetro = 0; //number of years to decrement for retrospective model run
    int iSeed = -1; //default random number generator seed
    random_number_generator rng(iSeed);//random number generator
    int iSimDataSeed = 0;
    random_number_generator rngSimData(-1);//random number generator for data simulation
    //debug flags
    int debugModelConfig     = 0;
    int debugModelDatasets   = 0;
    int debugModelParamsInfo = 0;
    int debugModelParams     = 0;
    int debugModelOptions    = 0;
    int debugDATA_SECTION    = 0;
    int debugPARAMS_SECTION  = 0;
    int debugPRELIM_CALCS    = 0;
    int debugPROC_SECTION    = 0;
    int debugREPORT_SECTION  = 0;
    int ctrDebugParams = 0;//PROCEDURE_SECTION call counter value to start debugging output
    int showActiveParams = 0;    
    int debugRunModel    = 0;    
    int debugObjFun      = 0;
    int debugOFL         = 0;
    int debugMCMC        = 0;
    //note: consider using std::bitset to implement debug functionality
    int dbgCalcProcs = 10;
    int dbgObjFun = 20;
    int dbgNLLs   = 25;
    int dbgPriors = tcsam::dbgPriors;
    int dbgPopDy  = 70;
    int dbgApply  = 80;
    int dbgDevs   = 90;
    int dbgAll    = tcsam::dbgAll;
    int phsItsRewgt  = 1000;//min phase to calculate effective weights for size comps
    int maxItsRewgt = 0;    //maximum number of terations for re-weighting
    int numItsRewgt = 0;    //number of re-weighting iterations completed
    int nSXs    = tcsam::nSXs;
    int MALE    = tcsam::MALE;
    int FEMALE  = tcsam::FEMALE;
    int ALL_SXs = tcsam::ALL_SXs;
    int nMSs     = tcsam::nMSs;
    int IMMATURE = tcsam::IMMATURE;
    int MATURE   = tcsam::MATURE;
    int ALL_MSs  = tcsam::ALL_MSs;
    int nSCs      = tcsam::nSCs;
    int NEW_SHELL = tcsam::NEW_SHELL;
    int OLD_SHELL = tcsam::OLD_SHELL;
    int ALL_SCs   = tcsam::ALL_SCs;
    double smlVal = 0.00001;//small value to keep things > 0
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <tcsam02.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
    {
        adstring msg = "#Starting "+model+" (ver "+modVer+") Code";
        PRINT2B1(msg)
        PRINT2B1("#Starting DATA_SECTION")
    }
    int on = 0;
    int flg = 0;
    PRINT2B1("#------Reading command line options---------")
    //configFile
    fnConfigFile = "tcsam02.ModelConfig.dat";//default model config filename
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-configFile"))>-1) {
        fnConfigFile = ad_comm::argv[on+1];
        rpt::echo<<"#config file changed to '"<<fnConfigFile<<"'"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg=1;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-pin"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        if (std::fstream(fnPin)){
            ad_comm::change_pinfile_name(fnPin);
            rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        } else {
            rpt::echo<<"#Initial parameter values from pin file: tcsam02.pin"<<endl;
        }
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-binp"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-ainp"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //mceval phase is on
    if ((option_match(ad_comm::argc,ad_comm::argv,"-mceval")>-1)){
        mcevalOn = 1;
        rpt::echo<<"#mceval is ON."<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //parameter input file for running NUTS MCMC 
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-mcpin"))>-1) {
        usePin = 2;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values for running NUTS MCMC from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //runAlt
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-runAlt"))>-1) {
        runAlt=1;
        rpt::echo<<"#run alternative population dynamics functions"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //calcOFL
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-calcOFL"))>-1) {
        doOFL=1;
        rpt::echo<<"#OFL calculations turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    // Add a section for calcOFL for Operating model?? 
    //if ((on=option_match(ad_comm::argc,ad_comm::argv,"-calcOFL_OpMod"))>-1) {
      //  doOFL=1;
        //rpt::echo<<"#OFL Op Mod calculations turned ON"<<endl;
        //rpt::echo<<"#-------------------------------------------"<<endl;
        //flg = 1;
    }
    //calcTAC
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-calcTAC"))>-1) {
        doTAC = atoi(ad_comm::argv[on+1]);
        doOFL = 1;
        rpt::echo<<"#Calculating TAC using harvest control rule "<<doTAC<<endl;
        rpt::echo<<"#OFL calculations also turned ON "<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //calcTAC_OpMod
    //if ((on=option_match(ad_comm::argc,ad_comm::argv,"-calcTAC_OpMod"))>-1) {
      //  doTAC = atoi(ad_comm::argv[on+1]);
        //doOFL = 1;
        //rpt::echo<<"#Calculating TAC for Op Mod using harvest control rule "<<doTAC<<endl;
        //rpt::echo<<"#OFL calculations also turned ON "<<endl;
        //rpt::echo<<"#-------------------------------------------"<<endl;
        //flg = 1;
    //}
  
    
    //run dynamic B0 calculations after final phase
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-calcDynB0"))>-1) {
        doDynB0 = 1;
        rpt::echo<<"#Running dynamic B0 calculations after final phase."<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //doRetro
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-doRetro"))>-1) {
        doRetro=1;
        cout<<"#doRetro turned ON"<<endl;
        rpt::echo<<"#doRetro turned ON"<<endl;
        if (on+1<argc) {
            yRetro=atoi(ad_comm::argv[on+1]);
            cout<<"#Retrospective model run using yRetro = "<<yRetro<<endl;
            rpt::echo<<"#Retrospective model run using yRetro = "<<yRetro<<endl;
        }
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //fitSimData
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-fitSimData"))>-1) {
        fitSimData=1;
        if (on+1<argc) {
            iSimDataSeed=atoi(ad_comm::argv[on+1]);
        } else {
            cout<<"-------------------------------------------"<<endl;
            cout<<"Enter random number seed (0 -> deterministic) for data simulation: ";
            cin>>iSimDataSeed;
        }
        if (iSimDataSeed) rng.reinitialize(iSimDataSeed);
        rpt::echo<<"#Simulating data to fit using "<<iSimDataSeed<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //ctrDebugParams
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-ctrDebugParams"))>-1) {
        ctrDebugParams=atoi(ad_comm::argv[on+1]);
        rpt::echo<<"Starting debugParams at counter "<<ctrDebugParams<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //jitter
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-jitter"))>-1) {
        jitter=1;
        iSeed=(long)start;
        if ((on=option_match(ad_comm::argc,ad_comm::argv,"-iSeed"))>-1) {
            if (on+1<argc) {
                iSeed=atoi(ad_comm::argv[on+1]);
            }
        } 
        rng.reinitialize(iSeed);
        rpt::echo<<"#Jittering for initial parameter values turned ON "<<endl;
        rpt::echo<<iSeed<<"  #iSeed"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        ofstream fs("jitterInfo.dat");
        fs<<"seed = "<<iSeed<<endl;
        fs.close();
        flg = 1;
    }
    //resample
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-resample"))>-1) {
        resample=1;
        rpt::echo<<"#Resampling for initial parameter values turned ON "<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //mseMode
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-mseMode"))>-1) {
        mseMode=1;
        adstring type = ad_comm::argv[on+1];
        if (type=="mseOpModMode"){
            mseOpModMode = 1;
            rpt::echo<<"#MSE operating model mode turned ON"<<endl;
            iSeed=(long)start;
            if ((on=option_match(ad_comm::argc,ad_comm::argv,"-iSeed"))>-1) {
                if (on+1<argc) {
                    iSeed=atoi(ad_comm::argv[on+1]);
                }
            } 
            rng.reinitialize(iSeed);
            rpt::echo<<tb<<iSeed<<"  #iSeed used for random recruitment"<<endl;
        } else if (type=="mseEstModMode"){
            mseEstModMode = 1;
            rpt::echo<<"#MSE estimation model mode turned ON"<<endl;
        }
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelConfig
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelConfig"))>-1) {
        debugModelConfig=1;
        rpt::echo<<"#debugModelConfig turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelParams
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelParams"))>-1) {
        debugModelParams=1;
        cout<<"#debugModelParams turned ON"<<endl;
        rpt::echo<<"#debugModelParams turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugDATA_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugDATA_SECTION"))>-1) {
        debugDATA_SECTION=1;
        rpt::echo<<"#debugDATA_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPARAMS_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPARAMS_SECTION"))>-1) {
        debugPARAMS_SECTION=1;
        rpt::echo<<"#debugPARAMS_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPRELIM_CALCS
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPRELIM_CALCS"))>-1) {
        debugPRELIM_CALCS=1;
        rpt::echo<<"debugPRELIM_CALCS turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPROC_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPROC_SECTION"))>-1) {
        debugPROC_SECTION=1;
        rpt::echo<<"#debugPROC_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugREPORT_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugREPORT_SECTION"))>-1) {
        debugREPORT_SECTION=1;
        rpt::echo<<"#debugREPORT_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugOFL
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugOFL"))>-1) {
        debugOFL=1;
        rpt::echo<<"#debugOFL turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelDatasets
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelDatasets"))>-1) {
        debugModelDatasets=1;
        rpt::echo<<"#debugModelDatasets turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelParamsInfo
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelParamsInfo"))>-1) {
        debugModelParamsInfo=1;
        rpt::echo<<"#debugModelParamsInfo turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelOptions
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelOptions"))>-1) {
        debugModelOptions=1;
        rpt::echo<<"#debugModelOptions turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugRunModel
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugRunModel")>-1) {
        debugRunModel=1;
        rpt::echo<<"#debugRunModel turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugObjFun
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugObjFun")>-1) {
        debugObjFun=1;
        rpt::echo<<"#debugObjFun turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debuMCMC
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugMCMC")>-1) {
        debugMCMC=1;
        rpt::echo<<"#debugMCMC turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //showActiveParams
    if (option_match(ad_comm::argc,ad_comm::argv,"-showActiveParams")>-1) {
        showActiveParams=1;
        rpt::echo<<"#showActiveParams turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    if (mseOpModMode&&mseEstModMode){
        PRINT2B1("");
        PRINT2B1("--------ERROR!-----------")
        PRINT2B1("mseOpModMode and mseEstModMode cannot both be 'on'.")
        PRINT2B1("Terminating model run!!")
        PRINT2B1("--------ERROR!-----------")
        PRINT2B1("")
        exit(1);
    }
    PRINT2B1("#-----------------------------------")
    PRINT2B1("#-----------------------------------")
    PRINT2B2("#--Reading configuration file ",fnConfigFile)
    ad_comm::change_datafile_name(fnConfigFile);
    ptrMC = new ModelConfiguration();
    ptrMC->read(*(ad_comm::global_datafile));
    PRINT2B1("#--Finished reading configuration file")
    
    mnYr   = ptrMC->mnYr;
    mxYr   = ptrMC->mxYr;
    if (doRetro){mxYr = mxYr-yRetro; ptrMC->setMaxModelYear(mxYr);}
    if (jitter)   {ptrMC->jitter   = 1;}
    if (resample) {ptrMC->resample = 1;}
    
    rpt::echo<<"#------------------ModelConfiguration-----------------"<<endl;
    rpt::echo<<(*ptrMC);
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#----finished model configuration---"<<endl;
    rpt::echo<<"#-----------------------------------"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelConfiguration-----------------"<<endl;
        cout<<(*ptrMC);
        cout<<"#-----------------------------------"<<endl;
        cout<<"#----finished model configuration---"<<endl;
        cout<<"#-----------------------------------"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
    
    mxYrp1 = mxYr+1;
    nFsh   = ptrMC->nFsh;
    nSrv   = ptrMC->nSrv;
    nZBs   = ptrMC->nZBs;
  zBs.allocate(1,nZBs);
zBs  = ptrMC->zMidPts;
    PRINT2B1("#-----------------------------------")
    PRINT2B2("#Reading parameters info file ",ptrMC->fnMPI)
    if (debugModelParamsInfo) ModelParametersInfo::debug=1;
    ptrMPI = new ModelParametersInfo(*ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMPI);
    ptrMPI->read(*(ad_comm::global_datafile));
    if (debugModelParamsInfo) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelParamsInfo;
        if (debugModelParamsInfo<0) exit(1);
        ModelParametersInfo::debug=debugModelParamsInfo;
    }
    PRINT2B1("#----finished reading model parameters info---")
    PRINT2B1("#-----------------------------------")
    PRINT2B2("#--Reading datasets file ",ptrMC->fnMDS);
    if (debugModelDatasets) {
        BioData::debug=1;
        FleetData::debug=1;
    }
    ptrMDS = new ModelDatasets(ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMDS);
    ptrMDS->read(*(ad_comm::global_datafile));
    if (debugModelDatasets) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelDatasets;
        if (debugModelDatasets<0) exit(1);
        ModelDatasets::debug=debugModelDatasets;
        BioData::debug=debugModelDatasets;
        FleetData::debug=debugModelDatasets;
    }
    PRINT2B1("#----finished model datasets---")
    if (debugDATA_SECTION){
        cout<<"#------------------ModelDatasets-----------------"<<endl;
        cout<<(*ptrMDS);
        cout<<"#----finished model datasets---"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading datasets file again to create SimMDS object '"<<ptrMC->fnMDS<<"'"<<endl;
    ptrSimMDS = new ModelDatasets(ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMDS);
    ptrSimMDS->read(*(ad_comm::global_datafile));
    PRINT2B1("#-----------------------------------")
    PRINT2B2("#--Reading model options file ",ptrMC->fnMOs)
    if (debugModelOptions) ModelOptions::debug=1;
    ptrMOs = new ModelOptions(*ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMOs);
    ptrMOs->read(*(ad_comm::global_datafile));
    if (debugModelOptions) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelOptions;
        if (debugModelOptions<0) exit(1);
        ModelOptions::debug=debugModelOptions;
    }
    doTAC=ptrMOs->HCR;
    PRINT2B1("#--finished reading model options file---")
    rpt::echo<<"#------------------ModelOptions-----------------"<<endl;
    rpt::echo<<(*ptrMOs);
    rpt::echo<<"#----finished model options---"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelOptions-----------------"<<endl;
        cout<<(*ptrMOs);
        cout<<"#----finished model options---"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
    //check commandline options for iterative reweighting overrides
    //min phase in which to calculate effective weights for size compositions
    phsItsRewgt = ptrMOs->phsIterativeReweighting;
    maxItsRewgt  = ptrMOs->maxIterations;
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-phsIterativeReweighing"))>-1) {
        if (on+1<argc) {
            phsItsRewgt=atoi(ad_comm::argv[on+1]);
        }
        flg = 1;
    }
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-maxIterationsReweighting"))>-1) {
        if (on+1<argc) {
            phsItsRewgt=atoi(ad_comm::argv[on+1]);
        }
        flg = 1;
    }
    rpt::echo<<phsItsRewgt<<tb<<"#phsIterativeReweighing"<<endl;
    rpt::echo<<maxItsRewgt<<tb<<"#maxIterationsReweighting"<<endl;
    rpt::echo<<"#-------------------------------------------"<<endl;
  mapD2MFsh.allocate(1,nFsh);
  mapM2DFsh.allocate(1,nFsh);
    {
     int idx;
     for (int f=1;f<=nFsh;f++){
         idx = wts::which(ptrMDS->ppFsh[f-1]->name,ptrMC->lblsFsh);
         if (idx<1){
             cout<<"\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
             cout<<"Error specifying fishery names and labels in data file and config file."<<endl;
             cout<<"Incorrect fishery name in data file is '"<<ptrMDS->ppFsh[f-1]<<"'"<<endl;
             cout<<"Please fix names in files!!"<<endl;
             exit(-1);
         }
         mapD2MFsh(f)   = idx;//map from fishery data object f to model fishery idx
         mapM2DFsh(idx) = f;  //map from model fishery idx to fishery data object f
     }
     PRINT2B2("model fisheries map to fishery data objects: ",mapM2DFsh)
    }
  mapD2MSrv.allocate(1,nSrv);
  mapM2DSrv.allocate(1,nSrv);
    {
     int idx;
     for (int v=1;v<=nSrv;v++){
         idx = wts::which(ptrMDS->ppSrv[v-1]->name,ptrMC->lblsSrv);
         if (idx<1){
             cout<<"\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
             cout<<"Error specifying survey names and labels in data file and config file."<<endl;
             cout<<"Incorrect survey name in data file is '"<<ptrMDS->ppSrv[v-1]<<"'"<<endl;
             cout<<"Please fix names in files!!"<<endl;
             exit(-1);
         }
         mapD2MSrv(v)   = idx;//map from survey data object v to model survey idx
         mapM2DSrv(idx) = v;  //map from model survey idx to survey data object v
     }
     PRINT2B2("model surveys map to survey data objects: ",mapM2DSrv)
    }
  mapD2MChd.allocate(1,ptrMDS->nCHD);
    {
     int idx;
     for (int v=1;v<=ptrMDS->nCHD;v++){
         idx = wts::which(ptrMDS->ppCHD[v-1]->survey,ptrMC->lblsSrv);
         if (idx<1){
             cout<<"\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
             cout<<"Error specifying survey names and labels in CH data file and config file."<<endl;
             cout<<"Incorrect survey name in CH data file is '"<<ptrMDS->ppCHD[v-1]<<"'"<<endl;
             cout<<"Please fix names in files!!"<<endl;
             exit(-1);
         }
         mapD2MChd(v)   = idx;//map from chela height dataset object v to model survey index
         ptrMDS->ppCHD[v-1]->calcSizeBinIndices(ptrMC->zCutPts);
     }
     PRINT2B2("chela height datasets map to model surveys: ",mapD2MChd)
    }
  phsLnR.allocate();
  lbLnR.allocate();
  ubLnR.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pLnR,npLnR,lbLnR,ubLnR,phsLnR,rpt::echo);
  phsRCV.allocate();
  lbRCV.allocate();
  ubRCV.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pRCV,npRCV,lbRCV,ubRCV,phsRCV,rpt::echo);
  phsRX.allocate();
  lbRX.allocate();
  ubRX.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pRX,npRX,lbRX,ubRX,phsRX,rpt::echo);
  phsRa.allocate();
  lbRa.allocate();
  ubRa.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pRa,npRa,lbRa,ubRa,phsRa,rpt::echo);
  phsRb.allocate();
  lbRb.allocate();
  ubRb.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pRb,npRb,lbRb,ubRb,phsRb,rpt::echo);
  mniDevsLnR.allocate();
  mxiDevsLnR.allocate();
  idxsDevsLnR.allocate();
  lbDevsLnR.allocate();
  ubDevsLnR.allocate();
  phsDevsLnR.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pDevsLnR,npDevsLnR,mniDevsLnR,mxiDevsLnR,idxsDevsLnR,lbDevsLnR,ubDevsLnR,phsDevsLnR,rpt::echo);
    if (mseOpModMode) {
        phsLnR = -1;
        phsRCV = -1;
        phsRX  = -1;
        phsRa  = -1;
        phsRb  = -1;
        phsDevsLnR = -1;
    }
  phsM.allocate();
  lbM.allocate();
  ubM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pM,npM,lbM,ubM,phsM,rpt::echo);
  phsDM1.allocate();
  lbDM1.allocate();
  ubDM1.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pDM1,npDM1,lbDM1,ubDM1,phsDM1,rpt::echo);
  phsDM2.allocate();
  lbDM2.allocate();
  ubDM2.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pDM2,npDM2,lbDM2,ubDM2,phsDM2,rpt::echo);
  phsDM3.allocate();
  lbDM3.allocate();
  ubDM3.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pDM3,npDM3,lbDM3,ubDM3,phsDM3,rpt::echo);
  phsDM4.allocate();
  lbDM4.allocate();
  ubDM4.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pDM4,npDM4,lbDM4,ubDM4,phsDM4,rpt::echo);
zMref = ptrMPI->ptrNM->zRef;
    if (mseOpModMode) {
        phsM = -1;
        phsDM1 = -1;
        phsDM2 = -1;
        phsDM3 = -1;
        phsDM4 = -1;
    }
  mniLgtPrMat.allocate();
  mxiLgtPrMat.allocate();
  idxsLgtPrMat.allocate();
  lbLgtPrMat.allocate();
  ubLgtPrMat.allocate();
  phsLgtPrMat.allocate();
tcsam::setParameterInfo(ptrMPI->ptrM2M->pvLgtPrM2M,npLgtPrMat,mniLgtPrMat,mxiLgtPrMat,idxsLgtPrMat,lbLgtPrMat,ubLgtPrMat,phsLgtPrMat,rpt::echo);
    if (mseOpModMode) {
        phsLgtPrMat = -1;
    }
  phsGrA.allocate();
  lbGrA.allocate();
  ubGrA.allocate();
tcsam::setParameterInfo(ptrMPI->ptrGrw->pGrA,npGrA,lbGrA,ubGrA,phsGrA,rpt::echo);
  phsGrB.allocate();
  lbGrB.allocate();
  ubGrB.allocate();
tcsam::setParameterInfo(ptrMPI->ptrGrw->pGrB,npGrB,lbGrB,ubGrB,phsGrB,rpt::echo);
  phsGrBeta.allocate();
  lbGrBeta.allocate();
  ubGrBeta.allocate();
tcsam::setParameterInfo(ptrMPI->ptrGrw->pGrBeta,npGrBeta,lbGrBeta,ubGrBeta,phsGrBeta,rpt::echo);
    if (mseOpModMode) {
        phsGrA = -1;
        phsGrB = -1;
        phsGrBeta = -1;
    }
  phsS1.allocate();
  lbS1.allocate();
  ubS1.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS1,npS1,lbS1,ubS1,phsS1,rpt::echo);
  phsS2.allocate();
  lbS2.allocate();
  ubS2.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS2,npS2,lbS2,ubS2,phsS2,rpt::echo);
  phsS3.allocate();
  lbS3.allocate();
  ubS3.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS3,npS3,lbS3,ubS3,phsS3,rpt::echo);
  phsS4.allocate();
  lbS4.allocate();
  ubS4.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS4,npS4,lbS4,ubS4,phsS4,rpt::echo);
  phsS5.allocate();
  lbS5.allocate();
  ubS5.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS5,npS5,lbS5,ubS5,phsS5,rpt::echo);
  phsS6.allocate();
  lbS6.allocate();
  ubS6.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS6,npS6,lbS6,ubS6,phsS6,rpt::echo);
  mniDevsS1.allocate();
  mxiDevsS1.allocate();
  idxsDevsS1.allocate();
  lbDevsS1.allocate();
  ubDevsS1.allocate();
  phsDevsS1.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS1,npDevsS1,mniDevsS1,mxiDevsS1,idxsDevsS1,lbDevsS1,ubDevsS1,phsDevsS1,rpt::echo);
  mniDevsS2.allocate();
  mxiDevsS2.allocate();
  idxsDevsS2.allocate();
  lbDevsS2.allocate();
  ubDevsS2.allocate();
  phsDevsS2.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS2,npDevsS2,mniDevsS2,mxiDevsS2,idxsDevsS2,lbDevsS2,ubDevsS2,phsDevsS2,rpt::echo);
  mniDevsS3.allocate();
  mxiDevsS3.allocate();
  idxsDevsS3.allocate();
  lbDevsS3.allocate();
  ubDevsS3.allocate();
  phsDevsS3.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS3,npDevsS3,mniDevsS3,mxiDevsS3,idxsDevsS3,lbDevsS3,ubDevsS3,phsDevsS3,rpt::echo);
  mniDevsS4.allocate();
  mxiDevsS4.allocate();
  idxsDevsS4.allocate();
  lbDevsS4.allocate();
  ubDevsS4.allocate();
  phsDevsS4.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS4,npDevsS4,mniDevsS4,mxiDevsS4,idxsDevsS4,lbDevsS4,ubDevsS4,phsDevsS4,rpt::echo);
  mniDevsS5.allocate();
  mxiDevsS5.allocate();
  idxsDevsS5.allocate();
  lbDevsS5.allocate();
  ubDevsS5.allocate();
  phsDevsS5.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS5,npDevsS5,mniDevsS5,mxiDevsS5,idxsDevsS5,lbDevsS5,ubDevsS5,phsDevsS5,rpt::echo);
  mniDevsS6.allocate();
  mxiDevsS6.allocate();
  idxsDevsS6.allocate();
  lbDevsS6.allocate();
  ubDevsS6.allocate();
  phsDevsS6.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS6,npDevsS6,mniDevsS6,mxiDevsS6,idxsDevsS6,lbDevsS6,ubDevsS6,phsDevsS6,rpt::echo);
  mniNPSel.allocate();
  mxiNPSel.allocate();
  idxsNPSel.allocate();
  lbNPSel.allocate();
  ubNPSel.allocate();
  phsNPSel.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pvNPSel,npNPSel,mniNPSel,mxiNPSel,idxsNPSel,lbNPSel,ubNPSel,phsNPSel,rpt::echo);
    if (mseOpModMode) {
        phsS1 = -1;
        phsS2 = -1;
        phsS3 = -1;
        phsS4 = -1;
        phsS5 = -1;
        phsS6 = -1;
        phsDevsS1 = -1;
        phsDevsS2 = -1;
        phsDevsS3 = -1;
        phsDevsS4 = -1;
        phsDevsS5 = -1;
        phsDevsS6 = -1;
        phsNPSel  = -1;
    }
  phsHM.allocate();
  lbHM.allocate();
  ubHM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pHM,npHM,lbHM,ubHM,phsHM,rpt::echo);
  phsLnC.allocate();
  lbLnC.allocate();
  ubLnC.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnC,npLnC,lbLnC,ubLnC,phsLnC,rpt::echo);
  phsDC1.allocate();
  lbDC1.allocate();
  ubDC1.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pDC1,npDC1,lbDC1,ubDC1,phsDC1,rpt::echo);
  phsDC2.allocate();
  lbDC2.allocate();
  ubDC2.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pDC2,npDC2,lbDC2,ubDC2,phsDC2,rpt::echo);
  phsDC3.allocate();
  lbDC3.allocate();
  ubDC3.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pDC3,npDC3,lbDC3,ubDC3,phsDC3,rpt::echo);
  phsDC4.allocate();
  lbDC4.allocate();
  ubDC4.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pDC4,npDC4,lbDC4,ubDC4,phsDC4,rpt::echo);
  mniDevsLnC.allocate();
  mxiDevsLnC.allocate();
  idxsDevsLnC.allocate();
  lbDevsLnC.allocate();
  ubDevsLnC.allocate();
  phsDevsLnC.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pDevsLnC,npDevsLnC,mniDevsLnC,mxiDevsLnC,idxsDevsLnC,lbDevsLnC,ubDevsLnC,phsDevsLnC,rpt::echo);
  phsLnEffX.allocate();
  lbLnEffX.allocate();
  ubLnEffX.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnEffX,npLnEffX,lbLnEffX,ubLnEffX,phsLnEffX,rpt::echo);
  phsLgtRet.allocate();
  lbLgtRet.allocate();
  ubLgtRet.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLgtRet,npLgtRet,lbLgtRet,ubLgtRet,phsLgtRet,rpt::echo);
    if (mseOpModMode) {
        phsHM = -1;
        phsLnC = -1;
        phsDC1 = -1;
        phsDC2 = -1;
        phsDC3 = -1;
        phsDC4 = -1;
        phsDevsLnC = -1;
        phsLnEffX = -1;
        phsLgtRet = -1;
    }
  phsQ.allocate();
  lbQ.allocate();
  ubQ.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pQ,npQ,lbQ,ubQ,phsQ,rpt::echo);
  phsDQ1.allocate();
  lbDQ1.allocate();
  ubDQ1.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pDQ1,npDQ1,lbDQ1,ubDQ1,phsDQ1,rpt::echo);
  phsDQ2.allocate();
  lbDQ2.allocate();
  ubDQ2.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pDQ2,npDQ2,lbDQ2,ubDQ2,phsDQ2,rpt::echo);
  phsDQ3.allocate();
  lbDQ3.allocate();
  ubDQ3.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pDQ3,npDQ3,lbDQ3,ubDQ3,phsDQ3,rpt::echo);
  phsDQ4.allocate();
  lbDQ4.allocate();
  ubDQ4.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pDQ4,npDQ4,lbDQ4,ubDQ4,phsDQ4,rpt::echo);
    if (mseOpModMode) {
        phsQ = -1;
        phsDQ1 = -1;
        phsDQ2 = -1;
        phsDQ3 = -1;
        phsDQ4 = -1;
    }
  phsMSE_LnC.allocate();
  lbMSE_LnC.allocate();
  ubMSE_LnC.allocate();
tcsam::setParameterInfo(ptrMPI->ptrMSE->pMSE_LnC,npMSE_LnC,lbMSE_LnC,ubMSE_LnC,phsMSE_LnC,rpt::echo);
if (!mseOpModMode) {phsMSE_LnC = -1;}
if ( mseOpModMode) {phsMSE_LnC =  1;}
  dtF_y.allocate(mnYr,mxYr);
dtF_y = ptrMDS->ptrBio->fshTiming_y(mnYr,mxYr);
  dtM_y.allocate(mnYr,mxYr);
dtM_y = ptrMDS->ptrBio->fshTiming_y(mnYr,mxYr);
nEASs = ptrMOs->ptrEffXtrapScenarios->ptrEffAvgScenarios->nAvgs;
  yrsAvgEff_ny.allocate(1,nEASs,mnYr,mxYr);
  eff_ny.allocate(1,nEASs,mnYr,mxYr);
  avgEff_n.allocate(1,nEASs);
nCRASs = ptrMOs->ptrEffXtrapScenarios->ptrCapRateAvgScenarios->nAvgs;
  obsEff_nxmsy.allocate(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr);
npcRec = ptrMPI->ptrRec->nPCs;
npcNM = ptrMPI->ptrNM->nPCs;
npcM2M = ptrMPI->ptrM2M->nPCs;
npcGrw = ptrMPI->ptrGrw->nPCs;
npcSel = ptrMPI->ptrSel->nPCs;
npcFsh = ptrMPI->ptrFsh->nPCs;
npcSrv = ptrMPI->ptrSrv->nPCs;
  nDevsLnR_c.allocate(1,npcRec);
  zGrA_xy.allocate(1,nSXs,mnYr,mxYr);
  zGrB_xy.allocate(1,nSXs,mnYr,mxYr);
  idxDevsLnC_fy.allocate(1,nFsh,mnYr,mxYr);
yDms = ptrMC->dimYrsToR;//years (mny:mxy)
xDms = ptrMC->dimSXsToR;//sex
mDms = ptrMC->dimMSsToR;//maturity
sDms = ptrMC->dimSCsToR;//shell condition
fDms = ptrMC->dimFshToR;//fisheries
vDms = ptrMC->dimSrvToR;//surveys
ypDms = ptrMC->dimYrsP1ToR;//years (mny:asy)
zbDms = ptrMC->dimZBsToR;//size bin midpoints
zpDms = ptrMC->dimZPsToR;//size bin midpoints (alternative)
zcDms = ptrMC->dimZCsToR;//size bin cuptoints
  hasF_fy.allocate(1,nFsh,mnYr,mxYr);
ctrProcCalls        = 0;
ctrProcCallsInPhase = 0;
    if (mseOpModMode){
        PRINT2B1("#--Creating ptrOMI")
        ptrOMI = new MSE_OpModInfo(ptrMC);
        ad_comm::change_datafile_name("OpModStateFile.txt");
        PRINT2B1("#--Reading OpModStateFile.txt")
        (*ad_comm::global_datafile)>>(*ptrOMI);
        PRINT2B1("#--Finished reading OpModStateFile.txt")
        prjR = 0.0;
        PRINT2B1("#--Reading TAC.txt")
        ad_comm::change_datafile_name("TAC.txt");
        (*ad_comm::global_datafile)>>inpTAC;
        (*ad_comm::global_datafile)>>inpOFL;
        PRINT2B1("#--Finished reading TAC.txt")
        rpt::echo<<"TAC, OFL is: "<<inpTAC<<tb<<inpOFL<<endl;
        cout<<"TAC, OFL is: "<<inpTAC<<tb<<inpOFL<<endl;
    }
PRINT2B1("#finished DATA_SECTION")
}

void model_parameters::initializationfunction(void)
{
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
PRINT2B1("#Starting PARAMETER_SECTION")
  pLnR.allocate(1,npLnR,lbLnR,ubLnR,phsLnR,"pLnR");
cout<<"pLnR = "<<pLnR<<endl;
  pRCV.allocate(1,npRCV,lbRCV,ubRCV,phsRCV,"pRCV");
cout<<"pRCV = "<<pRCV<<endl;
  pRX.allocate(1,npRX,lbRX,ubRX,phsRX,"pRX");
cout<<"pRX = "<<pRX<<endl;
  pRa.allocate(1,npRa,lbRa,ubRa,phsRa,"pRa");
cout<<"pRa = "<<pRa<<endl;
  pRb.allocate(1,npRb,lbRb,ubRb,phsRb,"pRb");
cout<<"pRb = "<<pRb<<endl;
  pDevsLnR.allocate(1,npDevsLnR,mniDevsLnR,mxiDevsLnR,lbDevsLnR,ubDevsLnR,phsDevsLnR,"pDevsLnR");
for (int p=1;p<=npDevsLnR;p++) cout<<"pDevsLnR["<<p<<"] = "<<pDevsLnR[p]<<endl;
  devsLnR.allocate(1,npDevsLnR,mniDevsLnR,mxiDevsLnR+1,"devsLnR");
  #ifndef NO_AD_INITIALIZE
    devsLnR.initialize();
  #endif
cout<<"got past recruitment parameters"<<endl;
  pM.allocate(1,npM,lbM,ubM,phsM,"pM");
cout<<"pM = "<<pM<<endl;
  pDM1.allocate(1,npDM1,lbDM1,ubDM1,phsDM1,"pDM1");
cout<<"pDM1 = "<<pDM1<<endl;
  pDM2.allocate(1,npDM2,lbDM2,ubDM2,phsDM2,"pDM2");
cout<<"pDM2 = "<<pDM2<<endl;
  pDM3.allocate(1,npDM3,lbDM3,ubDM3,phsDM3,"pDM3");
cout<<"pDM3 = "<<pDM3<<endl;
  pDM4.allocate(1,npDM4,lbDM4,ubDM4,phsDM4,"pDM4");
cout<<"pDM4 = "<<pDM4<<endl;
cout<<"got past natural mortality parameters"<<endl;
  pGrA.allocate(1,npGrA,lbGrA,ubGrA,phsGrA,"pGrA");
cout<<"pGrA = "<<pGrA<<endl;
  pGrB.allocate(1,npGrB,lbGrB,ubGrB,phsGrB,"pGrB");
cout<<"pGrB = "<<pGrB<<endl;
  pGrBeta.allocate(1,npGrBeta,lbGrBeta,ubGrBeta,phsGrBeta,"pGrBeta");
cout<<"pGrBeta = "<<pGrBeta<<endl;
cout<<"got past growth parameters"<<endl;
  pvLgtPrM2M.allocate(1,npLgtPrMat,mniLgtPrMat,mxiLgtPrMat,lbLgtPrMat,ubLgtPrMat,phsLgtPrMat,"pvLgtPrM2M");
for (int p=1;p<=npLgtPrMat;p++) cout<<"pvLgtPrM2M["<<p<<"] = "<<pvLgtPrM2M[p]<<endl;
cout<<"got past maturity parameters"<<endl;
  pS1.allocate(1,npS1,lbS1,ubS1,phsS1,"pS1");
cout<<"pS1 = "<<pS1<<endl;
  pS2.allocate(1,npS2,lbS2,ubS2,phsS2,"pS2");
cout<<"pS2 = "<<pS2<<endl;
  pS3.allocate(1,npS3,lbS3,ubS3,phsS3,"pS3");
cout<<"pS3 = "<<pS3<<endl;
  pS4.allocate(1,npS4,lbS4,ubS4,phsS4,"pS4");
cout<<"pS4 = "<<pS4<<endl;
  pS5.allocate(1,npS5,lbS5,ubS5,phsS5,"pS5");
cout<<"pS5 = "<<pS5<<endl;
  pS6.allocate(1,npS6,lbS6,ubS6,phsS6,"pS6");
cout<<"pS6 = "<<pS6<<endl;
  pDevsS1.allocate(1,npDevsS1,mniDevsS1,mxiDevsS1,lbDevsS1,ubDevsS1,phsDevsS1,"pDevsS1");
for (int p=1;p<=npDevsS1;p++) cout<<"pDevsS1["<<p<<"] = "<<pDevsS1[p]<<endl;
  pDevsS2.allocate(1,npDevsS2,mniDevsS2,mxiDevsS2,lbDevsS2,ubDevsS2,phsDevsS2,"pDevsS2");
for (int p=1;p<=npDevsS2;p++) cout<<"pDevsS2["<<p<<"] = "<<pDevsS2[p]<<endl;
  pDevsS3.allocate(1,npDevsS3,mniDevsS3,mxiDevsS3,lbDevsS3,ubDevsS3,phsDevsS3,"pDevsS3");
for (int p=1;p<=npDevsS3;p++) cout<<"pDevsS3["<<p<<"] = "<<pDevsS3[p]<<endl;
  pDevsS4.allocate(1,npDevsS4,mniDevsS4,mxiDevsS4,lbDevsS4,ubDevsS4,phsDevsS4,"pDevsS4");
for (int p=1;p<=npDevsS4;p++) cout<<"pDevsS4["<<p<<"] = "<<pDevsS4[p]<<endl;
  pDevsS5.allocate(1,npDevsS5,mniDevsS5,mxiDevsS5,lbDevsS5,ubDevsS5,phsDevsS5,"pDevsS5");
for (int p=1;p<=npDevsS5;p++) cout<<"pDevsS5["<<p<<"] = "<<pDevsS5[p]<<endl;
  pDevsS6.allocate(1,npDevsS6,mniDevsS6,mxiDevsS6,lbDevsS6,ubDevsS6,phsDevsS6,"pDevsS6");
for (int p=1;p<=npDevsS6;p++) cout<<"pDevsS6["<<p<<"] = "<<pDevsS6[p]<<endl;
  devsS1.allocate(1,npDevsS1,mniDevsS1,mxiDevsS1+1,"devsS1");
  #ifndef NO_AD_INITIALIZE
    devsS1.initialize();
  #endif
  devsS2.allocate(1,npDevsS2,mniDevsS2,mxiDevsS2+1,"devsS2");
  #ifndef NO_AD_INITIALIZE
    devsS2.initialize();
  #endif
  devsS3.allocate(1,npDevsS3,mniDevsS3,mxiDevsS3+1,"devsS3");
  #ifndef NO_AD_INITIALIZE
    devsS3.initialize();
  #endif
  devsS4.allocate(1,npDevsS4,mniDevsS4,mxiDevsS4+1,"devsS4");
  #ifndef NO_AD_INITIALIZE
    devsS4.initialize();
  #endif
  devsS5.allocate(1,npDevsS5,mniDevsS5,mxiDevsS5+1,"devsS5");
  #ifndef NO_AD_INITIALIZE
    devsS5.initialize();
  #endif
  devsS6.allocate(1,npDevsS6,mniDevsS6,mxiDevsS6+1,"devsS6");
  #ifndef NO_AD_INITIALIZE
    devsS6.initialize();
  #endif
  pvNPSel.allocate(1,npNPSel,mniNPSel,mxiNPSel,lbNPSel,ubNPSel,phsNPSel,"pvNPSel");
for (int p=1;p<=npNPSel;p++) cout<<"pvNPSel["<<p<<"] = "<<pvNPSel[p]<<endl;
cout<<"got past selectivity parameters"<<endl;
  pHM.allocate(1,npHM,lbHM,ubHM,phsHM,"pHM");
cout<<"pHM = "<<pHM<<endl;
  pLnC.allocate(1,npLnC,lbLnC,ubLnC,phsLnC,"pLnC");
cout<<"pLnC = "<<pLnC<<endl;
  pDC1.allocate(1,npDC1,lbDC1,ubDC1,phsDC1,"pDC1");
cout<<"pDC1 = "<<pDC1<<endl;
  pDC2.allocate(1,npDC2,lbDC2,ubDC2,phsDC2,"pDC2");
cout<<"pDC2 = "<<pDC2<<endl;
  pDC3.allocate(1,npDC3,lbDC3,ubDC3,phsDC3,"pDC3");
cout<<"pDC3 = "<<pDC3<<endl;
  pDC4.allocate(1,npDC4,lbDC4,ubDC4,phsDC4,"pDC4");
cout<<"pDC4 = "<<pDC4<<endl;
  pDevsLnC.allocate(1,npDevsLnC,mniDevsLnC,mxiDevsLnC,lbDevsLnC,ubDevsLnC,phsDevsLnC,"pDevsLnC");
for (int p=1;p<=npDevsLnC;p++) cout<<"pDevsLnC["<<p<<"] = "<<pDevsLnC[p]<<endl;
cout<<npLnEffX<<tb<<lbLnEffX<<tb<<ubLnEffX<<tb<<phsLnEffX<<endl;
  pLnEffX.allocate(1,npLnEffX,lbLnEffX,ubLnEffX,phsLnEffX,"pLnEffX");
cout<<"pLnEffX = "<<pLnEffX<<endl;
cout<<npLgtRet<<tb<<lbLgtRet<<tb<<ubLgtRet<<tb<<phsLgtRet<<endl;
  pLgtRet.allocate(1,npLgtRet,lbLgtRet,ubLgtRet,phsLgtRet,"pLgtRet");
cout<<"pLgtRet = "<<pLgtRet<<endl;
  devsLnC.allocate(1,npDevsLnC,mniDevsLnC,mxiDevsLnC+1,"devsLnC");
  #ifndef NO_AD_INITIALIZE
    devsLnC.initialize();
  #endif
cout<<"got past capture rate parameters"<<endl;
  pQ.allocate(1,npQ,lbQ,ubQ,phsQ,"pQ");
cout<<"pQ = "<<pQ<<endl;
  pDQ1.allocate(1,npDQ1,lbDQ1,ubDQ1,phsDQ1,"pDQ1");
cout<<"pDQ1 = "<<pDQ1<<endl;
  pDQ2.allocate(1,npDQ2,lbDQ2,ubDQ2,phsDQ2,"pDQ2");
cout<<"pDQ2 = "<<pDQ2<<endl;
  pDQ3.allocate(1,npDQ3,lbDQ3,ubDQ3,phsDQ3,"pDQ3");
cout<<"pDQ3 = "<<pDQ3<<endl;
  pDQ4.allocate(1,npDQ4,lbDQ4,ubDQ4,phsDQ4,"pDQ4");
cout<<"pDQ4 = "<<pDQ4<<endl;
cout<<"got past survey catchability parameters"<<endl;
  pMSE_LnC.allocate(1,npMSE_LnC,lbMSE_LnC,ubMSE_LnC,phsMSE_LnC,"pMSE_LnC");
cout<<"pMSE_LnC = "<<pMSE_LnC<<endl;
cout<<"got past MSE parameters"<<endl;
  objFun.allocate("objFun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  spB_yx.allocate(mnYr,mxYr,1,nSXs,"spB_yx");
  #ifndef NO_AD_INITIALIZE
    spB_yx.initialize();
  #endif
  n_yxmsz.allocate(mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"n_yxmsz");
  #ifndef NO_AD_INITIALIZE
    n_yxmsz.initialize();
  #endif
  nmN_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"nmN_yxmsz");
  #ifndef NO_AD_INITIALIZE
    nmN_yxmsz.initialize();
  #endif
  tmN_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"tmN_yxmsz");
  #ifndef NO_AD_INITIALIZE
    tmN_yxmsz.initialize();
  #endif
  initMnR.allocate("initMnR");
  #ifndef NO_AD_INITIALIZE
  initMnR.initialize();
  #endif
  R_y.allocate(mnYr,mxYr,"R_y");
  #ifndef NO_AD_INITIALIZE
    R_y.initialize();
  #endif
  Rx_c.allocate(1,npcRec,"Rx_c");
  #ifndef NO_AD_INITIALIZE
    Rx_c.initialize();
  #endif
  R_yx.allocate(mnYr,mxYr,1,nSXs,"R_yx");
  #ifndef NO_AD_INITIALIZE
    R_yx.initialize();
  #endif
  R_cz.allocate(1,npcRec,1,nZBs,"R_cz");
  #ifndef NO_AD_INITIALIZE
    R_cz.initialize();
  #endif
  R_yz.allocate(mnYr,mxYr,1,nZBs,"R_yz");
  #ifndef NO_AD_INITIALIZE
    R_yz.initialize();
  #endif
  R_yxz.allocate(mnYr,mxYr,1,nSXs,1,nZBs,"R_yxz");
  #ifndef NO_AD_INITIALIZE
    R_yxz.initialize();
  #endif
  stdvDevsLnR_c.allocate(1,npcRec,"stdvDevsLnR_c");
  #ifndef NO_AD_INITIALIZE
    stdvDevsLnR_c.initialize();
  #endif
  devsLnR_cy.allocate(1,npcRec,mnYr,mxYr,"devsLnR_cy");
  #ifndef NO_AD_INITIALIZE
    devsLnR_cy.initialize();
  #endif
  zscrDevsLnR_cy.allocate(1,npcRec,mnYr,mxYr,"zscrDevsLnR_cy");
  #ifndef NO_AD_INITIALIZE
    zscrDevsLnR_cy.initialize();
  #endif
  M_c.allocate(1,npcNM,"M_c");
  #ifndef NO_AD_INITIALIZE
    M_c.initialize();
  #endif
  M_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"M_yxmsz");
  #ifndef NO_AD_INITIALIZE
    M_yxmsz.initialize();
  #endif
  prM2M_cz.allocate(1,npcM2M,1,nZBs,"prM2M_cz");
  #ifndef NO_AD_INITIALIZE
    prM2M_cz.initialize();
  #endif
  prM2M_yxz.allocate(mnYr,mxYr,1,nSXs,1,nZBs,"prM2M_yxz");
  #ifndef NO_AD_INITIALIZE
    prM2M_yxz.initialize();
  #endif
  mnGrZ_cz.allocate(1,npcGrw,1,nZBs,"mnGrZ_cz");
  #ifndef NO_AD_INITIALIZE
    mnGrZ_cz.initialize();
  #endif
  prGr_czz.allocate(1,npcGrw,1,nZBs,1,nZBs,"prGr_czz");
  #ifndef NO_AD_INITIALIZE
    prGr_czz.initialize();
  #endif
  mnGrZ_yxsz.allocate(mnYr,mxYr,1,nSXs,1,nSCs,1,nZBs,"mnGrZ_yxsz");
  #ifndef NO_AD_INITIALIZE
    mnGrZ_yxsz.initialize();
  #endif
  prGr_yxszz.allocate(mnYr,mxYr,1,nSXs,1,nSCs,1,nZBs,1,nZBs,"prGr_yxszz");
  #ifndef NO_AD_INITIALIZE
    prGr_yxszz.initialize();
  #endif
  grA_xy.allocate(1,nSXs,mnYr,mxYr+1,"grA_xy");
  #ifndef NO_AD_INITIALIZE
    grA_xy.initialize();
  #endif
  grB_xy.allocate(1,nSXs,mnYr,mxYr+1,"grB_xy");
  #ifndef NO_AD_INITIALIZE
    grB_xy.initialize();
  #endif
  grBeta_xy.allocate(1,nSXs,mnYr,mxYr+1,"grBeta_xy");
  #ifndef NO_AD_INITIALIZE
    grBeta_xy.initialize();
  #endif
  npSel_cz.allocate(1,npNPSel,1,nZBs,"npSel_cz");
  #ifndef NO_AD_INITIALIZE
    npSel_cz.initialize();
  #endif
  sel_cz.allocate(1,npcSel,1,nZBs,"sel_cz");
  #ifndef NO_AD_INITIALIZE
    sel_cz.initialize();
  #endif
  sel_cyz.allocate(1,npcSel,mnYr,mxYr+1,1,nZBs,"sel_cyz");
  #ifndef NO_AD_INITIALIZE
    sel_cyz.initialize();
  #endif
  avgFc_nxms.allocate(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,"avgFc_nxms");
  #ifndef NO_AD_INITIALIZE
    avgFc_nxms.initialize();
  #endif
  avgFc2Eff_nxms.allocate(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,"avgFc2Eff_nxms");
  #ifndef NO_AD_INITIALIZE
    avgFc2Eff_nxms.initialize();
  #endif
  obsFc_nxmsy.allocate(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr,"obsFc_nxmsy");
  #ifndef NO_AD_INITIALIZE
    obsFc_nxmsy.initialize();
  #endif
  prdFc_nxmsy.allocate(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr,"prdFc_nxmsy");
  #ifndef NO_AD_INITIALIZE
    prdFc_nxmsy.initialize();
  #endif
  prdEff_nxmsy.allocate(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr,"prdEff_nxmsy");
  #ifndef NO_AD_INITIALIZE
    prdEff_nxmsy.initialize();
  #endif
  zscrEffX_nxmsy.allocate(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr,"zscrEffX_nxmsy");
  #ifndef NO_AD_INITIALIZE
    zscrEffX_nxmsy.initialize();
  #endif
  nllEffX_nxms.allocate(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,"nllEffX_nxms");
  #ifndef NO_AD_INITIALIZE
    nllEffX_nxms.initialize();
  #endif
  hmF_fy.allocate(1,nFsh,mnYr,mxYr,"hmF_fy");
  #ifndef NO_AD_INITIALIZE
    hmF_fy.initialize();
  #endif
  dvsLnC_fy.allocate(1,nFsh,mnYr,mxYr,"dvsLnC_fy");
  #ifndef NO_AD_INITIALIZE
    dvsLnC_fy.initialize();
  #endif
  cpF_fyxms.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,"cpF_fyxms");
  #ifndef NO_AD_INITIALIZE
    cpF_fyxms.initialize();
  #endif
  sel_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"sel_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    sel_fyxmsz.initialize();
  #endif
  ret_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"ret_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    ret_fyxmsz.initialize();
  #endif
  cpF_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"cpF_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    cpF_fyxmsz.initialize();
  #endif
  rmF_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"rmF_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    rmF_fyxmsz.initialize();
  #endif
  dmF_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"dmF_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    dmF_fyxmsz.initialize();
  #endif
  tmF_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"tmF_yxmsz");
  #ifndef NO_AD_INITIALIZE
    tmF_yxmsz.initialize();
  #endif
  cpN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"cpN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    cpN_fyxmsz.initialize();
  #endif
  dsN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"dsN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    dsN_fyxmsz.initialize();
  #endif
  rmN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"rmN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    rmN_fyxmsz.initialize();
  #endif
  dmN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"dmN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    dmN_fyxmsz.initialize();
  #endif
cout<<"got here 12"<<endl;
  cpB_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"cpB_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    cpB_fyxmsz.initialize();
  #endif
  dsB_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"dsB_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    dsB_fyxmsz.initialize();
  #endif
  rmB_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"rmB_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    rmB_fyxmsz.initialize();
  #endif
  dmB_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"dmB_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    dmB_fyxmsz.initialize();
  #endif
cout<<"got here 13"<<endl;
  mb_vyx.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,"mb_vyx");
  #ifndef NO_AD_INITIALIZE
    mb_vyx.initialize();
  #endif
  q_vyxms.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,"q_vyxms");
  #ifndef NO_AD_INITIALIZE
    q_vyxms.initialize();
  #endif
  a_vyxmsz.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"a_vyxmsz");
  #ifndef NO_AD_INITIALIZE
    a_vyxmsz.initialize();
  #endif
  s_vyxmsz.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"s_vyxmsz");
  #ifndef NO_AD_INITIALIZE
    s_vyxmsz.initialize();
  #endif
  q_vyxmsz.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"q_vyxmsz");
  #ifndef NO_AD_INITIALIZE
    q_vyxmsz.initialize();
  #endif
  n_vyxmsz.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"n_vyxmsz");
  #ifndef NO_AD_INITIALIZE
    n_vyxmsz.initialize();
  #endif
  b_vyxmsz.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"b_vyxmsz");
  #ifndef NO_AD_INITIALIZE
    b_vyxmsz.initialize();
  #endif
  prj_hmF_f.allocate(1,nFsh,"prj_hmF_f");
  #ifndef NO_AD_INITIALIZE
    prj_hmF_f.initialize();
  #endif
  prj_capF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs,"prj_capF_fmsz");
  #ifndef NO_AD_INITIALIZE
    prj_capF_fmsz.initialize();
  #endif
  prj_retF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs,"prj_retF_fmsz");
  #ifndef NO_AD_INITIALIZE
    prj_retF_fmsz.initialize();
  #endif
  prj_selF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs,"prj_selF_fmsz");
  #ifndef NO_AD_INITIALIZE
    prj_selF_fmsz.initialize();
  #endif
  prj_spB_x.allocate(1,nSXs,"prj_spB_x");
  #ifndef NO_AD_INITIALIZE
    prj_spB_x.initialize();
  #endif
  prj_n_xmsz.allocate(1,nSXs,1,nMSs,1,nSCs,1,nZBs,"prj_n_xmsz");
  #ifndef NO_AD_INITIALIZE
    prj_n_xmsz.initialize();
  #endif
  prj_cpN_fxmsz.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"prj_cpN_fxmsz");
  #ifndef NO_AD_INITIALIZE
    prj_cpN_fxmsz.initialize();
  #endif
  prj_rmN_fxmsz.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"prj_rmN_fxmsz");
  #ifndef NO_AD_INITIALIZE
    prj_rmN_fxmsz.initialize();
  #endif
  prj_dmN_fxmsz.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"prj_dmN_fxmsz");
  #ifndef NO_AD_INITIALIZE
    prj_dmN_fxmsz.initialize();
  #endif
  prj_dsN_fxmsz.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"prj_dsN_fxmsz");
  #ifndef NO_AD_INITIALIZE
    prj_dsN_fxmsz.initialize();
  #endif
  prjRetCatchMortBio_fx.allocate(1,nFsh,1,nSXs,"prjRetCatchMortBio_fx");
  #ifndef NO_AD_INITIALIZE
    prjRetCatchMortBio_fx.initialize();
  #endif
  prjDscCatchMortBio_fx.allocate(1,nFsh,1,nSXs,"prjDscCatchMortBio_fx");
  #ifndef NO_AD_INITIALIZE
    prjDscCatchMortBio_fx.initialize();
  #endif
  prjTotCatchMortBio_fx.allocate(1,nFsh,1,nSXs,"prjTotCatchMortBio_fx");
  #ifndef NO_AD_INITIALIZE
    prjTotCatchMortBio_fx.initialize();
  #endif
  prj_n_vxmsz.allocate(1,nSrv,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"prj_n_vxmsz");
  #ifndef NO_AD_INITIALIZE
    prj_n_vxmsz.initialize();
  #endif
  fPenRecDevs.allocate(1,npDevsLnR,"fPenRecDevs");
  #ifndef NO_AD_INITIALIZE
    fPenRecDevs.initialize();
  #endif
  fPenSmoothLgtPrMat.allocate(1,npLgtPrMat,"fPenSmoothLgtPrMat");
  #ifndef NO_AD_INITIALIZE
    fPenSmoothLgtPrMat.initialize();
  #endif
  fPenNonDecLgtPrMat.allocate(1,npLgtPrMat,"fPenNonDecLgtPrMat");
  #ifndef NO_AD_INITIALIZE
    fPenNonDecLgtPrMat.initialize();
  #endif
  fPenSmoothNPSel.allocate(1,npNPSel,"fPenSmoothNPSel");
  #ifndef NO_AD_INITIALIZE
    fPenSmoothNPSel.initialize();
  #endif
  fPenDevsS1.allocate(1,npDevsS1,"fPenDevsS1");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS1.initialize();
  #endif
  fPenDevsS2.allocate(1,npDevsS2,"fPenDevsS2");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS2.initialize();
  #endif
  fPenDevsS3.allocate(1,npDevsS3,"fPenDevsS3");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS3.initialize();
  #endif
  fPenDevsS4.allocate(1,npDevsS4,"fPenDevsS4");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS4.initialize();
  #endif
  fPenDevsS5.allocate(1,npDevsS5,"fPenDevsS5");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS5.initialize();
  #endif
  fPenDevsS6.allocate(1,npDevsS6,"fPenDevsS6");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS6.initialize();
  #endif
  fPenDevsLnC.allocate(1,npDevsLnC,"fPenDevsLnC");
  #ifndef NO_AD_INITIALIZE
    fPenDevsLnC.initialize();
  #endif
  R_z.allocate(1,nZBs,"R_z");
  #ifndef NO_AD_INITIALIZE
    R_z.initialize();
  #endif
  S1_msz.allocate(1,nMSs,1,nSCs,1,nZBs,"S1_msz");
  #ifndef NO_AD_INITIALIZE
    S1_msz.initialize();
  #endif
  Th_sz.allocate(1,nSCs,1,nZBs,"Th_sz");
  #ifndef NO_AD_INITIALIZE
    Th_sz.initialize();
  #endif
  T_szz.allocate(1,nSCs,1,nZBs,1,nZBs,"T_szz");
  #ifndef NO_AD_INITIALIZE
    T_szz.initialize();
  #endif
  S2_msz.allocate(1,nMSs,1,nSCs,1,nZBs,"S2_msz");
  #ifndef NO_AD_INITIALIZE
    S2_msz.initialize();
  #endif
  n_xmsz.allocate(1,nSXs,1,nMSs,1,nSCs,1,nZBs,"n_xmsz");
  #ifndef NO_AD_INITIALIZE
    n_xmsz.initialize();
  #endif
  sdrLnR_y.allocate(mnYr,mxYr,"sdrLnR_y");
  sdrSpB_xy.allocate(1,nSXs,mnYr,mxYr,"sdrSpB_xy");
  lkMMB.allocate("lkMMB");
  lkQM.allocate("lkQM");
  lkQF.allocate("lkQF");
  lkNMIM.allocate("lkNMIM");
  lkNMIF.allocate("lkNMIF");
  lkNMMM.allocate("lkNMMM");
  lkNMMF.allocate("lkNMMF");
PRINT2B1("#finished PARAMETER_SECTION")
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
    PRINT2B1("#Starting PRELIMINARY_CALCS_SECTION")
    int debug=1;
    
    PRINT2B1("PRELIMINARY_CALCS: ALL MODEL RUNS")
    //create population projection objects
    PRINT2B1("--Creating population projection objects")
    pPDI = new PopDyInfo(nZBs);         //  generic population dynamics info
    pCDI = new CatchInfo(nZBs,nFsh);    //  generic catch info
    pPPr = new PopProjector(pPDI,pCDI); //  generic population projector
    PRINT2B1("--Created population projection objects")
    //set initial values for all parameters
    if (usePin) {
        PRINT2B1("NOTE: setting initial values for parameters using pin file")
    } else {
        PRINT2B1("NOTE: setting initial values for parameters using MPI")
    }
    setInitVals(0,rpt::echo);
    
    //calculate average effort for fisheries over specified time periods and 
    //allocate associated arrays
    PRINT2B1(" ")
    PRINT2B1("--calculating average effort")
    rpt::echo<<"mapD2MFsh = "<<mapD2MFsh<<endl;    
    rpt::echo<<"nEASs     = "<<nEASs<<endl;    
    if (nEASs){
        EffAvgScenarios* ptrEASs = ptrMOs->ptrEffXtrapScenarios->ptrEffAvgScenarios;
        for (int n=1;n<=nEASs;n++){//effort averaging scenarios
            EffAvgScenario* ptrEAS = ptrEASs->ppEASs[n-1];
            int fm = ptrEAS->f;    //index for fishery associated with this averaging scenario
            int fd = mapM2DFsh(fm);//index for corresponding fishery data object
            rpt::echo<<"n = "<<n<<". fm = "<<fm<<". fd = "<<fd<<endl;
            if (!ptrMDS->ppFsh[fd-1]->ptrEff){
                cout<<"---------------------------------------------------"<<endl;
                cout<<"No effort data given for "<<ptrMC->lblsFsh[fm]<<","<<endl;
                cout<<"but effort averaging requested! Aborting..."<<endl;
                exit(-1);
            }
            rpt::echo<<"fishery has effort"<<endl;
            //extract years for effort averaging, keeping model limits in mind
            yrsAvgEff_ny(n).deallocate();
            yrsAvgEff_ny(n) = wts::extractVector(mnYr,mxYr,ptrEAS->ptrIB->getFwdIndexVector());
            //assign values for effort averaging
            eff_ny(n).deallocate();
            eff_ny(n) = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y(yrsAvgEff_ny(n));
            rpt::echo<<"eff_ny(n)       = "<<eff_ny(n)<<endl;
            rpt::echo<<"yrsAvgEff_ny(n) = "<<yrsAvgEff_ny(n)<<endl;
            avgEff_n(n) = sum(eff_ny(n))/yrsAvgEff_ny(n).size();
            rpt::echo<<"avgEff_n(n) = "<<avgEff_n(n)<<endl;
        }//n
        rpt::echo<<"eff_ny       = "<<endl<<eff_ny<<endl;
        rpt::echo<<"yrsAvgEff_ny = "<<endl<<yrsAvgEff_ny<<endl;
        rpt::echo<<"avgEff_n     = "<<avgEff_n<<endl;
        rpt::echo<<"finished calculating average effort scenarios"<<endl;
        cout<<"finished calculating average effort scenarios"<<endl;
        PRINT2B1(" ")
        PRINT2B1("--starting allocation of arrays for average capture rates and effort extrapolation")
        obsEff_nxmsy.initialize();
        for (int n=1;n<=nCRASs;n++){
            CapRateAvgScenario* ptrCRAS = ptrMOs->ptrEffXtrapScenarios->ptrCapRateAvgScenarios->ppCRASs[n-1];
            int idEAS = ptrCRAS->idEffAvgInfo;
            int mnx, mxx, mnm, mxm, mns, mxs;
            mnx = mxx = ptrCRAS->x;//sex index
            if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
            mnm = mxm = ptrCRAS->m;//maturity index
            if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
            mns = mxs = ptrCRAS->s;//shell index
            if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
            for (int x=mnx;x<=mxx;x++){
                for (int m=mnm;m<=mxm;m++){
                    for (int s=mns;s<=mxs;s++) {
                        for (int iy=yrsAvgEff_ny(idEAS).indexmin();iy<=yrsAvgEff_ny(idEAS).indexmax();iy++)
                            obsEff_nxmsy(n,x,m,s)(yrsAvgEff_ny(idEAS,iy)) = eff_ny(idEAS,iy);
                        rpt::echo<<"n = "<<n<<". idEAS = "<<idEAS<<endl;
                        rpt::echo<<"yrsAvgEff_ny(idEAS) = "<<yrsAvgEff_ny(idEAS)<<endl;
                        rpt::echo<<"eff_ny(idEAS)       ="<<eff_ny(idEAS)<<endl;
                        rpt::echo<<"obsEff_nxmsy("<<n<<cc<<x<<cc<<m<<cc<<s<<") = "<<obsEff_nxmsy(n,x,m,s)<<endl;
                    }//s
                }//m
            }//x
        }//n
        PRINT2B1("--finished allocation of arrays for average capture rates and effort extrapolation")
    } else {
        PRINT2B1("--NO effort averaging scenarios defined!")
    }
    if (!mseMode){
        PRINT2B1("--PRELIMINARY_CALCS: NON-MSE MODEL RUNS ONLY")
        {
            PRINT2B1("writing effective MPI after setInitVals")
            ofstream os; os.open("effectiveMPI.dat", ios::trunc);
            os.precision(12);
            ptrMPI->setToWriteVectorInitialValues(true);
            os<<(*ptrMPI)<<endl;
            os.close();
            PRINT2B1("finished writing effective MPI after setInitVals")
        }
        PRINT2B1("testing setAllDevs()")
        setAllDevs(tcsam::dbgAll,rpt::echo);
        PRINT2B1("finished testing setAllDevs()")
        {
            PRINT2B1("writing data to R")
            ofstream os; os.open("ModelData.R", ios::trunc);
            os.precision(12);
            ReportToR_Data(os,0,cout);
            PRINT2B1("finished writing data to R")
        }
        {
            PRINT2B1("writing parameters info to R")
            ofstream os; os.open("ModelParametersInfo.R", ios::trunc);
            os.precision(12);
            ptrMPI->writeToR(os);
            os.close();
            PRINT2B1("finished writing parameters info to R")
            //write initial parameter values to csv
            PRINT2B1("writing parameters info to csv")
            ofstream os1("tcsam02.params.all.init.csv", ios::trunc);
            os1.precision(12);
            writeParameters(os1,0,0);//all parameters
            os1.close();
            ofstream os2("tcsam02.params.active.init.csv", ios::trunc);
            os2.precision(12);
            writeParameters(os2,0,1);//only parameters that will be active (i.e., phase>0)
            os2.close();
            PRINT2B1("finished writing parameters info to csv")
        }
        if (!mcevalOn) {
            //this section runs for an "ordinary" model run
            int dbgLevel = 0; //set at dbgCalcProcs+1 to print to terminal 
            PRINT2B1("testing calcRecruitment():")
            calcRecruitment(dbgLevel,rpt::echo);
            PRINT2B1("testing calcNatMort():")
            calcNatMort(dbgLevel,rpt::echo);
            PRINT2B1("testing calcGrowth():")
            calcGrowth(dbgLevel,rpt::echo);
            PRINT2B1("testing calcPrM2M():")
            calcPrM2M(dbgCalcProcs+1,rpt::echo);
            PRINT2B1("testing calcSelectivities():")
            calcSelectivities(dbgLevel,rpt::echo);
            PRINT2B1("testing calcFisheryFs():")
            calcFisheryFs(dbgCalcProcs+1,rpt::echo);
            PRINT2B1("testing calcSurveyQs():")
            calcSurveyQs(dbgLevel,cout);
            if (!runAlt){
                PRINT2B1("testing runPopDyMod():")
                runPopDyMod(dbgLevel,cout);
                rpt::echo<<"n_yxm:"<<endl;
                for (int y=mnYr;y<=(mxYr+1);y++){
                    for (int x=1;x<=nSXs;x++){
                        for (int m=1;m<=nMSs;m++){
                           rpt::echo<<y<<cc;
                           rpt::echo<<tcsam::getSexType(x)<<cc;
                           rpt::echo<<tcsam::getMaturityType(m)<<cc;
                           rpt::echo<<sum(n_yxmsz(y,x,m))<<endl;
                        }
                    }
                }
                rpt::echo<<"n_yxmsz:"<<endl;
                for (int y=mnYr;y<=(mxYr+1);y++){
                    for (int x=1;x<=nSXs;x++){
                        for (int m=1;m<=nMSs;m++){
                            for (int s=1;s<=nSCs;s++){
                               rpt::echo<<y<<cc;
                               rpt::echo<<tcsam::getSexType(x)<<cc;
                               rpt::echo<<tcsam::getMaturityType(m)<<cc;
                               rpt::echo<<tcsam::getShellType(s)<<cc;
                               rpt::echo<<n_yxmsz(y,x,m,s)<<endl;
                            }
                        }
                    }
                }
            } else {
                PRINT2B1("testing runAltPopDyMod():")
                runAltPopDyMod(dbgLevel,cout);
                rpt::echo<<"n_yxm:"<<endl;
                for (int y=mnYr;y<=(mxYr+1);y++){
                    for (int x=1;x<=nSXs;x++){
                        for (int m=1;m<=nMSs;m++){
                           rpt::echo<<y<<cc;
                           rpt::echo<<tcsam::getSexType(x)<<cc;
                           rpt::echo<<tcsam::getMaturityType(m)<<cc;
                           rpt::echo<<sum(n_yxmsz(y,x,m))<<endl;
                        }
                    }
                }
                rpt::echo<<"n_yxmsz:"<<endl;
                for (int y=mnYr;y<=(mxYr+1);y++){
                    for (int x=1;x<=nSXs;x++){
                        for (int m=1;m<=nMSs;m++){
                            for (int s=1;s<=nSCs;s++){
                               rpt::echo<<y<<cc;
                               rpt::echo<<tcsam::getSexType(x)<<cc;
                               rpt::echo<<tcsam::getMaturityType(m)<<cc;
                               rpt::echo<<tcsam::getShellType(s)<<cc;
                               rpt::echo<<n_yxmsz(y,x,m,s)<<endl;
                            }
                        }
                    }
                }
            }
            if (doOFL&&debugOFL){
                PRINT2B1("Testing OFL calculations")
                ofstream echoOFL; echoOFL.open("calcOFL.init.txt", ios::trunc);
                echoOFL.precision(12);
                echoOFL<<"----Testing calcOFL()"<<endl;
                calcOFL(mxYr+1,debugOFL,echoOFL);//updates ptrOFLResults
                ptrOFLResults->writeCSVHeader(echoOFL); echoOFL<<endl;
                ptrOFLResults->writeToCSV(echoOFL);     echoOFL<<endl;
                echoOFL<<"----Finished testing calcOFL()!"<<endl;
                echoOFL.close();
                PRINT2B1("Finished testing OFL calculations!")
            }
            if (fitSimData){
                int dbgLevel = 0;
                PRINT2B1("creating sim data to fit in model")
                createSimData(dbgLevel,rpt::echo,iSimDataSeed,ptrMDS);//stochastic if iSimDataSeed<>0
                {
                    PRINT2B1("re-writing data to R")
                    ofstream echo1; echo1.open("ModelData.R", ios::trunc);
                    echo1.precision(12);
                    ReportToR_Data(echo1,0,cout);
                }
            }
            {
                PRINT2B1("--Testing calcObjFun()")
                if (!runAlt) runPopDyMod(0,cout); else runAltPopDyMod(0,cout);
                calcObjFun(dbgAll,rpt::echo);
                PRINT2B2("--Finished testing calcObjFun(): ",objFun)
            }
            {
                //write objective function components only
                PRINT2B1("Writing model fits to R")
                ofstream os0("tcsam02.ModelFits.init.R", ios::trunc);
                os0.precision(12);
                ReportToR_ModelFits(os0,-1.0,0,cout);
                os0.close();
                PRINT2B1("Finished writing model fits to R")
            }
            {
                PRINT2B1("writing model sim data to file")
                int dbgLevel = 0;
                createSimData(dbgLevel,rpt::echo,0,ptrSimMDS);//deterministic
                ofstream echo1; echo1.open("tcsam02.SimData.init.dat", ios::trunc);
                echo1.precision(12);
                writeSimData(echo1,0,rpt::echo,ptrSimMDS);
                echo1.close();
                PRINT2B1("finished writing model sim data to file")
            }
            {
                //must do this last because call to calcDynB0 in ReportToR
                //sets fishing mortality to 0 across all fleets
                //and re-runs population model
                PRINT2B1("writing initial model report to R")
                ofstream echo1; echo1.open("tcsam02.init.rep", ios::trunc);
                echo1.precision(12);
                ReportToR(echo1,-1.0,1,cout);
                echo1.close();
                PRINT2B1("finished writing model report to R")
            }        
        } else {
            writeMCMCHeader();
            PRINT2B1("MCEVAL is on")
        }//if mcevalOn
        PRINT2B1("")
        PRINT2B2("obj fun = ",objFun)
    } else {
         if (mseOpModMode){
            PRINT2B1("PRELIMINARY_CALCS: MSE OpModMode")
            dvector vLnR_y = log(ptrOMI->R_y(1982,mxYr));
            cout<<"vLnR_y = "<<vLnR_y<<endl;
            double mn = mean(vLnR_y);
            double sd = sqrt(wts::variance(vLnR_y));
            prjR = mfexp(wts::drawSampleNormal(rng, mn, sd));
            PRINT2B2("mean recruitment: ",mn);
            PRINT2B2("stdv recruitment: ",sd);
            PRINT2B2("expected total recruitment: ",mfexp(mn+square(sd)/2.0));
            PRINT2B2("projected total recruitment: ",prjR);
            if (inpTAC>0.0){
                dvariable mseCapF = mfexp(pMSE_LnC[1]);
                projectPopForTAC(mseCapF,dbgAll,rpt::echo);
                calcObjFunForTAC(dbgAll,rpt::echo);
            } else {
                PRINT2B1("#--TAC is 0: directed fishery is closed!!");
                PRINT2B1("#--RUNNING OpMod without fitting for TAC.")
                finishOpModMode();
                PRINT2B1("#--Exiting model without fitting for TAC.")
                exit(0);
            }
        } else if (mseEstModMode){
             PRINT2B1("PRELIMINARY_CALCS: MSE EstModMode")
        }
    }
    PRINT2B1("#finished PRELIMINARY_CALCS_SECTION")
    PRINT2B1("#----------------------------------")
    
}

void model_parameters::userfunction(void)
{
  objFun =0.0;
    int dbg = 0; //dbgAll;
    ctrProcCalls++;       //increment procedure section calls counter
    ctrProcCallsInPhase++;//increment in-phase procedure section calls counter
    if (dbg>=dbgObjFun){
        PRINT2B1(" ")
        PRINT2B1("--PROCEDURE_SECTION----------------")
    }
    if (mc_phase()&&(option_match(ad_comm::argc,ad_comm::argv,"-nuts")>0)){
        adstring msg;
        PRINT2B1("--Running NUTS MCMC-------")
        msg = "---phase="+str(current_phase())+". ctrProcCalls = "+str(ctrProcCalls)+tb+str(ctrProcCallsInPhase);
        PRINT2B1(msg)
        msg = "number of active params = "+str(initial_params::nvarcalc())+tb+str(initial_params::num_active_calc());
        PRINT2B1(msg);
        rpt::echo<<"number of active params = "<<initial_params::nvarcalc()<<cc<<initial_params::num_active_calc()<<endl;
        //writeParameters(std::cout,0,1);
        writeParameters(rpt::echo,0,1);
    }
    if (ctrDebugParams&&(ctrProcCalls>=ctrDebugParams)){
        adstring msg;
        PRINT2B1("--writing parameters to file for debugging")
        msg = "---phase="+str(current_phase())+". ctrProcCalls = "+str(ctrProcCalls)+tb+str(ctrProcCallsInPhase);
        PRINT2B1(msg)
        writeParameters(std::cout,0,1);
        writeParameters(rpt::echo,0,1);
        adstring fn = "tcsam02.Debug."+str(current_phase())+"."+str(ctrProcCalls)+".rep";
        ofstream os; os.open((char*) fn, ios::trunc);
        os.precision(12);
        PRINT2B1("--writing report to R file for debugging")
        ReportToR(os,-1.0,1,cout);
        PRINT2B1("--finished writing report to R file for debugging")
    }
    if (checkParams(0,std::cout)){
        adstring msg;
        PRINT2B1("--checking parameters for debugging")
        msg = "---phase="+str(current_phase())+". ctrProcCalls = "+str(ctrProcCalls)+tb+str(ctrProcCallsInPhase);
        PRINT2B1(msg)
        checkParams(1,std::cout);
        checkParams(1,rpt::echo);
        writeParameters(std::cout,0,1);
        writeParameters(rpt::echo,0,1);
        PRINT2B1("--exiting...")
        ad_exit(-1);
    }
    if (mseOpModMode){
        dvariable mseCapF = mfexp(pMSE_LnC[1]);//multiplier on capture rate in directed fishery
        projectPopForTAC(mseCapF,0,cout);
        calcObjFunForTAC(1000,cout);
    } else {
        if (!runAlt) runPopDyMod(0,cout); else runAltPopDyMod(0,cout);
        calcObjFun(dbg,rpt::echo);
    }
    if ((!mseOpModMode)&&(ctrProcCallsInPhase==1)){
        //write objective function components only
        adstring fn = "tcsam02.ModelFits."+itoa(current_phase(),10)+"-"+str(ctrProcCallsInPhase)+"."+"R";
        ofstream os0(fn, ios::trunc);
        os0.precision(12);
        ReportToR_ModelFits(os0,-1.0,0,cout);
        os0.close();
        PRINT2B2("just after ReportToR_ModelFits. obj fun = ",objFun)
        PRINT2B1(" ")            
    }
    //assign values to likelihood profile numbers
    lkMMB = spB_yx(mxYr,MALE); 
    lkQM  = q_vyxms(1,mxYr+1,  MALE,MATURE,NEW_SHELL); 
    lkQF  = q_vyxms(1,mxYr+1,FEMALE,MATURE,NEW_SHELL); 
    lkNMIM = M_yxmsz(mxYr,  MALE,IMMATURE,NEW_SHELL)(1); 
    lkNMIF = M_yxmsz(mxYr,FEMALE,IMMATURE,NEW_SHELL)(1);
    lkNMMM = M_yxmsz(mxYr,  MALE,MATURE,NEW_SHELL)(1); 
    lkNMMF = M_yxmsz(mxYr,FEMALE,MATURE,NEW_SHELL)(1);
    if (sd_phase()){
        sdrLnR_y = log(R_y);
        for (int x=1;x<=nSXs;x++){
            for (int y=mnYr; y<=mxYr; y++){
                sdrSpB_xy(x,y) = spB_yx(y,x);
            }
        }
    }
    if (mceval_phase()){
        updateMPI(0, cout);
        writeMCMCtoR(mcmc);
    }
    if (dbg>=dbgObjFun) {PRINT2B1("--END PROCEDURE_SECTION----------------")}
#ifdef DEBUG
  std::cout << "DEBUG: Total gradient stack used is " << gradient_structure::get()->GRAD_STACK1->total() << " out of " << gradient_structure::get_GRADSTACK_BUFFER_SIZE() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->GRAD_LIST->total_addresses() << " out of " << gradient_structure::get_MAX_DLINKS() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->ARR_LIST1->get_max_last_offset() << " out of " << gradient_structure::get_ARRAY_MEMBLOCK_SIZE() << std::endl;;
#endif
}

void model_parameters::runAltPopDyMod(int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting runAltPopDyMod()"<<endl;
    //initialize population model
    initAltPopDyMod(debug, cout);
    //run population model
    for (int y=mnYr;y<=mxYr;y++){
        doSurveys(y,debug,cout);
       if (debug>=dbgPopDy) cout<<"--year = "<<y<<endl;        
         for (int x=1;x<=nSXs;x++){
            if (debug>=dbgPopDy) cout<<"----sex = "<<x<<endl;
            runAltPopDyModOneYear(y,x,debug,cout);  
            if (debug>=dbgPopDy) {cout<<"------R_z    = "<<endl; cout<<R_yxz(y,x)<<endl;}
        }
        if (debug>=dbgPopDy) cout<<endl;
    }
    doSurveys(mxYr+1,debug,cout);//do final surveys
    if (debug>=dbgPopDy) cout<<"finished runAltPopDyMod()"<<endl;
}

void model_parameters::initAltPopDyMod(int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting initAltPopDyMod()"<<endl;
    spB_yx.initialize();   //annual spawning biomass
    n_yxmsz.initialize();  //annual population abundance
    nmN_yxmsz.initialize();//number of crab killed by natural mortality
    tmN_yxmsz.initialize();//total number of crab killed 
    tmF_yxmsz.initialize();//total fishing mortality rate 
    setAllDevs(debug,cout);//set devs vectors
    calcRecruitment(debug,cout);//calculate recruitment
    calcNatMort(debug,cout);    //calculate natural mortality rates
    calcGrowth(debug,cout);     //calculate growth transition matrices
    calcPrM2M(debug,cout);      //calculate maturity ogives
    calcSelectivities(debug,cout); //calculate selectivity functions
    calcFisheryFs(debug,cout);     //calculate fishery F's
    calcSurveyQs(debug,cout);      //calculate survey Q's
    if (ptrMOs->optInitNatZ==0){
        //will build up population from recruitment (like TCSAM2013)
        //do nothing, because n_yxmsz has already been initialized to 0
    } else if (ptrMOs->optInitNatZ==1){
        //use equilibrium calculation to set initial n-at-z (like gmacs)
        //assumes no fishing occurs before model start
        calcEqNatZF100(initMnR,mnYr,debug,cout);//calculate n_xmsz
        n_yxmsz(mnYr) = n_xmsz;
    } else {
        cout<<"Unrecognized option for initial n-at-z: "<<ptrMOs->optInitNatZ<<endl;
        cout<<"Terminating!"<<endl;
        exit(-1);
    }
    if (debug>=dbgPopDy) cout<<"finished initAltPopDyMod()"<<endl;
}

void model_parameters::runAltPopDyModOneYear(int y, int x, int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"Starting runAltPopDyModOneYear("<<y<<cc<<x<<")"<<endl;
   //1. Update population rates, based on year and sex
    pPDI->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(x);//assign weight-at-length
    pPDI->R_z   = R_yz(y);                   //relative recruitment-at-size
    pPDI->M_msz = M_yxmsz(y,x);              //rates of natural mortality
    pPDI->T_szz = prGr_yxszz(y,x);           //growth transition matrices
    for (int s=1;s<=nSCs;s++) pPDI->Th_sz(s) = prM2M_yxz(y,x);//pr(terminal molt|size)
    //if (debug) cout<<"updated pPDI for sex "<<x<<" in "<<y<<endl;
    //2. Update fishery conditions based on year and sex
    //TODO: move declarations to PARAMETER_SECTION
    dvar_vector hmF_f(1,nFsh);//handling mortality
    dvar4_array capF_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);//fishery capture rates
    dvar4_array retF_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);//fishery retention functions
    dvar4_array selF_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);//fishery selectivity functions
    hmF_f.initialize();
    capF_fmsz.initialize();
    retF_fmsz.initialize();
    selF_fmsz.initialize();
    for (int f=1;f<=nFsh;f++){
        hmF_f(f) = hmF_fy(f,y);
            for (int m=1;m<=nMSs;m++){
                            capF_fmsz(f,m) = cpF_fyxmsz(f,y,x,m);
                            retF_fmsz(f,m) = ret_fyxmsz(f,y,x,m);
                            selF_fmsz(f,m) = sel_fyxmsz(f,y,x,m);
            }//m
    }//f
    //update CatchInfo objects
    pCDI->setCaptureRates(capF_fmsz);
    pCDI->setRetentionFcns(retF_fmsz);
    pCDI->setSelectivityFcns(selF_fmsz);
    pCDI->setHandlingMortality(hmF_f);
    if (debug) cout<<"updated pCDI for sex "<<x<<" in "<<y<<endl;
    //3. run PopProjector
    pPPr->dtF = dtF_y(y); //time at which fishery occurs
    pPPr->dtM = dtM_y(y); //time at which mating occurs
    dvariable dirF = -1.0;//don't change scale on directed fishery F's
    n_yxmsz(y+1,x)       = pPPr->project(dirF,n_yxmsz(y,x),cout);
    n_yxmsz(y+1,x,IMMATURE,NEW_SHELL) += R_yxz(y,x);
    spB_yx(y,x)          = pPPr->getMatureBiomassAtMating();
    dvar4_array cpN_fmsz = pPPr->getFisheriesCaptureAbundance();
    dvar4_array rmN_fmsz = pPPr->getRetainedCatchMortality();
    dvar4_array dmN_fmsz = pPPr->getDiscardCatchMortality();
    for (int f=1;f<=nFsh;f++){
        for (int m=1;m<=nMSs;m++){
            cpN_fyxmsz(f,y,x,m) = cpN_fmsz(f,m);                          //capture abundance
            rmN_fyxmsz(f,y,x,m) = rmN_fmsz(f,m);                          //retained mortality abundance
            dmN_fyxmsz(f,y,x,m) = dmN_fmsz(f,m);                          //discard mortality abundance
            dsN_fyxmsz(f,y,x,m) = cpN_fyxmsz(f,y,x,m)-rmN_fyxmsz(f,y,x,m);//total discard abundance
        }//m
    }//f
    if (debug) cout<<"ran pPPr for sex "<<x<<" in "<<y<<endl;
    if (debug>=dbgPopDy) cout<<"finished runAltPopDyModOneYear("<<y<<cc<<x<<")"<<endl;
}

void model_parameters::projectPopForTAC(dvariable& mseCapF, int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"projectPopForTAC()"<<endl;
    prjRetCatchMortBio_fx.initialize();
    prjDscCatchMortBio_fx.initialize();
    prjTotCatchMortBio_fx.initialize();
    prj_cpN_fxmsz.initialize(); //capture abundance
    prj_rmN_fxmsz.initialize(); //retained mortality abundance
    prj_dmN_fxmsz.initialize(); //discard mortality abundance
    prj_dsN_fxmsz.initialize(); //total discard abundance
    dvariable maxCapF = 0.0;
    for (int x=1;x<=nSXs;x++){   
        //1. Update population rates, based on year and sex
        pPDI->w_mz  = ptrOMI->wAtZ_xmz(x);//assign weight-at-length
        pPDI->R_z   = ptrOMI->R_z;                    //relative recruitment-at-size
        pPDI->M_msz = ptrOMI->M_xmsz(x);               //rates of natural mortality
        pPDI->T_szz = ptrOMI->prGr_xszz(x);            //growth transition matrices
        //pr(terminal molt|size)
        for (int s=1;s<=nSCs;s++) pPDI->Th_sz(s) = ptrOMI->prM2M_xz(x);
        //if (debug) cout<<"updated pPDI for sex "<<x<<" in "<<y<<endl;
        //2. Update fishery conditions based on year and sex
        prj_hmF_f = ptrOMI->hmF_f;
        prj_capF_fmsz.initialize();
        prj_retF_fmsz.initialize();
        prj_selF_fmsz.initialize();
        for (int f=1;f<=nFsh;f++){
            for (int m=1;m<=nMSs;m++){
                prj_capF_fmsz(f,m) = ptrOMI->cpF_fxmsz(f,x,m);
                prj_retF_fmsz(f,m) = ptrOMI->ret_fxmsz(f,x,m);
                prj_selF_fmsz(f,m) = ptrOMI->sel_fxmsz(f,x,m);
            }//m
        }//f
        //determine proposed rates for directed fishery that will
        //result in TAC being taken
        prj_capF_fmsz(1) = mseCapF*prj_selF_fmsz(1);
        //update CatchInfo objects
        pCDI->setHandlingMortality(prj_hmF_f);
        pCDI->setCaptureRates(prj_capF_fmsz);
        pCDI->setRetentionFcns(prj_retF_fmsz);
        pCDI->setSelectivityFcns(prj_selF_fmsz);
        if (debug) cout<<"updated pCDI for sex "<<x<<endl;
        //3. run PopProjector
        pPPr->dtF = ptrOMI->dtF; //time at which fishery occurs
        pPPr->dtM = ptrOMI->dtM; //time at which mating occurs
        dvariable dirF = -1.0;//don't change scale on directed fishery F's
        dvar3_array n_msz   = ptrOMI->n_xmsz(x);//sex "x" population abundance at start of year
        prj_n_xmsz(x)       = pPPr->project(dirF,n_msz,cout);
        prj_n_xmsz(x,IMMATURE,NEW_SHELL) += prjR*ptrOMI->R_x(x)*ptrOMI->R_z;
        prj_spB_x(x)         = pPPr->getMatureBiomassAtMating();
        dvar4_array cpN_fmsz = pPPr->getFisheriesCaptureAbundance();
        dvar4_array rmN_fmsz = pPPr->getRetainedCatchMortality();
        dvar4_array dmN_fmsz = pPPr->getDiscardCatchMortality();
        for (int f=1;f<=nFsh;f++){
            for (int m=1;m<=nMSs;m++){
                prj_cpN_fxmsz(f,x,m) = cpN_fmsz(f,m); //capture abundance
                prj_rmN_fxmsz(f,x,m) = rmN_fmsz(f,m); //retained mortality abundance
                prj_dmN_fxmsz(f,x,m) = dmN_fmsz(f,m); //discard mortality abundance
                //total discard abundance
                prj_dsN_fxmsz(f,x,m) = prj_cpN_fxmsz(f,x,m)-prj_rmN_fxmsz(f,x,m);
                for (int s=1;s<=nSCs;s++){
                    //retained catch mortality (biomass)
                    prjRetCatchMortBio_fx(f,x) += prj_rmN_fxmsz(f,x,m,s)*ptrOMI->wAtZ_xmz(x,m);
                    //discarded catch mortality (biomass)
                    prjDscCatchMortBio_fx(f,x) += prj_dmN_fxmsz(f,x,m,s)*ptrOMI->wAtZ_xmz(x,m);
                }//s
            }//m
            prjTotCatchMortBio_fx(f,x) = prjRetCatchMortBio_fx(f,x) + prjDscCatchMortBio_fx(f,x);
        }//f
        if (debug) cout<<"ran pPPr for sex "<<x<<endl;
    }//x
    if (debug>=dbgPopDy) cout<<"finished projectPopForTAC()"<<endl;
}

void model_parameters::projectPopForZeroTAC(int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"projectPopForZeroTAC()"<<endl;
    prjRetCatchMortBio_fx.initialize();
    prjDscCatchMortBio_fx.initialize();
    prjTotCatchMortBio_fx.initialize();
    prj_cpN_fxmsz.initialize(); //capture abundance
    prj_rmN_fxmsz.initialize(); //retained mortality abundance
    prj_dmN_fxmsz.initialize(); //discard mortality abundance
    prj_dsN_fxmsz.initialize(); //total discard abundance
    dvariable maxCapF = 0.0;
    for (int x=1;x<=nSXs;x++){   
        //1. Update population rates, based on year and sex
        pPDI->w_mz  = ptrOMI->wAtZ_xmz(x);   //assign weight-at-length
        pPDI->R_z   = ptrOMI->R_z;           //relative recruitment-at-size
        pPDI->M_msz = ptrOMI->M_xmsz(x);     //rates of natural mortality
        pPDI->T_szz = ptrOMI->prGr_xszz(x);  //growth transition matrices
        //pr(terminal molt|size)
        for (int s=1;s<=nSCs;s++) pPDI->Th_sz(s) = ptrOMI->prM2M_xz(x);
        //if (debug) cout<<"updated pPDI for sex "<<x<<" in "<<y<<endl;
        //2. Update fishery conditions based on year and sex
        prj_hmF_f = ptrOMI->hmF_f;
        prj_capF_fmsz.initialize();
        prj_retF_fmsz.initialize();
        prj_selF_fmsz.initialize();
        for (int f=1;f<=nFsh;f++){
            for (int m=1;m<=nMSs;m++){
                prj_capF_fmsz(f,m) = ptrOMI->cpF_fxmsz(f,x,m);
                prj_retF_fmsz(f,m) = ptrOMI->ret_fxmsz(f,x,m);
                prj_selF_fmsz(f,m) = ptrOMI->sel_fxmsz(f,x,m);
            }//m
        }//f
        //directed fishery is closed,so set capture rates to 0
        prj_capF_fmsz(1) = 0.0*prj_selF_fmsz(1);
        //update CatchInfo objects
        pCDI->setHandlingMortality(prj_hmF_f);
        pCDI->setCaptureRates(prj_capF_fmsz);
        pCDI->setRetentionFcns(prj_retF_fmsz);
        pCDI->setSelectivityFcns(prj_selF_fmsz);
        if (debug) cout<<"updated pCDI for sex "<<x<<endl;
        //3. run PopProjector
        pPPr->dtF = ptrOMI->dtF; //time at which fishery occurs
        pPPr->dtM = ptrOMI->dtM; //time at which mating occurs
        dvariable dirF = -1.0;//don't change scale on directed fishery F's
        dvar3_array n_msz   = ptrOMI->n_xmsz(x);//sex "x" population abundance at start of year
        prj_n_xmsz(x)       = pPPr->project(dirF,n_msz,cout);
        prj_n_xmsz(x,IMMATURE,NEW_SHELL) += prjR*ptrOMI->R_x(x)*ptrOMI->R_z;
        prj_spB_x(x)         = pPPr->getMatureBiomassAtMating();
        dvar4_array cpN_fmsz = pPPr->getFisheriesCaptureAbundance();
        dvar4_array rmN_fmsz = pPPr->getRetainedCatchMortality();
        dvar4_array dmN_fmsz = pPPr->getDiscardCatchMortality();
        for (int f=1;f<=nFsh;f++){
            for (int m=1;m<=nMSs;m++){
                prj_cpN_fxmsz(f,x,m) = cpN_fmsz(f,m); //capture abundance
                prj_rmN_fxmsz(f,x,m) = rmN_fmsz(f,m); //retained mortality abundance
                prj_dmN_fxmsz(f,x,m) = dmN_fmsz(f,m); //discard mortality abundance
                //total discard abundance
                prj_dsN_fxmsz(f,x,m) = prj_cpN_fxmsz(f,x,m)-prj_rmN_fxmsz(f,x,m);
                for (int s=1;s<=nSCs;s++){
                    //retained catch mortality (biomass)
                    prjRetCatchMortBio_fx(f,x) += prj_rmN_fxmsz(f,x,m,s)*ptrOMI->wAtZ_xmz(x,m);
                    //discarded catch mortality (biomass)
                    prjDscCatchMortBio_fx(f,x) += prj_dmN_fxmsz(f,x,m,s)*ptrOMI->wAtZ_xmz(x,m);
                }//s
            }//m
            prjTotCatchMortBio_fx(f,x) = prjRetCatchMortBio_fx(f,x) + prjDscCatchMortBio_fx(f,x);
        }//f
        if (debug) cout<<"ran pPPr for sex "<<x<<endl;
    }//x
    if (debug>=dbgPopDy) cout<<"finished projectPopForZeroTAC()"<<endl;
}

void model_parameters::runPopDyMod(int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting runPopDyMod()"<<endl;
    //initialize population model
    initPopDyMod(0, cout);
    //run population model
    for (int y=mnYr;y<=mxYr;y++){
        doSurveys(y,0,cout);
        runPopDyModOneYear(y,debug,cout);        
    }
    doSurveys(mxYr+1,0,cout);//do final surveys
    if (debug>=dbgPopDy) cout<<"finished runPopDyMod()"<<endl;
}

void model_parameters::initPopDyMod(int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting initPopDyMod()"<<endl;
    spB_yx.initialize();
    n_yxmsz.initialize();
    nmN_yxmsz.initialize();
    tmN_yxmsz.initialize();
    tmF_yxmsz.initialize();
    setAllDevs(debug,cout);//set devs vectors
    calcRecruitment(debug,cout);//calculate recruitment
    calcNatMort(debug,cout);    //calculate natural mortality rates
    calcGrowth(debug,cout);     //calculate growth transition matrices
    calcPrM2M(debug,cout);      //calculate maturity ogives
    calcSelectivities(debug,cout); //calculate selectivity functions
    calcFisheryFs(debug,cout);     //calculate fishery F's
    calcSurveyQs(debug,cout);      //calculate survey Q's
    if (ptrMOs->optInitNatZ==0){
        //will build up population from recruitment (like TCSAM2013)
        //do nothing, because n_yxmsz has already been initialized to 0
    } else if (ptrMOs->optInitNatZ==1){
        //use equilibrium calculation to set initial n-at-z (like gmacs)
        //assumes no fishing occurs before model start
        calcEqNatZF100(initMnR,mnYr,debug,cout);//calculate n_xmsz
        n_yxmsz(mnYr) = n_xmsz;
    } else {
        cout<<"Unrecognized option for initial n-at-z: "<<ptrMOs->optInitNatZ<<endl;
        cout<<"Terminating!"<<endl;
        exit(-1);
    }
    if (debug>=dbgPopDy) cout<<"finished initPopDyMod()"<<endl;
}

void model_parameters::runPopDyModOneYear(int y, int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"Starting runPopDyModOneYear("<<y<<")"<<endl;
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n2_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n3_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n4_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n5_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    n1_xmsz.initialize();
    n2_xmsz.initialize();
    n3_xmsz.initialize();
    n4_xmsz.initialize();
    n5_xmsz.initialize();
    if (dtF_y(y)<=dtM_y(y)){//fishery occurs BEFORE molting/growth/maturity
        if (debug>=dbgPopDy) cout<<"Fishery occurs BEFORE molting/growth/maturity"<<endl;
        //apply natural mortality before fisheries
        n1_xmsz = applyNatMort(n_yxmsz(y),y,dtF_y(y),debug,cout);
        //conduct fisheries
        n2_xmsz = applyFshMort(n1_xmsz,y,debug,cout);
        //apply natural mortality from fisheries to molting/growth/maturity
        if (dtF_y(y)==dtM_y(y)) {
            n3_xmsz = n2_xmsz;
        } else {
            n3_xmsz = applyNatMort(n2_xmsz,y,dtM_y(y)-dtF_y(y),debug,cout);
        }
        //calc mature (spawning) biomass at time of mating, but BEFORE growth/maturity (TODO: does this make sense??)
        spB_yx(y) = calcSpB(n3_xmsz,y,debug,cout);
        //apply molting, growth and maturation
        n4_xmsz = applyMGM(n3_xmsz,y,debug,cout);
        //apply natural mortality to end of year
        if (dtM_y(y)==1.0) {
            n5_xmsz = n4_xmsz;
        } else {
            n5_xmsz = applyNatMort(n4_xmsz,y,1.0-dtM_y(y),debug,cout);
        }
    } else {              //fishery occurs AFTER molting/growth/maturity
        if (debug>=dbgPopDy) cout<<"Fishery occurs AFTER molting/growth/maturity"<<endl;
        //apply natural mortality before molting/growth/maturity
        n1_xmsz = applyNatMort(n_yxmsz(y),y,dtM_y(y),debug,cout);
        //calc mature (spawning) biomass at time of mating (TODO: does this make sense??)
        spB_yx(y) = calcSpB(n1_xmsz,y,debug,cout);
        //apply molting, growth and maturation
        n2_xmsz = applyMGM(n1_xmsz,y,debug,cout);
        //apply natural mortality from molting/growth/maturity to fisheries
        if (dtM_y(y)==dtF_y(y)) {
            n3_xmsz = n2_xmsz;
        } else {
            n3_xmsz = applyNatMort(n2_xmsz,y,dtF_y(y)-dtM_y(y),debug,cout);
        }
        //conduct fisheries
        n4_xmsz = applyFshMort(n3_xmsz,y,debug,cout);
        //apply natural mortality to end of year
        if (dtF_y(y)==1.0) {
            n5_xmsz = n4_xmsz;
        } else {
            n5_xmsz = applyNatMort(n4_xmsz,y,1.0-dtF_y(y),debug,cout);
        }
    }
    //advance surviving individuals to next year
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                n_yxmsz(y+1,x,m,s) = n5_xmsz(x,m,s);
            }
        }
    }
    //add in recruits (NOTE: R_y(y) here corresponds to R_y(y+1) in TCSAM2013)
    for (int x=1;x<=nSXs;x++) n_yxmsz(y+1,x,IMMATURE,NEW_SHELL) += R_yxz(y,x);
    if (debug>=dbgPopDy){
        cout<<"----year = "<<y<<endl;
        for (int x=1;x<=nSXs;x++){
            cout<<"----sex = "<<x<<endl;
            cout<<"------n0_msz = "<<endl; wts::print(n_yxmsz(y,x),cout,1);
            cout<<"------n1_msz = "<<endl; wts::print(n1_xmsz(x),cout,1);
            for (int f=1;f<=nFsh;f++){
                cout<<"------fishery = "<<f<<endl;
                cout<<"--------hmF_f = "<<hmF_fy(f,y)<<endl;
                cout<<"--------cpF_msz = "<<endl; wts::print(cpF_fyxmsz(f,y,x),cout,1);
                cout<<"--------rmF_msz = "<<endl; wts::print(rmF_fyxmsz(f,y,x),cout,1);
                cout<<"--------dmF_msz = "<<endl; wts::print(dmF_fyxmsz(f,y,x),cout,1);
            }//f
            cout<<"------n2_msz = "<<endl; wts::print(n2_xmsz(x),cout,1);
            cout<<"------n3_msz = "<<endl; wts::print(n3_xmsz(x),cout,1);
            cout<<"------n4_msz = "<<endl; wts::print(n4_xmsz(x),cout,1);
            cout<<"------n5_msz = "<<endl; wts::print(n5_xmsz(x),cout,1);
            cout<<"------R_z    = "<<endl; cout<<R_yxz(y,x)<<endl;
        }
        cout<<endl;
    }
    if (debug>=dbgPopDy) cout<<"finished runPopDyModOneYear("<<y<<")"<<endl;
}

void model_parameters::doSurveys(int y,int debug,ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting doSurveys("<<y<<")"<<endl;
    for (int v=1;v<=nSrv;v++){
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++){
                    n_vyxmsz(v,y,x,m,s) = elem_prod(q_vyxmsz(v,y,x,m,s),n_yxmsz(y,x,m,s));
                }
            }
        }
    }
    for (int v=1;v<=nSrv;v++){
        mb_vyx(v,y) = calcSpB(n_vyxmsz(v,y),y,debug,cout);
    }
    if (debug>=dbgPopDy) cout<<"finished doSurveys("<<y<<")"<<endl;
}

dvar_vector model_parameters::calcSpB(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting calcSpB("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar_vector spb(1,nSXs); spb.initialize();
    for (int x=1;x<=nSXs;x++){
        for (int s=1;s<=nSCs;s++) spb(x) += n0_xmsz(x,MATURE,s)*ptrMDS->ptrBio->wAtZ_xmz(x,MATURE);//dot product here
    }
    if (debug>dbgApply) cout<<"finished calcSpB("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return spb;
}

dvar4_array model_parameters::applyNatMort(dvar4_array& n0_xmsz, int y, double dt, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting applyNatMort("<<y<<cc<<dt<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                n1_xmsz(x,m,s) = elem_prod(mfexp(-M_yxmsz(y,x,m,s)*dt),n0_xmsz(x,m,s));//survivors
                nmN_yxmsz(y,x,m,s) += n0_xmsz(x,m,s)-n1_xmsz(x,m,s); //natural mortality
                tmN_yxmsz(y,x,m,s) += n0_xmsz(x,m,s)-n1_xmsz(x,m,s); //natural mortality
            }
        }
    }
    if (debug>dbgApply) cout<<"finished applyNatMort("<<y<<cc<<dt<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
}

dvar4_array model_parameters::applyFshMort(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting applyFshMort("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvector     tdF_z(1,nZBs);//for use in calculating fishing rate components
    dvar_vector tm_z(1,nZBs); //total mortality (numbers) by size
    dvar_vector tvF_z(1,nZBs);//for use in calculating fishing rate components
    dvar_vector tfF_z(1,nZBs);//for use in calculating fishing rate components
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);//numbers surviving fisheries
    n1_xmsz.initialize();
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                //tmF_yxmsz(y,x,m,s).initialize();//total fishing mortality rate
                for (int f=1;f<=nFsh;f++) tmF_yxmsz(y,x,m,s) += rmF_fyxmsz(f,y,x,m,s)+dmF_fyxmsz(f,y,x,m,s);
                n1_xmsz(x,m,s) = elem_prod(mfexp(-tmF_yxmsz(y,x,m,s)),n0_xmsz(x,m,s));//numbers surviving all fisheries
                tm_z = n0_xmsz(x,m,s)-n1_xmsz(x,m,s);                                 //numbers killed by all fisheries
                tmN_yxmsz(y,x,m,s) += tm_z;            //add in numbers killed by all fisheries to total killed
                //calculate fishing rate components (need to ensure NOT dividing by 0)
                tdF_z = value(tmF_yxmsz(y,x,m,s));
                tvF_z = elem_prod(1-wts::isEQ(tdF_z,0.0),tmF_yxmsz(y,x,m,s)) + wts::isEQ(tdF_z,0.0);//= tmF_yxmsz(y,x,m,s,z) if tmF_yxmsz(y,x,m,s,z) > 0, else = 1
                tfF_z = elem_div(1.0-mfexp(-tmF_yxmsz(y,x,m,s)),tvF_z);//= (1-exp(-tmF_yxmsz(y,x,m,s,z)))/tmF_yxmsz(y,x,m,s,z) if tmF_yxmsz(y,x,m,s,z) > 0, else = 1
                for (int f=1;f<=nFsh;f++){                   
                    cpN_fyxmsz(f,y,x,m,s) = elem_prod(elem_prod(cpF_fyxmsz(f,y,x,m,s),tfF_z),n0_xmsz(x,m,s)); //numbers captured in fishery f
                    rmN_fyxmsz(f,y,x,m,s) = elem_prod(elem_prod(rmF_fyxmsz(f,y,x,m,s),tfF_z),n0_xmsz(x,m,s)); //retained mortality in fishery f (numbers)
                    dmN_fyxmsz(f,y,x,m,s) = elem_prod(elem_prod(dmF_fyxmsz(f,y,x,m,s),tfF_z),n0_xmsz(x,m,s)); //discards mortality in fishery f (numbers)
                    dsN_fyxmsz(f,y,x,m,s) = cpN_fyxmsz(f,y,x,m,s)-rmN_fyxmsz(f,y,x,m,s);//discarded catch (NOT mortality) in fishery f (numbers)                    
                }
            }
        }
    }
    if (debug>dbgApply) cout<<"finished applyFshMort("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
}

dvar4_array model_parameters::applyMGM(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting applyMGM("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    n1_xmsz.initialize();
    for (int x=1;x<=nSXs;x++){
        dvar_vector np_z = prGr_yxszz(y,x,NEW_SHELL)*n0_xmsz(x,IMMATURE,NEW_SHELL);
        n1_xmsz(x,IMMATURE,NEW_SHELL) = elem_prod(1.0-prM2M_yxz(y,x),np_z);
        n1_xmsz(x,IMMATURE,OLD_SHELL) = 0.0;
        n1_xmsz(x,MATURE,NEW_SHELL)   = elem_prod(    prM2M_yxz(y,x),np_z);
        n1_xmsz(x,MATURE,OLD_SHELL)   = n0_xmsz(x,MATURE,NEW_SHELL)+n0_xmsz(x,MATURE,OLD_SHELL);
    }
    if (debug>dbgApply) cout<<"finished applyNatMGM("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
}

void model_parameters::calcRecruitment(int debug, ostream& cout)
{
    if (debug>dbgCalcProcs) cout<<"starting calcRecruitment()"<<endl;
    RecruitmentInfo* ptrRI = ptrMPI->ptrRec;
    R_y.initialize();
    Rx_c.initialize();
    R_yx.initialize();
    R_cz.initialize();
    R_yz.initialize();
    R_yxz.initialize();
    nDevsLnR_c.initialize();
    stdvDevsLnR_c.initialize();
    devsLnR_cy.initialize();
    zscrDevsLnR_cy.initialize();
    int k; int y;
    dvector dzs = zBs+(zBs[2]-zBs[1])/2.0-zBs[1];
    dvar_vector ptLnR  = ptrRI->pLnR->calcArithScaleVals(pLnR);
    dvar_vector ptRCV  = ptrRI->pRCV->calcArithScaleVals(pRCV);
    dvar_vector ptRX   = ptrRI->pRX->calcArithScaleVals(pRX);
    dvar_vector ptRa   = ptrRI->pRa->calcArithScaleVals(pRa);
    dvar_vector ptRb   = ptrRI->pRb->calcArithScaleVals(pRb);
    for (int pc=1;pc<=ptrRI->nPCs;pc++){
        ivector pids = ptrRI->getPCIDs(pc);
        k=ptrRI->nIVs+1;//first parameter variable column in ParameterComnbinations
        dvariable mnLnR = ptLnR(pids[k++]);
        dvariable cvR   = ptRCV(pids[k++]);
        dvariable xR    = ptRX(pids[k++]);
        dvariable aR    = ptRa(pids[k++]);
        dvariable bR    = ptRb(pids[k++]);
        if (debug>dbgCalcProcs){
            cout<<"pids  = "<<pids<<endl;
            cout<<"mnLnR = "<<mnLnR<<endl;
            cout<<"cvR   = "<<cvR<<endl;
            cout<<"xR    = "<<xR<<endl;
            cout<<"aR    = "<<aR<<endl;
            cout<<"bR    = "<<bR<<endl;
        }
        int useDevs = pids[k]; k++;
        dvariable mnR;   //mean recruitment
        dvariable varLnR;//ln-scale variance in recruitment
        dvar_vector dvsLnR;
        ivector idxDevsLnR;
        varLnR = log(1.0+square(cvR));    //ln-scale variance
        stdvDevsLnR_c(pc) = sqrt(varLnR); //ln-scale std dev
        mnR    = mfexp(mnLnR+varLnR/2.0); //mean recruitment
        if (useDevs) {
            dvsLnR     = devsLnR(useDevs);
            idxDevsLnR = idxsDevsLnR(useDevs);
            if (debug>dbgCalcProcs) {
                cout<<"lims(dvsLnR) = "<<dvsLnR.indexmin()<<cc<<dvsLnR.indexmax()<<endl;
                cout<<"idx(dvsLnR) = "<<idxDevsLnR<<endl;
                cout<<"dvsLnR = "<<dvsLnR<<endl;
            }
        }
        Rx_c(pc) = xR;
        R_cz(pc) = elem_prod(pow(dzs,(aR/bR)-1.0),mfexp(-dzs/bR));
        R_cz(pc) /= sum(R_cz(pc));//normalize to sum to 1
        imatrix idxs = ptrRI->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);
            if (y==mnYr) initMnR = mnR;
            if ((mnYr<=y)&&(y<=mxYr)){
                if (debug>dbgCalcProcs+10) cout<<"y,i = "<<y<<tb<<idxDevsLnR(y)<<endl;
                if (useDevs){
                    R_y(y) = mfexp(mnLnR+dvsLnR[idxDevsLnR[y]]);
                } else {
                    R_y(y) = mnR;
                }
                if (debug>dbgCalcProcs+10) cout<<"R_y(y)="<<R_y(y)<<tb;
                if (MALE==nSXs){
                    R_yx(y,MALE) = 1.0;//only tracking males
                } else {
                    R_yx(y,MALE)   = Rx_c(pc);
                    R_yx(y,FEMALE) = 1.0-R_yx(y,MALE);
                    if (debug>dbgCalcProcs+10) cout<<R_yx(y,MALE)<<endl;
                }
                R_yz(y) = R_cz(pc);
                if (debug>dbgCalcProcs+10) cout<<"R_yz(y)="<<R_yz(y)<<endl;
                for (int x=1;x<=nSXs;x++) R_yxz(y,x) = R_y(y)*R_yx(y,x)*R_yz(y);
                nDevsLnR_c(pc)++;//increment count of devs for this pc
                devsLnR_cy(pc,y) = dvsLnR(idxDevsLnR(y));//rec devs by pc and year
            } else {
                if (debug>dbgCalcProcs) cout<<"skipping y,i = "<<y<<tb<<idxDevsLnR(y)<<endl;
            }
        }//idx
        zscrDevsLnR_cy(pc) = devsLnR_cy(pc)/stdvDevsLnR_c(pc);//standardized zscores (assuming mean=0)
    }//pc
    if (debug>dbgCalcProcs) {
        cout<<"R_y = "<<R_y<<endl;
        cout<<"R_yx(MALE) = "<<column(R_yx,MALE)<<endl;
        cout<<"R_yz  = "<<endl<<R_yz<<endl;
        cout<<"nDevs = "<<nDevsLnR_c<<endl;
        cout<<"sigR  = "<<stdvDevsLnR_c<<endl;
        cout<<"devs  = "<<devsLnR_cy<<endl;
        cout<<"zscr  = "<<zscrDevsLnR_cy<<endl;
        cout<<"finished calcRecruitment()"<<endl;
    }
}

void model_parameters::calcNatMort(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcNatMort()"<<endl;
    NaturalMortalityInfo* ptrNM = ptrMPI->ptrNM;
    dvariable lnM;
    dvar_vector M_z(1,nZBs);
    M_c.initialize();
    M_yxmsz.initialize();
    int y; int mnx; int mxx; int mnm; int mxm; int mns; int mxs;
    dvar_vector ptM  = ptrNM->pM->calcArithScaleVals(pM);
    dvar_vector ptDM1  = ptrNM->pDM1->calcArithScaleVals(pDM1);
    dvar_vector ptDM2  = ptrNM->pDM2->calcArithScaleVals(pDM2);
    dvar_vector ptDM3  = ptrNM->pDM3->calcArithScaleVals(pDM3);
    dvar_vector ptDM4  = ptrNM->pDM4->calcArithScaleVals(pDM4);
    for (int pc=1;pc<=ptrNM->nPCs;pc++){
        if (debug>dbgCalcProcs) cout<<"pc = "<<pc<<endl;
        lnM.initialize();
        ivector pids = ptrNM->getPCIDs(pc);
        int k=ptrNM->nIVs+1;//1st parameter variable column
        if (debug>dbgCalcProcs) cout<<"pids = "<<pids(k,pids.indexmax())<<endl;
        if (ptrMOs->optParamNM==0){
            //add in base (arithmetic-scale) natural mortality
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptM["<<pids[k]<<"]: "<<ptM(pids[k])<<endl;
                lnM += log(ptM(pids[k]));
            } k++;
            //add in ln-scale offset 1
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptDM1["<<pids[k]<<"]: "<<ptDM1(pids[k])<<endl;
                lnM += ptDM1(pids[k]);
            } k++;
            //add in ln-scale offset 2
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptDM2["<<pids[k]<<"]: "<<ptDM2(pids[k])<<endl;
                lnM += ptDM2(pids[k]);
            } k++;
            //add in ln-scale offset 3
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptDM3["<<pids[k]<<"]: "<<ptDM3(pids[k])<<endl;
                lnM += ptDM3(pids[k]);
            } k++;
            //add in ln-scale offset 4
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptDM4["<<pids[k]<<"]: "<<ptDM4(pids[k])<<endl;
                lnM += ptDM4(pids[k]);
            }  k++; //advance k to zScaling in pids
        } else if (ptrMOs->optParamNM==1){
            //TCSAM2013 parameterization: arithmetic scale multipliers to base
            //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptM["<<pids[k]<<"]: "<<ptM(pids[k])<<endl;
                lnM += log(ptM(pids[k]));
            } k++;
            //multiply offset 1 (for immature crab)
            //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Multiplying by offset ptDM1["<<pids[k]<<"]: "<<ptDM1(pids[k])<<endl;
                lnM += log(ptDM1(pids[k]));//add on ln-scale
            } k++;
            //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Multiplying by offset ptDM2["<<pids[k]<<"]: "<<ptDM2(pids[k])<<endl;
                lnM += log(ptDM2(pids[k]));
            } k++;
            //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Multiplying by offset ptDM3["<<pids[k]<<"]: "<<ptDM3(pids[k])<<endl;
                lnM += log(ptDM3(pids[k]));
            } k++;
            //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Multiplying by offset ptDM4["<<pids[k]<<"]: "<<ptDM4(pids[k])<<endl;
                lnM += log(ptDM4(pids[k]));
            } k++;
        }//optParamNM
        //convert from ln-scale to arithmetic scale
        M_c(pc) = mfexp(lnM);
        if (debug>dbgCalcProcs) cout<<"lnM= "<<lnM<<tb<<"M_c = "<<M_c(pc)<<endl;
        //add in size-scaling, if requested
        M_z.initialize();
        //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
        if (pids[k]&&(current_phase()>=pids[k])) {
            if (debug>dbgCalcProcs) cout<<"adding size scaling"<<endl;
            M_z = M_c(pc)*(zMref/zBs);//factor in size dependence
        } else {
            if (debug>dbgCalcProcs) cout<<"not adding size scaling"<<endl;
            M_z = M_c(pc);//no size dependence
        }
        if (debug>dbgCalcProcs) cout<<"finished scaling"<<endl;
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrNM->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);//year
            if ((mnYr<=y)&&(y<=mxYr)){
                mnx = mxx = idxs(idx,2);//sex index
                if (mnx==tcsam::ALL_SXs){mnx = 1; mxx = tcsam::nSXs;}
                mnm = mxm = idxs(idx,3);//maturity state
                if (mnm==tcsam::ALL_MSs){mnm = 1; mxm = tcsam::nMSs;}
                mns = mxs = idxs(idx,4);//shell condition
                if (mns==tcsam::ALL_SCs){mns = 1; mxs = tcsam::nSCs;}
                for (int x=mnx;x<=mxx;x++){
                    for (int m=mnm;m<=mxm;m++){
                        for (int s=mns;s<=mxs;s++) M_yxmsz(y,x,m,s) = 1.0*M_z;
                    }
                }
            }
        }
    }//loop over pcs
    if (debug>dbgCalcProcs) {
        for (int y=mnYr;y<=mxYr;y++){
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) cout<<"M_yxmsz("<<y<<tb<<x<<tb<<m<<tb<<s<<")="<<M_yxmsz(y,x,m,s)<<endl;
                }
            }
        }
        cout<<"Finished calcNatMort()"<<endl;
    }
}

void model_parameters::calcPrM2M(int debug, ostream& cout)
{
    if (debug>dbgCalcProcs) cout<<"starting calcPrM2M()"<<endl;
    Molt2MaturityInfo* ptrM2MI = ptrMPI->ptrM2M;
    prM2M_cz.initialize();
    prM2M_yxz.initialize();
    int k; int y; int mnx; int mxx;
    for (int pc=1;pc<=ptrM2MI->nPCs;pc++){
        ivector pids = ptrM2MI->getPCIDs(pc);
        k=ptrM2MI->nIVs+1;//first parameter variable column in ParameterComnbinations
        BoundedVectorInfo* pBVI = (*ptrM2MI->pvLgtPrM2M)[pids[k]];
        dvar_vector lgtPrM2M = pBVI->calcArithScaleVals(pvLgtPrM2M(pids[k++]));
        ivector idxf = pBVI->getFwdIndices();//map indices to model size bins
        int imn = idxf.indexmin();
        int imx = idxf.indexmax();
        if (debug>dbgCalcProcs){
            cout<<"pc = "<<pc<<". mn = "<<imn<<", mx = "<<imx<<endl;
            cout<<"lgtPrM2M = "<<lgtPrM2M<<endl;
            cout<<"fwd indices = "<<idxf<<endl;
            ivector idxr = pBVI->getRevIndices();
            cout<<"rev indices min, max = "<<idxr.indexmin()<<cc<<idxr.indexmax()<<endl;
            cout<<"rev indices = "<<idxr<<endl;
        }
        prM2M_cz(pc) = 1.0;                                          //1 for size classes larger than max estimated
        if (idxf[1]>1) prM2M_cz(pc)(1,idxf[1]-1) = 0.0;              //0 for size classes smaller than min estimated
        dvar_vector prM2Mp = 1.0/(1.0+mfexp(-lgtPrM2M));             //logistic otherwise
        for (int i=imn;i<=imx;i++) prM2M_cz(pc)(idxf[i]) = prM2Mp[i];//assign to model size bins         
        if (debug>dbgCalcProcs){
            cout<<"prM2M = "<<prM2M_cz(pc)<<endl;
        }
        imatrix idxs = ptrM2MI->getModelIndices(pc);
        if (debug>dbgCalcProcs) cout<<"molt-to-maturity indices"<<endl<<idxs<<endl;
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);
            if ((mnYr<=y)&&(y<=mxYr)){
                mnx = mxx = idxs(idx,2);//sex index
                if (mnx==tcsam::ALL_SXs){mnx = 1; mxx = tcsam::nSXs;}
                if (debug>dbgCalcProcs) cout<<"y = "<<y<<tb<<"sex = "<<tcsam::getSexType(mnx)<<": "<<tcsam::getSexType(mnx)<<endl;
                for (int x=mnx;x<=mxx;x++) prM2M_yxz(y,x) = prM2M_cz(pc);//note: this change made a difference, but not sure why!
            }
        }
    }
    if (debug>dbgCalcProcs) {
        for (int y=mnYr;y<=mxYr;y++){
            for (int x=1;x<=nSXs;x++){
                cout<<"prM2M_yxz("<<y<<tb<<x<<")="<<prM2M_yxz(y,x)<<endl;
            }
        }
        cout<<"finished calcPrM2M()"<<endl;
    }
}

void model_parameters::calcGrowth(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcGrowth()"<<endl;
    GrowthInfo* ptrGrw = ptrMPI->ptrGrw;
    dvariable grA;
    dvariable grB;
    dvariable grBeta;
    double zGrA = 0.0;
    double zGrB = 0.0;
    mnGrZ_cz.initialize();
    prGr_czz.initialize();
    mnGrZ_yxsz.initialize();
    prGr_yxszz.initialize();
    grA_xy.initialize();
    grB_xy.initialize();
    grBeta_xy.initialize();
    zGrA_xy.initialize();
    zGrB_xy.initialize();
    dvar_matrix prGr_zz(1,nZBs,1,nZBs);
    int y; int mnx; int mxx;
    dvar_vector ptGrA    = ptrGrw->pGrA->calcArithScaleVals(pGrA);
    dvar_vector ptGrB    = ptrGrw->pGrB->calcArithScaleVals(pGrB);
    dvar_vector ptGrBeta = ptrGrw->pGrBeta->calcArithScaleVals(pGrBeta);
    if (debug>dbgCalcProcs){
        cout<<"pGrA:  "<<pGrA<<endl;
        cout<<"ptGrA: "<<ptGrA<<endl;
        cout<<"pGrB:  "<<pGrB<<endl;
        cout<<"ptGrB: "<<ptGrB<<endl;
        cout<<"pGrBeta:  "<<pGrBeta<<endl;
        cout<<"ptGrBeta: "<<ptGrBeta<<endl;
    }
    for (int pc=1;pc<=ptrGrw->nPCs;pc++){
        ivector pids = ptrGrw->getPCIDs(pc);
        int k=ptrGrw->nIVs+1;//1st parameter column
        grA = ptGrA(pids[k]); k++; //"a" coefficient for mean growth
        grB = ptGrB(pids[k]); k++; //"b" coefficient for mean growth
        grBeta = ptGrBeta(pids[k]); k++; //scale factor for gamma function growth transition
        if (debug>dbgCalcProcs){
            cout<<"pc: "<<pc<<tb<<"grA:"<<tb<<grA<<". grB:"<<tb<<grB<<". grBeta:"<<grBeta<<endl;
        }
        //compute growth transition matrix for this pc
        dvar_vector mnZs(1,nZBs); mnZs.initialize();  //mean post-molt sizes with zBs as pre-molt sizes
        if (ptrMOs->optGrowthParam==0){
            //TCSAM2013 parameterization with ln-scale intercept, slope
            mnZs = mfexp(grA+grB*log(zBs));
        } else if (ptrMOs->optGrowthParam==1){
            //parameterization at min, max pre-molt sizes
            dvector pXDs = ptrGrw->getPCXDs(pc);
            zGrA = pXDs[1]; //pre-molt size corresponding to pGrA as mean post-molt size
            zGrB = pXDs[2]; //pre-molt size corresponding to pGrB as mean post-molt size
            if (debug>dbgCalcProcs) cout<<"growth parameterization 1. "<<"zGrA:"<<tb<<zGrA<<". zGrB:"<<tb<<zGrB<<endl;
            mnZs = grA*mfexp(log(grB/grA)/log(zGrB/zGrA)*log(zBs/zGrA));
            if (debug>dbgCalcProcs) cout<<"mnZs:"<<tb<<mnZs<<endl;
        } else if (ptrMOs->optGrowthParam==2){
            //parameterization at min pre-molt size and ln-scale slope
            dvector pXDs = ptrGrw->getPCXDs(pc);
            zGrA = pXDs[1]; //pre-molt size corresponding to pGrA as mean post-molt size
            if (debug>dbgCalcProcs) cout<<"growth parameterization 2. "<<"zGrA:"<<tb<<zGrA<<endl;
            mnZs = grA*mfexp(grB*log(zBs/zGrA));
            if (debug>dbgCalcProcs) cout<<"mnZs:"<<tb<<mnZs<<endl;
        } else {
            //throw error
            PRINT2B1(" ")
            PRINT2B1("#---------------------")
            PRINT2B2("Invalid growth parameterization option ",ptrMOs->optGrowthParam)
            PRINT2B1("Please select a valid option and re-run.")
            ad_exit(-1);
        }
        dvar_vector mnIs = mnZs - zBs;              //mean molt increments
        dvariable invBeta = 1.0/grBeta;             //inverse scale for gamma density function
        dvar_vector alIs = mnIs*invBeta;            //gamma density alpha (location) parameters
        dvar_vector mnpIs = mnZs - (zBs - 2.5);     //mean molt increment (adjusted to start of size bin)
        dvar_vector alpIs = mnpIs*invBeta;          //gamma density alpha (location) parameters
        //check all mean molt increments are > 0
        if (isnan(value(sum(sqrt(mnIs))))){
            ofstream os("GrowthReport."+str(current_phase())+"."+str(ctrProcCallsInPhase)+".dat");
            os.precision(12);
            std::cout<<"##Found negative growth increments!"<<endl;
            std::cout<<"jitter seed = "<<iSeed<<endl;
            std::cout<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
            os<<"##Found negative growth increments!"<<endl;
            os<<"jitter seed = "<<iSeed<<endl;
            os<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
            os<<"pc: "<<pc<<tb<<"grA:"<<tb<<grA<<" grB:"<<tb<<grB<<" grBeta:"<<grBeta<<endl;
            os<<"zBs   = "<<zBs<<endl;
            os<<"mnZs  = "<<mnZs<<endl;
            os<<"mnIs  = "<<mnIs<<endl;
            os<<"sdIs  = "<<sqrt(mnIs*grBeta)<<endl;
            os<<"alIs  = "<<alIs<<endl;
            os.close();
            exit(-1);
        }
        prGr_zz.initialize();
        if (ptrMOs->optGrowthPDF==0) {
            //old style (TCSAM2013)
            for (int z=1;z<=nZBs;z++){//pre-molt growth bin
                dvector dpZs =  zBs(z,nZBs) - (zBs(z)-2.5);//realized growth increments (note non-neg. growth only)
                dvar_vector prs  = elem_prod(pow(dpZs,alpIs(z)-1.0),exp(-dpZs/grBeta)); //pr(dZ|z): use exp like TCSAM2013
                int mxIZ = min(nZBs,z+10);
                prGr_zz(z)(z,mxIZ) = prs(z,mxIZ);
                prGr_zz(z) /= sum(prGr_zz(z));
            }//zs
            prGr_zz(nZBs,nZBs) = 1.0; //no growth from max size
        } else if (ptrMOs->optGrowthPDF==1){
            //using cumd_gamma function like gmacs
            for (int z=1;z<nZBs;z++){
                dvar_vector sclIs = (ptrMC->zCutPts(z+1,nZBs+1)-zBs(z))/grBeta;//scaled increments at size bin cut points
                dvar_vector cprs(z,nZBs); cprs.initialize();
                dvar_vector prs(z,nZBs); prs.initialize();
                cprs(z) = cumd_gamma(sclIs(z+1),alIs(z));
                prs(z)  = cprs(z);
                for (int zp=z+1;zp<=nZBs;zp++){
                    cprs(zp) = cumd_gamma(sclIs(zp+1),alIs(z));
                    prs(zp)  = cprs(zp)-cprs(zp-1);//cumulative pr from zCs(zp) to zCs(zp+1)
                }
                prs(nZBs)   += 1.0 - cprs(nZBs);//treat final size bin as accumulator
                //cout<<"prs indices: "<<prs.indexmin()<<"  "<<prs.indexmax()<<endl;
                //check sum = 1
                if (sfabs(1.0-sum(prs))>1.0e-10){
                    std::cout<<"##Errors in calculating growth transition matrix: sum NE 1"<<endl;
                    std::cout<<"jitter seed = "<<iSeed<<endl;
                    std::cout<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
                    ofstream os("GrowthReportSum1."+str(current_phase())+"."+str(ctrProcCallsInPhase)+".dat");
                    os.precision(12);
                    os<<"Errors in calculating growth transition matrix: sum NE 1"<<endl;
                    os<<"jitter seed = "<<iSeed<<endl;
                    os<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
                    os<<"pc: "<<pc<<tb<<"grA:"<<tb<<grA<<" grB:"<<tb<<grB<<" grBeta:"<<grBeta<<endl;
                    os<<"zBs   = "<<zBs<<endl;
                    os<<"mnZs  = "<<mnZs<<endl;
                    os<<"mnIs  = "<<mnIs<<endl;
                    os<<"sdIs  = "<<sqrt(mnIs*grBeta)<<endl;
                    os<<"alIs  = "<<alIs<<endl;
                    os<<"z = "<<z<<tb<<"zB = "<<zBs(z)<<endl;
                    os<<"cutpts = "<<ptrMC->zCutPts(z+1,nZBs+1)<<endl;
                    os<<"sclIs = "<<sclIs<<endl;
                    os<<"sum(prs) = "<<sum(prs)<<endl;
                    os<<"prs  = "<<prs<<endl;
                    os<<"cprs = "<<cprs<<endl;
                    os.close();
                    exit(-1);
                }
                if (prs.size()>11) prs(z+11,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
                prs = prs/sum(prs);//normalize to sum to 1
                if (debug) cout<<prs<<endl;
                prGr_zz(z)(z,nZBs) = prs;
            }//zs
            prGr_zz(nZBs,nZBs) = 1.0; //no growth from max size
        } else {
            cout<<"Unrecognized growth option: "<<ptrMOs->optGrowthPDF<<endl;
            cout<<"Terminating!"<<endl;
            exit(-1);
        }
        mnGrZ_cz(pc) = mnZs;        
        prGr_czz(pc) = trans(prGr_zz);//transpose so rows are post-molt (i.e., lefthand z index is post-molt, or "to") z's so n+ = prGr_zz*n
        if (isnan(value(sum(prGr_zz)))){
            ofstream os("GrowthReportPrZZ."+str(current_phase())+"."+str(ctrProcCallsInPhase)+".dat");
            os.precision(12);
            std::cout<<"##Found NaN in prGz_zz in calcGrowth!"<<endl;
            std::cout<<"jitter seed = "<<iSeed<<endl;
            std::cout<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
            os<<"##Found NaN in prGz_zz in calcGrowth!"<<endl;
            os<<"jitter seed = "<<iSeed<<endl;
            os<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
            os<<"pc: "<<pc<<tb<<"grA:"<<tb<<grA<<" grB:"<<tb<<grB<<" grBeta:"<<grBeta<<endl;
            os<<"zBs   = "<<zBs<<endl;
            os<<"mnZs  = "<<mnZs<<endl;
            os<<"mnIs  = "<<mnIs<<endl;
            os<<"sdIs  = "<<sqrt(mnIs*grBeta)<<endl;
            os<<"alIs  = "<<alIs<<endl;
            if (ptrMOs->optGrowthPDF==0) {
                //old style (TCSAM2013)
                for (int z=1;z<nZBs;z++){//pre-molt growth bin
                    dvector dZs =  zBs(z,nZBs) - zBs(z);//realized growth increments (note non-neg. growth only)
                    os<<"Zpre = "<<zBs(z)<<endl;
                    os<<"dZs: "<<dZs<<endl;
                    os<<"Zbs: "<<zBs(z,nZBs)<<endl;
                    //dvar_vector prs = elem_prod(pow(dZs,alZ(z)-1.0),mfexp(-dZs/grBeta)); //pr(dZ|z)
                    dvar_vector prs = wts::log_gamma_density(dZs,alIs(z),invBeta);
                    os<<"log(prs): "<<prs<<endl;
                    prs = mfexp(prs);//gamma pdf
                    os<<"prs     : "<<prs<<endl;
                    if (prs.size()>11) prs(z+11,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
                    os<<"normalization factor = "<<sum(prs)<<endl;
                    os<<"prs     : "<<prs<<endl;
                }
            } else if (ptrMOs->optGrowthPDF==1){
                //using cumd_gamma function like gmacs
                for (int z=1;z<nZBs;z++){
                    dvar_vector sclIs = (ptrMC->zCutPts(z+1,nZBs+1)-zBs(z))/grBeta;//scaled increments at size bin cut points
                    dvar_vector prs(z,nZBs);
                    prs(z) = cumd_gamma(sclIs(z+1),alIs(z));
                    for (int zp=z+1;zp<=nZBs;zp++){
                        prs(zp) = cumd_gamma(sclIs(zp+1),alIs(z))-prs(zp-1);//cumulative pr from zCs(zp) to zCs(zp+1)
                    }
                    prs(nZBs) += 1.0 - cumd_gamma(sclIs(nZBs+1),alIs(z));//treat final size bin as accumulator
                    cout<<"zB = "<<zBs(z)<<tb<<"sum(prs) = "<<sum(prs)<<endl;
                    cout<<"prs = "<<prs<<endl;
                    if (prs.size()>11) prs(z+11,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
                    os<<"normalization factor = "<<sum(prs)<<endl;
                    os<<"prs     : "<<prs<<endl;
                }//zs
            }//optGrowth
            d3_array val = value(prGr_czz);
            os<<"prGr_czz = "<<endl; wts::print(val, os, 1); os<<endl;
            os.close();
            exit(-1);
            //testNaNs(value(sum(prGr_zz)),"Calculating growth");
        }
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrGrw->getModelIndices(pc);
        if (debug) cout<<"growth indices for pc "<<pc<<endl<<idxs<<endl;
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1); //year index
            if ((mnYr<=y)&&(y<=mxYr)){
                mnx = mxx = idxs(idx,2);//sex index
                if (mnx==tcsam::ALL_SXs){mnx = 1; mxx = tcsam::nSXs;}
                if (debug>dbgCalcProcs) cout<<"y = "<<y<<tb<<"sex = "<<tcsam::getSexType(mnx)<<": "<<tcsam::getSexType(mnx)<<endl;
                for (int x=mnx;x<=mxx;x++){
                    grA_xy(x,y) = grA;
                    grB_xy(x,y) = grB;
                    zGrA_xy(x,y) = zGrA;
                    zGrB_xy(x,y) = zGrB;
                    grBeta_xy(x,y)  = grBeta;
                    for (int s=1;s<=nSCs;s++){
                        mnGrZ_yxsz(y,x,s) = mnGrZ_cz(pc);
                        prGr_yxszz(y,x,s) = prGr_czz(pc);
                        //for (int z=1;z<=nZBs;z++) prGr_yxszz(y,x,s,z) = prGr_czz(pc,z);
                    }//s
                }//x
            }
        }//idx
    }//pc
    //set values for mxYr+1 to those for mxYr
    for (int x=1;x<=nSXs;x++){
        grA_xy(x,mxYr+1) = grA_xy(x,mxYr);
        grB_xy(x,mxYr+1) = grB_xy(x,mxYr);
        grBeta_xy(x,mxYr+1) = grBeta_xy(x,mxYr);
    }
    if (debug>dbgCalcProcs) cout<<"finished calcGrowth()"<<endl;
}

void model_parameters::calcSelectivities(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcSelectivities()"<<endl;
    SelectivityInfo* ptrSel = ptrMPI->ptrSel;
    double fsZ;    //fully selected size
    int idSel;     //selectivity function id
    int idxFSZ = 1;//index for fsZ in pXDs vector below
    ivector mniSelDevs(1,6);//min indices of devs vectors
    ivector mxiSelDevs(1,6);//max indices of devs vectors
    dvar_vector params(1,6);//vector for number_vector params
    dvar_vector paramsp(1,6);//vector for number_vector params + dev offsets
    npSel_cz.initialize();//nonparameteric selectivities
    sel_cz.initialize();//selectivities w/out deviations
    sel_cyz.initialize();//selectivity array
    int y;
    dvar_vector ptS1 = ptrSel->pS1->calcArithScaleVals(pS1);
    dvar_vector ptS2 = ptrSel->pS2->calcArithScaleVals(pS2);
    dvar_vector ptS3 = ptrSel->pS3->calcArithScaleVals(pS3);
    dvar_vector ptS4 = ptrSel->pS4->calcArithScaleVals(pS4);
    dvar_vector ptS5 = ptrSel->pS5->calcArithScaleVals(pS5);
    dvar_vector ptS6 = ptrSel->pS6->calcArithScaleVals(pS6);
    for (int pc=1;pc<=ptrSel->nPCs;pc++){
        ivector pids = ptrSel->getPCIDs(pc);
        dvector pXDs = ptrSel->getPCXDs(pc);
        fsZ   = pXDs[idxFSZ];
        idSel = pids[ptrSel->nIVs+ptrSel->nPVs+idxFSZ+1];
        if (debug) cout<<"pc = "<<pc<<tb<<"idSel = "<<idSel<<tb<<SelFcns::getSelFcnID(idSel)<<endl;
        if (idSel!=SelFcns::ID_NONPARAMETRIC){
            //extract the number parameters
            params.initialize();
            int k=ptrSel->nIVs+1;//1st parameter variable column
            if (pids[k]) {params[1] = ptS1(pids[k]);}   k++;
            if (pids[k]) {params[2] = ptS2(pids[k]);}   k++;
            if (pids[k]) {params[3] = ptS3(pids[k]);}   k++;
            if (pids[k]) {params[4] = ptS4(pids[k]);}   k++;
            if (pids[k]) {params[5] = ptS5(pids[k]);}   k++;
            if (pids[k]) {params[6] = ptS6(pids[k]);}   k++;
            if (debug>dbgCalcProcs) {
                cout<<"pc: "<<pc<<tb<<"pids = "<<pids<<endl;
                cout<<tb<<"params:"<<tb<<params<<endl;
            }
            int useDevsS1=pids[k++];
            dvar_vector dvsS1; ivector idxDevsS1;
            if (useDevsS1){
                dvsS1 = devsS1(useDevsS1);
                idxDevsS1 = idxsDevsS1(useDevsS1);
                if (debug>dbgCalcProcs){
                    cout<<"idx(dvsS1) = "<<idxDevsS1<<endl;
                    cout<<"dvsS1      = "<<dvsS1<<endl;
                }
            }
            int useDevsS2=pids[k++];
            dvar_vector dvsS2; ivector idxDevsS2;
            if (useDevsS2){
                dvsS2 = devsS2(useDevsS2);
                idxDevsS2 = idxsDevsS2(useDevsS2);
                if (debug>dbgCalcProcs){
                    cout<<"idx(dvsS2) = "<<idxDevsS2<<endl;
                    cout<<"dvsS2      = "<<dvsS2<<endl;
                }
            }
            int useDevsS3=pids[k++];
            dvar_vector dvsS3; ivector idxDevsS3;
            if (useDevsS3){
                dvsS3 = devsS3(useDevsS3);
                idxDevsS3 = idxsDevsS3(useDevsS3);
                if (debug>dbgCalcProcs){
                    cout<<"idx(dvsS3) = "<<idxDevsS3<<endl;
                    cout<<"dvsS3      = "<<dvsS3<<endl;
                }
            }
            int useDevsS4=pids[k++];
            dvar_vector dvsS4; ivector idxDevsS4;
            if (useDevsS4){
                dvsS4 = devsS4(useDevsS4);
                idxDevsS4 = idxsDevsS4(useDevsS4);
                if (debug>dbgCalcProcs){
                    cout<<"idx(dvsS4) = "<<idxDevsS4<<endl;
                    cout<<"dvsS4      = "<<dvsS4<<endl;
                }
            }
            int useDevsS5=pids[k++];
            dvar_vector dvsS5; ivector idxDevsS5;
            if (useDevsS5){
                dvsS5 = devsS5(useDevsS5);
                idxDevsS5 = idxsDevsS5(useDevsS5);
                if (debug>dbgCalcProcs){
                    cout<<"idx(dvsS5) = "<<idxDevsS5<<endl;
                    cout<<"dvsS5      = "<<dvsS5<<endl;
                }
            }
            int useDevsS6=pids[k++];
            dvar_vector dvsS6; ivector idxDevsS6;
            if (useDevsS6){
                dvsS6 = devsS6(useDevsS6);
                idxDevsS6 = idxsDevsS6(useDevsS6);
                if (debug>dbgCalcProcs){
                    cout<<"idx(dvsS6) = "<<idxDevsS6<<endl;
                    cout<<"dvsS6      = "<<dvsS6<<endl;
                }
            }
            if (debug>dbgCalcProcs) cout<<tb<<"fsZ: "<<fsZ<<tb<<"idSel"<<tb<<idSel<<tb<<SelFcns::getSelFcnID(idSel)<<tb<<params<<endl;
            //calc selectivity function WITHOUT any annual devs
            sel_cz(pc) = SelFcns::calcSelFcn(idSel, zBs, params, fsZ);
            if (debug>dbgCalcProcs) cout<<tb<<"pc = "<<pc<<tb<<"sel_cz: "<<sel_cz(pc)<<endl;
            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrSel->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                y = idxs(idx,1);//year
                if ((mnYr<=y)&&(y<=mxYr+1)){
                    paramsp = params;//set paramsp equal to base params
                    k=ptrSel->nIVs+1+6;//1st devs vector variable column
                    if (useDevsS1){
                        if (idxDevsS1[y]){
                            if (debug>dbgCalcProcs) cout<<tb<<idx<<tb<<y<<tb<<useDevsS1<<tb<<idxDevsS1[y]<<tb<<paramsp[1]<<tb<<devsS1(useDevsS1,idxDevsS1[y])<<endl;
                            paramsp[1] += devsS1(useDevsS1,idxDevsS1[y]);
                        }
                    }
                    if (useDevsS2) if (idxDevsS2[y]){paramsp[2] += devsS2(useDevsS2,idxDevsS2[y]);}
                    if (useDevsS3) if (idxDevsS3[y]){paramsp[3] += devsS3(useDevsS3,idxDevsS3[y]);}
                    if (useDevsS4) if (idxDevsS4[y]){paramsp[4] += devsS4(useDevsS4,idxDevsS4[y]);}
                    if (useDevsS5) if (idxDevsS5[y]){paramsp[5] += devsS5(useDevsS5,idxDevsS5[y]);}
                    if (useDevsS6) if (idxDevsS6[y]){paramsp[6] += devsS6(useDevsS6,idxDevsS6[y]);}
                    sel_cyz(pc,y) = SelFcns::calcSelFcn(idSel, zBs, paramsp, fsZ);
                    if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"paramsp: "<<paramsp<<tb<<"sel: "<<sel_cyz(pc,y)<<endl;
                } else {
                    if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"y outside model range--skipping year!"<<endl;
                }
            }//idx
        } else if (idSel==SelFcns::ID_NONPARAMETRIC){
            //calculate nonparametric selectivity function
            int pcNP = pids[ptrSel->nIVs+ptrSel->nPVs];
            dvar_vector nonParParams = pvNPSel(pcNP);
            int idZ = (int) fsZ;
            npSel_cz(pcNP) = SelFcns::calcSelFcn(idSel, zBs, nonParParams, idZ);
            sel_cz(pc)     = npSel_cz(pcNP);
            if (debug>dbgCalcProcs) cout<<tb<<"pcNP = "<<pcNP<<tb<<"pc = "<<pc<<tb<<"sel_cz: "<<sel_cz(pc)<<endl;
            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrSel->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                y = idxs(idx,1);//year
                if ((mnYr<=y)&&(y<=mxYr+1)) {
                    sel_cyz(pc,y) = npSel_cz(pcNP);
                    if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<endl<<"nonParParams: "<<paramsp<<endl<<"sel: "<<sel_cyz(pc,y)<<endl;
                } else {
                    if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"y outside model range--skipping year!"<<endl;
                }
            }
        }
    }//pc
    if (debug>dbgCalcProcs) cout<<"finished calcSelectivities()"<<endl;
}

void model_parameters::calcFisheryFs(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcFisheryFs()"<<endl;
    FisheriesInfo* ptrFsh = ptrMPI->ptrFsh;
    dvariable hm; //handling mortality
    dvariable lnC;//ln-scale capture rate
    dvariable arC;//arithmetic-scale capture rate
    dvsLnC_fy.initialize();
    for (int f=1;f<=nFsh;f++) idxDevsLnC_fy(f) = -1;
    hasF_fy.initialize();   //flags indicating whether or not fishery occurs
    hmF_fy.initialize();    //handling mortality
    cpF_fyxms.initialize(); //fully-selected capture rate
    sel_fyxmsz.initialize();//selectivity functions
    ret_fyxmsz.initialize();//retention functions
    cpF_fyxmsz.initialize();//size-specific capture rate
    rmF_fyxmsz.initialize();//retention rate
    dmF_fyxmsz.initialize();//discard mortality rate
    //tmF_yxmsz.initialize(); //total mortality rate
    /*********************************************************\n
     * Fully-selected annual capture rates are calculated     \n
     * using 2 approaches:                                    \n
     *  1. directly from parameter values (if useEX=0 below)  \n
     *  2. based on effort extrapolated to                    \n
     *      capture rate        (if useEX>1 below)            \n
     * Consequently, calculating all the above quantities     \n
     * requires 2 passes through parameter combinations.      \n
    ***********************************************************/
    int idxEX = ptrFsh->idxUseEX;//index into pids below for flag to use effort extrapolation
    int y, f, k, mnx, mxx, mnm, mxm, mns, mxs; 
    int idSel, idRet, useDevs;
    //Pass 1: calculations based on parameter values
    if (debug>dbgCalcProcs) cout<<"starting pass 1"<<endl;
    for (int pc=1;pc<=ptrFsh->nPCs;pc++){
        ivector pids = ptrFsh->getPCIDs(pc);
        if (debug>dbgCalcProcs) cout<<"pc: "<<pc<<tb<<"pids: "<<pids<<endl;
        int useEX = pids[ptrFsh->idxUseEX];//flag to use effort extrapolation
        if (!useEX){//calculate capture rates from parameters
            lnC.initialize();
            arC.initialize();
             //get handling mortality (default to 1)
            hm = 1.0;
            k = FisheriesInfo::idxHM;
            if (pids[k]) {hm = pHM(pids[k]);}
            //set base (ln-scale) capture rate
            k = FisheriesInfo::idxLnC;
            if (pids[k]) {lnC += pLnC(pids[k]);}
            //add in offset 1 (temporal, perhaps)
            k = FisheriesInfo::idxDC1;
            if (pids[k]) {lnC += pDC1(pids[k]);}
            //add in offset 2 (for females, perhaps)
            k = FisheriesInfo::idxDC2;
            if (pids[k]) {lnC += pDC2(pids[k]);}
            //add in offset 3 (for immature crab, perhaps)
            k = FisheriesInfo::idxDC3;
            if (pids[k]) {lnC += pDC3(pids[k]);}
            //add in offset 4 (for immature females, perhaps)
            k = FisheriesInfo::idxDC4;
            if (pids[k]) {lnC += pDC4(pids[k]);}
            //extract devs vector
            k = FisheriesInfo::idxDevsLnC;
            useDevs = pids[k];
            dvar_vector dvsLnC;             
            ivector idxDevsLnC;
            if (useDevs) {
                dvsLnC     = devsLnC(useDevs);
                idxDevsLnC = idxsDevsLnC(useDevs);
            } else {
                arC = mfexp(lnC);
            }
            //extract logistic scale retention fraction for old shell crab
            k = FisheriesInfo::idxLgtRet;
            dvariable retFrac = 1.0;
            if (pids[k]) {
                retFrac = 1.0/(1.0+mfexp(-pLgtRet[pids[k]]));
                if (debug>dbgCalcProcs) {
                    cout<<"pc: "<<pc<<". retFrac = "<<retFrac<<endl;
                }
            } k++;
            idSel = pids[FisheriesInfo::idxSelFcn];//selectivity function id
            idRet = pids[FisheriesInfo::idxRetFcn];//retention function id
            //convert from ln-scale to arithmetic scale
            if (debug>dbgCalcProcs){
                cout<<"pc: "<<pc<<". idSel = "<<idSel<<". idRet = "<<idRet<<". lnC = "<<lnC<<". retFrac = "<<retFrac<<endl;
                if (useDevs) {
                    cout<<tb<<tb<<"dvsLnC["<<dvsLnC.indexmin()<<cc<<dvsLnC.indexmax()<<"] = "<<dvsLnC<<endl;
                } else {
                    cout<<tb<<tb<<"arC:"<<arC<<endl;
                }
            }
            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrFsh->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                f = idxs(idx,1);//fishery
                y = idxs(idx,2);//year
                if ((mnYr<=y)&&(y<=mxYr)){
                    hasF_fy(f,y) = 1;//flag indicating occurrence of fishery in year y
                    hmF_fy(f,y) = hm;//save discard mortality rate
                    mnx = mxx = idxs(idx,3);//sex index
                    if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
                    mnm = mxm = idxs(idx,4);//maturity index
                    if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
                    mns = mxs = idxs(idx,5);//shell condition index
                    if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
                    for (int x=mnx;x<=mxx;x++){
                        if (debug>dbgCalcProcs) cout<<"f,y,x,useDevs = "<<f<<cc<<y<<cc<<x<<cc<<useDevs<<endl;
                        if (useDevs) {
                            idxDevsLnC_fy(f,y) = idxDevsLnC[y];
                            dvsLnC_fy(f,y)     = dvsLnC[idxDevsLnC[y]];
                            arC = mfexp(lnC+dvsLnC[idxDevsLnC[y]]);//recalculate arC w/ devs
                            if (debug>dbgCalcProcs) {
                                cout<<"idxDevsLnC[y],dvsLnC[idxDevsLnC[y]], arC: "<<idxDevsLnC[y]<<tb<<dvsLnC[idxDevsLnC[y]]<<tb<<arC<<endl;
                            }
                        }
                        for (int m=mnm;m<=mxm;m++){
                            cpF_fyxms(f,y,x,m)  = arC; //fully-selected capture rate (independent of shell condition)
                            for (int s=mns;s<=mxs;s++){
                                sel_fyxmsz(f,y,x,m,s) = sel_cyz(idSel,y);          //selectivity
                                cpF_fyxmsz(f,y,x,m,s) = cpF_fyxms(f,y,x,m,s)*sel_cyz(idSel,y);//size-specific capture rate
                                if (idRet){//fishery has retention
                                    ret_fyxmsz(f,y,x,m,s) = retFrac*sel_cyz(idRet,y);      //retention curves
                                    rmF_fyxmsz(f,y,x,m,s) = elem_prod(ret_fyxmsz(f,y,x,m,s),         cpF_fyxmsz(f,y,x,m,s));//retention mortality rate
                                    dmF_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-ret_fyxmsz(f,y,x,m,s)),cpF_fyxmsz(f,y,x,m,s));//discard mortality rate
                                } else {//discard only
                                    dmF_fyxmsz(f,y,x,m,s) = hm*cpF_fyxmsz(f,y,x,m,s);//discard mortality rate
                                }
                            }//s
                        }//m
                    }//x
                }//(mnYr<=y)&&(y<=mxYr)
            }//idx
        }//useEX=FALSE
    }//pc
    if (debug>dbgCalcProcs) cout<<"finished pass 1"<<endl;
    //calculate ratio of average capture rate to effort
    if (debug>dbgCalcProcs) cout<<"--calculating average capture rates"<<endl;
    avgFc_nxms.initialize();
    avgFc2Eff_nxms.initialize();
    obsFc_nxmsy.initialize();
    prdFc_nxmsy.initialize();
    prdEff_nxmsy.initialize();
    CapRateAvgScenarios* pCRASs = ptrMOs->ptrEffXtrapScenarios->ptrCapRateAvgScenarios;
    int nCRASs = pCRASs->nAvgs;
    for (int n=1;n<=nCRASs;n++){//capture rate averaging scenarios
        CapRateAvgScenario* ptrCRAS = pCRASs->ppCRASs[n-1];
        int idEAS = ptrCRAS->idEffAvgInfo;//index to associated average effort
        int idPar = ptrCRAS->idParam;     //index to associated extrapolation parameter
        int f     = ptrCRAS->f;  //fishery
        int fd    = mapM2DFsh(f);//index of corresponding fishery data object
        dvector eff_y = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y;//corresponding effort time series
        mnx = mxx = ptrCRAS->x;//sex index
        if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
        mnm = mxm = ptrCRAS->m;//maturity index
        if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
        mns = mxs = ptrCRAS->s;//shell index
        if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
        for (int x=mnx;x<=mxx;x++){
            for (int m=mnm;m<=mxm;m++){
                for (int s=mns;s<=mxs;s++){
                    if (debug>dbgCalcProcs)cout<<"capture rate averaging for n="<<n<<"idEAS="<<idEAS<<cc<<" avgEff_n(idEAS)="<<avgEff_n(idEAS)<<cc
                                               <<" f="<<f<<cc<<" x="<<x<<cc<<" m="<<m<<cc<<" s="<<s<<endl;
                    //loop over years to extract "observed" F's and calculate averages
                    for (int iy=yrsAvgEff_ny(idEAS).indexmin();iy<=yrsAvgEff_ny(idEAS).indexmax();iy++){
                        //cout<<iy<<cc<<yrsAvgEff_ny(idEAS,iy)<<endl;
                        obsFc_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy)) = cpF_fyxms(f,yrsAvgEff_ny(idEAS,iy),x,m,s);
                    }
                    avgFc_nxms(n,x,m,s) = sum(obsFc_nxmsy(n,x,m,s))/yrsAvgEff_ny(idEAS).size();
                    avgFc2Eff_nxms(n,x,m,s) = avgFc_nxms(n,x,m,s)/avgEff_n(idEAS);
                    //loop over year again to calculate "predicted" F's based on effort
                    for (int iy=yrsAvgEff_ny(idEAS).indexmin();iy<=yrsAvgEff_ny(idEAS).indexmax();iy++){
                        //fully-selected capture rate
                        if (idPar==0){
                            //extrapolation based on effort extrapolation ratio only
                            prdFc_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy)) =                      
                                                          avgFc2Eff_nxms(n,x,m,s)*eff_y(yrsAvgEff_ny(idEAS,iy)); 
                            //predicted effort from "observed" capture rates
                            prdEff_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy)) = 
                                    obsFc_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy))/avgFc2Eff_nxms(n,x,m,s);
                        } else {
                            //extrapolation based on effort extrapolation parameters, as well as ratio
                            prdFc_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy)) = 
                                    mfexp(pLnEffX(idPar))*avgFc2Eff_nxms(n,x,m,s)*eff_y(yrsAvgEff_ny(idEAS,iy));
                            //predicted effort from "observed" capture rates
                            prdEff_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy)) = 
                                    obsFc_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy))/(mfexp(pLnEffX(idPar))*avgFc2Eff_nxms(n,x,m,s));
                        }
                    }
                    if (debug>dbgCalcProcs){
                        cout<<"yrsAvgEff_ny(idEAS)     = "<<yrsAvgEff_ny(idEAS)<<endl;
                        cout<<"avgFc_nxms(n,x,m,s)     = "<<avgFc_nxms(n,x,m,s)<<endl;
                        cout<<"avgFc2Eff_nxms(n,x,m,s) = "<<avgFc2Eff_nxms(n,x,m,s)<<endl;
                        cout<<"obsFc_nxmsy(n,x,m,s)    = "<<obsFc_nxmsy(n,x,m,s)<<endl;
                        cout<<"prdFc_nxmsy(n,x,m,s)    = "<<prdFc_nxmsy(n,x,m,s)<<endl;
                        cout<<"obsEff_nxmsy(n,x,m,s)   = "<<obsEff_nxmsy(n,x,m,s)<<endl;
                        cout<<"prdEff_nxmsy(n,x,m,s)   = "<<prdEff_nxmsy(n,x,m,s)<<endl;
                    }
                }//s
            }//m
        }//x
    }//n - capture rate averaging
    if (debug>dbgCalcProcs) cout<<"calculated avgFc2Eff_nxms"<<endl;
    //Pass 2: calculations based on effort and effort extrapolation parameters or average effort:capture rate ratios
    if (debug>dbgCalcProcs) cout<<"Starting pass 2"<<endl;
    int fd; double eff;
    for (int pc=1;pc<=ptrFsh->nPCs;pc++){
        ivector pids = ptrFsh->getPCIDs(pc);
        int useEX   = pids[FisheriesInfo::idxUseEX];  //flag to use direct effort extrapolation (+ index to EX scenario [i.e., "n"])
        int idPar   = pids[FisheriesInfo::idxLnEffX]; //index to parameters for effort extrapolation
        if (useEX){//calculate capture rates from effort extrapolation
            //get handling mortality (default to 1)
            hm = 1.0;
            if (pids[ptrFsh->idxHM]) {hm = pHM(pids[ptrFsh->idxHM]);}
            //extract logistic scale retention fraction for old shell crab
            dvariable retFrac = 1.0;
            k = FisheriesInfo::idxLgtRet;
            if (pids[k]) {
                retFrac = 1.0/(1.0+mfexp(-pLgtRet[pids[k]]));
                if (debug>dbgCalcProcs) {
                    cout<<"pc: "<<pc<<". retFrac = "<<retFrac<<endl;
                }
            }
            idSel = pids[FisheriesInfo::idxSelFcn];//selectivity function id
            idRet = pids[FisheriesInfo::idxRetFcn];//retention function id
            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrFsh->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                f = idxs(idx,1);//fishery
                y = idxs(idx,2);//year
                if ((mnYr<=y)&&(y<=mxYr)){
                    hasF_fy(f,y) = 1;//flag indicating occurrence of fishery in year y
                    hmF_fy(f,y) = hm;//save discard mortality rate
                    fd = mapM2DFsh(f);//index of corresponding fishery data object
                    eff = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y(y);
                    mnx = mxx = idxs(idx,3);//sex index
                    if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
                    mnm = mxm = idxs(idx,4);//maturity index
                    if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
                    mns = mxs = idxs(idx,5);//shell index
                    if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
                    if (debug>dbgCalcProcs) cout<<"f, y n, x, m, s, pLnEffX avgFc2Eff, eff cpF  = "<<endl;
                    for (int x=mnx;x<=mxx;x++){
                        for (int m=mnm;m<=mxm;m++){
                            for (int s=mns;s<=mxs;s++){
                                //fully-selected capture rate
                                if (idPar==0){
                                        //extrapolation based on effort extrapolation ratio only
                                    //prdFc_nxmsy(useEX,x,m,s,y) = avgFc2Eff_nxms(useEX,x,m,s)*eff; 
                                    //cpF_fyxms(f,y,x,m,s)   = prdFc_nxmsy(useEX,x,m,s,y); 
                                    cpF_fyxms(f,y,x,m,s) = avgFc2Eff_nxms(useEX,x,m,s)*eff; 
                                } else {
                                    //extrapolation based on effort extrapolation parameters, as well as ratio
                                    //prdFc_nxmsy(useEX,x,m,s,y) = mfexp(pLnEffX(idPar))*avgFc2Eff_nxms(useEX,x,m,s)*eff;
                                    //cpF_fyxms(f,y,x,m,s) = prdFc_nxmsy(useEX,x,m,s,y);
                                    cpF_fyxms(f,y,x,m,s) = mfexp(pLnEffX(idPar))*avgFc2Eff_nxms(useEX,x,m,s)*eff;
                                }
                                if (debug>dbgCalcProcs) {
                                    if (idPar==0) cout<<f<<tb<<y<<useEX<<tb<<x<<tb<<0             <<avgFc2Eff_nxms(useEX,x,m,s)<<tb<<eff<<tb<<cpF_fyxms(f,y,x,m,s)<<endl;
                                    if (idPar)    cout<<f<<tb<<y<<useEX<<tb<<x<<tb<<pLnEffX(idPar)<<avgFc2Eff_nxms(useEX,x,m,s)<<tb<<eff<<tb<<cpF_fyxms(f,y,x,m,s)<<endl;
                                }
                                if (!debug) testNaNs(value(cpF_fyxms(f,y,x,m,s)),"calcFisheryFs: 2nd pass");
                                sel_fyxmsz(f,y,x,m,s) = sel_cyz(idSel,y);
                                cpF_fyxmsz(f,y,x,m,s) = cpF_fyxms(f,y,x,m,s)*sel_cyz(idSel,y);//size-specific capture rate
                                if (idRet){//fishery has retention
                                    ret_fyxmsz(f,y,x,m,s) = retFrac*sel_cyz(idRet,y);
                                    rmF_fyxmsz(f,y,x,m,s) = elem_prod(ret_fyxmsz(f,y,x,m,s),         cpF_fyxmsz(f,y,x,m,s));//retention mortality rate
                                    dmF_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-ret_fyxmsz(f,y,x,m,s)),cpF_fyxmsz(f,y,x,m,s));//discard mortality rate
                                } else {//discard only
                                    dmF_fyxmsz(f,y,x,m,s) = hm*cpF_fyxmsz(f,y,x,m,s);//discard mortality rate
                                }
                            }//s
                        }//m
                    }//x
                }//(mnYr<=y)&&(y<=mxYr)
            }//idx
        }//useEX=TRUE
    }//pc
    if (debug>dbgCalcProcs) cout<<"Finished pass 2"<<endl;
    if (debug>dbgCalcProcs) {
        for (int f=1;f<=nFsh;f++){
            cout<<"cpF_fyxmsz("<<f<<",mnYr:mxYr,  MALE,MATURE,NEW SHELL) = ";
            for (int y=mnYr;y<=mxYr;y++) {cout<<cpF_fyxmsz(f,y,  MALE,MATURE,NEW_SHELL)<<tb;} cout<<endl;
            cout<<"cpF_fyxmsz("<<f<<",mnYr:mxYr,FEMALE,MATURE,NEW SHELL) = ";
            for (int y=mnYr;y<=mxYr;y++) {cout<<cpF_fyxmsz(f,y,FEMALE,MATURE,NEW_SHELL)<<tb;} cout<<endl;
        }
        cout<<"finished calcFisheryFs()"<<endl;
    }
}

void model_parameters::calcSurveyQs(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcSurveyQs()"<<endl;
    SurveysInfo* ptrSrv = ptrMPI->ptrSrv;
    dvariable lnQ;//ln-scale Q's
    dvariable arQ;//arithmetic-scale Q's
    q_vyxms.initialize();
    a_vyxmsz.initialize();
    q_vyxmsz.initialize();
    int y; int v; int mnx; int mxx; int mnm; int mxm; int mns; int mxs;
    int idAvl; int idSel;
    dvar_vector ptQ  = ptrSrv->pQ->calcArithScaleVals(pQ);
    for (int pc=1;pc<=ptrSrv->nPCs;pc++){
        lnQ.initialize();
        arQ.initialize();
        ivector pids = ptrSrv->getPCIDs(pc);
        if (debug>dbgCalcProcs) cout<<"pc: "<<pc<<tb<<"pids = "<<pids<<endl;
        int k=ptrSrv->nIVs+1;//1st parameter variable column
        //add in base catchability (e.g., for mature male)
        if (pids[k]) {lnQ += log(ptQ(pids[k]));} k++;
        //add in ln-scale offset 1            (e.g., for time period)
        if (pids[k]) {lnQ += pDQ1(pids[k]);} k++;
        //add in ln-scale offset 2            (e.g., for females)
        if (pids[k]) {lnQ += pDQ2(pids[k]);} k++;
        //add in ln-scale offset 3            (e.g., for immature crab)
        if (pids[k]) {lnQ += pDQ3(pids[k]);} k++;
        //add in ln-scale offset 4            (e.g., for immature females)
        if (pids[k]) {lnQ += pDQ4(pids[k]);} k++; 
        idAvl = pids[k++];//availability function id
        idSel = pids[k++];//selectivity function id
        //convert from ln-scale to arithmetic scale
        arQ = mfexp(lnQ);
        if (debug>dbgCalcProcs){
            cout<<"lnQ:"<<lnQ<<tb<<"arQ:"<<endl<<arQ<<endl;
        }
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrSrv->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            v = idxs(idx,1);//survey
            y = idxs(idx,2);//year
            if ((mnYr<=y)&&(y<=mxYrp1)){
                mnx = mxx = idxs(idx,3);//sex index
                if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
                mnm = mxm = idxs(idx,4);//maturity state index
                if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
                mns = mxs = idxs(idx,5);//shell condition index
                if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
                for (int x=mnx;x<=mxx;x++){
                    for (int m=mnm;m<=mxm;m++){
                        q_vyxms(v,y,x,m) = arQ;
                        for (int s=mns;s<=mxs;s++){
                            a_vyxmsz(v,y,x,m,s) = 1.0;//default: availability = 1
                            if (idAvl>0) a_vyxmsz(v,y,x,m,s) = sel_cyz(idAvl,y);
                            if (idSel>0) s_vyxmsz(v,y,x,m,s) = sel_cyz(idSel,y);
                            q_vyxmsz(v,y,x,m,s) = arQ*elem_prod(a_vyxmsz(v,y,x,m,s),s_vyxmsz(v,y,x,m,s));
                        }//s
                    }//m
                }//x
            }//(mnYr<=y)&&(y<=mxYrp1)
        }//idx
    }//pc
    if (debug>dbgCalcProcs) cout<<"finished calcSurveyQs()"<<endl;
}

void model_parameters::calcEqNatZF100(dvariable& R, int yr, int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting void calcEqNatZF100(R,yr)"<<endl;
    n_xmsz.initialize();//equilibrium n-at-z
    for (int x=1;x<=nSXs;x++){
        S1_msz.initialize(); //survival until molting/mating
        Th_sz.initialize();  //pr(molt to maturity|pre-molt size, molt)
        T_szz.initialize();  //growth matrices (indep. of molt to maturity)
        S2_msz.initialize(); //survival after molting/mating
        R_z.initialize();    //recruitment size distribution
        R_z = R*R_yx(yr,x)*R_yz(yr);//initial mean recruitment by size
        for (int s=1;s<=nSCs;s++){
            Th_sz(s) = prM2M_yxz(yr,x); //pr(molt to maturity|pre-molt size, molt)
            for (int z=1;z<=nZBs;z++) T_szz(s,z) = prGr_yxszz(yr,x,s,z);//growth matrices
            for (int m=1;m<=nMSs;m++){ 
                S1_msz(m,s) = mfexp(-M_yxmsz(yr,x,m,s)*dtM_y(yr));      //survival until molting/growth/mating
                S2_msz(m,s) = mfexp(-M_yxmsz(yr,x,m,s)*(1.0-dtM_y(yr)));//survival after molting/growth/mating
            }//m
        }//s
        n_xmsz(x) = calcEqNatZ(R_z, S1_msz, Th_sz, T_szz, S2_msz, debug, cout);
    }
    if (debug>=dbgPopDy) cout<<"finished void calcEqNatZF100(R,yr)"<<endl;
}

dvar3_array model_parameters::calcEqNatZ(dvar_vector& R_z,dvar3_array& S1_msz, dvar_matrix& Th_sz, dvar3_array& T_szz, dvar3_array& S2_msz, int debug, ostream& cout)
{
    RETURN_ARRAYS_INCREMENT();
    if (debug>=dbgPopDy) cout<<"starting dvar3_array calcEqNatZ(...)"<<endl;
    //the equilibrium solution
    dvar3_array n_msz(1,nMSs,1,nSCs,1,nZBs); n_msz.initialize();
    //create an identity matrix
    dmatrix I = identity_matrix(1,nZBs);
    //--calc the state transition matrices
    int i = IMMATURE; 
    int m =   MATURE;
    int n = NEW_SHELL;
    int o = OLD_SHELL;
    //immature new shell crab
    dvar_matrix S2_in = wts::diag(S2_msz(i,n)); //pr(survival|size) for immature new shell crab after molting occurs
    dvar_matrix Tr_in = T_szz(n);               //pr(Z_post|Z_pre, molt) for immature new shell crab pre-terminal molt
    dvar_matrix Th_in = wts::diag(Th_sz(n));    //pr(molt to maturity|pre-molt size,new shell, molting)
    dvar_matrix Ph_in = identity_matrix(1,nZBs);//pr(molt|pre-molt size, new shell) [assumed that all immature crab molt]
    dvar_matrix S1_in = wts::diag(S1_msz(i,n)); //pr(survival|size) for immature new shell crab before molting occurs
    //immature old shell crab [shouldn't be any of these]
    dvar_matrix S2_io = wts::diag(S2_msz(i,o)); //pr(survival|size) for immature old shell crab after molting occurs (but they didn't molt)
    dvar_matrix Tr_io = T_szz(o);               //pr(Z_post|Z_pre, molt) for immature old shell crab pre-terminal molt
    dvar_matrix Th_io = wts::diag(Th_sz(o));    //pr(molt to maturity|pre-molt size,old shell, molting)
    dvar_matrix Ph_io = identity_matrix(1,nZBs);//pr(molt|pre-molt size, old shell) [assumed all immature crab molt]
    dvar_matrix S1_io = wts::diag(S1_msz(i,o)); //pr(survival|size) for immature old shell crab before molting occurs
    //mature new shell crab
    dvar_matrix Tr_mn = T_szz(n);               //pr(Z_post|Z_pre, molt) for immature new shell crab undergoing terminal molt (same as non-terminal molt)
    dvar_matrix Tr_mo = T_szz(o);               //pr(Z_post|Z_pre, molt) for immature old shell crab undergoing terminal molt (same as non-terminal molt)
    dvar_matrix S2_mn = wts::diag(S2_msz(m,n)); //pr(survival|size) for mature new shell crab after molting occurs
    dvar_matrix S1_mn = wts::diag(S1_msz(m,n)); //pr(survival|size) for mature new shell crab before molting occurs (but they won't molt)
    //mature old shell crab
    dvar_matrix S2_mo = wts::diag(S2_msz(m,o)); //pr(survival|size) for mature old shell crab after molting occurs (but they didn't molt)
    dvar_matrix S1_mo = wts::diag(S1_msz(m,o)); //pr(survival|size) for mature old shell crab before molting occurs (but they won't molt)
    //full state transition matrices
    dvar_matrix lA = S2_in * (I-Th_in) * Tr_in * Ph_in * S1_in;//imm, new -> imm, new
    dvar_matrix lB = S2_in * (I-Th_io) * Tr_io * Ph_io * S1_io;//imm, old -> imm, new
    dvar_matrix lC = S2_io * (I-Ph_in) * S1_in;                //imm, new -> imm, old
    dvar_matrix lD = S2_io * (I-Ph_io) * S1_io;                //imm, old -> imm, old
    dvar_matrix lE = S2_mn * Th_in * Tr_mn * Ph_in * S1_in;    //imm, new -> mat, new (terminal molt)
    dvar_matrix lF = S2_mn * Th_io * Tr_mo * Ph_io * S1_io;    //imm, old -> mat, new (terminal molt)
    dvar_matrix lG = S2_mo * S1_mn;                            //mat, new -> mat, old
    dvar_matrix lH = S2_mo * S1_mo;                            //mat, old -> mat, old
    //--done calculating transition matrices
    //calculate inverses of matrix quantities
    dvar_matrix iM1 = inv(I - lD);
    dvar_matrix iM2 = inv(I - lA - lB * iM1 * lC);
    dvar_matrix iM3 = inv(I - lH);
    //the equilibrium solution is
    n_msz(i,n) = iM2 * R_z;                         //immature, new shell
    n_msz(i,o) = iM1 * lC * n_msz(i,n);             //immature, old shell
    n_msz(m,n) = lE * n_msz(i,n) + lF * n_msz(i,o); //  mature, new shell
    n_msz(m,o) = iM3 * lG * n_msz(m,n);             //  mature, old shell
    if (debug>=dbgPopDy) cout<<"finished dvar3_array calcEqNatZ(...)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return(n_msz);
}

d5_array model_parameters::calcCohortProgression(int yr, int debug, ostream& cout)
{
    if (debug) {
        cout<<endl<<endl<<"#------------------------"<<endl;
        cout<<"starting calcCohortProgression(yr,debug,cout)"<<endl;
        cout<<"year for progression = "<<yr<<endl;
    }
    int nyp = 20;//number of years to track cohorts
    //1. set initial cohort abundance
    int nzp = 3;//number of size bins to initially distribute cohort over
    dvar5_array n_yxmsz(0,nyp,1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    n_yxmsz.initialize();
    n_yxmsz(0,  MALE,IMMATURE,NEW_SHELL)(1,nzp) = R_yz(yr)(1,nzp)/sum(R_yz(yr)(1,nzp));//set initial abundance-at-size
    n_yxmsz(0,FEMALE,IMMATURE,NEW_SHELL)(1,nzp) = R_yz(yr)(1,nzp)/sum(R_yz(yr)(1,nzp));//set initial abundance-at-size
    //2. Determine population rates, based on yr
    double dtF = dtF_y(yr);//time at which fisheries occur
    double dtM = dtM_y(yr);//time at which mating occurs
    PopDyInfo* pPIM = new PopDyInfo(nZBs);//  males info
    pPIM->R_z   = R_yz(yr);
    pPIM->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(MALE);
    pPIM->M_msz = M_yxmsz(yr,MALE);
    pPIM->T_szz = prGr_yxszz(yr,MALE);
    for (int s=1;s<=nSCs;s++) pPIM->Th_sz(s) = prM2M_yxz(yr,MALE);
    PopDyInfo* pPIF = new PopDyInfo(nZBs);//females info
    pPIF->R_z   = R_yz(yr);
    pPIF->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(FEMALE);
    pPIF->M_msz = M_yxmsz(yr,FEMALE);
    pPIF->T_szz = prGr_yxszz(yr,FEMALE);
    for (int s=1;s<=nSCs;s++) pPIF->Th_sz(s) = prM2M_yxz(yr,FEMALE);
    if (debug) cout<<"calculated pPIM, pPIF."<<endl;
    //3. Determine fishery conditions based on yr
    int avgPeriodYrs = 1;  
    //assumption here is that ALL fisheries EXCEPT the first are bycatch fisheries
    //a. Calculate average handling mortality, retention curves and capture rates
    int ny;   //number of years fishery is active
    dvar_vector avgHM_f(1,nFsh);
    avgHM_f.initialize();
    for (int f=1;f<=nFsh;f++){
        ny = 0;
        for (int y=yr-avgPeriodYrs+1;y<=yr;y++){
            ny         += hasF_fy(f,y);
            avgHM_f(f) += hmF_fy(f,y);
        }
        avgHM_f(f) /= wts::max(1.0,1.0*ny);
    }
    if (debug) cout<<"avgHm_f = "<<avgHM_f<<endl;
    dvar5_array avgCapF_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery capture rates
    dvar5_array avgRFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery retention functions
    dvar5_array avgSFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery selectivity functions
    avgCapF_xfmsz.initialize();
    avgRFcn_xfmsz.initialize();
    avgSFcn_xfmsz.initialize();
    for (int f=1;f<=nFsh;f++){
        avgPeriodYrs = 1;
        if (debug) cout<<"avgPeriodYrs("<<f<<") = "<<avgPeriodYrs<<endl;
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++) {
                    for (int z=1;z<=nZBs;z++){
                        ny = 0;
                        for (int y=(yr-avgPeriodYrs+1);y<=yr;y++) {
                            //if (debug) cout<<"y = "<<y<<endl;
                            ny += hasF_fy(f,y);
                            avgCapF_xfmsz(x,f,m,s,z) += cpF_fyxmsz(f,y,x,m,s,z);
                            avgRFcn_xfmsz(x,f,m,s,z) += ret_fyxmsz(f,y,x,m,s,z);
                            avgSFcn_xfmsz(x,f,m,s,z) += sel_fyxmsz(f,y,x,m,s,z);
                        }//y
                        double fac = 1.0/max(1.0,1.0*ny);
                        avgCapF_xfmsz(x,f,m,s,z) *= fac;
                        avgRFcn_xfmsz(x,f,m,s,z) *= fac;
                        avgSFcn_xfmsz(x,f,m,s,z) *= fac;
                    }//z
                }//s
            }//m
        }//x
    }//f
    if (debug){
        for (int f=1;f<=nFsh;f++){
            cout<<"avgCapF_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
            cout<<"avgCapF_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
            cout<<"avgRFcn_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
            cout<<"avgRFcn_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
            cout<<"avgSFcn_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgSFcn_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
            cout<<"avgSFcn_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgSFcn_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
        }
    }
    CatchInfo* pCIM = new CatchInfo(nZBs,nFsh);//male catch info
    pCIM->setCaptureRates(avgCapF_xfmsz(MALE));
    pCIM->setRetentionFcns(avgRFcn_xfmsz(MALE));
    pCIM->setHandlingMortality(avgHM_f);
    dvariable maxCapF = pCIM->findMaxTargetCaptureRate(cout);
    if (debug) cout<<"maxCapF = "<<maxCapF<<endl;
    CatchInfo* pCIF = new CatchInfo(nZBs,nFsh);//female catch info
    pCIF->setCaptureRates(avgCapF_xfmsz(FEMALE));
    pCIF->setRetentionFcns(avgRFcn_xfmsz(FEMALE));
    pCIF->setHandlingMortality(avgHM_f);
    pCIF->maxF = maxCapF;//need to set this for females
    //4. Create PopProjectors
    PopProjector* pPPM = new PopProjector(pPIM,pCIM);
    pPPM->dtF = dtF;
    pPPM->dtM = dtM;
    PopProjector* pPPF = new PopProjector(pPIF,pCIF);
    if (debug) cout<<"created pPPs."<<endl;
    pPPF->dtF = dtF;
    pPPF->dtM = dtM;
    //5. Create multi-year population projectors
    MultiYearPopProjector* pMYPPM = new MultiYearPopProjector(pPPM);
    MultiYearPopProjector* pMYPPF = new MultiYearPopProjector(pPPF);
    if (debug) cout<<"created pMPPs."<<endl;
    //6.project with no recruitment
    pMYPPM->project(nyp,0.0,maxCapF,n_yxmsz(0,  MALE),cout);
    pMYPPF->project(nyp,0.0,maxCapF,n_yxmsz(0,FEMALE),cout);
    if (debug) cout<<"projected cohort."<<endl;
    for (int y=0;y<=nyp;y++){
        n_yxmsz(y,  MALE) = pMYPPM->n_ymsz(y);
        n_yxmsz(y,FEMALE) = pMYPPF->n_ymsz(y);
    }
    d5_array vn_yxmsz = wts::value(n_yxmsz);
    cout<<"finished calcCohortProgression(yr,debug,cout)"<<endl<<endl<<endl;
    return vn_yxmsz;
}

void model_parameters::calcOFL(int yr, int debug, ostream& cout)
{
    if (debug) {
        cout<<endl<<endl<<"#------------------------"<<endl;
        cout<<"starting calcOFL(yr,debug,cout)"<<endl;
        cout<<"year for projection = "<<yr<<endl;
    }
    //1. get initial population for "upcoming" year, yr
    dvar4_array n_xmsz = n_yxmsz(yr);
    if (debug) {cout<<"  males_msz:"<<endl; wts::print(n_xmsz(  MALE),cout,1);}
    if (debug) {cout<<"females_msz:"<<endl; wts::print(n_xmsz(FEMALE),cout,1);}
    //2. set yr back one year to get population rates, etc., 
    //   from year prior to projection year
    yr = yr-1;//don't have pop rates, etc. for projection year
    if (debug) cout<<"year for pop rates = "<<yr<<endl;
    //3. Determine mean recruitment 
    //   1981 here corresponds to 1982 in TCSAM2013, the year recruitment enters
    //   the model population.
    dvar_vector avgRec_x(1,nSXs);
    if (debug) cout<<"R dims: "<<R_y.indexmin()<<cc<<R_y.indexmax()<<endl;
    for (int x=1;x<=nSXs;x++) 
        avgRec_x(x)= mean(elem_prod(R_y(1981,yr),column(R_yx,x)(1981,yr)));
    if (debug) {
        cout<<"R_y(  1981:"<<yr<<")      = "<<R_y(1981,yr)<<endl;
        cout<<"R_yx((1981:"<<yr<<",MALE) = "<<column(R_yx,MALE)(1981,yr)<<endl;
        cout<<"Average recruitment = "<<avgRec_x<<endl;
    }
    //4. Determine population rates for next year, using yr
    double dtF = dtF_y(yr);//time at which fisheries occur
    double dtM = dtM_y(yr);//time at which mating occurs
    PopDyInfo* pPIM = new PopDyInfo(nZBs);//  males info
    pPIM->R_z   = R_yz(yr);
    pPIM->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(MALE);
    pPIM->M_msz = M_yxmsz(yr,MALE);
    pPIM->T_szz = prGr_yxszz(yr,MALE);
    for (int s=1;s<=nSCs;s++) pPIM->Th_sz(s) = prM2M_yxz(yr,MALE);
    PopDyInfo* pPIF = new PopDyInfo(nZBs);//females info
    pPIF->R_z   = R_yz(yr);
    pPIF->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(FEMALE);
    pPIF->M_msz = M_yxmsz(yr,FEMALE);
    pPIF->T_szz = prGr_yxszz(yr,FEMALE);
    for (int s=1;s<=nSCs;s++) pPIF->Th_sz(s) = prM2M_yxz(yr,FEMALE);
    if (debug) cout<<"calculated pPIM, pPIF."<<endl;
    //5. Determine fishery conditions for next year based on averages for recent years
        int oflAvgPeriodYrs = 5;  //TODO: this should be an input
        //assumption here is that ALL fisheries EXCEPT the first are bycatch fisheries
        //a. Calculate average handling mortality, retention curves and capture rates
        int ny;   //number of years fishery is active
        dvar_vector avgHM_f(1,nFsh);
        avgHM_f.initialize();
        for (int f=1;f<=nFsh;f++){
            ny = 0;
            for (int y=yr-oflAvgPeriodYrs+1;y<=yr;y++){
                ny         += hasF_fy(f,y);
                avgHM_f(f) += hmF_fy(f,y);
            }
            avgHM_f(f) /= wts::max(1.0,1.0*ny);
        }
        if (debug) cout<<"avgHm_f = "<<avgHM_f<<endl;
        dvar4_array avgCapF_xfms(1,nSXs,1,nFsh,1,nMSs,1,nSCs);//averaged max fishery capture rates
        dvar5_array avgCapF_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery capture rates
        dvar5_array avgRFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery retention functions
        dvar5_array avgSFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery selectivity functions
        avgCapF_xfms.initialize();
        avgCapF_xfmsz.initialize();
        avgRFcn_xfmsz.initialize();
        avgSFcn_xfmsz.initialize();
        for (int f=1;f<=nFsh;f++){
            oflAvgPeriodYrs = ptrMOs->oflNumYrsForAvgCapRate(f);
            if (debug) cout<<"oflAvgPeriodYrs("<<f<<") = "<<oflAvgPeriodYrs<<endl;
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) {
                        for (int z=1;z<=nZBs;z++){
                            ny = 0;
                            for (int y=(yr-oflAvgPeriodYrs+1);y<=yr;y++) {
                                //if (debug) cout<<"y = "<<y<<endl;
                                ny += hasF_fy(f,y);
                                avgCapF_xfms(x,f,m,s) += cpF_fyxms(f,y,x,m,s);
                                avgCapF_xfmsz(x,f,m,s,z) += cpF_fyxmsz(f,y,x,m,s,z);
                                avgRFcn_xfmsz(x,f,m,s,z) += ret_fyxmsz(f,y,x,m,s,z);
                                avgSFcn_xfmsz(x,f,m,s,z) += sel_fyxmsz(f,y,x,m,s,z);
                            }//y
                            double fac = 1.0/max(1.0,1.0*ny);
                            avgCapF_xfms(x,f,m,s) *= fac;
                            avgCapF_xfmsz(x,f,m,s,z) *= fac;
                            avgRFcn_xfmsz(x,f,m,s,z) *= fac;
                            avgSFcn_xfmsz(x,f,m,s,z) *= fac;
                        }//z
                    }//s
                }//m
            }//x
        }//f
        if (debug){
            for (int f=1;f<=nFsh;f++){
                cout<<"avgCapF_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgCapF_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgRFcn_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgRFcn_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgSFcn_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgSFcn_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgSFcn_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgSFcn_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
            }
        }
        //determine average capture rates
        for (int f=1;f<=nFsh;f++){
            if (ptrMOs->optOFLAvgCapRate(f)==0){
                //use averaged selectivity functions, max capture rates
                if (debug) cout<<"Using average selectivity functions for fishery "<<f<<" for OFL calculations"<<endl;
                for (int x=1;x<=nSXs;x++){
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++) avgCapF_xfmsz(x,f,m,s) = avgCapF_xfms(x,f,m,s)*avgSFcn_xfmsz(x,f,m,s);
                    }//m
                }//x
            } else {
                if (debug) cout<<"Using size-specific average capture rates for fishery "<<f<<" for OFL calculations"<<endl;
                //use size-specific averaged capture rates
                //do nothing
            }
        }//f
        CatchInfo* pCIM = new CatchInfo(nZBs,nFsh);//male catch info
        pCIM->setCaptureRates(avgCapF_xfmsz(MALE));
        pCIM->setCaptureRates(avgCapF_xfms(MALE));
        pCIM->setSelectivityFcns(avgSFcn_xfmsz(MALE));
        pCIM->setRetentionFcns(avgRFcn_xfmsz(MALE));
        pCIM->setHandlingMortality(avgHM_f);
        dvariable maxCapF = pCIM->findMaxTargetCaptureRate(cout);
        if (debug) cout<<"maxCapF = "<<maxCapF<<endl;
        CatchInfo* pCIF = new CatchInfo(nZBs,nFsh);//female catch info
        pCIF->setCaptureRates(avgCapF_xfmsz(FEMALE));
        pCIF->setCaptureRates(avgCapF_xfms(FEMALE));
        pCIF->setSelectivityFcns(avgSFcn_xfmsz(FEMALE));
        pCIF->setRetentionFcns(avgRFcn_xfmsz(FEMALE));
        pCIF->setHandlingMortality(avgHM_f);
        pCIF->maxF = maxCapF;//need to set this for females
    //6. Create PopProjectors
        PopProjector* pPPM = new PopProjector(pPIM,pCIM);
        pPPM->dtF = dtF;
        pPPM->dtM = dtM;
        PopProjector* pPPF = new PopProjector(pPIF,pCIF);
        pPPF->dtF = dtF;
        pPPF->dtM = dtM;
        if (debug) cout<<"created pPPs."<<endl;
    //7. Create Equilibrium_Calculators
        Equilibrium_Calculator* pECM = new Equilibrium_Calculator(pPPM);
        Equilibrium_Calculator* pECF = new Equilibrium_Calculator(pPPF);
        if (debug) cout<<"created pECs."<<endl;
    //8. Define OFL_Calculator   
        OFL_Calculator*  pOC;
        if (debug) cout<<"declared pOC."<<endl;
    //9. Determine TIER LEVEL, define Tier_Calculators, calculate OFL
        int tier = 3;
        if (tier==3){
            //5. Determine Fmsy and Bmsy
            Tier3_Calculator* pT3CM = new Tier3_Calculator(0.35,pECM);
            Tier3_Calculator* pT3CF = new Tier3_Calculator(0.35,pECF);
            if (debug) cout<<"created pT3Cs."<<endl;
            pOC = new OFL_Calculator(pT3CM,pT3CF);
            if (debug) {
                cout<<"created pOC."<<endl;
                OFL_Calculator::debug=1;
                Tier3_Calculator::debug=1;
                Equilibrium_Calculator::debug=0;
                cout<<"Calculating ptrOFLResults"<<endl;
            }
            ptrOFLResults = pOC->calcOFLResults(avgRec_x,n_xmsz,cout);
            if (debug) {
                cout<<"calculated ptrOFLResults->"<<endl;
                ptrOFLResults->writeCSVHeader(cout); cout<<endl;
                ptrOFLResults->writeToCSV(cout); cout<<endl;
                ptrOFLResults->writeToR(cout,ptrMC,"oflResults",0); cout<<endl;
                OFL_Calculator::debug=0;
                Tier3_Calculator::debug=0;
                Equilibrium_Calculator::debug=0;
            }
        }//Tier 3 calculation
    if (debug) {
        int n = 100;
        MultiYearPopProjector* pMYPPM = new MultiYearPopProjector(pPPM);
        MultiYearPopProjector* pMYPPF = new MultiYearPopProjector(pPPF);
        pMYPPM->projectUnFished(n,avgRec_x(  MALE),n_xmsz(  MALE),cout);
        dvariable myPP_B0 = pMYPPM->matBio_y(n);
        pMYPPM->project(n,avgRec_x(  MALE),ptrOFLResults->Fmsy,n_xmsz(  MALE),cout);
        pMYPPF->project(n,avgRec_x(FEMALE),ptrOFLResults->Fmsy,n_xmsz(FEMALE),cout);
        dvariable myPP_Bmsy = pMYPPM->matBio_y(n);
        dvariable myPP_prjB = pMYPPM->matBio_y(1);
        dvariable myPP_MSY  = pMYPPM->totCM_y(n)+pMYPPF->totCM_y(n);
        dvariable myPP_OFL  = pMYPPM->totCM_y(1)+pMYPPF->totCM_y(1);
        cout<<"#------------------------"<<endl;
        cout<<"MYPP B0     = "<<myPP_B0  <<". B0     = "<<ptrOFLResults->B0  <<endl;
        cout<<"MYPP Bmsy   = "<<myPP_Bmsy<<". Bmsy   = "<<ptrOFLResults->Bmsy<<endl;
        cout<<"MYPP prjMMB = "<<myPP_prjB<<". prjMMB = "<<ptrOFLResults->prjB<<endl;
        cout<<"MYPP MSY    = "<<myPP_MSY <<". MSY    = "<<ptrOFLResults->MSY <<endl;
        cout<<"MYPP OFL    = "<<myPP_OFL <<". OFL    = "<<ptrOFLResults->OFL <<endl;
        cout<<"#------------------------"<<endl;
        cout<<"finished calcOFL(yr,debug,cout)"<<endl<<endl<<endl;
    }
}

void model_parameters::calcOFL_OpMod(int debug, ostream& cout)
{
    if (debug) {
        cout<<endl<<endl<<"#------------------------"<<endl;
        cout<<"starting calcOFL_OpMod(debug,cout)"<<endl;
        //cout<<"year for projection = "<<yr<<endl;
    }
    //1. get initial population -- unlike calcOFL there is no year component
    //dvar4_array n_xmsz = n_yxmsz(yr); Get rid of this line, as we don't want years
    if (debug) {cout<<"  males_msz:"<<endl; wts::print(prj_n_xmsz(  MALE),cout,1);}
    if (debug) {cout<<"females_msz:"<<endl; wts::print(prj_n_xmsz(FEMALE),cout,1);}
    //2. set yr back one year to get population rates, etc.,
        // -----CHECK THIS, the Op model is the projection, so yr=yr 
    //yr = yr;//don't have pop rates, etc. for projection year --- KEEP YEAR AT YEAR?
    //if (debug) cout<<"year for pop rates = "<<yr<<endl;
    //3. Determine mean recruitment --Keep this?
        //NOT SURE IF THIS SHOULD BE ALTERED 
    //   1981 here corresponds to 1982 in TCSAM2013, the year recruitment enters
    //   the model population.
    dvar_vector avgRec_x(1,nSXs); 
    if (debug) cout<<"R dims: "<<R_y.indexmin()<<cc<<R_y.indexmax()<<endl;
    int yr = R_y.indexmax();
    for (int x=1;x<=nSXs;x++) 
        avgRec_x(x)= mean(ptrOMI->R_y(1981,yr))*ptrOMI->R_x(x);
    if (debug) {
        cout<<"Average recruitment = "<<avgRec_x<<endl;
    }
    //4. Identify population rates
    double dtF = ptrOMI->dtF;//time at which fisheries occur
    double dtM = ptrOMI->dtM;//time at which mating occurs
    PopDyInfo* pPIM = new PopDyInfo(nZBs);//  males info
    pPIM->R_z = ptrOMI->R_z; //Recruitment by size
    //pPIM->R_x = ptrOMI->R_x; // Recruitment by sex 
    pPIM->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(MALE); // changed to just sex
    pPIM->M_msz = ptrOMI->M_xmsz(MALE);           // changed to just sex, added pointer
    pPIM->T_szz = ptrOMI->prGr_xszz(MALE);        // changed to just sex, added pointer
    for (int s=1;s<=nSCs;s++) pPIM->Th_sz(s) = ptrOMI->prM2M_xz(MALE); // changed to just sex, added pointer
    PopDyInfo* pPIF = new PopDyInfo(nZBs);//females info
    pPIF->R_z = ptrOMI->R_z;  // = R_yz(yr);
    //pPIM->R_x = ptrOMI->R_x; // Recruitment by sex 
    pPIF->w_mz = ptrMDS->ptrBio->wAtZ_xmz(FEMALE);  // changed to just sex 
    pPIF->M_msz = ptrOMI->M_xmsz(FEMALE);           // changed to just sex, added pointer
    pPIF->T_szz = ptrOMI->prGr_xszz(FEMALE);        // changed to just sex, added pointer
    for (int s=1;s<=nSCs;s++) pPIF->Th_sz(s) = ptrOMI->prM2M_xz(FEMALE); // changed to just sex, added pointer
    if (debug) cout<<"calculated pPIM, pPIF."<<endl;
    //5. Determine fishery conditions for next year based on averages for recent years
      // NOT SURE WHAT TO DO HERE --MS
        int oflAvgPeriodYrs = 1;  //TODO: this should be an input
        //assumption here is that ALL fisheries EXCEPT the first are bycatch fisheries
        //a. Calculate average handling mortality, retention curves and capture rates
        int ny;   //number of years fishery is active
        dvar_vector avgHM_f(1,nFsh);
        avgHM_f=ptrOMI->hmF_f;
        dvar5_array avgCapF_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery capture rates
        dvar5_array avgRFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged retention function
        // Added for HCR7
         dvar5_array avgSFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery selectivity functions
        avgCapF_xfmsz.initialize();
        avgRFcn_xfmsz.initialize();
        avgSFcn_xfmsz.initialize(); // Added for HCR 7 
        for (int f=1;f<=nFsh;f++){ 
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++){ 
                        avgCapF_xfmsz(x,f,m,s) = ptrOMI->cpF_fxmsz(f,x,m,s); 
                        avgRFcn_xfmsz(x,f,m,s) = ptrOMI->ret_fxmsz(f,x,m,s);
                        avgSFcn_xfmsz(x,f,m,s) = ptrOMI->sel_fxmsz(f,x,m,s); // Added for HCR7
                    }//s
                }//m
            }//x
        }//f
        CatchInfo* pCIM = new CatchInfo(nZBs,nFsh);//male catch info
        pCIM->setCaptureRates(avgCapF_xfmsz(MALE));
        //pCIM->setCaptureRates(avgCapF_xfms(MALE));
        pCIM->setSelectivityFcns(avgSFcn_xfmsz(MALE));
        pCIM->setRetentionFcns(avgRFcn_xfmsz(MALE));
        pCIM->setHandlingMortality(avgHM_f);
        dvariable maxCapF = pCIM->findMaxTargetCaptureRate(cout);
        if (debug) cout<<"maxCapF = "<<maxCapF<<endl;
        CatchInfo* pCIF = new CatchInfo(nZBs,nFsh);//female catch info
        pCIF->setCaptureRates(avgCapF_xfmsz(FEMALE));
        //pCIF->setCaptureRates(avgCapF_xfms(FEMALE));
        pCIF->setSelectivityFcns(avgSFcn_xfmsz(FEMALE));
        pCIF->setRetentionFcns(avgRFcn_xfmsz(FEMALE));
        pCIF->setHandlingMortality(avgHM_f);
        pCIF->maxF = maxCapF;//need to set this for females
    //6. Create PopProjectors
        PopProjector* pPPM = new PopProjector(pPIM,pCIM);
        pPPM->dtF = dtF;
        pPPM->dtM = dtM;
        PopProjector* pPPF = new PopProjector(pPIF,pCIF);
        pPPF->dtF = dtF;
        pPPF->dtM = dtM;
        if (debug) cout<<"created pPPs."<<endl;
    //7. Create Equilibrium_Calculators
        Equilibrium_Calculator* pECM = new Equilibrium_Calculator(pPPM);
        Equilibrium_Calculator* pECF = new Equilibrium_Calculator(pPPF);
        if (debug) cout<<"created pECs."<<endl;
    //8. Define OFL_Calculator   
        OFL_Calculator*  pOC;
        if (debug) cout<<"declared pOC."<<endl;
    //9. Determine TIER LEVEL, define Tier_Calculators, calculate OFL
        int tier = 3;
        if (tier==3){
            //5. Determine Fmsy and Bmsy
            Tier3_Calculator* pT3CM = new Tier3_Calculator(0.35,pECM);
            Tier3_Calculator* pT3CF = new Tier3_Calculator(0.35,pECF);
            if (debug) cout<<"created pT3Cs."<<endl;
            pOC = new OFL_Calculator(pT3CM,pT3CF);
            if (debug) {
                cout<<"created pOC."<<endl;
                OFL_Calculator::debug=1;
                Tier3_Calculator::debug=1;
                Equilibrium_Calculator::debug=0;
                cout<<"Calculating ptrOFLResults"<<endl;
            }
            ptrOFLResults = pOC->calcOFLResults(avgRec_x,prj_n_xmsz,cout);
            if (debug) {
                cout<<"calculated ptrOFLResults->"<<endl;
                ptrOFLResults->writeCSVHeader(cout); cout<<endl;
                ptrOFLResults->writeToCSV(cout); cout<<endl;
                ptrOFLResults->writeToR(cout,ptrMC,"oflResults",0); cout<<endl;
                OFL_Calculator::debug=0;
                Tier3_Calculator::debug=0;
                Equilibrium_Calculator::debug=0;
            }
        }//Tier 3 calculation
        PRINT2B1("OFL Op model function done")
}

void model_parameters::calcPenalties(int debug, ostream& cout)
{
    if (debug>dbgObjFun) cout<<"Started calcPenalties()"<<endl;
    if (debug<0) cout<<"list("<<endl;//start list of penalties by category
    //growth-related penalties
    if (debug<0) cout<<tb<<"growth=list(negativeGrowth=list("<<endl;//start of growth penalties list
    GrowthInfo* ptrGrw = ptrMPI->ptrGrw;
    //calculate growth parameters on arithmetic scale
    dvar_vector ptGrA = ptrGrw->pGrA->calcArithScaleVals(pGrA);
    dvar_vector ptGrB = ptrGrw->pGrB->calcArithScaleVals(pGrB);
    //loop over parameter combinations
    dvariable grA, grB;
    dvector zBsp(1,nZBs+2);
    dvar_vector dZ(1,nZBs+2);
    zBsp(1,nZBs) = zBs;
    zBsp(nZBs+1) = ptrMOs->minGrowthCW;
    zBsp(nZBs+2) = ptrMOs->maxGrowthCW;
    for (int pc=1;pc<=(ptrGrw->nPCs-1);pc++){
        ivector pids = ptrGrw->getPCIDs(pc);
        int k=ptrGrw->nIVs+1;//1st parameter column
        grA = ptGrA(pids[k]); k++; //"a" coefficient for mean growth
        grB = ptGrB(pids[k]); k++; //"b" coefficient for mean growth
        //add objective function penalties to keep mean size increments positive
        dvariable pen; pen.initialize();
        dZ.initialize();
        if (ptrMOs->optGrowthParam==0){
            dZ = mfexp(grA+grB*log(zBsp)) - zBsp;
        } else if (ptrMOs->optGrowthParam==1){
            dvector pXDs = ptrGrw->getPCXDs(pc);
            double zGrA = pXDs[1];
            double zGrB = pXDs[2];
            dZ = grA*mfexp(log(grB/grA)/log(zGrB/zGrA)*log(zBsp/zGrA)) - zBsp;
        } else if (ptrMOs->optGrowthParam==2){
            dvector pXDs = ptrGrw->getPCXDs(pc);
            double zGrA = pXDs[1];
            double zGrB = pXDs[2];
            dZ = grA*mfexp(grB*log(zBsp/zGrA)) - zBsp;
        }
        posfun(dZ,ptrMOs->epsNegGrowth,pen);
        if (pen>0.0){
            rpt::echo<<"--Growth Increments Approaching 0 for pc = "<<pc<<". pen = "<<pen<<endl;
            rpt::echo<<"params = "<<grA<<cc<<grB<<endl;
            rpt::echo<<"zBsp = "<<zBsp<<endl;
            rpt::echo<<"dZ   = "<<dZ<<endl;
            rpt::echo<<"--------"<<endl;
        }
        objFun += ptrMOs->wgtNegGrowth*pen;
        if (debug<0) {
            cout<<tb<<tb<<tb<<"'"<<pc<<"'=list(wgt="<<ptrMOs->wgtNegGrowth<<cc<<"pen="<<pen<<cc<<"objfun="<<ptrMOs->wgtNegGrowth*pen<<"),"<<endl;
        }
    }
    {
        int pc = ptrGrw->nPCs;
        ivector pids = ptrGrw->getPCIDs(pc);
        int k=ptrGrw->nIVs+1;//1st parameter column
        grA = (*ptrMPI->ptrGrw->pGrA)[pids[k]]->calcArithScaleVal(pGrA(pids[k])); k++; //"a" coefficient for mean growth
        grB = (*ptrMPI->ptrGrw->pGrB)[pids[k]]->calcArithScaleVal(pGrB(pids[k])); k++; //"b" coefficient for mean growth
        //add objective function penalty to keep mean size increments positive
        dvariable pen; pen.initialize();
        dZ.initialize();
        if (ptrMOs->optGrowthParam==0){
            dZ = mfexp(grA+grB*log(zBsp)) - zBsp;
        } else if (ptrMOs->optGrowthParam==1){
            dvector pXDs = ptrGrw->getPCXDs(pc);
            double zGrA = pXDs[1];
            double zGrB = pXDs[2];
            dZ = grA*mfexp(log(grB/grA)/log(zGrB/zGrA)*log(zBsp/zGrA)) - zBsp;
        } else if (ptrMOs->optGrowthParam==2){
            dvector pXDs = ptrGrw->getPCXDs(pc);
            double zGrA = pXDs[1];
            double zGrB = pXDs[2];
            dZ = grA*mfexp(grB*log(zBsp/zGrA)) - zBsp;
        }
        posfun(dZ,ptrMOs->epsNegGrowth,pen);
        if (pen>0.0){
            rpt::echo<<"--Growth Increments Approaching 0 for pc = "<<pc<<". pen = "<<pen<<endl;
            rpt::echo<<"params = "<<grA<<cc<<grB<<endl;
            rpt::echo<<"zBsp = "<<zBsp<<endl;
            rpt::echo<<"dZ   = "<<dZ<<endl;
            rpt::echo<<"--------"<<endl;
        }
        objFun += ptrMOs->wgtNegGrowth*pen;
        if (debug<0) {
            cout<<tb<<tb<<tb<<"'"<<pc<<"'=list(wgt="<<ptrMOs->wgtNegGrowth<<cc<<"pen="<<pen<<cc<<"objfun="<<ptrMOs->wgtNegGrowth*pen<<")"<<endl;
        }
    }
    if (debug<0) cout<<tb<<tb<<"))"<<cc<<endl;//end of growth penalties list
    //maturity-related penalties
    if (debug<0) cout<<tb<<"maturity=list("<<endl;//start of maturity penalties list
    //smoothness penalties
    dvector penWgtSmthLgtPrMat = ptrMOs->wgtPenSmthPrM2M;
    if (ptrMOs->optPenSmthPrM2M==0){
        //smoothness penalties on maturity PARAMETERS (NOT maturity ogives)
        fPenSmoothLgtPrMat.initialize();
        if (debug<0) cout<<tb<<tb<<"smoothness=list(";//start of smoothness penalties list
        for (int i=1;i<npLgtPrMat;i++){
            dvar_vector v = 1.0*pvLgtPrM2M(i);
            fPenSmoothLgtPrMat(i) = 0.5*norm2(calc2ndDiffs(v));
            objFun += penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i);
            if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat(i)<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i)<<"),"<<endl;
        }
        {
            int i = npLgtPrMat;
            dvar_vector v = 1.0*pvLgtPrM2M(i);
            fPenSmoothLgtPrMat(i) = 0.5*norm2(calc2ndDiffs(v));
            objFun += penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i);
            if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat(i)<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i)<<")"<<endl;
        }
        if (debug<0) cout<<tb<<tb<<")"<<cc<<endl;//end of smoothness penalties list
    } else if (ptrMOs->optPenSmthPrM2M==1){
        //smoothness penalties on maturity OGIVES (NOT maturity parameters)
        fPenSmoothLgtPrMat.initialize();
        if (debug<0) cout<<tb<<tb<<"smoothness=list(";//start of smoothness penalties list
        for (int i=1;i<ptrMPI->ptrM2M->nPCs;i++){
            dvar_vector v = 1.0*prM2M_cz(i);
            fPenSmoothLgtPrMat(i) = 0.5*norm2(calc2ndDiffs(v));
            objFun += penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i);
            if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat(i)<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i)<<"),"<<endl;
        }
        {
            int i = ptrMPI->ptrM2M->nPCs;
            dvar_vector v = 1.0*prM2M_cz(i);
            fPenSmoothLgtPrMat(i) = 0.5*norm2(calc2ndDiffs(v));
            objFun += penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i);
            if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat(i)<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i)<<")"<<endl;
        }
        if (debug<0) cout<<tb<<tb<<")"<<cc<<endl;//end of smoothness penalties list
    }
    //penalties for decreasing maturity parameters/ogives
    dvector penWgtNonDecLgtPrMat = ptrMOs->wgtPenNonDecPrM2M;
    fPenNonDecLgtPrMat.initialize();
    if (debug<0) cout<<tb<<tb<<"nondecreasing=list(";//start of non-decreasing penalties list
    int np;
    if (ptrMOs->optPenNonDecPrM2M==0||ptrMOs->optPenNonDecPrM2M==1) np = npLgtPrMat;
    else if (ptrMOs->optPenNonDecPrM2M==2||ptrMOs->optPenNonDecPrM2M==3) np = ptrMPI->ptrM2M->nPCs;
    for (int i=1;i<np;i++){
        dvar_vector v; 
        if (ptrMOs->optPenNonDecPrM2M==0){
            v = calc1stDiffs(pvLgtPrM2M(i));
            for (int iv=v.indexmin();iv<=v.indexmax();iv++){
                posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
            } 
        } else if (ptrMOs->optPenNonDecPrM2M==1){
            v = calc1stDiffs(pvLgtPrM2M(i));
            fPenNonDecLgtPrMat(i) = sum(mfexp(-10.0*v));
        } else if (ptrMOs->optPenNonDecPrM2M==2){
            v = calc1stDiffs(prM2M_cz(i));
            for (int iv=v.indexmin();iv<=v.indexmax();iv++){
                posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
            } 
        } else if (ptrMOs->optPenNonDecPrM2M==3){
            v = calc1stDiffs(prM2M_cz(i));
            fPenNonDecLgtPrMat(i) = sum(mfexp(-10.0*v));
        }
        objFun += penWgtNonDecLgtPrMat(i)*fPenNonDecLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtNonDecLgtPrMat(i)<<cc<<"pen="<<fPenNonDecLgtPrMat(i)<<cc<<"objfun="<<penWgtNonDecLgtPrMat(i)*fPenNonDecLgtPrMat(i)<<"),";
    }
    {
        int i = np;
        dvar_vector v; 
        if (ptrMOs->optPenNonDecPrM2M==0){
            v = calc1stDiffs(pvLgtPrM2M(i));
            for (int iv=v.indexmin();iv<=v.indexmax();iv++){
                posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
            }
        } else if (ptrMOs->optPenNonDecPrM2M==1){
            v = calc1stDiffs(pvLgtPrM2M(i));
            fPenNonDecLgtPrMat(i) = sum(mfexp(-10.0*v));
        } else if (ptrMOs->optPenNonDecPrM2M==2){
            v = calc1stDiffs(prM2M_cz(i));
            for (int iv=v.indexmin();iv<=v.indexmax();iv++){
                posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
            } 
        } else if (ptrMOs->optPenNonDecPrM2M==3){
            v = calc1stDiffs(prM2M_cz(i));
            fPenNonDecLgtPrMat(i) = sum(mfexp(-10.0*v));
        }
        objFun += penWgtNonDecLgtPrMat(i)*fPenNonDecLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtNonDecLgtPrMat(i)<<cc<<"pen="<<fPenNonDecLgtPrMat(i)<<cc<<"objfun="<<penWgtNonDecLgtPrMat(i)*fPenNonDecLgtPrMat(i)<<")";
    }
    if (debug<0) cout<<tb<<tb<<")"<<endl;//end of non-decreasing penalties list    
    if (debug<0) cout<<tb<<"),";//end of maturity penalties list
    //penalties on nonparametric selectivity functions
    if (debug<0) cout<<tb<<"nonParSelFcns=list("<<endl;//start of nonparametric selectivity function penalties list
    //smoothness penalties
    fPenSmoothNPSel.initialize();
    if (ptrMPI->ptrSel->pvNPSel->getSize()>0){
        dvector penWgtSmthNPSel = ptrMOs->wgtPenSmthNPSel;
        if (ptrMOs->optPenSmthNPSel==0){
            //smoothness penalties on selectivity PARAMETERS (NOT selectivity functions)
            if (debug<0) cout<<tb<<tb<<"smoothness=list(";//start of smoothness penalties list
            for (int i=1;i<npNPSel;i++){
                dvar_vector v = 1.0*pvNPSel(i);
                fPenSmoothNPSel(i) = 0.5*norm2(calc2ndDiffs(v));
                objFun += penWgtSmthNPSel(i)*fPenSmoothNPSel(i);
                if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthNPSel(i)<<cc<<"pen="<<fPenSmoothNPSel(i)<<cc<<"objfun="<<penWgtSmthNPSel(i)*fPenSmoothNPSel(i)<<"),"<<endl;
            }
            {
                int i = npNPSel;
                dvar_vector v = 1.0*pvNPSel(i);
                fPenSmoothNPSel(i) = 0.5*norm2(calc2ndDiffs(v));
                objFun += penWgtSmthNPSel(i)*fPenSmoothNPSel(i);
                if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthNPSel(i)<<cc<<"pen="<<fPenSmoothNPSel(i)<<cc<<"objfun="<<penWgtSmthNPSel(i)*fPenSmoothNPSel(i)<<")"<<endl;
            }
            if (debug<0) cout<<tb<<tb<<")"<<cc<<endl;//end of smoothness penalties list
        } else if (ptrMOs->optPenSmthNPSel==1){
            //smoothness penalties on nonparametric selectivity curves (NOT parameter vectors)
            fPenSmoothNPSel.initialize();
            if (debug<0) cout<<tb<<tb<<"smoothness=list(";//start of smoothness penalties list
            for (int i=ptrMOs->wgtPenSmthNPSel.indexmin();i<ptrMOs->wgtPenSmthNPSel.indexmax();i++){
                dvar_vector v = 1.0*npSel_cz(i);
                fPenSmoothNPSel(i) = 0.5*norm2(calc2ndDiffs(v));
                objFun += penWgtSmthNPSel(i)*fPenSmoothNPSel(i);
                if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthNPSel(i)<<cc<<"pen="<<fPenSmoothNPSel(i)<<cc<<"objfun="<<penWgtSmthNPSel(i)*fPenSmoothNPSel(i)<<"),"<<endl;
            }
            {
                int i = ptrMOs->wgtPenSmthNPSel.indexmax();
                dvar_vector v = 1.0*npSel_cz(i);
                fPenSmoothNPSel(i) = 0.5*norm2(calc2ndDiffs(v));
                objFun += penWgtSmthNPSel(i)*fPenSmoothNPSel(i);
                if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthNPSel(i)<<cc<<"pen="<<fPenSmoothNPSel(i)<<cc<<"objfun="<<penWgtSmthNPSel(i)*fPenSmoothNPSel(i)<<")"<<endl;
            }
        }
        if (debug<0) cout<<tb<<tb<<")"<<cc<<endl;//end of smoothness penalties list
    }
    //penalties on sums of dev vectors to enforce sum-to-zero
    double penWgt = 0.0;
    if (current_phase()>=ptrMOs->phsSqSumDevsPen) penWgt = ptrMOs->wgtSqSumDevsPen;
    if (debug<0) rpt::echo<<"#----Check on sum(devs) in calcPenalties "<<current_phase()<<tb<<ctrProcCallsInPhase<<endl;
    if (debug<0) cout<<tb<<"devsSumSq=list("<<endl;//start of devs penalties list
    //recruitment devs
    if (ptrMPI->ptrRec->pDevsLnR->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsLnR = ";
            for (int i=1;i<=pDevsLnR.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsLnR[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsLnR=";
        calcDevsPenalties(debug,cout,penWgt,pDevsLnR,devsLnR);        
        if (debug<0) cout<<cc<<endl;
    }
    //S1 devs
    if (ptrMPI->ptrSel->pDevsS1->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS1 = ";
            for (int i=1;i<=pDevsS1.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS1[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS1=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS1,devsS1);
        if (debug<0) cout<<cc<<endl;
    }
    //S2 devs
    if (ptrMPI->ptrSel->pDevsS2->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS2 = ";
            for (int i=1;i<=pDevsS2.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS2[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS2=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS2,devsS2);
        if (debug<0) cout<<cc<<endl;
    }
    //S3 devs
    if (ptrMPI->ptrSel->pDevsS3->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS3 = ";
            for (int i=1;i<=pDevsS3.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS3[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS3=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS3,devsS3);
        if (debug<0) cout<<cc<<endl;
    }
    //S4 devs
    if (ptrMPI->ptrSel->pDevsS4->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS4 = ";
            for (int i=1;i<=pDevsS4.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS4[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS4=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS4,devsS4);
        if (debug<0) cout<<cc<<endl;
    }
    //S5 devs
    if (ptrMPI->ptrSel->pDevsS5->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS5 = ";
            for (int i=1;i<=pDevsS5.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS5[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS5=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS5,devsS5);
        if (debug<0) cout<<cc<<endl;
    }
    //S6 devs
    if (ptrMPI->ptrSel->pDevsS6->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS6 = ";
            for (int i=1;i<=pDevsS6.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS6[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS6=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS6,devsS6);
        if (debug<0) cout<<cc<<endl;
    }
    //capture rate devs
    if (ptrMPI->ptrFsh->pDevsLnC->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsLnC = ";
            for (int i=1;i<=pDevsLnC.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsLnC[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsLnC=";
        calcDevsPenalties(debug,cout,penWgt,pDevsLnC,devsLnC);
        if (debug<0) cout<<cc<<endl;
    }
    if (debug<0) rpt::echo<<"#------"<<endl<<endl;
    if (debug<0) cout<<tb<<"NULL)"<<endl;//end of devs penalties list
    if (debug<0) cout<<")"<<cc<<endl;//end of penalties list
    //Apply diminishing penalties to variances of F-devs
    {
        if (debug<0) cout<<"penFDevs=list("<<endl;
        double penWgt = 1.0/log(1.0+square(ptrMOs->cvFDevsPen));
        if (ptrMOs->phsDecrFDevsPen<=current_phase()){
            double scl = max(1.0,(double)(ptrMOs->phsZeroFDevsPen-ptrMOs->phsDecrFDevsPen));
            penWgt *= max(0.0,(double)(ptrMOs->phsZeroFDevsPen-current_phase()))/scl;
        }
        double effCV = std::numeric_limits<double>::infinity();
        if (penWgt>0) {
            effCV = sqrt(mfexp(1.0/penWgt)-1.0);
            if (debug<0) rpt::echo<<"phase: "<<current_phase()<<"; penWgt = "<<penWgt<<"; effCV = "<<effCV<<endl;
        } else {
            if (debug<0) rpt::echo<<"phase: "<<current_phase()<<"; penWgt = "<<penWgt<<"; effCV = Inf"<<endl;
        }
        for (int i=pDevsLnC.indexmin();i<=pDevsLnC.indexmax();i++){
            dvariable fpen = 0.5*norm2(pDevsLnC(i));
            objFun += penWgt*fpen;
            if (debug<0) {
                double rmse = sqrt(value(norm2(pDevsLnC(i)))/pDevsLnC(i).size());
                rpt::echo<<tb<<i<<": pen="<<fpen<<cc<<" objfun="<<penWgt*fpen<<endl;
                if (penWgt>0) {
                    cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"effCV="<<effCV<<cc
                            <<"pen="<<fpen<<cc<<"rmse="<<rmse<<cc<<"objfun="<<penWgt*fpen<<"),"<<endl;
                } else {
                    cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"effCV=Inf"    <<cc
                            <<"pen="<<fpen<<cc<<"rmse="<<rmse<<cc<<"objfun="<<penWgt*fpen<<"),"<<endl;
                }
            }
        }
        if (debug<0) cout<<tb<<tb<<"NULL)"<<endl;//end of penFDevs lists
    }
    if (!debug) testNaNs(value(objFun),"in calcPenalties()");
    if (debug>dbgObjFun) cout<<"Finished calcPenalties()"<<endl;
}

void model_parameters::calcDevsPenalties(int debug, ostream& cout, double penWgt, param_init_bounded_vector_vector& pDevs, dvar_matrix devs)
{
    dvariable fPen;
    if (debug<0) cout<<"list(";//start of list
    for (int i=pDevs.indexmin();i<=pDevs.indexmax();i++){
        if (pDevs(i).get_phase_start()){
            fPen.initialize();
            fPen = square(sum(devs(i)));
            objFun += penWgt*fPen;
            if (debug<0) {
                double rmse = sqrt(value(norm2(devs(i)))/devs(i).size());
                cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"pen="<<fPen<<cc<<"objfun="<<penWgt*fPen<<cc<<"val="<<sum(devs(i))<<cc<<
                                                  "rmse="<<rmse<<"),"<<endl;
            }
        }
    }
    if (!debug) testNaNs(value(objFun),"in calcDevsPenalties()");
    if (debug<0) cout<<tb<<tb<<"NULL)";//end of penalties list
}

dvar_vector model_parameters::calc1stDiffs(const dvar_vector& d)
{
    RETURN_ARRAYS_INCREMENT();
    dvar_vector r = first_difference(d);
    RETURN_ARRAYS_DECREMENT();
    return r;
}

dvar_vector model_parameters::calc2ndDiffs(const dvar_vector& d)
{
    RETURN_ARRAYS_INCREMENT();
    dvar_vector r = first_difference(first_difference(d));
    RETURN_ARRAYS_DECREMENT();
    return r;
}

void model_parameters::calcNLLs_Recruitment(int debug, ostream& cout)
{
    if (debug>dbgObjFun) cout<<"Starting calcNLLs_Recruitment"<<endl;
    dvariable nllRecDevs;
    if (debug<0) cout<<"list("<<endl;
    if (debug<0) cout<<tb<<"recDevs=list("<<endl;
    for (int pc=1;pc<npcRec;pc++){
        double nllWgtRecDevs = ptrMPI->ptrRec->xd(pc,1);
        nllRecDevs.initialize();
        nllRecDevs = 0.5*norm2(zscrDevsLnR_cy(pc));
        nllRecDevs += nDevsLnR_c(pc)*log(stdvDevsLnR_c(pc));
        objFun += nllWgtRecDevs*nllRecDevs;
        if (debug<0){
            double rmse = sqrt(value(norm2(devsLnR_cy(pc)))/nDevsLnR_c(pc));
            cout<<tb<<tb<<"'"<<pc<<"'=list(type='normal',wgt="<<nllWgtRecDevs<<cc<<"nll="<<nllRecDevs<<cc
                    <<"objfun="<<nllWgtRecDevs*nllRecDevs<<cc
                    <<"n="<<nDevsLnR_c(pc)<<cc<<"rmse="<<rmse<<cc<<"sigmaR="<<stdvDevsLnR_c(pc)<<cc<<endl;
            cout<<"devs="; wts::writeToR(cout,value(devsLnR_cy(pc))); cout<<cc<<endl;
            cout<<"zscrs="; wts::writeToR(cout,value(zscrDevsLnR_cy(pc))); cout<<")"<<cc<<endl;
        }
    }//pc
    {
        int pc = npcRec;
        double nllWgtRecDevs = ptrMPI->ptrRec->xd(pc,1);
        nllRecDevs.initialize();
        nllRecDevs = 0.5*norm2(zscrDevsLnR_cy(pc));
        nllRecDevs += nDevsLnR_c(pc)*log(stdvDevsLnR_c(pc));
        objFun += nllWgtRecDevs*nllRecDevs;
        if (debug<0){
            double rmse = sqrt(value(norm2(devsLnR_cy(pc)))/nDevsLnR_c(pc));
            cout<<tb<<tb<<"'"<<pc<<"'=list(type='normal',wgt="<<nllWgtRecDevs<<cc<<"nll="<<nllRecDevs<<cc
                    <<"objfun="<<nllWgtRecDevs*nllRecDevs<<cc
                    <<"n="<<nDevsLnR_c(pc)<<cc<<"rmse="<<rmse<<cc<<"sigmaR="<<stdvDevsLnR_c(pc)<<cc<<endl;
            cout<<"devs="; wts::writeToR(cout,value(devsLnR_cy(pc))); cout<<cc<<endl;
            cout<<"zscrs="; wts::writeToR(cout,value(zscrDevsLnR_cy(pc))); cout<<")"<<endl;
        }
    }//pc
    if (debug<0) cout<<tb<<")";//recDevs
    if (debug<0) cout<<")";
    if (!debug) testNaNs(value(objFun),"in calcNLLs_Recruitment()");
    if (debug>dbgObjFun) cout<<"Finished calcNLLs_Recruitment"<<endl;
}

void model_parameters::calcObjFunForTAC(int debug, ostream& cout)
{
    if ((debug>=dbgObjFun)||(debug<0)) {PRINT2B1("----Starting calcObjFunForTAC()")}
    int k = 0;
    dvector objFunV(0,20);
    //reset objective function
    objFun.initialize();                     objFunV[k++] = value(objFun);
    //Fit OpMod catch to TAC
    dvariable prdTAC = 0.0;
    dvariable prdTot = 0.0;
    for (int f=1;f<=nFsh;f++){
        prdTAC += sum(prjRetCatchMortBio_fx(f));
        prdTot += sum(prjTotCatchMortBio_fx(f));
    }
    dvariable nllForTAC = square(inpTAC-prdTAC);
    dvariable nllForOFL = 0.0; 
    posfun(inpOFL-prdTot,1.0e-2,nllForOFL);
    nllForOFL *= 100.0;
    objFun += (nllForTAC + nllForOFL);   objFunV[k++] = value(objFun);
    if ((debug>=dbgObjFun)||(debug<0)){
        int k=0;
        PRINT2B2("proc call          = ",ctrProcCalls)
        PRINT2B2("proc call in phase = ",ctrProcCallsInPhase)
        PRINT2B2("pMSE_LnC =",pMSE_LnC[1]);
        PRINT2B2("mseCapF  =",mfexp(pMSE_LnC[1]));
        PRINT2B2("after initialization. objFun =",objFunV[k++])
        PRINT2B2("inpTAC    = ",inpTAC)    
        PRINT2B2("prdTAC    = ",prdTAC)    
        PRINT2B2("nllForTAC = ",nllForTAC)    
        PRINT2B2("inpOFL    = ",inpOFL)    
        PRINT2B2("prdTot    = ",prdTot)    
        PRINT2B2("nllForOFL = ",nllForOFL)    
        PRINT2B2("after MSE OpMod call objFun =",value(objFun))
        PRINT2B1("----Finished calcObjFunForTAC()")
    }
}

void model_parameters::calcObjFun(int debug, ostream& cout)
{
    if ((debug>=dbgObjFun)||(debug<0)) {PRINT2B1("----Starting calcObjFun()")}
    int k = 0;
    dvector objFunV(0,20);
    //reset objective function
    objFun.initialize();                     objFunV[k++] = value(objFun);
    //objective function penalties
    calcPenalties(debug,cout);               objFunV[k++] = value(objFun);
    //prior likelihoods
    calcAllPriors(debug,cout);               objFunV[k++] = value(objFun);
    //recruitment component
    calcNLLs_Recruitment(debug,cout);        objFunV[k++] = value(objFun);    
    //effort extrapolation
    calcNLLs_ExtrapolatedEffort(debug,cout); objFunV[k++] = value(objFun);
    //data components
    calcNLLs_Fisheries(debug,cout);       objFunV[k++] = value(objFun);
    calcNLLs_Surveys(debug,cout);         objFunV[k++] = value(objFun);
    calcNLLs_GrowthData(debug,cout);      objFunV[k++] = value(objFun);
    calcNLLs_ChelaHeightData(debug,cout); objFunV[k++] = value(objFun);
    if ((debug>=dbgObjFun)||(debug<0)){
        PRINT2B2("proc call          = ",ctrProcCalls)
        PRINT2B2("proc call in phase = ",ctrProcCallsInPhase)
        PRINT2B2("after initialization. objFun =",objFunV[0])
        PRINT2B2("after calcPenalties.  objFun =",objFunV[1])
        PRINT2B2("after calcAllPriors.  objFun =",objFunV[2])
        PRINT2B2("after calcNLLs_Recruitment.        objFun =",objFunV[3])
        PRINT2B2("after calcNLLs_ExtrapolatedEffort. objFun =",objFunV[4])
        PRINT2B2("after calcNLLs_Fisheries.       objFun =",objFunV[5])
        PRINT2B2("after calcNLLs_Surveys.         objFun =",objFunV[6])
        PRINT2B2("after calcNLLs_GrowthData.      objFun =",objFunV[7])
        PRINT2B2("after calcNLLs_ChelaHeightData. objFun =",objFunV[8])    
        PRINT2B1("----Finished calcObjFun()")
    }
}

void model_parameters::calcNLLs_ChelaHeightData(int debug, ostream& cout)
{
    if(debug>dbgObjFun) cout<<"Starting calcNLLs_ChelaHeightData()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int i=0;i<ptrMDS->nCHD;i++){
        ChelaHeightData* pCHD = ptrMDS->ppCHD[i];
        if ((pCHD->llWgt>0.0)||(debug<0)){
            if (debug<0) cout<<"`"<<pCHD->name<<"`=list("<<endl;
            int nObs = pCHD->nObs;
            //cout<<"nObs= "<<nObs<<tb;
            if (nObs>0) {
                /* index of model survey corresponding to dataset */
                int v = mapD2MChd(i+1);//
                //cout<<"v= "<<v<<endl;
                /* likelihood multiplier for this dataset */
                double wgt = pCHD->llWgt;
                /* year corresponding to observed fractions */
                ivector y_n = pCHD->obsYear_n;
                /* sample sizes from observations */
                dvector ss_n = pCHD->obsSS_n;
               /* model size bin index for size corresponding to observation */
                ivector obsIZ = pCHD->obsSizeBinIndex_n;
                /* observed fractions of new shell mature crab at size */
                dvector obsPM = 1.0*pCHD->obsPrMat_n;
                dvar_vector modPM(1,nObs);  modPM.initialize();
                dvar_vector nlls_n(1,nObs); nlls_n.initialize();
                dvector zscrs_n(1,nObs);    zscrs_n.initialize();
                for (int n=1;n<=nObs;n++){
                    //cout<<"n="<<n<<tb<<"obsIZ="<<tb<<obsIZ(n)<<tb;
                    if ((obsIZ(n)>0)&(obsIZ(n)<=nZBs)){
                        //cout<<y_n(n)<<tb;
                        dvariable nMat = n_vyxmsz(v,y_n(n),MALE,MATURE,NEW_SHELL,obsIZ(n));
                        dvariable nTot = nMat + n_vyxmsz(v,y_n(n),MALE,IMMATURE,NEW_SHELL,obsIZ(n));
                        modPM(n) = nMat/nTot;
                        //cout<<nMat<<tb<<nTot<<tb<<modPM(n)<<tb;
                        if ((modPM(n)>0.0)&(modPM(n)<1.0)){
                            if (obsPM(n)>0.0) nlls_n(n) -= ss_n(n)*obsPM(n)*(log(modPM(n))-log(obsPM(n)));
                            if (obsPM(n)<1.0) nlls_n(n) -= ss_n(n)*(1.0-obsPM(n))*(log(1.0-modPM(n))-log(1.0-obsPM(n)));
                            //cout<<nlls_n(n)<<tb;
                            double modPMv = value(modPM(n));
                            zscrs_n(n) = (obsPM(n)-modPMv)/sqrt(modPMv*(1.0-modPMv)/ss_n(n));
                            //cout<<"zscrs_n(n)="<<zscrs_n(n)<<endl;
                        }
                    }
                }//loop over n
                dvariable nll = sum(nlls_n);
                objFun += wgt*nll;
                if (debug<0) {
                    cout<<"type='binomial',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl;
                    cout<<"y=";     wts::writeToR(cout,y_n);             cout<<cc<<endl;
                    cout<<"n=";     wts::writeToR(cout,ss_n);            cout<<cc<<endl;
                    cout<<"z=";     wts::writeToR(cout,pCHD->obsSize_n); cout<<cc<<endl;
                    cout<<"i=";     wts::writeToR(cout,obsIZ);           cout<<cc<<endl;
                    cout<<"obsPM="; wts::writeToR(cout,obsPM);           cout<<cc<<endl;
                    cout<<"modPM="; wts::writeToR(cout,value(modPM));    cout<<cc<<endl;
                    cout<<"nlls=";  wts::writeToR(cout,value(nlls_n));   cout<<cc<<endl;
                    cout<<"zscrs="; wts::writeToR(cout,zscrs_n);         cout<<cc<<endl;
                    cout<<"rmse="<<sqrt(norm2(zscrs_n)/zscrs_n.size())<<"),"<<endl;
                }
            }//nObs>0
        }//((pCHD->llWgt>0)||(debug<0))
    }//datasets (i)
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug>dbgObjFun) cout<<"finished calcNLLs_ChelaHeightData()"<<endl;
}

void model_parameters::calcNLLs_GrowthData(int debug, ostream& cout)
{
    if(debug>dbgObjFun) cout<<"Starting calcNLLs_GrowthData()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int i=0;i<ptrMDS->nGrw;i++){
        GrowthData* pGD = ptrMDS->ppGrw[i];
        if ((pGD->llWgt>0.0)||(debug<0)){
            if (debug<0) cout<<"`"<<ptrMDS->ppGrw[i]->name<<"`=list("<<endl;
            for (int x=1;x<=nSXs;x++){
                int nObs = pGD->nObs_x(x);
                if (nObs>0) {
                    double wgt = pGD->llWgt;
                    /* observation year */
                    ivector year_n = pGD->obsYears_xn(x);
                    /* pre-molt size, by observation */
                    dvar_vector zpre_n = pGD->inpData_xcn(x,2);
                    /* post-molt size, by observation */
                    dvar_vector zpst_n = ptrMDS->ppGrw[i]->inpData_xcn(x,3);
                    /* molt increment, by observation */
                    dvar_vector incZ_n = zpst_n - zpre_n;
                    /* mean post-molt size, by observation */
                    //dvar_vector mnZ_n = elem_prod(mfexp(grA_xy(x)(year_n)),pow(zpre_n,grB_xy(x)(year_n)));
                    dvar_vector mnZ_n(zpre_n.indexmin(),zpre_n.indexmax());
                    if (ptrMOs->optGrowthParam==0){
                        mnZ_n = mfexp(grA_xy(x)(year_n)+elem_prod(grB_xy(x)(year_n),log(zpre_n)));
                    } else if (ptrMOs->optGrowthParam==1){
                        mnZ_n = elem_prod(
                                    grA_xy(x)(year_n),
                                    mfexp(
                                        elem_prod(
                                            elem_div(log(elem_div( grB_xy(x)(year_n), grA_xy(x)(year_n))),
                                                     log(elem_div(zGrB_xy(x)(year_n),zGrA_xy(x)(year_n)))),
                                            log(elem_div(zpre_n,zGrA_xy(x)(year_n)))
                                        )
                                    )
                                );
                    } else if (ptrMOs->optGrowthParam==2){
                        mnZ_n = elem_prod(
                                    grA_xy(x)(year_n),
                                    mfexp(
                                        elem_prod(
                                            grB_xy(x)(year_n),
                                            log(elem_div(zpre_n,zGrA_xy(x)(year_n)))
                                        )
                                    )
                                );
                    } else {
                        //throw error
                        PRINT2B1(" ")
                        PRINT2B1("#---------------------")
                        PRINT2B2("Unknown growth parameterization option",ptrMOs->optGrowthParam)
                        PRINT2B1("Terminating model run. Please correct.")
                        ad_exit(-1);
                    }
                    /* multiplicative scale factor, by observation */
                    dvar_vector ibeta_n = 1.0/grBeta_xy(x)(year_n);
                    /* location factor, by observation */
                    dvar_vector alpha_n = elem_prod(mnZ_n-zpre_n,ibeta_n);
                    dvar_vector nlls_n(1,nObs); nlls_n.initialize();
                    nlls_n = -wts::log_gamma_density(incZ_n,alpha_n,ibeta_n);
                    dvariable nll = sum(nlls_n);
                    if (isnan(value(nll))){
                        dvar_vector zscrs = elem_div((zpst_n-mnZ_n),sqrt(elem_prod(mnZ_n,grBeta_xy(x)(year_n))));
                        ofstream os("GrowthData.NLLs.NanReport.dat");
                        os.precision(12);
                        os<<"phase = "<<current_phase()<<endl;
                        os<<"sex   = "<<tcsam::getSexType(x)<<endl;
                        os<<"nll   = "<<nll<<endl;
                        os<<"nObs  = "<<nObs<<endl;
                        os<<"year  grA   zGrA    grB    zGrB   zpre_n  zpst_n  mnZ_n   incZ   mnInc  ibeta_n alpha_n nll_n  zscr"<<endl;
                        for (int n=1;n<=nObs;n++){
                            os<<year_n(n)<<tb<<
                                    grA_xy(x)(year_n(n))<<tb<<grB_xy(x)(year_n(n))<<tb<<
                                    zGrA_xy(x)(year_n(n))<<tb<<zGrB_xy(x)(year_n(n))<<tb<<
                                    zpre_n(n)<<tb<<zpst_n(n)<<tb<<mnZ_n(n)<<tb<<
                                    zpst_n(n)-zpre_n(n)<<tb<<mnZ_n(n)-zpre_n(n)<<tb<<
                                    ibeta_n(n)<<tb<<alpha_n(n)<<tb<<nlls_n(n)<<tb<<zscrs(n)<<endl;
                        }
                        os.close();
                        exit(-1);
                    }
                    objFun += wgt*nll;
                    if (debug<0) {
                        dvar_vector zscrs = elem_div((zpst_n-mnZ_n),sqrt(elem_prod(mnZ_n,grBeta_xy(x)(year_n))));
                        cout<<tcsam::getSexType(x)<<"=list(type='gamma',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl;
                        cout<<"years="; wts::writeToR(cout,year_n);         cout<<cc<<endl;
                        cout<<"zPre=";  wts::writeToR(cout,value(zpre_n));  cout<<cc<<endl;
                        cout<<"zPst=";  wts::writeToR(cout,value(zpst_n));  cout<<cc<<endl;
                        cout<<"grA=";   wts::writeToR(cout,value(grA_xy(x)(year_n))); cout<<cc<<endl;
                        cout<<"grB=";   wts::writeToR(cout,value(grB_xy(x)(year_n))); cout<<cc<<endl;
                        cout<<"mnZ =";  wts::writeToR(cout,value(mnZ_n));   cout<<cc<<endl;
                        cout<<"ibeta="; wts::writeToR(cout,value(ibeta_n)); cout<<cc<<endl;
                        cout<<"alpha="; wts::writeToR(cout,value(alpha_n)); cout<<cc<<endl;
                        cout<<"nlls=";  wts::writeToR(cout,value(nlls_n));  cout<<cc<<endl;
                        cout<<"zscrs="; wts::writeToR(cout,value(zscrs));   cout<<cc<<endl;
                        cout<<"rmse="<<sqrt(value(norm2(zscrs))/zscrs.size())<<"),"<<endl;
                    }
                }//nObs>0
            }//x
            if (debug<0) cout<<"NULL),";
        }//(pGD->llWgt>0.0)||(debug<0)
    }//datasets
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug>dbgObjFun) cout<<"finished calcNLLs_GrowthData()"<<endl;
}

void model_parameters::testNaNs(double v, adstring str)
{
    if (isnan(v)){
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"----NaN detected: "<<str<<"---"<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        ofstream os("NaNReport.rep");
        os.precision(12);
        os<<"#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        os<<"#----NaN detected: "<<str<<"---"<<endl;
        os<<"#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        updateMPI(0,cout);
        os<<"nanRep=list(phase="<<current_phase()<<",pcCtr="<<ctrProcCalls<<cc<<",pcCtrInPhs="<<ctrProcCallsInPhase<<cc<<endl;
        ReportToR_Params(os,0,cout);         os<<","<<endl;
        ReportToR_ModelProcesses(os,0,cout); os<<","<<endl;
        ReportToR_ModelResults(os,0,cout);   os<<","<<endl;
        ReportToR_ModelFits(os,-1.0,-1,cout);     os<<endl;
        os<<")"<<endl;
        os.close();
        exit(-1);
    }
}

void model_parameters::calcNorm2NLL(double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNorm2NLL()"<<endl;
    int y;
    double rmse = 0.0; int cnt = 0;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    for (int i=1;i<=yrs.size();i++){
        y = yrs(i);
        if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
            zscr(y) = (obs[i]-mod[y]); cnt++;
        }
    }
    if (cnt>0) rmse = sqrt(value(norm2(zscr))/cnt);
    nll += 0.5*norm2(zscr);
    objFun += wgt*nll;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='norm2',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<cc<<endl;
        cout<<"rmse="<<rmse<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcNorm2NLL()"<<endl;
}

void model_parameters::calcNormalNLL(double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
   if (debug>=dbgAll) cout<<"Starting calcNormalNLL()"<<endl;
    int y;
    double rmse = 0.0; int cnt = 0;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    if (sum(stdv)>0){
        for (int i=1;i<=yrs.size();i++){
            y = yrs(i);
            if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
                zscr(y) = (obs[i]-mod[y])/stdv[i]; cnt++;
            }
        }
        if (cnt>0) rmse = sqrt(value(norm2(zscr))/cnt);
        nll += 0.5*norm2(zscr);
    }
    objFun += wgt*nll;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='normal',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<cc<<endl;
        cout<<"rmse="<<rmse<<")";
    }
   if (debug>=dbgAll) cout<<"Finished calcNormalNLL()"<<endl;
}

void model_parameters::calcLognormalNLL(double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcLognormalNLL()"<<endl;
    int y;
    dvariable nll = 0.0;
    double rmse = 0.0; int cnt = 0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    if (sum(stdv)>0){
        for (int i=1;i<=yrs.size();i++){
            y = yrs(i);
            if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
                zscr(y) = (log(obs[i]+smlVal)-log(mod[y]+smlVal))/stdv[i]; cnt++;
            }
        }
        if (cnt>0) {
            rmse = sqrt(value(norm2(zscr))/cnt);
            nll = 0.5*norm2(zscr);
            objFun += wgt*nll;
        }
    }
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='lognormal',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<cc<<endl;
        cout<<"rmse="<<rmse<<")";
    }
   if (debug>=dbgAll) cout<<"Finished calcLognormalNLL()"<<endl;
}

void model_parameters::calcNoneNLL(double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNoneNLL(agg catch)"<<endl;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='none',wgt="<<0<<cc<<"nll="<<0<<cc<<"objfun="<<0<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcNoneNLL(agg catch)"<<endl;
}

void model_parameters::calcMultinomialNLL(double wgt, dvar_vector& mod, dvector& obs, double& ss, int& yr, int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcMultinomialNLL()"<<endl;
    dvariable nll = -ss*(obs*(log(mod+smlVal)-log(obs+smlVal)));//note dot-product sums
    objFun += wgt*nll;
    if (debug<0){
        dvector vmod = value(mod);
        dvector nlls = -ss*(elem_prod(obs,log(vmod+smlVal)-log(obs+smlVal)));
        dvector zscrs = elem_div(obs-vmod,sqrt(elem_prod((vmod+smlVal),1.0-(vmod+smlVal))/ss));//pearson residuals
        double effN = 0.0;
        if ((ss>0)&&(norm2(obs-vmod)>0)) effN = (vmod*(1.0-vmod))/norm2(obs-vmod);
        cout<<"list(nll.type='multinomial',yr="<<yr<<cc<<"wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<"ss="<<ss<<cc<<"effN="<<effN<<cc<<endl; 
        adstring dzbs = "size=c("+ptrMC->csvZBs+")";
        cout<<"nlls=";  wts::writeToR(cout,nlls, dzbs); cout<<cc<<endl;
        cout<<"obs=";   wts::writeToR(cout,obs,  dzbs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,vmod, dzbs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,zscrs,dzbs); cout<<endl;
        cout<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcMultinomialNLL()"<<endl;
}

void model_parameters::calcNoneNLL(double wgt, dvar_vector& mod, dvector& obs, double& ss, int& yr, int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNoneNLL(size comps)"<<endl;
    if (debug<0){
        cout<<"list(nll.type='none',yr="<<yr<<cc<<"wgt="<<wgt<<cc<<"nll="<<0.0<<cc<<"objfun="<<0.0<<cc<<"ss="<<0<<cc<<"effN="<<0<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcNoneNLL(size comps)"<<endl;
}

void model_parameters::calcNLL(int llType, double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
    switch (llType){
        case tcsam::LL_NONE:
            calcNoneNLL(wgt,mod,obs,stdv,yrs,debug,cout);
            break;
        case tcsam::LL_LOGNORMAL:
            calcLognormalNLL(wgt,mod,obs,stdv,yrs,debug,cout);
            break;
        case tcsam::LL_NORMAL:
            calcNormalNLL(wgt,mod,obs,stdv,yrs,debug,cout);
            break;
        case tcsam::LL_NORM2:
            calcNorm2NLL(wgt,mod,obs,stdv,yrs,debug,cout);
            break;
        default:
            cout<<"Unrecognized likelihood type in calcNLL(1)"<<endl;
            cout<<"Input type was "<<llType<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
    }    
}

void model_parameters::calcNLL(int llType, double wgt, dvar_vector& mod, dvector& obs, double& ss, int& yr, int debug, ostream& cout)
{
    switch (llType){
        case tcsam::LL_NONE:
            calcNoneNLL(wgt,mod,obs,ss,yr,debug,cout);
            break;
        case tcsam::LL_MULTINOMIAL:
            calcMultinomialNLL(wgt,mod,obs,ss,yr,debug,cout);
            break;
        default:
            cout<<"Unrecognized likelihood type in calcNLL(2)"<<endl;
            cout<<"Input type was "<<llType<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
    }
}

void model_parameters::calcNLLs_AggregateCatch(AggregateCatchData* ptrAB, dvar5_array& mA_yxmsz, int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNLLs_AggregateCatch()"<<endl;
    int mny = mA_yxmsz.indexmin();
    int mxy = mA_yxmsz.indexmax();//may NOT be mxYr
    dvar_vector tAB_y(mny,mxy);
    int isBio = ptrAB->type==AggregateCatchData::KW_BIOMASS_DATA;
    if (debug>=dbgAll) cout<<"isBio="<<isBio<<tb<<"type="<<ptrAB->type<<endl;
    if (debug<0) cout<<"list(fit.type='"<<tcsam::getFitType(ptrAB->optFit)<<"',fits=list("<<endl;
    if (ptrAB->optFit==tcsam::FIT_BY_TOT){
        tAB_y.initialize();
        if (isBio){
            for (int x=1;x<=nSXs;x++) {
                for (int m=1;m<=nMSs;m++) {
                    if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                    for (int s=1;s<=nSCs;s++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//s
                }//m
            }//x
        } else {
            for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y));//sum over x,m,s,z
        }
        if (debug>=dbgAll) cout<<"FIT_BY_TOT: "<<tAB_y<<endl;
        if (debug<0) {
            cout<<"list(";
            cout<<"x="<<qt<<tcsam::getSexType(ALL_SXs)     <<qt<<cc;
            cout<<"m="<<qt<<tcsam::getMaturityType(ALL_MSs)<<qt<<cc;
            cout<<"s="<<qt<<tcsam::getShellType(ALL_SCs)   <<qt<<cc;
            cout<<"nll=";
        }
        calcNLL(ptrAB->llType, ptrAB->llWgt, tAB_y, ptrAB->C_xmsy(ALL_SXs,ALL_MSs,ALL_SCs), ptrAB->sd_xmsy(ALL_SXs,ALL_MSs,ALL_SCs), ptrAB->yrs, debug, cout);                
        if (debug<0) cout<<")";
        if (debug<0) cout<<")";
    } else if (ptrAB->optFit==tcsam::FIT_BY_X){
        for (int x=1;x<=nSXs;x++){
            tAB_y.initialize();
            if (isBio){
                for (int m=1;m<=nMSs;m++) {
                    if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                    for (int s=1;s<=nSCs;s++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//s
                }//m
            } else {
                for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y,x));//sum over m,s,z
            }
            if (debug>=dbgAll) cout<<"FIT_BY_X("<<x<<"): "<<tAB_y<<endl;
            if (debug<0) {
                cout<<"list(";
                cout<<"x="<<qt<<tcsam::getSexType(x)           <<qt<<cc;
                cout<<"m="<<qt<<tcsam::getMaturityType(ALL_MSs)<<qt<<cc;
                cout<<"s="<<qt<<tcsam::getShellType(ALL_SCs)   <<qt<<cc;
                cout<<"nll=";
            }
            calcNLL(ptrAB->llType, ptrAB->llWgt, tAB_y, ptrAB->C_xmsy(x,ALL_MSs,ALL_SCs), ptrAB->sd_xmsy(x,ALL_MSs,ALL_SCs), ptrAB->yrs, debug, cout); 
            if (debug<0) cout<<"),"<<endl;
        }//x
        if (debug<0) cout<<"NULL)";
    } else if ((ptrAB->optFit==tcsam::FIT_BY_XM)||(ptrAB->optFit==tcsam::FIT_BY_X_MATONLY)){
        int mnM = IMMATURE;
        if (ptrAB->optFit==tcsam::FIT_BY_X_MATONLY) mnM = MATURE;
        for (int x=1;x<=nSXs;x++){
            for (int m=mnM;m<=nMSs;m++){
                tAB_y.initialize();
                if (isBio){
                    if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                    for (int s=1;s<=nSCs;s++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//s
                } else {
                    for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y,x,m));//sum over s,z
                }
                if (debug<0) {
                    cout<<"list(";
                    cout<<"x="<<qt<<tcsam::getSexType(x)        <<qt<<cc;
                    cout<<"m="<<qt<<tcsam::getMaturityType(m)   <<qt<<cc;
                    cout<<"s="<<qt<<tcsam::getShellType(ALL_SCs)<<qt<<cc;
                    cout<<"nll=";
                }
                calcNLL(ptrAB->llType, ptrAB->llWgt, tAB_y, ptrAB->C_xmsy(x,m,ALL_SCs), ptrAB->sd_xmsy(x,m,ALL_SCs), ptrAB->yrs, debug, cout); 
                if (debug<0) cout<<"),"<<endl;
            }//m
        }//x
        if (debug<0) cout<<"NULL)";
    } else if (ptrAB->optFit==tcsam::FIT_BY_XS){
        for (int x=1;x<=nSXs;x++){
            for (int s=1;s<=nSCs;s++){
                tAB_y.initialize();
                if (isBio){
                    for (int m=1;m<=nMSs;m++) {
                        if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//m
                } else {
                    for (int m=1;m<=nMSs;m++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += sum(mA_yxmsz(y,x,m,s));//sum over m,z
                        }//y
                    }//m
                }
                if (debug<0) {
                    cout<<"list(";
                    cout<<"x="<<qt<<tcsam::getSexType(x)           <<qt<<cc;
                    cout<<"m="<<qt<<tcsam::getMaturityType(ALL_MSs)<<qt<<cc;
                    cout<<"s="<<qt<<tcsam::getShellType(s)         <<qt<<cc;
                    cout<<"nll=";
                }
                calcNLL(ptrAB->llType, ptrAB->llWgt, tAB_y, ptrAB->C_xmsy(x,ALL_MSs,s), ptrAB->sd_xmsy(x,ALL_MSs,s), ptrAB->yrs, debug, cout); 
                if (debug<0) cout<<"),"<<endl;
            }//s
        }//x
        if (debug<0) cout<<"NULL)";
    } else if (ptrAB->optFit==tcsam::FIT_BY_XMS){
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++){
                    tAB_y.initialize();
                    if (isBio){
                        if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                        for (int y=mny;y<=mxy;y++) tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                    } else {
                        for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y,x,m,s));//sum over z
                    }
                    if (debug<0) {
                        cout<<"list(";
                        cout<<"x="<<qt<<tcsam::getSexType(x)     <<qt<<cc;
                        cout<<"m="<<qt<<tcsam::getMaturityType(m)<<qt<<cc;
                        cout<<"s="<<qt<<tcsam::getShellType(s)   <<qt<<cc;
                        cout<<"nll=";
                    }
                    calcNLL(ptrAB->llType, ptrAB->llWgt, tAB_y, ptrAB->C_xmsy(x,m,s), ptrAB->sd_xmsy(x,m,s), ptrAB->yrs, debug, cout); 
                    if (debug<0) cout<<"),"<<endl;
                }//s
            }//m
        }//x
        if (debug<0) cout<<"NULL)";
    } else {
        std::cout<<"Calling calcNLLs_AggregateCatch with invalid fit option."<<endl;
        std::cout<<"Invalid fit option was '"<<tcsam::getFitType(ptrAB->optFit)<<qt<<endl;
        std::cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug<0) cout<<")";
    if (debug>=dbgAll){
        cout<<"Finished calcNLLs_AggregateCatch()"<<endl;
    }
}

d5_array model_parameters::calcNLLs_CatchNatZ(SizeFrequencyData* ptrZFD, dvar5_array& mA_yxmsz, int debug, ostream& cout)
{
    if (debug>=dbgNLLs) cout<<"Starting calcNLLs_CatchNatZ()"<<endl;
    d5_array effWgtComps_xmsyn;
    if (ptrZFD->optFit==tcsam::FIT_NONE) return effWgtComps_xmsyn;
    ivector yrs = ptrZFD->yrs;
    int y;
    double ss;
    dvariable nT;
    int mny = mA_yxmsz.indexmin();
    int mxy = mA_yxmsz.indexmax();//may NOT be mxYr
    dvector     oP_z;//observed size comp.
    dvar_vector mP_z;//model size comp.
    if (debug<0) cout<<"list("<<endl;
    if (ptrZFD->optFit==tcsam::FIT_BY_TOT){
        effWgtComps_xmsyn.allocate(tcsam::ALL_SXs,tcsam::ALL_SXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                ss = 0;
                nT = sum(mA_yxmsz(y));//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=ALL_SXs;x++){
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                    }
                    if (sum(oP_z)>0) oP_z /= sum(oP_z);
                    if (debug>=dbgNLLs){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    for (int x=1;x<=nSXs;x++){
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++)mP_z += mA_yxmsz(y,x,m,s);
                        }
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                    if (debug<0) {
                        cout<<"'"<<y<<"'=list(";
                        cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                        cout<<"y="<<y<<cc;
                        cout<<"x='"<<tcsam::getSexType(ALL_SXs)<<"'"<<cc;
                        cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                        cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                        cout<<"fit=";
                    }
                    calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,y,debug,cout);
                    effWgtComps_xmsyn(tcsam::ALL_SXs,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                    if (debug<0) cout<<")"<<cc<<endl;
                }//value(nT)>0
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_TOT
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_X){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                        if (sum(oP_z)>0) oP_z /= sum(oP_z);
                        if (debug>=dbgNLLs){
                            cout<<"ss = "<<ss<<endl;
                            cout<<"oP_Z = "<<oP_z<<endl;
                        }
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++) mP_z += mA_yxmsz(y,x,m,s);
                        }
                        mP_z /= nT;//normalize model size comp
                        if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                        if (debug<0) {
                            cout<<"'"<<y<<"'=list(";
                            cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                            cout<<"y="<<y<<cc;
                            cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                            cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                            cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                            cout<<"fit=";
                        }
                        calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,y,debug,cout);
                        effWgtComps_xmsyn(x,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//nT>0
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_X
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XM){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               1,nMSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++){
                        ss = 0;
                        nT = sum(mA_yxmsz(y,x,m));//=0 if not calculated
                        if (value(nT)>0){
                            oP_z.initialize();//observed size comp.
                            mP_z.initialize();//model size comp.
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            if (sum(oP_z)>0) oP_z /= sum(oP_z);
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            for (int s=1;s<=nSCs;s++) mP_z += mA_yxmsz(y,x,m,s);
                            mP_z /= nT;//normalize model size comp
                            if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,y,debug,cout);
                            effWgtComps_xmsyn(x,m,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//nT>0
                    }//m
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XM
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XS){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               1,nSCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                for (int x=1;x<=nSXs;x++) {
                    for (int s=1;s<=nSCs;s++){
                        ss = 0;
                        nT.initialize();
                        for (int m=1;m<=nMSs;m++) nT += sum(mA_yxmsz(y,x,m,s));//=0 if not calculated
                        if (value(nT)>0){
                            oP_z.initialize();//observed size comp.
                            mP_z.initialize();//model size comp.
                            for (int m=1;m<=ALL_MSs;m++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            if (sum(oP_z)>0) oP_z /= sum(oP_z);
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            for (int m=1;m<=nMSs;m++) mP_z += mA_yxmsz(y,x,m,s);
                            mP_z /= nT;//normalize model size comp
                            if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,y,debug,cout);
                            effWgtComps_xmsyn(x,tcsam::ALL_MSs,s,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//nT>0
                    }//s
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XS
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XMS){
        effWgtComps_xmsyn.allocate(1,nSXs,1,nMSs,1,nSCs,mny,mxy,1,3);
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++) {
                            ss = 0;
                            nT = sum(mA_yxmsz(y,x,m,s));//=0 if not calculated
                            if (value(nT)>0){
                                oP_z.initialize();//observed size comp.
                                mP_z.initialize();//model size comp.                            
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                                if (sum(oP_z)>0) oP_z /= sum(oP_z);
                                if (debug>=dbgNLLs){
                                    cout<<"ss = "<<ss<<endl;
                                    cout<<"oP_Z = "<<oP_z<<endl;
                                }
                                mP_z += mA_yxmsz(y,x,m,s);
                                mP_z /= nT;//normalize model size comp
                                if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                                if (debug<0) {
                                    cout<<"'"<<y<<"'=list(";
                                    cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                    cout<<"y="<<y<<cc;
                                    cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                    cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                    cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                    cout<<"fit=";
                                }
                                calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,y,debug,cout);
                                effWgtComps_xmsyn(x,m,s,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                                if (debug<0) cout<<")"<<cc<<endl;
                            }//nT>0
                        }//s
                    }//m
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XMS
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XE){
        effWgtComps_xmsyn.allocate(tcsam::ALL_SXs,tcsam::ALL_SXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nSXs*nZBs);
        mP_z.allocate(1,nSXs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                ss = 0;
                nT = sum(mA_yxmsz(y));//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=nSXs;x++) {
                        int mnz = 1+(x-1)*nZBs;
                        int mxz = x*nZBs;
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++) mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                        }
                    }//x
                    if (sum(oP_z)>0) oP_z /= sum(oP_z);
                    if (debug>=dbgNLLs){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                    for (int x=1;x<=nSXs;x++) {
                        int mnz = 1+(x-1)*nZBs;
                        int mxz = x*nZBs;
                        dvar_vector mPt = mP_z(mnz,mxz);
                        dvector oPt = oP_z(mnz,mxz);
                        if (debug<0) {
                            cout<<"'"<<y<<"'=list(";
                            cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                            cout<<"y="<<y<<cc;
                            cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                            cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                            cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                            cout<<"fit=";
                        }
                        calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//x
                    effWgtComps_xmsyn(tcsam::ALL_SXs,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                }//nT>0
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XE
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_X_ME){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nMSs*nZBs);
        mP_z.allocate(1,nMSs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs;
                            int mxz = m*nZBs;
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }//s
                            for (int s=1;s<=nSCs;s++) {
                                mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                            }//s
                        }//m
                        if (sum(oP_z)>0) oP_z /= sum(oP_z);
                        if (debug>=dbgNLLs){
                            cout<<"ss = "<<ss<<endl;
                            cout<<"oP_Z = "<<oP_z<<endl;
                        }
                        mP_z /= nT;//normalize model size comp
                        if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs;
                            int mxz = m*nZBs;
                            dvar_vector mPt = mP_z(mnz,mxz);
                            dvector oPt = oP_z(mnz,mxz);
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//m
                        effWgtComps_xmsyn(x,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                    }//nT>0
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_X_ME
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_X_SE){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nSCs*nZBs);
        mP_z.allocate(1,nSCs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int s=1;s<=nSCs;s++) {
                            int mnz = 1+(s-1)*nZBs;
                            int mxz = s*nZBs;
                            for (int m=1;m<=ALL_MSs;m++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }//m
                            for (int m=1;m<=nMSs;m++) {
                                mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                            }//m
                        }//s
                        if (sum(oP_z)>0) oP_z /= sum(oP_z);
                        if (debug>=dbgNLLs){
                            cout<<"ss = "<<ss<<endl;
                            cout<<"oP_Z = "<<oP_z<<endl;
                        }
                        mP_z /= nT;//normalize model size comp
                        if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                        for (int s=1;s<=nSCs;s++) {
                            int mnz = 1+(s-1)*nZBs;
                            int mxz = s*nZBs;
                            dvar_vector mPt = mP_z(mnz,mxz);
                            dvector oPt = oP_z(mnz,mxz);
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//s
                        effWgtComps_xmsyn(x,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                    }//nT>0
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_X_SE
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_XME){
        effWgtComps_xmsyn.allocate(tcsam::ALL_SXs,tcsam::ALL_SXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nSXs*nMSs*nZBs);
        mP_z.allocate(1,nSXs*nMSs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                ss = 0;
                nT = sum(mA_yxmsz(y));//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=nSXs;x++) {
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs+(x-1)*nMSs*nZBs;
                            int mxz = mnz+nZBs-1;
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            if (m<=nMSs) {for (int s=1;s<=nSCs;s++) mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);}
                        }//m
                    }//x
                    if (sum(oP_z)>0) oP_z /= sum(oP_z);
                    if (debug>=dbgNLLs){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                    for (int x=1;x<=nSXs;x++) {
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs+(x-1)*nMSs*nZBs;
                            int mxz = mnz+nZBs-1;
                            dvar_vector mPt = mP_z(mnz,mxz);
                            dvector oPt = oP_z(mnz,mxz);
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//m
                    }//x
                    effWgtComps_xmsyn(tcsam::ALL_SXs,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                }//nT>0
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XME
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XM_SE){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               1,nMSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nSCs*nZBs);
        mP_z.allocate(1,nSCs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++) {
                        ss = 0;
                        nT = sum(mA_yxmsz(y,x,m));//=0 if not calculated
                        if (value(nT)>0){
                            oP_z.initialize();//observed size comp.
                            mP_z.initialize();//model size comp.
                            for (int s=1;s<=nSCs;s++) {
                                int mnz = 1+(s-1)*nZBs;
                                int mxz = s*nZBs;
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                                mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                            }//s
                            if (sum(oP_z)>0) oP_z /= sum(oP_z);
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            mP_z /= nT;//normalize model size comp
                            if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                            for (int s=1;s<=nSCs;s++) {
                                int mnz = 1+(s-1)*nZBs;
                                int mxz = s*nZBs;
                                dvar_vector mPt = mP_z(mnz,mxz);
                                dvector oPt = oP_z(mnz,mxz);
                                if (debug<0) {
                                    cout<<"'"<<y<<"'=list(";
                                    cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                    cout<<"y="<<y<<cc;
                                    cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                    cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                    cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                    cout<<"fit=";
                                }
                                calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                                if (debug<0) cout<<")"<<cc<<endl;
                            }//s
                            effWgtComps_xmsyn(x,m,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                        }//nT>0
                    }//m
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XM_SE
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_X_MSE){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nMSs*nSCs*nZBs);
        mP_z.allocate(1,nMSs*nSCs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++) {
                                int mnz = 1+(s-1)*nZBs+(m-1)*nSCs*nZBs;
                                int mxz =       s*nZBs+(m-1)*nSCs*nZBs;
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                                mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                            }//s
                        }//m
                        if (sum(oP_z)>0) oP_z /= sum(oP_z);
                        if (debug>=dbgNLLs){
                            cout<<"ss = "<<ss<<endl;
                            cout<<"oP_Z = "<<oP_z<<endl;
                        }
                        mP_z /= nT;//normalize model size comp
                        if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++) {
                                int mnz = 1+(s-1)*nZBs+(m-1)*nSCs*nZBs;
                                int mxz =       s*nZBs+(m-1)*nSCs*nZBs;
                                dvar_vector mPt = mP_z(mnz,mxz);
                                dvector oPt = oP_z(mnz,mxz);
                                if (debug<0) {
                                    cout<<"'"<<y<<"'=list(";
                                    cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                    cout<<"y="<<y<<cc;
                                    cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                    cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                    cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                    cout<<"fit=";
                                }
                                calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                                if (debug<0) cout<<")"<<cc<<endl;
                            }//s
                        }//m
                        effWgtComps_xmsyn(x,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                    }//nT>0
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_X_MSE
    } else 
    {
        std::cout<<"Calling calcNLLs_CatchNatZ with invalid fit option."<<endl;
        std::cout<<"Invalid fit option was '"<<tcsam::getFitType(ptrZFD->optFit)<<qt<<endl;
        std::cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug<0) cout<<"NULL)";
    if (debug>=dbgNLLs) cout<<"Finished calcNLLs_CatchNatZ()"<<endl;
    return effWgtComps_xmsyn;
}

dvector model_parameters::calcEffWgtComponents(double ss, dvector& obs, dvar_vector& mod, int debug, ostream& cout)
{
    //if (debug) cout<<"Starting calcEffWgtComponents(...)"<<endl;
    dvector effWgts(1,3); effWgts = 0.0;
    dvector vmd = value(mod);
    //set counter for valid comps
    effWgts(1) = 1.0;
    //McAllister-Ianelli E_y/N_y ala Punt, 2017
    effWgts(2) = ((vmd*(1.0-vmd))/norm2(obs-vmd))/ss;
    //Francis weights ala Punt, 2017
    int N = obs.size();
    dvector L(1,N); L.fill_seqadd(1.0,1.0);
    double obsMnL = L*obs;
    double prdMnL = L*vmd;
    double stdMnL = sqrt(vmd*elem_prod(L-prdMnL,L-prdMnL)/N);
    effWgts(3) = (obsMnL-prdMnL)/stdMnL;
    //if (debug) cout<<"Starting calcEffWgtComponents(...)"<<endl;
    return effWgts;
}

void model_parameters::calcNLLs_Fisheries(int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNLLs_Fisheries()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int f=1;f<=nFsh;f++){
        if (debug>=dbgAll) cout<<"calculating NLLs for fishery "<<ptrMC->lblsFsh[f]<<endl;
        if (debug<0) cout<<"`"<<ptrMC->lblsFsh[f]<<"`=list("<<endl;
        int fd = mapM2DFsh(f);//get index for fishery data corresponding to model fishery f
        FleetData* ptrObs = ptrMDS->ppFsh[fd-1];
        if (ptrObs->hasRCD){//retained catch data
            if (debug<0) cout<<"retained.catch=list("<<endl;
            if (ptrObs->ptrRCD->hasN && ptrObs->ptrRCD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---retained catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrRCD->ptrN,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasB && ptrObs->ptrRCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---retained catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrRCD->ptrB,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasZFD && ptrObs->ptrRCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---retained catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                smlVal = 0.001;//match to TCSAM2013
                calcNLLs_CatchNatZ(ptrObs->ptrRCD->ptrZFD,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (ptrObs->hasTCD){//observed total catch data
            if (debug<0) cout<<"total.catch=list("<<endl;
            if (ptrObs->ptrTCD->hasN && ptrObs->ptrTCD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---total catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrTCD->ptrN,cpN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasB && ptrObs->ptrTCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---total catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrTCD->ptrB,cpN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasZFD && ptrObs->ptrTCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---total catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                smlVal = 0.001;//match to TCSAM2013
                calcNLLs_CatchNatZ(ptrObs->ptrTCD->ptrZFD,cpN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (ptrObs->hasDCD){//observed discard catch data
            if (debug<0) cout<<"discard.catch=list("<<endl;
            if (ptrObs->ptrDCD->hasN && ptrObs->ptrDCD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---discard catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrDCD->ptrN,dsN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasB && ptrObs->ptrDCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---discard catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrDCD->ptrB,dsN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasZFD && ptrObs->ptrDCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---discard catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                smlVal = 0.001;//match to TCSAM2013
                calcNLLs_CatchNatZ(ptrObs->ptrDCD->ptrZFD,dsN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (debug<0) cout<<"NULL),"<<endl;
    }//fisheries
    if (debug<0) cout<<"NULL)"<<endl;
    if (!debug) testNaNs(value(objFun),"in calcNLLs_Fisheries()");
    if (debug>=dbgAll) cout<<"Finished calcNLLs_Fisheries()"<<endl;
}

void model_parameters::calcNLLs_ExtrapolatedEffort(int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNLLs_ExtrapolatedEffort()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    CapRateAvgScenarios* ptrCRASs = ptrMOs->ptrEffXtrapScenarios->ptrCapRateAvgScenarios;
    double stdv = sqrt(log(1.0+square(0.1)));//TODO: check this!!
    int mnx, mxx, mnm, mxm, mns, mxs;
    nllEffX_nxms.initialize();
    zscrEffX_nxmsy.initialize();
    for (int n=1;n<=ptrCRASs->nAvgs;n++){
        CapRateAvgScenario* ptrCRAS = ptrCRASs->ppCRASs[n-1];
        int idEAS = ptrCRAS->idEffAvgInfo;//index to associated average effort
        int f     = ptrCRAS->f;//fishery
        mnx = mxx = ptrCRAS->x;//sex index
        if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
        mnm = mxm = ptrCRAS->m;//maturity index
        if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
        mns = mxs = ptrCRAS->s;//shell index
        if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
        for (int x=mnx;x<=mxx;x++){
            for (int m=mnm;m<=mxm;m++){
                for (int s=mns;s<=mxs;s++){
                    zscrEffX_nxmsy(n,x,m,s) = (log(obsEff_nxmsy(n,x,m,s)+smlVal)-log(prdEff_nxmsy(n,x,m,s)+smlVal))/stdv;
                    nllEffX_nxms(n,x,m,s)   = norm2(zscrEffX_nxmsy(n,x,m,s));
                    if (debug>dbgAll){
                        cout<<"n    x   m    s      stdv   wgt     nll"<<endl;
                        cout<<n<<tb<<x<<tb<<m<<tb<<s<<tb<<stdv<<tb<<ptrCRAS->llWgt<<tb<<nllEffX_nxms(n,x,m,s)<<endl;
                        cout<<"obsEff_nxmsy(n,x,m,s) = "<<obsEff_nxmsy(n,x,m,s)<<endl;
                        cout<<"prdEff_nxmsy(n,x,m,s) = "<<prdEff_nxmsy(n,x,m,s)<<endl;
                        cout<<"zscores               = "<<zscrEffX_nxmsy(n,x,m,s)<<endl;
                    }
                    if ((ptrCRAS->idParam)&&(active(pLnEffX(ptrCRAS->idParam)))){
                        objFun += (ptrCRAS->llWgt)*nllEffX_nxms(n,x,m,s);
                    } else {
                        objFun += 0.0;//add nothing, for now
                    }
                }//s
            }//m
        }//x
        if (debug<0){
            cout<<"`"<<n<<"`=list(n="<<n<<",f='"<<ptrMC->lblsFsh(f)<<"',stdv="<<stdv<<",wgt="<<ptrCRAS->llWgt<<cc<<endl;
            cout<<"obsFc=";   wts::writeToR(cout,obsFc_nxmsy(n),   xDms,mDms,sDms,yDms);cout<<cc<<endl;
            cout<<"prdFc=";   wts::writeToR(cout,prdFc_nxmsy(n),   xDms,mDms,sDms,yDms);cout<<cc<<endl;
            cout<<"obsEff=";  wts::writeToR(cout,obsEff_nxmsy(n),  xDms,mDms,sDms,yDms);cout<<cc<<endl;
            cout<<"prdEff=";  wts::writeToR(cout,prdEff_nxmsy(n),  xDms,mDms,sDms,yDms);cout<<cc<<endl;
            cout<<"zscores="; wts::writeToR(cout,zscrEffX_nxmsy(n),xDms,mDms,sDms,yDms);cout<<cc<<endl;
            cout<<"nlls=";    wts::writeToR(cout,nllEffX_nxms(n),  xDms,mDms,sDms);cout<<endl;
            cout<<")"<<cc<<endl;
        }
    }//n
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug>=dbgAll) cout<<"Finished calcNLLs_ExtrapolatedEffort()"<<endl;
}

void model_parameters::calcNLLs_Surveys(int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNLLs_Surveys()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int v=1;v<=nSrv;v++){
        if (debug>=dbgAll) cout<<"calculating NLLs for survey "<<ptrMC->lblsSrv[v]<<endl;
        if (debug<0) cout<<"`"<<ptrMC->lblsSrv[v]<<"`=list("<<endl;
        FleetData* ptrObs = ptrMDS->ppSrv[v-1];
        if (ptrObs->hasICD){//index catch data
            if (debug<0) cout<<"index.catch=list("<<endl;
            if (ptrObs->ptrICD->hasN && ptrObs->ptrICD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---index catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                smlVal = 0.000001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrICD->ptrN,n_vyxmsz(v),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrICD->hasB && ptrObs->ptrICD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---index catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                smlVal = 0.000001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrICD->ptrB,n_vyxmsz(v),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrICD->hasZFD && ptrObs->ptrICD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---index catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                smlVal = 0.001;//match to TCSAM2013
                calcNLLs_CatchNatZ(ptrObs->ptrICD->ptrZFD,n_vyxmsz(v),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (debug<0) cout<<"NULL),"<<endl;
    }//surveys loop
    if (debug<0) cout<<"NULL)"<<endl;
    if (!debug) testNaNs(value(objFun),"in calcNLLs_Surveys()");
    if (debug>=dbgAll) cout<<"Finished calcNLLs_Surveys()"<<endl;
}

void model_parameters::calcAllPriors(int debug, ostream& cout)
{
    if (debug>=dbgPriors) cout<<"Starting calcAllPriors()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    //recruitment parameters
    if (debug<0) cout<<tb<<"recruitment=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnR,pLnR,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pRCV,pRCV,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pRX,pRX,  debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pRa, pRa, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pRb, pRb, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pDevsLnR,devsLnR,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    //natural mortality parameters
    if (debug<0) cout<<tb<<"'natural mortality'=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pM,   pM,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pDM1, pDM1, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pDM2, pDM2, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pDM3, pDM3, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pDM4, pDM4, debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    //growth parameters
    if (debug<0) cout<<tb<<"growth=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGrw->pGrA,   pGrA,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGrw->pGrB,   pGrB,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGrw->pGrBeta,pGrBeta,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    //maturity parameters
    if (debug<0) cout<<tb<<"maturity=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrM2M->pvLgtPrM2M,pvLgtPrM2M,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    //selectivity parameters
    if (debug<0) cout<<tb<<"'selectivity functions'=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS1,pS1,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS2,pS2,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS3,pS3,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS4,pS4,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS5,pS5,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS6,pS6,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS1,devsS1,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS2,devsS2,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS3,devsS3,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS4,devsS4,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS5,devsS5,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS6,devsS6,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    //fishing mortality parameters
    if (debug<0) cout<<tb<<"fisheries=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnC, pLnC, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDC1, pDC1, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDC2, pDC2, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDC3, pDC3, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDC4, pDC4, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDevsLnC,devsLnC,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    //survey catchability parameters
    if (debug<0) cout<<tb<<"surveys=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pQ,   pQ,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pDQ1, pDQ1, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pDQ2, pDQ2, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pDQ3, pDQ3, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pDQ4, pDQ4, debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<endl;
    if (debug<0) cout<<")"<<endl;
    if (!debug) testNaNs(value(objFun),"in calcAllPriors()");
    if (debug>=dbgPriors) cout<<"Finished calcAllPriors()"<<endl;
}

void model_parameters::setInitVals(int debug, ostream& os)
{
    //recruitment parameters
    setInitVals(ptrMPI->ptrRec->pLnR, pLnR, usePin, debug, os);
    setInitVals(ptrMPI->ptrRec->pRCV, pRCV, usePin, debug, os);
    setInitVals(ptrMPI->ptrRec->pRX,  pRX,  usePin, debug, os);
    setInitVals(ptrMPI->ptrRec->pRa,  pRa,  usePin, debug, os);
    setInitVals(ptrMPI->ptrRec->pRb,  pRb,  usePin, debug, os);
    setInitVals(ptrMPI->ptrRec->pDevsLnR,pDevsLnR,usePin, debug, os);
    //natural mortality parameters
    setInitVals(ptrMPI->ptrNM->pM,   pM,   usePin, debug, os);
    setInitVals(ptrMPI->ptrNM->pDM1, pDM1, usePin, debug, os);
    setInitVals(ptrMPI->ptrNM->pDM2, pDM2, usePin, debug, os);
    setInitVals(ptrMPI->ptrNM->pDM3, pDM3, usePin, debug, os);
    setInitVals(ptrMPI->ptrNM->pDM4, pDM4, usePin, debug, os);
    //growth parameters
    setInitVals(ptrMPI->ptrGrw->pGrA,   pGrA,   usePin, debug, os);
    setInitVals(ptrMPI->ptrGrw->pGrB,   pGrB,   usePin, debug, os);
    setInitVals(ptrMPI->ptrGrw->pGrBeta,pGrBeta,usePin, debug, os);
    //maturity parameters
    setInitVals(ptrMPI->ptrM2M->pvLgtPrM2M,pvLgtPrM2M,usePin, debug, os);
    //selectivity parameters
    setInitVals(ptrMPI->ptrSel->pS1, pS1,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pS2, pS2,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pS3, pS3,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pS4, pS4,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pS5, pS5,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pS6, pS6,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS1, pDevsS1,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS2, pDevsS2,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS3, pDevsS3,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS4, pDevsS4,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS5, pDevsS5,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS6, pDevsS6,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pvNPSel, pvNPSel,usePin, debug, os);
    //fully-selected fishing capture rate parameters
    setInitVals(ptrMPI->ptrFsh->pHM,  pHM,  usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pLnC, pLnC, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pDC1, pDC1, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pDC2, pDC2, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pDC3, pDC3, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pDC4, pDC4, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pDevsLnC,pDevsLnC,usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pLnEffX, pLnEffX, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pLgtRet, pLgtRet, usePin, debug, os);
    //survey catchability parameters
    setInitVals(ptrMPI->ptrSrv->pQ,   pQ,   usePin, debug, os);
    setInitVals(ptrMPI->ptrSrv->pDQ1, pDQ1, usePin, debug, os);
    setInitVals(ptrMPI->ptrSrv->pDQ2, pDQ2, usePin, debug, os);
    setInitVals(ptrMPI->ptrSrv->pDQ3, pDQ3, usePin, debug, os);
    setInitVals(ptrMPI->ptrSrv->pDQ4, pDQ4, usePin, debug, os);
    //MSE parameters    
    setInitVals(ptrMPI->ptrMSE->pMSE_LnC, pMSE_LnC, usePin, debug, os);
}

int model_parameters::checkParams(int debug, ostream& os)
{
    int res = 0;
    //recruitment parameters
    res += checkParams(ptrMPI->ptrRec->pLnR, pLnR, debug,os);
    res += checkParams(ptrMPI->ptrRec->pRCV, pRCV, debug,os);
    res += checkParams(ptrMPI->ptrRec->pRX,  pRX,  debug,os);
    res += checkParams(ptrMPI->ptrRec->pRa,  pRa,  debug,os);
    res += checkParams(ptrMPI->ptrRec->pRb,  pRb,  debug,os);
    res += checkParams(ptrMPI->ptrRec->pDevsLnR,pDevsLnR,debug,os);
    //natural mortality parameters
    res += checkParams(ptrMPI->ptrNM->pM, pM, debug,os);
    res += checkParams(ptrMPI->ptrNM->pDM1, pDM1, debug,os);
    res += checkParams(ptrMPI->ptrNM->pDM2, pDM2, debug,os);
    res += checkParams(ptrMPI->ptrNM->pDM3, pDM3, debug,os);
    res += checkParams(ptrMPI->ptrNM->pDM4, pDM4, debug,os);
    //growth parameters
    res += checkParams(ptrMPI->ptrGrw->pGrA,   pGrA,   debug,os);
    res += checkParams(ptrMPI->ptrGrw->pGrB,   pGrB,   debug,os);
    res += checkParams(ptrMPI->ptrGrw->pGrBeta,pGrBeta,debug,os);
    //maturity parameters
    res += checkParams(ptrMPI->ptrM2M->pvLgtPrM2M,pvLgtPrM2M,debug,os);
    //selectivity parameters
    res += checkParams(ptrMPI->ptrSel->pS1, pS1,debug,os);
    res += checkParams(ptrMPI->ptrSel->pS2, pS2,debug,os);
    res += checkParams(ptrMPI->ptrSel->pS3, pS3,debug,os);
    res += checkParams(ptrMPI->ptrSel->pS4, pS4,debug,os);
    res += checkParams(ptrMPI->ptrSel->pS5, pS5,debug,os);
    res += checkParams(ptrMPI->ptrSel->pS6, pS6,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS1, pDevsS1,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS2, pDevsS2,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS3, pDevsS3,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS4, pDevsS4,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS5, pDevsS5,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS6, pDevsS6,debug,os);
    res += checkParams(ptrMPI->ptrSel->pvNPSel, pvNPSel,debug,os);
    //fully-selected fishing capture rate parameters
    res += checkParams(ptrMPI->ptrFsh->pHM,     pHM,     debug,os);
    res += checkParams(ptrMPI->ptrFsh->pLnC,    pLnC,    debug,os);
    res += checkParams(ptrMPI->ptrFsh->pDC1,    pDC1,    debug,os);
    res += checkParams(ptrMPI->ptrFsh->pDC2,    pDC2,    debug,os);
    res += checkParams(ptrMPI->ptrFsh->pDC3,    pDC3,    debug,os);
    res += checkParams(ptrMPI->ptrFsh->pDC4,    pDC4,    debug,os);
    res += checkParams(ptrMPI->ptrFsh->pDevsLnC,pDevsLnC,debug,os);
    res += checkParams(ptrMPI->ptrFsh->pLnEffX, pLnEffX, debug,os);
    res += checkParams(ptrMPI->ptrFsh->pLgtRet, pLgtRet, debug,os);
    //survey catchability parameters
    res += checkParams(ptrMPI->ptrSrv->pQ,   pQ,   debug,os);
    res += checkParams(ptrMPI->ptrSrv->pDQ1, pDQ1, debug,os);
    res += checkParams(ptrMPI->ptrSrv->pDQ2, pDQ2, debug,os);
    res += checkParams(ptrMPI->ptrSrv->pDQ3, pDQ3, debug,os);
    res += checkParams(ptrMPI->ptrSrv->pDQ4, pDQ4, debug,os);
    //MSE parameters
    res += checkParams(ptrMPI->ptrMSE->pMSE_LnC,pMSE_LnC,debug,os);
    return res;
}

int model_parameters::checkParams(NumberVectorInfo* pI, param_init_number_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting checkParams(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    int r = 0;
    if (np){
        for (int i=1;i<=np;i++) {
            if (isnan(value(p(i)))){
                r++;
                if (debug) cout<<"NaN detected in checkParams() for "<<p(1).label()<<"["<<i<<"]"<<endl;
            }
        }
    }
    if (debug>=dbgAll) {
        std::cout<<"Finished checkParams(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    }
    return r;
}

int model_parameters::checkParams(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting checkParams(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    int r = 0;
    if (np){
        for (int i=1;i<=np;i++) {
            if (isnan(value(p(i)))){
                r++;
                if (debug) cout<<"NaN detected in checkParams() for "<<p(1).label()<<"["<<i<<"]"<<endl;
            }
        }
    }
    if (debug>=dbgAll) {
        std::cout<<"Finished checkParams(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    }
    return r;
}

int model_parameters::checkParams(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting checkParams(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    int r = 0;
    if (np){
        for (int i=1;i<=np;i++) {
            for (int j=p(i).indexmin();j<=p(i).indexmax();j++) {
                if (isnan(value(p(i,j)))){
                    r++;
                    if (debug) cout<<"NaN detected in checkParams() for "<<p(1).label()<<"["<<i<<cc<<j<<"]"<<endl;
                }
            }
        }
    }
    if (debug>=dbgAll) {
        std::cout<<"Finished checkParams(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }
    return r;
}

int model_parameters::checkParams(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting checkParams(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    int r = 0;
    if (np){
        for (int i=1;i<=np;i++) {
            for (int j=p(i).indexmin();j<=p(i).indexmax();j++) {
                if (isnan(value(p(i,j)))){
                    r++;
                    if (debug) cout<<"NaN detected in checkParams() for "<<p(1).label()<<"["<<i<<cc<<j<<"]"<<endl;
                }
            }
        }
    }
    if (debug>=dbgAll) {
        std::cout<<"Finished checkParams(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }
    return r;
}

void model_parameters::writeMCMCHeader(void)
{
    mcmc.open((char*)(fnMCMC),ofstream::out|ofstream::trunc);
    ctrMCMC = 0;
    mcmc<<"mcmc<-list();"<<endl;
    mcmc.close();
}

void model_parameters::writeMCMCtoR(ostream& mcmc,NumberVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
}

void model_parameters::writeMCMCtoR(ostream& mcmc,BoundedNumberVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
}

void model_parameters::writeMCMCtoR(ostream& mcmc,BoundedVectorVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
}

void model_parameters::writeMCMCtoR(ostream& mcmc,DevsVectorVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
}

void model_parameters::writeMCMCtoR(ofstream& mcmc)
{
    mcmc.open((char *) fnMCMC, ofstream::out|ofstream::app);
    ctrMCMC+=1;
    std::cout<<"writing mcmc iteration "<<ctrMCMC<<endl;
    mcmc<<"mcmc[["<<ctrMCMC<<"]]<-list(objFun="<<objFun<<cc<<endl;
    //write parameter values
        //recruitment values
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnR);     mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pRCV);     mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pRX);      mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pRa);      mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pRb);      mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pDevsLnR); mcmc<<cc<<endl;
        //natural mortality parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pM);   mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pDM1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pDM2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pDM3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pDM4); mcmc<<cc<<endl;
        //growth parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrGrw->pGrA);    mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrGrw->pGrB);    mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrGrw->pGrBeta); mcmc<<cc<<endl;
        //maturity parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrM2M->pvLgtPrM2M); mcmc<<cc<<endl;
        //selectivity parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS4); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS5); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS6); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS4); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS5); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS6); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pvNPSel); mcmc<<cc<<endl;
        //fully-selected fishing capture rate parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pHM);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnC); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDC1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDC2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDC3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDC4); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDevsLnC); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnEffX);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLgtRet);  mcmc<<cc<<endl;
        //survey catchability parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pQ);    mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pDQ1);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pDQ2);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pDQ3);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pDQ4);  mcmc<<cc<<endl;
        //write other quantities
        mcmc<<"R_y="; wts::writeToR(mcmc,value(R_y)); mcmc<<cc<<endl;
        ivector bnds = wts::getBounds(spB_yx);
        mcmc<<"MB_xy="; wts::writeToR(mcmc,trans(value(spB_yx)),xDms,yDms); 
        if (doOFL){
            mcmc<<cc<<endl;
            calcOFL(mxYr+1,0,cout);//updates oflresults
            ptrOFLResults->writeToR(mcmc,ptrMC,"ptrOFLResults",0);//mcm<<cc<<endl;
        }
    mcmc<<");"<<endl;
    mcmc.close();
}

void model_parameters::createSimData(int debug, ostream& cout, int iSimDataSeed, ModelDatasets* ptrSim)
{
    if (debug)cout<<"simulating model results as data"<<endl;
    d6_array vn_vyxmsz = wts::value(n_vyxmsz);
    d6_array vcN_fyxmsz = wts::value(cpN_fyxmsz);
    d6_array vrmN_fyxmsz = wts::value(rmN_fyxmsz);
    for (int f=1;f<=nFsh;f++) {
        if (debug) cout<<"fishery f: "<<f<<endl;
        (ptrSim->ppFsh[f-1])->replaceFisheryCatchData(iSimDataSeed,rngSimData,vcN_fyxmsz(f),vrmN_fyxmsz(f),ptrSim->ptrBio->wAtZ_xmz);
    }
    for (int v=1;v<=nSrv;v++) {
        if (debug) cout<<"survey "<<v<<endl;
        (ptrSim->ppSrv[v-1])->replaceIndexCatchData(iSimDataSeed,rngSimData,vn_vyxmsz(v),ptrSim->ptrBio->wAtZ_xmz);
    }
    if (debug) cout<<"finished simulating model results as data"<<endl;
}

void model_parameters::writeSimData(ostream& os, int debug, ostream& cout, ModelDatasets* ptrSim)
{
    if (debug)cout<<"writing model results as data"<<endl;
    for (int v=1;v<=nSrv;v++) {
        os<<"#------------------------------------------------------------"<<endl;
        os<<(*(ptrSim->ppSrv[v-1]))<<endl;
    }
    //     cout<<4<<endl;
    for (int f=1;f<=nFsh;f++) {
        os<<"#------------------------------------------------------------"<<endl;
        os<<(*(ptrSim->ppFsh[f-1]))<<endl;
    }
    if (debug) cout<<"finished writing model results as data"<<endl;
}

void model_parameters::setInitVals(NumberVectorInfo* pI, param_init_number_vector& p, int usePin, int debug, ostream& os)
{
    if (debug>=dbgAll) os<<"Starting setInitVals(NumberVectorInfo* pI, param_init_number_vector& p, usePin, debug, os) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        //parameters have been defined
        dvector aovls = 1.0*pI->getInitVals();             //original initial values on arithmetic scales from parameter info
        dvector povls = 1.0*pI->getInitValsOnParamScales();//original initial values on parameter scales from parameter info
        dvector afvls = 1.0*pI->getInitVals();             //final initial values on arithmetic scales from parameter info
        if (usePin) {
            //use values from pinfile assigned to p as initial values
            os<<"Using pin file to set initial values for "<<p(1).label()<<endl;
            pI->setInitValsFromParamVals(p);//update initial values in pI to those from pinfile
            afvls = 1.0*pI->getInitVals();
            os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
            for (int i=1;i<=np;i++) os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
        } else {
            //use values from pI as initial values
            for (int i=1;i<=np;i++) {
                p(i) = povls(i);  //assign original initial value from parameter info
                NumberInfo* ptrI = (*pI)[i];
                if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                    //assign final initial value based on resampling prior pdf
                    os<<"Using resampling to set initial values for "<<p(i).label()<<endl;
                    afvls(i) = ptrI->drawInitVal(rng,ptrMC->vif); //get resampled initial value on arithmetic scale
                    p(i) = ptrI->calcParamScaleVal(afvls(i));     //calc initial parameter value
                    ptrI->setInitVal(afvls(i));                  //update ptrI
                    os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
                    os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
                } else {
                    os<<"Using MPI to set initial values for "<<p(i).label()<<endl;
                    os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
                    os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
                }
            }
        }
    } else {
        //parameters have not been defined
        os<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        os<<"Finished setInitVals(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    }
}

void model_parameters::setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int usePin, int debug, ostream& os)
{
    if (debug>=dbgAll) os<<"Starting setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        //parameters have been defined
        dvector aovls = 1.0*pI->getInitVals();             //original initial values on arithmetic scales from parameter info
        dvector povls = 1.0*pI->getInitValsOnParamScales();//original initial values on parameter scales from parameter info
        dvector afvls = 1.0*pI->getInitVals();             //original initial values on arithmetic scales from parameter info
        if (usePin) {
            //use values from pinfile assigned to p as initial values
            os<<"Using pin file to set initial values for "<<p(1).label()<<endl;
            pI->setInitValsFromParamVals(p);//update initial values in pI to those from pinfile
            afvls = 1.0*pI->getInitVals();
            os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
            for (int i=1;i<=np;i++) os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
        } else {
            //use values from pI as initial values
            for (int i=1;i<=np;i++) {
                p(i) = povls(i);  //assign original initial value from parameter info
                BoundedNumberInfo* ptrI = (*pI)[i];
                if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                    //assign final initial value based on jittering
                    os<<"Using jittering to set initial values for "<<p(i).label()<<endl;
                    afvls(i) = ptrI->jitterInitVal(rng,ptrMC->jitFrac);//get jittered initial value on arithmetic scale
                    p(i) = ptrI->calcParamScaleVal(afvls(i));          //calc initial parameter value
                    ptrI->setInitVal(afvls(i));                       //update ptrI
                    os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
                    os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
                } else 
                if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                    //assign final initial value based on resampling prior pdf
                    os<<"Using resampling to set initial values for "<<p(i).label()<<endl;
                    afvls(i) = ptrI->drawInitVal(rng,ptrMC->vif); //get resampled initial value on arithmetic scale
                    p(i) = ptrI->calcParamScaleVal(afvls(i));     //calc initial parameter value
                    ptrI->setInitVal(afvls(i));                  //update ptrI
                    os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
                    os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
                } else {
                    os<<"Using MPI to set initial values for "<<p(i).label()<<endl;
                    os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
                    os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
                }
            }
        }
    } else {
        os<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    }
}

void model_parameters::setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int usePin, int debug, ostream& os)
{
    //debug=dbgAll;
    if (debug>=dbgAll) os<<"Starting setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        if (usePin){
            //use values from pinfile assigned to p as initial values
            for (int i=1;i<=np;i++) {
                os<<"Using pin file to set initial values for "<<p(i).label()<<endl;
                BoundedVectorInfo* ptrI = (*pI)[i];
                dvector pnvls = value(p(i));                       //original initial values on parameter scale from pin file
                dvector aovls = 1.0*ptrI->getInitVals();            //original initial values on arithmetic scale from parameter info
                dvector povls = 1.0*ptrI->getInitValsOnParamScale();//original initial values on parameter scale from parameter info
                ptrI->setInitValsFromParamVals(p(i));               //set final initial values on arithmetic scale for parameter info
                os<<tb<<"pinfile    inits : "<<pnvls<<endl;
                os<<tb<<"orig arith inits : "<<aovls<<endl;
                os<<tb<<"orig param inits : "<<povls<<endl;
                os<<tb<<"final arith inits: "<<ptrI->getInitVals()<<endl;//final initial values on arithmetic scale from parameter info
                os<<tb<<"final param inits: "<<value(p(i))<<endl;        //final initial values on parameter scale from parameter info
            }
        } else {
            //use values based on pI as initial values for p
            for (int i=1;i<=np;i++) {
                os<<"Using MPI to set initial values for "<<p(i).label()<<endl;
                dvector pnvls = value(p(i));                            //original initial values on parameter scale from pin file
                if (debug>=dbgAll) os<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<pnvls.indexmin()<<tb<<pnvls.indexmax()<<endl;
                os<<tb<<"pinfile    inits : "<<pnvls<<endl;
                BoundedVectorInfo* ptrI = static_cast<DevsVectorInfo*>((*pI)[i]);
                dvector aovls = 1.0*ptrI->getInitVals();            //initial values on arithmetic scale from parameter info
                os<<tb<<"orig arith inits : "<<aovls<<endl;
                dvector povls = 1.0*ptrI->getInitValsOnParamScale();//initial values on parameter scales from parameter info
                os<<tb<<"orig param inits : "<<povls<<endl;
                p(i)=povls;//set initial values on parameter scales from parameter info
                if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                    os<<"Using jittering to set initial values for "<<p(i).label()<<" : "<<endl;
                    if (debug>=dbgAll) os<<"Using jittering to set initial values for "<<p(i).label()<<" : "<<endl;
                    dvector afvls = ptrI->jitterInitVals(rng,ptrMC->jitFrac);//get jittered values on arithmetic scale
                    ptrI->setInitVals(afvls);                               //set jittered values as initial values on arithmetic scale
                    p(i)=ptrI->getInitValsOnParamScale();                   //get jittered values on param scale as initial parameter values
                } else
                if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                    os<<"Using resampling to set initial values for "<<p(i).label()<<" : "<<endl;
                    if (debug>=dbgAll) os<<"Using resampling to set initial values for "<<p(i).label()<<" : "<<endl;
                    dvector afvls = ptrI->drawInitVals(rng,ptrMC->vif);//get resampled values on arithmetic scale
                    ptrI->setInitVals(afvls);                          //set resampled values as initial values on arithmetic scale
                    p(i)=ptrI->getInitValsOnParamScale();              //get resampled values on param scale as initial parameter values
                } else {
                    os<<"Using MPI to set initial values for "<<p(i).label()<<endl;
                }
                os<<tb<<"final arith inits: "<<ptrI->getInitVals()<<endl;//final initial values on arithmetic scale from parameter info
                os<<tb<<"final param inits: "<<value(p(i))<<endl;        //final initial values on parameter scale from parameter info
            }
        }
    } else {
        os<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        os<<"Finished setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }
}

void model_parameters::setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int usePin, int debug, ostream& os)
{
    //debug=dbgAll;
    if (debug>=dbgAll) os<<"Starting setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        if (usePin){
            //use values from pinfile assigned to p as initial values
            for (int i=1;i<=np;i++) {
                os<<"Using pin file to set initial values for "<<p(i).label()<<endl;
                DevsVectorInfo* ptrI = static_cast<DevsVectorInfo*>((*pI)[i]);
                dvector pnvls = value(p(i));                       //original initial values on parameter scale from pin file
                dvector aovls = 1.0*ptrI->getInitVals();            //original initial values on arithmetic scale from parameter info
                dvector povls = 1.0*ptrI->getInitValsOnParamScale();//original initial values on parameter scale from parameter info
                ptrI->setInitValsFromParamVals(p(i));               //set final initial values on arithmetic scale for parameter info
                os<<tb<<"pinfile    inits : "<<pnvls<<endl;
                os<<tb<<"orig arith inits : "<<aovls<<endl;
                os<<tb<<"orig param inits : "<<povls<<endl;
                os<<tb<<"final arith inits: "<<ptrI->getInitVals()<<endl;//final initial values on arithmetic scale from parameter info
                os<<tb<<"final param inits: "<<value(p(i))<<endl;        //final initial values on parameter scale from parameter info
            }
        } else {
            //use values based on pI as initial values for p
            for (int i=1;i<=np;i++) {
                dvector pnvls = value(p(i));                            //original initial values on parameter scale from pin file
                dvector aovls = 1.0*(*pI)[i]->getInitVals();            //initial values on arithmetic scale from parameter info
                dvector povls = 1.0*(*pI)[i]->getInitValsOnParamScale();//initial values on parameter scales from parameter info
                if (debug>=dbgAll) os<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<pnvls.indexmin()<<tb<<pnvls.indexmax()<<endl;
                p(i)=povls;//set initial values on parameter scales from parameter info
                DevsVectorInfo* ptrI = (*pI)[i];
                if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                    os<<"Using jittering to set initial values for "<<p(i).label()<<" : "<<endl;
                    if (debug>=dbgAll) os<<"Using jittering to set initial values for "<<p(i).label()<<" : "<<endl;
                    dvector afvls = ptrI->jitterInitVals(rng,ptrMC->jitFrac);//get jittered values on arithmetic scale
                    ptrI->setInitVals(afvls);                               //set jittered values as initial values on arithmetic scale
                    p(i)=ptrI->getInitValsOnParamScale();                   //get jittered values on param scale as initial parameter values
                } else
                if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                    os<<"Using resampling to set initial values for "<<p(i).label()<<" : "<<endl;
                    if (debug>=dbgAll) os<<"Using resampling to set initial values for "<<p(i).label()<<" : "<<endl;
                    dvector afvls = ptrI->drawInitVals(rng,ptrMC->vif);//get resampled values on arithmetic scale
                    ptrI->setInitVals(afvls);                          //set resampled values as initial values on arithmetic scale
                    p(i)=ptrI->getInitValsOnParamScale();              //get resampled values on param scale as initial parameter values
                } else {
                    os<<"Using MPI to set initial values for "<<p(i).label()<<endl;
                }
                os<<tb<<"pinfile    inits : "<<pnvls<<endl;
                os<<tb<<"orig arith inits : "<<aovls<<endl;
                os<<tb<<"orig param inits : "<<povls<<endl;
                os<<tb<<"final arith inits: "<<ptrI->getInitVals()<<endl;//final initial values on arithmetic scale from parameter info
                os<<tb<<"final param inits: "<<value(p(i))<<endl;        //final initial values on parameter scale from parameter info
            }
        }
    } else {
        os<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        os<<"Finished setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }
}

void model_parameters::setAllDevs(int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"starting setAllDevs()"<<endl;
    if (debug>=dbgAll) cout<<"setDevs() for pLnR"<<endl;
    tcsam::setDevs(devsLnR, pDevsLnR, ptrMPI->ptrRec->pDevsLnR,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS1"<<endl;
    tcsam::setDevs(devsS1, pDevsS1, ptrMPI->ptrSel->pDevsS1, debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS2"<<endl;
    tcsam::setDevs(devsS2, pDevsS2, ptrMPI->ptrSel->pDevsS2,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS3"<<endl;
    tcsam::setDevs(devsS3, pDevsS3, ptrMPI->ptrSel->pDevsS3,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS4"<<endl;
    tcsam::setDevs(devsS4, pDevsS4, ptrMPI->ptrSel->pDevsS4,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS5"<<endl;
    tcsam::setDevs(devsS5, pDevsS5, ptrMPI->ptrSel->pDevsS5,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS6"<<endl;
    tcsam::setDevs(devsS6, pDevsS6, ptrMPI->ptrSel->pDevsS6,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsLnC"<<endl;
    tcsam::setDevs(devsLnC, pDevsLnC, ptrMPI->ptrFsh->pDevsLnC,debug,cout);
    if (debug>=dbgAll) cout<<"finished setAllDevs()"<<endl;
}

void model_parameters::ReportToR_Data(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_Data(...)"<<endl;
    ptrMDS->writeToR(os,"data",0);
    if (debug) cout<<"Finished ReportToR_Data(...)"<<endl;
}

void model_parameters::ReportToR_Params(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_Params(...)"<<endl;
    ptrMPI->writeToR(os);
    if (debug) cout<<"Finished ReportToR_Params(...)"<<endl;
}

void model_parameters::ReportToR_ModelProcesses(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_ModelProcesses(...)"<<endl;
    os<<"mp=list("<<endl;
        os<<"M_c      ="; wts::writeToR(os,value(M_c),     adstring("pc=1:"+str(npcNM )));      os<<cc<<endl;
        os<<"prM2M_cz ="; wts::writeToR(os,value(prM2M_cz),adstring("pc=1:"+str(npcM2M)),zbDms);os<<cc<<endl;
        os<<"sel_cz   ="; wts::writeToR(os,value(sel_cz),  adstring("pc=1:"+str(npcSel)),zbDms);os<<cc<<endl;
        os<<"M_yxmsz  ="; wts::writeToR(os,wts::value(M_yxmsz),   yDms,xDms,mDms,sDms,zbDms);   os<<cc<<endl;
        os<<"prM2M_yxz =";wts::writeToR(os,     value(prM2M_yxz), yDms,xDms,zbDms);             os<<cc<<endl;
        os<<"T_list=list("<<endl;
            os<<"mnZAM_cz   ="; wts::writeToR(os,value(mnGrZ_cz),adstring("pc=1:"+str(npcGrw )),zbDms);       os<<cc<<endl;
            os<<"T_czz      ="; wts::writeToR(os,value(prGr_czz),adstring("pc=1:"+str(npcGrw )),zbDms,zpDms); os<<cc<<endl;
            os<<"mnZAM_yxsz =";  wts::writeToR(os,wts::value(mnGrZ_yxsz),yDms,xDms,sDms,zbDms);              os<<cc<<endl;
            if (0){//TDODO: develop option to output
                os<<"T_yxszz    =";  wts::writeToR(os,wts::value(prGr_yxszz),yDms,xDms,sDms,zbDms,zpDms);    os<<endl;
            } else {
                os<<"T_yxszz    =NULL"<<endl;
            }
        os<<")"<<cc<<endl;
        os<<"R_list=list("<<endl;
            os<<"R_y  ="; wts::writeToR(os,value(R_y), yDms);                                os<<cc<<endl;
            os<<"R_yx ="; wts::writeToR(os,value(R_yx),yDms,xDms);                           os<<cc<<endl;
            os<<"R_yz ="; wts::writeToR(os,value(R_yz),yDms,zbDms);                          os<<cc<<endl;
            os<<"Rx_c ="; wts::writeToR(os,value(Rx_c),adstring("pc=1:"+str(npcRec)));       os<<cc<<endl;
            os<<"R_cz ="; wts::writeToR(os,value(R_cz),adstring("pc=1:"+str(npcRec)),zbDms); os<<endl;
        os<<")"<<cc<<endl;
        if (nEASs){
        os<<"EffX_list=list("<<endl;
            os<<"effAvgScenarios="; ptrMOs->ptrEffXtrapScenarios->ptrEffAvgScenarios->writeToR(os); os<<cc<<endl;
            os<<"capRateScenarios="; ptrMOs->ptrEffXtrapScenarios->ptrCapRateAvgScenarios->writeToR(os); os<<cc<<endl;
            os<<"avgEff=list("<<endl;
            for (int n=eff_ny.indexmin();n<=eff_ny.indexmax();n++){
                os<<"`"<<n<<"`=list(avgEff="<<avgEff_n(n)<<cc;
                os<<"eff="; wts::writeToR(os,eff_ny(n),ptrMOs->ptrEffXtrapScenarios->ptrEffAvgScenarios->ppEASs[n-1]->getYDimsForR()); os<<"), "<<endl;
            }
            os<<"NULL)"<<endl;
        os<<")"<<cc<<endl;
        } else {
            os<<"EffX_list=NULL"<<cc<<endl;
        }
        d6_array tmF_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);
        for (int f=1;f<=nFsh;f++){
            for (int y=mnYr;y<=mxYr;y++){
                for (int x=1;x<=nSXs;x++){
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++){
                            tmF_fyxmsz(f,y,x,m,s) = value(rmF_fyxmsz(f,y,x,m,s)+dmF_fyxmsz(f,y,x,m,s));
                        }
                    }
                }
            }
        }
        os<<"F_list=list("<<endl;
            os<<"hm_fy     ="; wts::writeToR(os,     value(hmF_fy),    fDms,yDms);                     os<<cc<<endl;
            os<<"cpF_fyxms ="; wts::writeToR(os,wts::value(cpF_fyxms ),fDms,yDms,xDms,mDms,sDms);      os<<cc<<endl;
            os<<"sel_fyxmsz="; wts::writeToR(os,wts::value(sel_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"ret_fyxmsz="; wts::writeToR(os,wts::value(ret_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"tmF_yxmsz ="; wts::writeToR(os,wts::value(tmF_yxmsz),      yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"cpF_fyxmsz="; wts::writeToR(os,wts::value(cpF_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"rmF_fyxmsz="; wts::writeToR(os,wts::value(rmF_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dmF_fyxmsz="; wts::writeToR(os,wts::value(dmF_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"tmF_fyxmsz="; wts::writeToR(os,            tmF_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")"<<cc<<endl;
        os<<"S_list=list("<<endl;
            os<<"sel_vyxmsz="; wts::writeToR(os,wts::value(s_vyxmsz),vDms,ypDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"Q_vyxms   ="; wts::writeToR(os,wts::value(q_vyxms), vDms,ypDms,xDms,mDms,sDms);       os<<cc<<endl;
            os<<"Q_vyxmsz  ="; wts::writeToR(os,wts::value(q_vyxmsz),vDms,ypDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")";
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelProcesses(...)"<<endl;
}

void model_parameters::ReportToR_ModelResults(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_ModelResults(...)"<<endl;
    d3_array ones(1,nSXs,1,nMSs,1,nZBs);//weighting for simple sums
    for (int x=1;x<=nSXs;x++) ones(x) = 1.0;
    //population numbers, biomass
    d5_array vn_yxmsz = wts::value(n_yxmsz);
    d5_array vb_yxmsz = tcsam::calcBiomass(vn_yxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    //population male maturity ogives
    dmatrix vPM_yz(mnYr,mxYrp1,1,nZBs);
    for (int y=mnYr;y<=mxYrp1;y++){
        vPM_yz(y) = elem_div(vn_yxmsz(y,MALE,MATURE,NEW_SHELL),
                            (1.0e-5)+vn_yxmsz(y,MALE,IMMATURE,NEW_SHELL)+vn_yxmsz(y,MALE,MATURE,NEW_SHELL));
    }
    //numbers, biomass captured (NOT mortality)
    d6_array vcpN_fyxmsz = wts::value(cpN_fyxmsz);
    d6_array vcpB_fyxmsz = tcsam::calcBiomass(vcpN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    //numbers, biomass discards (NOT mortality)
    d6_array vdsN_fyxmsz = wts::value(dsN_fyxmsz);
    d6_array vdsB_fyxmsz = tcsam::calcBiomass(vdsN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    //numbers, biomass retained (mortality)
    d6_array vrmN_fyxmsz = wts::value(rmN_fyxmsz);
    d6_array vrmB_fyxmsz = tcsam::calcBiomass(vrmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    //numbers, biomass discard mortality
    d6_array vdmN_fyxmsz = wts::value(dmN_fyxmsz);
    d6_array vdmB_fyxmsz = tcsam::calcBiomass(vdmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    //survey numbers, biomass
    d6_array vn_vyxmsz = wts::value(n_vyxmsz);
    d6_array vb_vyxmsz = tcsam::calcBiomass(vn_vyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    d3_array vPM_vyz(1,nSrv,mnYr,mxYrp1,1,nZBs);
    for (int v=1;v<=nSrv;v++){
        for (int y=mnYr;y<=mxYrp1;y++){
            vPM_vyz(v,y) = elem_div(vn_vyxmsz(v,y,MALE,MATURE,NEW_SHELL),
                                    (1.0e-5)+vn_vyxmsz(v,y,MALE,IMMATURE,NEW_SHELL)+vn_vyxmsz(v,y,MALE,MATURE,NEW_SHELL));
        }
    }
    os<<"mr=list("<<endl;
        os<<"iN_xmsz ="; wts::writeToR(os,vn_yxmsz(mnYr),xDms,mDms,sDms,zbDms); os<<cc<<endl;
        os<<"P_list=list("<<endl;
            os<<"MB_yx    ="; wts::writeToR(os,value(spB_yx), yDms,xDms);                        os<<cc<<endl;
            os<<"N_yxmsz  ="; wts::writeToR(os,vn_yxmsz,ypDms,xDms,mDms,sDms,zbDms);             os<<cc<<endl;
            os<<"B_yxmsz  ="; wts::writeToR(os,vb_yxmsz,ypDms,xDms,mDms,sDms,zbDms);             os<<cc<<endl;
            os<<"prM_yz   ="; wts::writeToR(os,vPM_yz,ypDms,zbDms);                              os<<cc<<endl;
            os<<"nmN_yxmsz="; wts::writeToR(os,wts::value(nmN_yxmsz),yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"tmN_yxmsz="; wts::writeToR(os,wts::value(tmN_yxmsz),yDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")"<<cc<<endl;    
        os<<"F_list=list("<<endl;
            os<<"cpN_fyxmsz="; wts::writeToR(os,vcpN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dsN_fyxmsz="; wts::writeToR(os,vdsN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"rmN_fyxmsz="; wts::writeToR(os,vrmN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dmN_fyxmsz="; wts::writeToR(os,vdmN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"cpB_fyxmsz="; wts::writeToR(os,vcpB_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dsB_fyxmsz="; wts::writeToR(os,vdsB_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"rmB_fyxmsz="; wts::writeToR(os,vrmB_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dmB_fyxmsz="; wts::writeToR(os,vdmB_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")"<<cc<<endl;
        os<<"S_list=list("<<endl;
           os<<"MB_vyx  ="; wts::writeToR(os,value(mb_vyx),vDms,ypDms,xDms);                 os<<cc<<endl;
           os<<"N_vyxmsz="; wts::writeToR(os,    vn_vyxmsz,vDms,ypDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
           os<<"B_vyxmsz="; wts::writeToR(os,    vb_vyxmsz,vDms,ypDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
           os<<"prM_vyz ="; wts::writeToR(os,    vPM_vyz,vDms,ypDms,zbDms);                  os<<endl;
       os<<")";
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelResults(...)"<<endl;
}

void model_parameters::ReportToR_CohortProgression(ostream& os,  int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_CohortProgression(...)"<<endl;
    d5_array n_yxmsz = calcCohortProgression(mxYr,1,cout);
    int ny = n_yxmsz.indexmax();
    os<<"cohortprogression=list("<<endl;
        os<<"n_yxmsz="; wts::writeToR(os,n_yxmsz,adstring("y=0:"+str(ny)),xDms,mDms,sDms,zbDms); os<<endl;
    os<<")";
    if (debug) cout<<"Finished ReportToR_CohortProgression(...)"<<endl;
}

void model_parameters::ReportToR_ModelFits(ostream& os, double maxGrad, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_ModelFits(...)"<<endl;
    os<<"model.fits=list("<<endl;
        os<<tb<<"objfun="<<value(objFun)<<cc<<"maxGrad="<<maxGrad<<cc<<endl;
        //recalc objective function components and and write results to os
        objFun.initialize();
        os<<tb<<"penalties="; calcPenalties(-1,os);      os<<cc<<endl;
        os<<tb<<"#end of penalties"<<endl;
        os<<tb<<"priors=";    calcAllPriors(-1,os);      os<<cc<<endl;
        os<<tb<<"#end of priors"<<endl;
        os<<tb<<"components=list("<<endl;
            os<<tb<<tb<<"recruitment="; calcNLLs_Recruitment(-1,os); os<<endl;
        os<<tb<<")"<<cc<<endl;
        os<<tb<<"#end of components"<<endl;
        os<<tb<<"fisheries=";  calcNLLs_Fisheries(-1,os);          os<<cc<<endl; 
        os<<tb<<"#end of fisheries"<<endl;
        os<<tb<<"surveys=";    calcNLLs_Surveys(-1,os);            os<<cc<<endl;  
        os<<tb<<"#end of surveys"<<endl;
        os<<tb<<"growthdata="; calcNLLs_GrowthData(-1,os);         os<<cc<<endl;
        os<<tb<<"#end of growthdata"<<endl;
        os<<tb<<"maturitydata="; calcNLLs_ChelaHeightData(-1,os);  os<<cc<<endl;
        os<<tb<<"#end of maturitydata"<<endl;
        os<<tb<<"effortdata="; calcNLLs_ExtrapolatedEffort(-1,os); os<<endl;
        os<<tb<<"#end of effortdata"<<endl;
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelFits(...)"<<endl;
}

void model_parameters::updateMPI(int debug, ostream& cout)
{
    if (debug) cout<<"Starting updateMPI(...)"<<endl;
    //recruitment parameters
    if (debug) cout<<"starting recruitment parameters"<<endl;
    ptrMPI->ptrRec->pLnR->setFinalValsFromParamVals(pLnR);
    ptrMPI->ptrRec->pRCV->setFinalValsFromParamVals(pRCV);
    ptrMPI->ptrRec->pRX->setFinalValsFromParamVals(pRX);
    ptrMPI->ptrRec->pRa->setFinalValsFromParamVals(pRa);
    ptrMPI->ptrRec->pRb->setFinalValsFromParamVals(pRb);
    if (debug) cout<<"setting final vals for pDevsLnR"<<endl;
    for (int p=1;p<=ptrMPI->ptrRec->pDevsLnR->getSize();p++) 
        (*ptrMPI->ptrRec->pDevsLnR)[p]->setFinalValsFromParamVals(pDevsLnR(p));
    if (debug) cout<<"finished recruitment parameters"<<endl;
    //natural mortality parameters
    if (debug) cout<<"starting natural mortality parameters"<<endl;
    ptrMPI->ptrNM->pM->setFinalValsFromParamVals(pM);
    ptrMPI->ptrNM->pDM1->setFinalValsFromParamVals(pDM1);
    ptrMPI->ptrNM->pDM2->setFinalValsFromParamVals(pDM2);
    ptrMPI->ptrNM->pDM3->setFinalValsFromParamVals(pDM3);
    ptrMPI->ptrNM->pDM4->setFinalValsFromParamVals(pDM4);
    if (debug) cout<<"finished natural mortality parameters"<<endl;
    //growth parameters
    if (debug) cout<<"starting growth parameters"<<endl;
    ptrMPI->ptrGrw->pGrA->setFinalValsFromParamVals(pGrA);
    ptrMPI->ptrGrw->pGrB->setFinalValsFromParamVals(pGrB);
    ptrMPI->ptrGrw->pGrBeta->setFinalValsFromParamVals(pGrBeta);
    if (debug) cout<<"finished growth parameters"<<endl;
    //maturity parameters
    //cout<<"setting final vals for pvLgtPrM2M"<<endl;
    if (debug) cout<<"starting prM2M parameters"<<endl;
    for (int p=1;p<=ptrMPI->ptrRec->pDevsLnR->getSize();p++) 
        (*ptrMPI->ptrM2M->pvLgtPrM2M)[p]->setFinalValsFromParamVals(pvLgtPrM2M(p));
    if (debug) cout<<"finished prM2M parameters"<<endl;
    //selectivity parameters
    if (debug) cout<<"starting selectivities"<<endl;
    ptrMPI->ptrSel->pS1->setFinalValsFromParamVals(pS1);
    ptrMPI->ptrSel->pS2->setFinalValsFromParamVals(pS2);
    ptrMPI->ptrSel->pS3->setFinalValsFromParamVals(pS3);
    ptrMPI->ptrSel->pS4->setFinalValsFromParamVals(pS4);
    ptrMPI->ptrSel->pS5->setFinalValsFromParamVals(pS5);
    ptrMPI->ptrSel->pS6->setFinalValsFromParamVals(pS6);
    if (debug) cout<<"setting final vals for pDevsS1"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS1->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS1)[p]->setFinalValsFromParamVals(pDevsS1(p));
    if (debug) cout<<"setting final vals for pDevsS2"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS2->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS2)[p]->setFinalValsFromParamVals(pDevsS2(p));
    if (debug) cout<<"setting final vals for pDevsS3"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS3->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS3)[p]->setFinalValsFromParamVals(pDevsS3(p));
    if (debug) cout<<"setting final vals for pDevsS4"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS4->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS4)[p]->setFinalValsFromParamVals(pDevsS4(p));
    if (debug) cout<<"setting final vals for pDevsS5"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS5->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS5)[p]->setFinalValsFromParamVals(pDevsS5(p));
    if (debug) cout<<"setting final vals for pDevsS6"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS6->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS6)[p]->setFinalValsFromParamVals(pDevsS6(p));
    if (debug) cout<<"setting final vals for pvNPSel"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pvNPSel->getSize();p++) 
        (*ptrMPI->ptrSel->pvNPSel)[p]->setFinalValsFromParamVals(pvNPSel(p));
    if (debug) cout<<"finished selectivities"<<endl;
    //fully-selected fishing capture rate parameters
    if (debug) cout<<"starting fisheries parameters"<<endl;
    ptrMPI->ptrFsh->pHM->setFinalValsFromParamVals(pHM);
    ptrMPI->ptrFsh->pLnC->setFinalValsFromParamVals(pLnC);
    ptrMPI->ptrFsh->pDC1->setFinalValsFromParamVals(pDC1);
    ptrMPI->ptrFsh->pDC2->setFinalValsFromParamVals(pDC2);
    ptrMPI->ptrFsh->pDC3->setFinalValsFromParamVals(pDC3);
    ptrMPI->ptrFsh->pDC4->setFinalValsFromParamVals(pDC4);
    if (debug) cout<<"setting final vals for pDevsLnC"<<endl;
    for (int p=1;p<=ptrMPI->ptrFsh->pDevsLnC->getSize();p++) 
        (*ptrMPI->ptrFsh->pDevsLnC)[p]->setFinalValsFromParamVals(pDevsLnC(p));
    ptrMPI->ptrFsh->pLnEffX->setFinalValsFromParamVals(pLnEffX);
    ptrMPI->ptrFsh->pLgtRet->setFinalValsFromParamVals(pLgtRet);
    if (debug) cout<<"finished fisheries parameters"<<endl;
    //survey catchability parameters
    if (debug) cout<<"starting surveys parameters"<<endl;
    ptrMPI->ptrSrv->pQ->setFinalValsFromParamVals(pQ);
    ptrMPI->ptrSrv->pDQ1->setFinalValsFromParamVals(pDQ1);
    ptrMPI->ptrSrv->pDQ2->setFinalValsFromParamVals(pDQ2);
    ptrMPI->ptrSrv->pDQ3->setFinalValsFromParamVals(pDQ3);
    ptrMPI->ptrSrv->pDQ4->setFinalValsFromParamVals(pDQ4);
    if (debug) cout<<"finished surveys parameters"<<endl;
    //MSE-related parameters
    ptrMPI->ptrMSE->pMSE_LnC->setFinalValsFromParamVals(pMSE_LnC);
    if (debug) cout<<"finished MSE parameters"<<endl;
    if (debug) cout<<"Finished updateMPI(...)"<<endl;
}

void model_parameters::ReportToR(ostream& os, double maxGrad, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR(...)"<<endl;
    updateMPI(debug,cout);
    os<<"res=list("<<endl;
        if (jitter) os<<"jitter="<<iSeed<<cc<<endl;
        os<<"objFun="<<value(objFun)<<cc<<"maxGrad="<<maxGrad<<cc<<endl;
        //model configuration
        ptrMC->writeToR(os,"mc",0); os<<","<<endl;
        os<<"#end of mc"<<endl;
        //model data
        ptrMDS->writeToR(os,"data",0); os<<","<<endl;
        os<<"#end of data"<<endl;
        //parameter values
        ReportToR_Params(os,debug,cout); os<<","<<endl;
        os<<"#end of params"<<endl;
        //model processes
        ReportToR_ModelProcesses(os,debug,cout); os<<","<<endl;
        os<<"#end of modelprocesses"<<endl;
        //model results
        ReportToR_ModelResults(os,debug,cout); os<<","<<endl;
        os<<"#end of modelresults"<<endl;
        //model fit quantities
        ReportToR_ModelFits(os,maxGrad,debug,cout); os<<","<<endl;
        os<<"#end of modelfits"<<endl;
        //simulated model data
        createSimData(debug, cout, 0, ptrSimMDS);//deterministic
        ptrSimMDS->writeToR(os,"sim.data",0); os<<","<<endl;
        os<<"#end of sim.data"<<endl;
        //cohort projections
        ReportToR_CohortProgression(os,debug,cout);
        os<<"#end of cohortprogression"<<endl;
        //do OFL calculations
        if (doOFL&&last_phase()){
            cout<<"ReportToR: starting OFL calculations"<<endl;
            ofstream echoOFL; echoOFL.open("calcOFL.final.txt", ios::trunc);
            echoOFL.precision(12);
            calcOFL(mxYr+1,1,echoOFL);//updates ptrOFLResults
            ptrOFLResults->writeCSVHeader(echoOFL); echoOFL<<endl;
            ptrOFLResults->writeToCSV(echoOFL);     echoOFL<<endl;
            echoOFL.close();
            os<<","<<endl;
            ptrOFLResults->writeToR(os,ptrMC,"ptrOFLResults",0);
            os<<"#end of ptrOFLResults"<<endl;
            cout<<"ReportToR: finished OFL calculations"<<endl;
        }
        //do dynamic B0 calculations
        if (doDynB0&&last_phase()){
            cout<<"ReportToR: starting dynamic B0 calculations"<<endl;
            os<<","<<endl;
            calcDynB0(-1,os);
            os<<"#end of dynamic B0 results"<<endl;
            cout<<"ReportToR: finished dynamic B0 calculations"<<endl;
        }
    os<<")"<<endl;
    if (debug) cout<<"Finished ReportToR(...)"<<endl;
}

void model_parameters::writeParameters(ostream& os,int toR, int willBeActive)
{
    adstring ctg1, ctg2;
    if (!toR) tcsam::writeCSVHeaderForParametersFile(os);
    //recruitment parameters
    ctg1="population processes";
    ctg2="recruitment";
    tcsam::writeParameters(os,pLnR,ctg1,ctg2,ptrMPI->ptrRec->pLnR,toR,willBeActive);      
    tcsam::writeParameters(os,pRCV,ctg1,ctg2,ptrMPI->ptrRec->pRCV,toR,willBeActive);      
    tcsam::writeParameters(os,pRX, ctg1,ctg2,ptrMPI->ptrRec->pRX, toR,willBeActive);      
    tcsam::writeParameters(os,pRa, ctg1,ctg2,ptrMPI->ptrRec->pRa, toR,willBeActive);      
    tcsam::writeParameters(os,pRb, ctg1,ctg2,ptrMPI->ptrRec->pRb, toR,willBeActive);      
    tcsam::writeParameters(os,pDevsLnR,ctg1,ctg2,ptrMPI->ptrRec->pDevsLnR,toR,willBeActive);      
    //natural mortality parameters
    ctg1="population processes";
    ctg2="natural mortality";
    tcsam::writeParameters(os,pM,ctg1,ctg2,ptrMPI->ptrNM->pM,toR,willBeActive);      
    tcsam::writeParameters(os,pDM1,ctg1,ctg2,ptrMPI->ptrNM->pDM1,toR,willBeActive);      
    tcsam::writeParameters(os,pDM2,ctg1,ctg2,ptrMPI->ptrNM->pDM2,toR,willBeActive);      
    tcsam::writeParameters(os,pDM3,ctg1,ctg2,ptrMPI->ptrNM->pDM3,toR,willBeActive);      
    tcsam::writeParameters(os,pDM4,ctg1,ctg2,ptrMPI->ptrNM->pDM4,toR,willBeActive);      
    //growth parameters
    ctg1="population processes";
    ctg2="growth";
    tcsam::writeParameters(os,pGrA,ctg1,ctg2,ptrMPI->ptrGrw->pGrA,toR,willBeActive);      
    tcsam::writeParameters(os,pGrB,ctg1,ctg2,ptrMPI->ptrGrw->pGrB,toR,willBeActive);      
    tcsam::writeParameters(os,pGrBeta,ctg1,ctg2,ptrMPI->ptrGrw->pGrBeta,toR,willBeActive);      
    //maturity parameters
    ctg1="population processes";
    ctg2="maturity";
    tcsam::writeParameters(os,pvLgtPrM2M,ctg1,ctg2,ptrMPI->ptrM2M->pvLgtPrM2M,toR,willBeActive);      
    //selectivity parameters
    ctg1="selectivity";
    ctg2="selectivity";
    tcsam::writeParameters(os,pS1,ctg1,ctg2,ptrMPI->ptrSel->pS1,toR,willBeActive);      
    tcsam::writeParameters(os,pS2,ctg1,ctg2,ptrMPI->ptrSel->pS2,toR,willBeActive);      
    tcsam::writeParameters(os,pS3,ctg1,ctg2,ptrMPI->ptrSel->pS3,toR,willBeActive);      
    tcsam::writeParameters(os,pS4,ctg1,ctg2,ptrMPI->ptrSel->pS4,toR,willBeActive);      
    tcsam::writeParameters(os,pS5,ctg1,ctg2,ptrMPI->ptrSel->pS5,toR,willBeActive);      
    tcsam::writeParameters(os,pS6,ctg1,ctg2,ptrMPI->ptrSel->pS6,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS1,ctg1,ctg2,ptrMPI->ptrSel->pDevsS1,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS2,ctg1,ctg2,ptrMPI->ptrSel->pDevsS2,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS3,ctg1,ctg2,ptrMPI->ptrSel->pDevsS3,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS4,ctg1,ctg2,ptrMPI->ptrSel->pDevsS4,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS5,ctg1,ctg2,ptrMPI->ptrSel->pDevsS5,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS6,ctg1,ctg2,ptrMPI->ptrSel->pDevsS6,toR,willBeActive);      
    tcsam::writeParameters(os,pvNPSel,ctg1,ctg2,ptrMPI->ptrSel->pvNPSel,toR,willBeActive);      
    //fishery parameters
    ctg1="fisheries";
    ctg2="fisheries";
    tcsam::writeParameters(os,pHM, ctg1,ctg2,ptrMPI->ptrFsh->pHM, toR,willBeActive);      
    tcsam::writeParameters(os,pLnC,ctg1,ctg2,ptrMPI->ptrFsh->pLnC,toR,willBeActive);      
    tcsam::writeParameters(os,pDC1,ctg1,ctg2,ptrMPI->ptrFsh->pDC1,toR,willBeActive);      
    tcsam::writeParameters(os,pDC2,ctg1,ctg2,ptrMPI->ptrFsh->pDC2,toR,willBeActive);      
    tcsam::writeParameters(os,pDC3,ctg1,ctg2,ptrMPI->ptrFsh->pDC3,toR,willBeActive);      
    tcsam::writeParameters(os,pDC4,ctg1,ctg2,ptrMPI->ptrFsh->pDC4,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsLnC,ctg1,ctg2,ptrMPI->ptrFsh->pDevsLnC,toR,willBeActive);      
    tcsam::writeParameters(os,pLnEffX, ctg1,ctg2,ptrMPI->ptrFsh->pLnEffX, toR,willBeActive);      
    tcsam::writeParameters(os,pLgtRet, ctg1,ctg2,ptrMPI->ptrFsh->pLgtRet, toR,willBeActive);      
    //survey parameters
    ctg1="surveys";
    ctg2="surveys";
    tcsam::writeParameters(os,pQ,ctg1,ctg2,ptrMPI->ptrSrv->pQ,toR,willBeActive);      
    tcsam::writeParameters(os,pDQ1,ctg1,ctg2,ptrMPI->ptrSrv->pDQ1,toR,willBeActive);      
    tcsam::writeParameters(os,pDQ2,ctg1,ctg2,ptrMPI->ptrSrv->pDQ2,toR,willBeActive);      
    tcsam::writeParameters(os,pDQ3,ctg1,ctg2,ptrMPI->ptrSrv->pDQ3,toR,willBeActive);      
    tcsam::writeParameters(os,pDQ4,ctg1,ctg2,ptrMPI->ptrSrv->pDQ4,toR,willBeActive);
    //MSE directed fishery
    if (mseOpModMode){
        ctg1="MSE-related";
        ctg2="ln-scale directed fishery capture rate";
        tcsam::writeParameters(os,pMSE_LnC,ctg1,ctg2,ptrMPI->ptrMSE->pMSE_LnC,toR,willBeActive);
    }
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
    PRINT2B1(" ")
    PRINT2B1("--REPORT_SECTION---------")
    //max gradient
    double maxGrad = max(fabs(gradients));
    //write active parameters to rpt::echo
    PRINT2B2("Finished phase ",current_phase())
    {
        //write final gradients to file
        ofstream os0("tcsam02.Gradients."+itoa(current_phase(),10)+".csv", ios::trunc);
        for (int i=gradients.indexmin();i<=gradients.indexmax();i++){
            os0<<i<<cc<<gradients[i]<<endl;
        }
        os0.close();
    }
    if (!mseOpModMode) {
        //write objective function components only
        ofstream os0("tcsam02.ModelFits."+itoa(current_phase(),10)+".R", ios::trunc);
        os0.precision(12);
        ReportToR_ModelFits(os0,maxGrad,0,cout);
        os0.close();
    }
    if (last_phase()) {
        PRINT2B1("--ReportToR: last phase ")
        PRINT2B2("----last phase objFun =",value(objFun))
        if (!mseOpModMode){
            //write report as R file
            report.precision(12);
            ReportToR(report,maxGrad,1,rpt::echo);
            //write parameter values to csv
            ofstream os1("tcsam02.params.all.final.csv", ios::trunc);
            os1.precision(12);
            writeParameters(os1,0,0);
            os1.close();
            //write parameter values to csv
            ofstream os2("tcsam02.params.active.final.csv", ios::trunc);
            os2.precision(12);
            writeParameters(os2,0,1);
            os2.close();
            if (option_match(ad_comm::argc,ad_comm::argv,"-jitter")>-1) {
                ofstream fs("jitterInfo.csv");
                fs.precision(20);
                fs<<"seed"<<cc<<"objfun"<<cc<<"maxGrad"<<cc<<"MMB";
                if (doOFL) fs<<cc<<"B0"<<cc<<"Bmsy"<<cc<<"Fmsy"<<cc<<"OFL"<<cc<<"curB";
                fs<<endl;
                fs<<iSeed<<cc<<value(objFun)<<cc<<maxGrad<<cc<<spB_yx(mxYr,MALE);
                if (doOFL) fs<<cc<<ptrOFLResults->B0<<cc<<ptrOFLResults->Bmsy<<cc<<ptrOFLResults->Fmsy<<cc<<ptrOFLResults->OFL<<cc<<ptrOFLResults->curB;
                fs<<endl;
                fs.close();
            }
        } else if (mseOpModMode){
            PRINT2B1("TODO: CREATE A REPORT FILE FOR OpModMode")
        } 
    }
    PRINT2B1("--FINISHED REPORT_SECTION---------")
}

void model_parameters::between_phases_calculations(void)
{
    PRINT2B1(" ")
    PRINT2B1("#--BETWEEN_PHASES_SECTION---------------------")
    adstring msg = "#----Starting phase "+str(current_phase())+" of "+str(initial_params::max_number_phases);
    PRINT2B1(msg)
    if (mseOpModMode){
        dvariable mseCapF = mfexp(pMSE_LnC[1]);
        projectPopForTAC(mseCapF,0,cout);
    } else {
        if (!runAlt) runPopDyMod(0,cout); else runAltPopDyMod(0,cout);
        calcObjFun(dbgObjFun,cout);
        for (int n=1;n<=npLnEffX;n++){
            rpt::echo<<"#----is pLnEffX["<<n<<"] active? "<<active(pLnEffX[n])<<endl;
        }
        if ((current_phase()>=phsItsRewgt)){
            PRINT2B1("#--Calculating effective weights for size compositions")
            //note that the following modifies the objective function value
            ofstream os; 
            if ((current_phase()==phsItsRewgt)){
                os.open("effectiveWeights.R", ios::trunc);
                os<<"effWgts=list("<<endl;
            } else {
                os.open("effectiveWeights.R", ios::app);
            }
            os<<"#--BETWEEN_PHASES_SECTION: Calculating effective weights for size compositions at start of "<<str(current_phase())<<endl;
            os<<"effWgts."<<current_phase()-phsItsRewgt+1<<"=list("<<endl;
            os<<"surveys="  <<endl; calcWeightsForSurveySizeComps( -1,os); os<<","<<endl;
            os<<"fisheries="<<endl; calcWeightsForFisherySizeComps(-1,os); os<<endl;
            os<<"),"<<endl;
            os<<"#--BETWEEN_PHASES_SECTION: Finished calculating effective weights for size compositions"<<endl<<endl;
            if ((ptrMOs->optIterativeReweighting)&&(numItsRewgt<maxItsRewgt)){
                PRINT2B1("#--Re-weighting size compositions using reWeightXXXSizeComps")
                os<<"#--BETWEEN_PHASES_SECTION: Re-weighting size compositions using reWeightXXXSizeComps"<<endl;
                //reWeightSurveySizeComps(-1,os);
                reWeightFisherySizeComps(-1,os);
                os<<"#--BETWEEN_PHASES_SECTION: Finished r-weighting size compositions using reWeightXXXSizeComps"<<endl<<endl;
                PRINT2B1("#--Finished re-weighting size compositions using reWeightXXXSizeComps")
                numItsRewgt++;
                calcObjFun(dbgObjFun,cout);
            }
            os.close();
        }
    }
    ctrProcCallsInPhase=0;//reset in-phase counter
        
    PRINT2B1("#---END BETWEEN_PHASES_SECTION")
}

void model_parameters::reWeightSurveySizeComps(int debug,ostream& cout)
{
    if (debug) cout<<"#--Starting reWeightSurveySizeComps()"<<endl;
    if (ptrMOs->optIterativeReweighting){
        for (int v=1;v<=nSrv;v++){
            if (debug) cout<<"#--reweighting size comps for survey "<<ptrMC->lblsSrv[v]<<endl;
            int vd = mapM2DSrv(v);//get index for survey data corresponding to model survey v
            FleetData* ptrObs = ptrMDS->ppSrv[vd-1];
            if (ptrObs->hasICD){//index catch data
                 if (ptrObs->ptrICD->hasZFD && ptrObs->ptrICD->ptrZFD->optFit){
                    if (debug) cout<<"#---index catch size frequencies"<<endl;
                    ptrObs->ptrICD->ptrZFD->applyReWeightingFactors();
                }
            }
        }//surveys loop
    }
    if (debug) cout<<"#--Finished reWeightSurveySizeComps()"<<endl<<endl;
}

void model_parameters::reWeightFisherySizeComps(int debug,ostream& cout)
{
    if (debug) cout<<"#--Starting reWeightsForFisherySizeComps()"<<endl;
    if (ptrMOs->optIterativeReweighting){
        for (int f=1;f<=nFsh;f++){
            if (debug) cout<<"#---reweighting size comps for fishery "<<ptrMC->lblsFsh[f]<<endl;
            int fd = mapM2DFsh(f);//get index for fishery data corresponding to model fishery f
            FleetData* ptrObs = ptrMDS->ppFsh[fd-1];
            if (ptrObs->hasRCD){//retained catch data
                if (ptrObs->ptrRCD->hasZFD && ptrObs->ptrRCD->ptrZFD->optFit){
                    if (debug) cout<<"#----retained catch size frequencies"<<endl;
                    ptrObs->ptrRCD->ptrZFD->applyReWeightingFactors();
                }
            }
            if (ptrObs->hasTCD){//observed total catch data
                 if (ptrObs->ptrTCD->hasZFD && ptrObs->ptrTCD->ptrZFD->optFit){
                    if (debug) cout<<"#----total catch size frequencies"<<endl;
                    ptrObs->ptrTCD->ptrZFD->applyReWeightingFactors();
                }
            }
            if (ptrObs->hasDCD){//observed discard catch data
                if (ptrObs->ptrDCD->hasZFD && ptrObs->ptrDCD->ptrZFD->optFit){
                    if (debug) cout<<"#----discard catch size frequencies"<<endl;
                    ptrObs->ptrDCD->ptrZFD->applyReWeightingFactors();
                }
            }
        }//fisheries
    }
    if (debug) cout<<"#--Finished reWeightFisherySizeComps()"<<endl;
}

void model_parameters::calcWeightsForSurveySizeComps(int debug, ostream& cout)
{
    if (debug) cout<<"#--Starting calcWeightsForSurveySizeComps()"<<endl;
    int opt = ptrMOs->optIterativeReweighting;
    adstring nDims = "c('N','McAllister-Ianelli','Francis')";
    if (debug<0) cout<<"list("<<endl;
    for (int v=1;v<=nSrv;v++){
        if (debug) cout<<"#--calculating size comps weights for survey "<<ptrMC->lblsSrv[v]<<endl;
        if (debug<0) cout<<"`"<<ptrMC->lblsSrv[v]<<"`=list("<<endl;
        int vd = mapM2DSrv(v);//get index for survey data corresponding to model survey v
        FleetData* ptrObs = ptrMDS->ppSrv[vd-1];
        if (ptrObs->hasICD){//index catch data
            if (debug<0) cout<<"index.catch=list("<<endl;
            if (ptrObs->ptrICD->hasZFD && ptrObs->ptrICD->ptrZFD->optFit){
                smlVal = 0.001;//match to TCSAM2013
                d5_array effWgtComps_xmsyn = calcNLLs_CatchNatZ(ptrObs->ptrICD->ptrZFD,n_vyxmsz(v),0,cout);//don't debug
                if (effWgtComps_xmsyn.size()>0){
                    if (debug>0) cout<<"#--Index catch size comps weights using "<<tcsam::getFitType(ptrObs->ptrICD->ptrZFD->optFit)<<endl;
                    d4_array effWgts_nxms = calcEffWgts(effWgtComps_xmsyn,debug,cout);
                    if (debug>0) {cout<<"effWgts_nxms ="<<endl; wts::print(effWgts_nxms,cout,0); cout<<endl;}
                    if (debug<0) {
                        ivector bnds = wts::getBounds(effWgts_nxms);
                        adstring x = tcsamDims::getSXsForR(bnds(3),bnds(4));
                        adstring m = tcsamDims::getMSsForR(bnds(5),bnds(6));
                        adstring s = tcsamDims::getSCsForR(bnds(7),bnds(8));
                        cout<<"effWgts="<<endl; wts::writeToR(cout, effWgts_nxms,nDims,x,m,s); cout<<cc<<endl;
                    }
                    if (opt) {
                        cout<<"#----Calculating reweighting size comps in call to calcWeightsForSurveySizeComps!!"<<endl;
                        if (debug<0) {
                            cout<<"opt="<<opt<<cc<<endl;
                            cout<<"applied="<<endl;
                        }
                        ptrObs->ptrICD->ptrZFD->calcReWeightingFactors(effWgts_nxms(opt),debug,cout);
                    }
                }//if effWgtComps_xmsyn.size()>0
            }
            if (debug<0) cout<<")"<<endl;
        }//ptrObs->hasICD
        if (debug<0) cout<<"),"<<endl;
    }//surveys loop
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug) cout<<"#--Finished calcWeightsForSurveySizeComps()"<<endl;
}

void model_parameters::calcWeightsForFisherySizeComps(int debug, ostream& cout)
{
    if (debug) cout<<"#--Starting calcWeightsForFisherySizeComps()"<<endl;
    int opt = ptrMOs->optIterativeReweighting;
    adstring nDims = "c('N','McAllister-Ianelli','Francis')";
    if (debug<0) cout<<"list("<<endl;
    for (int f=1;f<=nFsh;f++){
        if (debug) cout<<"#--calculating effective weights for fishery "<<ptrMC->lblsFsh[f]<<endl;
        if (debug<0) cout<<"`"<<ptrMC->lblsFsh[f]<<"`=list("<<endl;
        int fd = mapM2DFsh(f);//get index for fishery data corresponding to model fishery f
        FleetData* ptrObs = ptrMDS->ppFsh[fd-1];
        if (ptrObs->hasRCD){//retained catch data
            if (debug<0) cout<<"retained.catch=list("<<endl;
            if (ptrObs->ptrRCD->hasZFD && ptrObs->ptrRCD->ptrZFD->optFit){
                if (debug) cout<<"#---retained catch size frequencies"<<endl;
                smlVal = 0.001;//match to TCSAM2013
                d5_array effWgtComps_xmsyn = calcNLLs_CatchNatZ(ptrObs->ptrRCD->ptrZFD,rmN_fyxmsz(f),0,cout);
                if (effWgtComps_xmsyn.size()>0){
                    if (debug>0) cout<<"#--Retained catch size comps weights using "<<tcsam::getFitType(ptrObs->ptrRCD->ptrZFD->optFit)<<endl;
                    d4_array effWgts_nxms = calcEffWgts(effWgtComps_xmsyn,debug,cout);
                    if (debug>0) {cout<<"effWgts_nxms ="<<endl; wts::print(effWgts_nxms,cout,0); cout<<endl;}
                    if (debug<0) {
                        ivector bnds = wts::getBounds(effWgts_nxms);
                        adstring x = tcsamDims::getSXsForR(bnds(3),bnds(4));
                        adstring m = tcsamDims::getMSsForR(bnds(5),bnds(6));
                        adstring s = tcsamDims::getSCsForR(bnds(7),bnds(8));
                        cout<<"effWgts="<<endl; wts::writeToR(cout, effWgts_nxms,nDims,x,m,s); cout<<cc<<endl;
                    }
                    if (opt) {
                        cout<<"#----Calculating reweighting size comps in call to calcWeightsForFisherySizeComps!!"<<endl;
                        if (debug<0) {
                            cout<<"opt="<<opt<<cc<<endl;
                            cout<<"applied="<<endl;
                        }
                        ptrObs->ptrRCD->ptrZFD->calcReWeightingFactors(effWgts_nxms(opt),debug,cout);
                    }
                }//if effWgtComps_xmsyn.size()>0
            }
            if (debug<0) cout<<"),"<<endl;
        }//ptrObs->hasRCD
        if (ptrObs->hasTCD){//observed total catch data
            if (debug<0) cout<<"total.catch=list("<<endl;
            if (ptrObs->ptrTCD->hasZFD && ptrObs->ptrTCD->ptrZFD->optFit){
                if (debug) cout<<"#---total catch size frequencies"<<endl;
                smlVal = 0.001;//match to TCSAM2013
                d5_array effWgtComps_xmsyn = calcNLLs_CatchNatZ(ptrObs->ptrTCD->ptrZFD,cpN_fyxmsz(f),0,cout);
                if (effWgtComps_xmsyn.size()>0){
                    if (debug>0) cout<<"#--Total catch size comps weights using "<<tcsam::getFitType(ptrObs->ptrTCD->ptrZFD->optFit)<<endl;
                    d4_array effWgts_nxms = calcEffWgts(effWgtComps_xmsyn,debug,cout);
                    if (debug>0) {cout<<"effWgts_nxms ="<<endl; wts::print(effWgts_nxms,cout,0); cout<<endl;}
                    if (debug<0) {
                        ivector bnds = wts::getBounds(effWgts_nxms);
                        adstring x = tcsamDims::getSXsForR(bnds(3),bnds(4));
                        adstring m = tcsamDims::getMSsForR(bnds(5),bnds(6));
                        adstring s = tcsamDims::getSCsForR(bnds(7),bnds(8));
                        cout<<"effWgts="<<endl; wts::writeToR(cout, effWgts_nxms,nDims,x,m,s); cout<<cc<<endl;
                    }
                    if (opt) {
                        cout<<"#----Calculating reweighting size comps in call to calcWeightsForFisherySizeComps!!"<<endl;
                        if (debug<0) {
                            cout<<"opt="<<opt<<cc<<endl;
                            cout<<"applied="<<endl;
                        }
                        ptrObs->ptrTCD->ptrZFD->calcReWeightingFactors(effWgts_nxms(opt),debug,cout);
                    }
                }//if effWgtComps_xmsyn.size()>0
            }
            if (debug<0) cout<<"),"<<endl;
        }//ptrObs->hasTCD
        if (ptrObs->hasDCD){//observed discard catch data
            if (debug<0) cout<<"discard.catch=list("<<endl;
            if (ptrObs->ptrDCD->hasZFD && ptrObs->ptrDCD->ptrZFD->optFit){
                if (debug) cout<<"#---discard catch size frequencies"<<endl;
                smlVal = 0.001;//match to TCSAM2013
                d5_array effWgtComps_xmsyn = calcNLLs_CatchNatZ(ptrObs->ptrDCD->ptrZFD,dsN_fyxmsz(f),0,cout);
                if (effWgtComps_xmsyn.size()>0){
                    if (debug>0) cout<<"#--Discard catch size comps weights using "<<tcsam::getFitType(ptrObs->ptrDCD->ptrZFD->optFit)<<endl;
                    d4_array effWgts_nxms = calcEffWgts(effWgtComps_xmsyn,debug,cout);
                    if (debug>0) {cout<<"effWgts_nxms ="<<endl; wts::print(effWgts_nxms,cout,0); cout<<endl;}
                    if (debug<0) {
                        ivector bnds = wts::getBounds(effWgts_nxms);
                        adstring x = tcsamDims::getSXsForR(bnds(3),bnds(4));
                        adstring m = tcsamDims::getMSsForR(bnds(5),bnds(6));
                        adstring s = tcsamDims::getSCsForR(bnds(7),bnds(8));
                        cout<<"effWgts="<<endl; wts::writeToR(cout, effWgts_nxms,nDims,x,m,s); cout<<cc<<endl;
                    }
                    if (opt) {
                        cout<<"#----Calculating reweighting size comps in call to calcWeightsForFisherySizeComps!!"<<endl;
                        if (debug<0) {
                            cout<<"opt="<<opt<<cc<<endl;
                            cout<<"applied="<<endl;
                        }
                        ptrObs->ptrDCD->ptrZFD->calcReWeightingFactors(effWgts_nxms(opt),debug,cout);
                    }
                }//if effWgtComps_xmsyn.size()>0
            }
            if (debug<0) cout<<"),"<<endl;
        }//ptrObs->hasDCD
        if (debug<0) cout<<"NULL),"<<endl;
    }//fisheries
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug) cout<<"#--Finished calcWeightsForFisherySizeComps()"<<endl;
}

d4_array model_parameters::calcEffWgts(d5_array& effWgtComps,int debug, ostream& cout)
{
    if (debug>0) cout<<"Starting calcEffWgts()"<<endl;
    ivector d = wts::getBounds(effWgtComps);
    d4_array effWgts_nxms(0,2,d[1],d[2],d[3],d[4],d[5],d[6]);
    //calculate McAllister-Ianelli weights using the harmonic mean (1.B in Punt, 2017)
    for (int x=d[1];x<=d[2];x++){
        for (int m=d[3];m<=d[4];m++){
            for (int s=d[5];s<=d[6];s++){
                double N = sum(column(effWgtComps(x,m,s),1));
                double harmn = 0.0;
                for (int y=effWgtComps(x,m,s).indexmin();y<=effWgtComps(x,m,s).indexmax();y++){
                    if (effWgtComps(x,m,s,y,1)>0.0) harmn += 1.0/effWgtComps(x,m,s,y,2);
                }
                effWgts_nxms(0,x,m,s) = N;      //number of actual size comps
                if (N>0){//to avoid nan's
                    effWgts_nxms(1,x,m,s) = N/harmn;//harmonic mean  of annual Mc-I tuning weights 
                } else effWgts_nxms(1,x,m,s) = 1.0;
            }
        }
    }
    //calculate Francis weights  (1.C in Punt, 2017))
    for (int x=d[1];x<=d[2];x++){
        for (int m=d[3];m<=d[4];m++){
            for (int s=d[5];s<=d[6];s++){
                double N = sum(column(effWgtComps(x,m,s),1));   //number of actual size comps
                int Np   = column(effWgtComps(x,m,s),1).size(); //total number of years
                dvector effWgt_y = column(effWgtComps(x,m,s),3);//z-scores for Francis tuning weight
                effWgts_nxms(0,x,m,s) = N;
                if (N>0){//to avoid nan's
                    effWgts_nxms(2,x,m,s) = 1.0/(wts::variance(effWgt_y)*Np/N);//Francis tuning weight
                } else effWgts_nxms(2,x,m,s) = 1.0;
            }
        }
    }
    if (debug>0) wts::print(effWgts_nxms,cout,1);
    if (debug>0) cout<<"Finished calcEffWgts()"<<endl;
    return effWgts_nxms;
}

void model_parameters::save_params(void)
{
    adstring fn = "tcsam02."+str(current_phase())+"."+str(ctrProcCallsInPhase)+".par";
    ofstream os1(fn, ios::trunc);
    os1.precision(12);
    initial_params::save(os1, 12);
    os1.close();
}

void model_parameters::calcDynB0(int debug, ostream& cout)
{
    if (debug>0) cout<<"starting calcDynB0()"<<endl;
    //open file for output
    adstring fn = "DynamicB0.R";
    ofstream os(fn, ios::trunc);
    os.precision(12);
    //write MMB from final phase to R
    os<<"res=list("<<endl;
    os<<"MB_yx    ="; wts::writeToR(os,value(spB_yx), yDms,xDms); os<<cc<<endl;
    if (debug<0){
        cout<<"res=list("<<endl;
        cout<<"MB_yx    ="; wts::writeToR(cout,value(spB_yx), yDms,xDms); cout<<cc<<endl;
    }
    //initialize population model
    if (!runAlt) initPopDyMod(0, cout); else initAltPopDyMod(0,cout);
    //reset fishery F's to 0
    hasF_fy.initialize();   //flags indicating whether or not fishery occurs
    hmF_fy.initialize();    //handling mortality
    cpF_fyxms.initialize(); //fully-selected capture rate
    cpF_fyxmsz.initialize();//size-specific capture rate
    rmF_fyxmsz.initialize();//retained mortality rate
    dmF_fyxmsz.initialize();//discard mortality rate
    tmF_yxmsz.initialize(); //total mortality rate
    //run population model
    if (!runAlt){ 
        for (int y=mnYr;y<=mxYr;y++) runPopDyModOneYear(y,0,cout);
    } else {
        for (int y=mnYr;y<=mxYr;y++) 
            for (int x=1;x<=nSXs;x++) runAltPopDyModOneYear(y,x,0,cout);
    }
    //write MMB from dynamic B0 calculations to R
    os<<"dB0_yx    ="; wts::writeToR(os,value(spB_yx), yDms,xDms); os<<endl;
    os<<")"<<endl;
    if (debug<0){
        cout<<"dB0_yx    ="; wts::writeToR(cout,value(spB_yx), yDms,xDms); cout<<endl;
        cout<<")"<<endl;
    }
    os.close();
    if (debug>0) cout<<"finished calcDynB0()"<<endl;
}

void model_parameters::writeStateForOpMod(int y,ostream& os)
{
    if (!mseOpModMode){
        os<<"#--OpMod state for projecting year "<<y+1<<"--"<<endl;
        os<<mnYr    <<tb<<"#start year for recruitment"  <<endl;
        os<<y+1     <<tb<<"#year for projection"<<endl;
        os<<dtF_y(y)<<tb<<"#dtF"<<endl;
        os<<dtM_y(y)<<tb<<"#dtM"<<endl;
        os<<"#wAtZ_xmz:"<<endl; wts::print(ptrMDS->ptrBio->wAtZ_xmz,os,1); os<<endl;
        os<<"#R_y:"<<endl<<R_y(mnYr,y)<<endl;
        os<<"#R_x:"<<endl<<R_yx(y)<<endl;
        os<<"#R_z:"<<endl<<R_yz(y)<<endl;
        os<<"#M_xmsz:"   <<endl; wts::print(M_yxmsz(y),os,1);    os<<endl;
        os<<"#prGr_xszz:"<<endl; wts::print(prGr_yxszz(y),os,1); os<<endl;
        os<<"#prM2M_xz:" <<endl; wts::print(prM2M_yxz(y),os,1);  os<<endl;
        os<<"#hmF_f:"<<endl<<column(hmF_fy,y)<<endl;
        os<<"#cpF_fxmsz:"<<endl;
        for (int f=1;f<=nFsh;f++) {os<<"  #"<<f<<endl; wts::print(cpF_fyxmsz(f,y),os,1); os<<endl;}
        os<<"#ret_fxmsz:"<<endl;
        for (int f=1;f<=nFsh;f++) {os<<"  #"<<f<<endl; wts::print(ret_fyxmsz(f,y),os,1); os<<endl;}
        os<<"#sel_fxmsz:"<<endl;
        for (int f=1;f<=nFsh;f++) {os<<"  #"<<f<<endl; wts::print(sel_fyxmsz(f,y),os,1); os<<endl;}
        os<<"#q_vxmsz:"<<endl;
        for (int v=1;v<=nSrv;v++) {os<<"  #"<<v<<endl; wts::print(q_vyxmsz(v,y),os,1); os<<endl;}
        os<<"#n_xmsz (pop state at start of year to be projected):"<<endl; wts::print(n_yxmsz(y+1),os,1); os<<endl;
        // Add section in for Mature Biomass Average calculations
        //PRINT2B1("Adding code for mature biomass averaging in state file")
        os<<"#MB_ave:"<<endl; wts::print(spB_yx,os,1); os<<endl;
    } else {
        os<<"#--OpMod state for projecting year "<<ptrOMI->mxYr+1<<"--"<<endl;
        os<<ptrOMI->mnYr  <<tb<<"#start year for recruitment"<<endl;
        os<<ptrOMI->mxYr+1<<tb<<"#year for projection"       <<endl;
        os<<ptrOMI->dtF<<tb<<"#dtF"<<endl;
        os<<ptrOMI->dtM<<tb<<"#dtM"<<endl;
        os<<"#wAtZ_xmz:"<<endl; wts::print(ptrOMI->wAtZ_xmz,os,1); os<<endl; 
        os<<"#R_y:"<<endl<<ptrOMI->R_y(mnYr,mxYr)<<tb<<prjR<<endl;
        os<<"#R_x:"<<endl<<ptrOMI->R_x<<endl;
        os<<"#R_z:"<<endl<<ptrOMI->R_z<<endl;
        os<<"#M_xmsz:"   <<endl; wts::print(ptrOMI->M_xmsz,   os,1); os<<endl;
        os<<"#prGr_xszz:"<<endl; wts::print(ptrOMI->prGr_xszz,os,1); os<<endl;
        os<<"#prM2M_xz:" <<endl; wts::print(ptrOMI->prM2M_xz, os,1); os<<endl;
        os<<"#hmF_f:"    <<endl<<ptrOMI->hmF_f<<endl;
        os<<"#cpF_fxmsz:"<<endl; wts::print(ptrOMI->cpF_fxmsz,os,1); os<<endl;
        os<<"#ret_fxmsz:"<<endl; wts::print(ptrOMI->ret_fxmsz,os,1); os<<endl;
        os<<"#sel_fxmsz:"<<endl; wts::print(ptrOMI->sel_fxmsz,os,1); os<<endl;
        os<<"#q_vxmsz:"  <<endl; wts::print(ptrOMI->q_vxmsz,  os,1); os<<endl;
        os<<"#n_xmsz (pop state at start of year to be projected):"<<endl;         
          wts::print(prj_n_xmsz,os,1); os<<endl;
        // Add section in for Mature Biomass Average calculations
       // PRINT2B1("Adding code for mature biomass averaging in state file part2")
        os<<"#MB:"<<endl;wts::print(ptrOMI->spB_yx,os,1); os<<endl;
    }
}

void model_parameters::writeEstModPinFile(int closed, ostream& os)
{
    updateMPI(0,cout);
    ptrMPI->addNextYearToInfo(closed);
    os<<"#####---EstMod Pin File"<<endl;
    ptrMPI->writePin(os);
}

int model_parameters::calcTAC(int hcr, double OFL)
{
    int closed = 1;
    double TAC = 0.0;
    adstring info;
    if (hcr==1){ 
      // Identify Biomass 
        dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
        PRINT2B2("MMB_estmod=",MMB)
        PRINT2B2("MFB_estmod=",MFB)
      // Identify Biomass Average
        ivector perm(1,2); perm[1]=2;perm[2]=1; //transposing to move from sex and year to yr and sex 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        //double aveMFB = mean(vspB_xy(FEMALE)(1982,2017)); // test
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        PRINT2B2("aveMFB_estmod=",aveMFB)
      //TAC Call
        TAC = HarvestStrategies::HCR1_FemaleRamp(MFB, aveMFB, MMB);
        info = "#--HCR1: MFB = "+str(MFB)+cc+"aveMFB = "+str(aveMFB)+cc+"ratio = "+str(MFB/aveMFB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==2){  
    //Identify Biomass
        dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        PRINT2B2("MMB_estmod=",MMB)
    //Identify Biomass Average
        ivector perm(1,2); perm[1]=2;perm[2]=1; //what does this line do?
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));// HCRx,copy in for HCR2 in model options
        PRINT2B2("aveMMB_estmod=",aveMMB)
        PRINT2B1("#----Establish exploitation rate rampID")
        int rampID = ptrMOs->HCR2_rampID; // HAVE BUCK CHECK THIS rampID 
        PRINT2B1("#----Testing rampID")
        PRINT2B2("rampID=",rampID)
        TAC = HarvestStrategies::HCR2_MaleRamp(MMB, aveMMB, rampID);
        info = "#--HCR2: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==22){         
        PRINT2B1("Starting HCR22")
        //Identify survey Biomass for the year 
        d6_array vn_vyxmsz = wts::value(n_vyxmsz);
        d6_array vb_vyxmsz = tcsam::calcBiomass(vn_vyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
        //PRINT2B2("array of biomass=", vb_vyxmsz)
        double MMB = 0;
            for (int s = 1; s<=2; s++){
                for(int z=1;z<=nZBs;z++){
                    MMB += vb_vyxmsz(1,mxYr,MALE,MATURE,s,z);
                   }//z is size
                }//s is shell
        PRINT2B2("MMB=", MMB)
        //Indentify Long Term Average 
        dmatrix vspB_yx = value(spB_yx);
       // PRINT2B2("SpawnB=", vspB_yx)
        ivector perm(1,2); perm[1]=2;perm[2]=1; //what does this line do?
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));// HCRx,copy in for HCR2 in model options
        PRINT2B2("aveMMB_estmod=",aveMMB)
        TAC = HarvestStrategies::HCR22_MaleRamp_SurvEst(MMB, aveMMB);
        info = "#--HCR22: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==23){      
        // Define MMB 
        d3_array n_msz     = value(this->n_vyxmsz(1,mxYr,MALE)); // ESTIMATED male survey abundance in final year
        dmatrix w_mz       = value(ptrMDS->ptrBio->wAtZ_xmz(MALE));   //weight at size
        dvector MMB_z(20,32);MMB_z.initialize();
        for (int m=1;m<=nMSs;m++){                       //loops over maturity state
            for (int s=1;s<=nMSs;s++){                   //loops over shell condition
                MMB_z += elem_prod(n_msz(m,s)(20,32),w_mz(m)(20,32));
            }
        }                                
        double MMB = sum(MMB_z);
        // Define MMBAve
        dmatrix vspB_yx = value(spB_yx);
        ivector perm(1,2); perm[1]=2;perm[2]=1; //what does this line do?
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));// HCRx,copy in for HCR2 in model options
        //int rampID = ptrMOs->HCR2_rampID; // HAVE BUCK CHECK THIS rampID 
        TAC = HarvestStrategies::HCR23_MaleRamp_ModSurvEst(MMB, aveMMB);
        info = "#--HCR23: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==3){                             
       //double buffer = ptrMOs->HCR3_buffer; 
       // HARD CODED BUFFER DUE TO ODD ESTIMATION OF VERY HIGH BUFFER, ASK BUCK HOW TO FIX 
       PRINT2B1("-- BUFFER HAS BEEN HARD CODED AT 0.20!!!!!!!!! --")
       double buffer = 0.20;
       PRINT2B2 ("buffer=", buffer)
       PRINT2B2 ("OFL=", OFL)
        TAC = HarvestStrategies::HCR3_ABC(OFL, buffer);
        info = "#--HCR3: buffer = "+str(buffer)+cc+"OFL = "+str(OFL)+cc+"TAC = "+str(TAC);
    }
    if (hcr==4){                             
    // Define Biomass 
       dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
   // Define Average Biomass
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
         info = "#--HCR4: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB);
    }
    if (hcr==41){                             
    // Define Biomass 
       dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
   // Define Average Biomass
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
         info = "#--HCR41: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB);
    }
    if (hcr==42){                             
    // Define Biomass 
       dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
   // Define Average Biomass
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
         info = "#--HCR42: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB);
    }
    if (hcr==43){                             
    // Define Biomass 
       dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
   // Define Average Biomass
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
         info = "#--HCR43: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB);
    }
    if (hcr==5){                             
    // Define Biomass
       dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
    // Define Average Biomass
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR5_FemaleBlocks(MFB, aveMFB, MMB, aveMMB);
        info = "#--HCR5: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB);
    }
    if (hcr==6){                             
    // Define Biomass 
        dmatrix vspB_yx = value(spB_yx); 
        ivector perm(1,2); perm[1]=2;perm[2]=1; // creating a 2d vector perm <- 2,1 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx); 
        //d3_array weights = ptrMDS->ptrBio->wAtZ_xmz; // array from pointer to a pointer with weights by sex, maturity, size   
        dvector weights = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE);
        double newshell = value(sum(n_yxmsz(mxYr,MALE,MATURE, NEW_SHELL)(20,32))); //n_yxmsz is 5d array by year, sex, maturity, stage, size
        double total = 0.0;
        for (int s = 1; s<=nSCs; s++){
             total += value(sum(n_yxmsz(mxYr, MALE, MATURE, s)(20,32))); 
            }
        double propNS = newshell/total;
        dvector abundELM(20,32);
        abundELM.initialize();
        for (int s = 1; s<=nSCs; s++){
            abundELM += value(n_yxmsz(mxYr, MALE, MATURE, s)(20,32)); 
            }
        double xpRate = ptrMOs->HCR6_xpRate; // exploitation rate
        double sOS = ptrMOs->HCR6_sOS; // old shell selectivity
        TAC = HarvestStrategies::HCR6_ELM(propNS,abundELM, weights, sOS, xpRate);
        info = "#--HCR6: xpRate = "+str(xpRate)+cc+"sOs = "+str(sOS)+cc+"TAC = "+str(TAC); 
    }
    if (hcr==7){  
        double Fmsy        = value(ptrOFLResults->Fmsy);
        d3_array selF_msz  = value(ptrOFLResults->pCIM->selF_fmsz(1));//pull out selectivity for directed fishery
        d3_array M_msz     = value(ptrOFLResults->pPDIM->M_msz);      //natural mortality
        d3_array n_msz     = value(this->n_vyxmsz(1,mxYr,MALE));      //ESTIMATED male survey abundance in final year
        dmatrix w_mz       = value(ptrMDS->ptrBio->wAtZ_xmz(MALE));   //weight at size
        dvector cpB_z(20,32); cpB_z.initialize();
        for (int m=1;m<=nMSs;m++){                       //loops over maturity state
            for (int s=1;s<=nSCs;s++){                   //loops over shell condition
                cpB_z += elem_prod(
                    elem_prod(
                      elem_prod(n_msz(m,s)(20,32),  exp(-0.625*M_msz(m,s)(20,32))),
                      (1-exp(-Fmsy* selF_msz(m,s)(20,32)))),  
                        w_mz(m)(20,32));       
            }
        }
        double CWmsy = sum(cpB_z);
        //now calculate TAC using the HCR
        // Identify Biomass
        dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
        // Identify Average Biomass 
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR7_StatusQuo(MFB,aveMFB, MMB, aveMMB, CWmsy); 
        // Half TAC Rule
        if(rmN_fyxmsz(1, mxYr, MALE, MATURE, 1, (20,32)) == 0) TAC = 0.5*TAC; // Half TAC rule VERIFY THIS DOES WHAT IT'S SUPPOSED TO
        info = "#--HCR7: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB)+cc+"CWmsy = "+str(CWmsy);
    }
            // CAP TAC at 50% of ELMB
        dvector weights = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE);
        double newshell = value(sum(n_yxmsz(mxYr,MALE,MATURE, NEW_SHELL)(20,32))); //n_yxmsz is 5d array by year, sex, maturity, stage, size
        double total = 0.0;
        for (int s = 1; s<=nSCs; s++){
             total += value(sum(n_yxmsz(mxYr, MALE, MATURE, s)(20,32))); 
            }
        double propNS = newshell/total;
        dvector abundELM(20,32);
        abundELM.initialize();
        for (int s = 1; s<=nSCs; s++){
            abundELM += value(n_yxmsz(mxYr, MALE, MATURE, s)(20,32)); 
            }
        PRINT2B1("#GET ELM")
        dvector ELM(20,32);
        ELM.initialize();
        ELM = (propNS*abundELM)+(0.40*(1-propNS)*abundELM);
        PRINT2B2("#ELM=", ELM) 
        PRINT2B1("#GET ELMB")
        double ELMB = 0;
        //Abundance*weight for Exploitable Legal Male Biomass
            for(int i = ELM.indexmin(); i<=ELM.indexmax(); ++i){ // check this indexing format     
                ELMB += ELM[i]*weights[i];
            }
        PRINT2B2("#ELMB=", ELMB) 
        double maxTAC = 0;
        double maxTAC30 = 0;
        maxTAC = 0.5*ELMB;
        PRINT2B2("#maxTAC=", maxTAC) 
        maxTAC30 = 0.3*ELMB; // SET MAX TAC at 30% for Third Dimmer iteration hcr 43
    if(hcr==3){ TAC = TAC; //uncap TAC
    }else if(hcr==43){ if (TAC>maxTAC30) TAC=maxTAC30; // SET MAX TAC at 30% ELM for HCR 43
        PRINT2B2("#maxTAC30=", maxTAC30)
    }else{if (TAC>maxTAC) TAC=maxTAC;} // SET MAX TAC
    if (TAC>0.0) closed=0;  // If there is a TAC the fishery is not closed 
    //--save TAC and OFL to file for OpMod to read
    adstring fn = "TAC_"+str(mxYr+1)+".txt";
    ofstream os; os.open(fn, ios::trunc);
    os<<"#---TAC, OFL for MSE OpMod"<<endl;
    os<<"# TAC      OFL"<<endl;
    os<<TAC<<tb<<OFL<<endl;
    os<<endl;
    os<<info<<endl;
    os.close();
    return(closed);
}

double model_parameters::repTAC(int hcr, double OFL)
{
    double TAC = 0.0;
    adstring info;
    if (hcr==1){       
       // ID Biomass 
        dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        PRINT2B2("MMB_estmodTAC=", MMB)
        double MFB = vspB_yx(mxYr,FEMALE);
        PRINT2B2("MFB_estmodTAC=", MFB)
       //ID AveBiomass
        ivector perm(1,2); perm[1]=2;perm[2]=1; //transposing to move from sex and year to yr and sex 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        //double aveMFB = mean(vspB_xy(FEMALE)(1982,2017));
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        PRINT2B2("aveMFB_estmodTAC=", aveMFB)
       // TAC 
        TAC = HarvestStrategies::HCR1_FemaleRamp(MFB, aveMFB, MMB);
        PRINT2B2("tac=", TAC)
        info = "#--HCR1: MFB = "+str(MFB)+cc+"aveMFB = "+str(aveMFB)+cc+"ratio = "+str(MFB/aveMFB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==2){     
        //Identify Biomass
        dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        PRINT2B2("MMB_estmod=",MMB)
        //Identify Biomass Average
        ivector perm(1,2); perm[1]=2;perm[2]=1; //what does this line do?
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));// HCRx,copy in for HCR2 in model options
        PRINT2B2("aveMMB_estmod=",aveMMB)
        PRINT2B1("#----Establish exploitation rate rampID")
        int rampID = ptrMOs->HCR2_rampID; // HAVE BUCK CHECK THIS rampID 
        PRINT2B1("#----Testing rampID")
        PRINT2B2("rampID=",rampID)
        TAC = HarvestStrategies::HCR2_MaleRamp(MMB, aveMMB, rampID);
        info = "#--HCR2: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==22){      
             PRINT2B1("Starting HCR22")
        //Identify survey Biomass for the year 
        d6_array vn_vyxmsz = wts::value(n_vyxmsz);
        d6_array vb_vyxmsz = tcsam::calcBiomass(vn_vyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
        //PRINT2B2("array of biomass=", vb_vyxmsz)
        double MMB = 0;
            for (int s = 1; s<=2; s++){
                for(int z=1;z<=nZBs;z++){
                    MMB += vb_vyxmsz(1,mxYr,MALE,MATURE,s,z);
                   }//z is size
                }//s is shell
        //PRINT2B2("MMB=", MMB)
        //Indentify Long Term Average 
        dmatrix vspB_yx = value(spB_yx);
        //PRINT2B2("SpawnB=", vspB_yx)
        ivector perm(1,2); perm[1]=2;perm[2]=1; //what does this line do?
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));// HCRx,copy in for HCR2 in model options
        PRINT2B2("aveMMB_estmod=",aveMMB)
        TAC = HarvestStrategies::HCR22_MaleRamp_SurvEst(MMB, aveMMB);
        info = "#--HCR22: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
 }
    if (hcr==23){    
        // Define MMB 
        d3_array n_msz     = value(this->n_vyxmsz(1,mxYr,MALE)); // ESTIMATED male survey abundance in final year
        dmatrix w_mz       = value(ptrMDS->ptrBio->wAtZ_xmz(MALE));   //weight at size
        dvector MMB_z(20,32);MMB_z.initialize();
        for (int m=1;m<=nMSs;m++){                       //loops over maturity state
            for (int s=1;s<=nMSs;s++){                   //loops over shell condition
                MMB_z += elem_prod(n_msz(m,s)(20,32),w_mz(m)(20,32));
        }
        }                                
        double MMB = sum(MMB_z);
        // Define MMBAve
        dmatrix vspB_yx = value(spB_yx);
        ivector perm(1,2); perm[1]=2;perm[2]=1; //what does this line do?
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));// HCRx,copy in for HCR2 in model options
        //int rampID = ptrMOs->HCR2_rampID; // HAVE BUCK CHECK THIS rampID 
        TAC = HarvestStrategies::HCR23_MaleRamp_ModSurvEst(MMB, aveMMB);
        info = "#--HCR23: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==3){  
       // ID buffer from model options 
       //double buffer = ptrMOs->HCR3_buffer; 
       // HARD CODED BUFFER DUE TO BUFFER at 0 or ODD ESTIMATION OF VERY HIGH BUFFER, ASK BUCK HOW TO FIX 
       PRINT2B1("-- BUFFER HAS BEEN HARD CODED AT 0.20!!!!!!!!! --")
       double buffer = 0.20;
       PRINT2B2 ("buffer=", buffer)
       //does the OFL function need to be called?
       PRINT2B2 ("OFL=", OFL)
        TAC = HarvestStrategies::HCR3_ABC(OFL, buffer);
        info = "#--HCR3: buffer = "+str(buffer)+cc+"OFL = "+str(OFL)+cc+"TAC = "+str(TAC);
        }
    if (hcr==4){                             
       dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR4: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
    }
    if (hcr==41){                             
    // Define Biomass 
       dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
   // Define Average Biomass
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
         info = "#--HCR41: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB);
    }
    if (hcr==42){                             
    // Define Biomass 
       dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
   // Define Average Biomass
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
         info = "#--HCR42: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB);
    }
    if (hcr==42){                             
    // Define Biomass 
       dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
   // Define Average Biomass
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
         info = "#--HCR43: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB);
    }
    if (hcr==5){                             
       dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR5_FemaleBlocks(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR5: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
    }
    if (hcr==6){                             
        dmatrix vspB_yx = value(spB_yx); 
        ivector perm(1,2); perm[1]=2;perm[2]=1; // creating a 2d vector perm <- 2,1 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx); 
        //d3_array weights = ptrMDS->ptrBio->wAtZ_xmz; // array from pointer to a pointer with weights by sex, maturity, size   
        dvector weights = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE);
        double newshell = value(sum(n_yxmsz(mxYr,MALE,MATURE, NEW_SHELL)(20,32))); //n_yxmsz is 5d array by year, sex, maturity, stage, size
        double total = 0.0;
        for (int s = 1; s<=nSCs; s++){
             total += value(sum(n_yxmsz(mxYr, MALE, MATURE, s)(20,32))); 
            }
        double propNS = newshell/total;
        dvector abundELM(20,32);
        abundELM.initialize();
        for (int s = 1; s<=nSCs; s++){
            abundELM += value(n_yxmsz(mxYr, MALE, MATURE, s)(20,32)); 
            }
        double xpRate = ptrMOs->HCR6_xpRate; // exploitation rate
        double sOS = ptrMOs->HCR6_sOS; // old shell selectivity
        TAC = HarvestStrategies::HCR6_ELM(propNS,abundELM, weights, sOS, xpRate);
        info = "#--HCR6: xpRate = "+str(xpRate)+cc+"sOs = "+str(sOS)+cc+"TAC = "+str(TAC); 
    }
    if (hcr==7){  
        double Fmsy        = value(ptrOFLResults->Fmsy);
        d3_array selF_msz  = value(ptrOFLResults->pCIM->selF_fmsz(1));//pull out selectivity for directed fishery
        d3_array M_msz     = value(ptrOFLResults->pPDIM->M_msz);      //natural mortality
        d3_array n_msz     = value(this->n_vyxmsz(1,mxYr,MALE));      //ESTIMATED male survey abundance in final year
        dmatrix w_mz       = value(ptrMDS->ptrBio->wAtZ_xmz(MALE));   //weight at size
        dvector cpB_z(20,32); cpB_z.initialize();
        for (int m=1;m<=nMSs;m++){                       //loops over maturity state
            for (int s=1;s<=nSCs;s++){                   //loops over shell condition
                cpB_z += elem_prod(
                    elem_prod(
                      elem_prod(n_msz(m,s)(20,32),  exp(-0.625*M_msz(m,s)(20,32))),
                      (1-exp(-Fmsy* selF_msz(m,s)(20,32)))),  
                        w_mz(m)(20,32)); 
                                                PRINT2B2("#CWmsy_check1= ", cpB_z)
                                                PRINT2B2("#selF_check1= ", selF_msz)
                                                PRINT2B2("#M_check1= ", M_msz)
                                                PRINT2B2("#Num_check1=", n_msz)
                                                PRINT2B2("#weight_mz=", w_mz)
                                                PRINT2B2("#Fmsy=", Fmsy)
            }                   
        }
        double CWmsy = sum(cpB_z);
        PRINT2B2("#CWmsy_check2= ", CWmsy)
        //now calculate TAC using the HCR
        dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR7_StatusQuo(MFB,aveMFB, MMB, aveMMB, CWmsy); 
        // Half TAC Rule
        if(rmN_fyxmsz(1, mxYr, MALE, MATURE, 1, (20,32)) == 0) TAC = 0.5*TAC; // Half TAC rule VERIFY THIS DOES WHAT IT'S SUPPOSED TO           
        info = "#--HCR7: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB)+cc+"CWmsy = "+str(CWmsy);
    }
        // CAP TAC at 50% of ELMB
        dvector weights = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE);
        double newshell = value(sum(n_yxmsz(mxYr,MALE,MATURE, NEW_SHELL)(20,32))); //n_yxmsz is 5d array by year, sex, maturity, stage, size
        double total = 0.0;
        for (int s = 1; s<=nSCs; s++){
             total += value(sum(n_yxmsz(mxYr, MALE, MATURE, s)(20,32))); 
            }
        double propNS = newshell/total;
        dvector abundELM(20,32);
        abundELM.initialize();
        for (int s = 1; s<=nSCs; s++){
            abundELM += value(n_yxmsz(mxYr, MALE, MATURE, s)(20,32)); 
            }
        PRINT2B1("#GET ELM")
        dvector ELM(20,32);
        ELM.initialize();
        ELM = (propNS*abundELM)+(0.40*(1-propNS)*abundELM);
        PRINT2B2("#ELM=", ELM) 
        PRINT2B1("#GET ELMB")
        double ELMB = 0;
        //Abundance*weight for Exploitable Legal Male Biomass
            for(int i = ELM.indexmin(); i<=ELM.indexmax(); ++i){ // check this indexing format     
                ELMB += ELM[i]*weights[i];
            }
        PRINT2B2("#ELMB=", ELMB) 
        double maxTAC = 0;
        double maxTAC30 = 0;
        maxTAC = 0.5*ELMB;
        PRINT2B2("#maxTAC=", maxTAC) 
        maxTAC30 = 0.3*ELMB; // SET MAX TAC at 30% for Third Dimmer iteration hcr 43
    if(hcr==3){ TAC = TAC; //uncap TAC
    }else if(hcr==43){ if (TAC>maxTAC30) TAC=maxTAC30; // SET MAX TAC at 30% ELM for HCR 43
        PRINT2B2("#maxTAC30=", maxTAC30)
    }else{if (TAC>maxTAC) TAC=maxTAC;} // SET MAX TAC
    //--save TAC and OFL to file for OpMod to read
    //adstring fn = "TAC_"+str(mxYr+1)+".txt";
    //ofstream os; os.open(fn, ios::trunc);
    //os<<"#---TAC, OFL for MSE OpMod"<<endl;
    //os<<"# TAC      OFL"<<endl;
    //os<<TAC<<tb<<OFL<<endl;
    //os<<endl;
    //os<<info<<endl;
    //os.close();
    return(TAC);
}

int model_parameters::calcTAC_OpMod(int hcr, double OFL)
{
    int closed = 1;
    double TAC = 0.0;
    adstring info;
    if (hcr==1){ 
    //Calculate MMB and MFB from projections 
        double MMB = value(prj_spB_x(MALE));
        PRINT2B2("MMB_opmodTAC=", MMB)
        double MFB = value(prj_spB_x(FEMALE));
        PRINT2B2("MFB_opmodTAC=", MFB)
     //Calculate ave MFB 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
       // PRINT2B2("SpB_yx", vspB_yx);
        ivector perm(1,2); perm[1]=2;perm[2]=1; //transposing to move from sex and year to yr and sex 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        PRINT2B2("aveMFB_opmodTAC=", aveMFB)
        TAC = HarvestStrategies::HCR1_FemaleRamp(MFB, aveMFB, MMB);
        double ratio = MFB/aveMFB;
        info = "#--HCR1: MFB = "+str(MFB)+cc+"aveMFB_test = "+str(aveMFB)+cc+"ratio = "+str(ratio)+cc+"TAC = "+str(TAC);
    }
    if (hcr==2){  
        // Define projected MMB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        PRINT2B2("MMB_Opmod=",MMB)
        //Define Average MMB
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        //PRINT2B2("opmod_spB =", vspB_yx) 
        //PRINT2B2("opmod_spBMale1 =", vspB_yx(ptrMOs->HCR_avgMinYr, MALE)) 
        //PRINT2B2("opmod_spBMale2 =", vspB_yx(ptrMOs->HCR_avgMaxYr, MALE)) 
        //PRINT2B2("opmod_spBFemale1 =", vspB_yx(ptrMOs->HCR_avgMinYr,FEMALE)) 
        //PRINT2B2("opmod_spBFemale2 =", vspB_yx(ptrMOs->HCR_avgMaxYr,FEMALE)) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        PRINT2B2("aveMMB_Opmod=",aveMMB)
        // Ramp for exploitation 
        PRINT2B1("#----Establish exploitation rate rampID")
        int rampID = ptrMOs->HCR2_rampID; // HAVE BUCK CHECK THIS rampID 
        PRINT2B1("#testing rampID")
        PRINT2B2("rampID_opMod=",rampID)
        TAC = HarvestStrategies::HCR2_MaleRamp(MMB, aveMMB, rampID);
        info = "#--HCR2: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==22){  
        // Define projected MMB
        d3_array n_msz     = value(this->prj_n_vxmsz(1,MALE)); // ESTIMATED male survey abundance in final year
        dmatrix w_mz       = value(ptrMDS->ptrBio->wAtZ_xmz(MALE));   //weight at size
        dvector MMB_z(20,32);MMB_z.initialize();
        for (int m=1;m<=nMSs;m++){                       //loops over maturity state
            for (int s=1;s<=nMSs;s++){                   //loops over shell condition
                MMB_z += elem_prod(n_msz(m,s)(20,32),w_mz(m)(20,32));
            }
        }                                
        double MMB = sum(MMB_z);
        //Define Average MMB
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        PRINT2B2("opmod_spB =", vspB_yx) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        PRINT2B2("aveMMB_Opmod=",aveMMB)
        // HCRx,copy in for HCR2 in model options
        //int rampID = ptrMOs->HCR2_rampID; // HAVE BUCK CHECK THIS rampID 
        TAC = HarvestStrategies::HCR22_MaleRamp_SurvEst(MMB, aveMMB);
        info = "#--HCR22: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
        }
    if (hcr==23){  
        // Define projected MMB
        d3_array n_msz     = value(this->prj_n_vxmsz(1,MALE)); // ESTIMATED male survey abundance in final year
        dmatrix w_mz       = value(ptrMDS->ptrBio->wAtZ_xmz(MALE));   //weight at size
        dvector MMB_z(20,32);MMB_z.initialize();
        for (int m=1;m<=nMSs;m++){                       //loops over maturity state
            for (int s=1;s<=nMSs;s++){                   //loops over shell condition
                MMB_z += elem_prod(n_msz(m,s)(20,32),w_mz(m)(20,32));
            }
        }                                
        double MMB = sum(MMB_z);
        //Define Average MMB
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        PRINT2B2("opmod_spB =", vspB_yx) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(1982,2017));
        //int rampID = ptrMOs->HCR2_rampID; // HAVE BUCK CHECK THIS rampID 
        TAC = HarvestStrategies::HCR23_MaleRamp_ModSurvEst(MMB, aveMMB);
        info = "#--HCR23: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==3){                             
       // ID buffer from Model options 
       //double buffer = ptrMOs->HCR3_buffer; 
       // HARD CODED BUFFER DUE TO BUFFER at 0 or ODD ESTIMATION OF VERY HIGH BUFFER, ASK BUCK HOW TO FIX 
       PRINT2B1("-- BUFFER HAS BEEN HARD CODED AT 0.20!!!!!!!!! --")
       double buffer = 0.20;
       PRINT2B2("Buffer_opMod=",buffer)
       //does OFL function need to be called?
       PRINT2B2("OFL_opMod=", OFL)
        TAC = HarvestStrategies::HCR3_ABC(OFL, buffer); // ASK ABOUT THIS LINE, OP OFL?
        info = "#--HCR3: buffer = "+str(buffer)+cc+"OFL = "+str(OFL)+cc+"TAC = "+str(TAC);
    }
    if (hcr==4){                             
        // ID biomass for males and females for  OP MOD using prj_spB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        // Averages for Male and Female Biomass 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        PRINT2B2("opmod_spB =", vspB_yx) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR4: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
        }
    if (hcr==41){                             
        // ID biomass for males and females for  OP MOD using prj_spB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        // Averages for Male and Female Biomass 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        PRINT2B2("opmod_spB =", vspB_yx) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR41: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
        }
    if (hcr==42){                             
        // ID biomass for males and females for  OP MOD using prj_spB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        // Averages for Male and Female Biomass 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        PRINT2B2("opmod_spB =", vspB_yx) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR43: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
        }
    if (hcr==43){                             
        // ID biomass for males and females for  OP MOD using prj_spB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        // Averages for Male and Female Biomass 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        PRINT2B2("opmod_spB =", vspB_yx) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR43: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
        }
    if (hcr==5){      
        // ID biomass for males and females for  OP MOD using prj_spB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        // Averages for Male and Female Biomass 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        PRINT2B2("opmod_spB =", vspB_yx) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR5_FemaleBlocks(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR5: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
    }
    if (hcr==6){                             
        //d3_array weights = ptrMDS->ptrBio->wAtZ_xmz; // array from pointer to a pointer with weights by sex, maturity, size   
        dvector weights = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE);
        //Calc only newshell mature males 
        double newshell = value(sum(prj_n_xmsz(MALE,MATURE, NEW_SHELL)(20,32))); // 
        //Calculate total mature males 
        double total = 0.0;
        for (int s = 1; s<=nSCs; s++){
             total += value(sum(prj_n_xmsz(MALE, MATURE, s)(20,32))); //WHY SUM
            }
        double propNS = newshell/total;
        dvector abundELM(20,32);
        abundELM.initialize();
        for (int s = 1; s<=nSCs; s++){
            abundELM += value(prj_n_xmsz(MALE, MATURE, s)(20,32)); // WHY NO SUM
            }
        double xpRate = ptrMOs->HCR6_xpRate; // exploitation rate
        double sOS = ptrMOs->HCR6_sOS; // old shell selectivity
        PRINT2B2("xpRate=", xpRate)
        PRINT2B2("sOs=", sOS)
        TAC = HarvestStrategies::HCR6_ELM(propNS,abundELM, weights, sOS, xpRate);
        info = "#--HCR6: xpRate = "+str(xpRate)+cc+"sOs = "+str(sOS)+cc+"TAC = "+str(TAC); 
    }
    if (hcr==7){  
        double Fmsy        = value(ptrOFLResults->Fmsy);
        d3_array selF_msz  = value(ptrOFLResults->pCIM->selF_fmsz(1));//pull out selectivity for directed fishery
        d3_array M_msz     = value(ptrOFLResults->pPDIM->M_msz);      //natural mortality
        d3_array n_msz     = value(this->prj_n_vxmsz(1,MALE));      //ESTIMATED male survey abundance in final year
        dmatrix w_mz       = value(ptrMDS->ptrBio->wAtZ_xmz(MALE));   //weight at size
        dvector cpB_z(20,32); cpB_z.initialize();
        for (int m=1;m<=nMSs;m++){                       //loops over maturity state
            for (int s=1;s<=nSCs;s++){                   //loops over shell condition
                cpB_z += elem_prod(
                    elem_prod(
                      elem_prod(n_msz(m,s)(20,32),  exp(-0.625*M_msz(m,s)(20,32))),
                      (1-exp(-Fmsy* selF_msz(m,s)(20,32)))),  
                        w_mz(m)(20,32));        
            }
        }
        double CWmsy = sum(cpB_z);
        PRINT2B2("CWmsy=", CWmsy)
        //now calculate TAC using the HCR
        // Biomass for projectected year 
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        PRINT2B2("opmod_spB =", vspB_yx) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR7_StatusQuo(MFB,aveMFB, MMB, aveMMB, CWmsy); 
         //----------------------------------    
        // Half TAC Rule CHECK ME
        if(prj_rmN_fxmsz(1, MALE, MATURE, 1, (20,32)) == 0) TAC = 0.5*TAC; // Half TAC rule VERIFY THIS DOES WHAT IT'S SUPPOSED TO           
        info = "#--HCR7: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB)+cc+"CWmsy = "+str(CWmsy);
    }
    //  SET TAC CAP FOR HARVEST to 50% ELMB
        PRINT2B1("TAC CAP at 50% ELMB")
        dvector weights = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE);
        //Calc only newshell mature males 
        double newshell = value(sum(prj_n_xmsz(MALE,MATURE, NEW_SHELL)(20,32))); // 
        //Calculate total mature males 
        double total = 0.0;
        for (int s = 1; s<=nSCs; s++){
             total += value(sum(prj_n_xmsz(MALE, MATURE, s)(20,32))); //WHY SUM
            }
        double propNS = newshell/total;
        dvector abundELM(20,32);
        abundELM.initialize();
        for (int s = 1; s<=nSCs; s++){
            abundELM += value(prj_n_xmsz(MALE, MATURE, s)(20,32)); // WHY NO SUM
            }
        PRINT2B1("#GET ELM")
        dvector ELM(20,32);
        ELM.initialize();
        ELM = (propNS*abundELM)+(0.40*(1-propNS)*abundELM);
        PRINT2B2("#ELM=", ELM) 
        PRINT2B1("#GET ELMB")
        double ELMB = 0;
        //Abundance*weight for Exploitable Legal Male Biomass
            for(int i = ELM.indexmin(); i<=ELM.indexmax(); ++i){ // check this indexing format     
                ELMB += ELM[i]*weights[i];
            }
        PRINT2B2("#ELMB=", ELMB) 
        double maxTAC = 0;
        double maxTAC30 = 0;
        maxTAC = 0.5*ELMB;
        PRINT2B2("#maxTAC=", maxTAC) 
        maxTAC30 = 0.3*ELMB; // SET MAX TAC at 30% for Third Dimmer iteration hcr 43
    if(hcr==3){ TAC = TAC; //uncap TAC
    }else if(hcr==43){ if (TAC>maxTAC30) TAC=maxTAC30; // SET MAX TAC at 30% ELM for HCR 43
        PRINT2B2("#maxTAC30=", maxTAC30)
    }else{if (TAC>maxTAC) TAC=maxTAC;} // SET MAX TAC
    if (TAC>0.0) closed=0;  // If there is a TAC the fishery is not closed 
    //--save TAC and OFL to file for OpMod to read
    adstring fn = "TAC_"+str(mxYr+1)+".txt";
    ofstream os; os.open(fn, ios::trunc);
    os<<"#---TAC, OFL for MSE OpMod"<<endl;
    os<<"# TAC      OFL"<<endl;
    os<<TAC<<tb<<OFL<<endl;
    os<<endl;
    os<<info<<endl;
    os.close();
    return(closed);
}

double model_parameters::repTAC_OpMod(int hcr, double OFL)
{
    int closed = 1;
    double TAC = 0.0;
    adstring info;
    if (hcr==1){ 
   //Calculate MMB and MFB from projections 
        double MMB = value(prj_spB_x(MALE));
        PRINT2B2("MMB_opmodTAC=", MMB)
        double MFB = value(prj_spB_x(FEMALE));
        PRINT2B2("MFB_opmodTAC=", MFB)
     //Calculate ave MFB 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        //PRINT2B2("SpB_yx (RepTAC)", vspB_yx);
        ivector perm(1,2); perm[1]=2;perm[2]=1; //transposing to move from sex and year to yr and sex 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR1_FemaleRamp(MFB, aveMFB, MMB);
        double ratio = MFB/aveMFB;
        info = "#--HCR1: MFB = "+str(MFB)+cc+"aveMFB_test = "+str(aveMFB)+cc+"ratio = "+str(ratio)+cc+"TAC = "+str(TAC);
        }
    if (hcr==2){  
        // Define projected MMB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        PRINT2B2("MMB_Opmod=",MMB)
        //Define Average MMB
      //Calculate ave MFB 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        //PRINT2B2("SpB_yx (RepTAC)", vspB_yx);
        ivector perm(1,2); perm[1]=2;perm[2]=1; //transposing to move from sex and year to yr and sex 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        // Ramp for exploitation 
        PRINT2B1("#----Establish exploitation rate rampID")
        int rampID = ptrMOs->HCR2_rampID; // HAVE BUCK CHECK THIS rampID 
        PRINT2B1("#Testing rampID")
        PRINT2B2("rampID_opMod=",rampID)
        TAC = HarvestStrategies::HCR2_MaleRamp(MMB, aveMMB, rampID);
        info = "#--HCR2: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==22){  
        // Define projected MMB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        //Define Average MMB
        //Calculate ave MFB 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        //PRINT2B2("SpB_yx (RepTAC)", vspB_yx);
        ivector perm(1,2); perm[1]=2;perm[2]=1; //transposing to move from sex and year to yr and sex 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR22_MaleRamp_SurvEst(MMB, aveMMB);
        info = "#--HCR2: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==23){  
        // Define projected MMB
        d3_array n_msz     = value(this->prj_n_vxmsz(1,MALE)); // ESTIMATED male survey abundance in final year
        dmatrix w_mz       = value(ptrMDS->ptrBio->wAtZ_xmz(MALE));   //weight at size
        dvector MMB_z(20,32);MMB_z.initialize();
        for (int m=1;m<=nMSs;m++){                       //loops over maturity state
            for (int s=1;s<=nMSs;s++){                   //loops over shell condition
                MMB_z += elem_prod(n_msz(m,s)(20,32),w_mz(m)(20,32));
            }
        }                                
        double MMB = sum(MMB_z);
        //Define Average MMB
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        //PRINT2B2("SpB_yx (RepTAC)", vspB_yx);
        ivector perm(1,2); perm[1]=2;perm[2]=1; //transposing to move from sex and year to yr and sex 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        //int rampID = ptrMOs->HCR2_rampID; // HAVE BUCK CHECK THIS rampID 
        TAC = HarvestStrategies::HCR23_MaleRamp_ModSurvEst(MMB, aveMMB);
        info = "#--HCR2: MMB = "+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"ratio = "+str(MMB/aveMMB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==3){                             
        // ID buffer from Model options 
        //double buffer = ptrMOs->HCR3_buffer; 
        // HARD CODED BUFFER DUE TO BUFFER at 0 or ODD ESTIMATION OF VERY HIGH BUFFER, ASK BUCK HOW TO FIX 
       PRINT2B1("-- BUFFER HAS BEEN HARD CODED AT 0.20!!!!!!!!! --")
       double buffer = 0.20;
       PRINT2B2("Buffer_opMod=",buffer)
       //does OFL function need to be called?
       PRINT2B2("OFL_opMod=", OFL)
        TAC = HarvestStrategies::HCR3_ABC(OFL, buffer); // ASK ABOUT THIS LINE, OP OFL?
        info = "#--HCR3: buffer = "+str(buffer)+cc+"OFL = "+str(OFL)+cc+"TAC = "+str(TAC);
       }
    if (hcr==4){                             
        // ID biomass for males and females for  OP MOD using prj_spB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        // Averages for Male and Female Biomass 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        //PRINT2B2("SpB_yx (RepTAC)", vspB_yx);
        ivector perm(1,2); perm[1]=2;perm[2]=1; //transposing to move from sex and year to yr and sex 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR4: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
          }
    if (hcr==41){                             
        // ID biomass for males and females for  OP MOD using prj_spB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        // Averages for Male and Female Biomass 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        PRINT2B2("opmod_spB =", vspB_yx) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR41: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
        }
    if (hcr==42){                             
        // ID biomass for males and females for  OP MOD using prj_spB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        // Averages for Male and Female Biomass 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        PRINT2B2("opmod_spB =", vspB_yx) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR42: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
        }
    if (hcr==43){                             
        // ID biomass for males and females for  OP MOD using prj_spB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        // Averages for Male and Female Biomass 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        PRINT2B2("opmod_spB =", vspB_yx) 
        ivector perm(1,2); perm[1]=2;perm[2]=1; //permutates year and sex in matrix structure
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR4_FemaleDimmer(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR43: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
        }
    if (hcr==5){      
        // ID biomass for males and females for  OP MOD using prj_spB
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        // Averages for Male and Female Biomass 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        //PRINT2B2("SpB_yx (RepTAC)", vspB_yx);
        ivector perm(1,2); perm[1]=2;perm[2]=1; //transposing to move from sex and year to yr and sex 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR5_FemaleBlocks(MFB, aveMFB, MMB, aveMMB);
        info ="#--HCR5: MFB = "+str(MFB)+cc+"aveMFB= "+str(aveMFB)+cc+"MMB ="+str(MMB)+cc+"aveMMB = "+str(aveMMB)+cc+"TAC = "+str(TAC); 
    }
    if (hcr==6){                             
        //d3_array weights = ptrMDS->ptrBio->wAtZ_xmz; // array from pointer to a pointer with weights by sex, maturity, size   
        dvector weights = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE);
        //Calc only newshell mature males 
        double newshell = value(sum(prj_n_xmsz(MALE,MATURE, NEW_SHELL)(20,32))); // 
        //Calculate total mature males 
        double total = 0.0;
        for (int s = 1; s<=nSCs; s++){
             total += value(sum(prj_n_xmsz(MALE, MATURE, s)(20,32))); //WHY SUM
            }
        double propNS = newshell/total;
        dvector abundELM(20,32);
        abundELM.initialize();
        for (int s = 1; s<=nSCs; s++){
            abundELM += value(prj_n_xmsz(MALE, MATURE, s)(20,32)); // WHY NO SUM
            }
        double xpRate = ptrMOs->HCR6_xpRate; // exploitation rate
        double sOS = ptrMOs->HCR6_sOS; // old shell selectivity
        TAC = HarvestStrategies::HCR6_ELM(propNS,abundELM, weights, sOS, xpRate);
        info = "#--HCR6: xpRate = "+str(xpRate)+cc+"sOs = "+str(sOS)+cc+"TAC = "+str(TAC); 
    }
    // Harvest Control Rule 7--"Status Quo" 
    if (hcr==7){  
        double Fmsy        = value(ptrOFLResults->Fmsy);
        d3_array selF_msz  = value(ptrOFLResults->pCIM->selF_fmsz(1));//pull out selectivity for directed fishery
        //PRINT2B2("selF_msz=", selF_msz)
        d3_array M_msz     = value(ptrOFLResults->pPDIM->M_msz);      //natural mortality
        //PRINT2B2("M_msz=", M_msz)
        d3_array n_msz     = value(this->prj_n_vxmsz(1,MALE));      //ESTIMATED male survey abundance in final year
        //PRINT2B2("n_msz=", n_msz)
        dmatrix w_mz       = value(ptrMDS->ptrBio->wAtZ_xmz(MALE));   //weight at size
        dvector cpB_z(20,32); cpB_z.initialize();
        for (int m=1;m<=nMSs;m++){                       //loops over maturity state
            for (int s=1;s<=nMSs;s++){                   //loops over shell condition
                cpB_z += elem_prod(
                    elem_prod(
                      elem_prod(n_msz(m,s)(20,32),  exp(-0.625*M_msz(m,s)(20,32))),
                      (1-exp(-Fmsy* selF_msz(m,s)(20,32)))),  
                        w_mz(m)(20,32));     
                                                PRINT2B2("#CWmsy_check1= ", cpB_z)
                                                PRINT2B2("#selF_check1= ", selF_msz)
                                                PRINT2B2("#M_check1= ", M_msz)
                                                PRINT2B2("#Num_check1=", n_msz)
                                                PRINT2B2("#weight_mz=", w_mz)
                                                PRINT2B2("#Fmsy=", Fmsy)
            }
        }
        double CWmsy = sum(cpB_z);
        PRINT2B2("CWmsy=",CWmsy)
        //now calculate TAC using the HCR
        // Biomass for projectected year 
        dvector vspB_x = value(prj_spB_x);
        double MMB = vspB_x(MALE);
        double MFB = vspB_x(FEMALE);
        // Mature Biomass Averages 
        dmatrix vspB_yx = value(ptrOMI->spB_yx);
        //PRINT2B2("SpB_yx (RepTAC)", vspB_yx);
        ivector perm(1,2); perm[1]=2;perm[2]=1; //transposing to move from sex and year to yr and sex 
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        double aveMMB = mean(vspB_xy(MALE)(ptrMOs->HCR_avgMinYr,ptrMOs->HCR_avgMaxYr));
        TAC = HarvestStrategies::HCR7_StatusQuo(MFB,aveMFB, MMB, aveMMB, CWmsy); 
    //----------------------------------    
        // Half TAC Rule CHECK ME
        if(prj_rmN_fxmsz(1, MALE, MATURE, 1, (20,32)) == 0) TAC = 0.5*TAC; // Half TAC rule VERIFY THIS DOES WHAT IT'S SUPPOSED TO           
    //------------------------------------ 
    info = "#--HCR7: MMB = "+str(MMB)+cc+"MFB = "+str(MFB)+cc+"aveMMB = "+str(aveMMB)+cc+"aveMFB = "+str(aveMFB)+cc+"CWmsy = "+str(CWmsy);
    }
    PRINT2B1("TAC CAP at 50% ELMB")
        dvector weights = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE);
        //Calc only newshell mature males 
        double newshell = value(sum(prj_n_xmsz(MALE,MATURE, NEW_SHELL)(20,32))); // 
        //Calculate total mature males 
        double total = 0.0;
        for (int s = 1; s<=nSCs; s++){
             total += value(sum(prj_n_xmsz(MALE, MATURE, s)(20,32))); //WHY SUM
            }
        double propNS = newshell/total;
        dvector abundELM(20,32);
        abundELM.initialize();
        for (int s = 1; s<=nSCs; s++){
            abundELM += value(prj_n_xmsz(MALE, MATURE, s)(20,32)); // WHY NO SUM
            }
        PRINT2B1("#GET ELM")
        dvector ELM(20,32);
        ELM.initialize();
        ELM = (propNS*abundELM)+(0.40*(1-propNS)*abundELM);
        //PRINT2B2("#ELM=", ELM) 
        PRINT2B1("#GET ELMB")
        double ELMB = 0;
        //Abundance*weight for Exploitable Legal Male Biomass
            for(int i = ELM.indexmin(); i<=ELM.indexmax(); ++i){ // check this indexing format     
                ELMB += ELM[i]*weights[i];
            }
        PRINT2B2("#ELMB=", ELMB) 
        double maxTAC = 0;
        double maxTAC30 = 0;
        maxTAC = 0.5*ELMB;
        maxTAC30 = 0.3*ELMB; // SET MAX TAC at 30% for Third Dimmer iteration hcr 43
        PRINT2B2("#maxTAC=", maxTAC) 
    if(hcr==3){ TAC = TAC; //uncap TAC
    }else if(hcr==43){ if (TAC>maxTAC30) TAC=maxTAC30; // SET MAX TAC at 30% ELM for HCR 43
    }else{if (TAC>maxTAC) TAC=maxTAC;} // SET MAX TAC
    if (TAC>0.0) closed=0;  // If there is a TAC the fishery is not closed 
    //--save TAC and OFL to file for OpMod to read
    adstring fn = "TAC_"+str(mxYr+1)+".txt";
    ofstream os; os.open(fn, ios::trunc);
    os<<"#---TAC, OFL for MSE OpMod"<<endl;
    os<<"# TAC      OFL"<<endl;
    os<<TAC<<tb<<OFL<<endl;
    os<<endl;
    os<<info<<endl;
    os.close();
    return(TAC);
}

void model_parameters::finishOpModMode(void)
{
    PRINT2B1("#----MSE OpModMode: Recalculating population dynamics")
    if (inpTAC>0.0){
        dvariable mseCapF = mfexp(pMSE_LnC[1]);//exponentiated version of estimated param (by opmod)
        projectPopForTAC(mseCapF,0,cout);
        calcObjFunForTAC(dbgObjFun,cout); //recalc
        PRINT2B2("#--Final obj fun = ",objFun)
    } else {
        projectPopForZeroTAC(0,cout); //if TAC 0 project w/o dir fishery
    }
    //--define parent folder for subsequent files
    ptrMC->mxYr   = mxYr+1;
    ptrMC->asYr   = ptrMC->asYr+1;
    adstring nxtOpModFolder  = "../"+str(ptrMC->asYr)+".OpMod";
    adstring nxtEstModFolder = "../"+str(ptrMC->asYr)+".EstMod";
    //--write model configuration file for OpMod in mxYr+1
    adstring nwMC = wts::concatenateFilePaths(nxtOpModFolder,
                                              "OpMod.Configuration.inp");
    PRINT2B2("writing MC to ",nwMC)
    ptrMC->write(nwMC);
    //--write model datasets file and dataset files for estimation model in mxYr+1
    //--bio data
    {
        adstring nwF = wts::concatenateFilePaths(nxtEstModFolder,ptrMDS->fnBioData);
        ofstream ofs; ofs.open(nwF,ios::trunc);
        ptrMDS->ptrBio->write(ofs);
        ofs.close();
    }
    //--fishery data
    for (int i=1;i<=ptrMDS->nFsh;i++){
        if ((i!=1)||(inpTAC>0.0)){
            //update fishery i (skipping directed fishery if TAC=0)
            adstring name  = ptrMDS->ppFsh[i-1]->name;
            //update fishery data for year mxYr+1
            FleetData::debug=0;
            ptrMDS->ppFsh[i-1]->addFisheryCatchData(ptrMC->mxYr,
                                                    prj_cpN_fxmsz(i),
                                                    prj_rmN_fxmsz(i),
                                                    ptrMDS->ptrBio->wAtZ_xmz,
                                                    0.0,
                                                    0.0,
                                                    rng);
        }
        //write new fishery data
        adstring nwF = wts::concatenateFilePaths(nxtEstModFolder,ptrMDS->fnsFisheryData[i]);
        ofstream ofs; ofs.open(nwF,ios::trunc);
        ptrMDS->ppFsh[i-1]->write(ofs);
        ofs.close();
        FleetData::debug=0;
    }
    //--survey data
    //----calculate survey data for year mxYr+1
    cout<<"Calculating survey data for for year "<<ptrMC->mxYr+1<<endl;
    prj_n_vxmsz.initialize();
    for (int v=1;v<=nSrv;v++){ //survey # (only one in this case) 
        for (int x=1;x<=nSXs;x++){ 
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++){
                    //NOTE: using survey characteristics from final year of "real" data
                    prj_n_vxmsz(v,x,m,s) = elem_prod(ptrOMI->q_vxmsz(v,x,m,s),prj_n_xmsz(x,m,s));
                }
            }
        }
    }
    cout<<"Calculated survey data for for year "<<ptrMC->mxYr+1<<endl;
    //----write survey data files
    for (int i=1;i<=ptrMDS->nSrv;i++){
        adstring name  = ptrMDS->ppSrv[i-1]->name;
        //update survey data for year mxYr+1
        FleetData::debug=0;
        ptrMDS->ppSrv[i-1]->addIndexCatchData(ptrMC->mxYr+1,
                                              prj_n_vxmsz(i),
                                              ptrMDS->ptrBio->wAtZ_xmz,
                                              0.0,
                                              0.0,
                                              rng);
        //write new survey data
        adstring nwF = wts::concatenateFilePaths(nxtEstModFolder,ptrMDS->fnsSurveyData[i]);
        ofstream ofs; ofs.open(nwF,ios::trunc);
        ptrMDS->ppSrv[i-1]->write(ofs);
        ofs.close();
        FleetData::debug=0;
    }
    //--growth data
    for (int i=1;i<=ptrMDS->nGrw;i++){
        adstring name  = ptrMDS->ppGrw[i-1]->name;
        //update growth data for year mxYr+1
        //NOTHING TO UPDATE
        //write new growth data
        adstring nwF = wts::concatenateFilePaths(nxtEstModFolder,ptrMDS->fnsGrowthData[i]);
        ofstream ofs; ofs.open(nwF,ios::trunc);
        ptrMDS->ppGrw[i-1]->write(ofs);
        ofs.close();
    }
    //--chela height data
    for (int i=1;i<=ptrMDS->nCHD;i++){
        adstring name  = ptrMDS->ppCHD[i-1]->name;
        //update chela height data for year mxYr+1
        //NOTHING TO UPDATE
        //write new chela height data
        adstring nwF = wts::concatenateFilePaths(nxtEstModFolder,ptrMDS->fnsChelaHeightData[i]);
        ofstream ofs; ofs.open(nwF,ios::trunc);
        ptrMDS->ppCHD[i-1]->write(ofs);
        ofs.close();
    }
    adstring nwMDS = wts::concatenateFilePaths(nxtEstModFolder,ptrMC->fnMDS);
    rpt::echo<<"writing MDS to '"<<nwMDS<<"'"<<endl;
    ptrMDS->write(nwMDS);
    //write OpMod "state" to file
    {
    adstring nwF = "OpModStateFile_"+str(mxYr+2)+".txt";//writing out for July 1, mxYr+2 (yes, this is correct!)
    ofstream ofs; ofs.open(nwF,ios::trunc);
    writeStateForOpMod(0,ofs);
    ofs.close();
    }
    //SEE if FIXES TAC 
    //--calculate OFL
        {
           cout<<"#----Starting OFL calculations for Operating Model"<<endl;
           ofstream echoOFL_OpMod; echoOFL_OpMod.open("calcOFL_OpMod.final.txt", ios::trunc);//changed name?
           echoOFL_OpMod.precision(12);
           calcOFL_OpMod(1,echoOFL_OpMod);//updates ptrOFLResults
           ptrOFLResults->writeCSVHeader(echoOFL_OpMod); echoOFL_OpMod<<endl;
           ptrOFLResults->writeToCSV(echoOFL_OpMod);     echoOFL_OpMod<<endl;
           echoOFL_OpMod.close();
           cout<<"#----Finished OFL calculations"<<endl;
        }
    //--calculate TAC for upcoming year using harvest control rule
            cout<<"#----Starting TAC calculations"<<endl;
            calcTAC_OpMod(doTAC, value(ptrOFLResults->OFL));
            //int closed calcTAC(doTAC, value(ptrOFLResults->OFL));
            //PRINT2B2("TAC for opMod =", calcTAC(doTAC, value(ptrOFLResults->OFL)))
    //-- performance metrics for Operating Model 
        {
        // MALE AND FEMALE BIOMASS 
            dvector weightsMALE = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE); //male weight matrix
            dvector weightsFEM = ptrMDS->ptrBio->wAtZ_xmz(FEMALE, MATURE);//female weight matrix
            //PRINT2B2("Male weights", weightsMALE)
            //PRINT2B2("Female weights", weightsFEM)
            dvector abundELM(20,32);
            abundELM.initialize();
            double ELMB = 0; // exploitable legal male biomass
            dvector MMA(1,32);
            MMA.initialize();
            dvector MFA(1,32);
            MFA.initialize();
            double MMB=0; // Mature Male biomass 
            double MFB=0; // Mature Female biomass 
            //look over shell condition 1 and 2 (new and old) 
            for (int s = 1; s<=2; s++){
                //for(int z=1;z<=nZBs;z++){
                    abundELM += value(prj_n_xmsz(MALE, MATURE, s)(20,32));
                    MMA +=value(prj_n_xmsz( MALE, MATURE, s));
                    MFA +=value(prj_n_xmsz(FEMALE, MATURE, s));
                   //}//z = number of size bins
                }//s = shell condition
            //Abundance*weight for Exploitable Legal Male Biomass
            for(int i = abundELM.indexmin(); i<=abundELM.indexmax(); ++i){ // check this indexing format     
                ELMB += abundELM[i]*weightsMALE[i];
            }
            //Abundance*weight for Mature Male Biomass
             for(int i = MMA.indexmin(); i<=MMA.indexmax(); ++i){ // check this indexing format     
                MMB += MMA[i]*weightsMALE[i];
            }
            //Abundance*weight for Mature Female Biomass
            for(int i = MFA.indexmin(); i<=MFA.indexmax(); ++i){ // check this indexing format     
                MFB += MFA[i]*weightsFEM[i];
            }
            //testing
            //PRINT2B2("MMA", MMA)
            //PRINT2B2("MFA", MFA)
            //PRINT2B2("ELMA", abundELM)
            // NARROW DOWN CATCH AND DISCARDS TO MALES AND FEMALES MATURE AND IMMATURE BY WEIGHT
              // CATCH 
                PRINT2B1("Calculating CATCH Biomass for OP MODEL")
                // loop over fishery, shells, for sizes
                d5_array vrmN_fxmsz = wts::value(prj_rmN_fxmsz);
                d5_array vrmB_fxmsz = tcsam::calcBiomass(vrmN_fxmsz,ptrMDS->ptrBio->wAtZ_xmz);
                //PRINT2B2("Catch", vrmB_fxmsz)
                PRINT2B1(" YOU MADE IT TO CATCH")
               double MFCB = 0;
               double IFCB = 0;
               dvector ELMC_temp(20,32);
               ELMC_temp.initialize(); 
               double ELMC = 0;
               double MMCB = 0;
               double IMCB = 0;
                    for(int f=1;f<=nFsh;f++){          //Look over all fisheries 
                        for(int s = 1; s<=2; s++){          //look over shell condition 1 and 2 (new and old)  
                            for(int z=1;z<=nZBs;z++){   // Look over all size bins 
                        // Female Catch Biomass
                        MFCB += vrmB_fxmsz(f,FEMALE,MATURE,s,z);
                        IFCB +=vrmB_fxmsz(f,FEMALE,IMMATURE,s,z);
                        //Male Catch Biomass
                        MMCB +=vrmB_fxmsz(f,MALE,MATURE,s,z);
                        IMCB +=vrmB_fxmsz(f,MALE,IMMATURE,s,z);                        
                    } //z
                   }//s
                }//f
                //Exploitable Legal Male Catches 
                     for(int f=1;f<=nFsh;f++){          //Look over all fisheries 
                        for(int s = 1; s<=2; s++){          //look over shell condition 1 and 2 (new and old)  
                            ELMC_temp += vrmB_fxmsz(f,MALE,MATURE,s)(20,32);
                        } // shell
                    } // fishery
                //PRINT2B2("ELMC_vec_OPMod=", ELMC_temp)
                 for(int i = ELMC_temp.indexmin(); i<=ELMC_temp.indexmax(); ++i){ // check this indexing format     
                ELMC += ELMC_temp[i];
            }
                //PRINT2B2("ELMC_OPMod=", ELMC)
                //PRINT2B2("IFCatch", IFCB)
                //PRINT2B2("MMCatch", MMCB)
                //PRINT2B2("IMCatch", IMCB)
                // !!!!! DISCARDS !!!!!!
                    PRINT2B1("Calculating DISCARD Biomass for OP MODEL") //EDIT
                // loop over fishery, shells, for sizes
                d5_array vdmN_fxmsz = wts::value(prj_dmN_fxmsz);
                d5_array vdmB_fxmsz = tcsam::calcBiomass(vdmN_fxmsz,ptrMDS->ptrBio->wAtZ_xmz);
                //PRINT2B2("Discards", vdmB_fxmsz)
                PRINT2B1(" YOU MADE IT TO DISCARD Biomass")
               double MFDB = 0; // Mature Female Discard Biomass
               double IFDB = 0; // Immature Female Dicard Biomass 
               dvector ELMD_temp(20,32); // Exploitable Legal Male Discard Biomass
               ELMD_temp.initialize();
               double ELMD = 0; //
               double MMDB = 0; // Mature Male Discard Biomass 
               double IMDB = 0; // Immature Male Discard Biomass
                    for(int f=1;f<=nFsh;f++){          //Look over all fisheries 
                        for(int s = 1; s<=2; s++){          //look over shell condition 1 and 2 (new and old)  
                            for(int z=1;z<=nZBs;z++){   // Look over all size bins 
                        // Female Discard Biomass
                        MFDB += vdmB_fxmsz(f,FEMALE,MATURE,s,z);
                        IFDB +=vdmB_fxmsz(f,FEMALE,IMMATURE,s,z);
                        //Male Discard Biomass
                        MMDB +=vdmB_fxmsz(f,MALE,MATURE,s,z);
                        IMDB +=vdmB_fxmsz(f,MALE,IMMATURE,s,z);
                    } //z
                   }//s
                }//f
                 //Exploitable Legal Male Discards
                     for(int f=1;f<=nFsh;f++){          //Look over all fisheries 
                        for(int s = 1; s<=2; s++){          //look over shell condition 1 and 2 (new and old)  
                            ELMD_temp += vdmB_fxmsz(f,MALE,MATURE,s)(20,32);
                        } // shell
                    } // fishery
                //PRINT2B2("ELMD_vec_OPMod=", ELMD_temp)
                 for(int i = ELMD_temp.indexmin(); i<=ELMD_temp.indexmax(); ++i){ // check this indexing format     
                ELMD += ELMD_temp[i];
            }
                //PRINT2B2("ELMD_vec_OPMod=", ELMD_temp)
                //PRINT2B2("MFDiscards", MFDB)
                //PRINT2B2("IFDiscards", IFDB)
                //PRINT2B2("MMDiscards", MMDB)
                //PRINT2B2("IMDiscards", IMDB)
              //  PRINT2B2("ELMD_OPMod=", ELMD)
            // RECRUITMENT (1 value)
            double RecAve = mean(ptrOMI->R_y(1981,mxYr))*ptrOMI->R_x(MALE);
            //dvariable test = mean(ptrOMI->R_y(1981,mxYr));
            //PRINT2B2("TestRecAve", test)
            //PRINT2B2("RecAve", RecAve)
            ///////////////////////////////////////
            // Add performance metric for ELM as defined by the State of Alaska with Selectivity of OS animals
            //////////////////////////////////////
            PRINT2B1( " Calculating ELM as defined by the State of AK, with soS (OPMOD) ")
        dvector weights = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE);
        //Calc only newshell mature males 
        double newshell = value(sum(prj_n_xmsz(MALE,MATURE, NEW_SHELL)(20,32))); // 
        //Calculate total mature males 
        double total = 0.0;
        for (int s = 1; s<=nSCs; s++){
             total += value(sum(prj_n_xmsz(MALE, MATURE, s)(20,32))); //WHY SUM
            }
        double propNS = newshell/total;
        dvector abundELM_State(20,32);
        abundELM_State.initialize();
        for (int s = 1; s<=nSCs; s++){
            abundELM_State += value(prj_n_xmsz(MALE, MATURE, s)(20,32)); // WHY NO SUM
            }
        PRINT2B1("#GET ELM")
        dvector ELM_State(20,32);
        ELM_State.initialize();
        ELM_State = (propNS*abundELM_State)+(0.40*(1-propNS)*abundELM_State);
        //PRINT2B2("#ELM=", ELM_State) 
        PRINT2B1("#GET ELMB")
        double ELMB_State = 0;
        //Abundance*weight for Exploitable Legal Male Biomass
            for(int i = ELM_State.indexmin(); i<=ELM_State.indexmax(); ++i){ // check this indexing format     
                ELMB_State += ELM_State[i]*weights[i];
            }
        //PRINT2B2("#ELMB=", ELMB_State) 
            ///////////////////////////////////
            PRINT2B1("writing performance metrics OP MODEL")
            cout<<"#---writing performance metrics OP MODEL "<<ptrMC->mxYr+1<<endl;
            adstring perfMetricsOPMod = "perfMetricsOPMod_"+str(mxYr+1)+".txt"; // Add perfMet somewhere earlier in function?
            ofstream os; os.open("perfMetricsOPMod.txt", ios::trunc);       // ::trunc or app?
            os.precision(6); //number of decimal places in file
            os<<"#--Performance Metrics for Operating Model--"<<endl;
            os<<"perfMetricsOPMod.final=list("<<endl; 
            os<<"TACset=";os<<calcTAC_OpMod(doTAC, value(ptrOFLResults->OFL));os<<cc<<endl;
            //PRINT2B2("TAC = ", repTAC_OpMod(doTAC, value(ptrOFLResults->OFL)))
            //PRINT2B2("TAC input = ", doTAC)
            //PRINT2B2("OFL = ", ptrOFLResults->OFL)
            os<<"TAC=";os<<repTAC_OpMod(doTAC, value(ptrOFLResults->OFL)); os<<cc<<endl; // showing 0 value
            os<<"OFL=";os<<ptrOFLResults->OFL; os<<cc<<endl; // OFL
            os<<"B0=";os<<value(ptrOFLResults->B0);os<<cc<<endl; // B0
            os<<"Bmsy=";os<<value(ptrOFLResults->Bmsy);os<<cc<<endl; // Bmsy
            os<<"Fmsy=";os<<value(ptrOFLResults->Fmsy);os<<cc<<endl; // Fmsy
            os<<"Fofl=";os<<value(ptrOFLResults->Fofl);os<<cc<<endl; // Fofl
            os<<"MMB=";os<<MMB; os<<cc<<endl; // Mature Male Biomass 
            os<<"MFB=";os<<MFB; os<<cc<<endl; // Mature Female Biomass 
            os<<"ELMB=";os<<ELMB; os<<cc<<endl; // 5"< males 
            os<<"ELMB_State=";os<<ELMB_State; os<<cc<<endl; // ELMB state definition
            //os<<"AveRec=";wts::writeToR(os,mean(ptrOMI->R_y(1981,mxYr))*ptrOMI->R_x);os<<cc<<endl;
            //PRINT2B2("aveRec=", mean(ptrOMI->R_y(1981,mxYr))*ptrOMI->R_x)
            // Recruitment (1 value)
            os<<"RecAve=";os<<RecAve; os<<cc<<endl;
            os<<"Rec=";os<<ptrOMI->R_y(mxYr);os<<cc<<endl;
            // Catch Biomass (By sex and maturity)
            os<<"MFCB="; os<<MFCB; os<<cc<<endl;
            os<<"IFCB="; os<<IFCB; os<<cc<<endl;
            os<<"ELMC="; os<<ELMC; os<<cc<<endl;
            os<<"MMCB="; os<<MMCB; os<<cc<<endl;
            os<<"IMCB="; os<<IMCB; os<<cc<<endl;
            // Discard Biomass (By sex and maturity)
            os<<"MFDB="; os<<MFDB; os<<cc<<endl;
            os<<"IFDB="; os<<IFDB; os<<cc<<endl;
            os<<"ELMD="; os<<ELMD; os<<cc<<endl;
            os<<"MMDB="; os<<MMDB; os<<cc<<endl;
            os<<"IMDB="; os<<IMDB; os<<endl;
            //os<<"DISCARDS="; wts::writeToR(os,wts::value(prj_dmN_fxmsz),fDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            //os<<"CATCH="; wts::writeToR(os,wts::value(prj_rmN_fxmsz),fDms,xDms,mDms,sDms,zbDms); os<<endl;
            os<<")"<<endl;
            os<<"#--Finished performance metrics OP MODEL FINAL_PHASE--"<<endl;
            //os<<")"<<endl;
            os.close();
            PRINT2B1("finished writing performance metrics OP MODEL")
        }
}

void model_parameters::final_calcs()
{
    PRINT2B1(" ")
    PRINT2B1("#--Starting FINAL_SECTION")
    if (!mseMode){
        PRINT2B1("#--mseMode is OFF") //not running OP or EST model
        if (mcevalOn) {
            PRINT2B1("#--mceval is ON")
            PRINT2B1("#----Closing mcmc file")
            mcmc.open((char*)(fnMCMC),ios::app);
            mcmc.precision(12);
            //mcmc<<"NULL)"<<endl;
            mcmc.close();
            PRINT2B1(" ")
        }
        if (!mcevalOn){
            PRINT2B1("#--mceval is OFF")
            PRINT2B2("obj fun = ",objFun)
            PRINT2B1(" ")
            {
                PRINT2B1("#----Writing sim data to file")
                ofstream echo1; echo1.open("ModelSimData.dat", ios::trunc);
                echo1.precision(12);
                writeSimData(echo1,0,cout,ptrSimMDS);
                PRINT2B2("obj fun = ",objFun)
                PRINT2B1(" ")
            }
            {
                PRINT2B1("#----Calculating final effective weights for size compositions")
                PRINT2B1("#----Note that the value of objFun is not valid now!!")
                //note that this modifies the value of the objective function!
                ofstream os; os.open("effectiveWeights.R", ios::app);
                os<<"#--Calculating effective weights in FINAL_PHASE--"<<endl;
                os<<"effWgts.final=list("<<endl;
                os<<"surveys=";   calcWeightsForSurveySizeComps( -1,os); os<<","<<endl;
                os<<"fisheries="; calcWeightsForFisherySizeComps(-1,os); os<<endl;
                os<<")"<<endl;
                os<<"#--Finished calculating effective weights in FINAL_PHASE--"<<endl;
                os<<")"<<endl;
                os.close();
                PRINT2B2("#--obj fun = ",objFun)
                PRINT2B1(" ")
            }
            if (doDynB0>0){
                PRINT2B1("#----Calculating dynamic B0")
                PRINT2B1("#----Note that the value of objFun is not valid now!!")
                calcDynB0(1,rpt::echo);
                PRINT2B2("#--obj fun = ",objFun)
                PRINT2B1(" ")
            }
            PRINT2B1("#----Recalculating final objective function value")
            if (!runAlt) runPopDyMod(0,cout); else runAltPopDyMod(0,cout);
            calcObjFun(dbgObjFun,cout);
            PRINT2B2("#--Final obj fun = ",objFun)
            //do OFL calculations
            if (doOFL){
                PRINT2B1("#----Starting OFL calculations")
                ofstream echoOFL; echoOFL.open("calcOFL.final.txt", ios::trunc);
                echoOFL.precision(12);
                calcOFL(mxYr+1,1,echoOFL);//updates ptrOFLResults
                ptrOFLResults->writeCSVHeader(echoOFL); echoOFL<<endl;
                ptrOFLResults->writeToCSV(echoOFL);     echoOFL<<endl;
                echoOFL.close();
                PRINT2B1("#----Finished OFL calculations")
            }
            //do TAC calculations
            PRINT2B2("#----doTAC = ",doTAC)
            int closed = 0;
            if (doTAC>0) {
                PRINT2B1("#----Starting TAC calculations")
                closed = calcTAC(doTAC, value(ptrOFLResults->OFL));
                PRINT2B1("#----Finished TAC calculations")
            }
            //write state for operating model
            {
                PRINT2B1("#----Writing OpMod state file")
                adstring nwF = "OpModStateFile_"+str(mxYr+1)+".txt";
                //nwF = wts::concatenateFilePaths(parent,nwF);
                ofstream ofs; ofs.open(nwF,ios::trunc);
                writeStateForOpMod(mxYr,ofs);
                ofs.close();
                PRINT2B1("#----Finished writing OpMod state file")
            }
            //write MPI for EstMod for upcoming year
            //----devs for upcoming year are zero
            {
                ptrMPI->addNextYearToInfo(closed);
                ptrMPI->setToWriteVectorInitialValues(false);
                adstring nwF = "EstMod.ParametersInfo."+str(mxYr+1)+".inp";
                ofstream ofs; ofs.open(nwF, ios::trunc);
                ptrMPI->write(ofs);
                ofs.close();
            }
            //--write pin for EstMod for upcoming year
            {
                adstring nwF = "EstModPinFile_"+str(mxYr+1)+".txt";
                ofstream ofs; ofs.open(nwF, ios::trunc);
                ptrMPI->writePin(ofs);
                ofs.close();
            }
        }
    } else if (mseEstModMode){
        //running in estModMode
        PRINT2B1("#----MSE EstModMode: Recalculating population dynamics")
        if (!runAlt) runPopDyMod(0,cout); else runAltPopDyMod(0,cout);
        calcObjFun(dbgObjFun,cout);
        PRINT2B2("#--Final obj fun = ",objFun)
        //--calculate OFL
        {
            cout<<"#----Starting OFL calculations"<<endl;
            ofstream echoOFL; echoOFL.open("calcOFL.final.txt", ios::trunc);
            echoOFL.precision(12);
            calcOFL(mxYr+1,1,echoOFL);//updates ptrOFLResults
            ptrOFLResults->writeCSVHeader(echoOFL); echoOFL<<endl;
            ptrOFLResults->writeToCSV(echoOFL);     echoOFL<<endl;
            echoOFL.close();
            cout<<"#----Finished OFL calculations"<<endl;
        }
        //--calculate TAC for upcoming year using harvest control rule
        int closed = calcTAC(doTAC, value(ptrOFLResults->OFL));
        //write MPI for EstMod for upcoming year
        //----devs for upcoming year are zero
        {
            ptrMPI->addNextYearToInfo(closed);
            ptrMPI->setToWriteVectorInitialValues(false);
            adstring nwF = "EstMod.ParametersInfo."+str(mxYr+1)+".inp";
            ofstream ofs; ofs.open(nwF, ios::trunc);
            ptrMPI->write(ofs);
            ofs.close();
        }
        //--write pin for EstMod for upcoming year
        {
            adstring nwF = "EstModPinFile_"+str(mxYr+1)+".txt";
            ofstream ofs; ofs.open(nwF, ios::trunc);
            ptrMPI->writePin(ofs);
            ofs.close();
        }
        // WRITE OFstream for TAC,OFL, biomass performance metrics for ESTIMATION MODEL 
        {
         // IDENTIFY MMB, MFB, and ELMB for Per metrics 
            dvector weightsMALE = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE); //male weight matrix
            dvector weightsFEM = ptrMDS->ptrBio->wAtZ_xmz(FEMALE, MATURE);//female weight matrix 
            dvector abundELM(20,32);
            abundELM.initialize();
            double ELMB = 0; // exploitable legal male biomass
            dvector MMA(1,32);
            MMA.initialize();
            dvector MFA(1,32);
            MFA.initialize();
            double MMB=0; // Mature Male biomass 
            double MFB=0; // Mature Female biomass 
            //look over shell condition 1 and 2 (new and old) 
            for (int s = 1; s<=2; s++){
                //for(int z=1;z<=nZBs;z++){
                    abundELM += value(n_yxmsz(mxYr, MALE, MATURE, s)(20,32));
                    MMA +=value(n_yxmsz(mxYr, MALE, MATURE, s));
                    MFA +=value(n_yxmsz(mxYr, FEMALE, MATURE, s));
                   //}//z
                }//s
            //PRINT2B2("ELMA",abundELM)
            //PRINT2B2("MMA", MMA)
            //PRINT2B2("MFA", MFA)
            //Abundance*weight for Exploitable Legal Male Biomass
            for(int i = abundELM.indexmin(); i<=abundELM.indexmax(); ++i){ // check this indexing format     
                ELMB += abundELM[i]*weightsMALE[i];
            }
            //Abundance*weight for Mature Male Biomass
             for(int i = MMA.indexmin(); i<=MMA.indexmax(); ++i){ // check this indexing format     
                MMB += MMA[i]*weightsMALE[i];
            }
            //Abundance*weight for Mature Female Biomass
            for(int i = MFA.indexmin(); i<=MFA.indexmax(); ++i){ // check this indexing format     
                MFB += MFA[i]*weightsFEM[i];
            }
             // NARROW DOWN CATCH AND DISCARDS TO MALES AND FEMALES MATURE AND IMMATURE BY WEIGHT
                // CATCH 
                PRINT2B1("Calculating CATCH Biomass for EST MODEL")
                // loop over fishery, shells, for sizes
                d6_array vrmN_fyxmsz = wts::value(rmN_fyxmsz);
                d6_array vrmB_fyxmsz = tcsam::calcBiomass(vrmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
                //PRINT2B2("Catch", vrmB_fyxmsz)
                PRINT2B1(" YOU MADE IT TO CATCH")
               double MFCB = 0; // Mature female catch biomass
               double IFCB = 0; // Immature female catch biomass
               dvector ELMC_temp(20,32); // ELM catch biomass vector
               ELMC_temp.initialize();
               double ELMC = 0; // ELM catch biomass
               double MMCB = 0; // Mature male catch biomass
               double IMCB = 0; // Immature male catch biomass
                    for(int f=1;f<=nFsh;f++){          //Look over all fisheries 
                        for(int s = 1; s<=2; s++){          //look over shell condition 1 and 2 (new and old)  
                            for(int z=1;z<=nZBs;z++){   // Look over all size bins 
                        // Female Catch Biomass
                        MFCB += vrmB_fyxmsz(f,mxYr,FEMALE,MATURE,s,z);
                        IFCB +=vrmB_fyxmsz(f,mxYr,FEMALE,IMMATURE,s,z);
                        //Male Catch Biomass
                        MMCB +=vrmB_fyxmsz(f,mxYr,MALE,MATURE,s,z);
                        IMCB +=vrmB_fyxmsz(f,mxYr,MALE,IMMATURE,s,z);
                    } //z
                   }//s
                }//f
                //Exploitable Legal Male Catches 
                     for(int f=1;f<=nFsh;f++){          //Look over all fisheries 
                        for(int s = 1; s<=2; s++){          //look over shell condition 1 and 2 (new and old)  
                            ELMC_temp += vrmB_fyxmsz(f,mxYr,MALE,MATURE,s)(20,32);
                        } // shell
                    } // fishery
                //PRINT2B2("ELMC_vec_EstMod=", ELMC_temp)
                 for(int i = ELMC_temp.indexmin(); i<=ELMC_temp.indexmax(); ++i){ // check this indexing format     
                ELMC += ELMC_temp[i];
            }
             //PRINT2B2("ELMC_EstMod=", ELMC)
                //PRINT2B2("IFCatch", IFCB)
                //PRINT2B2("MMCatch", MMCB)
                //PRINT2B2("IMCatch", IMCB)
                // !!!!! DISCARDS !!!!!!
                    PRINT2B1("Calculating DISCARD Biomass for EST MODEL") //EDIT
                // loop over fishery, shells, for sizes
                d6_array vdmN_fyxmsz = wts::value(dmN_fyxmsz);
                d6_array vdmB_fyxmsz = tcsam::calcBiomass(vdmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
                //PRINT2B2("Discards", vdmB_fyxmsz)
                PRINT2B1(" YOU MADE IT TO DISCARD Biomass")
               double MFDB = 0; // Mature Female Discard Biomass
               double IFDB = 0; // Immature Female Dicard Biomass 
               dvector ELMD_temp(20,32); // ELM Discard Biomass 
               ELMD_temp.initialize();
               double ELMD =0;
               double MMDB = 0; // Mature Male Discard Biomass 
               double IMDB = 0; // Immature Male Discard Biomass
                    for(int f=1;f<=nFsh;f++){          //Look over all fisheries 
                        for(int s = 1; s<=2; s++){          //look over shell condition 1 and 2 (new and old)  
                            for(int z=1;z<=nZBs;z++){   // Look over all size bins 
                        // Female Discard Biomass
                        MFDB += vdmB_fyxmsz(f,mxYr,FEMALE,MATURE,s,z);
                        IFDB +=vdmB_fyxmsz(f,mxYr,FEMALE,IMMATURE,s,z);
                        //Male Discard Biomass
                        MMDB +=vdmB_fyxmsz(f,mxYr,MALE,MATURE,s,z);
                        IMDB +=vdmB_fyxmsz(f,mxYr,MALE,IMMATURE,s,z);
                    } //z                      
                   }//s
                }//f
               // PRINT2B2("MFDiscards", MFDB)
               // PRINT2B2("IFDiscards", IFDB)
               // PRINT2B2("MMDiscards", MMDB)
               // PRINT2B2("IMDiscards", IMDB)
               //Exploitable Legal Male Discards
                     for(int f=1;f<=nFsh;f++){          //Look over all fisheries 
                        for(int s = 1; s<=2; s++){          //look over shell condition 1 and 2 (new and old)  
                            ELMD_temp += vdmB_fyxmsz(f,mxYr,MALE,MATURE,s)(20,32);
                        } // shell
                    } // fishery
                //PRINT2B2("ELMD_vec_EstMod=", ELMD_temp)
                 for(int i = ELMD_temp.indexmin(); i<=ELMD_temp.indexmax(); ++i){ // check this indexing format     
                ELMD += ELMD_temp[i];
            }
            //PRINT2B2("ELMD=", ELMD)
            // RECRUITMENT (1 value)
            prevariable RecAve = mean(R_y(1981,mxYr))*R_yx(mxYr, MALE);
            //PRINT2B2("RecAve", RecAve)
    /////////////////////////////////////////////
   // Exploitable legal males as defined by the State of Alaska with Selectivity of old shell animals
    ////////////////////////////////////////////////
            PRINT2B1( " Calculating ELM as defined by the State of AK, with soS")
            dvector weights = ptrMDS->ptrBio->wAtZ_xmz(MALE, MATURE);
            double newshell = value(sum(n_yxmsz(mxYr,MALE,MATURE, NEW_SHELL)(20,32))); //n_yxmsz is 5d array by year, sex, maturity, stage, size
            double total = 0.0;
            for (int s = 1; s<=nSCs; s++){
                total += value(sum(n_yxmsz(mxYr, MALE, MATURE, s)(20,32))); 
                }
            double propNS = newshell/total;
            dvector abundELM_State(20,32);
            abundELM_State.initialize();
            for (int s = 1; s<=nSCs; s++){
                abundELM_State += value(n_yxmsz(mxYr, MALE, MATURE, s)(20,32)); 
            }
            PRINT2B1("#GET ELM")
            dvector ELM_State(20,32);
            ELM_State.initialize();
            ELM_State = (propNS*abundELM)+(0.40*(1-propNS)*abundELM);
            //PRINT2B2("#ELM=", ELM_State) 
            PRINT2B1("#GET ELMB")
            double ELMB_State = 0;
        //Abundance*weight for Exploitable Legal Male Biomass
            for(int i = ELM_State.indexmin(); i<=ELM_State.indexmax(); ++i){ // check this indexing format     
                ELMB_State += ELM_State[i]*weights[i];
            }
            //PRINT2B2("#ELMB=", ELMB_State) 
        ////////////////////////////////////
            PRINT2B1("writing performance metrics EST MODEL") //EDIT
            adstring perfMetricsEstMod = "perfMetricsEstMod_"+str(mxYr+1)+".txt"; // Add perfMet somewhere earlier in function?
            ofstream os; os.open(perfMetricsEstMod, ios::trunc);       // ::trunc or app? Eventually appended but FIX IT FIRST
            os<<"#--Performance Metrics for Estimation Model--1--"<<endl;
            os<<"perfMetricsEstMod.final=list("<<endl;                  //Check line 
            os<<"TACset=";os<<calcTAC(doTAC,value(ptrOFLResults->OFL));os<<cc<<endl;
            os<<"TAC=";os<<repTAC(doTAC,value(ptrOFLResults->OFL));os<<cc<<endl;
            os<<"OFL=";os<<value(ptrOFLResults->OFL);os<<cc<<endl; // OFL
            os<<"B0=";os<<value(ptrOFLResults->B0);os<<cc<<endl; // B0
            os<<"Bmsy=";os<<value(ptrOFLResults->Bmsy);os<<cc<<endl; // Bmsy
            os<<"Fmsy=";os<<value(ptrOFLResults->Fmsy);os<<cc<<endl; // Bmsy
            os<<"Fofl=";os<<value(ptrOFLResults->Fofl);os<<cc<<endl; // Bmsy
            os<<"MMB="; os<<MMB; os<<cc<<endl; //
            os<<"MFB="; os<<MFB; os<<cc<<endl; // 
            os<<"ELMB=";os<<ELMB;os<<cc<<endl;
            os<<"ELMB_State=";os<<ELMB_State;os<<cc<<endl;
            os<<"AveRec=";os<<RecAve;os<<cc<<endl;
            //os<<"AveRec_test="; wts::writeToR(os,value(mean(R_y(1981,mxYr))*R_yx(mxYr)));os<<cc<<endl;
            //PRINT2B2("aveRec=", mean(R_y(1981,mxYr))*R_yx(mxYr))
            os<<"Rec=";os<<R_y(mxYr);os<<cc<<endl;
            // Catch Biomass (By sex and maturity)
            os<<"MFCB="; os<<MFCB; os<<cc<<endl;
            os<<"IFCB="; os<<IFCB; os<<cc<<endl;
            os<<"ELMC="; os<<ELMC; os<<cc<<endl;
            os<<"MMCB="; os<<MMCB; os<<cc<<endl;
            os<<"IMCB="; os<<IMCB; os<<cc<<endl;
            // Discard Biomass (By sex and maturity)
            os<<"MFDB="; os<<MFDB; os<<cc<<endl;
            os<<"IFDB="; os<<IFDB; os<<cc<<endl;
            os<<"ELMD="; os<<ELMD; os<<cc<<endl;
            os<<"MMDB="; os<<MMDB; os<<cc<<endl;
            os<<"IMDB="; os<<IMDB; os<<endl;
            //os<<"DISCARDS="; wts::writeToR(os,wts::value(dmN_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            //os<<"CATCH="; wts::writeToR(os,wts::value(rmN_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<endl;
            os<<")"<<endl;
            os<<"#--Finished TAC,OFL, biomass performance metrics EST MODEL FINAL_PHASE--"<<endl;
            //os<<")"<<endl;
            os.close();
            PRINT2B1("finished writing TAC,OFL, biomass, performance metrics EST MODEL")
        }
    } else if (mseOpModMode){
        finishOpModMode();
    }//mseOpModMode
    int hour,minute,second;
    double elapsed_time;
    time(&finish); 
    elapsed_time = difftime(finish,start);
    hour   = (int) long(elapsed_time)/3600;
    minute = (int) long(elapsed_time)%3600/60;
    second = (int) (long(elapsed_time)%3600)%60;
    PRINT2B1("")
    PRINT2B1("#------------------")
    PRINT2B2("Starting time : ",ctime(&start))
    PRINT2B2("Finishing time: ",ctime(&finish))
    adstring hms = "This run took: "+str(hour)+" hours, "+str(minute)+" minutes, "+str(second)+" seconds.";
    PRINT2B1(hms)
    PRINT2B1("#--Finished FINAL_SECTION")
    PRINT2B1("#------------------------")
    PRINT2B1("#------------------------")
    PRINT2B1("")
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{5000,5000,5000,5000,5000,5000,10000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{0.5,0.1,.01,.001,1e-4,1e-5,1e-6}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 2147000000; //must be smaller than 2,147,483,647
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(40000000); // this may be incorrect in the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1500000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(7000);
  time(&start);
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint = defaults::iprint;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
