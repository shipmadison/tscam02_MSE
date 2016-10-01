/* 
 * File:   ModelData.hpp
 * Author: william.stockhausen
 *
 * Created on 2014-02-11.
 * 
 * 20140930: 1. Added recLag to BioData.
 * 20150209: 1. Removed prMature_xz, frMature_xsz, cvMnMxZ_xc from BioData.
 * 20150121: 1. Combined FisheryData, SurveyData classes into a FleetData class
 *           2. Updated some documentation.
 */

#ifndef MODELDATA_HPP
#define MODELDATA_HPP

//**********************************************************************
//  Includes
//      AggregateCatchData
//      EffortData
//      SizeFrequencyData
//      CatchData
//      BioData
//      FleetData
//      ModelDatasets
//**********************************************************************
class ModelConfiguration; //forward definition
class IndexRange;

//--------------------------------------------------------------------------------
//          AggregateCatchData
// Change notes:
//  2014-12-08: 1. changed input data row format for so columns are now
//                  year, sex, maturity, shell_condition, value, cv
//              2. changed inpC_yc dimensions from yc to yxmsc to match new input format
//              3. changed C_xy, cv_xy, sd_xy dimensions from xy to xmsy (and renamed them accordingly)
//
//--------------------------------------------------------------------------------
    class AggregateCatchData {
    public:
        static int debug;
        const static adstring KW_ABUNDANCE_DATA;//keyword indicating abundance (numbers) data
        const static adstring KW_BIOMASS_DATA;  //keyword indicating biomass (weight) data
    protected:
        d5_array inpC_xmsyc; //input aggregate catch data (value, cv by xmsy)
    public:
        adstring type;  //type (abundance, biomass) of data (matches one of the keywords above)
        int optFit;     //objective function fitting option
        int llType;     //likelihood function type
        int ny;         //number of data rows for aggregate catch
        adstring units; //units for aggregate catch data
        ivector yrs;    //years for aggregate catch data
        wts::adstring_matrix factors;//factor combinations for input numbers-at-size
        d4_array C_xmsy;   //aggregate catch by sex, maturity, shell_condition, year (converted from units to THOUSANDS of crab or MT))
        d4_array cv_xmsy;  //aggregate catch cv's by sex, maturity, shell_condition, year
        d4_array sd_xmsy;  //aggregate catch stdv's by sex, maturity, shell_condition, year
        d4_array nlls_xmsy;//negative log-likelihood components
        
    public:
        AggregateCatchData(){}
        ~AggregateCatchData(){}
        /**
         * Replace catch data C_xmsy with new data. 
         * Also modifies inpC_xmsyc to reflect new data.
         * Error-related quantities remain the same.
         * 
         * @param iSeed - flag (!=0) to add random noise
         * @param rng - random number generator
         * @param newC_xmsy - d4_array with new catch data
         */
        void replaceCatchData(int iSeed,random_number_generator& rng,d4_array& newC_xmsy);
        void read(cifstream & is);//read file in ADMB format
        void write(ostream & os); //write object to file in ADMB format
        void writeToR(ostream& os, std::string nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, AggregateCatchData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   AggregateCatchData & obj){obj.write(os); return os;}
    
    protected:
        /**
         * Aggregate catch data C_xmsy, cv_xmsy, sd_xmsy over summary indices.
         */
        void aggregateData(void);

};

//--------------------------------------------------------------------------------
//          EffortData
//--------------------------------------------------------------------------------
    class EffortData {
    public:
        /* flag to print debugging info */
        static int debug;
        /* keyword indicating effort data */
        const static adstring KW_EFFORT_DATA;
    public:
        /* number of years of effort data */
        int ny;
        /* interval to average effort/fishing mortality */
        IndexRange* ptrAvgIR; 
        /* units for potlifts */
        adstring units;   
        /* input effort data (year)x(year,potlifts) */
        dmatrix  inpEff_yc;   
        /* vector of years w/ effort data */
        dvector yrs;
        /* effort data vector corresponding to \code{yrs} */
        dvector eff_y;        
    public:
        /**
         * Constructor.
         */
        EffortData(){ptrAvgIR=0;}
        /**
         * Destructor.
         */
        ~EffortData();
        /**
         * Read input data in ADMB format from a file stream
         * 
         * @param is - input file stream
         */
        void read(cifstream & is);//read file in ADMB format
        /**
         * Write data to an output stream in ADMB format
         * 
         * @param os output stream
         */
        void write(ostream & os); //write object to file in ADMB format
        /**
         * Write data to an output stream as an R-formatted list object
         * 
         * @param os output stream
         */
        void writeToR(ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into an EffortData object.
         */
        friend cifstream& operator >>(cifstream & is, EffortData & obj){obj.read(is); return is;}
        /**
         * Operators to write data to an output stream in ADMB format from an EffortData object.
         */
        friend ostream&   operator <<(ostream & os,   EffortData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//          SizeFrequencyData
//--------------------------------------------------------------------------------
    class SizeFrequencyData {
    public:
        static int debug;
        /* keyword indicating size frequency data */
        const static adstring KW_SIZEFREQUENCY_DATA;
    private:
        /* factor combinations for input numbers-at-size */
        wts::adstring_matrix factors;
        /* input numbers-at-size data (sex,maturity state,shell condition,year,year+sample_size+nAtZ) */
        d5_array inpNatZ_xmsyc;       
    public:
        /* objective function fitting option */
        int optFit; 
        /* likelihood function type */
        int llType; 
        
        /* number of size bin cut pts */
        int nZCs;    
        /* vector of cut points for size bins (1:nZCs) */
        dvector zCs; 
        /* vector of size bin centers (1:nZCs-1) */
        dvector zBs; 
        
        /* number of years of size frequency data */
        int ny;         
        /* units for numbers-at-size data */
        adstring units; 
        /* years of size frequency data */
        ivector  yrs;       
        /* sample sizes for size frequency data */
        d4_array ss_xmsy;   
        /* raw size frequency data */
        d5_array NatZ_xmsyz;
        /* normalized size frequency data (sums to 1 over xmsz for each y) */
        d5_array PatZ_xmsyz;
    public:
        /**
         * Constructor.
         */
        SizeFrequencyData(){}
        /**
         * Destructor.
         */
        ~SizeFrequencyData(){}
        /**
         * Replace catch-at-size data NatZ_xmsyz with new data. 
         * Also modifies inpNatZ_xmsyc to reflect new data.
         * Error-related quantities remain the same.
         * 
         * @param iSeed - flag (!=0) to add random noise
         * @param rng - random number generator
         * @param newNatZ_yxmsz - d5_array with new numbers-at-size data
         */
        void replaceSizeFrequencyData(int iSeed,random_number_generator& rng,d5_array& newNatZ_xmsyz);
        /**
         * Save the negative log-likelihoods from a model fit (values only).
         * 
         * @param nlls
         */
        void saveNLLs(dvar_matrix& nlls);
        /**
         * Read input data in ADMB format from a file stream
         * 
         * @param is - input file stream
         */
        void read(cifstream & is);//read file in ADMB format
        /**
         * Write data to an output stream in ADMB format
         * 
         * @param os output stream
         */
        void write(ostream & os); //write object to file in ADMB format
        /**
         * Write data to an output stream as an R-formatted list object
         * 
         * @param os output stream
         */
        void writeToR(ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into an SizeFrequencyData object.
         */
        friend cifstream& operator >>(cifstream & is, SizeFrequencyData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   SizeFrequencyData & obj){obj.write(os); return os;}
    public:
        /**
         * Calculate normalized size compositions PatZ_xmsyz based on NatZ_xmsyz.
         */
        void normalize(void);
    };

//--------------------------------------------------------------------------------
//          BioData
//--------------------------------------------------------------------------------
    class BioData {
    public:
        static int debug;
        const static adstring KW_BIO_DATA;
    public:
        int nZBins;           //number of size bin cut pts
        dvector zBins;        //size bins
        adstring unitsWatZ;   //units for weight-at-size
        d3_array wAtZ_xmz;    //weight at size (in kg)
//        dmatrix prMature_xz;  //probability of maturity at size by sex
//        d3_array frMature_xsz;//fraction mature at size by sex, shell condition
//        dmatrix cvMnMxZ_xc;   //cv of min size, max size by sex
        
        int recLag;         //recruitment lag (in years)
        double fshTimingTypical;//typical midpoint of fishery seasons
        double matTimingTypical;//typical timing of mating
        int nyAtypicalT;        //number of years for atypical fishery season midpoints, time-at-mating
        dvector fshTiming_y;    //timing of midpoint of fisheries seasons
        dvector matTiming_y;    //timing of mating
    protected:
        dmatrix timing_yc;//midpoint of fishery season, time-at-mating by year
    public:
        BioData(){}
        ~BioData(){}
        void read(cifstream & is);//read file in ADMB format
        void write(std::ostream & os); //write object to file in ADMB format
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into a BioData object.
         */
        friend cifstream& operator >>(cifstream & is, BioData & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   BioData & obj){obj.write(os); return os;}
    };

//------------------------------------------------------------------------------
    /**
     * CatchData class.
     */
    class CatchData {
    public:
        static int debug;
        const static adstring KW_CATCH_DATA;
    protected:
        dmatrix inpN_yc;       //input catch abundance data (year,female abundance,female cv,male abundance, male cv, total abundance, cv_total)
        dmatrix inpB_yc;       //input catch biomass   data (year,female biomass,female cv,male biomass, male cv, total biomass, cv_total)
        //TODO: need sample sizes (tows/potlifts, non-zero tows/potlifts, num individ.s by xms, calculated ss)
        d5_array inpNatZ_xmsyc;//input numbers-at-size data (sex,maturity,shell,year) x (year,sample_size,nAtZ)
    public:
        adstring type;    //data type (e.g. survey or fishery)
        adstring name;    //data source name
        int hasN;         //abundance data flag
        AggregateCatchData* ptrN;//pointer to aggregate abundance data
        int hasB;         //biomass data flag
        AggregateCatchData* ptrB;//pointer to aggregate biomass data
        int hasZFD;       //numbers-at-size data flag
        SizeFrequencyData* ptrZFD;//pointer to numbers-at-size data
        
    public:
        /**
         * Constructor.
         */
        CatchData();
        /**
         * Destructor
         */
        ~CatchData();
        /**
         * Replace catch data based on newNatZ_yxmsz.
         * 
         * @param iSeed - flag (!=0) to add random noise
         * @param rng - random number generator
         * @param newNatZ_yxmsz - d5_array of catch-at-size by sex/maturity/shell condition/year
         * @param wAtZ_xmz - d3_arrray of weight-at-size by sex/maturity
         */
        virtual void replaceCatchData(int iSeed,random_number_generator& rng,d5_array& newNatZ_yxmsz, d3_array& wAtZ_xmz);
        virtual void read(cifstream & is);//read file in ADMB format
        virtual void write(ostream & os); //write object to file in ADMB format
        virtual void writeToR(ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into a CatchData object.
         */
        friend cifstream& operator >>(cifstream & is, CatchData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   CatchData & obj){obj.write(os); return os;}
    };

//------------------------------------------------------------------------------    
    /**
     * Fleet data class.
     */
    class FleetData {
    public:
        static int debug;
        const static adstring KW_FISHERY;
        const static adstring KW_SURVEY;
    public:
        adstring name;     //fleet name
        adstring type;     //fleet type (survey or fishery)
        int hasICD;        //flag indicating observed index (survey) catch data
        CatchData* ptrICD; //pointer to CatchData object for index survey) catch data
        int hasRCD;        //flag indicating retained catch data
        CatchData* ptrRCD; //pointer to CatchData object for retained catch
        int hasDCD;        //flag indicating observed discard catch data
        CatchData* ptrDCD; //pointer to CatchData object for observed discards
        int hasTCD;        //flag indicating observed total catch data
        CatchData* ptrTCD; //pointer to CatchData object for observed total catch
        int hasEff;        //flag indicating effort data
        EffortData* ptrEff;//pointer to effort data
    public:
        /**
         * Constructor.
         */
        FleetData();
        /**
         * Destructor.
         */
        ~FleetData();
        /**
         * Replace existing catch data with new values.
         * 
         * @param iSeed - flag (!=0) to add random noise
         * @param rng - random number generator
         * @param newNatZ_yxmsz - catch data array
         * @param wAtZ_xmz - weight-at-size array
         */
        void replaceIndexCatchData(int iSeed,random_number_generator& rng,d5_array& newNatZ_yxmsz, d3_array& wAtZ_xmz);
        /**
         * Replace existing fishery catch (retained, discarded, total) data with new values.
         * 
         * @param iSeed - flag (!=0) to add random noise
         * @param rng - random number generator
         * @param newCatZ_yxmsz - total catch data array
         * @param newRatZ_yxmsz - retained catch data array
         * @param wAtZ_xmz      - weight-at-size array
         */
        void replaceFisheryCatchData(int iSeed,random_number_generator& rng,d5_array& newCatZ_yxmsz,d5_array& newRatZ_yxmsz,d3_array& wAtZ_xmz);
        /**
         * Read input file in ADMB format from input stream.
         * 
         * @param is - input stream
         */
        void read(cifstream & is);//read file in ADMB format
        /**
         * Write the fleet data in ADMB format (i.e., as an input file)
         * 
         * @param os - output stream to write to
         */
        void write(ostream & os); //write object to file in ADMB format
        void writeToR(ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into a FleetData object.
         */
        friend cifstream& operator >>(cifstream & is, FleetData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   FleetData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//         ModelDatasets
//--------------------------------------------------------------------------------
    class ModelDatasets {
    public:
        static int debug;
    public:
        ModelConfiguration* pMC;//pointer to ModelConfiguration object
        adstring fnBioData;//bio data file name
        BioData* ptrBio;   //pointer to bio dataset object
        
        int nFsh;//number of fishery datasets to read
        adstring_array fnsFisheryData;//fishery data file names       
        FleetData**    ppFsh;         //pointer to array of pointers to fishery dataset objects
        
        int nSrv;//number of survey datasets to read
        adstring_array fnsSurveyData;//survey data files names
        FleetData**    ppSrv;        //pointer to array of pointers to survey dataset objects
    public:
        ModelDatasets(ModelConfiguration* ptrMC);
        ~ModelDatasets();
        void read(cifstream & is);//read file in ADMB format
        void write(std::ostream & os); //write object to file in ADMB format
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into a ModelDatasets object.
         */
        friend cifstream& operator >>(cifstream & is, ModelDatasets & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   ModelDatasets & obj){obj.write(os); return os;}
    };
#endif  /* MODELDATA_HPP */

