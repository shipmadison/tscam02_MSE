/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HarvestStrategies.hpp
 * Author: william.stockhausen
 *
 * Created on September 26, 2018, 10:35 AM
 * 2018-9-28 Updated Directories in NetBeans Project
 * 2018-10-1 Meeting with Andre   
 * 2018-10-8 initial changes to header file 
 * 2018-10-15 Meeting with Hayden Corbin, first rule tested
 */


#ifndef HARVESTSTRATEGIES_HPP
#define HARVESTSTRATEGIES_HPP

class HarvestStrategies {
    public:
        static int debug;
        static const double MMB_AveMilLb; //dummy values for testing
        static const double MFB_AveMilLb; // dummy values for testing 
        static const double z; // dummy values for testing, number of ELM size classes 
        
    public:
        /**
         * HCR 1: female ramp.
         * 
         * @param MFB - Mature female biomass for the year 
         * @param aveMFB - averaged MFB from 1982-2016
         * @param MMB - Mature Male Biomass for the Year 
         * 
         * @return TAC = exploitation rate time Mature Male Biomass 
         */
        double static HCR1_FemaleRamp(double MFB, double aveMFB, double MMB);
        
        /**
         * HCR 2. 3 ramps.
         * 
         * @param MMB Mature Male Biomass in millions of lbs for the year 
         * @param aveMMB constant, MMB average from 1982-2016
         * @param rampID identifies the level of exploitation rate, low, moderate, or aggressive 
         * 
         * @return 
         */
        double static HCR2_MaleRamp(double MMB, double aveMMB, int rampID);
        
        /**
         * HCR 22. Survey Estimates of MMB, Slope set to Middle Ramp (10/75)
         * 
         * @param MMB Mature Male Biomass in millions of lbs for the year 
         * @param aveMMB constant, MMB average from 1982-2016
         * @param rampID identifies the level of exploitation rate, low, moderate, or aggressive 
         * 
         * @return 
         */
        
        double static HCR22_MaleRamp_SurvEst(double MMB, double aveMMB);
        
         /**
         * HCR 23. Model Survey Estimates of MMB, Slope set to Middle Ramp (10/75)
         * 
         * @param MMB Mature Male Biomass in millions of lbs for the year 
         * @param aveMMB constant, MMB average from 1982-2016
         * @param rampID identifies the level of exploitation rate, low, moderate, or aggressive 
         * 
         * @return 
         */
        
        double static HCR23_MaleRamp_ModSurvEst(double MMB, double aveMMB);
        
        /**
         * HCR 3. Takes OFL from operating model, uses pre-determined buffer.
         * 
         * @param OFL from the model
         * @param buffer <- current buffer on OFL (20% currently?)
         * 
         * @return 
         */
        double static HCR3_ABC(double OFL, double buffer);
        
        /**
         * HCR 4. Scaled with females--> moving from 5% to 20% exploitation. 
         * 
         * @param MFB
         * @param aveMFB
         * @param MMB
         * @param aveMMB
         * @return 
         */
        double static HCR4_FemaleDimmer(double MFB, double aveMFB, double MMB, double aveMMB);
        
        /**
         * HCR 41. Scaled with females--> moving from 10% to 20% exploitation. 
         * 
         * @param MFB
         * @param aveMFB
         * @param MMB
         * @param aveMMB
         * @return 
         */
        double static HCR41_FemaleDimmer(double MFB, double aveMFB, double MMB, double aveMMB);
        
        /**
         * HCR 42. Scaled with females--> moving from 10% to 22.5% exploitation. 
         * 
         * @param MFB
         * @param aveMFB
         * @param MMB
         * @param aveMMB
         * @return 
         */
        double static HCR42_FemaleDimmer(double MFB, double aveMFB, double MMB, double aveMMB);
        
        /**
         * HCR 43 --> Written in for if statement for Max TAC at 30% ELM
         * 
         * @param MFB
         * @param aveMFB
         * @param MMB
         * @param aveMMB
         * @return 
         */
        double static HCR43_FemaleDimmer(double MFB, double aveMFB, double MMB, double aveMMB);
        
        /**
         * HCR 5. Scaled with females --> blocks.
         * 
         * @param MFB
         * @param aveMFB
         * @param MMB
         * @param aveMMB
         * 
         * @return 
         */
        double static HCR5_FemaleBlocks(double MFB, double aveMFB, double MMB, double aveMMB);
    
        /**
         * HCR 6. ELM. 
         * 
         * @param propNS proportion new shell
         * @param abunELM
         * @param weights
         * @param sOS
         * @param xpRate 30%, 40%, or 50%
         * 
         * @return 
         */
        double static HCR6_ELM(double propNS, dvector abunELM, dvector weights, double sOS, double xpRate );
        
        /**
         * HCR 7. Status Quo
         * @param MFB
         * @param aveMFB
         * @param MMB
         * @param aveMMB
         * @param abunELM
         * @param weights
         * @param FISHERY_SELECTIVITY Fix this, find the selectivity(sex,size) in tpl
         * @return 
         */
        
        double static HCR7_StatusQuo(double MFB, double aveMFB, double MMB, double aveMMB, double CWmsy);


private:
        //class constructor
        HarvestStrategies(){};
        //class destructor
       ~HarvestStrategies(){};
};



#endif /* HARVESTSTRATEGIES_HPP */

