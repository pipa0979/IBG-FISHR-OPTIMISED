//#include "ErrorFinderManager.hpp"
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <algorithm>
#include "ErrorFinderManager.hpp"
#include "ErrorCalculator.hpp"

#define GetCurrentDir getcwd




ErrorFinderManager::ErrorFinderManager():
					WINDOW(50),
                    MIN_SNP(30),GAP(1), MA_SNP_END(50), TRUESNP( 600 ),
                    MIN_CM(0.4), MA_ERR_THRESHOLD_START(0.08),
                    MA_ERR_THRESHOLD_END(0.08),PCT_ERR_THRESHOLD( 0.90 ),
                    HO_THRESHOLD( 0.98 ), TRUECM( 6 ),PIELENGTH( 3 ),
                    ISMOL( false), COUNTGAPERR( false ),MA_THRESHOLD(0.8),EMPIRICAL_MA_RESULT(-1.0), EMPIRICAL_PIE_RESULT(-1.0),EXTENDSNP(0),T_ERROR(false),T_TRIM(false),CERRMESSAGE(""),SITUATION_NO(0)
{
}

template <typename T>
std::string NumberToString ( T Number )
{
	std::ostringstream ss;
	ss << Number;
	return ss.str();
}


void ErrorFinderManager::performConsolidation(int argc,char *argv[])
{
   bool goodParam=true;
   bool thresholdError = false; //if the user supplies both -empricial-ma-threshold and -ma-threshold, that is an error. This will be used to detect such an error.
   bool pieThresholdError = false;
   char currentPath[FILENAME_MAX];
   if(!GetCurrentDir(currentPath,sizeof(currentPath)))
	   {
			std::cerr << "Error reading current directory" << std::endl;
			return;
	   }

   for(int i=1;i<argc;i++)
		{

			/*Code for expanding window*/
			if (strcmp(argv[i],"-extendSNP")==0 && i<argc-1)
				{
					EXTENDSNP=atoi(argv[++i]);
					//cout<<"entered snp = "<<EXTENDSNP<<endl;
				}

			else if(strcmp(argv[i],"-bmatch")==0&&i<argc-1)
				{
					BMATCHFILE=std::string(argv[++i]);
				}
			else if(strcmp(argv[i],"-bmid")==0&&i<argc-1)
				{
					BMIDFILE=std::string(argv[++i]);
				}
			else  if(strcmp(argv[i],"-bsid")==0&&i<argc-1)
				{
					BSIDFILE=std::string(argv[++i]);
				}
			 else if(strcmp(argv[i],"-reduced")==0&&i<argc-2)
				{
					MIN_SNP=atoi(argv[++i]);
					MIN_CM=atof(argv[++i]);
				}
			else  if(strcmp(argv[i],"-ped-file")==0&&i<argc-1)
				{
					PEDFILE=std::string(argv[++i]);
				}
				else  if(strcmp(argv[i],"-holdout-ped")==0&&i<argc-1)
				{
					HPEDFILE=std::string(argv[++i]);
				}
				else  if(strcmp(argv[i],"-holdout-map")==0&&i<argc-1)
				{
					HMAPFILE=std::string(argv[++i]);
				}
			else  if(strcmp(argv[i],"-window")==0&&i<argc-1)
				{
					WINDOW=atoi(argv[++i]);
				}
			else  if(strcmp(argv[i],"-holdout-threshold")==0&&i<argc-1)
				{
					HO_THRESHOLD=atof(argv[++i]);
				}
			else  if(strcmp(argv[i],"-trueCM")==0&&i<argc-1)
				{
					TRUECM=atof(argv[++i]);
				}
			else  if( strcmp( argv[i],"-trueSNP" )==0 && i < argc-1 )
				{
					TRUESNP=atoi(argv[++i]);
				}
			else  if(strcmp(argv[i],"-holdout-missing")==0&&i<argc-1)
				{
					HO_MISSING= std::string(argv[++i]);
				}
			else  if(strcmp(argv[i],"-gap")==0&&i<argc-1)
				{
					GAP=atoi(argv[++i]);
				}
			else  if(strcmp(argv[i],"-ma-snp")==0&&i<argc-1)
				{
					MA_SNP_END=atoi(argv[++i]);
				}
			else  if(strcmp(argv[i],"-pct-err-threshold")==0&&i<argc-1)
				{
					 if(pieThresholdError == true)
						{
							std::cerr << "ERROR: You have supplied both -emp-pie-threshold and -pct-err-threshold parameters, but only one is allowed. Please try again." << std::endl;
							exit(1);
						}
					 PCT_ERR_THRESHOLD=atof(argv[++i]);
					 pieThresholdError = true;
				}
			else  if(strcmp(argv[i],"-emp-pie-threshold")==0&&i<argc-1)
				{
					 if(pieThresholdError == true)
						 {
							 std::cerr << "ERROR: You have supplied both -emp-pie-threshold and -pct-err-threshold parameters, but only one is allowed. Please try again." << std::endl;
							 exit(1);
						 }
					 EMPIRICAL_PIE_RESULT=atof(argv[++i]);
					 pieThresholdError = true;
				}

			else  if(strcmp(argv[i],"-output.type")==0&&i<argc-1)
				{
					 OPTION=std::string(argv[++i]);
					 //cout<<OPTION<<endl;
				}
			else  if(strcmp(argv[i], "-snpfile") ==0&&i<argc-1)
				{
					 SNPWEIGHTFILE = std::string(argv[++i]);
				}
			else  if(strcmp(argv[i],"-log.file")==0&&i<argc-1)
				{
					 LOGFILE=std::string(argv[++i]);
				}
			else if(strcmp(argv[i],"-ma-threshold")==0&&i<argc-1)//adding new -ma-threshold argument
				{
					if(thresholdError == true)
						{ //user has already supplied an empirical-ma-threshold, so exit the program with an error message
							std::cerr << "ERROR: You have supplied both -empirical-ma-threshold and -ma-threshold parameters, but only one is allowed. Exiting program."<< std::endl;
							exit(1);
						}
					MA_THRESHOLD=atof(argv[++i]);
					thresholdError = true;
				}
			else if(strcmp(argv[i],"-empirical-ma-threshold")==0 && i<argc-1)
				{
						if(thresholdError == true)
							{
								std::cerr << "ERROR: You have supplied both -empirical-ma-threshold and -ma-threshold parameters, but only one is allowed. Exiting program."<< std::endl;
								exit(1);
							}
						EMPIRICAL_MA_RESULT = atof(argv[++i]); //use the user supplied empirical ma threshold, instead of calculating it via true ibd segments
						thresholdError = true;
				}
			else  if(strcmp(argv[i],"-PIE.dist.length")==0&&i<argc-1)
				{
					std::string MOL=std::string(argv[++i]);
					if( MOL.compare( "MOL" ) ==0 )
						{
							ISMOL = true;
						}
					else
						{
							PIELENGTH = atof( MOL.c_str() );//redundant string oeration
						}
				}
			else  if(strcmp(argv[i],"-count.gap.errors")==0&&i<argc-1)
				{
					std::string option=std::string(argv[++i]);
					if( option.compare( "TRUE" ) ==0 )
						 {
							COUNTGAPERR = true;
						 }
				}

			else
				{
					wrongParam += " " + std::string(argv[i]);
					goodParam=false;
				}
		}	//for loop ends


	if((!goodParam)||BMATCHFILE.compare("")==0||BSIDFILE.compare("")==0||BMIDFILE.compare("")==0||PEDFILE.compare("")==0)
		{
			displayError( argv[0] );
			return;
		}
	if( OPTION.compare( "" ) == 0 )
		{
			std::cerr<< " please provide a valid output.type option " <<std::endl;
			exit( -1 );
		}
	if( LOGFILE.compare( "" ) == 0 )
		{
			std::cerr<< " default log file name is FISH " <<std::endl;
			LOGFILE = "FISH";
		}
	eCalculator.createLogFile( LOGFILE  );
	eCalculator.countGapErrors( COUNTGAPERR );
	time_t startTime;
   // cout<<"helloe5"<<endl;
	time (&startTime);
	std::string str_head = "****************************************************\n";
	str_head += "****************************************************\n";
	str_head += "FISHR LOG FILE INFORMATION\n\n";
	std::string str1 = " The program started at: " + std::string( ctime ( &startTime ) );
	std::string str = str_head + "Program working directory was: " + std::string(currentPath) +
		 " \nProgram version was: " + std::string(argv[0]) +
		 " \nProgram options:\n-bmatch file: " + BMATCHFILE +
		 " \n-bmid file: " + BMIDFILE +
		 " \n-bsid file: " + BSIDFILE +
		 " \n-ped file: " + PEDFILE;
	if(HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) !=0)
		{
			str +=
			" \n-holdout ped file: " + HPEDFILE +
			" \n-holdoutmap file: " + HMAPFILE;
		}
	str += " \n-output type: " + OPTION +
	" \n- missing SNP representation in pedfile: " + HO_MISSING +
	" \n-log file: " + LOGFILE;
	str = str + " \nmin snp length : " + NumberToString( MIN_SNP )  +
				" \nmin cm length : " + NumberToString( MIN_CM  ) +
				" \ngap to consolidate : " + NumberToString( GAP ) +
				" \nmoving averages window size : " + NumberToString(  WINDOW ) +
				" \ndiscard ends to calculate pct err : " + NumberToString( MA_SNP_END );

	if(EMPIRICAL_PIE_RESULT < 0.0)
		{
			  str += " \nuser supplied percentage error threshold: " + NumberToString( PCT_ERR_THRESHOLD ) +
			  "\nuser did not supply an empirical-pie-threshold: NA";
		}
	else
		{
			  str += "\nuser did not supply a percentage error threshold: NA";
			  str += "\nuser supplied empirical-pie-threshold: " + NumberToString( EMPIRICAL_PIE_RESULT );
		}
	if(EMPIRICAL_MA_RESULT < 0.0)
		{
			  str += "\nuser supplied ma threshold: " + NumberToString(MA_THRESHOLD) +
			  "\nuser did not supply an empirical-ma-threshold: NA";
		}
	else
		{
			  str += "\nuser did not supply an ma-threshold: NA";
			  str += "\nuser supplied empirical-ma-threshold: " + NumberToString(EMPIRICAL_MA_RESULT);
		}
	//"\n****************************************************\n";

	eCalculator.log( str );


	initiateErrorFinder_PRE();
	int pers_count=eCalculator.getNoOfPersons();
	//consolidator.person_count = pers_count;
	SITUATION_NO = 0;
	//falls in none of the 4 categories (explained in the header file), if 0 throw error and exit
	/*
	 * 1 indicates HPED HMAP: TRUE,ISMOL: TRUE
	 * 2 indicates HPED HMAP: TRUE,ISMOL: FALSE
	 * 3 indicates HPED HMAP: FALSE,ISMOL: TRUE
	 * 4 indicates HPED HMAP: FALSE,ISMOL: FALSE
	*/
	initiateErrorFinder_CURRENT(pers_count);
if ((SITUATION_NO == 1)	||	(SITUATION_NO == 3))	//findTruePctErrors
		//consolidator.findTruePctErrors( eCalculator, MA_SNP_END, true, WINDOW,MA_THRESHOLD, EMPIRICAL_MA_RESULT);
{

    if(T_TRIM )
    {
  	  std::cerr<<"entering into HoldOut Mode Matching"<<std::endl;
    }
    else
	{
	  if( OPTION.compare( "Error3") == 0 )
	  {
		std::cerr<< " Error: You have provided option:Error3 " << std::endl
			 <<" you can use Error3 only if you provided"
			 <<" hold out ped and map file, program with not output anything" << std::endl;
		exit( -1 );
	  }
	}
    ErrorCalculator &e_obj = eCalculator;
    int ma_snp_ends = MA_SNP_END;
    //bool holdOut,
    int window = WINDOW;
    float ma_threshold =MA_THRESHOLD;
    float empirical_ma_threshold = EMPIRICAL_MA_RESULT;

    //ErrorCalculator& e_obj=eCalculator;
    //int window=WINDOW;
	//int ma_snp_ends=MA_SNP_END;
	//float ma_threshold=MA_THRESHOLD;
	int min_snp=MIN_SNP;
	float min_cm=MIN_CM;
	float per_err_threshold=PCT_ERR_THRESHOLD;
	std::string option=OPTION;
	float hThreshold=HO_THRESHOLD;
	//bool holdOut,
	float empirical_threshold=EMPIRICAL_MA_RESULT;
	float empirical_pie_threshold=EMPIRICAL_PIE_RESULT;
	int extendSnp=EXTENDSNP;

	int removed1 =0, removed2 = 0, removed3 = 0, removed4 = 0;
	int not_removed = 0;
	int total_count = consolidator.global_initial;
	//bool wrongOption = false;
	float per_err_threshold1;
	  if(empirical_pie_threshold >= 0.0){
	    per_err_threshold1 = empirical_pie_threshold;
	  } else {
	    per_err_threshold1 = consolidator.getPctErrThreshold( per_err_threshold );
	  }
	  std::stringstream sstr;
	    sstr << std::fixed << std::setprecision(10) << per_err_threshold1;
	    std::string per_err_value = sstr.str();
	    consolidator.emp_pie_thresh_str = "empirical pie threshold is : " + per_err_value  + " \n";
	    float hThreshold1 = 0;
	    if( T_TRIM )
	    {
	      hThreshold1 = consolidator.getHoldOutThreshold( hThreshold );
	    }

	    per_err_threshold = per_err_threshold1;
	    hThreshold = hThreshold1;

	    bool ifcontd = false;	//whether the simplePCT loop with continue

	   for(int i=0;i<consolidator.person_count;i++)
	   {
	      for(int j=i;j<consolidator.person_count;j++)
	      {
	    	  for (int l=0;l<consolidator.m_matches[i][j].size();l++ ) //for(int l=0;l<consolidator.m_trueMatches[i][j].size();l++)
	          {
	    		  total_count++;
				  if (l<consolidator.m_trueMatches[i][j].size()	)
					  {
					  	  if (consolidator.m_trueMatches[i][j][l].end == -1	)
							  {
					  		  	  ifcontd = true;
							  }
					  	  else
							  {
								  ifcontd = false;
							  }
					  }
				  else
					  {
						  ifcontd = true;
					  }

/*	              if(consolidator.m_trueMatches[i][j][l].end==-1)
	              {
	                  continue;
	              }*/
					if (!ifcontd)
						{
						int t1_ = consolidator.m_trueMatches[ i ][ j ][ l ].start +
													  ( consolidator.m_trueMatches[ i ][ j ][ l ].end -
															  consolidator.m_trueMatches[ i ][ j ][ l ].start ) * 0.25;
									int t2_ = consolidator.m_trueMatches[ i ][ j ][ l ].end -
											   ( consolidator.m_trueMatches[ i ][ j ][ l ].end -
													   consolidator.m_trueMatches[ i ][ j ][ l ].start ) * 0.25;
							     //now we have the positions of the first and last 25% of the truly ibd SH
							     //all that's left to do is to pass them into the moving averages function, and obtain the max ma
							     //then store that in a vector, sort them, and find the xth percentile of that vector. That will be
							     //the ma that we use later
							     //for that "finalErrors" parameters, need to get the number of errors along the truly IBD SH first...
							              std::vector<std::vector<int> > trueErrors=e_obj.checkErrors( i, j, t1_, t2_);
							              std::vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );

							     //handles MA calculations
							     std::vector<float> av;
							     float current_max;
							     if(empirical_ma_threshold < 0.0){
							     av = e_obj.getTrueMovingAverages(finalTrueErrors,t1_,t2_,window);
							              current_max = av[0];
							              for(int q = 1; q < av.size(); q++){
							                   if(av[q] > current_max){
							                           current_max = av[q];
							                   }
							              }
							              e_obj.addMaxAverage(current_max);
							     }

							     //
							              int temp1_ = consolidator.m_trueMatches[ i ][ j ][ l ].start +
							                          ( consolidator.m_trueMatches[ i ][ j ][ l ].end -
							                        		  consolidator.m_trueMatches[ i ][ j ][ l ].start ) * 0.15; //Should probably stop doing this
							              int temp2_ = consolidator.m_trueMatches[ i ][ j ][ l ].end -
							                           ( consolidator.m_trueMatches[ i ][ j ][ l ].end -
							                        		   consolidator.m_trueMatches[ i ][ j ][ l ].start ) * 0.15;
							              int start =0, end =0, fend = ( temp2_ -temp1_ )  ;

							 //since we are using MOL at this point, this will pick out a random SH from the set of non-truly IBD SH
							 //and use that length to define the region over which we find PIE. Unless you are changing something with MOL,
							 //don't ever read this next block
							                  int randPers1, randPers2, pos;
							                  randPers1 = std::rand() % consolidator.person_count;
							                  randPers2 = std::rand() % consolidator.person_count;
							                  if( randPers1 > randPers2 )
							                  {
							                     randPers1 = randPers1 + randPers2;
							                     randPers2 = randPers1 - randPers2;
							                     randPers1 = randPers1 - randPers2;
							                  }
							                  while( consolidator.m_matches[ randPers1 ][ randPers2 ].size() <= 0 )
							                  {
							                    randPers1 = std::rand() % consolidator.person_count;
							                    randPers2 = std::rand() % consolidator.person_count;
							                   if( randPers1 > randPers2 )
							                   {
							                      randPers1 = randPers1 + randPers2;
							                      randPers2 = randPers1 - randPers2;
							                      randPers1 = randPers1 - randPers2;
							                   }

							                  }
							                  pos = std::rand() % consolidator.m_matches[ randPers1 ][ randPers2 ].size();
							                  int len = consolidator.m_matches[ randPers1 ][ randPers2 ][ pos ].end
							                            - consolidator.m_matches[ randPers1 ][ randPers2 ][ pos ].start;
							                  if( len >= fend || len <= 0)
							                  {
							                      continue;
							                  }
							                  //temp1 = temp1;//Why?
							                  temp2_ = temp1_ + len;
							 //end crazy MOL stuff
							                  std::vector<std::vector<int> > errors=e_obj.checkErrors( i, j, temp1_, temp2_);

							                  std::vector<int>finalErrors=e_obj.getFinalErrors( errors );
							//                  float per_err = e_obj.getThreshold(finalErrors,temp1,temp2,ma_snp_ends );
							                  float per_err = e_obj.getThreshold(finalErrors,temp1_,temp2_);//overload
							                  consolidator.m_errors.push_back( per_err );
							                 if( T_ERROR  )
							                 {
							                        float oppHom = ( e_obj.getOppHomThreshold( i, j, temp1_, temp2_ ) ) / ( temp2_ -temp1_ );
							                        consolidator.m_holdOutErrors.push_back( oppHom );
							                 }

						}
			        if(consolidator.m_matches[i][j][l].end==-1)
			        {
			          continue;
			        }
			        int temp1=consolidator.m_matches[i][j][l].start;
			                //cout<<"temp1 start begin= "<<temp1<<endl;

			        int temp2=consolidator.m_matches[i][j][l].end;
			                //cout<<"temp2 end begin= "<<temp2<<endl;


			        if (extendSnp != 0)
			        {

			                /*<piyush1>*/
			                if(temp1-extendSnp <0)
			                                             {
			                                             	temp1=0;
			                                             	//cout<<"New value of temp1= "<<temp1<<endl;
			                                             }
			                                             else
			                                             {
			                                             	temp1=consolidator.m_matches[i][j][l].start-extendSnp;
			                                             	//cout<<"New value of temp1= "<<temp1<<endl;
			                                             }

			                                             if(temp2+extendSnp > 4443)// change this constant
			                                             {
			                                             	temp2=consolidator.m_matches[i][j][l].end;
			                                             	//cout<<"New value of temp2= "<<temp2<<endl;
			                                             }
			                                             else
			                                             {
			                                             	temp2=consolidator.m_matches[i][j][l].end+extendSnp;
			                                             	//cout<<"New value of temp2= "<<temp2<<endl;
			                                             }


			                /*till here*/
			                                             /*cout<<"temp1 start after= "<<temp1<<endl;
			                                             cout<<"temp2 end after= "<<temp2<<endl;*/


			                                             //cout<<"perform trim temp1= "<<temp1<<endl;
			                                             //cout<<"perform trim temp2= "<<temp2<<endl;


			        }

			                int pers1 = i, pers2 = j;
			                if( option.compare( "ErrorRandom1" ) == 0 || option.compare( "ErrorRandom2" ) == 0 || option.compare( "ErrorRandom3" ) == 0 )
			                {
			                  pers1 = std::rand() % e_obj.getNoOfPersons();
			                  pers2 = std::rand() % e_obj.getNoOfPersons();
			                  if( pers1 > pers2 )
			                  {
			                    pers1 = pers1 + pers2;
			                    pers2 = pers1 - pers2;
			                    pers1 = pers1 - pers2;
			                  }
			                }

			                std::vector<std::vector<int> > errors=e_obj.checkErrors(pers1, pers2, temp1, temp2);
			                std::vector<int>finalErrors=e_obj.getFinalErrors(errors);//<piyush for errors>

			                //cout<<"finalErrors size= "<<finalErrors.size()<<endl;
			                /*Inject implied error at start/end of SH here*/
			                std::vector<int>::iterator it;
			                it = finalErrors.begin(); //go to the start of the vector
			                if(finalErrors[0] != 1){
			                  finalErrors.insert(it,1); //inject an error at position 1, if not already there
			                }
			                /*End inject implied error section*/

			                std::vector<int>trimPositions;
			                std::vector<float>movingAverages;
			                float threshold;
			                if( (e_obj.isInitialCmDrop(temp1,temp2,min_cm)) || ((temp2-temp1) < min_snp) ){ //initial drop. Don't calculate MA
			                  trimPositions.push_back(temp1);
			                  trimPositions.push_back(temp2);
			                  trimPositions.push_back(1);
			                }else{

			                  movingAverages = e_obj.getMovingAverages(finalErrors,temp1,temp2,window,extendSnp);//<piyush> get moving averages are calculated from this part
			                  if(empirical_threshold < 0.0){
			                    threshold = e_obj.getCutoff();
			                  } else {
			                    threshold = empirical_threshold;
			                  }
			                  trimPositions = e_obj.getTrimPositions(movingAverages,temp1,temp2,threshold,min_cm);

			                }
			                //-----------------

			                int beforeTrimStart = temp1;
			                int beforeTrimEnd = temp2;
			                consolidator.m_matches[i][j][l].end = temp2 = temp1+trimPositions[1];
			                consolidator.m_matches[i][j][l].start = temp1 = temp1+trimPositions[0];
			                int del0 = trimPositions[0];
			                int del1 = trimPositions[1];
			                float per_err = e_obj.getThreshold(finalErrors,del0,del1,ma_snp_ends);

			                //add new weighted option
			                /*
			                 For this new option, we only output SH that are not dropped. So, the output is finalOutput + weighted column.
			                */
			                if( (option.compare("weightedOutput") == 0) || (option.compare("weightedOutputBP") == 0) ){
			                  int snp1 = 0, snp2 = 0, hlength = 0;
			                  float noOfOppHom = 0;
			                  if( T_TRIM )
			                  {
			                    snp1 = e_obj.getNewSnp( temp1 );
			                    snp2 = e_obj.getNewSnp( temp2 );
			                    hlength = snp2 - snp1;
			                    if( hlength <= 0 )
			                    {
			                      hlength = 1;
			                    }
			                    noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, consolidator.m_matches[i][j][l].start, consolidator.m_matches[i][j][l].end );
			                  }
			                  if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){
			                	  consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    removed4++;
			                    continue;
			                  }
			                  if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ //removed2 a tpos.size of 3 indicates trimming due ot cM
			                    removed2++;
			                    consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    continue;
			                  }
			                  if( per_err > per_err_threshold){
			                    removed1++;
			                    continue;
			                  }
			                  if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength ){
			                    removed3++;
			                    consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    continue;
			                  } //removed3

			                  not_removed++;
			                  consolidator.m_weighted_sh.push_back(Weighted_SH(temp1,temp2,i,j)); //build the vector of SH that passed
			                  continue;
			                }//end weghtedOutput
			                /*Add new finalErrorsOutput*/
			                if( (option.compare("finalErrorsOutput") == 0) ){
			                  int snp1 = 0, snp2 = 0, hlength = 0;
			                  float noOfOppHom = 0;

			                  if( T_TRIM )
			                  {
			                    snp1 = e_obj.getNewSnp( temp1 );
			                    snp2 = e_obj.getNewSnp( temp2 );
			                    hlength = snp2 - snp1;
			                    if( hlength <= 0 )
			                    {
			                      hlength = 1;
			                    }
			                    noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, consolidator.m_matches[i][j][l].start, consolidator.m_matches[i][j][l].end );
			                  }


			                  if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){
			                    std::vector<float>movingAverages;
			                    temp1 = beforeTrimStart;
			                    temp2 = beforeTrimEnd;
			                    consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    removed4++;
			                    continue;
			                  }

			                  if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ //removed2 a tpos.size of 3 indicates trimming due ot cM
			                    removed2++;
			                    consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    continue;
			                  }
			                  if( per_err > per_err_threshold){
			                    removed1++;
			                    continue;
			                  }

			                  if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength ){
			                    removed3++;
			                    consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    continue;
			                  } //removed3
			                  not_removed++;
			                  e_obj.finalErrorsOutput(i,j,temp1,temp2,min_cm,per_err);// <piyush>this is where the final output is written is called at
			                  continue;
			                }//end finalErrorsOutput
			                if( (option.compare("FullPlusDropped") == 0) ){
			                  int snp1 = 0, snp2 = 0, hlength = 0;
			                  float noOfOppHom = 0;

			                  if( T_TRIM )
			                  {
			                    snp1 = e_obj.getNewSnp( temp1 );
			                    snp2 = e_obj.getNewSnp( temp2 );
			                    hlength = snp2 - snp1;
			                    if( hlength <= 0 )
			                    {
			                      hlength = 1;
			                    }
			                    noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, consolidator.m_matches[i][j][l].start, consolidator.m_matches[i][j][l].end );
			                  }


			                  if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){
			                    std::vector<float>movingAverages;
			                    temp1 = beforeTrimStart;
			                    temp2 = beforeTrimEnd;
			                    e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,1);//standardize the error codes
			                    consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    removed4++;
			                    continue;
			                  }

			                  if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ //removed2 a tpos.size of 3 indicates trimming due ot cM
			                    e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,2);
			                    removed2++;
			                    consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    continue;
			                  }
			                  if( per_err > per_err_threshold){
			                    e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,3);
			                    removed1++;
			                    continue;
			                  }

			                  if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength ){
			                    e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,4);
			                    removed3++;
			                    consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    continue;
			                  } //removed3
			                  not_removed++;
			                  e_obj.finalOutPut(i,j,temp1,temp2,min_cm);
			                  continue;
			                } //end FullPlusDropped

			                //Calculate Error1
			                if( (option.compare("Error1") == 0 ) || (option.compare("ErrorRandom1") == 0) || (option.compare("Error") == 0) ){

			                  if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){ //dropped before trimming
			                    //don't bother printing out ma for this one. But go back and change it so that it doesn't actually calc it
			                    std::vector<float>movingAverages;//null
			                    //trying something special in this case. This can be removed once idrops aren't being trimmed
			                    //test code
			                    temp1 = beforeTrimStart;
			                    temp2 = beforeTrimEnd;
			                    //
			                    e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,1);
			                    consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    removed4++; //seems ok
			                    continue;
			                  }
			                  if( (( temp2-temp1 ) < min_snp) || ((trimPositions.size() == 3) && (trimPositions[2] == 2) ) ) //dropped after trimming
			                  {
			                    e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,2);
			                    ++removed2;
			                    consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    continue;
			                  }
			                  if( per_err > per_err_threshold ) //dropped due to pie
			                  {
			                    e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,3);
			                    ++removed1;
			                    consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                    continue;
			                  }
			                  not_removed++;
			                  e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,0);//no drop
			                  continue;
			                }//end error1

			                int snp1 = 0, snp2 = 0, hlength = 0;
			                float noOfOppHom = 0;
			                if( T_TRIM )
			                {
			                  snp1 = e_obj.getNewSnp( temp1 );
			                  snp2 = e_obj.getNewSnp( temp2 );
			                  hlength = snp2 - snp1;
			                  if( hlength <= 0 )
			                  {
			                    hlength = 1;
			                  }
			                  noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, consolidator.m_matches[i][j][l].start, consolidator.m_matches[i][j][l].end );
			                }
			                //update drop order 2/26/14
			                if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) )
			                {
			                  ++removed2;
			                  consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                  continue;
			                }

			                if( per_err > per_err_threshold )
			                {

			                  ++removed1;
			                  consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                  continue;
			                }
			                //probably not removed?
			                not_removed++;
			                if( option.compare("MovingAverages")==0 ) //make this ma2
			                {
			                  if( T_TRIM)
			                  {
			                    e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp,min_cm,movingAverages,trimPositions,per_err, noOfOppHom, hlength );
			                  }
			                  else
			                  {
			                    e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm,movingAverages, trimPositions,per_err );
			                  }
			                  continue;
			                }

			                if(option.compare("Error2")==0 || option.compare( "ErrorRandom2" ) == 0)
			                {
			                  if( T_TRIM)
			                  {
			                    e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
			                  }
			                  else
			                  {
			                    e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err);
			                  }
			                  continue;
			                }
			                if ( T_TRIM && hThreshold < ( noOfOppHom ) / hlength )
			                {
			                  ++removed3;
			                  consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			                  continue;
			                }
			                if( option.compare("Error3")==0 || option.compare( "ErrorRandom3" ) == 0  )
			                {
			                  e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
			                }

	     //handle moving averages calculation
	     //

	          }//inner for loop closed
	       }//middle for loop  closed
	   }//outer for loop closed

		std::vector<float>maxes;
		float cutoff = empirical_ma_threshold;
		if(empirical_ma_threshold < 0.0){
		maxes = e_obj.getMaxAverages();
		std::sort(maxes.begin(),maxes.end());
		e_obj.setMaxAverage(maxes);
		cutoff = e_obj.getXthPercentile(ma_threshold);
		}
		e_obj.setCutoff(cutoff);//set the actual threshold to be used when calculating MA in all other SH
		//
		std::sort( consolidator.m_errors.begin(), consolidator.m_errors.end() );
		std::sort( consolidator.m_holdOutErrors.begin(), consolidator.m_holdOutErrors.end() );
		std::string str =  " \n No of elements in error check are: "
						+ NumberToString( consolidator.m_errors.size() );
		str  = str + " \n No of elements in hold  error check are: "
						+ NumberToString( consolidator.m_holdOutErrors.size() );
		e_obj.log( str );
		  if( option.compare("weightedOutput") == 0 ){
		    float snp_average_count = 0.0;
		    int start_position;
		    int end_position;
		    int genome_length;
		    if(consolidator.isUserSuppliedWeights()){ //the user has supplied their own weights.
		      //in this case, the min and max values correspond to the number of lines in the input file,
		      //since each line represents a snp. So the min is always 0, and the max is always the number of lines-1.
		      start_position = 0;
		      end_position = consolidator.user_supplied_snp_weights.size() - 1;
		    }else {
		      start_position = consolidator.find_genome_min();
		      end_position = consolidator.find_genome_max();
		    }//end else
		    genome_length = (end_position - start_position)+1;
		    consolidator.genome_vector.resize(genome_length,0);
		    if(consolidator.isUserSuppliedWeights()){
		      for(int i = 0; i < consolidator.user_supplied_snp_weights.size(); i++){
		    	  consolidator.update_genome(i,consolidator.user_supplied_snp_weights[i]);
		      }
		    }else{
		      /*This next for loop adds one to each snp in a SH. Bypass it if the user gives a files of weights*/
		      for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
		    	  consolidator.update_genome(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2);
		      }
		    }
		    //this part is next...will probably need to add stuff to that weighted object...
		    snp_average_count = consolidator.average_snp_count();
		    for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
		    	consolidator.m_weighted_sh[i].snp_weight = consolidator.update_snp_weight(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2);
		    }
		    for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
		    	consolidator.m_weighted_sh[i].final_weight = ( snp_average_count / (consolidator.m_weighted_sh[i].snp_weight));
		      e_obj.weightedOutput(consolidator.m_weighted_sh[i].per1, consolidator.m_weighted_sh[i].per2, consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2, consolidator.m_weighted_sh[i].final_weight);
		    }
		  }
		  if (option.compare("weightedOutputBP") == 0){

		  //begin new test code section here: Dec 4th 2014
		  int genome_length = e_obj.getGenomeBPLength();
		  float adjusted_genome_length = genome_length / 1000.0; //L using kbp for now
		  int genome_min = e_obj.getMinimumBP(); std::cout<<"genome_min= "<<genome_min<<std::endl;
		  int genome_max = e_obj.getMaximumBP(); std::cout<<"genome_max= "<<genome_max<<std::endl;
		  int genome_size_snps = (consolidator.find_genome_max() - consolidator.find_genome_min())+1; //used for genome_vector
		  float wprime_numerator = 0.0;  //This is Ci / L
		  float total_sh_length_sum = 0.0;
		  float w2prime_denominator = 0.0;

		  consolidator.genome_vector.resize(genome_size_snps,0); //resize and zero out the genome. shit that needs to be snps.

		  //update all of the snp counts in the genome. This looks fine.
		  for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
			  consolidator.update_genome(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2);
		  }

		  //calculate the w' numerator by summing up all of the snp counts and dividing by the genome length.
		  //WARNING: This can cause wprime_numerator to overflow. Currently using kbp units to avoid this, but
		  //this needs to be addressed.
		  for(int i = 0; i < consolidator.genome_vector.size(); i++){
		    wprime_numerator += consolidator.genome_vector[i] / adjusted_genome_length;
		  }

		  //Calculate w' for each SH.
		  for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
		    float wprime_denominator = 0.0;
		    consolidator.m_weighted_sh[i].mbp_length = (e_obj.getSHBPLength(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2)/1000.0);
		    wprime_denominator = consolidator.get_snps_over_range(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2, consolidator.m_weighted_sh[i].mbp_length);
		    consolidator.m_weighted_sh[i].wprime = wprime_numerator / wprime_denominator;
		  }

		  //This is the total length of all SH. This can probably overflow as well...ugh.
		  for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
		    total_sh_length_sum += consolidator.m_weighted_sh[i].mbp_length;
		  }

		  //Calculate the w2prime denominator - this value is a constant
		  for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
		    float temp = consolidator.m_weighted_sh[i].mbp_length * consolidator.m_weighted_sh[i].wprime;
		    w2prime_denominator += temp / total_sh_length_sum;
		  }

		  //Calculate and output w2' for each SH
		  for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
			  consolidator.m_weighted_sh[i].w2prime = (consolidator.m_weighted_sh[i].wprime) / w2prime_denominator;
		    e_obj.weightedOutput(consolidator.m_weighted_sh[i].per1, consolidator.m_weighted_sh[i].per2, consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2, consolidator.m_weighted_sh[i].w2prime);
		  }
		}

		  /*End weighted output*/



		  /*END TESTING AREA DEC 4th 2014*/
		  /*****************************
		  ******************************/
		  consolidator.ma_drop_str = "No of matches removed due to length of trimming by moving averages: " + NumberToString( removed2 );
		  consolidator.pie_drop_str = "No of matches removed due to percentage error: " + NumberToString( removed1 );
		  if(T_TRIM){
		  //  str = str+ " \n No of matches removed due hold out ped file checking: "+ NumberToString( removed3 );
		  }
		  //begin log output
		  std::string parameter_string_1 = "\n\n**********Parameters used in program**********\n";
		  e_obj.log(parameter_string_1);
		  e_obj.log(consolidator.emp_ma_thresh_str); //keep
		  e_obj.log(consolidator.emp_pie_thresh_str);//keep
		  parameter_string_1 = "**********************************************\n\n";
		  e_obj.log(parameter_string_1);
		  std::string total_count_str = "The total number of SH in the input file was: " + NumberToString(total_count);
		  e_obj.log(total_count_str);
		  e_obj.log(consolidator.consolidated_str);
		  e_obj.log(consolidator.initial_drop_str);
		  //  e_obj.log(ibg_str);
		  e_obj.log(consolidator.ma_drop_str);
		  e_obj.log(consolidator.pie_drop_str);
		  consolidator.final_sh_str = "Total number of SH that were not dropped is: " + NumberToString(not_removed);
		  e_obj.log(consolidator.final_sh_str);







}//ending of findTruePctErrors & error




















else
	// for (	(situation_no == 2)	&&	(situation_no == 4)		)
{
			//consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT );




	if (T_TRIM)
	{
		std::cerr<<"entering into HoldOut Mode Matching"<<std::endl;
	}
	else
	{
		 if( OPTION.compare( "Error3") == 0 )
		            {
		            	std::cerr<< " Error: You have provided option:Error3 " << std::endl
		                   <<" you can use Error3 only if you provided"
		                   <<" hold out ped and map file, program with not output anything" << std::endl;
		              exit( -1 );
		            }
	}


ErrorCalculator &e_obj= eCalculator;
float PIElength = PIELENGTH;
int window=WINDOW;
float ma_threshold = MA_THRESHOLD;
float empirical_ma_threshold =EMPIRICAL_MA_RESULT;
/*float current_max=0.0;*/
int ma_snp_ends=MA_SNP_END;
//float ma_threshold=MA_THRESHOLD;
int min_snp=MIN_SNP;
float min_cm=MIN_CM;
float per_err_threshold=PCT_ERR_THRESHOLD;
std::string option=OPTION;
float hThreshold=HO_THRESHOLD;
//bool holdOut,
float empirical_threshold=EMPIRICAL_MA_RESULT;
float empirical_pie_threshold=EMPIRICAL_PIE_RESULT;
int extendSnp=EXTENDSNP;
int removed1 =0, removed2 = 0, removed3 = 0, removed4 = 0;
int not_removed = 0;
int total_count = consolidator.global_initial;
//bool wrongOption = false;
float per_err_threshold1;

if(empirical_pie_threshold >= 0.0){
per_err_threshold1 = empirical_pie_threshold;
} else {
per_err_threshold1 = consolidator.getPctErrThreshold( per_err_threshold );
}
std::stringstream sstr;
sstr << std::fixed << std::setprecision(10) << per_err_threshold1;
std::string per_err_value = sstr.str();
consolidator.emp_pie_thresh_str = "empirical pie threshold is : " + per_err_value  + " \n";
float hThreshold1 = 0;
if( T_TRIM )
{
hThreshold1 = consolidator.getHoldOutThreshold( hThreshold );
}

per_err_threshold = per_err_threshold1;
hThreshold = hThreshold1;



/*std::cout<<"sddsdsdsdsd"<<consolidator.person_count;
std::cout<<std::endl;
exit(0);*/

bool ifcontd = false;	//whether the simplePCT loop with continue

for(int i=0;i<consolidator.person_count;i++)
{
	for(int j=i;j<consolidator.person_count;j++)
		{
		  for (int l=0;l<consolidator.m_matches[i][j].size();l++ ) // for(int l=0;l<consolidator.m_trueMatches[i][j].size();l++)
			{
			  total_count++;
			  if (l<consolidator.m_trueMatches[i][j].size()	)
				  {
				  	  if (consolidator.m_trueMatches[i][j][l].end == -1	)
						  {
				  		  	  ifcontd = true;
						  }
				  	  else
						  {
							  ifcontd = false;
						  }
				  }
			  else
				  {
					  ifcontd = true;
				  }
/*

			  if(consolidator.m_trueMatches[i][j][l].end==-1)
			  {
				  continue;
			  }
*/

			//-------------------------------------------------------------------------------------------------

		if (!ifcontd)
		{

			int t1_ = consolidator.m_trueMatches[ i ][ j ][ l ].start +
						  ( consolidator.m_trueMatches[ i ][ j ][ l ].end -
								  consolidator.m_trueMatches[ i ][ j ][ l ].start ) * 0.25;
			int t2_ = consolidator.m_trueMatches[ i ][ j ][ l ].end -
			   ( consolidator.m_trueMatches[ i ][ j ][ l ].end -
					   consolidator.m_trueMatches[ i ][ j ][ l ].start ) * 0.25;

			std::vector<std::vector<int> > trueErrors=e_obj.checkErrors( i, j, t1_, t2_);
			std::vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );
			//This section handles finding the maximum moving averages amongst trulyIBD segments
			std::vector<float> av;
			float current_max;
			if(empirical_ma_threshold < 0.0){
			 av = e_obj.getTrueMovingAverages(finalTrueErrors,t1_,t2_,window);//<true moving avg>
			 current_max = av[0];
			 for(int q = 1; q < av.size(); q++){
				   if(av[q] > current_max){
						   current_max = av[q];
				   }
			  }
			 e_obj.addMaxAverage(current_max);
			}
			//------------------------------------------------------------------------------

			int temp1_ = consolidator.m_trueMatches[i][j][l].start;
			int temp2_ = consolidator.m_trueMatches[i][j][l].end;
			float startCM = e_obj.getCMDistance( temp1_ );
			float endCM = e_obj.getCMDistance( temp2_ );
			float mid1CM = startCM + ( endCM - startCM ) / 2 - PIElength / 2;
			float mid2CM = startCM + ( endCM - startCM ) / 2 + PIElength / 2;
			while( e_obj.getCMDistance( temp1_ ) <= mid1CM || e_obj.getCMDistance( temp2_ ) >=mid2CM )
			 {
				if( e_obj.getCMDistance( temp1_ ) <= mid1CM )
				{
				  ++temp1_;
				}
				if( e_obj.getCMDistance( temp2_ ) >=mid2CM )
				{
				  --temp2_;
				}
			 }
			 std::vector<std::vector<int> > errors_=e_obj.checkErrors( i, j, temp1_, temp2_);
			 std::vector<int>finalErrors_=e_obj.getFinalErrors( errors_ );
			//                  float per_err = e_obj.getThreshold(finalErrors,temp1, temp2, 0 );
			 float per_err_ = e_obj.getThreshold(finalErrors_,temp1_,temp2_); //overload!
			 consolidator.m_errors.push_back( per_err_ );
			 if( T_ERROR  )
			 {
					float oppHom =
								 ( e_obj.getOppHomThreshold( i, j,
										 temp1_,
										 temp2_ ) ) / ( temp2_ -temp1_ );
					consolidator.m_holdOutErrors.push_back( oppHom );
			 }
		}
//----------------------------------------------------------------
		        if(consolidator.m_matches[i][j][l].end==-1)
		        {
		          continue;
		        }





			 int temp1=consolidator.m_matches[i][j][l].start;
			        //cout<<"temp1 start begin= "<<temp1<<endl;

			int temp2=consolidator.m_matches[i][j][l].end;
			//cout<<"temp2 end begin= "<<temp2<<endl;
			if (extendSnp != 0)
			{

			        /*<piyush1>*/
			        if(temp1-extendSnp <0)
			                                     {
			                                     	temp1=0;
			                                     	//cout<<"New value of temp1= "<<temp1<<endl;
			                                     }
			                                     else
			                                     {
			                                     	temp1=consolidator.m_matches[i][j][l].start-extendSnp;
			                                     	//cout<<"New value of temp1= "<<temp1<<endl;
			                                     }

			                                     if(temp2+extendSnp > 4443)// change this constant
			                                     {
			                                     	temp2=consolidator.m_matches[i][j][l].end;
			                                     	//cout<<"New value of temp2= "<<temp2<<endl;
			                                     }
			                                     else
			                                     {
			                                     	temp2=consolidator.m_matches[i][j][l].end+extendSnp;
			                                     	//cout<<"New value of temp2= "<<temp2<<endl;
			                                     }


						/*till here*/
						/*cout<<"temp1 start after= "<<temp1<<endl;
						cout<<"temp2 end after= "<<temp2<<endl;*/


						//cout<<"perform trim temp1= "<<temp1<<endl;
						//cout<<"perform trim temp2= "<<temp2<<endl;
			}
			int pers1 = i, pers2 = j;
			        if( option.compare( "ErrorRandom1" ) == 0 || option.compare( "ErrorRandom2" ) == 0 || option.compare( "ErrorRandom3" ) == 0 )
			        {
			          pers1 = std::rand() % e_obj.getNoOfPersons();
			          pers2 = std::rand() % e_obj.getNoOfPersons();
			          if( pers1 > pers2 )
			          {
			            pers1 = pers1 + pers2;
			            pers2 = pers1 - pers2;
			            pers1 = pers1 - pers2;
			          }
			        }

			        std::vector<std::vector<int> > errors=e_obj.checkErrors(pers1, pers2, temp1, temp2);
			        std::vector<int>finalErrors=e_obj.getFinalErrors(errors);//<piyush for errors>

			        //cout<<"finalErrors size= "<<finalErrors.size()<<endl;
			        /*Inject implied error at start/end of SH here*/
			        std::vector<int>::iterator it;
			        it = finalErrors.begin(); //go to the start of the vector
			        if(finalErrors[0] != 1){
			          finalErrors.insert(it,1); //inject an error at position 1, if not already there
			        }
			        /*End inject implied error section*/

			        std::vector<int>trimPositions;
			        std::vector<float>movingAverages;
			        float threshold;
			        if( (e_obj.isInitialCmDrop(temp1,temp2,min_cm)) || ((temp2-temp1) < min_snp) ){ //initial drop. Don't calculate MA
			          trimPositions.push_back(temp1);
			          trimPositions.push_back(temp2);
			          trimPositions.push_back(1);
			        }else{

			          movingAverages = e_obj.getMovingAverages(finalErrors,temp1,temp2,window,extendSnp);//<piyush> get moving averages are calculated from this part
			          if(empirical_threshold < 0.0){
			            threshold = e_obj.getCutoff();
			          } else {
			            threshold = empirical_threshold;
			          }
			          trimPositions = e_obj.getTrimPositions(movingAverages,temp1,temp2,threshold,min_cm);

			        }
			        //-----------------

			        int beforeTrimStart = temp1;
			        int beforeTrimEnd = temp2;
			        consolidator.m_matches[i][j][l].end = temp2 = temp1+trimPositions[1];
			        consolidator.m_matches[i][j][l].start = temp1 = temp1+trimPositions[0];
			        int del0 = trimPositions[0];
			        int del1 = trimPositions[1];
			        float per_err = e_obj.getThreshold(finalErrors,del0,del1,ma_snp_ends);

			        //add new weighted option
			        /*
			         For this new option, we only output SH that are not dropped. So, the output is finalOutput + weighted column.
			        */
			        if( (option.compare("weightedOutput") == 0) || (option.compare("weightedOutputBP") == 0) ){
			          int snp1 = 0, snp2 = 0, hlength = 0;
			          float noOfOppHom = 0;
			          if( T_TRIM )
			          {
			            snp1 = e_obj.getNewSnp( temp1 );
			            snp2 = e_obj.getNewSnp( temp2 );
			            hlength = snp2 - snp1;
			            if( hlength <= 0 )
			            {
			              hlength = 1;
			            }
			            noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, consolidator.m_matches[i][j][l].start, consolidator.m_matches[i][j][l].end );
			          }
			          if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){
			        	  consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            removed4++;
			            continue;
			          }
			          if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ //removed2 a tpos.size of 3 indicates trimming due ot cM
			            removed2++;
			            consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            continue;
			          }
			          if( per_err > per_err_threshold){
			            removed1++;
			            continue;
			          }
			          if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength ){
			            removed3++;
			            consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            continue;
			          } //removed3

			          not_removed++;
			          consolidator.m_weighted_sh.push_back(Weighted_SH(temp1,temp2,i,j)); //build the vector of SH that passed
			          continue;
			        }//end weghtedOutput
			        /*Add new finalErrorsOutput*/
			        if( (option.compare("finalErrorsOutput") == 0) ){
			          int snp1 = 0, snp2 = 0, hlength = 0;
			          float noOfOppHom = 0;

			          if( T_TRIM )
			          {
			            snp1 = e_obj.getNewSnp( temp1 );
			            snp2 = e_obj.getNewSnp( temp2 );
			            hlength = snp2 - snp1;
			            if( hlength <= 0 )
			            {
			              hlength = 1;
			            }
			            noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, consolidator.m_matches[i][j][l].start, consolidator.m_matches[i][j][l].end );
			          }


			          if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){
			            std::vector<float>movingAverages;
			            temp1 = beforeTrimStart;
			            temp2 = beforeTrimEnd;
			            consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            removed4++;
			            continue;
			          }

			          if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ //removed2 a tpos.size of 3 indicates trimming due ot cM
			            removed2++;
			            consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            continue;
			          }
			          if( per_err > per_err_threshold){
			            removed1++;
			            continue;
			          }

			          if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength ){
			            removed3++;
			            consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            continue;
			          } //removed3
			          not_removed++;
			          e_obj.finalErrorsOutput(i,j,temp1,temp2,min_cm,per_err);// <piyush>this is where the final output is written is called at
			          continue;
			        }//end finalErrorsOutput
			        if( (option.compare("FullPlusDropped") == 0) ){
			          int snp1 = 0, snp2 = 0, hlength = 0;
			          float noOfOppHom = 0;

			          if( T_TRIM )
			          {
			            snp1 = e_obj.getNewSnp( temp1 );
			            snp2 = e_obj.getNewSnp( temp2 );
			            hlength = snp2 - snp1;
			            if( hlength <= 0 )
			            {
			              hlength = 1;
			            }
			            noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, consolidator.m_matches[i][j][l].start, consolidator.m_matches[i][j][l].end );
			          }


			          if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){
			            std::vector<float>movingAverages;
			            temp1 = beforeTrimStart;
			            temp2 = beforeTrimEnd;
			            e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,1);//standardize the error codes
			            consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            removed4++;
			            continue;
			          }

			          if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ //removed2 a tpos.size of 3 indicates trimming due ot cM
			            e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,2);
			            removed2++;
			            consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            continue;
			          }
			          if( per_err > per_err_threshold){
			            e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,3);
			            removed1++;
			            continue;
			          }

			          if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength ){
			            e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,4);
			            removed3++;
			            consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            continue;
			          } //removed3
			          not_removed++;
			          e_obj.finalOutPut(i,j,temp1,temp2,min_cm);
			          continue;
			        } //end FullPlusDropped

			        //Calculate Error1
			        if( (option.compare("Error1") == 0 ) || (option.compare("ErrorRandom1") == 0) || (option.compare("Error") == 0) ){

			          if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){ //dropped before trimming
			            //don't bother printing out ma for this one. But go back and change it so that it doesn't actually calc it
			            std::vector<float>movingAverages;//null
			            //trying something special in this case. This can be removed once idrops aren't being trimmed
			            //test code
			            temp1 = beforeTrimStart;
			            temp2 = beforeTrimEnd;
			            //
			            e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,1);
			            consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            removed4++; //seems ok
			            continue;
			          }
			          if( (( temp2-temp1 ) < min_snp) || ((trimPositions.size() == 3) && (trimPositions[2] == 2) ) ) //dropped after trimming
			          {
			            e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,2);
			            ++removed2;
			            consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            continue;
			          }
			          if( per_err > per_err_threshold ) //dropped due to pie
			          {
			            e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,3);
			            ++removed1;
			            consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			            continue;
			          }
			          not_removed++;
			          e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,0);//no drop
			          continue;
			        }//end error1

			        int snp1 = 0, snp2 = 0, hlength = 0;
			        float noOfOppHom = 0;
			        if( T_TRIM )
			        {
			          snp1 = e_obj.getNewSnp( temp1 );
			          snp2 = e_obj.getNewSnp( temp2 );
			          hlength = snp2 - snp1;
			          if( hlength <= 0 )
			          {
			            hlength = 1;
			          }
			          noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, consolidator.m_matches[i][j][l].start, consolidator.m_matches[i][j][l].end );
			        }
			        //update drop order 2/26/14
			        if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) )
			        {
			          ++removed2;
			          consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			          continue;
			        }

			        if( per_err > per_err_threshold )
			        {

			          ++removed1;
			          consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			          continue;
			        }
			        //probably not removed?
			        not_removed++;
			        if( option.compare("MovingAverages")==0 ) //make this ma2
			        {
			          if( T_TRIM)
			          {
			            e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp,min_cm,movingAverages,trimPositions,per_err, noOfOppHom, hlength );
			          }
			          else
			          {
			            e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm,movingAverages, trimPositions,per_err );
			          }
			          continue;
			        }

			        if(option.compare("Error2")==0 || option.compare( "ErrorRandom2" ) == 0)
			        {
			          if( T_TRIM)
			          {
			            e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
			          }
			          else
			          {
			            e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err);
			          }
			          continue;
			        }
			        if ( T_TRIM && hThreshold < ( noOfOppHom ) / hlength )
			        {
			          ++removed3;
			          consolidator.m_matches[i][j][l].start= consolidator.m_matches[i][j][l].end=-1;
			          continue;
			        }
			        if( option.compare("Error3")==0 || option.compare( "ErrorRandom3" ) == 0  )
			        {
			          e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
			        }

		}
	}
}

	   //this section actually handles the sorting of the max averages, and the setting of the user supplied percentile.
	   std::vector<float>maxes;
	   float cutoff = empirical_ma_threshold; //assume the user wanted to supply a value. This value will be overwritten shortly if they did not.
	   if(empirical_ma_threshold < 0.0){
	   maxes = e_obj.getMaxAverages();
	   std::sort(maxes.begin(),maxes.end());
	   e_obj.setMaxAverage(maxes);
	   cutoff = e_obj.getXthPercentile(ma_threshold); //<-make that an actual user input value
	   }

	   e_obj.setCutoff(cutoff);//set the actual threshold to be used when calculating MA in all other SH

	   if(empirical_ma_threshold < 0.0){
	   //std::string outt = "\n User supplied ma-threshold is: " + NumberToString(ma_threshold);
	   //outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
	   //e_obj.log(outt);
		   consolidator.ma_thresh_str = "User supplied ma-threshold is: " + NumberToString(ma_threshold);
		   consolidator.emp_ma_thresh_str = "Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
	   } else {
	   //std::string outt = "\n User supplied empirical-ma-threshold is: " + NumberToString(cutoff);
	   //outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
		   consolidator.emp_ma_thresh_str = "Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
	   //e_obj.log(outt);
	   }
	   //----------------------------------------


	   std::sort( consolidator.m_errors.begin(), consolidator.m_errors.end() );
	   std::sort( consolidator.m_holdOutErrors.begin(), consolidator.m_holdOutErrors.end() );
	   /*std::string str =  " \n No of segments deemed to be IBD for finding empirical error threshold "
	                      + NumberToString( m_errors.size() );*/
	   consolidator.ibg_str = "No of segments deemed to be IBD for finding empirical error threshold "
	                      + NumberToString( consolidator.m_errors.size() );
	   /*str  = str + " \n No of segments deemed to be IBD for finding empirical error threshold in hold out are "
	                      + NumberToString( m_holdOutErrors.size() );*/
	        //e_obj.log( str );
	   if( option.compare("weightedOutput") == 0 ){
	       float snp_average_count = 0.0;
	       int start_position;
	       int end_position;
	       int genome_length;
	       if(consolidator.isUserSuppliedWeights()){ //the user has supplied their own weights.
	         //in this case, the min and max values correspond to the number of lines in the input file,
	         //since each line represents a snp. So the min is always 0, and the max is always the number of lines-1.
	         start_position = 0;
	         end_position = consolidator.user_supplied_snp_weights.size() - 1;
	       }else {
	         start_position = consolidator.find_genome_min();
	         end_position = consolidator.find_genome_max();
	       }//end else
	       genome_length = (end_position - start_position)+1;
	       consolidator.genome_vector.resize(genome_length,0);
	       if(consolidator.isUserSuppliedWeights()){
	         for(int i = 0; i < consolidator.user_supplied_snp_weights.size(); i++){
	        	 consolidator.update_genome(i,consolidator.user_supplied_snp_weights[i]);
	         }
	       }else{
	         /*This next for loop adds one to each snp in a SH. Bypass it if the user gives a files of weights*/
	         for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
	        	 consolidator.update_genome(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2);
	         }
	       }
	       //this part is next...will probably need to add stuff to that weighted object...
	       snp_average_count = consolidator.average_snp_count();
	       for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
	    	   consolidator.m_weighted_sh[i].snp_weight = consolidator.update_snp_weight(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2);
	       }
	       for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
	    	   consolidator.m_weighted_sh[i].final_weight = ( snp_average_count / (consolidator.m_weighted_sh[i].snp_weight));
	         e_obj.weightedOutput(consolidator.m_weighted_sh[i].per1, consolidator.m_weighted_sh[i].per2, consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2, consolidator.m_weighted_sh[i].final_weight);
	       }
	     }

	     if (option.compare("weightedOutputBP") == 0){

	     //begin new test code section here: Dec 4th 2014
	     int genome_length = e_obj.getGenomeBPLength();
	     float adjusted_genome_length = genome_length / 1000.0; //L using kbp for now
	     int genome_min = e_obj.getMinimumBP(); std::cout<<"genome_min= "<<genome_min<<std::endl;
	     int genome_max = e_obj.getMaximumBP(); std::cout<<"genome_max= "<<genome_max<<std::endl;
	     int genome_size_snps = (consolidator.find_genome_max() - consolidator.find_genome_min())+1; //used for genome_vector
	     float wprime_numerator = 0.0;  //This is Ci / L
	     float total_sh_length_sum = 0.0;
	     float w2prime_denominator = 0.0;

	     consolidator.genome_vector.resize(genome_size_snps,0); //resize and zero out the genome. shit that needs to be snps.

	     //update all of the snp counts in the genome. This looks fine.
	     for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
	    	 consolidator.update_genome(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2);
	     }

	     //calculate the w' numerator by summing up all of the snp counts and dividing by the genome length.
	     //WARNING: This can cause wprime_numerator to overflow. Currently using kbp units to avoid this, but
	     //this needs to be addressed.
	     for(int i = 0; i < consolidator.genome_vector.size(); i++){
	       wprime_numerator += consolidator.genome_vector[i] / adjusted_genome_length;
	     }

	     //Calculate w' for each SH.
	     for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
	       float wprime_denominator = 0.0;
	       consolidator.m_weighted_sh[i].mbp_length = (e_obj.getSHBPLength(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2)/1000.0);
	       wprime_denominator = consolidator.get_snps_over_range(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2, consolidator.m_weighted_sh[i].mbp_length);
	       consolidator.m_weighted_sh[i].wprime = wprime_numerator / wprime_denominator;
	     }

	     //This is the total length of all SH. This can probably overflow as well...ugh.
	     for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
	       total_sh_length_sum += consolidator.m_weighted_sh[i].mbp_length;
	     }

	     //Calculate the w2prime denominator - this value is a constant
	     for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
	       float temp = consolidator.m_weighted_sh[i].mbp_length * consolidator.m_weighted_sh[i].wprime;
	       w2prime_denominator += temp / total_sh_length_sum;
	     }

	     //Calculate and output w2' for each SH
	     for(int i = 0; i < consolidator.m_weighted_sh.size(); i++){
	    	 consolidator.m_weighted_sh[i].w2prime = (consolidator.m_weighted_sh[i].wprime) / w2prime_denominator;
	       e_obj.weightedOutput(consolidator.m_weighted_sh[i].per1, consolidator.m_weighted_sh[i].per2, consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2, consolidator.m_weighted_sh[i].w2prime);
	     }
	   }

	     /*End weighted output*/



	     /*END TESTING AREA DEC 4th 2014*/
	     /*****************************
	     ******************************/
	     consolidator.ma_drop_str = "No of matches removed due to length of trimming by moving averages: " + NumberToString( removed2 );
	     consolidator.pie_drop_str = "No of matches removed due to percentage error: " + NumberToString( removed1 );
	     if(T_TRIM){
	     //  str = str+ " \n No of matches removed due hold out ped file checking: "+ NumberToString( removed3 );
	     }
	     //begin log output
	     std::string parameter_string_1 = "\n\n**********Parameters used in program**********\n";
	     e_obj.log(parameter_string_1);
	     e_obj.log(consolidator.emp_ma_thresh_str); //keep
	     e_obj.log(consolidator.emp_pie_thresh_str);//keep
	     parameter_string_1 = "**********************************************\n\n";
	     e_obj.log(parameter_string_1);
	     std::string total_count_str = "The total number of SH in the input file was: " + NumberToString(total_count);
	     e_obj.log(total_count_str);
	     e_obj.log(consolidator.consolidated_str);
	     e_obj.log(consolidator.initial_drop_str);
	     //  e_obj.log(ibg_str);
	     e_obj.log(consolidator.ma_drop_str);
	     e_obj.log(consolidator.pie_drop_str);
	     consolidator.final_sh_str = "Total number of SH that were not dropped is: " + NumberToString(not_removed);
	     e_obj.log(consolidator.final_sh_str);














}////ending of findTrueSimplePctErrors & error





































	//initiateErrorFinder(); //but bypass true calculations if empirical threshold is supplied


/*
	if( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 )
		{
			//std::cerr<<"entering into HoldOut Mode Matching"<<std::endl;
			           consolidator.performTrim(eCalculator,
							  WINDOW,MA_SNP_END,
							 MA_ERR_THRESHOLD_START,
							 MA_ERR_THRESHOLD_END,MIN_SNP,
							 MIN_CM,PCT_ERR_THRESHOLD, OPTION, HO_THRESHOLD, true  );
			consolidator.performTrim(eCalculator, WINDOW, MA_SNP_END,MA_THRESHOLD,MIN_SNP,MIN_CM,PCT_ERR_THRESHOLD,OPTION,HO_THRESHOLD,true,EMPIRICAL_MA_RESULT,EMPIRICAL_PIE_RESULT,EXTENDSNP);//<piyush> added the param EXTENDSNP for calculating moving window
			//std::cerr<<" Main Trim operation has completed "<< std::endl;
			//std::cerr<< " Hold out trim has completed" <<std::endl;
			if( (OPTION.compare( "finalOutput" ) == 0) || (OPTION.compare( "Full" ) == 0 ) )
				{
					consolidator.finalOutPut( eCalculator, MIN_CM, MIN_SNP );
				}
		}
	else
		{
			if( (OPTION.compare( "finalOutput" ) == 0) || (OPTION.compare( "Full" ) == 0 ) )
				{
					//cerr <<"DEBUG: ENTERING CONSOLIDATOR FINAL OUTPUT" << endl;
						consolidator.finalOutPut( eCalculator, MIN_CM, MIN_SNP );
					//cerr <<"DEBUG: EXITING CONSOLIDATOR FINAL OUTPUT" << endl;
				}
		}

*/

if (T_TRIM)
{
	std::cerr<<" Main Trim operation has completed "<< std::endl;
	std::cerr<< " Hold out trim has completed" <<std::endl;
	if( (OPTION.compare( "finalOutput" ) == 0) || (OPTION.compare( "Full" ) == 0 ) )
	{
		consolidator.finalOutPut( eCalculator, MIN_CM, MIN_SNP );
	}
}
else
{
	if( (OPTION.compare( "finalOutput" ) == 0) || (OPTION.compare( "Full" ) == 0 ) )
	{
	//cerr <<"DEBUG: ENTERING CONSOLIDATOR FINAL OUTPUT" << endl;
		consolidator.finalOutPut( eCalculator, MIN_CM, MIN_SNP );
	//cerr <<"DEBUG: EXITING CONSOLIDATOR FINAL OUTPUT" << endl;
	}

}




//---------------------------------------------------------------------------------


	time_t endTime;
	time (&endTime);
	/*Provide the log file with information about start and end times*/
	std::string str2 = " The program ended at: " +
	std::string( ctime ( &endTime ) );
	//remove trailing newline, for readability.
	str1.erase(std::remove(str1.begin(),str1.end(),'\n'),str1.end());
	str2.erase(std::remove(str2.begin(),str2.end(),'\n'),str2.end());
	str2 = str2 +  "           Total time ( in seconds): " + NumberToString( ( endTime - startTime ) );

	//add start time
	str2 = "\n\n" + str1 + "        " + str2;
	//eCalculator.log( str1 );
	eCalculator.log( str2 );

}




void ErrorFinderManager::initiateErrorFinder_PRE()
{
	eCalculator.readBmidFile(BMIDFILE);
	std::cerr<<"Reading bmid file completed"<<std::endl;
	eCalculator.readBsidFile(BSIDFILE);
	std::cerr<<"Reading bsid file completed"<<std::endl;
	eCalculator.readPedFile(PEDFILE, HO_MISSING);
	std::cerr<<"Reading ped file completed"<<std::endl;
}

void ErrorFinderManager::initiateErrorFinder_CURRENT(int pers_count)
{
	consolidator.readMatches(BMATCHFILE, pers_count, eCalculator, TRUESNP, TRUECM, EXTENDSNP,PEDFILE  );
	std::cerr<<"Reading bmatch file completed"<<std::endl;
	if( !(SNPWEIGHTFILE.empty()))
		{ //verify that this check actually works
			consolidator.readUserSuppliedSnpWeights(SNPWEIGHTFILE);
		}
	consolidator.performConsolidation(eCalculator,GAP,MIN_SNP,MIN_CM,EXTENDSNP);
	std::cerr<<"Consolidation completed"<<std::endl;

	if(( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 ) &&	ISMOL)
		{
			SITUATION_NO	=	1;
			T_ERROR	=	true;
			//T_TRIM = true;
			eCalculator.changeMapFile( HMAPFILE );
			eCalculator.readHPedFile( HPEDFILE, HO_MISSING );
			std::cerr<< " new map and ped File has been read" <<std::endl;
			std::cerr<< " calculating true percentage errors" <<std::endl;
			//consolidator.findTruePctErrors( eCalculator, MA_SNP_END, T_VAL, WINDOW,MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count, TRUESNP, TRUECM);//T_VAL = true
			CERRMESSAGE	=	" true percentage errors calculated ";
		}
	else if(( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 ) &&	!ISMOL)
		{
			SITUATION_NO	=	2;
			T_ERROR	=	true;
			//T_TRIM = true;
			eCalculator.changeMapFile( HMAPFILE );
			eCalculator.readHPedFile( HPEDFILE, HO_MISSING );
			std::cerr<< " new map and ped File has been read" <<std::endl;
			std::cerr<< " calculating true percentage errors" <<std::endl;
			//consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, T_VAL, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );//T_VAL = true
			CERRMESSAGE	=	" true percentage errors calculated ";
		}
	else if(	!( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 ) &&	ISMOL)
			{
				SITUATION_NO	=	3;
				T_ERROR	=	false;
				//T_TRIM = false;
				std::cerr<< " calculating true percentage errors" <<std::endl;
				//consolidator.findTruePctErrors( eCalculator, MA_SNP_END, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );
				CERRMESSAGE	=	" true hold out percentage errors calculated ";
			}
	else if(	!( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 ) &&	!ISMOL)
			{
				SITUATION_NO	=	4;
				T_ERROR	=	false;
				//T_TRIM = false;
				std::cerr<< " calculating true percentage errors" <<std::endl;
				//consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );
				CERRMESSAGE	=	" true hold out percentage errors calculated ";
			}
		else
			{

			}

		if (( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 ))
			{
				T_TRIM = true;
			}
		else
			{
				T_TRIM = false;
			}

		if (!SITUATION_NO)
			{
				std::cerr<<"Wrong Situation:HPED HMAP, ISMOL Error "<<std::endl;
				exit(0);
			}


}

/*
void ErrorFinderManager::initiateErrorFinder()
{

	//int pers_count=eCalculator.getNoOfPersons();
	//consolidator.readMatches(BMATCHFILE, pers_count, eCalculator, TRUESNP, TRUECM, EXTENDSNP,PEDFILE  );
	//std::cerr<<"Reading bmatch file completed"<<std::endl;

	if( !(SNPWEIGHTFILE.empty()))
		{ //verify that this check actually works
			consolidator.readUserSuppliedSnpWeights(SNPWEIGHTFILE);
		}
	//testing to see what happens when consolidation is removed
	//consolidator.performConsolidation(eCalculator,GAP,MIN_SNP,MIN_CM,EXTENDSNP);
	//std::cerr<<"Consolidation completed"<<std::endl;
	if( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 )
		{
			eCalculator.changeMapFile( HMAPFILE );
			eCalculator.readHPedFile( HPEDFILE, HO_MISSING );
			std::cerr<< " new map and ped File has read" <<std::endl;
			std::cerr<< " calculating true percentage errors" <<std::endl;
		if( ISMOL )
			{
				consolidator.findTruePctErrors( eCalculator, MA_SNP_END, true, WINDOW,MA_THRESHOLD, EMPIRICAL_MA_RESULT);
			}
		else
			{
				consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, true, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT );
			}
			std::cerr<< " true percentage errors calculated "<<std::endl;
		}
	else
		{
			std::cerr<< " calculating true percentage errors" <<std::endl;
			if( ISMOL )
				{
					consolidator.findTruePctErrors( eCalculator, MA_SNP_END, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT );
				}
			else
				{
					consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT );//<piyush>
				}
			std::cerr<< " true hold out percentage errors calculated "<<std::endl;

		}

}
*/

void ErrorFinderManager::displayError(std::string argv)
{
	std::cerr<<"these parameters are not allowed "<<wrongParam<<std::endl;
	std::cerr << "Usage: " << argv << " -bmatch [BMATCH FILE]  -bsid [BSID FILE] -bmid [BMID FILE] -reduced [min_snp] [min_cm] "
             <<" -ped-file [ped file] -window [window width to calculate moving averages] "
             <<" -gap [max gap to consolidate two matches]"
               <<" -pct-err-threshold [max percentage of errors in a match after the trim] OR -emp-pie-threshold"
               <<" -ma-threshold [specifies percentile to be drawn from trulyIBD data for MA calculations] OR -empirical-ma-threshold"
               <<" Note that if both -emp-pie-threshold and empirical-ma-threshold are supplied, then -trueSNP and -trueCM will be ignored"
              <<"-output.type [ must provide any of these. it can be "
               << "MovingAverages  or Error1 or Error2 or Error3 or ErrorRandom1 "
                << "or ErrorRandom2 or Error3 or ErrorRandom3 or Full "
                <<  "look at the description about how these works in wiki ]"
              << "(optional) -holdout-ped [new ped file path] -holdout-map [new map file] "
              << "-holdout-threshold [threshold to drop a match with new ped file ]"
              << " -holdout-missing [missing value representation in new ped file] "
              << " -log.file [log file name]"
              << " -trueCM [ true match maximum cm length] "
              << " - trueSNP [ true match SNP length]"
              << " -PIE.dist.length [ can be MOL or any cm distance length "
              << "please refer wiki for more details on how to use this option"
              << "-count.gap.errors [ TRUE or FALSE to include gap errors in errors count ]"
            << std::endl;


}
/*
void ErrorFinderManager::initiateErrorFinder()
{
        eCalculator.readBmidFile(BMIDFILE);
        std::cerr<<"Reading bmid file completed"<<std::endl;
        eCalculator.readBsidFile(BSIDFILE);
        std::cerr<<"Reading bsid file completed"<<std::endl;
        eCalculator.readPedFile(PEDFILE, HO_MISSING);
        std::cerr<<"Reading ped file completed"<<std::endl;
        int pers_count=eCalculator.getNoOfPersons();
        consolidator.readMatches(BMATCHFILE, pers_count, eCalculator, TRUESNP, TRUECM, EXTENDSNP,PEDFILE  );
        std::cerr<<"Reading bmatch file completed"<<std::endl;
        if( !(SNPWEIGHTFILE.empty())){ //verify that this check actually works
          consolidator.readUserSuppliedSnpWeights(SNPWEIGHTFILE);
        }
//testing to see what happens when consolidation is removed
        consolidator.performConsolidation(eCalculator,GAP,MIN_SNP,MIN_CM,EXTENDSNP);
        std::cerr<<"Consolidation completed"<<std::endl;
        if( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 )
       {
        	eCalculator.changeMapFile( HMAPFILE );
        	eCalculator.readHPedFile( HPEDFILE, HO_MISSING );
        	std::cerr<< " new map and ped File has read" <<std::endl;
        	std::cerr<< " calculating true percentage errors" <<std::endl;
        if( ISMOL )
        {
          consolidator.findTruePctErrors( eCalculator, MA_SNP_END, true, WINDOW,MA_THRESHOLD, EMPIRICAL_MA_RESULT);

        }
        else
        {
          consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, true, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT );
        }
        std::cerr<< " true percentage errors calculated "<<std::endl;
       }
       else
       {
    	   std::cerr<< " calculating true percentage errors" <<std::endl;
          if( ISMOL )
          {
            consolidator.findTruePctErrors( eCalculator, MA_SNP_END, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT );
          }
          else
          {
            consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT );//<piyush>
          }
          std::cerr<< " true hold out percentage errors calculated "<<std::endl;
	
      }   
  
}*/

