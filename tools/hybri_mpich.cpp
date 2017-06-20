// and command for run mpirun -np 4  ./a.out A 0.4 0
//for rahaman /home/soumen/softwares/mpich-3.1.4/bin/mpic++ -std=c++1y hybri_boost.cpp -o hybri
// for desktop /usr/bin/mpic++ -std=c++1y -o hybri hybri_mpich.cpp 
// for rahamn  CC -std=c++1y -o hybri hybri_mpich.cpp 
/*  Created on 10/04/16
*   ITS for 2D square lattice hybridization function for IHM. 
*   Warning:: Please specify the Delta, t, t2, beta, GRID for your problem. please make sure they are correct one.
*   Warning:: if you want to use the previous U hybridization finction the set the previous U last iteration sigma to SigmaA.out.0 and SigmaB.out.0
*   Veriable:: A/B(subtitle type ),mu(chemical_potential),it(current iteration no) pass these argument with executable.
*/

#include <chrono>        // header for time
#include<iostream>
#include <math.h>       /* cos */
#include <vector>
#include<fstream>       // header file for input and output
#include<complex>
#include<string>
#include <iterator>
#include <sstream>      // header for converting int to string
#include <mpi.h>
#include <stdio.h>
//#include <boost/mpi/environment.hpp>
//#include <boost/mpi/communicator.hpp>
//#include <boost/mpi.hpp>
#define PI 3.14159265

//namespace mpi = boost::mpi;

class HYBRI {
  //private variable define
  long double t,t2,delta,beta;
  const static int dim = 2,N= 1025,GRID = 100;
  long double cos_[2*GRID] = {};
  

  public:
  //public veraible
  std::complex<long double> DAAU[N] = {};
  std::complex<long double> DAAD[N] = {};
  std::complex<long double> DBBU[N] = {};
  std::complex<long double> DBBD[N] = {};
  std::complex<long double> sigmaAU[N] = {};
  std::complex<long double> sigmaAD[N] = {};
  std::complex<long double> sigmaBU[N] = {};
  std::complex<long double> sigmaBD[N] = {};
  long double wn[N] = {};

  //constructor
   HYBRI(long double t_, long double t2_, long double delta_, long double beta_){
	t = t_;
	t2 = t2_;
	delta = delta_;
	beta = beta_;

	for (int i=0; i<N; i++){
		wn[i] = (2*i+1)*PI/beta;
		} //for

        for(int i = 0; i<2*GRID; i++){
	 	cos_[i] = cos( (i-GRID)*PI/GRID);
		}

	}//HYBRI()

  //Calculate Hybridization for A sublattice
  void hybridizationA(long double mu,int start, int end){
        std::cout<<"start----:"<<start<<"  end-------"<<end<<std::endl;
	for (int n =start; n<end; n++){   		//loop for matsubara frequency
           auto alpaAU =std::complex<long double>(0.0,wn[n]) + mu - delta - sigmaAU[n];
	   auto alpaBU = std::complex<long double>(0.0,wn[n]) + mu + delta - sigmaBU[n];
           auto alpaAD = std::complex<long double>(0.0,wn[n]) + mu - delta - sigmaAD[n];
	   auto alpaBD = std::complex<long double>(0.0,wn[n]) + mu + delta - sigmaBD[n];
	   std::complex<long double> GU = (0.0,0.0);
	   std::complex<long double> GD = (0.0,0.0);

	   for (int j =0; j<2*GRID; j++){                //look fo kx    
                  for (int k =0; k<2*GRID; k++){         //loop for ky
		     auto l = 4*t2*cos_[j]*cos_[k];
		     auto m = 2*t*(cos_[j] + cos_[k]);
		     GU = GU + 	(alpaBU + l)/((alpaBU + l)*(alpaAU + l) - m*m);
 		     GD = GD + 	(alpaBD + l)/((alpaBD + l)*(alpaAD + l) - m*m);		 
                  }//for (int k; k<2*GRID; k++)
           } // for (int j; j<2*GRID; j++)
          
          //normalise the grren function
          long double normalise = 1.0/(4.0*GRID*GRID);
	  GU = GU*normalise;
	  GD = GD*normalise;

        // Delta function for wach matsubara frequency
	//std::cout<<"GU:"<<n<<GU<<std::endl;
	DAAU[n] = alpaAU - std::complex<long double>(1.0,0.0)/GU;
	DAAD[n] = alpaAD - std::complex<long double>(1.0,0.0)/GD;
	}//for (int n; n<N; n++)\

	/*
        MPI_Barrier(MPI_COMM_WORLD) ;
        for (int i = 0; i<N; i++){
        long double sumr;
        long double re = std::real(DAAD[i]);
        mpi::reduce(c, re, sumr, std::plus<long double>(),0);
        boost::mpi::reduce(world, hf.DAAD[i], hf.DAAD[i], std::plus<std::complex<long double>>(),0);
         }*/
      
	}//hybridizationA()

   //Module for writing hybr. fn. for A sub. lat. in Delta.inp 
   void write_hybA(){
	std::ofstream DeltaA;
        DeltaA.open("Delta_mpich.inp"); 
	 for (int i =0; i<N; i++){
		std::real(DAAU[i]);
		DeltaA<<wn[i]<<"   "<<std::real(DAAU[i])<<"   "<<std::imag(DAAU[i])<<"   "<<std::real(DAAD[i])<<"   "<<std::imag(DAAD[i])<<std::endl;
		}
	DeltaA.close();
	}//write_hybA()

   //Calculate Hybridization forB sublattice
   void hybridizationB(long double mu, int start, int end){

	for (int n =start; n<end; n++){                          //matsubara frequency
           auto alpaAU = std::complex<long double>(0.0,wn[n]) + mu - delta - sigmaAU[n];
	   auto alpaBU = std::complex<long double>(0.0,wn[n]) + mu + delta - sigmaBU[n];
           auto alpaAD = std::complex<long double>(0.0,wn[n]) + mu - delta - sigmaAD[n];
	   auto alpaBD = std::complex<long double>(0.0,wn[n]) + mu + delta - sigmaBD[n];
	   std::complex<long double> GU = (0.0,0.0);
	   std::complex<long double> GD = (0.0,0.0);

	   for (int j =0; j<2*GRID; j++){                        //kx sum
                  for (int k =0; k<2*GRID; k++){                 //ky sum
		     auto l = 4*t2*cos_[j]*cos_[k];
		     auto m = 2*t*(cos_[j] + cos_[k]);
		     GU = GU + 	(alpaAU + l)/((alpaBU + l)*(alpaAU + l) - m*m);
 		     GD = GD + 	(alpaAD + l)/((alpaBD + l)*(alpaAD + l) - m*m);		 
                  }//for (int k; k<2*GRID; k++)
           } // for (int j; j<2*GRID; j++)
          
          long double normalise = 1.0/(4.0*GRID*GRID);
	  GU = GU*normalise;
	  GD = GD*normalise;

	  //std::cout<<"GU:"<<n<<GU<<std::endl;
	DBBU[n] = alpaBU - std::complex<long double>(1.0,0.0)/GU;
	DBBD[n] = alpaBD - std::complex<long double>(1.0,0.0)/GD;
	}//for (int n; n<N; n++)
	}//hybridizationA()
  
   //Module for writing hybr. fn. for B sub. lat. in Delta.inp 
   void write_hybB(){
	std::ofstream DeltaB;
        DeltaB.open("Delta.inp"); 
	 for (int i =0; i<N; i++){
		DeltaB<<wn[i]<<"   "<<std::real(DBBU[i])<<"   "<<std::imag(DBBU[i])<<"   "<<std::real(DBBD[i])<<"   "<<std::imag(DBBD[i])<<std::endl;
		}
	DeltaB.close();
	}//write_hybB()
  
};//class HYBRI
    
 
//##############################################
//        MAIN FUNCTION
//############################################## 

int main(int argc , char * argv[]){ // A/B(subtitle type ),mu(chemical_potential),it(current iteration no) pass these argument with executable
 
  //variable define
  //mpi::environment env;
  //mpi::communicator world;

  long double t,t2,delta,U,mu,beta;
  int dim,GRID,it;
  int start, end,rank,size; // these variable are in a given processor from where to where execution will go.

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //initialising Veriables
  const int N = 1025;
  t = 0.5;
  t2 = 0.0;
  delta = 0.5;
  beta = 50.0;

  // intialising the requency range for different processor;
  //size = world.size();
  int step = int(N/size);
  start = rank*step;
  if ( rank == size-1) end = N;
  else end = (rank+1)*step ;

  // simulation starting time.
  auto start_time = std::chrono::high_resolution_clock::now();

  //creating class object
  HYBRI hf(t,t2,delta,beta);

  //taking the chemical potential and current iteration no. through argument
  std::stringstream convert_mu(argv[2]);
  std::stringstream convert_it(argv[3]);
  convert_mu >> mu;
  convert_it >> it; //it is current iteration no in DMFT needed to choose which self energy to take.

  std::string fileA; //veriable for SigA.out.it-1
  std::string fileB; //veriable for SigB.out.it-1
  double long wn; // this veraible will not be used naywhere it is only because to read verable.
  if (*argv[1] == 'A'){
        if (rank==0) std::cout<<"A sublatice hybridization is being created."<<std::endl;
       
       //choosing SigA and SigB from previous iteration or previos U value
       if (it ==0){
         fileA = "SigA.out."+std::to_string(0);  //from previous U value:: save Prev U Sig in Sig?.out.0
	 fileB = "SigB.out."+std::to_string(0);
	}
       else{
         fileA = "SigA.out."+std::to_string(it-1);  // from previous iteration
	 fileB = "SigB.out."+std::to_string(it-1);
	}

       //READING SigA.out.it-1
       std::ifstream SigmaA;
       SigmaA.open(fileA);
       long double SAUR,SAUI,SADR,SADI;
       if (SigmaA.is_open())
	{
		if (rank==0) std::cout<<"Reading Data From File::"<<fileA<<std::endl;
		std::string line;
		getline(SigmaA,line);
		for (int i = 0; i<N; i++){
			SigmaA>>wn>>SAUR>>SAUI>>SADR>>SADI;
                        //std::cout<<wn<<SAUR<<SAUI<<SADR<<SADI<<std::endl;
                        hf.sigmaAU[i] = std::complex<long double>(SAUR,SAUI);
			hf.sigmaAD[i] = std::complex<long double>(SADR,SADI);
			}//for (int i = 0; i<1024; i++)
	}//if (SigmaA.is_open())
       else{
	 if (rank==0) std::cout<<"file:"<<fileA<<"in not in current location called for it:"<<it<<std::endl;
	exit(1);
	}

       //READING SigB.out.it-1
       std::ifstream SigmaB;
       SigmaB.open(fileB);
       long double SBUR,SBUI,SBDR,SBDI;
       if (SigmaB.is_open())
	{       
		
	if (rank==0)	std::cout<<"Reading Data From File::"<<fileB<<std::endl;
		std::string line;
		getline(SigmaB,line);
		for (int i = 0; i<N; i++){
			SigmaB>>wn>>SBUR>>SBUI>>SBDR>>SBDI;
                        //std::cout<<wn<<SAUR<<SAUI<<SADR<<SADI<<std::endl;
                        hf.sigmaBU[i] = std::complex<long double>(SBUR,SBUI);
			hf.sigmaBD[i] = std::complex<long double>(SBDR,SBDI);
			}//for (int i = 0; i<1024; i++)
	}//if (SigmaB.is_open())
	else{
	 if (rank==0) std::cout<<"file:"<<fileB<<"in not in current location called for it:"<<it<<std::endl;
	exit(1);
	}

       // calling for Hybridization function calculation
       hf.hybridizationA(mu,start,end); 
       MPI_Barrier(MPI_COMM_WORLD) ;
       
       
       //MPI massage transfer
       for (int i = 0; i<N; i++){
       long double sumrU,sumrD,sumiU,sumiD;
       long double reU= std::real(hf.DAAU[i]);
       long double imU= std::imag(hf.DAAU[i]);
       long double reD= std::real(hf.DAAD[i]);
       long double imD= std::imag(hf.DAAD[i]);	
       MPI_Reduce(&reU, &sumrU, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       MPI_Reduce(&imU, &sumiU, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       MPI_Reduce(&reD, &sumrD, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       MPI_Reduce(&imD, &sumiD, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

       if (rank == 0) {
        hf.DAAU[i] = std::complex<long double>(sumrU,sumiU);
	hf.DAAD[i] = std::complex<long double>(sumrD,sumiD);
           }//(world.rank() == 0)
       } //for (int i = 0; i<N; i++)
      
       // calling for creating Delta.inp
       MPI_Barrier(MPI_COMM_WORLD) ;
       if (rank==0){
	 std::cout<<"writing Delta.inp for A"<<std::endl;
	 hf.write_hybA();
	}
       
     }//if (*argv[1] == 'A')

   if (*argv[1] == 'B'){
        if (rank==0) std::cout<<"B sublatice hybridization is being created."<<std::endl;
       
       //choosing SigA and SigB from previous iteration or previos U value
       if (it ==0){
         fileA = "SigA.out."+std::to_string(0); //from previous U value:: save Prev U Sig in Sig?.out.0
	 fileB = "SigB.out."+std::to_string(0);
	}
       else{
         fileA = "SigA.out."+std::to_string(it); // from previous iteration
	 fileB = "SigB.out."+std::to_string(it-1);
	}

       //READING SigA.out.it-1
       std::ifstream SigmaA;
       SigmaA.open(fileA);
       long double SAUR,SAUI,SADR,SADI;
       if (SigmaA.is_open())
	{
		 if (rank==0) std::cout<<"Reading Data From File::"<<fileA<<std::endl;
		std::string line;
		getline(SigmaA,line);
		for (int i = 0; i<N; i++){
			SigmaA>>wn>>SAUR>>SAUI>>SADR>>SADI;
                        //std::cout<<wn<<SAUR<<SAUI<<SADR<<SADI<<std::endl;
                        hf.sigmaAU[i] = std::complex<long double>(SAUR,SAUI);
			hf.sigmaAD[i] = std::complex<long double>(SADR,SADI);
			}//for (int i = 0; i<1024; i++)
	}//if (SigmaA.is_open())
	else{
	 if (rank==0) std::cout<<"file:"<<fileA<<"in not in current location called for it:"<<it<<std::endl;
	exit(1);
	}

       //READING SigB.out.it-1
       std::ifstream SigmaB;
       SigmaB.open(fileB);
       long double SBUR,SBUI,SBDR,SBDI;
       if (SigmaB.is_open())
	{       
		
		 if (rank==0) std::cout<<"Reading Data From File::"<<fileB<<std::endl;
		std::string line;
		getline(SigmaB,line);
		for (int i = 0; i<N; i++){
			SigmaB>>wn>>SBUR>>SBUI>>SBDR>>SBDI;
                        //std::cout<<wn<<SAUR<<SAUI<<SADR<<SADI<<std::endl;
                        hf.sigmaBU[i] = std::complex<long double>(SBUR,SBUI);
			hf.sigmaBD[i] = std::complex<long double>(SBDR,SBDI);
			}//for (int i = 0; i<1024; i++)
	}//if (SigmaB.is_open())
	else{
	 if (rank==0) std::cout<<"file:"<<fileB<<"in not in current location called for it:"<<it<<std::endl;
	exit(1);
	}
	
     // calling for Hybridization function calculation
       hf.hybridizationB(mu,start,end); 

      MPI_Barrier(MPI_COMM_WORLD) ;
       //MPI massage transfer
       for (int i = 0; i<N; i++){
       long double sumrU,sumrD,sumiU,sumiD;
       long double reU= std::real(hf.DBBU[i]);
       long double imU= std::imag(hf.DBBU[i]);
       long double reD= std::real(hf.DBBD[i]);
       long double imD= std::imag(hf.DBBD[i]);	
       // MPI_Allreduce( in, out, count, MPI_INT, MPI_SUM, MPI_COMM_WORLD ); where count is array size, if its a number then 1
       //MPI_Reduce(&local_sum, &global_sum, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       //mpi::reduce(world, reU, sumrU, std::plus<long double>(),0);
       //mpi::reduce(world, imU, sumiU, std::plus<long double>(),0);
       //mpi::reduce(world, reD, sumrD, std::plus<long double>(),0);
       //mpi::reduce(world, imD, sumiD, std::plus<long double>(),0);

       MPI_Reduce(&reU, &sumrU, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       MPI_Reduce(&imU, &sumiU, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       MPI_Reduce(&reD, &sumrD, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       MPI_Reduce(&imD, &sumiD, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

       if (rank == 0) {
        hf.DBBU[i] = std::complex<long double>(sumrU,sumiU);
	hf.DBBD[i] = std::complex<long double>(sumrD,sumiD);
           }//(world.rank() == 0)
       } //for (int i = 0; i<N; i++)

      //printing delta in file
      if (rank==0) {
	std::cout<<"writing Delta.inp for B"<<std::endl;
        hf.write_hybB();
	   }

	}//if (*argv[1] == 'B')
  
  auto end_time = std::chrono::high_resolution_clock::now();
  
 //if(world.rank() == 0) std::cout<<"time taken in:"<<size<<"procesor is:"<<start<<":"<<end<<"second"<<std::endl;
 if (rank==0) std::cout << "tame taken" << std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count() << std::endl ;
  MPI_Finalize();
}//main
  
//TODO check how one would do the mpi::reduce for array of complex nomber.   

