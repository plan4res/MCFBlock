/*--------------------------------------------------------------------------*/
/*----------------------------- File main.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small main() for testing MCFBlock. It just creates one and loads it from a
 * file, then prints it back to a file. It also de-serializes it,
 * serialize a copy out of it, and finally de-serialize the copy.
 *
 * Little more than a compilation check.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>

#include "MCFBlock.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- Main -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 if( argc < 2 ) {
  cerr << "Usage: " << argv[ 0 ]
       << " DIMACS_in [ DIMACS_out netCDF_out netCDF_out_2 ]" << endl;
  return( 1 );
  }

 // open input file in DIMACS format
 ifstream ProbFile( argv[ 1 ] );
 if( ! ProbFile.is_open() ) {
  cerr << "Error: cannot open file " << argv[ 1 ] << endl;
  return( 1 );
  }

 // load the MCFBlock from the file
 auto MCFB = dynamic_cast< MCFBlock * >( Block::new_Block( "MCFBlock" ) );
 ProbFile >> *MCFB;
 ProbFile.close();

 if( argc == 2 ) {  // just print a little and quit
  MCFB->set_verbosity( Block::low );
  cout << *MCFB;
  delete MCFB;
  return( 0 );
  }

 // open output DIMACS file
 ofstream OutFile( argv[ 2 ] );
 if( ! OutFile.is_open() ) {
  cerr << "Error: cannot open file " << argv[ 2 ] << endl;
  delete MCFB;
  return( 1 );
  }

 // save the MCFBLock there
 MCFB->set_verbosity( Block::complete );
 OutFile << *MCFB;

 if( argc == 3 ) {
  delete MCFB;
  return( 0 );
  }

 // serialize the MCFBlock (in a netCDF BlockFile)
 MCFB->Block::serialize( argv[ 3 ] , eBlockFile );

 // for better or for worse, delete the MCFBlock;
 delete MCFB;

 if( argc == 4 )
  return( 0 );

 // de-serialize the MCFBlock from the same file
 MCFB = dynamic_cast< MCFBlock * >( Block::deserialize( argv[ 3 ] ) );

 // now re-serialize it on a different netCDF BlockFile
 MCFB->Block::serialize( argv[ 4 ] , eBlockFile );

 // delete it for good
 delete MCFB;

 // all done
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File main.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
