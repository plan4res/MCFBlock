/*--------------------------------------------------------------------------*/
/*--------------------------- File MCFSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MCFSolver class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/* If any of the following macros is defined, then the corresponding
 * :MCFClass solver is included and the corresponding version of
 * MCFSolver<> is defined:
 *
 * - HAVE_CSCL2      for the CS2 class
 *
 * - HAVE_CPLEX      for the MCFCplex class
 *
 * - HAVE_MFSMX      for the MCFSimplex class
 *
 * - HAVE_MFZIB      for the MCFZIB class
 *
 * - HAVE_RELAX      for the RelaxIV class
 *
 * - HAVE_SPTRE      for the SPTree class; note that SPTree cannot solve
 *                   most MCF instances, except those with SPT structure.
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFSolver.h"

#ifdef HAVE_CSCL2
 #include "CS2.h"
#endif

#ifdef HAVE_CPLEX
 #include "MCFCplex.h"
#endif

#ifdef HAVE_MFSMX
 #include "MCFSimplex.h"
#endif

#ifdef HAVE_MFZIB
 #include "MCFZIB.h"
#endif

#ifdef HAVE_RELAX
 #include "RelaxIV.h"
#endif

#ifdef HAVE_SPTRE
 #include "SPTree.h"
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register the various MCFSolver< * > to the Solver factory

#ifdef HAVE_CSCL2
 SMSpp_insert_in_factory_cpp_0_t( MCFSolver< CS2 > );
#endif

#ifdef HAVE_CPLEX
 SMSpp_insert_in_factory_cpp_0_t( MCFSolver< MCFCplex > );
#endif

#ifdef HAVE_MFSMX
 SMSpp_insert_in_factory_cpp_0_t( MCFSolver< MCFSimplex > );
#endif

#ifdef HAVE_MFZIB
 SMSpp_insert_in_factory_cpp_0_t( MCFSolver< MCFZIB > );
#endif

#ifdef HAVE_RELAX
 SMSpp_insert_in_factory_cpp_0_t( MCFSolver< RelaxIV > );
#endif

#ifdef HAVE_SPTRE
 SMSpp_insert_in_factory_cpp_0_t( MCFSolver< SPTree > );
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register MCFSolverState to the State factory

SMSpp_insert_in_factory_cpp_0( MCFSolverState );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the various static maps

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* Managing parameters for MCFSimplex---------------------------------------*/
/*
 * MCFSimplex has the following extra parameters:
 *
 * - kAlgPrimal     parameter to set algorithm (Primal/Dual):
 * - kAlgPricing    parameter to set algorithm of pricing
 * - kNumCandList   parameter to set the number of candidate list for
 *                  Candidate List Pivot method
 * - kHotListSize   parameter to set the size of Hot List
 *
 * (plus, actually, a few about Quadratic MCF that are not relevent here).
 * These are all "int" parameters, hence the "double" versions only issue the
 * method of the base CDASolver class, and therefore need not be defined. */

#ifdef HAVE_MFSMX

/*--------------------------------------------------------------------------*/

template<>
const std::vector< int > MCFSolver< MCFSimplex >::Solver_2_MCFClass_int = {
 MCFClass::kMaxIter ,        // intMaxIter
 -1 ,                        // intMaxThread
 -1 ,                        // intEverykIt
 -1 ,                        // intMaxSol
 -1 ,                        // intLogVerb
 -1 ,                        // intMaxDSol
 MCFClass::kReopt ,          // intLastParCDAS
 MCFSimplex::kAlgPrimal ,
 MCFSimplex::kAlgPricing ,
 MCFSimplex::kNumCandList ,
 MCFSimplex::kHotListSize
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<>
const std::vector< int > MCFSolver< MCFSimplex >::Solver_2_MCFClass_dbl = {
 MCFClass::kMaxTime ,       // dblMaxTime
 -1 ,                       // dblEveryTTm
 -1 ,                       // dblRelAcc
 MCFClass::kEpsFlw ,        // dblAbsAcc
 -1 ,                       // dblUpCutOff
 -1 ,                       // dblLwCutOff
 -1 ,                       // dblRAccSol
 -1 ,                       // dblAAccSol
 -1 ,                       // dblFAccSol
 -1 ,                       // dblRAccDSol
 MCFClass::kEpsCst ,        // dblAAccDSol
 -1                         // dblFAccDSol
 };

/*--------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< MCFSimplex >::get_num_int_par( void ) const
{
 return( CDASolver::get_num_int_par() + 5 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< MCFSimplex >::get_num_dbl_par( void ) const { }

----------------------------------------------------------------------------*/

template<>
int MCFSolver< MCFSimplex >::get_dflt_int_par( idx_type par ) const
{
 static const std::array< int , 5 > my_dflt_int_par = { MCFClass::kYes ,
		  MCFClass::kYes , MCFSimplex::kCandidateListPivot , 0 , 0 };

 return( par >= intLastParCDAS ? my_dflt_int_par[ par - intLastParCDAS ]
	                       : CDASolver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
double MCFSolver< MCFSimplex >::get_dflt_dbl_par( const idx_type par )
 const { }

----------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< MCFSimplex >::int_par_str2idx(
					     const std::string & name ) const
{
 if( name == "kReopt" )
  return( intLastParCDAS );
 if( name == "kAlgPrimal" )
  return( intLastParCDAS + 1 );
 if( name == "kAlgPricing" )
  return( intLastParCDAS + 2 );
 if( name == "kNumCandList" )
  return( intLastParCDAS + 3 );
 if( name == "kHotListSize" )
  return( intLastParCDAS + 4 );

 return( CDASolver::dbl_par_str2idx( name ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< MCFSimplex >::dbl_par_str2idx(
					const std::string & name ) const { }

----------------------------------------------------------------------------*/

template<>
const std::string & MCFSolver< MCFSimplex >::int_par_idx2str( idx_type idx )
 const {
 static const std::array< std::string , 5 > my_int_pars_str = {
  "kReopt" , "kAlgPrimal" , "kAlgPricing" , "kNumCandList" , "kHotListSize"
  };

 return( idx >= intLastParCDAS ? my_int_pars_str[ idx - intLastParCDAS ]
	                       : CDASolver::int_par_idx2str( idx ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
const std::string & MCFSolver< MCFSimplex >::dbl_par_idx2str( idx_type idx )
 const { }
 */

#endif

/*--------------------------------------------------------------------------*/
/* Managing parameters for RelaxIV -----------------------------------------*/
/*
 * RelaxIV has the following extra parameters:
 *
 * - kAuction     the auction/shortest paths initialization is used
 *
 * These are all "int" parameters, hence the "double" versions only issue the
 * method of the base CDASolver class, and therefore need not be defined. */

#ifdef HAVE_RELAX

/*--------------------------------------------------------------------------*/

template<>
const std::vector< int > MCFSolver< RelaxIV >::Solver_2_MCFClass_int = {
 MCFClass::kMaxIter ,        // intMaxIter
 -1 ,                        // intMaxThread
 -1 ,                        // intEverykIt
 -1 ,                        // intMaxSol
 -1 ,                        // intLogVerb
 -1 ,                        // intMaxDSol
 MCFClass::kReopt ,          // intLastParCDAS
 RelaxIV::kAuction
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<>
const std::vector< int > MCFSolver< RelaxIV >::Solver_2_MCFClass_dbl = {
 MCFClass::kMaxTime ,       // dblMaxTime
 -1 ,                       // dblEveryTTm
 -1 ,                       // dblRelAcc
 MCFClass::kEpsFlw ,        // dblAbsAcc
 -1 ,                       // dblUpCutOff
 -1 ,                       // dblLwCutOff
 -1 ,                       // dblRAccSol
 -1 ,                       // dblAAccSol
 -1 ,                       // dblFAccSol
 -1 ,                       // dblRAccDSol
 MCFClass::kEpsCst ,        // dblAAccDSol
 -1                         // dblFAccDSol
 };

/*--------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< RelaxIV >::get_num_int_par( void ) const
{
 return( CDASolver::get_num_int_par() + 2 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< RelaxIV >::get_num_dbl_par( void ) const { }

----------------------------------------------------------------------------*/

template<>
int MCFSolver< RelaxIV >::get_dflt_int_par( const idx_type par ) const
{
 static const std::array< int , 2 > my_dflt_int_par = { MCFClass::kYes ,
						        MCFClass::kYes };

 return( par >= intLastParCDAS ? my_dflt_int_par[ par - intLastParCDAS ]
	                       : CDASolver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
double MCFSolver< RelaxIV >::get_dflt_dbl_par( idx_type par ) const { }

----------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< RelaxIV >::int_par_str2idx(
					     const std::string & name ) const
{
 if( name == "kReopt" )
  return( intLastParCDAS );
 if( name == "kAuction" )
  return( intLastParCDAS + 1 );

 return( CDASolver::dbl_par_str2idx( name ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< RelaxIV >::dbl_par_str2idx(
					 const std::string & name ) const { }

----------------------------------------------------------------------------*/

template<>
const std::string & MCFSolver< RelaxIV >::int_par_idx2str( idx_type idx )
 const {
 static const std::array< std::string , 2 > my_int_pars_str = { "kReopt" ,
								"kAuction" };

 return( idx >= intLastParCDAS ? my_int_pars_str[ idx - intLastParCDAS ]
	                       : CDASolver::int_par_idx2str( idx ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
const std::string & MCFSolver< RelaxIV >::dbl_par_idx2str( idx_type idx )
 const { }
 */

#endif

/*--------------------------------------------------------------------------*/
/* Managing parameters for CPLEX -------------------------------------------*/

#ifdef HAVE_CPLEX

/*--------------------------------------------------------------------------*/

template<>
const std::vector< int > MCFSolver< MCFCplex >::Solver_2_MCFClass_int = {
 MCFClass::kMaxIter ,        // intMaxIter
 -1 ,                        // intMaxThread
 -1 ,                        // intEverykIt
 -1 ,                        // intMaxSol
 -1 ,                        // intLogVerb
 -1 ,                        // intMaxDSol
 MCFClass::kReopt ,          // intLastParCDAS
 MCFCplex::kQPMethod
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<>
const std::vector< int > MCFSolver< MCFCplex >::Solver_2_MCFClass_dbl =
{
 MCFClass::kMaxTime ,       // dblMaxTime
 -1 ,                       // dblEveryTTm
 -1 ,                       // dblRelAcc
 MCFClass::kEpsFlw ,        // dblAbsAcc
 -1 ,                       // dblUpCutOff
 -1 ,                       // dblLwCutOff
 -1 ,                       // dblRAccSol
 -1 ,                       // dblAAccSol
 -1 ,                       // dblFAccSol
 -1 ,                       // dblRAccDSol
 MCFClass::kEpsCst ,        // dblAAccDSol
 -1                         // dblFAccDSol
 };

/*--------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< MCFCplex >::get_num_int_par( void ) const {
 return( CDASolver::get_num_int_par() + 2 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< MCFCplex >::get_num_dbl_par( void ) const { }

----------------------------------------------------------------------------*/

template<>
int MCFSolver< MCFCplex >::get_dflt_int_par( idx_type par ) const {
 static const std::array< int , 2 > my_dflt_int_par = { MCFClass::kYes ,
                                                        MCFClass::kYes };
 return( par >= intLastParCDAS ?
         my_dflt_int_par[ par - intLastParCDAS ] :
         CDASolver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 template<>
 double MCFSolver< MCFCplex >::get_dflt_dbl_par( idx_type par ) const { }

----------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< MCFCplex >::int_par_str2idx(
					   const std::string & name ) const {
 if( name == "kReopt" )
  return( intLastParCDAS );
 if( name == "kQPMethod" )
  return( kQPMethod + 1 );

 return( CDASolver::dbl_par_str2idx( name ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< MCFCplex >::dbl_par_str2idx(
                                         const std::string & name ) const { }

----------------------------------------------------------------------------*/

template<>
const std::string & MCFSolver< MCFCplex >::int_par_idx2str( idx_type idx )
 const {
 static const std::array< std::string , 2 > my_int_pars_str = { "kReopt" ,
                                                                "kQPMethod" };
 return( idx >= intLastParCDAS ?
         my_int_pars_str[ idx - intLastParCDAS ] :
         CDASolver::int_par_idx2str( idx ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
const std::string & MCFSolver< MCFCplex >::dbl_par_idx2str( idx_type idx )
 const { }
*/

#endif

/*--------------------------------------------------------------------------*/
/* Managing parameters for SPTree ------------------------------------------*/

#ifdef HAVE_SPTRE

/*--------------------------------------------------------------------------*/

template<>
const std::vector< int > MCFSolver< SPTree >::Solver_2_MCFClass_int = {
 MCFClass::kMaxIter,        // intMaxIter
 -1,                        // intMaxSol
 -1,                        // intLogVerb
 -1,                        // intMaxDSol
 MCFClass::kReopt,          // intLastParCDAS
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<>
const std::vector< int > MCFSolver< SPTree >::Solver_2_MCFClass_dbl = {
 MCFClass::kMaxTime,        // dblMaxTime
 -1,                        // dblRelAcc
 MCFClass::kEpsFlw,         // dblAbsAcc
 -1,                        // dblUpCutOff
 -1,                        // dblLwCutOff
 -1,                        // dblRAccSol
 -1,                        // dblAAccSol
 -1,                        // dblFAccSol
 -1,                        // dblRAccDSol
 MCFClass::kEpsCst,         // dblAAccDSol
 -1                         // dblFAccDSol
 };

/*--------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< SPTree >::get_num_int_par( void ) const {
 return( CDASolver::get_num_int_par() + 1 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< SPTree >::get_num_dbl_par( void ) const { }

----------------------------------------------------------------------------*/

template<>
int MCFSolver< SPTree >::get_dflt_int_par( idx_type par ) const {
 static const std::array< int , 1 > my_dflt_int_par = { MCFClass::kYes };

 return( par >= intLastParCDAS ?
         my_dflt_int_par[ par - intLastParCDAS ] :
         CDASolver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
double MCFSolver< SPTree >::get_dflt_dbl_par( idx_type par ) const { }

----------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< SPTree >::int_par_str2idx(
					  const std::string & name ) const {
 if( name == "kReopt" )
  return( intLastParCDAS );

 return( CDASolver::dbl_par_str2idx( name ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< SPTree >::dbl_par_str2idx(
                                         const std::string & name ) const { }

----------------------------------------------------------------------------*/

template<>
const std::string & MCFSolver< SPTree >::int_par_idx2str( idx_type idx )
 const {
 static const std::array< std::string , 1 > my_int_pars_str = { "kReopt" };

 return( idx >= intLastParCDAS ?
         my_int_pars_str[ idx - intLastParCDAS ] :
         CDASolver::int_par_idx2str( idx ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 template<>
 const std::string & MCFSolver< SPTree >::dbl_par_idx2str( idx_type idx )
 const { }
*/

#endif

/*--------------------------------------------------------------------------*/
/*----------------------- End File MCFSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
