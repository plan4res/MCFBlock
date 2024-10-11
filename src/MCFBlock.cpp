/*--------------------------------------------------------------------------*/
/*-------------------------- File MCFBlock.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MCFBlock class.
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
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFBlock.h"

#include <iomanip>

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 #define CHECK_DS 0
 /* Perform long and costly checks on the data structures representing the
  * abstract and the physical representations agree. */
#else
 #define CHECK_DS 0
 // never change this
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- TYPES -----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;
using c_Index = Block::c_Index;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset;

using FNumber = MCFBlock::FNumber;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr auto dNAN = std::numeric_limits< double >::quiet_NaN();

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

// in the DIMACS format, comment lines start with 'c'

static std::istream & eatDMXcomments( std::istream& is )
{
 for(;;) {
  is >> std::ws;  // skip whitespaces
  if( is.peek() == is.widen( 'c' ) )
   // a comment: skip the rest of line and move to next
   is.ignore( std::numeric_limits< std::streamsize >::max() , is.widen( '\n' ) );
  else
   break;
  }

 return( is );
 }

/*--------------------------------------------------------------------------*/

static void print_UB( std::ostream & os , FNumber ub )
{
 if( ub == Inf< FNumber >() )
  os << "+Inf";
 else
  os << ub;
 }

/*--------------------------------------------------------------------------*/

static FNumber read_UB( std::istream & iStrm )
{
 iStrm >> eatcomments;
 int c = iStrm.peek();
 if( ! iStrm )
  throw( std::invalid_argument( "error reading the input stream" ) );
  
 if( ( c != 'I' ) && ( c != 'i' ) ) {
  FNumber res;
  iStrm >> res;
  if( ! iStrm )
   throw( std::invalid_argument( "error reading the input stream" ) );
  return( res );
  }

 do { c = iStrm.get(); c = iStrm.peek();
      if( ! iStrm )
       throw( std::invalid_argument( "error reading the input stream" ) );

  } while( ( c != iStrm.widen( ' ' ) ) &&
	   ( c != iStrm.widen( '\n' ) ) &&
	   ( c != iStrm.widen( '\t' ) ) );

 return( Inf< FNumber >() );
 }

/*--------------------------------------------------------------------------*/
// returns the number of elements where two vectors differ

template< typename T >
static Index countdiff( T beg , T end , T cmp )
{
 Index ndiff = 0;
 for( ; beg != end ; )
  if( *(beg++) != *(cmp++) )
   ndiff++;

 return( ndiff );
 }

/*--------------------------------------------------------------------------*/
// returns true if two vectors differ, one of them being given as a base
// vector and a subset of indices

template< typename T >
static bool is_equal( std::vector< T > & vec , c_Subset & nms ,
		      typename std::vector< T >::const_iterator cmp ,
		      Index n_max )
{
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *(cmp++) )
   return( false );
  }

 return( true );
 }

/*--------------------------------------------------------------------------*/
// returns the number of elements where two vectors differ, one of them
// being given as a base vector and a subset of indices

template< typename T >
static Index countdiff( std::vector< T > & vec , c_Subset & nms ,
			typename std::vector< T >::const_iterator cmp ,
			Index n_max )
{
 Index ndiff = 0;
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *(cmp++) )
   ndiff++;
  }

 return( ndiff );
 }

/*--------------------------------------------------------------------------*/
// copys one vector to a given subset of another

template< typename T >
static void copyidx( std::vector< T > & vec , c_Subset & nms ,
		     typename std::vector< T >::const_iterator cpy )
{
 for( auto nm : nms )
  vec[ nm ] = *(cpy++);
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register MCFBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( MCFBlock );

// register MCFSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( MCFSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF MCFBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void MCFBlock::load( Index n , Index m , c_Subset & pEn , c_Subset & pSn ,
		     c_Vec_FNumber & pU , c_Vec_CNumber & pC ,
		     c_Vec_FNumber & pB , Index dn , Index dm ,
		     Index mdn , Index mdm )
{
 // sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( pSn.size() < m )
  throw( std::invalid_argument( "pSn too small" ) );

 if( pEn.size() < m )
  throw( std::invalid_argument( "pEn too small" ) );

 if( ( pC.size() > 0 ) && ( pC.size() < m ) )
  throw( std::invalid_argument( "pC nonempty but too small" ) );

 if( ( pU.size() > 0 ) && ( pU.size() < m ) )
  throw( std::invalid_argument( "pU nonempty but too small" ) );

 if( ( pB.size() > 0 ) && ( pB.size() < n ) )
  throw( std::invalid_argument( "pB nonempty but too small" ) );

 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( MaxNNodes || get_MaxNArcs() )
  guts_of_destructor();
		   
 // copy over problem data - - - - - - - - - - - - - - - - - - - - - - - - - -

 NNodes = n;
 NArcs = m;
 MaxNNodes = NNodes + ( mdn > dn ? mdn - dn : 0 );
 c_Index MaxNArcs = NArcs + ( mdm > dm ? mdm - dm : 0 );
 NStaticNodes = dn > n ? 0 : n - dn;
 NStaticArcs = dm > NArcs ? 0 : NArcs - dm;

 SN.resize( MaxNArcs );
 if( std::any_of( pSn.begin() , pSn.begin() + m ,
		  [ n ]( c_Index sn ) { return( ( sn < 1 ) || ( sn > n ) ); }
		  ) )
  throw( std::invalid_argument( "wrong starting node" ) );
 std::copy( pSn.begin() , pSn.begin() + m , SN.begin() );

 EN.resize( MaxNArcs );
 if( std::any_of( pEn.begin() , pEn.begin() + m ,
		  [ n ]( c_Index en ) { return( ( en < 1 ) || ( en > n ) ); }
		  ) )
  throw( std::invalid_argument( "wrong ending node" ) );
 std::copy( pEn.begin() , pEn.begin() + m , EN.begin() );

 if( pC.empty() )
  C.assign( MaxNArcs , 0 );
 else {
  C.resize( MaxNArcs );
  std::copy( pC.begin() , pC.begin() + m , C.begin() );
  }

 if( std::any_of( pU.begin() , pU.begin() + m ,
		  []( c_FNumber ui ) { return( ui < Inf< FNumber >() ); } ) ) {
  U.resize( MaxNArcs );
  std::copy( pU.begin() , pU.begin() + m , U.begin() );
  }
 else
  U.clear();

 if( std::any_of( pB.begin() , pB.begin() + n ,
		  []( c_FNumber bi ) { return( bi != 0 ); } ) ) {
  B.resize( MaxNNodes );
  std::copy( pB.begin() , pB.begin() + n , B.begin() );
  }
 else
  B.clear();

 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // the arcs whose cost is infinite have to be closed
 // in addition the cost has to be set to 0 - - - - - - - - - - - - - - - - -

 for( Index j = 0 ; j < C.size() ; ++j )
  if( C[ j ] >= Inf< CNumber >() ) {
   close_arc( j , eNoMod , eNoMod );
   C[ j ] = 0;
   }

 // throw Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( MCFBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::load( std::istream & input , char frmt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( MaxNNodes || get_MaxNArcs() )
  guts_of_destructor();

 // read first non-comment line - - - - - - - - - - - - - - - - - - - - - - -

 char c;
 if( ! ( input >> eatDMXcomments >> c ) )
  throw( std::invalid_argument( "error reading the input stream" ) );

 if( c != 'p' )
  throw( std::invalid_argument( "format error in the input stream" ) );

 input >> eatDMXcomments;
 input.ignore( 3 , ' ' );  // skip "min"

 if( ! ( input >> eatDMXcomments >> NNodes ) )
  throw( std::invalid_argument( "LoadDMX: error reading number of nodes" ) );

 Index tm;
 if( ! ( input >> eatDMXcomments >> NArcs ) )
  throw( std::invalid_argument( "LoadDMX: error reading number of arcs" ) );

 // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SN.resize( NArcs );
 EN.resize( NArcs );
 C.assign( NArcs , 0 );
 U.assign( NArcs , Inf< FNumber >() );
 B.assign( NNodes , 0 );

 NStaticNodes = MaxNNodes = NNodes;
 NStaticArcs = NArcs;

 // read problem data - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index i = 0;  // arc counter
 for(;;) {
  if( ! ( input >> eatDMXcomments >> c ) )  // read next descriptor
   break;                                   // if none, end

  switch( c ) {
   case( 'n' ):  // description of a node
    Index j;
    if( ! ( input >> j ) )
     throw( std::invalid_argument( "error reading node name" ) );

    if( ( j < 1 ) || ( j > NNodes ) )
     throw( std::invalid_argument( "invalid node name" ) );

    FNumber Dfctj;
    if( ! ( input >> Dfctj ) )
     throw( std::invalid_argument( "error reading deficit" ) );

    B[ j - 1 ] -= Dfctj;
    break;

   case( 'a' ):  // description of an arc
    if( i == NArcs )
     throw( std::invalid_argument( "too many arc descriptors" ) );

    if( ! ( input >> SN[ i ] ) )
     throw( std::invalid_argument( "error reading start node" ) );

    if( ( SN[ i ] < 1 ) || ( SN[ i ] > NNodes ) )
     throw( std::invalid_argument( "invalid start node" ) );

    if( ! ( input >> EN[ i ] ) )
     throw( std::invalid_argument( "error reading end node" ) );

    if( ( EN[ i ] < 1 ) || ( EN[ i ] > NNodes ) )
     throw( std::invalid_argument( "LoadDMX: invalid end node" ) );

    if( SN[ i ] == EN[ i ] )
     throw( std::invalid_argument( "self-loops not permitted" ) );

    FNumber LB;
    if( ! ( input >> LB ) )
     throw( std::invalid_argument( "error reading lower bound" ) );

    U[ i ] = read_UB( input );

    if( ! ( input >> C[ i ] ) )
     throw( std::invalid_argument( "error reading arc cost" ) );

    if( U[ i ] < LB )
     throw( std::invalid_argument( "lower bound > upper bound" ) );

    if( LB > 0 ) {
     if( U[ i ] < Inf< MCFBlock::FNumber >() )
      U[ i ] -= LB;
     B[ SN[ i ] - 1 ] += LB;
     B[ EN[ i ] - 1 ] -= LB;
     }
    i++;
    break; 

   default:  // invalid code- - - - - - - - - - - - - - - - - - - - - - - - -
    throw( std::invalid_argument( "invalid DMX code" ) );

   }  // end( switch( c ) )
  }  // end( for( ever ) )

 if( i < NArcs )
  throw( std::invalid_argument( "too few arc descriptors" ) );

 f_cond_lower = dNAN;  // reset conditional bounds

 // simplify out the deta structures- - - - - - - - - - - - - - - - - - - - -

 if( std::all_of( B.begin() , B.end() ,
		  []( c_FNumber bi ) { return( bi == 0 ); } ) )
  B.clear();

 if( std::all_of( U.begin() , U.end() ,
		  []( c_FNumber ui ) { return( ui == Inf< FNumber >() ); } ) )
  U.clear();
 
 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // the arcs whose cost is infinite have to be closed
 // in addition the cost has to be set to 0 - - - - - - - - - - - - - - - - -

 for( Index j = 0 ; j < C.size() ; ++j )
  if( C[ j ] >= Inf< double >() ) {
   close_arc( j , eNoMod , eNoMod );
   C[ j ] = 0;
   }

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( MCFBlock::load( istream ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::deserialize( const netCDF::NcGroup & group )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( MaxNNodes || get_MaxNArcs() )
  guts_of_destructor();
		   
 // read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 netCDF::NcDim nn = group.getDim( "NNodes" );
 if( nn.isNull() )
  throw( std::logic_error( "NNodes dimension is required" ) );
 NNodes = nn.getSize();

 netCDF::NcDim na = group.getDim( "NArcs" );
 if( na.isNull() )
  throw( std::logic_error( "NArcs dimension is required" ) );
 NArcs = na.getSize();

 Index DynNNodes = 0;
 netCDF::NcDim dn = group.getDim( "DynNNodes" );
 if( dn.isNull() )
  NStaticNodes = NNodes;
 else {
  DynNNodes = dn.getSize();
  NStaticNodes = DynNNodes > NNodes ? 0 : NNodes - DynNNodes;
  }

 Index DynNArcs = 0;
 netCDF::NcDim dm = group.getDim( "DynNArcs" );
 if( dm.isNull() )
  NStaticArcs = NArcs;
 else {
  DynNArcs = dm.getSize();
  NStaticArcs = DynNArcs > NArcs ? 0 : NArcs - DynNArcs;
  }

 MaxNNodes = NNodes;
 netCDF::NcDim mdn = group.getDim( "MaxDynNNodes" );
 if( ( ! mdn.isNull() ) && ( mdn.getSize() > DynNNodes ) )
  MaxNNodes += mdn.getSize() - DynNNodes;

 Index MaxNArcs = NArcs;
 netCDF::NcDim mdm = group.getDim( "MaxDynNArcs" );
 if( ( ! mdm.isNull() ) && ( mdm.getSize() > DynNArcs ) )
  MaxNArcs += mdm.getSize() - DynNArcs;
 
 netCDF::NcVar sn = group.getVar( "SN" );
 if( sn.isNull() )
  throw( std::logic_error( "Starting Nodes not found" ) );

 SN.resize( MaxNArcs );

 sn.getVar( SN.data() );

 netCDF::NcVar en = group.getVar( "EN" );
 if( en.isNull() )
  throw( std::logic_error( "Ending Nodes not found" ) );

 EN.resize( MaxNArcs );
 en.getVar( EN.data() );

 netCDF::NcVar cst = group.getVar( "C" );
 if( ! cst.isNull() ) {
  C.resize( MaxNArcs );
  cst.getVar( C.data() );
  }
 else
  C.assign( MaxNArcs , 0 );

 netCDF::NcVar cap = group.getVar( "U" );
 if( ! cap.isNull() ) {
  U.resize( MaxNArcs );
  cap.getVar( U.data() );
  if( std::all_of( U.begin() , U.begin() + NArcs ,
		   []( c_FNumber ui ) { return( ui == Inf< FNumber >() ); } ) )
   U.clear();
  }

 netCDF::NcVar dfc = group.getVar( "B" );
 if( ! dfc.isNull() ) {
  B.resize( MaxNNodes );
  std::vector< size_t > countn = { NNodes };
  dfc.getVar( B.data() );
  if( std::all_of( B.begin() , B.begin() + NNodes ,
		   []( c_FNumber bi ) { return( bi == 0 ); } ) )
   B.clear();
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // the arcs whose cost is infinite have to be closed
 // in addition the cost has to be set to 0 - - - - - - - - - - - - - - - - -

 for( Index j = 0 ; j < C.size() ; ++j )
  if( C[ j ] >= Inf< CNumber >() ) {
   close_arc( j , eNoMod , eNoMod );
   C[ j ] = 0;
   }

 // call the method of Block- - - - - - - - - - - - - - - - - - - - - - - - -
 // inside this the NBModification, the "nuclear option",  is issued

 Block::deserialize( group );

 }  // end( MCFBlock::deserialize )

/*--------------------------------------------------------------------------*/

void MCFBlock::generate_abstract_variables( Configuration *stvv )
{
 if( AR & HasVar )  // the variables are there already
  return;           // nothing to do

 if( HasStaticX() ) {
  x.resize( get_NStaticArcs() );
  for( auto & var : x )
   var.is_positive( true , eNoBlck );

  add_static_variable( x );
  }

 if( MayHaveDynX() ) {
  dx.resize( get_NArcs() - get_NStaticArcs() );
  for( auto & var : dx )
   var.is_positive( true , eNoBlck );

  add_dynamic_variable( dx );
  } 

 AR |= HasVar;

 }  // end( MCFBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void MCFBlock::generate_abstract_constraints( Configuration *stcc )
{
 if( ! ( AR & HasFlw ) ) {
  // count number of nonzeroes in each constraint, i.e., #FS( i ) + #BS( i )
  Subset count( get_NNodes() );
 
  for( Index i = 0 ; i < get_NStaticArcs() ; ++i ) {
   count[ SN[ i ] - 1 ]++;
   count[ EN[ i ] - 1 ]++;
   }

  for( Index i = NStaticArcs ; i < get_NArcs() ; ++i )
   if( ! is_deleted( i ) ) {
    count[ SN[ i ] - 1 ]++;
    count[ EN[ i ] - 1 ]++;
    }

  // initialize the vectors of coefficients, and reset count[]
  std::vector< LinearFunction::v_coeff_pair > coeffs( get_NNodes() );

  for( Index i = 0 ; i < get_NNodes() ; ++i ) {
   coeffs[ i ].resize( count[ i ] );
   count[ i ] = 0;
   }

  // construct the vector of coefficients, static phase
  Index i = 0;
  for( ; i < get_NStaticArcs() ; ++i ) {
   coeffs[ SN[ i ] - 1 ][ count[ SN[ i ] - 1 ]++ ] =
                                    std::make_pair( &x[ i ] , double( -1 ) );
   coeffs[ EN[ i ] - 1 ][ count[ EN[ i ] - 1 ]++ ] =
                                    std::make_pair( &x[ i ] , double( 1 ) );
   }

  // construct the vector of coefficients, dynamic phase
  if( MayHaveDynX() )
   for( auto dxi = dx.begin() ; i < get_NArcs() ; ++i , ++dxi )
    if( ! is_deleted( i ) ) {
     coeffs[ SN[ i ] - 1 ][ count[ SN[ i ] - 1 ]++ ] =
                                    std::make_pair( &(*dxi) , double( -1 ) );
     coeffs[ EN[ i ] - 1 ][ count[ EN[ i ] - 1 ]++ ] =
                                    std::make_pair( &(*dxi) , double( 1 ) );
     }

  // generate the node-arc incidence matrix - - - - - - - - - - - - - - - - -
  // each constraint is an equality, i.e., LHS = RHS = B[ i ]

  // static part
  if( HasStaticE() ) {
   E.resize( get_NStaticNodes() );

   for( Index i = 0 ; i < get_NStaticNodes() ; ++i ) {
    E[ i ].set_both( B.empty() ? 0 : B[ i ] );
    E[ i ].set_function( new LinearFunction( std::move( coeffs[ i ] ) , 0 ) );
    }

   add_static_constraint( E );
   }

  // dynamic part
  if( MayHaveDynE() ) {
   dE.resize( get_NNodes() - get_NStaticNodes() );

   Index i = get_NStaticNodes();
   for( auto & cnst : dE ) {
    cnst.set_both( B.empty() ? 0 : B[ i ] );
    cnst.set_function( new LinearFunction( std::move( coeffs[ i++ ] ) , 0 ) );
    }

   add_dynamic_constraint( dE );
   }

  AR |= HasFlw;
  }

 // generate the bound constraints- - - - - - - - - - - - - - - - - - - - - -

 if( AR & HasBnd )  // bound constraints there already
  return;           // nothing to do

 if( U.empty() ) {
  // if upper bounds are not there and the Configuration says so, the
  // LB0Constraint are not constructed

  auto tstcc = dynamic_cast< SimpleConfiguration< int > * >( stcc );

  if( ( ! tstcc ) && f_BlockConfig &&
      f_BlockConfig->f_static_constraints_Configuration )
   tstcc = dynamic_cast< SimpleConfiguration< int > * >(
                         f_BlockConfig->f_static_constraints_Configuration );
  if( tstcc && ( tstcc->f_value != 0 ) )
   return;
  }

 // static part
 if( HasStaticX() ) {
  UB.resize( get_NStaticArcs() );
  for( Index i = 0 ; i < get_NStaticArcs() ; ++i ) {
   UB[ i ].set_variable( & x[ i ] , eNoBlck );
   UB[ i ].set_rhs( U[ i ] , eNoBlck );
   }

  add_static_constraint( UB );
  }

 // dynamic part
 if( MayHaveDynX() ) {
  dUB.resize( get_NArcs() - get_NStaticArcs() );

  auto dxi = dx.begin();
  auto ui = U.begin() + get_NStaticArcs();
  for( auto & cnst : dUB ) {
   cnst.set_variable( &(*(dxi++)) , eNoBlck );
   cnst.set_rhs( *(ui++) , eNoBlck );
   }

  add_dynamic_constraint( dUB );
  }

 AR |= HasBnd;

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void MCFBlock::generate_objective( Configuration *objc )
{
 if( AR & HasObj )  // the objective is there already
  return;           // cowardly (and silently) return

 // initialize objective function - - - - - - - - - - - - - - - - - - - - - -

 LinearFunction::v_coeff_pair p( get_NArcs() );

 // construct a "dense" LinearFunction- - - - - - - - - - - - - - - - - - - -

 Index i = 0;
 auto Cit = C.begin();

 // static part
 if( HasStaticX() )
  for( ; i < get_NStaticArcs() ; ++i ) {
   p[ i ].first = &x[ i ];
   auto ci = *(Cit++);
   p[ i ].second = std::isnan( ci ) ? 0 : ci;
   }

 // dynamic part
 if( HasDynamicX() ) {
  auto dxi = dx.begin();
  for( ; i < get_NArcs() ; ++i ) {
   p[ i ].first = &(*(dxi++));
   auto ci = *(Cit++);
   p[ i ].second = std::isnan( ci ) ? 0 : ci;
   }
  }

 // ensure no Modification is issued: this may happen in case a MCFBlock
 // is re-loaded, so that set_objective( c ) had already been called
 c.set_function( new LinearFunction( std::move( p ) , 0 ) , eNoMod );

 set_objective( & c , eNoMod );

 AR |= HasObj;

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/

bool MCFBlock::flow_feasible( c_FNumber feps , bool useabstract )
{
 if( useabstract && ( AR & HasFlw ) ) {
  // do it using the abstract representation, if possible - - - - - - - - - -

  if( ! RowConstraint::is_feasible( E , feps ) )   // static part
   return( false );

  if( ! RowConstraint::is_feasible( dE , feps ) )  // dynamic part
   return( false );
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  Index i = 0;
  Vec_FNumber tB = B;

  // static part
  for( ; i < get_NStaticArcs() ; ++i ) {
   c_FNumber xi = x[ i ].get_value();
   tB[ SN[ i ] - 1 ] += xi;
   tB[ EN[ i ] - 1 ] -= xi;
   }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();
   for( ; i < get_NArcs() ; ++i , ++dxi )
    if( ! is_deleted( i ) ) {
     c_FNumber xi = dxi->get_value();
     tB[ SN[ i ] - 1 ] += xi;
     tB[ EN[ i ] - 1 ] -= xi;
     }
   }

  for( Index i = 0 ; i < get_NNodes() ; ++i ) {
   c_FNumber slck = B[ i ] == 0 ? std::abs( tB[ i ] )
                                : std::abs( tB[ i ] / B[ i ] );
   if( slck > feps )
    return( false );
   }
  }

 return( true );

 }  // end( MCFBlock::flow_feasible )

/*--------------------------------------------------------------------------*/

bool MCFBlock::bound_feasible( c_FNumber feps , bool useabstract )
{
 if( useabstract && ( AR & ( HasFlw | HasVar ) ) ) {
  // do it using the abstract representation, if possible - - - - - - - - - -

  // static part
  if( HasStaticX() ) {
   if( UB.empty() ) {
    if( ! ColVariable::is_feasible( x , feps ) )
     return( false );
    }
   else
    if( ! RowConstraint::is_feasible( UB , feps ) )
     return( false );
   }

  // dynamic part
  if( HasDynamicX() ) {
   if( dUB.empty() ) {
    if( ! ColVariable::is_feasible( dx , feps ) )
     return( false );
    }
   else
    if( ! RowConstraint::is_feasible( dUB , feps ) )
     return( false );
   }
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -
  Index i = 0;

  // static part
  for( ; i < get_NStaticArcs() ; ++i ) {
   c_FNumber Ui = get_U( i );
   c_FNumber xi = x[ i ].get_value();
   if( Ui >= Inf< FNumber >() ) {
    if( xi < - feps )
     return( false );
    }
   else {
    c_FNumber slck = Ui == 0 ? std::abs( xi ) :
                               std::max( - xi , xi - Ui ) / std::abs( Ui );
    if( slck > feps )
     return( false );
    }
   }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();
   for(  ; i < get_NArcs() ; ++i ) {
    c_FNumber Ui = get_U( i );
    c_FNumber xi = (*(dxi++)).get_value();
    if( Ui >= Inf< FNumber >() ) {
     if( xi < - feps )
      return( false );
     }
    else {
     c_FNumber slck = Ui == 0 ? std::abs( xi ) :
                                std::max( - xi , xi - Ui ) / std::abs( Ui );
     if( slck > feps )
      return( false );
     }
    }
   }
  }

 return( true );

 }  // end( MCFBlock::bound_feasible )

/*--------------------------------------------------------------------------*/

bool MCFBlock::dual_feasible( c_CNumber ceps , bool useabstract )
{
 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  if( ( ! ( AR & HasFlw ) ) || ( ! ( AR & HasObj ) ) )
   throw( std::logic_error(
	 "abstract representation not there in dual_feasible( , true )" ) );

  auto obj = static_cast< FRealObjective * >( get_objective() );
  assert( obj );
  auto lfo = get_lfo();

  for( auto & pi : (*lfo).get_v_var() ) {
   auto xi = pi.first;
   if( xi->is_fixed() )  // closed or deleted arc
    continue;            // need not be checked
   auto RCi = pi.second;
   for( Index j = 0 ; j < xi->get_num_active() ; ++j ) {
    auto ci = xi->get_active( j );
    if( auto rci = dynamic_cast< FRowConstraint * >( ci ) ) {
     auto lrci = get_lfc( rci );
     RCi -= rci->get_dual() * lrci->get_coefficient( lrci->is_active( xi ) );
     }
    else {
     auto bci = dynamic_cast< BoxConstraint * >( ci );
     assert( bci );
     RCi -= bci->get_dual();
     }
    }

   RCi = std::abs( RCi );
   if( pi.second != 0 )
    RCi /= pi.second;

   if( RCi > ceps )
    return( false );
   }
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  Vec_CNumber RC;
  get_rc( RC.begin() );
  Vec_CNumber Pi;
  get_pi( Pi.begin() );

  for( Index i = 0 ; i < get_NArcs() ; ++i ) {
   if( is_closed( i ) || is_deleted( i ) )
    continue;

   c_CNumber Ci = get_C( i );
   c_CNumber RCi = Ci + Pi[ SN[ i ] - 1 ] - Pi[ EN[ i ] - 1 ];
   c_CNumber df = std::abs( RCi - RC[ i ] );
   c_CNumber mx = std::max( std::abs( Ci ) , df );
   if( mx == 0 ) {
    if( df > ceps )
     return( false );
    }
   else
    if( df > ceps * mx )
     return( false );
   }
  }

 return( true );

 }  // end( MCFBlock::dual_feasible )

/*--------------------------------------------------------------------------*/

bool MCFBlock::complementary_slackness( c_CNumber ceps , c_FNumber feps ,
					bool useabstract )
{
 Vec_CNumber RC;
 get_rc( RC.begin() );

 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  if( ( ! ( AR & HasVar ) ) || ( ! ( AR & HasObj ) ) )
   throw( std::logic_error(
    "abstract representation not there in complementary_slackness(( , true )"
			   ) );

  auto obj = static_cast< FRealObjective * >( get_objective() );
  assert( obj );
  auto lfo = get_lfo();
  Index i = 0;

  // static part
  if( HasStaticX() ) {
   if( UB.empty() ) {
    for( ; i < get_NStaticArcs() ; ++i ) {
     if( x[ i ].is_fixed() )  // arc closed or deleted
      continue;               // need not be checked
     c_CNumber Ci = lfo->get_coefficient( i );
     CNumber RCi = RC[ i ];
     if( Ci )
      RCi /= Ci;
     if( ( x[ i ].get_value() > feps ) && ( RCi < - ceps ) )
      return( false );
     }
    }
   else
    for( ; i < get_NStaticArcs() ; ++i ) {
     if( x[ i ].is_fixed() )  // arc closed or deleted
      continue;               // need not be checked
     c_CNumber Ci = lfo->get_coefficient( i );
     CNumber RCi = RC[ i ];
     if( Ci )
      RCi /= Ci;
     c_FNumber xiv = x[ i ].get_value();
     c_FNumber UBi = UB[ i ].get_rhs();
     if( UBi >= Inf< RowConstraint::RHSValue >() ) {
      if( ( xiv > feps ) && ( RCi < - ceps ) )
       return( false );
      }
     else {
      c_FNumber sfeps = ( UBi == 0 ? feps : feps * UBi );
      if( ( ( xiv > sfeps ) && ( RCi < - ceps ) ) ||
	  ( ( UBi - xiv > sfeps ) && ( RCi > ceps ) ) )
       return( false );
      }
     }
    }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();

   if( dUB.empty() ) {
    for( ; i < get_NStaticArcs() ; ++i , ++dxi)
     if( ! dxi->is_fixed() ) {  // neither closed nor deleted
      c_CNumber Ci = lfo->get_coefficient( i );
      CNumber RCi = RC[ i ];
      if( Ci )
       RCi /= Ci;
      if( ( dxi->get_value() > feps ) && ( RCi < - ceps ) )
       return( false );
      }
     }
   else {
    auto dubi = dUB.begin();

    for( ; i < get_NArcs() ; ++i , ++dxi , ++dubi )
     if( ! dxi->is_fixed() ) {  // neither closed nor deleted
      c_CNumber Ci = lfo->get_coefficient( i );
      CNumber RCi = RC[ i ];
      if( Ci )
       RCi /= Ci;
      c_FNumber dxiv = dxi->get_value();
      c_FNumber UBi = dubi->get_rhs();
      if( UBi >= Inf< RowConstraint::RHSValue >() ) {
       if( ( dxiv > feps ) && ( RCi < - ceps ) )
	return( false );
       }
      else {
       c_FNumber sfeps = ( UBi == 0 ? feps : feps * UBi );
       if( ( ( dxiv > sfeps ) && ( RCi < - ceps ) ) ||
	   ( ( UBi - dxiv > sfeps ) && ( RCi > ceps ) ) )
	return( false );
       }
      }
    }
   }
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -
  Index i = 0;

  // static part
  for(  ; i < get_NStaticArcs() ; ++i )
   if( ( ! std::isnan( C[ i ] ) ) && ( ! x[ i ].is_fixed() ) ) {
    // neither closed nor deleted
    c_CNumber Ci = get_C( i );
    CNumber RCi = RC[ i ];
    if( Ci != 0 )
     RCi /= C[ i ];

    c_FNumber Ui = get_U( i );
    c_FNumber xiv = x[ i ].get_value();

    if( Ui >= Inf< FNumber >() ) {
     if( ( xiv > feps ) && ( RCi < - ceps ) )
      return( false );
     }
    else {
     c_FNumber sfeps = ( Ui == 0 ? feps : feps * Ui );
     if( ( ( xiv > sfeps ) && ( RCi < - ceps ) ) ||
	 ( ( Ui - xiv > sfeps ) && ( RCi > ceps ) ) )
      return( false );
     }
    }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();

   for( ; i < get_NArcs() ; ++i , ++dxi )
    if( ( ! std::isnan( C[ i ] ) ) && ( ! dxi->is_fixed() ) ) {
     // neither closed nor deleted
     c_FNumber dxiv = dxi->get_value();
     c_CNumber Ci = get_C( i );
     CNumber RCi = RC[ i ];
     if( Ci != 0 )
      RCi /= C[ i ];

     c_FNumber Ui = get_U( i );
     if( Ui >= Inf< FNumber >() ) {
      if( ( dxiv > feps ) && ( RCi < - ceps ) )
       return( false );
      }
     else {
      c_FNumber sfeps = ( Ui == 0 ? feps : feps * Ui );
      if( ( ( dxiv > sfeps ) && ( RCi < - ceps ) ) ||
	  ( ( Ui - dxiv > sfeps ) && ( RCi > ceps ) ) )
       return( false );
      }
     }
   }
  }

 return( true );

 }  // end( MCFBlock::complementary_slackness )

/*--------------------------------------------------------------------------*/

bool MCFBlock::is_feasible( bool useabstract , Configuration *fsbc )
{
 FNumber eps = 0;
 auto tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >( fsbc );

 if( ( ! tfsbc ) && f_BlockConfig &&
     f_BlockConfig->f_is_feasible_Configuration )
  tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
                        f_BlockConfig->f_is_feasible_Configuration );
 if( tfsbc )
  eps = tfsbc->f_value;

 return( flow_feasible( eps , useabstract ) &&
	 bound_feasible( eps , useabstract ) );

 }  // end( MCFBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool MCFBlock::is_optimal( bool useabstract , Configuration *optc )
{
 CNumber ceps = 0;
 FNumber feps = 0;
 if( optc ) {
  if( auto toptc =
      dynamic_cast< SimpleConfiguration< std::pair< CNumber , FNumber > > * >(
								    optc ) ) {
   ceps = toptc->f_value.first;
   feps = toptc->f_value.second;
   }
  else {
   auto ttoptc = dynamic_cast< SimpleConfiguration< CNumber > * >( optc );

   if( ( ! ttoptc ) && f_BlockConfig &&
       f_BlockConfig->f_is_optimal_Configuration )
    ttoptc = dynamic_cast< SimpleConfiguration< CNumber > * >(
                           f_BlockConfig->f_is_optimal_Configuration );
   if( ttoptc )
    ceps = ttoptc->f_value;

   if( f_BlockConfig && f_BlockConfig->f_is_feasible_Configuration ) {
    auto fsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
                              f_BlockConfig->f_is_feasible_Configuration );
    if( fsbc )
     feps = fsbc->f_value;
    }
   }
  }
 else
  if( f_BlockConfig ) {
   if( f_BlockConfig->f_is_optimal_Configuration )
    if( auto csbc = dynamic_cast< SimpleConfiguration< CNumber > * >(
                              f_BlockConfig->f_is_optimal_Configuration ) )
     ceps = csbc->f_value;

   if( f_BlockConfig->f_is_feasible_Configuration )
    if( auto fsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
                              f_BlockConfig->f_is_feasible_Configuration ) )
     feps = fsbc->f_value;
   }

 return( flow_feasible( feps , useabstract ) &&
	 bound_feasible( feps , useabstract ) &&
	 dual_feasible( ceps , useabstract ) &&
	 complementary_slackness( ceps , feps , useabstract ) );

 }  //  end( MCFBlock::is_optimal )

/*--------------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/

Block * MCFBlock::get_R3_Block( Configuration *r3bc , Block * base  ,
				Block * father )
{
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 MCFBlock *MCFB;
 if( base ) {
  MCFB = dynamic_cast< MCFBlock * >( base );
  if( ! MCFB )
   throw( std::invalid_argument( "base is not a MCFBlock" ) );
  }
 else
  MCFB = new MCFBlock( father );

 MCFB->load( get_NNodes() , get_NArcs() , EN , SN , U , C , B ,
	     get_NNodes() - get_NStaticNodes() ,
	     get_NArcs() - get_NStaticArcs() ,
	     get_MaxNNodes() - get_NStaticNodes() ,
	     get_MaxNArcs() - get_NStaticArcs() );
 
 return( MCFB );

 }  // end( MCFBlock::get_R3_Block )

/*--------------------------------------------------------------------------*/

void MCFBlock::map_back_solution( Block *R3B , Configuration *r3bc ,
				               Configuration *solc )
{
 // process Configuration - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto MCFB = dynamic_cast< MCFBlock * >( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a MCFBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 int wsol = 0;
 auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc );

 if( ( ! tsolc ) && f_BlockConfig && f_BlockConfig->f_solution_Configuration )
   tsolc = dynamic_cast< SimpleConfiguration< int > * >(
                                    f_BlockConfig->f_solution_Configuration );
 if( tsolc )
  wsol = tsolc->f_value;

 // if required, map back primal solution - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 2 ) && ( AR & HasVar ) ) {  // ... if any
  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static flow size" ) );

  // static part
  if( HasStaticX() )
   for( auto xi = x.begin() , r3bxi = MCFB->x.begin() ; xi != x.end() ;
	++xi , ++r3bxi )
    if( ! xi->is_fixed() )
     xi->set_value( r3bxi->get_value() );
 
  // dynamic part
  // note that if MCFB->dx is longer than this->dx the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( HasDynamicX() ) {
   auto dxi = dx.begin();
   for( auto r3bdxi = MCFB->dx.begin() ;
	( dxi != dx.end() ) && ( r3bdxi != MCFB->dx.end() ) ;
	++dxi , ++r3bdxi )
    if( ! dxi->is_fixed() )
     dxi->set_value( r3bdxi->get_value() );

   for( ; dxi != dx.end() ; ++dxi )
    if( ! dxi->is_fixed() )
     dxi->set_value();
   }
  }

 // if required, map back dual solution - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 1 ) && ( AR & HasFlw ) ) {  // ... if any

  // map back the potentials- - - - - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticNodes() != get_NStaticNodes() )
   throw( std::invalid_argument( "incompatible static potential size" ) );

  // static part
  if( HasStaticE() )
   for( auto ei = E.begin() , r3bei = MCFB->E.begin() ; ei != E.end() ; )
    (ei++)->set_dual( (r3bei++)->get_dual() );
 
  // dynamic part
  // note that if MCFB->dE is longer than this->dE the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( HasDynamicE() ) {
   auto dei = dE.begin();

   for( auto r3bdei = MCFB->dE.begin() ;
	( dei != dE.end() ) && ( r3bdei != MCFB->dE.end() ) ; )
    (dei++)->set_dual( (r3bdei++)->get_dual() );
 
   while( dei != dE.end() )
    (dei++)->set_dual( 0 );
   }

  // map back the reduced costs - - - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static reduced cost size" ) );

  // static part
  if( HasStaticX() && ( ! UB.empty() ) ) {
   if( MCFB->UB.empty() ) {
    for( auto & cnst : UB )
     cnst.set_dual();
    }
   else
    for( auto drci = UB.begin() , r3bdrci = MCFB->UB.begin() ;
	 drci != UB.end() ; )
      (drci++)->set_dual( (r3bdrci++)->get_dual() );
   }

  // dynamic part
  if( HasDynamicX() && ( ! dUB.empty() ) ) {
   if( MCFB->dUB.empty() ) {
    for( auto & cnst : dUB )
     cnst.set_dual();
    }
   else {
    auto drci = dUB.begin();
    for( auto r3bdrci = MCFB->dUB.begin() ;
	 ( drci != dUB.end() ) && ( r3bdrci != MCFB->dUB.end() ) ;
	 ++drci , ++r3bdrci )
     drci->set_dual( r3bdrci->get_dual() );
 
    while( drci != dUB.end() )
     (drci++)->set_dual();
    }
   }
  }
 }  // end( MCFBlock::map_back_solution )

/*--------------------------------------------------------------------------*/

void MCFBlock::map_forward_solution( Block *R3B , Configuration *r3bc ,
				                  Configuration *solc )
{
 // process Configuration - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto MCFB = dynamic_cast< MCFBlock * >( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a MCFBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 int wsol = 0;
 auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc );

 if( ( ! tsolc ) && f_BlockConfig && f_BlockConfig->f_solution_Configuration )
  tsolc = dynamic_cast< SimpleConfiguration< int > * >(
                                    f_BlockConfig->f_solution_Configuration );
 if( tsolc )
  wsol = tsolc->f_value;

 // if required, map forward primal solution- - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 1 ) && ( AR & HasFlw ) ) {  // ... if any
  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static flow size" ) );

  // static part
  if( MCFB->HasStaticX() )
   for( auto xi = x.begin() , r3bxi = MCFB->x.begin() ; xi != x.end() ;
	++xi , ++r3bxi )
    if( ! r3bxi->is_fixed() )
     r3bxi->set_value( xi->get_value() );
 
  // dynamic part
  // note that if this->dx is longer than MCFB->dx the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( MCFB->HasDynamicX() ) {
   auto r3bdxi = MCFB->dx.begin();
   for( auto dxi = dx.begin();
	( dxi != dx.end() ) && ( r3bdxi != MCFB->dx.end() ) ;
	++dxi , ++r3bdxi )
    if( ! r3bdxi->is_fixed() )
     r3bdxi->set_value( dxi->get_value() );

   for( ; r3bdxi != MCFB->dx.end() ; ++r3bdxi )
    if( ! r3bdxi->is_fixed() )
     r3bdxi->set_value();
   }
  }

 // if required, map forward dual solution- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 1 ) && ( AR & HasFlw ) ) {  // ... if any
  // map forward the potentials - - - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticNodes() != get_NStaticNodes() )
   throw( std::invalid_argument( "incompatible static potential size" ) );

  // static part
  if( MCFB->HasStaticE() )
   for( auto ei = E.begin() , r3bei = MCFB->E.begin() ; ei != E.end() ; )
    (r3bei++)->set_dual( (ei++)->get_dual() );
 
  // dynamic part
  // note that if this->dE is longer than MCFB->dE the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( MCFB->HasDynamicE() ) {
   auto r3bdei = MCFB->dE.begin();

   for( auto dei = dE.begin() ;
	( dei != dE.end() ) && ( r3bdei != MCFB->dE.end() ) ; )
    (r3bdei++)->set_dual( (dei++)->get_dual() );
 
   while( r3bdei != MCFB->dE.end() )
    (r3bdei++)->set_dual( 0 );
   }

  // map forward the reduced costs- - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static reduced cost size" ) );

  // static part
  if( MCFB->HasStaticX() && ( ! MCFB->UB.empty() ) ) {
   if( UB.empty() ) {
    for( auto & cnst : MCFB->UB )
     cnst.set_dual( 0 );
    }
   else
    for( auto drci = UB.begin() , r3bdrci = MCFB->UB.begin() ;
	 drci != UB.end() ; )
      (r3bdrci++)->set_dual( (drci++)->get_dual() );
   }

  // dynamic part
  if( MCFB->HasDynamicX() && ( ! MCFB->dUB.empty() ) ) {
   if( UB.empty() ) {
    for( auto & cnst : MCFB->dUB )
     cnst.set_dual( 0 );
    }
   else {
    auto r3bdrci = MCFB->dUB.begin();

    for( auto drci = dUB.begin() ;
	 ( drci != dUB.end() ) && ( r3bdrci != MCFB->dUB.end() ) ; )
     (r3bdrci++)->set_dual( (drci++)->get_dual() );
 
    while( r3bdrci != MCFB->dUB.end() )
     (r3bdrci++)->set_dual( 0 );
    }
   }
  }
 }  // end( MCFBlock::map_forward_solution )

/*--------------------------------------------------------------------------*/

bool MCFBlock::map_forward_Modification( Block * R3B , c_p_Mod mod ,
					 Configuration * r3bc ,
					 ModParam issuePMod ,
					 ModParam issueAMod )
{
 if( mod->concerns_Block() )  // an abstract Modification
  return( false );            // none of my business
 
 auto MCFB = dynamic_cast< MCFBlock * >( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a MCFBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 auto iPM = issuePMod;
 auto iPA = un_ModBlock( issueAMod );

 /* Use a Lambda to define a "guts" of the method that can be called
    recursively without having to pass "local globals". Note the trick of
    defining the std::function object and "passing" it to the lambda,
    which allows recursive calls. Note the need to explicitly capture
    "this" to use fields/methods of the class. */

 std::function< bool( c_p_Mod ) > guts_of_mfM;
 guts_of_mfM = [ this , & guts_of_mfM , & MCFB , & iPM , & iPA ]( c_p_Mod mod
								  ) {
  // process Modification- - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* This requires to patiently sift through the possible Modification types
     to find what this Modification exactly is, and call the appropriate
     method of either MCFB, for a "physical Modification", or of the "abstract
     representation" of MCFB for an "abstract Modification". */

  //!! std::cout << *mod << std::endl;
  
  // GroupModification - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( auto tmod = dynamic_cast< const GroupModification * >( mod ) ) {
   // open or nest the two channels
   iPM = make_par( par2mod( iPM ) , MCFB->open_channel( par2chnl( iPM ) ) );
   iPA = make_par( par2mod( iPA ) , MCFB->open_channel( par2chnl( iPA ) ) );

   bool ok = true;
   for( const auto & submod : tmod->sub_Modifications() )
    if( ! guts_of_mfM( submod.get() ) )
     ok = false;

   // close or un-nest the channels
   MCFB->close_channel( par2chnl( iPM ) );
   MCFB->close_channel( par2chnl( iPA ) );

   return( ok );
   }

  // MCFBlockRngdMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* Note: in the following we can assume that C, B and U are nonempty. This
     is because they can be empty only if they are so when the object is
     loaded. But if a Modification has been issued they are no longer empty
     (a Modification changing nothing from the "empty" state is not issued). */

  if( auto tmod = dynamic_cast< const MCFBlockRngdMod * >( mod ) ) {
   switch( tmod->type() ) {
    case( MCFBlockMod::eChgCost ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      MCFB->chg_cost( C[ tmod->rng().first ] , tmod->rng().first ,
		      iPM , iPA );
     else
      MCFB->chg_costs( C.begin() + tmod->rng().first , tmod->rng() ,
		       iPM , iPA );
     break;
    case( MCFBlockMod::eChgCaps ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      MCFB->chg_ucap( U.empty() ? Inf< FNumber >() : U[ tmod->rng().first ] ,
		      tmod->rng().first , iPM , iPA );
     else
      if( U.empty() ) {
       Vec_FNumber NCap( tmod->rng().second - tmod->rng().first );
       auto NCit = NCap.begin();
       for( Index i = tmod->rng().first ; i < tmod->rng().second ; ++i )
	*(NCit++) = U[ i ];

       MCFB->chg_ucaps( NCap.begin() , tmod->rng() , iPM , iPA );
       }
      else
       MCFB->chg_ucaps( U.begin() + tmod->rng().first , tmod->rng() ,
			iPM , iPA );
     break;
    case( MCFBlockMod::eChgDfct ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NNodes() ) ||
	  ( tmod->rng().second > MCFB->get_NNodes() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      MCFB->chg_dfct( B.empty() ? 0 : B[ tmod->rng().first ] ,
		      tmod->rng().first , iPM , iPA );
     else
      if( B.empty() ) {
       Vec_FNumber NDfct( tmod->rng().second - tmod->rng().first );
       auto NDit = NDfct.begin();
       for( Index i = tmod->rng().first ; i < tmod->rng().second ; ++i )
	*(NDit++) = B[ i ];

       MCFB->chg_ucaps( NDfct.begin() , tmod->rng() , iPM , iPA );
       }
      else
       MCFB->chg_dfcts( B.begin() + tmod->rng().first , tmod->rng() ,
			iPM , iPA );
     break;
    case( MCFBlockMod::eOpenArc ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      MCFB->open_arc( tmod->rng().first , iPM , iPA );
     else
      MCFB->open_arcs( tmod->rng() , iPM , iPA );
     break;
    case( MCFBlockMod::eCloseArc ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      MCFB->close_arc( tmod->rng().first , iPM , iPA );
     else
      MCFB->close_arcs( tmod->rng() , iPM , iPA );
     break;
    case( MCFBlockMod::eAddArc ):
     #ifndef NDEBUG
      if( tmod->rng().first > get_NArcs() )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
     if( MCFB->add_arc( get_SN( tmod->rng().first ) ,
			get_EN( tmod->rng().first ) ,
			get_C( tmod->rng().first ) ,
			get_U( tmod->rng().first ) , iPM , iPA )
	 != tmod->rng().first )
      throw( std::logic_error( "inconsistency between arc names" ) );       
     break;
    case( MCFBlockMod::eRmvArc ):
     #ifndef NDEBUG
      if( tmod->rng().first > MCFB->get_NArcs() )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
     MCFB->remove_arc( tmod->rng().second - 1 , iPM , iPA );
     break;
    default:
     throw( std::invalid_argument( "unknown MCFBlockRngdMod type" ) );
    }
   return( true );
   }

  // MCFBlockSbstMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* Note that tmod->nms() need be copied, since the chg_*() methods
   * *in principle* "consume" the names vector. This is actually not true
   * if MCFB will *not* issue a physical modification, which one may
   * actually know beforehand, but it has to be done anyway because the
   * MCFBlockSbstMod only provides read-only access to the vector.
   * However, tmod->nms() is guaranteed to be ordered. */

  if( auto tmod = dynamic_cast< const MCFBlockSbstMod * >( mod ) ) {
   switch( tmod->type() ) {
    case( MCFBlockMod::eChgCost ): {
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
     Vec_CNumber NCost( tmod->nms().size() );
     for( Index i = 0 ; i < NCost.size() ; i++ )
      NCost[ i ] = C[ tmod->nms()[ i ] ];

     MCFB->chg_costs( NCost.begin() , Subset( tmod->nms() ) , true ,
		      iPM , iPA );
     break;
     }
    case( MCFBlockMod::eChgCaps ): {
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
     Vec_FNumber NCap( tmod->nms().size() , Inf< FNumber >() );
     if( ! U.empty() )
      for( Index i = 0 ; i < NCap.size() ; i++ )
       NCap[ i ] = U[ tmod->nms()[ i ] ];

     MCFB->chg_ucaps( NCap.begin() , Subset( tmod->nms() ) , true ,
		      iPM , iPA );
     break;
     }
    case( MCFBlockMod::eChgDfct ): {
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NNodes() ) ||
	  ( tmod->nms().back() >= MCFB->get_NNodes() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
     Vec_FNumber NDfct( tmod->nms().size() , 0 );
     if( ! B.empty() )
      for( Index i = 0 ; i < NDfct.size() ; i++ )
       NDfct[ i ] = B[ tmod->nms()[ i ] ];

     MCFB->chg_dfcts( NDfct.begin() , Subset( tmod->nms() ) , true ,
		      iPM , iPA );
     break;
     }
    case( MCFBlockMod::eOpenArc ):
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
    MCFB->open_arcs( Subset( tmod->nms() ) , true , iPM , iPA );
     break;
    case( MCFBlockMod::eCloseArc ):
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible MCFBlock" ) );
     #endif
    MCFB->close_arcs( Subset( tmod->nms() ) , true , iPM , iPA );
     break;
    default:
     throw( std::invalid_argument( "unknown MCFBlockSbstMod type" ) );
    }
   return( true );
   }

  // NBModification- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // this is the "nuclear option": the MCFBlock has been re-loaded
  // one should check that the Block is this MCFBlock, but it cannot
  // be otherwise, can it?

  if( auto tmod = dynamic_cast< const NBModification * >( mod ) ) {
   MCFB->load( get_NNodes() , get_NArcs() , EN , SN , U , C , B ,
	       get_NNodes() - get_NStaticNodes() ,
	       get_NArcs() - get_NStaticArcs() ,
	       get_MaxNNodes() - get_NStaticNodes() ,
	       get_MaxNArcs() - get_NStaticArcs() );
   return( true );
   }

  return( false );

  };  // end( guts_of_mfM )- - - - - - - - - - - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // finally, call the "guts of"- - - - - - - - - - - - - - - - - - - - - - - -
 return( guts_of_mfM( mod ) );

 }  // end( MCFBlock::map_forward_Modification )

/*--------------------------------------------------------------------------*/

bool MCFBlock::map_back_Modification( Block *R3B , c_p_Mod mod ,
				      Configuration *r3bc ,
				      ModParam issuePMod ,
				      ModParam issueAMod )
{
 /* Fantastically dirty trick: because the two objects are copies, mapping
  * back a Modification to this from R3B is the same as mapping forward a
  * Modification from R3B to this. */

 auto MCFB = dynamic_cast< MCFBlock * >( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a MCFBlock" ) );

 return( MCFB->map_forward_Modification( this , mod , r3bc , issuePMod ,
					 issueAMod ) );

 }  // end( MCFBlock::map_back_Modification )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * MCFBlock::get_Solution( Configuration * solc , bool emptys )
{
 int wsol = 0;
 if( ( ! solc ) && f_BlockConfig )
  solc = f_BlockConfig->f_solution_Configuration;

 if( auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  wsol = tsolc->f_value;

 auto *sol = new MCFSolution();

 if( wsol != 2 )
  sol->v_x.resize( get_NArcs() );

 if( wsol != 1 )
  sol->v_pi.resize( get_NNodes() );

 if( ! emptys )
  sol->read( this );

 return( sol );

 }  // end( MCFBlock::get_Solution )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_x( Vec_FNumber_it FSol , Range rng ) const
{
 for( ; rng.first < std::min( rng.second , get_NStaticArcs() ) ; )
  *(FSol++) = x[ rng.first++ ].get_value();

 if( HasDynamicX() ) {
  auto dxi = dx.begin();
  for( ; rng.first++ < std::min( rng.second , get_NArcs() ) ; )
   *(FSol++) = (*(dxi++)).get_value();
  }
 }  // end( MCFBlock::get_x( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_x( Vec_FNumber_it FSol , c_Subset & nms ) const
{
 if( ! ( AR & HasVar ) )
  throw( std::logic_error( "flow Variable not available" ) );

 auto nmsi = nms.begin();

 if( HasDynamicX() ) {
  while( ( nmsi != nms.end() ) && ( *nmsi < get_NStaticArcs() ) )
   *(FSol++) = x[ *(nmsi++) ].get_value();

  if( nmsi == nms.end() )
   return;

  auto i = get_NStaticArcs();
  for( auto dxi = dx.begin() ; nmsi != nms.end() ; ++dxi )
   if( *nmsi == i++ ) {
    *(FSol++) = (*dxi).get_value();
    nmsi++;
    }
  }
 else
  while( nmsi != nms.end() ) 
   *(FSol++) = x[ *(nmsi++) ].get_value();

 }  // end( MCFBlock::get_x( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_pi( Vec_CNumber_it PSol , Range rng ) const
{
 if( ! ( AR & HasFlw ) )
  throw( std::logic_error( "potentials unavailable if Constraint aren't" ) );

 Index i = rng.first;
 for( ; i < std::min( rng.second , get_NStaticNodes() ) ; )
  *(PSol++) = E[ i++ ].get_dual();

 if( HasDynamicE() ) {
  auto dei = std::next( dE.begin() , rng.first >= get_NStaticNodes() ?
			             i - get_NStaticNodes() : 0 );
  for( ; i++ < std::min( rng.second , get_NNodes() ) ; )
   *(PSol++) = (*(dei++)).get_dual();
  }
 }  // end( MCFBlock::get_pi( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_pi( Vec_CNumber_it PSol , c_Subset & nms ) const
{
 if( ! ( AR & HasFlw ) )
  throw( std::logic_error( "potentials unavailable if Constraint aren't" ) );

 auto nmsi = nms.begin();

 if( HasDynamicE() ) {
  while( ( nmsi != nms.end() ) && ( *nmsi < get_NStaticNodes() ) )
   *(PSol++) = E[ *(nmsi++) ].get_dual();

  if( nmsi == nms.end() )
   return;

  auto i = get_NStaticNodes();
  for( auto dei = dE.begin() ; nmsi != nms.end() ; ++dei )
   if( *nmsi == i++ ) {
    *(PSol++) = (*dei).get_dual();
    nmsi++;
    }
  }
 else
  while( nmsi != nms.end() )
   *(PSol++) = E[ *(nmsi++) ].get_dual();

 }  // end( MCFBlock::get_pi( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_rc( Vec_CNumber_it RC , Range rng ) const
{
 if( ! ( AR & HasFlw ) )
  throw( std::logic_error( "reduced costs unavailable if Constraint aren't" )
	 );

 if( AR & HasBnd ) {
  Index i = rng.first;
  if( HasStaticX() )
   for( ; i < std::min( rng.second , get_NStaticArcs() ) ; )
    *(RC++) = UB[ i++ ].get_dual();

  if( HasDynamicX() ) {
   auto dubi = std::next( dUB.begin() , rng.first >= get_NStaticArcs() ?
			                i - get_NStaticArcs() : 0 );
   for( ; i++ < std::min( rng.second , get_NArcs() ) ; )
    *(RC++) = (*(dubi++)).get_dual();
   }
  }
 else
  for( ; rng.first < std::min( rng.second , get_NArcs() ) ; ++rng.first )
   *(RC++) = C[ rng.first ] + get_pi( SN[ rng.first ] - 1 )
                            - get_pi( EN[ rng.first ] - 1 );

 }  // end( MCFBlock::get_rc( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_rc( Vec_CNumber_it RC , c_Subset & nms ) const
{
 if( ! ( AR & HasFlw ) )
  throw( std::logic_error( "reduced costs unavailable if Constraint aren't" )
	 );

 auto nmsi = nms.begin();

 if( AR & HasBnd ) {
  if( HasDynamicX() ) {
   while( ( nmsi != nms.end() ) && ( *nmsi < get_NStaticArcs() ) )
     *(RC++) = UB[ *(nmsi++) ].get_dual();

   if( nmsi == nms.end() )
    return;

   auto i = get_NStaticArcs();
   for( auto dubi = dUB.begin() ; nmsi != nms.end() ; ++dubi )
    if( *nmsi == i++ ) {
     *(RC++) = (*dubi).get_dual();
     nmsi++;
     }
   }
  else
   while( nmsi != nms.end() ) 
    *(RC++) = UB[ *(nmsi++) ].get_dual();
  }
 else
  for( ; nmsi != nms.end() ; ++nmsi ) 
   *(RC++) = C[ *nmsi ] + get_pi( SN[ *nmsi ] - 1 )
                        - get_pi( EN[ *nmsi ] - 1 );
 
 }  // end( MCFBlock::get_rc( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::set_x( c_Vec_FNumber_it fstrt , Range rng )
{
 if( ! ( AR & HasVar ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 if( rng.second > get_NArcs() )
  rng.second = get_NArcs();

 Index i = rng.first;

 if( HasStaticX() )
  for( auto xi = x.begin() + i ;
       i < std::min( rng.second , get_NStaticArcs() ) ; ++i )
   (xi++)->set_value( *(fstrt++) );

 if( ( ! HasDynamicX() ) || ( rng.second <= get_NStaticArcs() ) )
  return;

 auto dxi = dx.begin();
 if( i > get_NStaticArcs() )
  dxi = std::next( dxi , i - get_NStaticArcs() );

 for( ; i < std::min( rng.second , get_NArcs() ) ; ++i )
  (dxi++)->set_value( *(fstrt++) );

 }  // end( MCFBlock::set_x( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::set_x( c_Vec_FNumber_it fstrt , c_Subset sbst )
{
 if( ! ( AR & HasVar ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 auto it = sbst.begin();

 if( HasStaticX() )
  while( ( it != sbst.end() ) && ( *it < get_NStaticArcs() ) )
   x[ *(it++) ].set_value( *(fstrt++) );

 if( ! HasDynamicX() )
  return;

 auto dxi = dx.begin();
 Index prev = get_NStaticArcs();

 while( it != sbst.end() ) {
  Index h = *(it++);
  dxi = std::next( dxi , h - prev );
  dxi->set_value( *(fstrt++) );
  prev = h;
  }
 }  // end( MCFBlock::set_x( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::set_pi( c_Vec_CNumber_it pstrt , Range rng )
{
 if( ! ( AR & HasFlw ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 if( rng.second > get_NNodes() )
  rng.second = get_NNodes();

 Index i = rng.first;

 if( HasStaticE() )
  for( auto ei = E.begin() + i ;
       i < std::min( rng.second , get_NStaticNodes() ) ; ++i )
   (ei++)->set_dual( *(pstrt++) );

 if( ( ! HasDynamicE() ) || ( rng.second <= get_NStaticNodes() ) )
  return;

 auto dei = dE.begin();
 if( i > get_NStaticNodes() )
  dei = std::next( dei , i - get_NStaticNodes() );

 for( ; i < std::min( rng.second , get_NNodes() ) ; ++i )
  (dei++)->set_dual( *(pstrt++) );

 }  // end( MCFBlock::set_pi( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::set_pi( c_Vec_CNumber_it pstrt , c_Subset sbst )
{
 if( ! ( AR & HasFlw ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 auto it = sbst.begin();

 if( HasStaticE() )
  while( ( it != sbst.end() ) && ( *it < get_NStaticNodes() ) )
   E[ *(it++) ].set_dual( *(pstrt++) );

 if( ! HasDynamicE() )
  return;

 auto dEi = dE.begin();
 Index prev = get_NStaticNodes();

 while( it != sbst.end() ) {
  Index h = *(it++);
  dEi = std::next( dEi , h - prev );
  dEi->set_dual( *(pstrt++) );
  prev = h;
  }
 }  // end( MCFBlock::set_pi( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::set_rc( c_Vec_CNumber_it rcstrt , Range rng )
{
 if( ! ( AR & HasBnd ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 if( rng.second > get_NArcs() )
  rng.second = get_NArcs();

 Index i = rng.first;

 if( HasStaticX() )
  for( auto ubi = UB.begin() + i ;
       i < std::min( rng.second , get_NStaticArcs() ) ; ++i )
   (ubi++)->set_dual( *(rcstrt++) );

 if( ( ! HasDynamicX() ) || ( rng.second <= get_NStaticArcs() ) )
  return;

 auto dubi = dUB.begin();
 if( i > get_NStaticArcs() )
  dubi = std::next( dubi , i - get_NStaticArcs() );

 for( ; i < std::min( rng.second , get_NArcs() ) ; ++i )
  (dubi++)->set_dual( *(rcstrt++) );

 }  // end( MCFBlock::set_rc( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::set_rc( c_Vec_CNumber_it rcstrt , c_Subset sbst )
{
 if( ! ( AR & HasBnd ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 auto it = sbst.begin();

 if( HasStaticX() )
  while( ( it != sbst.end() ) && ( *it < get_NStaticArcs() ) )
   UB[ *(it++) ].set_dual( *(rcstrt++) );

 if( ! HasDynamicX() )
  return;

 auto dubi = dUB.begin();
 Index prev = get_NStaticArcs();

 while( it != sbst.end() ) {
  Index h = *(it++);
  dubi = std::next( dubi , h - prev );
  dubi->set_dual( *(rcstrt++) );
  prev = h;
  }
 }  // end( MCFBlock::set_rc( subset ) )

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void MCFBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 //!! std::cout << *mod << std::endl;

 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );
 }

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE MCFBlock ---------*/
/*--------------------------------------------------------------------------*/

void MCFBlock::print( std::ostream  & output , char vlvl ) const
{
 if( vlvl != 'C' ) {  // non-complete version
  // only basic information 
  output << "MCFBlock with: " << NNodes << " nodes and " << SN.size()
	 << " arcs" << std::endl;

  if( ! vlvl ) {     // print the graph
   if( ! B.empty() )
    for( Index i = 0 ; i < get_NNodes() ; ++i )
     if( B[ i ] != 0 )
      output << "B[ " << i + 1 << " ] = " << B[ i ] << std::endl;

   for( Index i = 0 ; i < get_NArcs() ; ++i ) {
    output << "( " << SN[ i ] << " , " << EN[ i ] << " ): C = ";
    if( is_deleted( i ) )
     output << "1000000";
    else
     output << C[ i ];

    if( ! U.empty() ) {
     output << ", U = ";
     print_UB( output , U[ i ] );
     }

    output << std::endl;
    }
   }
  }
 else  {
  // print header in DIMACS standard format
  output << "p min " << get_NNodes() << " " << get_NArcs() << std::endl;
  output << std::setprecision( 16 );

  // print node descriptors in DIMACS standard format
  if( ! B.empty() )
   for( Index i = 0 ; i < get_NNodes() ; ++i )
    if( B[ i ] != 0 )
     output << "n\t" << i + 1 << "\t" << - B[ i ] << std::endl;

  // print arc descriptors in DIMACS standard format
  for( Index i = 0 ; i < get_NArcs() ; ++i ) {
   output << "a\t" << SN[ i ] << "\t" << EN[ i ] << "\t0\t";

   if( ( is_deleted( i ) ) || is_closed( i ) )
    output << "0";
   else {
    if( U.empty() )
     output << "+Inf";
    else
     print_UB( output , U[ i ] );
    }

   output << "\t";
   if( is_deleted( i ) )
    output << "1000000";
   else
    output << C[ i ];
   output << std::endl;
   }
  }
 }  // end( MCFBlock::print )

/*--------------------------------------------------------------------------*/

void MCFBlock::serialize( netCDF::NcGroup & group ) const
{
 // call the method of Block- - - - - - - - - - - - - - - - - - - - - - - - -

 Block::serialize( group );

 // now the MCFBlock data - - - - - - - - - - - - - - - - - - - - - - - - - -

 netCDF::NcDim nn = group.addDim( "NNodes" , get_NNodes() );
 netCDF::NcDim na = group.addDim( "NArcs" , get_NArcs() );

 if( get_NNodes() > get_NStaticNodes() )
  group.addDim( "DynNNodes" , get_NNodes() - get_NStaticNodes() );

 if( get_NArcs() > get_NStaticArcs() )
  group.addDim( "DynNArcs" , get_NArcs() - get_NStaticArcs() );

 if( get_MaxNNodes() > get_NStaticNodes() )
  group.addDim( "MaxDynNNodes" , get_MaxNNodes() - get_NStaticNodes() );

 if( get_MaxNArcs() > get_NStaticArcs() )
  group.addDim( "MaxDynNArcs" , get_MaxNArcs() - get_NStaticArcs() );

 ( group.addVar( "SN" , netCDF::NcUint64() , na ) ).putVar( SN.data() );

 ( group.addVar( "EN" , netCDF::NcUint64() , na ) ).putVar( EN.data() );

 ( group.addVar( "C" , netCDF::NcDouble() , na ) ).putVar( C.data() );

 if( ! U.empty() )
  ( group.addVar( "U" , netCDF::NcDouble() , na ) ).putVar( U.data() );

 if( ! B.empty() )
  ( group.addVar( "B" , netCDF::NcDouble() , nn ) ).putVar( B.data() );

 }  // end( MCFBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void MCFBlock::chg_costs( c_Vec_CNumber_it NCost , Range rng ,
			  ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 // check to see how many of the initial arcs are either deleted or not
 // really changing the costs
 while( ( std::isnan( C[ rng.first ] ) || ( *NCost == C[ rng.first ] ) )
	&& ( rng.first < rng.second ) ) {
  ++rng.first;
  ++NCost;
  }

 if( rng.second <= rng.first )  // nothing left to change
  return;                       // cowardly (and silently) return

 // check to see how many of the final arcs are either deleted or not
 // really changing the costs
 auto NCEit = NCost + ( rng.second - rng.first );
 while( ( ( std::isnan( C[ rng.second - 1 ] ) ) ||
	  ( *(--NCEit) == C[ rng.second - 1 ] ) )
	&& ( rng.first < rng.second ) )
  --rng.second;

 if( rng.second <= rng.first )  // nothing left to change
  return;                       // cowardly (and silently) return

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  Vec_CNumber NC( rng.second - rng.first );
  auto NCit = NC.begin();

  for( Index i = rng.first ; i < rng.second ; ++i , ++NCost )
   if( std::isnan( C[ i ] ) )  // arc is deleted
    *(NCit++) = 0;             // give it an "harmless" coefficient
   else                        // arc is there    
    *(NCit++) = C[ i ] = *NCost;

  get_lfo()->modify_coefficients( std::move( NC ) , rng ,
				  un_ModBlock( issueAMod ) );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   for( Index i = rng.first ; i < rng.second ; ++i , ++NCost )
    if( ! std::isnan( C[ i ] ) )
     C[ i ] = *NCost;

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
					      MCFBlockMod::eChgCost , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::chg_costs( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_costs( c_Vec_CNumber_it NCost , Subset && nms ,
			  bool ordered  ,
			  ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 // eliminate from NCost and nms the entries corresponding to either
 // deleted arcs or arcs whose cost actually does not change; meanwhile,
 // if nms is not ordered, order it
 Vec_CNumber NC;
 if( ordered ) {
  NC.resize( nms.size() );
  auto NCit = NC.begin();
  auto nmsit = nms.begin();
  for( auto i : nms ) {
   auto nci = *(NCost++);
   if( ( ! std::isnan( C[ i ] ) ) && ( nci != C[ i ] ) ) {
    *(nmsit++) = i;
    *(NCit++) = nci;
    }
   }
  nms.resize( std::distance( nms.begin() , nmsit ) );
  NC.resize( nms.size() );
  }
 else {
  using TP = std::pair< Index , CNumber >;
  std::vector< TP > pairs;
  pairs.reserve( nms.size() );
  for( auto i : nms ) {
   auto nci = *(NCost++);
   if( ( ! std::isnan( C[ i ] ) ) && ( nci != C[ i ] ) )
    pairs.push_back( std::make_pair( i , nci ) );
   }
  std::sort( pairs.begin() , pairs.end() ,
	     []( auto & a , auto & b ) { return( a.first < b.first ); } );
  NC.resize( pairs.size() );
  auto NCit = NC.begin();
  auto nmsit = nms.begin();
  for( auto & el : pairs ) {
   *(nmsit++) = el.first;
   *(NCit++) = el.second;
   }
  nms.resize( pairs.size() );
  }

 if( nms.empty() )  // nothing left to change
  return;           // cowardly (and silently) return

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  copyidx( C , nms , NC.begin() );

  // note that modify_coefficients owns both vectors, so two copies have
  // to be made
  get_lfo()->modify_coefficients( std::move( NC ) , Subset( nms ) , true ,
				  un_ModBlock( issueAMod ) );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( C , nms , NC.begin() );

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockSbstMod >( this ,
                                  MCFBlockMod::eChgCost , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::chg_costs( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_cost( CNumber NCost , Index arc , 
			 ModParam issueMod , ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( std::isnan( C[ arc ] ) || ( C[ arc ] == NCost ) )
  return;  // if the arc is deleted or the cost does not change, all done

 f_cond_lower = dNAN;  // reset conditional bounds

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  C[ arc ] = NCost;

  get_lfo()->modify_coefficient( arc , NCost , un_ModBlock( issueAMod ) );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   C[ arc ] = NCost;

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
			   MCFBlockMod::eChgCost , Range( arc , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::chg_cost )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_ucaps( c_Vec_FNumber_it NCap , Range rng ,
			  ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 if( U.empty() ) {
  if( std::all_of( NCap , NCap + ( rng.second - rng.first ) ,
		   []( c_FNumber cap ) { return( cap >= Inf< FNumber >() ); }
		   ) )
   return;

  U.assign( get_MaxNArcs() , Inf< FNumber >() );
  }

 // check to see how many of the initial arcs are either deleted or not
 // really changing the capacity
 while( ( std::isnan( C[ rng.first ] ) || ( *NCap == U[ rng.first ] ) )
	&& ( rng.first < rng.second ) ) {
  ++rng.first;
  ++NCap;
  }

 if( rng.second <= rng.first )  // nothing left to change
  return;                       // cowardly (and silently) return

 // check to see how many of the final arcs are either deleted or not
 // really changing the capacity
 auto NCEit = NCap + ( rng.second - rng.first );
 while( ( ( std::isnan( C[ rng.second - 1 ] ) ) ||
	  ( *(--NCEit) == U[ rng.second - 1 ] ) )
	&& ( rng.first < rng.second ) )
  --rng.second;

 if( rng.second <= rng.first )  // nothing left to change
  return;                       // cowardly (and silently) return

 f_cond_lower = dNAN;  // reset conditional bounds

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( ! ( AR & HasBnd ) )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , rng.second - rng.first );

  Index i = rng.first;

  // static part
  for( ; i < std::min( rng.second , get_NStaticArcs() ) ;  ++i , ++NCap )
   if( ( ! std::isnan( C[ i ] ) ) && ( U[ i ] != *NCap ) ) {
    U[ i ] = *NCap;
    UB[ i ].set_rhs( *NCap , ampar );
    }

  // dynamic part
  for( auto dubi = std::next( dUB.begin() ,
			      rng.first >= get_NStaticArcs() ?
			      i - get_NStaticArcs() : 0 ) ;
       i < rng.second ; ++i , ++NCap , ++dubi )
   if( ( ! std::isnan( C[ i ] ) ) && ( U[ i ] != *NCap ) ) {
    U[ i ] = *NCap;
    dubi->set_rhs( *NCap , ampar );
    }

  close_if_needed( ampar , rng.second - rng.first );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  // note that this also changes the capacity of deleted arcs, but since that
  // is never really used, it does not matter
  if( not_dry_run( issueMod ) )
   std::copy( NCap , NCap + ( rng.second - rng.first ) ,
	      U.begin() + rng.first );

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
				              MCFBlockMod::eChgCaps , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::chg_ucaps( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_ucaps( c_Vec_FNumber_it NCap , Subset && nms ,
			  bool ordered  ,
			  ModParam issueMod , ModParam issueAMod )
{
 if( U.empty() ) {
  if( std::all_of( NCap , NCap + nms.size() ,
		   []( c_FNumber cap ) { return( cap >= Inf< FNumber >() ); }
		   ) )
   return;

  U.assign( get_MaxNArcs() , Inf< FNumber >() );
  }

 // eliminate from NCap and nms the entries corresponding to either
 // deleted arcs or arcs whose capacity actually does not change;
 // meanwhile, if nms is not ordered, order it
 Vec_FNumber NC;
 if( ordered ) {
  NC.resize( nms.size() );
  auto NCit = NC.begin();
  auto nmsit = nms.begin();
  for( auto i : nms ) {
   auto nci = *(NCap++);
   if( ( ! std::isnan( C[ i ] ) ) && ( nci != U[ i ] ) ) {
    *(nmsit++) = i;
    *(NCit++) = nci;
    }
   }
  nms.resize( std::distance( nms.begin() , nmsit ) );
  }
 else {
  using TP = std::pair< Index , FNumber >;
  std::vector< TP > pairs;
  pairs.reserve( nms.size() );
  for( auto i : nms ) {
   auto nci = *(NCap++);
   if( ( ! std::isnan( C[ i ] ) ) && ( nci != U[ i ] ) )
    pairs.push_back( std::make_pair( i , nci ) );
   }
  std::sort( pairs.begin() , pairs.end() ,
	     []( auto & a , auto & b ) { return( a.first < b.first ); } );
  NC.resize( pairs.size() );
  auto NCit = NC.begin();
  auto nmsit = nms.begin();
  for( auto & el : pairs ) {
   *(nmsit++) = el.first;
   *(NCit++) = el.second;
   }
  nms.resize( pairs.size() );
  }

 if( nms.empty() )  // nothing left to change
  return;           // cowardly (and silently) return

 f_cond_lower = dNAN;  // reset conditional bounds

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( ! ( AR & HasBnd ) )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , nms.size() );

  // static part
  auto nit = nms.begin();
  auto NCit = NC.begin();
  for( ; ( nit != nms.end() ) && ( *nit < get_NStaticArcs() ) ;
       ++NCit , ++nit ) {
   if( ( ! std::isnan( C[ *nit ] ) ) && ( U[ *nit ] != *NCit ) ) {
    U[ *nit ] = *NCit;
    UB[ *nit ].set_rhs( *NCit , ampar );
    }
   }

  if( HasDynamicX() ) {
   auto dubi = dUB.begin();
   for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dubi )
    if( *nit == i ) {
     if( ( ! std::isnan( C[ i ] ) ) && ( U[ i ] != *NCit ) ) {
      U[ i ] = *NCit;
      dubi->set_rhs( *NCit , ampar );
      }
     ++nit;
     ++NCit;
     }
   }

  close_if_needed( ampar , nms.size() );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  // note that this also changes the capacity of deleted arcs, but since that
  // is never really used, it does not matter
  if( not_dry_run( issueMod ) )
   copyidx( U , nms , NC.begin() );

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockSbstMod >( this ,
                                  MCFBlockMod::eChgCaps , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::chg_ucaps( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_ucap( FNumber NCap , Index arc ,
			 ModParam issueMod , ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( U.empty() ) {
  if( NCap < Inf< FNumber >() )
   return;

  U.assign( get_MaxNArcs() , Inf< FNumber >() );
  }

 if( U[ arc ] == NCap )
  return;

 f_cond_lower = dNAN;  // reset conditional bounds

 if( not_dry_run( issueMod ) )
  U[ arc ] = NCap;  // only change the physical representation - - - - - - -

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change the abstract representation - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( ! ( AR & HasBnd ) )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  if( arc < get_NStaticArcs() )
   UB[ arc ].set_rhs( NCap , issueAMod );
  else
   std::next( dUB.begin() , arc - get_NStaticArcs()
	      )->set_rhs( NCap , un_ModBlock( issueAMod ) );
  }

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
			   MCFBlockMod::eChgCaps , Range( arc , arc + 1 )  ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::chg_ucap )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_dfcts( c_Vec_CNumber_it NDfct , Range rng ,
			  ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NNodes() );
 if( rng.second <= rng.first )  // nothing to change
  return;                 // cowardly (and silently) return

 if( B.empty() ) {
  if( std::all_of( NDfct , NDfct + ( rng.second - rng.first ) ,
		   []( c_FNumber dfct ) { return( dfct == 0 ); } ) )
   return;

  B.assign( get_MaxNNodes() , 0 );
  }

 c_Index ndiff = countdiff( NDfct , NDfct + ( rng.second - rng.first ) ,
			    B.cbegin() + rng.first );
 if( ! ndiff )
  return;

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , ndiff );

  Index i = rng.first;

  // static part
  for( ; i < std::min( rng.second , get_NStaticNodes() ) ;  ++i , ++NDfct )
   if( B[ i ] != *NDfct ) {
    B[ i ] = *NDfct;
    E[ i ].set_both( *NDfct , ampar );
    }

  // dynamic part
  for( auto dei = std::next( dE.begin() ,
			     rng.first >= get_NStaticNodes() ?
			     i - get_NStaticNodes() : 0 ) ;
       i < rng.second ; ++i , ++NDfct , ++dei )
   if( B[ i ] != *NDfct ) {
    B[ i ] = *NDfct;
    dei->set_both( *NDfct , ampar );
    }

  close_if_needed( ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NDfct , NDfct + ( rng.second - rng.first ) ,
	      B.begin() + rng.first );

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
				              MCFBlockMod::eChgDfct , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::chg_dfcts( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_dfcts( c_Vec_CNumber_it NDfct , Subset && nms ,
			  bool ordered ,
			  ModParam issueMod , ModParam issueAMod )
{
 if( B.empty() ) {
  if( std::all_of( NDfct , NDfct + nms.size() ,
		   []( c_FNumber dfct ) { return( dfct == 0 ); } ) )
   return;

  B.assign( get_MaxNNodes() , 0 );
  }

 Index ndiff = countdiff( B , nms , NDfct , get_NNodes() );
 if( ! ndiff )
  return;

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , ndiff );

  if( HasDynamicE() )
   if( ordered ) {
    // static part
    auto nit = nms.begin();
    for( ; ( nit != nms.end() ) && ( *nit < get_NStaticNodes() ) ;
	 ++NDfct , ++nit ) {
     if( B[ *nit ] != *NDfct ) {
      B[ *nit ] = *NDfct;
      E[ *nit ].set_both( *NDfct , ampar );
      }
     }

    // dynamic part
    auto dei = dE.begin();
    for( Index i = get_NStaticNodes() ; nit != nms.end() ; ++i , ++dei )
     if( *nit == i ) {
      if( B[ i ] != *NDfct ) {
       B[ i ] = *NDfct;
       dei->set_both( *NDfct , ampar );
       }
      nit++;
      NDfct++;
      }
    }
   else {
    // make a vector of pairs < arc index , new capacity >
    typedef std::pair< Index , FNumber > index_pair;
    std::vector< index_pair > pairs( nms.size() );
    for( Index i = 0 ; i < nms.size() ; ++i )
     pairs[ i ] = std::make_pair( nms[ i ] , *(NDfct++) );

    // sort the vector for increasing index
    std::sort( pairs.begin() , pairs.end() ,
	       []( index_pair i , index_pair j )
	       { return( i.first < j.first ); } );

    // static part
    auto pit = pairs.begin();
    for( ; ( pit != pairs.end() ) && ( pit->first < get_NStaticArcs() ) ;
	 ++pit )
     if( B[ pit->first ] != pit->second ) {
      B[ pit->first ] = pit->second;
      E[ pit->first ].set_both( pit->second , ampar );
      }

    // dynamic part
    auto dei = dE.begin();
    for( Index i = get_NStaticNodes() ; pit != pairs.end() ; ++i , ++dei )
     if( pit->first == i ) {
      if( B[ i ] != pit->second ) {
       B[ i ] = pit->second;
       dei->set_both( pit->second , ampar );
       }
      pit++;
      }
    }
  else
   for( auto nit = nms.begin() ; nit != nms.end() ; ++NDfct , ++nit ) {
    if( B[ *nit ] != *NDfct ) {
     B[ *nit ] = *NDfct;
     E[ *nit ].set_both( *NDfct , ampar );
     }
    }

  close_if_needed( ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( B , nms , NDfct );

 // TODO: eliminate from nms the "fake" changes

 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  // ensure the names are ordered even if they were not so originally
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  Block::add_Modification( std::make_shared< MCFBlockSbstMod >( this ,
                                 MCFBlockMod::eChgDfct , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
  }

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::chg_dfcts( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_dfct( FNumber NDfct , Index nde ,
			 ModParam issueMod , ModParam issueAMod )
{
 if( nde >= get_NNodes() )
  throw( std::invalid_argument( "invalid node name" ) );

 if( B.empty() ) {
  if( ! NDfct )
   return;

  B.assign( get_NNodes() , 0 );
  }

 if( B[ nde ] == NDfct )
  return;

 if( not_dry_run( issueMod ) )
  B[ nde ] = NDfct;  // change the physical representation- - - - - - - - - -

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change the abstract representation - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  if( nde < get_NStaticNodes() )
   E[ nde ].set_both( NDfct , issueAMod );
  else
   std::next( dE.begin() , nde - get_NStaticNodes()
	      )->set_both( NDfct , un_ModBlock( issueAMod ) );
  }
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
			   MCFBlockMod::eChgDfct , Range( nde , nde + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::chg_dfct )

/*--------------------------------------------------------------------------*/

void MCFBlock::close_arcs( Range rng ,
			   ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  std::vector< ColVariable * > toclose;
  toclose.reserve( rng.second - rng.first );
  Index i = rng.first;

  // static part
  for( ; i < std::min( rng.second , get_NStaticArcs() ) ; ++i )
   if( ! x[ i ].is_fixed() )
    toclose.push_back( & x[ i ] );

  // dynamic part
  if( ( rng.second > get_NStaticArcs() ) && HasDynamicX() )
   for( auto dxi = std::next( dx.begin() , i - get_NStaticArcs() ) ;
	i++ < rng.second ; ++dxi )
    if( ! dxi->is_fixed() )
     toclose.push_back( & (*dxi) );

  if( toclose.empty() )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , toclose.size() );

  for( auto vi : toclose ) {
   vi->set_value( 0 );
   vi->is_fixed( true , ampar );
   }

  close_if_needed( ampar , toclose.size() );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
				             MCFBlockMod::eCloseArc , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::close_arcs( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::close_arcs( Subset && nms , bool ordered  ,
			   ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )
  return;

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  std::vector< ColVariable * > toclose;
  toclose.reserve( nms.size() );

  // static part
  auto nit = nms.begin();
  for( ; ( nit != nms.end() ) && ( *nit < get_NStaticArcs() ) ; ++nit )
   if( ! x[ *nit ].is_fixed() )
    toclose.push_back( & x[ *nit ] );

  // dynamic part
  auto dxi = dx.begin();
  for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dxi )
   if( *nit == i ) {
    if( ! dxi->is_fixed() )
     toclose.push_back( & (*dxi) );
    ++nit;
    }

  if( toclose.empty() )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , toclose.size() );

  for( auto vi : toclose ) {
   vi->set_value( 0 );
   vi->is_fixed( true , ampar );
   }

  close_if_needed( ampar , toclose.size() );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockSbstMod >( this ,
                                 MCFBlockMod::eCloseArc , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::close_arcs( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::close_arc( Index arc ,
			  ModParam issueMod , ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  auto xa = i2p_x( arc );

  if( xa->is_fixed() )
   return;

  xa->set_value( 0 );

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  xa->is_fixed( true , un_ModBlock( issueAMod ) );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
			   MCFBlockMod::eCloseArc , Range( arc , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::close_arc )

/*--------------------------------------------------------------------------*/

void MCFBlock::open_arcs( Range rng ,
			  ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  std::vector< ColVariable * > toopen;
  toopen.reserve( rng.second - rng.first );
  Index i = rng.first;

   // static part
  for( ; i < std::min( rng.second , get_NStaticArcs() ) ; ++i )
   if( x[ i ].is_fixed() && ( ! std::isnan( C[ i ] ) ) )
    toopen.push_back( & x[ i ] );

  // dynamic part
  if( ( rng.second > get_NStaticArcs() ) && HasDynamicX() )
   for( auto dxi = std::next( dx.begin() , i - get_NStaticArcs() ) ;
	i < rng.second ; ++i , ++dxi )
    if( dxi->is_fixed() && ( ! std::isnan( C[ i ] ) ) )
     toopen.push_back( & (*dxi) );

  if( toopen.empty() )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , toopen.size() );

  for( auto vi : toopen )
   vi->is_fixed( false , ampar );

  close_if_needed( ampar , toopen.size() );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
				             MCFBlockMod::eOpenArc , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::open_arcs( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::open_arcs( Subset && nms , bool ordered  ,
			  ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )
  return;

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  std::vector< ColVariable * > toopen;
  toopen.reserve( nms.size() );

  // static part
  auto nit = nms.begin();
  for( ; ( nit != nms.end() ) && ( *nit < get_NStaticArcs() ) ; ++nit )
   if( x[ *nit ].is_fixed() && ( ! std::isnan( C[ *nit ] ) ) )
    toopen.push_back( & x[ *nit ] );

  // dynamic part
  auto dxi = dx.begin();
  for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dxi )
   if( *nit == i ) {
    if( dxi->is_fixed() && ( ! std::isnan( C[ i ] ) ) )
     toopen.push_back( & (*dxi) );
    ++nit;
    }

  if( toopen.empty() )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , toopen.size() );

  for( auto vi : toopen )
   vi->is_fixed( false , ampar );

  close_if_needed( ampar , toopen.size() );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockSbstMod >( this ,
                                  MCFBlockMod::eOpenArc , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::open_arcs( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::open_arc( Index arc ,
			 ModParam issueMod , ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  auto xa = i2p_x( arc );

  if( ( ! xa->is_fixed() ) || std::isnan( C[ arc ] ) )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  xa->is_fixed( false , un_ModBlock( issueAMod ) );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
			   MCFBlockMod::eOpenArc , Range( arc , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::open_arc )

/*--------------------------------------------------------------------------*/

MCFBlock::Index MCFBlock::add_arc( Index sn , Index en ,
				   CNumber cst , FNumber cap ,
				   ModParam issueMod , ModParam issueAMod )
{
 if( ( sn < 1 ) || ( sn > get_NNodes() ) )
  throw( std::invalid_argument(
			 "MCFBlock::add_arc: invalid starting node name" ) );
 
 if( ( en < 1 ) || ( en > get_NNodes() ) )
  throw( std::invalid_argument(
			   "MCFBlock::add_arc: invalid ending node name" ) );

 Index arc = get_NStaticArcs();
 while( ( arc < get_NArcs() ) && ( ! is_deleted( arc ) ) )
  ++arc;

 if( arc >= get_MaxNArcs() )
  return( Inf< Index >() );

 // change the abstract representation- - - - - - - - - - - - - - - - - - - -
 // in the meantime, if so instructed also issue abstract Modification(s)
 // note that this is *always* done, unless issueAMod says this is a dry
 // run, because at least the BlockModAdd corresponding to adding the
 // Variable, or unfixing it, is always issued since the Variable are
 // always present

 if( not_dry_run( issueAMod ) ) {
  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod ,
			       AR & ( HasFlw | HasObj ) ? 4 : 1 );
  ColVariable * nx;
  LB0Constraint * nUB = nullptr;
  if( arc == get_NArcs() ) {  // the new arc is physically constructed
   // create the new variable
   std::list< ColVariable > na;
   na.emplace_back( this , ColVariable::kNonNegative );
   nx = &(na.back());

   Block::add_dynamic_variables( dx , na , ampar );  // now add it

   // add the new coefficient in the objective
   if( AR & HasObj )
    get_lfo()->add_variable( nx , cst , ampar );

   if( ( cap < Inf< FNumber >() ) && ( AR & HasFlw ) && ( ! ( AR & HasBnd ) ) )
    throw( std::logic_error( "cannot set finite capacity" ) );

   if( AR & HasBnd ) {  // construct new arc capacity constraint
    std::list< LB0Constraint > nub;
    nub.emplace_back( this , nx );
    nUB = &(nub.back());

    Block::add_dynamic_constraints( dUB , nub , ampar );  // now add it
    }

   SN[ arc ] = EN[ arc ] = Inf< Index >();  // ensure sn and en change
   }
  else {  // the arc is just inserted in a previously deleted slot
   nx = i2p_x( arc );  // recover pointer to the flow Variable

   nx->is_fixed( false , ampar );   // un-fix the Variable

   // recover pointer to the bound Constraint (if any)
   if( AR & HasBnd )
    nUB = i2p_ub( arc );

   // change the cost coefficient in the objective (if any)
   if( AR & HasObj )
    get_lfo()->modify_coefficient( arc , cst , ampar );

   if( AR & HasFlw ) {
    // in the (very likely) case that the start (end) node is different
    // from before, remove the contribution to the previous start (end)
    // node flow constraint
    if( sn != SN[ arc ] ) {
     auto osn = get_lfc( i2p_e( SN[ arc ] - 1 ) );
     auto i = osn->is_active( nx );
     if( i >= osn->get_num_active_var() )
      throw( std::logic_error(
	       "created arc not present in previous sn flow constraint" ) );
     osn->remove_variable( i , ampar );
     }
    if( en != EN[ arc ] ) {
     auto oen = get_lfc( i2p_e( EN[ arc ] - 1 ) );
     auto i = oen->is_active( nx );
     if( i >= oen->get_num_active_var() )
      throw( std::logic_error(
	       "created arc not present in previous en flow constraint" ) );
     oen->remove_variable( i , ampar );
     }
    }
   }  // end( inserting in the middle )

  if( nUB )
   nUB->set_rhs( cap , ampar );  // set new arc capacity

  if( AR & HasFlw ) {
   // set contribution to flow constraint, unless in the (very unlikely)
   // case that the start (end) node is the same as before
   if( sn != SN[ arc ] )
    get_lfc( i2p_e( sn - 1 ) )->add_variable( nx , -1 , ampar );
   if( en != EN[ arc ] )
    get_lfc( i2p_e( en - 1 ) )->add_variable( nx ,  1 , ampar );
   }

  close_if_needed( ampar , AR & ( HasFlw | HasObj ) ? 4 : 1 );
  }

 // change the physical representation- - - - - - - - - - - - - - - - - - - -
 if( not_dry_run( issueMod ) ) {
  C[ arc ] = cst;  // set new arc cost

  if( U.empty() && ( cap < Inf< FNumber >() ) )
   U.assign( get_MaxNArcs() , Inf< FNumber >() );

  if( ! U.empty() )
   U[ arc ] = cap;  // set new arc capacity

  SN[ arc ] = sn;  // set Start Node
  EN[ arc ] = en;  // set End Node
  }

 if( arc == get_NArcs() )
  ++NArcs;  // increase arc count

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
			    MCFBlockMod::eAddArc , Range( arc , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 return( arc );

 }  // end( MCFBlock::add_arc )

/*--------------------------------------------------------------------------*/

void MCFBlock::remove_arc( Index arc ,
			   ModParam issueMod , ModParam issueAMod )
{
 if( ( arc < get_NStaticArcs() ) || ( arc >= get_NArcs() ) )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( is_deleted( arc ) )  // arc deleted already
  return;                 // nothing to do

 // change the abstract representation- - - - - - - - - - - - - - - - - - - -
 // in the meantime, if so instructed also issue abstract Modification(s)
 // note that this is *always* done, unless issueAMod says this is a dry
 // run, because at least the BlockModAD corresponding to deleting the
 // Variable(s), or fixing them, is always issued since the Variable are
 // always present

 if( not_dry_run( issueAMod ) ) {
  not_ModBlock( issueAMod );
  if( arc == get_NArcs() - 1 ) {  // removing the last arc(s)
   Index rmvdarcs = 1;  // how many arcs are removed in the end
   auto ampar = open_if_needed( issueAMod ,
				AR & ( HasFlw | HasObj ) ? 4 : 1 );
   auto ritdx = dx.rbegin();  // reverse iterator into dx

   // pointer to LinearFunction in the objective (if any)
   LinearFunction * lfo = nullptr;
   if( AR & HasObj )
    lfo = get_lfo();

   // scan from the end backwards, eliminate all deleted arcs
   for( Index ai = arc ; ritdx != dx.rend() ; ++rmvdarcs , --ai ) {
    auto sn = SN[ ai ]; sn--;
    auto en = EN[ ai ]; en--;
    auto rxi = &(*(ritdx++));

    // delete contribution to objective (if any)
    if( lfo )
     lfo->remove_variable( ai , ampar );

    // delete contribution to flow constraint (if any)
    if( AR & HasFlw ) {
     auto snc = get_lfc( i2p_e( sn ) );
     auto sni = snc->is_active( rxi );
     if( sni >= snc->get_num_active_var() )
      throw( std::logic_error( "x variable not active in flow constraint" ) );
     snc->remove_variable( sni , ampar );
     auto enc = get_lfc( i2p_e( en ) );
     auto eni = enc->is_active( rxi );
     if( eni >= enc->get_num_active_var() )
      throw( std::logic_error( "x variable not active in flow constraint" ) );
     enc->remove_variable( eni , ampar );
     }

    // all remaining dynamic arcs have been removed
    if( rmvdarcs >= get_NArcs() - get_NStaticArcs() )
     break;  // all done

    // the first remaining arc is not deleted
    if( ! is_deleted( get_NArcs() - rmvdarcs - 1 ) )
     break;  // all done
    }

   // define the range of removed stuff
   Range range( get_NArcs() - get_NStaticArcs() - rmvdarcs ,
		get_NArcs() - get_NStaticArcs() );
   
   // now actually remove and clear the UB Constraint(s) (if any)
   // do this before removing the flow Variable(s), so that if they are
   // processed in FIFO order it is seen before
   if( AR & HasBnd )
    Block::remove_dynamic_constraints( dUB , range , ampar );

   // now actually remove the flow Variable(s) (if any)
   Block::remove_dynamic_variables( dx , range , ampar );

   close_if_needed( ampar , AR & ( HasFlw | HasObj ) ? 4 : 1 );   
   }
  else {  // deleting one arc in the middle
   auto rx = i2p_x( arc );
   rx->set_value( 0 );                // set the Variable to 0
   rx->is_fixed( true , issueAMod );  // fix it
   }
  }
 else  // at the very least ensure the value is 0
  i2p_x( arc )->set_value( 0 );

 f_cond_lower = dNAN;  // reset conditional bounds

 // change the physical representation- - - - - - - - - - - - - - - - - - - -
 if( not_dry_run( issueMod ) ) {
  Index rmvdarcs = 1;  // how many arcs are removed in the end

  if( arc == NArcs - 1 ) {    // deleted last arc(s)
   for( --NArcs ; std::isnan( C[ NArcs - 1 ] ) ; --NArcs , ++rmvdarcs )
    ;  // decrease arc count
   }
  else                       // deleted one arc in the middle
   C[ arc ] = std::numeric_limits< CNumber >::quiet_NaN();

  if( issue_pmod( issueMod ) )  // issue "physical Modification"- - - - - - -
   Block::add_Modification( std::make_shared< MCFBlockRngdMod >( this ,
				    MCFBlockMod::eRmvArc ,
				    Range( arc - rmvdarcs + 1 , arc + 1 ) ) ,
			    Observer::par2chnl( issueMod ) );
  }

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( MCFBlock::remove_arc )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void MCFBlock::guts_of_destructor( void )
{
 // clear() all Constraint to ensure that they do not bother to un-register
 // themselves from Variable that are going to be deleted anyway

 // clear the bound constraints
 Constraint::clear( UB );   // static
 Constraint::clear( dUB );  // dynamic

 // clear the flow conservation constraints
 Constraint::clear( E );   // static
 Constraint::clear( dE );  // dynamic

 c.clear();  // clear the Objective

 // delete all Variable
 dx.clear();  // dynamic
 x.clear();   // static

 // explicitly reset all Constraint and Variable
 // this is done for the case where this method is called prior to re-loading
 // a new instance: if not, the new representation would be added to the
 // (no longer current) one
 reset_static_constraints();
 reset_static_variables();
 reset_dynamic_constraints();
 reset_dynamic_variables();
 reset_objective();

 AR = 0;

 }  // end( MCFBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/

void MCFBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since MCFBlock is a "leaf" Block (has no sub-Block), this
  * method does not have to deal with GroupModification since these are
  * produced by Block::add_Modification(), but this method is called
  * *before* that one is.
  *
  * As an important consequence,
  *
  *   THE STATE OF THE DATA STRUCTURE IN MCFBlock WHEN THIS METHOD IS
  *   EXECUTED IS PRECISELY THE ONE IN WHICH THE Modification WAS ISSUED:
  *   NO COMPLICATED OPERATIONS (Variable AND/OR Constraint BEING
  *   ADDED/REMOVED ...) CAN HAVE BEEN PERFORMED IN THE MEANTIME
  *
  * This assumption drastically simplifies some of the logic here.*/

 // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod ) ) {
  if( ! ( AR & HasObj ) )
   throw( std::invalid_argument( "Modification to non-constructed Objective"
				 ) );

  auto lfo = static_cast< LinearFunction * const >( tmod->function() );
  if( static_cast< LinearFunction * const >( c.get_function() ) != lfo )
   throw( std::invalid_argument( "Modification to non-Objective" ) );

  // note: in the following we can assume that the Range in tmod is
  //       precisely the one we have to use since no Variable can have
  //       been added or deleted, which saves *a lot* of trouble

  if( tmod->range().second == tmod->range().first + 1 )
   // changing one cost only
   chg_cost( lfo->get_coefficient( tmod->range().first ) ,
	     tmod->range().first , make_par( eNoBlck , chnl ) , eDryRun );
  else {                            // changing many costs at once
   Vec_CNumber NC( tmod->range().second - tmod->range().first );
   auto NCit = NC.begin();
   for( Index i = tmod->range().first ; i < tmod->range().second ; )
    *(NCit++) = lfo->get_coefficient( i++ );

   chg_costs( NC.begin() , tmod->range() ,
	      make_par( eNoBlck , chnl ) , eDryRun );
   }

  return;
  }

 // C05FunctionModLinSbst - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< C05FunctionModLinSbst * >( mod ) ) {
  if( ! ( AR & HasObj ) )
   throw( std::invalid_argument( "Modification to non-constructed Objective"
				 ) );

  auto lfo = static_cast< LinearFunction * const >( tmod->function() );
  if( static_cast< LinearFunction * const >( c.get_function() ) != lfo )
   throw( std::invalid_argument( "Modification to non-Objective" ) );

  // note: in the following we can assume that the Subset in tmod is
  //       precisely the one we have to use since no Variable can have
  //       been added or deleted, which saves *a lot* of trouble
  // note: chg_costs() owns subset, so a copy has to be made

  Vec_CNumber NC( tmod->subset().size() );
  auto NCit = NC.begin();
  for( auto i : tmod->subset() )
   *(NCit++) = lfo->get_coefficient( i++ );

  chg_costs( NC.begin() , Subset( tmod->subset() ) , true ,
	     make_par( eNoBlck , chnl ) , eDryRun );
  return;
  }

 // RowConstraintMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< RowConstraintMod * >( mod ) ) {
  if( ! ( AR & HasFlw ) )
   throw( std::invalid_argument(
			    "Modification to non-constructed Constraint" ) );

  if( tmod->type() == RowConstraintMod::eChgRHS ) {
   auto cp = dynamic_cast< LB0Constraint * const >( tmod->constraint() );
   if( ! cp )
    throw( std::invalid_argument( "invalid Modification to Constraint" ) );

   chg_ucap( cp->get_rhs() , p2i_ub( cp ) ,
	     make_par( eNoBlck , chnl ) , eDryRun );
   return;
   }

  if( tmod->type() == RowConstraintMod::eChgBTS ) {
   auto cp = static_cast< FRowConstraint * const >( tmod->constraint() );
   if( ! cp )
    throw( std::invalid_argument( "invalid Modification to Constraint" ) );

   chg_dfct( cp->get_rhs() , p2i_e( cp ) ,
	     make_par( eNoBlck , chnl ) , eDryRun );
   return;
   }

  throw( std::invalid_argument( "illegal Modification to Constraint" ) );
  }

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< VariableMod * >( mod ) ) {
  auto xi = dynamic_cast< const ColVariable * >( tmod->variable() );
  if( ! xi )
   throw( std::logic_error( "Modification to wrong type of Variable" ) );
  if( ( xi->get_type() != ColVariable::kNonNegative ) &&
      ( xi->get_type() != ColVariable::kNatural ) )
   throw( std::logic_error( "changing type of flow Variable not allowed" ) );
   
  auto i = p2i_x( xi );
  if( std::isnan( C[ i ] ) )
   throw( std::logic_error( "[un]fixing deleted arc not allowed" ) );
  if( xi->is_fixed() )
   close_arc( i , make_par( eNoBlck , chnl ) , eDryRun );
  else
   open_arc( i , make_par( eNoBlck , chnl ) , eDryRun );

  return;
  }

 throw( std::invalid_argument( "unsupported Modification to MCFBlock" ) );

 }  // end( MCFBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/

void MCFBlock::compute_conditional_bounds( void )
{
 f_cond_lower = f_cond_upper = 0;

 auto tC = C.begin();
 auto tU = U.begin();

 for( ; tC < C.end() ; ++tC , ++tU ) {
  if( *tC == 0 )
   continue;

  if( *tC < 0 ) {
   if( *tU == Inf< FNumber >() ) {
    f_cond_lower = -Inf< double >();
    break;
    }
   else
    f_cond_lower += *tC * (*tU);
   }
  else
   if( *tU == Inf< FNumber >() ) {
    f_cond_upper = Inf< double >();
    break;
    }
   else
    f_cond_upper += *tC * (*tU);
   }

 if( f_cond_lower > -Inf< double >() ) {
  for( ; tC < C.end() ; ++tC , ++tU )
   if( *tC < 0 ) {
    if( *tU == Inf< FNumber >() ) {
     f_cond_lower = -Inf< double >();
     break;
     }
    else
     f_cond_lower += *tC * (*tU);
    }
  }

 if( f_cond_upper < Inf< double >() ) {
  for( ; tC < C.end() ; ++tC , ++tU )
   if( *tC > 0 ) {
    if( *tU == Inf< FNumber >() ) {
     f_cond_upper = Inf< double >();
     break;
     }
    else
     f_cond_upper += *tC * (*tU);
    }
  }
 }  // end( MCFBlock::compute_conditional_bounds )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

void MCFBlock::CheckAbsVSPhys( void )
{
 // check that the (part that has actually been constructed of the) abstract
 // representation coincides with the physical representation

 // check variables, these are always there - - - - - - - - - - - - - - - - -
 if( x.size() != get_NStaticArcs() )
  std::cerr << "x.size() != NStaticArcs" << std::endl;

 if( dx.size() < get_NArcs() - get_NStaticArcs() )
  std::cerr << "dx.size() too small" << std::endl;

 // check flow constraints- - - - - - - - - - - - - - - - - - - - - - - - - -
 if( AR & HasFlw ) { 
  if( E.size() != get_NStaticNodes() )
   std::cerr << "E.size() != NStaticNodes" << std::endl;

  if( dE.size() < get_NNodes() - get_NStaticNodes() )
   std::cerr << "dE.size() too small" << std::endl;
  
  Subset NRIncid( get_NNodes() , 0 );

  Index objbnd = 0;          // number of actives between obj and bound
  if( AR & HasBnd )
   ++objbnd;
  if( AR & HasObj )
   ++objbnd;
  Index expnc = 2 + objbnd;  // total number of active stuff per variable

  // static arcs
  Index a = 0;
  for( auto xi = x.begin() ; a < get_NStaticArcs() ; ++a , ++xi ) {
   if( is_deleted( a ) ) {
    std::cerr << "static arc " << a << " deleted" << std::endl;
    continue;
    }
   
   if( xi->get_num_active() != expnc )
    std::cerr << "arc " << a << " active in " << xi->get_num_active()
	      << " constraints" << std::endl;

   Index sna = SN[ a ];
   if( ! sna )
    std::cerr << "SN[ " << a << " ] == 0" << std::endl;
   else
    --sna;
   if( sna >= get_NArcs() )
    std::cerr << "SN[ " << a << " ] == " << sna << " >= |A|" << std::endl;

   if( ( SN[ a ] > 0 ) && ( sna < get_NArcs() ) ) {
    ++NRIncid[ sna ];
    auto snc = get_lfc( i2p_e( sna ) );
    auto sni = snc->is_active( &(*xi) );
    if( sni >= snc->get_num_active_var() )
     std::cerr << "static arc " << a
	       << " absent in flow constraint for SN[ a ] == "
	       << sna << std::endl;
    }
    
   Index ena = EN[ a ];
   if( ! ena )
    std::cerr << "EN[ " << a << " ] == 0" << std::endl;
   else
    --ena;
   if( ena >= get_NArcs() )
    std::cerr << "EN[ " << a << " ] == " << ena << " >= |A|" << std::endl;

   if( ( EN[ a ] > 0 ) && ( ena < get_NArcs() ) ) {
    ++NRIncid[ ena ];
    auto enc = get_lfc( i2p_e( ena ) );
    auto eni = enc->is_active( &(*xi) );
    if( eni >= enc->get_num_active_var() )
     std::cerr << "static arc " << a
	       << " absent in flow constraint for EN[ a ] == "
	       << ena << std::endl;
    }
   }

  // dynamic arcs
  for( auto xi = dx.begin() ; a < get_NArcs() ; ++a , ++xi ) {
   if( xi->get_num_active() != expnc )
    std::cerr << "arc " << a << " active in " << xi->get_num_active()
	      << " constraints" << std::endl;

   Index sna = SN[ a ];
   if( ! sna )
    std::cerr << "SN[ " << a << " ] == 0 " << std::endl;
   else
    --sna;
   if( sna >= get_NArcs() )
    std::cerr << "SN[ " << a << " ] == " << sna << " >= |A!" << std::endl;

   if( ( SN[ a ] > 0 ) && ( sna < get_NArcs() ) ) {
    ++NRIncid[ sna ];
    auto snc = get_lfc( i2p_e( sna ) );
    auto sni = snc->is_active( &(*xi) );
    if( sni >= snc->get_num_active_var() )
     std::cerr << "static arc " << a
	       << " absent in flow constraint for SN[ a ] == "
	       << sna << std::endl;
    }

   Index ena = EN[ a ];
   if( ! ena )
    std::cerr << "EN[ " << a << " ] == 0 " << std::endl;
   else
    --ena;
   if( ena >= get_NArcs() )
    std::cerr << "EN[ " << a << " ] == " << ena << " >= |A!" << std::endl;

   if( ( EN[ a ] > 0 ) && ( ena < get_NArcs() ) ) {
    ++NRIncid[ ena ];
    auto enc = get_lfc( i2p_e( ena ) );
    auto eni = enc->is_active( &(*xi) );
    if( eni >= enc->get_num_active_var() )
     std::cerr << "static arc " << a
	       << " absent in flow constraint for EN[ a ] == "
	       << ena << std::endl;
    }
   }

  if( HasDynamicX() && is_deleted( get_NArcs() - 1 ) )
   std::cerr << "last dynamic arc is deleted" << std::endl;

  // static nodes
  Index n = 0;
  for( auto ni = E.begin() ; n < get_NStaticNodes() ; ++n , ++ni ) {
   auto lni = get_lfc( &(*ni) );
   if( lni->get_num_active_var() != NRIncid[ n ] )
    std::cerr << "active variables in static flow constraint " << n
	      << " == " << lni->get_num_active_var()
	      << " do not match with incident arcs " << NRIncid[ n ]
	      << std::endl;
   }
   
  // dynamic nodes
  for( auto ni = dE.begin() ; n < get_NNodes() ; ++n , ++ni ) {
   auto lni = get_lfc( &(*ni) );
   if( lni->get_num_active_var() != NRIncid[ n ] )
    std::cerr << "active variables in dynamic flow constraint " << n
	      << " == " << lni->get_num_active_var()
	      << " do not match with incident arcs " << NRIncid[ n ]
	      << std::endl;
   }
  }  // end( if( AR & HasFlw ) )

 // check bound constraints - - - - - - - - - - - - - - - - - - - - - - - - -
 if( AR & HasBnd ) {
  // static bounds
  Index a = 0;
  auto UBi = UB.begin();
  for( auto xi = x.begin() ; xi != x.end() ; ++a , ++xi , ++UBi ) {
   if( UBi->is_active( &(*xi) ) >= UBi->get_num_active_var() )
    std::cerr << "static arc " << a << " absent in bound constraint"
	      << std::endl;
   }

 // dynamic bounds
  auto dUBi = dUB.begin();
  for( auto xi = dx.begin() ; xi != dx.end() ; ++a , ++xi , ++dUBi ) {
   if( is_deleted( a ) )
    continue;

   if( dUBi->is_active( &(*xi) ) >= dUBi->get_num_active_var() )
    std::cerr << "dynamic arc " << a << " absent in bound constraint"
	      << std::endl;
   }
  }  // end( if( AR & HasBnd ) )

 // check objective - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( AR & HasObj ) {
  auto lfo = get_lfo();
  if( lfo->get_num_active_var() != get_NArcs() )
   std::cerr << "objective has " << lfo->get_num_active_var()
	     << " variables while |A| = " << get_NArcs() << std::endl;

  Index a = 0;
  for( auto & xi : x ) {
   if( lfo->is_active( & xi ) >= get_NArcs() )
    std::cerr << "static arc " << a << " absent from objective " << std::endl;
   ++a;
   }

  for( auto & xi : dx ) {
   if( lfo->is_active( & xi ) >= get_NArcs() )
    std::cerr << "dynamic arc " << a << " absent from objective "
	      << std::endl;
   ++a;
   }
  }  // end( if( AR & HasObj ) )
 }  // end( MCFBlock::CheckAbsVSPhys )

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS OF MCFSolution ------------------------*/
/*--------------------------------------------------------------------------*/

void MCFSolution::deserialize( const netCDF::NcGroup & group )
{
 netCDF::NcDim na = group.getDim( "NumArcs" );
 if( na.isNull() )
  v_x.clear();
 else {
  netCDF::NcVar fs = group.getVar( "FlowSolution" );
  if( fs.isNull() )
   v_x.clear();
  else {
   v_x.resize( na.getSize() );
   fs.getVar( v_x.data() );
   }
  }

 netCDF::NcDim nn = group.getDim( "NumNodes" );
 if( nn.isNull() )
  v_pi.clear();
 else {
  netCDF::NcVar ps = group.getVar( "Potentials" );
  if( ps.isNull() )
   v_pi.clear();
  else {
   v_pi.resize( nn.getSize() );
   ps.getVar( v_pi.data() );
   }
  }
 }  // end( MCFSolution::deserialize )

/*--------------------------------------------------------------------------*/

void MCFSolution::read( const Block * block )
{
 auto MCFB = dynamic_cast< const MCFBlock * >( block );
 if( ! MCFB )
  throw( std::invalid_argument( "block is not a MCFBlock" ) );

 // read flows- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_x.empty() ) {
  if( v_x.size() < MCFB->get_NArcs() )
   v_x.resize( MCFB->get_NArcs() );

  MCFB->get_x( v_x.begin() );
  }

 // read potentials - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( MCFB->E.empty() && MCFB->dE.empty() )  // no potentials available
  return;

 if( ! v_pi.empty() ) {
  if( v_pi.size() < MCFB->get_NNodes() )
   v_pi.resize( MCFB->get_NNodes() );

  MCFB->get_pi( v_pi.begin() );
  }
 }  // end( MCFSolution::read )

/*--------------------------------------------------------------------------*/

void MCFSolution::write( Block * block ) 
{
 auto MCFB = dynamic_cast< MCFBlock * >( block );
 if( ! MCFB )
  throw( std::invalid_argument( "block is not a MCFBlock" ) );

 // write flows - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_x.empty() ) {
  if( v_x.size() < MCFB->get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  MCFB->set_x( v_x.begin() );
  }

 // write potentials- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( v_pi.empty() )  // no potentials to write
  return;

 if( MCFB->E.empty() && MCFB->dE.empty() )  // no Constraint to write to
  return;

 if( v_pi.size() < MCFB->get_NStaticNodes() )
  throw( std::invalid_argument( "incompatible potential size" ) );

 MCFB->set_pi( v_pi.begin() );

 // write reduced costs (if any)- - - - - - - - - - - - - - - - - - - - - - -
  
 if( MCFB->UB.empty() && MCFB->dUB.empty() )  // no bounds to write to
  return;

 MCFBlock::Index i = 0;

 // static part
 for( auto ubi = MCFB->UB.begin() ; ubi != MCFB->UB.end() ; ++i )
  (ubi++)->set_dual( MCFB->get_C( i ) + v_pi[ MCFB->SN[ i ] - 1 ]
		                      - v_pi[ MCFB->EN[ i ] - 1 ] );
 // dynamic part
 for( auto dubi = MCFB->dUB.begin() ;
      ( dubi != MCFB->dUB.end() ) && ( i < MCFB->get_NArcs() ) ; ++i )
  (dubi++)->set_dual( MCFB->get_C( i ) + v_pi[ MCFB->SN[ i ] - 1 ]
		                       - v_pi[ MCFB->EN[ i ] - 1 ] );
 }  // end( MCFSolution::write )

/*--------------------------------------------------------------------------*/

void MCFSolution::serialize( netCDF::NcGroup & group ) const
{
 // always call the method of the base class first
 Solution::serialize( group );

 std::vector< size_t > startp = { 0 };

 if( ! v_x.empty() ) {
  netCDF::NcDim na = group.addDim( "NumArcs" , v_x.size() );

  std::vector< size_t > countpa = { v_x.size() };

  ( group.addVar( "FlowSolution" , netCDF::NcDouble() , na ) ).putVar(
					      startp , countpa , v_x.data() );
  }

 if( v_pi.empty() )
  return;

 netCDF::NcDim nn = group.addDim( "NumNodes" ,  v_pi.size() );
 std::vector< size_t > countpn = { v_pi.size() };
 ( group.addVar( "Potentials" , netCDF::NcDouble() , nn ) ).putVar(
					     startp , countpn , v_pi.data() );
 
 }  // end( MCFSolution::serialize )

/*--------------------------------------------------------------------------*/

MCFSolution * MCFSolution::scale( double factor ) const
{
 auto * sol = MCFSolution::clone( true );

 if( ! v_x.empty() )
  for( MCFBlock::Index i = 0 ; i < v_x.size() ; ++i )
   sol->v_x[ i ] = v_x[ i ] * factor;

 if( ! v_pi.empty() )
  for( MCFBlock::Index i = 0 ; i < v_pi.size() ; ++i )
   sol->v_pi[ i ] = v_pi[ i ] * factor;

 return( sol );

 }  // end( MCFSolution::scale )

/*--------------------------------------------------------------------------*/

void MCFSolution::sum( const Solution * solution , double multiplier )
{
 auto MCFS = dynamic_cast< const MCFSolution * >( solution );
 if( ! MCFS )
  throw( std::invalid_argument( "solution is not a MCFSolution" ) );

 if( ! v_x.empty() ) {
  if( v_x.size() != MCFS->v_x.size() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  for( MCFBlock::Index i = 0 ; i < v_x.size() ; ++i )
   v_x[ i ] += MCFS->v_x[ i ] * multiplier;
  }

 if( ! v_pi.empty() ) {
  if( v_pi.size() != MCFS->v_pi.size()  )
   throw( std::invalid_argument( "incompatible potential size" ) );

  for( MCFBlock::Index i = 0 ; i < v_pi.size() ; ++i )
   v_pi[ i ] += MCFS->v_pi[ i ] * multiplier;
  }
 }  // end( MCFSolution::sum )

/*--------------------------------------------------------------------------*/

MCFSolution * MCFSolution::clone( bool empty ) const
{
 auto *sol = new MCFSolution();

 if( empty ) {
  if( ! v_x.empty() )
   sol->v_x.resize( v_x.size() );

  if( ! v_pi.empty() )
   sol->v_pi.resize( v_pi.size() );
  }
 else {
  sol->v_x = v_x;
  sol->v_pi = v_pi;
  }

 return( sol );

 }  // end( MCFSolution::clone )

/*--------------------------------------------------------------------------*/
/*----------------------- End File MCFBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
