// tests-main.cpp
// #define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

int argc;
char** argv;
extern int GRIDNUM;
extern int globtestnum;
extern int SURFOPT;
extern double RADIUS;

std::string INPUTFILE;// protein file used for testing

int main( int argc, char* argv[] )
{
  Catch::Session session; // There must be exactly one instance
  
   // Some user variable you want to be able to set
  // Build a new parser on top of Catch's
  using namespace Catch::clara;
  auto cli 
    = session.cli() // Get Catch's composite command line parser
    | Opt( INPUTFILE, "file" ) // bind variable to a new option, with a hint string
        ["-pf"]["--file"]    // the option names it will respond to
        ("file name")        // description string for the help output
    | Opt( SURFOPT, "SURFOPT" )
        ["-s"]["--SURFOPT"]
        ("surface option")
    | Opt( GRIDNUM, "GRIDNUM" )
        ["-n"]["--GRIDNUM"]
        ("grid number")
    | Opt( RADIUS, "radius" )
        ["-R"]["--RADIUS"]
        ("radius")
    | Opt( globtestnum, "globtestnum" )
        ["-T"]["--globtestnum"]
        ("test");  
        
  // Now pass the new composite back to Catch so it uses that
  session.cli( cli ); 
  
  // Let Catch (using Clara) parse the command line
  int returnCode = session.applyCommandLine( argc, argv );
  if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

  return session.run();
}