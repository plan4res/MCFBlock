##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of MCFBlock / MCFSolver                                         #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by MCFBlock / MCFSolver, i.e., core SMS++ and           #
#   MCFClass. These are *not* copied into $(MCFBkINC): adding those -I       #
#   directives to the compile commands will have to done by whatever "main"  #
#   makefile is using this. Analogously, any external library and the        #
#   corresponding -L< libdirs > will have to be added to the final linking   #
#   command by  whatever "main" makefile is using this.                      #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by MCFBlock / MCFSolver.                                            #
#                                                                            #
#   The makefile defines internally (cf. MCFClssSlvr below) which            #
#   MCFSolver< :MCFClass > will be available.                                #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(libMCFClINC) = the -I$( MCFClass library )                     #
#           $(libMCFClOBJ) = the MCFClass library                            #
#           $(MCFBkSDR)    = the directory where the source is               #
#                                                                            #
#   Output: $(MCFBkOBJ)    = the final object(s) / library                   #
#           $(MCFBkH)      = the .h files to include                         #
#           $(MCFBkINC)    = the -I$( source directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# define the set of MCFSolver< :MCFClass > that will be available by
# uncommenting the -DHAVE_* below corresponding to the :MCFClass; see
# MCFSolver.h for details. note that, obviously, the :MCFClass selected
# here must have been compiled in the MCFClass library
MCFClssSlvr = -DHAVE_MFSMX -DHAVE_CPLEX -DHAVE_RELAX
# -DHAVE_CSCL2 -DHAVE_MFZIB -DHAVE_SPTRE

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

MCFBkOBJ = $(MCFBkSDR)/obj/MCFBlock.o $(MCFBkSDR)/obj/MCFSolver.o

MCFBkINC = -I$(MCFBkSDR)/include

MCFBkH   = $(MCFBkSDR)/include/MCFBlock.h $(MCFBkSDR)/include/MCFSolver.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MCFBkOBJ) $(MCFBkSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(MCFBkSDR)/obj/MCFBlock.o: $(MCFBkSDR)/src/MCFBlock.cpp \
	$(MCFBkSDR)/include/MCFBlock.h $(SMS++OBJ)
	$(CC) -c $(MCFBkSDR)/src/MCFBlock.cpp -o $@ \
	$(MCFBkINC) $(SMS++INC) $(SW)

$(MCFBkSDR)/obj/MCFSolver.o: $(MCFBkSDR)/src/MCFSolver.cpp $(MCFBkH) \
	$(SMS++OBJ) $(libMCFClOBJ) 
	$(CC) -c $(MCFBkSDR)/src/MCFSolver.cpp -o $@ \
	$(MCFBkINC) $(SMS++INC) $(libMCFClINC) $(MCFClssSlvr) $(SW)

########################## End of makefile ###################################
