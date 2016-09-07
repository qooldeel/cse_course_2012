CC = c++
INCLDIR = ..
INCL = -I$(INCLDIR)
ERR = -Wfatal-errors
CC_FLAGS = -Wall -O3 $(ERR) -pedantic 
## profiling flag: uncomment the following whenever you want to profile
PROF = #-pg #-funroll-loops

###  USAGE:
## 1.) replace FILENAME = with the name of the file you wanna to compile 
## 2.) go back to terminal and execute: make clean; make; ./file
##     where 'file' is the executable from EXEC 

## ========== TODO: replace FILENAME (without suffix) by specific name of file
##                  you're going to compile, e.g. vec, blatt6, etc.
##  CAUTION: BEWARE OF WHITE CHARS IN FILENAME!!! (otherwise your file might
##           not be found ;) ) 
FILENAME = expvsimpl
##============================================================================

SUFFIX=.cc

EXEC = $(FILENAME)

MAIN_PROG = $(FILENAME).o    
OBJS = $(MAIN_PROG)

LIBS = 


## mind the order here!!
$(EXEC): $(OBJS)
	@echo "LINKING & BUILDING EXECUTABLE..."
	$(CC)  $(CC_FLAGS) -o $(EXEC) $(OBJS) $(PROF) $(LIBS)
	@echo "Done! ;)"

$(OBJS):
	@echo "COMPILING..."
	$(CC) $(CC_FLAGS) $(INCL) $(PROF) -c $*$(SUFFIX)   


#proper clean up when invoking 'make clean'
.PHONY: clean
clean:
	rm -rf $(OBJS) $(EXEC) *.out