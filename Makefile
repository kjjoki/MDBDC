# Makefile for the multiobjective double bundle method MDBDC for constrained DC optimization

FF = gfortran
#FFLAGS = -fbounds-check -Wall
#OPEN =  
OPEN =  -fopenmp 
#RM = del   #in windows
RM = rm    #in unix/linux/mac

all: tmdbdc

tmdbdc: constants.o bundle2.o functions.o bundle1.o mdbdc.o tmdbdc.o plqdf1.o 
	$(FF) -o tmdbdc $(FFLAGS) $(OPEN) constants.o bundle2.o functions.o bundle1.o mdbdc.o tmdbdc.o plqdf1.o 

constants.mod: constants.o constants.f95
	$(FF) -c $(FFLAGS) $(OPEN) constants.f95
	
constants.o: constants.f95
	$(FF) -c $(FFLAGS) $(OPEN) constants.f95
		
bundle2.mod: constants.mod bundle2.o bundle2.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f95 
	
bundle2.o: constants.mod bundle2.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f95 	

functions.mod: constants.mod functions.o functions.f95
	$(FF) -c $(FFLAGS) $(OPEN) functions.f95 
	
functions.o: constants.mod functions.f95
	$(FF) -c $(FFLAGS) $(OPEN) functions.f95 
	
bundle1.mod: constants.mod functions.mod bundle1.o bundle1.f95 
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f95
	
bundle1.o: constants.mod functions.mod bundle1.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f95 

mdbdc.mod: constants.mod bundle2.mod functions.mod bundle1.mod mdbdc.o mdbdc.f95 
	$(FF) -c $(FFLAGS) $(OPEN) mdbdc.f95	 
	
mdbdc.o: constants.mod bundle2.mod functions.mod bundle1.mod mdbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) mdbdc.f95 
	
tmdbdc.o: constants.mod bundle2.mod functions.mod bundle1.mod mdbdc.mod tmdbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) tmdbdc.f95 

plqdf1.o: plqdf1.f
	$(FF) -c $(FFLAGS) $(OPEN) plqdf1.f 

clean:	
	rm tmdbdc constants.mod constants.o bundle1.mod bundle1.o bundle2.mod bundle2.o functions.mod functions.o mdbdc.mod mdbdc.o tmdbdc.o plqdf1.o  
	echo Clean done