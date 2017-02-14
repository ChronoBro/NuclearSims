
/****************************************************
This is my attempt(DH) at consolidating and improving
     the Nuclear Reaction Monte Carlo Simulations
     	 done at WUSTL
	 
	      created by,
	      	      Bob Charity (Much Love)

****************************************************/

.h files are now in /include
.cpp files are now in /src
.o files created are now in /objs (Feel free to change this Cole)
.root files should be put in /root (but like, w/e man)

To add a new "object(class)" go to the Makefile and edit $OBJECTS and add "objs/(name of object).o"
to the list and create (src\include)(name of object).(cpp\h)

If you want to simulate a new reaction create a new "target" in the Makefile, use example format of sim34Ca,
and edit the reaction specifics in src/(simName).cpp

The simulation will be stored in /sims with simName after you make (if in example format)

Feel free to edit the /src files but try not to break anything that currently works so we can consolidate.
The idea being we keep a consistent set of objects so all we have to do is create new sim files.

Please add as much as you want!

I think we should put this to github