 FC=gfortran-5  -Wall -O3 -w 
 FLAGS=-frecord-marker=4 -DVERBOSE=1 -x f95-cpp-input  
 FLAGS77 = -frecord-marker=4 
 OBJ = diskparams.o derived_types.o diskgravmjp.o user_forces.o user_initialisation.o mercury_MJP_1.8.o

potpert: $(OBJ)
	$(FC) $(LFLAGS) -o mercury $(OBJ)

clean:
	rm -f *.o *.mod

veryclean:
	rm -f *.o *.mod potpert mercury *.aei *.tmp *.dmp *.out *.aei *~

%.o:%.f90
	$(FC) $(FLAGS) -c $^ -o $@
%.o:%.f
	$(FC) $(FLAGS77) -c $^ -o $@

