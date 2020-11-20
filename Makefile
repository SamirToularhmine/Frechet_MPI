MPICXX = mpic++
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o)
EXEC=frechet

%.o: %.cpp
	$(MPICXX) -o $@ -c $< 

frechet: $(OBJ) 
	$(MPICXX) -o $@ $^

run: frechet
	mpirun -np 4 ./frechet traj1.txt traj2.txt

clean:
	@rm -rf *.o

cleanall: clean
	@rm -rf $(EXEC)

