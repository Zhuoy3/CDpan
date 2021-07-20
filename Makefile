OBJ1=./src/main.o ./src/cmdwrapper.o
OUT1=./bin/CDpan.bin
INCL=-I./incl/ -I./lib/
PYSO=-lstdc++
install: $(OBJ1)
	gcc -Wall -O2 ${INCL} -o ${OUT1} $(OBJ1) ${PYSO}

./src/main.o: ./src/main.c ./incl/*.h
	gcc -Wall -O2 ${INCL} -c -o $@ $<

./src/cmdwrapper.o: ./src/cmdwrapper.cpp ./incl/cmdwrapper.h
	g++ -Wall -O2 ${INCL} -I./lib/ -c -o $@ $<

.PHONY: clean
clean:
	-rm -f $(OBJ1)

.PHONY: unstall
unstall:
	-rm -f ${OUT1}
