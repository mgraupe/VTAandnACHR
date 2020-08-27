CXXFLAGS=-Wall -pg
OBJS = parameter.o compose.o nrutils.o motif.o nicotine.o

main: $(OBJS)
	g++ -fpermissive -pg -o nicotine $(OBJS)
	cp nicotine output/

.PHONY: clean
clean:
	rm -f *.o
	rm -f nicotine

.PHONY: data
data: 
	rm -f *.dat
