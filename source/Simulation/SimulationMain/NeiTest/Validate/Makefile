PREFLASH_PATH=/home/chris/Flash/PreFlash

include ${PREFLASH_PATH}/settings.make

CCFLAGS += -I ${PREFLASH_PATH}/include $(CCFLAGS_NOSHARED)
LDFLAGS += ${PREFLASH_PATH}/lib/libpreflash.a $(LDFLAGS_NOSHARED) -L/opt/mpich2/lib -lmpich -lpthread -lrt -lz


BIN_FILES = extract_data


all : $(BIN_FILES)

clean : clean-obj
	rm -f $(BIN_FILES)

clean-obj :
	rm -f *.o

extract_data : extract_data.o
	$(LD) extract_data.o -o $(BIN_FILES) $(LDFLAGS) $(LDFLAGS_HDF5)

extract_data.o : extract_data.cpp
	$(CXX) -c extract_data.cpp $(CCFLAGS) $(CCFLAGS_HDF5) -o extract_data.o
