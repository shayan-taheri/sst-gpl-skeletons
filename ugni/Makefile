TARGET := libugni.so
SRC := ugni.cc

#CXX :=    g++
#CC :=     gcc
CXX :=    libsst++
CC :=     libsstcc
CXXFLAGS := -fPIC
CPPFLAGS := -I. -I$(HOME)/Programs/install/sst-core/clang/include
LIBDIR :=  
PREFIX := $(shell $(CXX) --prefix)
LDFLAGS :=  -Wl,-rpath,$(PREFIX)/lib

OBJ := $(SRC:.cc=.o) 
OBJ := $(OBJ:.cpp=.o)
OBJ := $(OBJ:.c=.o)

.PHONY: clean install 

all: $(TARGET)

$(TARGET): $(OBJ) 
	$(CXX) -shared -o $@ $+ $(LDFLAGS) $(LIBS)  $(CXXFLAGS)

%.o: %.cc 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean: 
	rm -f $(TARGET) $(OBJ) 

install: $(TARGET) gni_pub.h pmi.h
	cp $< $(PREFIX)/lib
	cp gni_pub.h $(PREFIX)/include
	cp pmi.h $(PREFIX)/include

