TARGETS := libugni.la libibverbs.la
SRC := ugni.cc fakeverbs.cc

#CXX :=    g++
#CC :=     gcc
LIBTOOL = glibtool
CXX :=    libsst++
CC :=     libsstcc
CXXFLAGS := -fPIC
CPPFLAGS := -I. -I$(HOME)/Programs/install/sst-core/clang/include
LIBDIR :=  
PREFIX := $(shell $(CXX) --prefix)
LDFLAGS :=  -Wl,-rpath,$(PREFIX)/lib

.PHONY: clean install 

%.la: %.lo 
	$(LIBTOOL) --mode=link --tag=CXX $(CXX) -o $@ $+ $(LDFLAGS) $(LIBS) $(CXXFLAGS) -rpath $(PREFIX)/lib

%.lo: %.cc 
	$(LIBTOOL) --mode=compile --tag=CXX $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.cc: 

all: $(TARGETS)

clean: 
	rm -rf $(TARGETS) *.o *.lo .libs

install: $(TARGETS) gni_pub.h pmi.h
	$(LIBTOOL) --mode=install /usr/bin/install -c libugni.la $(PREFIX)/lib
	$(LIBTOOL) --mode=install /usr/bin/install -c libibverbs.la $(PREFIX)/lib
	-mkdir -p $(PREFIX)/include/infiniband
	cp verbs.h $(PREFIX)/include/infiniband
	cp gni_pub.h $(PREFIX)/include
	cp pmi.h $(PREFIX)/include

