TARGETS := libugni.la libibverbs.la libpmi.la

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

HEADERS = gni_pub.h pmi.h pmi_cray.h udreg_pub.h hugetlbfs.h

.PHONY: clean install 

%.la: %.lo 
	$(LIBTOOL) --mode=link --tag=CXX $(CXX) -o $@ $+ $(LDFLAGS) $(LIBS) $(CXXFLAGS) -rpath $(PREFIX)/lib

%.lo: %.cc $(HEADERS)
	$(LIBTOOL) --mode=compile --tag=CXX $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.lo: %.c 
	$(LIBTOOL) --mode=compile --tag=CC $(CC) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.cc: 

all: $(TARGETS)

clean: 
	rm -rf $(TARGETS) *.o *.lo .libs

install: $(TARGETS) $(HEADERS)
	$(LIBTOOL) --mode=install /usr/bin/install -c libugni.la $(PREFIX)/lib
	$(LIBTOOL) --mode=install /usr/bin/install -c libibverbs.la $(PREFIX)/lib
	$(LIBTOOL) --mode=install /usr/bin/install -c libpmi.la $(PREFIX)/lib
	-mkdir -p $(PREFIX)/include/infiniband
	cp verbs.h $(PREFIX)/include/infiniband
	cp $(HEADERS) $(PREFIX)/include

