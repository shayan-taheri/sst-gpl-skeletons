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

SED_PREFIX := $(subst /,\/,$(PREFIX))

install: $(TARGETS) $(HEADERS)
	echo $(SED_PREFIX)
	$(LIBTOOL) --mode=install /usr/bin/install -c libugni.la $(PREFIX)/lib
	$(LIBTOOL) --mode=install /usr/bin/install -c libibverbs.la $(PREFIX)/lib
	$(LIBTOOL) --mode=install /usr/bin/install -c libpmi.la $(PREFIX)/lib
	-mkdir -p $(PREFIX)/include/infiniband
	-mkdir -p $(PREFIX)/pkgconfig
	sed -e 's/replInc/$(SED_PREFIX)\/include/g' cray-ugni.pc.in | sed -e 's/replLib/$(SED_PREFIX)\/lib/g' > cray-ugni.pc
	sed -e 's/replInc/$(SED_PREFIX)\/include/g' cray-gni-headers.pc.in > cray-gni-headers.pc
	cp cray-ugni.pc $(PREFIX)/pkgconfig
	cp cray-gni-headers.pc $(PREFIX)/pkgconfig
	cp verbs.h $(PREFIX)/include/infiniband
	cp $(HEADERS) $(PREFIX)/include

