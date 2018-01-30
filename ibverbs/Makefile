TARGETS := librdmacm.la

LIBTOOL = glibtool
CXX :=    libsst++
CC :=     libsstcc
CXXFLAGS := -fPIC -O0 -g
CPPFLAGS := -I. -I$(HOME)/install/sst-macro-devel/include -I/usr/lib/sys
LIBDIR := $(HOME)/install/sst-macro-devel/lib 
PREFIX := $(HOME)/install/sst-ibverbs
LDFLAGS :=  -Wl,-rpath,$(LIBDIR)

HEADERS = rdma_cma.h

.PHONY: clean install 

%.la: %.lo 
	$(LIBTOOL) --mode=link --tag=CXX $(CXX) -o $@ $+ $(LDFLAGS) $(LIBS) $(CXXFLAGS) -rpath $(LIBDIR) 

%.lo: %.cc $(HEADERS)
	$(LIBTOOL) --mode=compile --tag=CXX $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.lo: %.c 
	$(LIBTOOL) --mode=compile --tag=CC $(CC) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.cc: 

all: $(TARGETS)

clean: 
	rm -rf $(TARGETS) *.o *.lo .libs


install: $(TARGETS) $(HEADERS)
	mkdir -p $(PREFIX)
	mkdir -p $(PREFIX)/include
	cp -r infiniband $(PREFIX)/include
	mkdir -p $(PREFIX)/include/rdma
	cp rdma_cma.h $(PREFIX)/include/rdma
	mkdir -p $(PREFIX)/lib
	$(LIBTOOL) --mode=install /usr/bin/install -c librdmacm.la $(LIBDIR)
	$(LIBTOOL) --finish $(PREFIX)/lib
