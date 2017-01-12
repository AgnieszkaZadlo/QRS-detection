MKOCTFILE ?= mkoctfile
SED       ?= sed

OCT_FILES = \
  __fwht__.oct \
  __ultrwin__.oct \
  cl2bp.oct \
  medfilt1.oct \
  remez.oct \
  sosfilt.oct \
  upfirdn.oct

CL2BP_OBJECTS = cl2bp.o cl2bp_lib.o

OCT_SOURCES = $(patsubst %.oct,%.cc,$(OCT_FILES))

all: $(OCT_FILES)

%.o: %.cc
	$(MKOCTFILE) -Wall -c $< -o $@

%.oct: %.cc
	$(MKOCTFILE) -Wall $< -o $@

cl2bp.oct: $(CL2BP_OBJECTS)
	$(MKOCTFILE) -Wall $^ -o $@

$(CL2BP_OBJECTS): cl2bp_lib.h

PKG_ADD PKG_DEL: $(OCT_SOURCES)
	$(SED) -n -e 's/.*$@: \(.*\)/\1/p' $^ > $@-t
	mv $@-t $@

clean:
	-rm -f *.o *.oct PKG_*

.PHONY: all clean
