FSRCS    = $(wildcard *.f)

LIB      = libpythia.so

FF       = gfortran

FFLAGS   = -fPIC -g -O -fno-automatic -fno-second-underscore \
            -ffixed-line-length-none -funroll-loops -Wunused

SHFLAGS  = -fPIC $(MFLAG) -shared

FOBJS    = $(FSRCS:%.f=%.o)

all: $(LIB)

$(LIB): $(FOBJS)
	@echo "" ; echo "Creating shared library $(LIB) ..."
	$(FF) $(SHFLAGS) $(FOBJS) -o $@

.f.o :
	@echo Compiling $@
	@$(FF) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(FOBJS) $(LIB)
