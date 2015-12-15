PYDIR   = $(PWD)/pylib
pylib   = $(PYDIR)/libpythia.so

CCSRCS  = $(wildcard ./src/*.cc)
DEPS    = $(wildcard ./include/*.hh)

EXE        = pythia.exe

CXX        = g++

CCFLAGS    = -g -Wall -fPIC -Wno-deprecated
CCFLAGS    += `root-config --cflags`
CCFLAGS    += -I./include -I$(ROOTSYS)/include

PYLIB      = -L$(PYDIR) -lpythia
ROOTLIB    = `root-config --glibs`
MYLIBS     = -Wl,-rpath, -lm -lgfortran
LIBS       = $(ROOTLIB) $(PYLIB) $(MYLIBS)

CCOBJS    = $(CCSRCS:%.cc=%.o)

all: $(pylib) $(EXE)


$(EXE): $(pylib) $(CCOBJS)
	$(CXX) $(CCFLAGS) -o $@ $(CCOBJS) $(LIBS)

$(pylib):
	for dir in $(PYDIR); do make -C $$dir; done;

#.cc.o :
%.o: %.cc $(DEPS)
	@echo Compiling $@
	@$(CXX) $(CCFLAGS) -c $< -o $@

clean:
	rm -f $(EXE) $(CCOBJS)

cleanall:
	rm -f $(EXE) $(CCOBJS)
	for dir in $(PYDIR); do make -C $$dir clean; done;

print:
	@echo $(DEPS)
