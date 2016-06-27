OBJECT = simplemodel.o simutrans.o Zeipelmodel.o Laramodel.o 
PROJECT = simutranscc
SIMPLEMODEL_SOURCES = simplemodel.h simplemodel.cpp Zeipelmodel.cpp Zeipelmodel.h Laramodel.cpp Laramodel.h GDmodel.h 
#TRANSITMODEL_INTERFACE = transitmodel.i
SIMPLEMODEL_INTERFACE = simplemodel.i Laramodel.i Zeipelmodel.i
SOURCES = $(SIMPLEMODEL_SOURCES)
INTERFACES = $(SIMPLEMODEL_INTERFACE)
WRAPPERS = $(INTERFACES:.i=_wrap.cxx)
PROXIES = $(INTERFACES:.i=.py)

.PHONY: all
all: $(WRAPPERS) $(SOURCES)
	./setup.py build_ext -i

%_wrap.cxx: %.i %.h ./numpy.i
	swig -c++ -python $<


.PHONY: cc
cc: $(PROJECT)

$(PROJECT) :$(OBJECT)
	$(CXX) -pg -o $@ $^

#.PHONY: test
#test: $(PROJECT_DEBUG)
#
#$(PROJECT_DEBUG) :$(OBJECT_DEBUG)
#	$(CXX) -pg -o $@ $^


.PHONY: clean
clean:
	$(RM) $(PROJECT)
	$(RM) $(PROJECT_DEBUG)
	$(RM) *.so *.pyc *_wrap.h
	$(RM) -r *.dSYM
	$(RM) -r build
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
	$(RM) *.o
	$(RM) .depend
	$(RM) *.out
	$(RM) *~
