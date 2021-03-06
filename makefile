


ROOTFLAGS   = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs)
QCDLOOPLIBS = -lgfortran ext/ql/libqcdloop.a ext/ff/libff.a
LOOPTOOLS   = ext/LoopTools-2.12/build/libooptools.a
GSLLIBS     = -lgsl -lgslcblas
BOOSTLIBS   = -lboost_system -lboost_filesystem

EXT_LIBS    = $(ROOTLIBS) $(QCDLOOPLIBS) $(GSLLIBS) $(BOOSTLIBS) -lLHAPDF

DEBUG-FLAGS = -g -O0 -DDEBUG -pedantic
PERF1-FLAGS = -m64 -O1 -march=native
PERF2-FLAGS = -m64 -Ofast -flto -march=native -funroll-loops


#DEP-PATH = .dep
LIB-PATH = lib
BIN-PATH = bin
CC-FLAGS = -Wall -std=c++11 -L/home/clemens/LHAPDF-6.1.4/include

DEBUG ?= 0
ifeq ($(DEBUG), 1)
    CC-FLAGS += $(DEBUG-FLAGS)
else 
    CC-FLAGS += $(PERF2-FLAGS)
endif


SRC-PATH = src
CC = g++
CPP = /usr/lib/cpp



####################################################################
###################### local libraries #############################
####################################################################

_LIBS_COMMON   =  Global.o Lorentz.o HistArray.o Functions_Shared.o ScalarIntegrals.o FileBrowser.o Cuts.o
_LIBS_PP_HX    =  $(_LIBS_COMMON) HiggsModel_.o Integrator_.o PhaseSpace_.o Observables_.o Integrands_pp_HX.o Functions_pp_HX.o

_LIBS_PP_TTX   = $(_LIBS_COMMON) HiggsModel_.o  Integrator_.o  PhaseSpace_.o  Observables_.o  Integrands_pp_ttX_.o  Functions_pp_ttX_V_.o  Functions_pp_ttX_ID_.o  Functions_pp_ttX_R_.o  Functions_pp_ttX_UID_.o 
_LIBS_PP_TTX_S = $(_LIBS_COMMON) HiggsModel_S.o Integrator_S.o PhaseSpace_S.o Observables_S.o Integrands_pp_ttX_S.o Functions_pp_ttX_V_S.o Functions_pp_ttX_ID_S.o Functions_pp_ttX_R_S.o Functions_pp_ttX_UID_S.o Functions_tDecay_S.o
####################################################################
###################### object files ################################
####################################################################

LIBS_COMMON    = $(patsubst %,$(LIB-PATH)/%,$(_LIBS_COMMON)) 
LIBS_PP_HX     = $(patsubst %,$(LIB-PATH)/%,$(_LIBS_PP_HX))
LIBS_PP_TTX    = $(patsubst %,$(LIB-PATH)/%,$(_LIBS_PP_TTX))
LIBS_PP_TTX_S  = $(patsubst %,$(LIB-PATH)/%,$(_LIBS_PP_TTX_S)) 


$(LIB-PATH)/%_.o: $(SRC-PATH)/%.cpp
	$(CC) -c -o $@  $(ROOTFLAGS)  $(CC-FLAGS) $^

$(LIB-PATH)/%_S.o: $(SRC-PATH)/%.cpp
	$(CC) -c -o $@  $(ROOTFLAGS)  $(CC-FLAGS) -DWITH_T_SPIN $^

$(LIB-PATH)/%.o: $(SRC-PATH)/%.cpp 
	$(CC) -c -o $@  $(ROOTFLAGS)  $(CC-FLAGS) $^

####################################################################
###################### executables #################################
####################################################################

.PHONY: exec
exec:    $(BIN-PATH)/Test  $(BIN-PATH)/Integrate_pp_ttX $(BIN-PATH)/Integrate_pp_ttX_S  $(BIN-PATH)/Integrate_pp_ttX_withTdecay  $(BIN-PATH)/Integrate_pp_HX  $(BIN-PATH)/EvalSI



$(BIN-PATH)/Integrate_pp_ttX: $(LIB-PATH)/Integrate_pp_ttX.o $(LIBS_COMMON) $(LIBS_PP_TTX) 
	$(CC) -o $@  $^ -L/usr/local/lib $(EXT_LIBS)

$(BIN-PATH)/Integrate_pp_ttX_S: $(LIB-PATH)/Integrate_pp_ttX_S.o $(LIBS_COMMON) $(LIBS_PP_TTX_S) 
	$(CC) -o $@  $^ -L/usr/local/lib $(EXT_LIBS)

$(BIN-PATH)/Integrate_pp_HX: $(LIB-PATH)/Integrate_pp_HX.o $(LIBS_PP_HX)
	$(CC) -o $@  $^ -L/usr/local/lib $(EXT_LIBS)

$(BIN-PATH)/Integrate_LO: $(LIB-PATH)/Integrate_LO.o $(LIBS_COMMON) $(LIBS_PP_TTX) 
	$(CC) -o $@  $^ -L/usr/local/lib $(EXT_LIBS)

$(BIN-PATH)/Test: $(LIB-PATH)/Test.o $(LIBS_COMMON) $(LIBS_PP_TTX) 
	$(CC) -o $@  $^ -L/usr/local/lib $(EXT_LIBS)

$(BIN-PATH)/Test_S: $(LIB-PATH)/Test_S.o $(LIBS_COMMON) $(LIBS_PP_TTX_S) 
	$(CC) -o $@  $^ -L/usr/local/lib $(EXT_LIBS)

$(BIN-PATH)/EvalSI: $(LIB-PATH)/EvalSI.o 
	$(CC) -o $@  $^ $(LIB-PATH)/Global.o $(QCDLOOPLIBS) $(GSLLIBS) 

$(BIN-PATH)/TestVEGAS: $(LIB-PATH)/TestVEGAS.o $(LIB-PATH)/pVEGAS.o
	$(CC) -o $@  $^ -lgsl -lgomp  

$(BIN-PATH)/makePlots: $(LIB-PATH)/makePlots.o
	$(CC) -o $@  $^ -L/usr/local/lib -lboost_system -lboost_filesystem $(ROOTLIBS)

$(BIN-PATH)/readDat: $(LIB-PATH)/readDat.o $(LIB-PATH)/FileBrowser.o
	$(CC) -o $@  $^ -L/usr/local/lib -lboost_system -lboost_filesystem $(ROOTLIBS)
####################################################################

### targets without reference to a file
.PHONY: doc
doc:
	doxygen Doxyfile

.PHONY: clean
clean:
	-rm -f $(BIN-PATH)/*
	-rm -f $(LIB-PATH)/*.o
	-rm -f $(LIB-PATH)/*.d
