
ROOTFLAGS   = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs)
QCDLOOPLIBS = ext/ql/libqcdloop.a ext/ff/libff.a  -lgfortran
LOOPTOOLS   = ext/LoopTools-2.12/build/libooptools.a
BOOSTLIBS   = /usr/lib64/libboost_program_options.so

DEBUG-FLAGS = -g -O0 -DDEBUG -pedantic
PERF1-FLAGS = -m64 -O1 -march=native
PERF2-FLAGS = -m64 -Ofast -flto -march=native -funroll-loops


DEP-PATH = .dep
LIB-PATH = lib
BIN-PATH = bin
CC-FLAGS = -Wall -std=c++11 -L/home/clemens/LHAPDF-6.1.4/include $(PERF1-FLAGS) -fopenmp
### -DDUMP_DIPOLE_PS 

SRC-PATH = src
CC = g++
CPP = /usr/lib/cpp


####################################################################
###################### local libraries #############################
####################################################################

_LIBS_COMMON   =  Global.o Lorentz.o HistArray.o Integrator.o Functions_Shared.o
_LIBS_PP_HX    =  $(_LIBS_COMMON) PhaseSpace_.o ScalarIntegrals_.o Integrands_pp_HX.o Functions_pp_HX.o

_LIBS_PP_TTX   = $(_LIBS_COMMON) HiggsModel_.o  PhaseSpace_.o  ScalarIntegrals_.o  Integrands_pp_ttX_.o  Functions_pp_ttX_V_.o  Functions_pp_ttX_ID_.o  Functions_pp_ttX_R_.o  Functions_pp_ttX_UID_.o 
_LIBS_PP_TTX_S = $(_LIBS_COMMON) HiggsModel_S.o PhaseSpace_S.o ScalarIntegrals_S.o Integrands_pp_ttX_S.o Functions_pp_ttX_V_S.o Functions_pp_ttX_ID_S.o Functions_pp_ttX_R_S.o Functions_pp_ttX_UID_S.o Functions_tDecay.o
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

exec:    $(BIN-PATH)/Test  $(BIN-PATH)/Integrate_pp_ttX $(BIN-PATH)/Integrate_pp_ttX_S  $(BIN-PATH)/Integrate_pp_ttX_withTdecay  $(BIN-PATH)/Integrate_pp_HX  $(BIN-PATH)/EvalSI



$(BIN-PATH)/Integrate_pp_ttX: $(LIB-PATH)/Integrate_pp_ttX.o $(LIBS_COMMON) $(LIBS_PP_TTX) 
	$(CC) -o $@  $^ -L/usr/local/lib $(QCDLOOPLIBS) $(LOOPTOOLS) -lgsl -lLHAPDF $(ROOTLIBS)

$(BIN-PATH)/Integrate_pp_ttX_S: $(LIB-PATH)/Integrate_pp_ttX_S.o $(LIBS_COMMON) $(LIBS_PP_TTX_S) 
	$(CC) -o $@  $^ -L/usr/local/lib $(QCDLOOPLIBS) $(LOOPTOOLS) -lgsl -lLHAPDF $(ROOTLIBS)

$(BIN-PATH)/Integrate_pp_HX: $(LIB-PATH)/Integrate_pp_HX.o $(LIBS_PP_HX)
	$(CC) -o $@  $^ -L/usr/local/lib $(QCDLOOPLIBS) $(LOOPTOOLS) -lgsl -lLHAPDF $(ROOTLIBS)

$(BIN-PATH)/Test: $(LIB-PATH)/Test.o $(LIBS_COMMON) $(LIBS_PP_TTX) 
	$(CC) -o $@  $^ -L/usr/local/lib $(QCDLOOPLIBS) $(LOOPTOOLS) -lgsl -lLHAPDF $(ROOTLIBS) $(BOOSTLIBS)

$(BIN-PATH)/Test_S: $(LIB-PATH)/Test_S.o $(LIBS_COMMON) $(LIBS_PP_TTX_S) 
	$(CC) -o $@  $^ -L/usr/local/lib $(QCDLOOPLIBS) $(LOOPTOOLS) -lgsl -lLHAPDF $(ROOTLIBS) $(BOOSTLIBS)

$(BIN-PATH)/EvalSI: $(LIB-PATH)/EvalSI.o 
	$(CC) -o $@  $^ $(LIB-PATH)/Global.o $(QCDLOOPLIBS) $(LOOPTOOLS) -lgsl

$(BIN-PATH)/TestVEGAS: $(LIB-PATH)/TestVEGAS.o $(LIB-PATH)/pVEGAS.o
	$(CC) -o $@  $^ -lgsl -lgomp  

####################################################################

### targets without reference to a file
.PHONY: clean all exec doc

doc:
	doxygen Doxyfile

clean:
	-rm -f $(BIN-PATH)/*
	-rm -f $(LIB-PATH)/*.o
	-rm -f $(LIB-PATH)/*.d
