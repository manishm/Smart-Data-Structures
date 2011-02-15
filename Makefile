CPPSRCS		= main.cpp cpp_framework.cpp FCBase.cpp LearningEngine.cpp rl_agent.cpp
REGRESS         = regress_intel64
TARGET		= main_intel64

CPP		= g++

INCFLAGS        = -Idata_structures -Iframework -Imonitor -Irl_engine -Irl_engine/Seldon-5.0.1
DEFFLAGS        = -DNDEBUG -DINTEL64 -D_REENTRANT
CPPFLAGS	= -MD -O3 -m64 $(INCFLAGS) $(DEFFLAGS) -lrt -pthread -llapack
LFLAGS		= -O3 -m64 -DNDEBUG -DINTEL64 -D_REENTRANT -lrt -pthread -llapack

OBJS		= $(CPPSRCS:.cpp=.o)
REGRESSOBJS     = $(subst main,regress,$(OBJS)) 
DEPS            = $(CPPSRCS:.cpp=.d) regress.d

all: $(TARGET)

main.o:
	$(CPP) $(CPPFLAGS) -c ./test/main.cpp

regress.o:
	$(CPP) $(CPPFLAGS) -c ./test/regress.cpp

cpp_framework.o:
	$(CPP) $(CPPFLAGS) -c ./framework/cpp_framework.cpp

FCBase.o:
	$(CPP) $(CPPFLAGS) -c ./data_structures/FCBase.cpp

LearningEngine.o:
	$(CPP) $(CPPFLAGS) -c ./data_structures/LearningEngine.cpp

rl_agent.o:
	$(CPP) $(CPPFLAGS) -c ./rl_engine/rl_agent.cpp

$(TARGET): $(OBJS)
	$(CPP) $(LFLAGS) $(OBJS) -o $(TARGET)

$(REGRESS): $(REGRESSOBJS)
	$(CPP) $(LFLAGS) $(REGRESSOBJS) -o $(REGRESS)

clean:
	rm -f $(DEPS) $(OBJS) $(REGRESSOBJS) $(TARGET) $(REGRESS)

-include $(DEPS)