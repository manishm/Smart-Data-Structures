CPPSRCS		= main.cpp cpp_framework.cpp LearningEngine.cpp rl_agent.cpp
TARGET		= main_intel64
REGRESS         = regress_intel64
TSP             = tsp_intel64

CPP		= g++

INCFLAGS        = -Idata_structures -Iframework -Imonitor -Irl_engine -Irl_engine/Seldon-5.0.1
TSPFLAGS        = 
DEFFLAGS        = -DINTEL64 -D_REENTRANT
CPPFLAGS	= -MD -O3 -m64 $(INCFLAGS) $(DEFFLAGS) -lrt -pthread -llapack
LFLAGS		= -O3 -m64 $(DEFFLAGS) -lrt -pthread -llapack

OBJS		= $(CPPSRCS:.cpp=.o)
REGRESSOBJS     = $(subst main,regress,$(OBJS)) 
TSPOBJS         = $(subst main,tsp,$(OBJS))
DEPS            = $(CPPSRCS:.cpp=.d) regress.d tsp.d

all: $(TARGET)

main.o:
	$(CPP) $(CPPFLAGS) -c ./test/main.cpp

regress.o:
	$(CPP) $(CPPFLAGS) -c ./test/regress.cpp

cpp_framework.o:
	$(CPP) $(CPPFLAGS) -c ./framework/cpp_framework.cpp

LearningEngine.o:
	$(CPP) $(CPPFLAGS) -c ./data_structures/LearningEngine.cpp

rl_agent.o:
	$(CPP) $(CPPFLAGS) -c ./rl_engine/rl_agent.cpp

$(TARGET): $(OBJS)
	$(CPP) $(LFLAGS) $(OBJS) -o $(TARGET)

$(REGRESS): $(REGRESSOBJS)
	$(CPP) $(LFLAGS) $(REGRESSOBJS) -o $(REGRESS)


tsp.o:
	$(CPP) $(CPPFLAGS) $(TSPFLAGS) -DOYAMAQUEUE -c ./test/tsp.cpp

$(TSP): $(TSPOBJS)
	rm -f tsp.o
	$(CPP) $(CPPFLAGS) $(TSPFLAGS) -DOYAMAQUEUE -c ./test/tsp.cpp
	$(CPP) $(LFLAGS) $(TSPFLAGS) -DOYAMAQUEUE $(TSPOBJS) -o $(TSP)
	mv $(TSP) $(TSP)_oyqueue
	rm -f tsp.o
	$(CPP) $(CPPFLAGS) $(TSPFLAGS) -DFCQUEUE -c ./test/tsp.cpp
	$(CPP) $(LFLAGS) $(TSPFLAGS) -DFCQUEUE $(TSPOBJS) -o $(TSP)
	mv $(TSP) $(TSP)_fcqueue
	rm -f tsp.o
	$(CPP) $(CPPFLAGS) $(TSPFLAGS) -DSMARTQUEUE -c ./test/tsp.cpp
	$(CPP) $(LFLAGS) $(TSPFLAGS) -DSMARTQUEUE $(TSPOBJS) -o $(TSP)
	mv $(TSP) $(TSP)_smartqueue
	rm -f tsp.o
	$(CPP) $(CPPFLAGS) $(TSPFLAGS) -DMUTEXQUEUE -c ./test/tsp.cpp
	$(CPP) $(LFLAGS) $(TSPFLAGS) -DMUTEXQUEUE $(TSPOBJS) -o $(TSP)
	mv $(TSP) $(TSP)_mutexqueue
	rm -f tsp.o
	$(CPP) $(CPPFLAGS) $(TSPFLAGS) -DBASKETSQUEUE -c ./test/tsp.cpp
	$(CPP) $(LFLAGS) $(TSPFLAGS) -DBASKETSQUEUE $(TSPOBJS) -o $(TSP)
	mv $(TSP) $(TSP)_basketsqueue
	rm -f tsp.o
	$(CPP) $(CPPFLAGS) $(TSPFLAGS) -DMSQUEUE -c ./test/tsp.cpp
	$(CPP) $(LFLAGS) $(TSPFLAGS) -DMSQUEUE $(TSPOBJS) -o $(TSP)
	mv $(TSP) $(TSP)_msqueue
	rm -f tsp.o
	$(CPP) $(CPPFLAGS) $(TSPFLAGS) -DFCPAIRHEAP -c ./test/tsp.cpp
	$(CPP) $(LFLAGS) $(TSPFLAGS) -DFCPAIRHEAP $(TSPOBJS) -o $(TSP)
	mv $(TSP) $(TSP)_fcpairheap
	rm -f tsp.o
	$(CPP) $(CPPFLAGS) $(TSPFLAGS) -DSMARTPAIRHEAP -c ./test/tsp.cpp
	$(CPP) $(LFLAGS) $(TSPFLAGS) -DSMARTPAIRHEAP $(TSPOBJS) -o $(TSP)
	mv $(TSP) $(TSP)_smartpairheap

clean:
	rm -f $(DEPS) $(OBJS) $(REGRESSOBJS) $(TSPOBJS) $(TARGET) $(REGRESS) $(TSP)
	rm -f $(TSP)_oyqueue $(TSP)_fcqueue $(TSP)_smartqueue $(TSP)_fcpairheap $(TSP)_smartpairheap $(TSP)_basketsqueue $(TSP)_mutexqueue $(TSP)_msqueue output* scratch* tmp

-include $(DEPS)