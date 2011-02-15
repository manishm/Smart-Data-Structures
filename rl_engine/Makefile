
INC= -ISeldon-5.0.1
#LIB= -lcblas -llapack
LIB= -llapack
ALL= test_all test_gauss test_perm test_discrete test_binary

all: $(ALL)

test_all: test_all.cpp rl_agent.o rl_agent_c.h distance.o
	g++ $(INC) -o test_all test_all.cpp rl_agent.o distance.o $(LIB)

test_gauss: test_gauss.cpp rl_agent.o rl_agent_c.h 
	g++ $(INC) -o test_gauss test_gauss.cpp rl_agent.o $(LIB)

test_perm: test_perm.cpp rl_agent.o rl_agent_c.h distance.o
	g++ $(INC) -o test_perm test_perm.cpp rl_agent.o distance.o $(LIB)

test_discrete: test_discrete.cpp rl_agent.o rl_agent_c.h distance.o
	g++ $(INC) -o test_discrete test_discrete.cpp rl_agent.o distance.o $(LIB)

test_binary: test_binary.cpp rl_agent.o rl_agent_c.h distance.o
	g++ $(INC) -o test_binary test_binary.cpp rl_agent.o distance.o $(LIB)

rl_agent.o: rl_agent.cpp rl_agent.h rl_agent_c.h rl_util.cpp
	g++ -c $(INC) rl_agent.cpp

distance.o: distance.cpp distance.h
	g++ -c distance.cpp

clean:
	rm -f *\~ $(ALL) core* *.o
