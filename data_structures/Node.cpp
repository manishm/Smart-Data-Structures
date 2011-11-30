#include "Node.h"

//g++ -DINTEL64 -I../framework Node.cpp ../framework/cpp_framework.cpp -o Node

using namespace std;

int main(int argc, char* argv[])
{
        const int threads = 14;
        const int els = 7;
        int total = els*threads*10;
        WorkQueueFreeList<FCIntPtr> list(threads, els);
        PtrNode<FCIntPtr>* res[threads][els]; 

        //test head wrapping
        cout << "starting head wrapping test" << endl;
        for(int t = 0; t < total; t+=els*threads) {
                for(int j = 0; j < threads; j++) {
		        for(int i = 0; i < els; i++) {  //make this 8 and expect a hang
                                res[j][i] = list.alloc(j);
                                res[j][i]->setkey(t+j*els+i);
		                res[j][i]->setvalue((long long int*) (t+j*els+i));
				cout << (long long int) (res[j][i]->getvalue()) << endl;
			}
			for(int i = 0; i < els; i++) {
			        list.dealloc(j, res[j][i]);
			}
		}
	}

        //sequential test of code for concurrency
        cout << "starting sequential test of code for concurrency" << endl;
        for(int i = 0; i < els; i++) 
                res[0][i] = list.alloc(0);
        for(int i = 0; i < els; i++) {
	        if ( 1 /* i != els-1*/ )   //set this and expect a hang
	                list.dealloc(1, res[0][i]);
	}
        for(int i = 0; i < els; i++) 
                res[0][i] = list.alloc(0);


        //parallel test




        return 0;
}


