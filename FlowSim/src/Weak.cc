#include "../hipr-3.5/hi_pr.h"


#include "RelationLP.h"
#include <vector>
//delete the following line in the furture
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <queue>
#include <iterator>

using namespace std;


string int2string(const int& number)
{
/*	great function, but I don't know why it needs much more time then the other.
ostringstream oss;
	oss << number;
	return oss.str();*/
	char s[255];
	sprintf(s, "%d", number);
	return s;
}
	
Relation::Relation(Sparse* r){
	debug = false;
	DIVIDE_AND_CONQUER = true;
	FLOW = true;
	
	INPUT_PRECISION = 1000;
	
	rmsm = r;
	
	non_zeros = rmsm->non_zeros;
	cols = rmsm->cols;
	row_starts = rmsm->row_starts;
	
	n_states = rmsm->n;
	
	//simulated_by_s2 = new int[n_states]; //indicates whether the state is simulated by s2.
	//simulating_s1 = new int[n_states]; //indicates whether the state simulates s1.
	suc_s1 = new int[n_states];  //indicates whether suc_s1[i] is a successor of the state s1. If yes the index is saved
	suc_s2 = new int[n_states];  //indicates whether suc_s2[i] is a successor of the state s2. If yes the index is saved
	
	relation_forward_index = new int[n_states+1];
	relation_backward_index = new int[n_states+1];
	
	in_rel = new int*[n_states];
	for(int i=0; i<n_states; i++)	in_rel[i] = new int[n_states];
	
	size_of_relation = 0;	
	relation_forward = NULL;
	relation_backward = NULL;
	
	InitialRelations(3);
}

	
Relation:: ~Relation(){
	if(relation_forward_index != NULL)
		delete [] relation_forward_index;
	if(relation_forward != NULL)
		delete [] relation_forward;
	
	if(relation_backward_index != NULL)
		delete [] relation_backward_index;
	if(relation_backward != NULL)
		delete [] relation_backward;

	for(int i=0; i<n_states; i++)
		if(in_rel[i] != NULL)
			delete [] in_rel[i];
	if(in_rel != NULL)
		delete [] in_rel;	
	if(simulated_by_s2 != NULL)
		delete [] simulated_by_s2;
	if(simulating_s1 != NULL)
		delete [] simulating_s1;
	if(suc_s1 != NULL)
		delete [] suc_s1;
	if(suc_s2 != NULL)
		delete [] suc_s2;	
	
	if(relation != NULL)
		delete [] relation;
}
	
	
//the parameter num_par is the number of partitions of the states space.
void
Relation::InitialRelations(int num_par){
	//the relation R.
	int par_size = (int)n_states/num_par;
	size_of_relation = par_size * par_size * (num_par - 1)
		+ (n_states - (num_par - 1) * par_size) * (n_states - (num_par -1 ) * par_size);
	
	
	//the following line should be changed to Array class
	relation = new Pair[size_of_relation];
	
	if(relation == NULL){
		if(debug) printf("Error: assignment memory, take measure! Function: InitialRelations, file Relation.cc");
	}
	
	int index=0;
	if(debug) printf("number of states =%d, number of partition = %d, initial number of transitions=%d\n",n_states, num_par, size_of_relation);
	//init the relation, for partition 1,2,3, ..., num_par-1
	for(int np=1; np<num_par; np++){
		for(int i = (np-1) * par_size; i< np * par_size; i++){
			for(int j = (np-1) * par_size; j < np * par_size; j++){
				Pair p;      p.x=i;
				p.y=j;
				relation[index] = p;
				index++;
			}
		}
	}
	
	//init the relation, for partition num_par, i.e., the last partition
	for(int i=(num_par -1 ) * par_size; i<n_states; i++){
		for(int j=(num_par -1) * par_size; j<n_states; j++){
			Pair p;
			p.x=i;      p.y=j;
			relation[index]=p;
			index++;
		}
	}
	
}
	
void 
Relation::ConstructRelationStructure(){
	for(int i=0; i<n_states; i++)
		for(int j=0; j<n_states; j++)
			in_rel[i][j] = 0;
	
	if(relation_forward != NULL)
		delete [] relation_forward;
	if(relation_backward != NULL)
		delete [] relation_backward;
	relation_forward = new int[size_of_relation];
	relation_backward = new int[size_of_relation];
	
	for(int i=0; i<size_of_relation; i++){ //initialize relation_forward, relation_backward and calculate in_rel.
		relation_forward[i]=0;
		relation_backward[i]=0;
		in_rel[relation[i].x][relation[i].y] = 1; // for more detail see the header file.
	}
	
	for(int i=0; i<n_states+1; i++){
		relation_forward_index[i]=0;
		relation_backward_index[i]=0;
	}
	
	/* Build the relation forward and backward index
	1 since relation_forward_index[0] = 0, and ...
	*/
	for(int i=0; i<size_of_relation; i++){
		relation_forward_index[relation[i].x+1]++;
		relation_backward_index[relation[i].y+1]++;
	}
	for(int i=0; i<n_states; i++){
		relation_forward_index[i+1] += relation_forward_index[i];
		relation_backward_index[i+1] += relation_backward_index[i];
	}	
	
	//Build the forward successors.
	int *help_index_forward = new int[n_states];
	int *help_index_backward = new int[n_states];
	for(int i=0; i<n_states; i++){
		help_index_forward[i]=0;
		help_index_backward[i]=0;
	}
	
	for(int i=0; i<size_of_relation; i++){
		int k1 = relation_forward_index[relation[i].x] 
			+ help_index_forward[relation[i].x];
		relation_forward[k1] = relation[i].y;
		help_index_forward[relation[i].x]++;
		
		int k2 = relation_backward_index[relation[i].y]
			+ help_index_backward[relation[i].y];
		relation_backward[k2] = relation[i].x;
		help_index_backward[relation[i].y]++;
	}
	
	
	//   for(int i=0; i<n_states; i++){
	//     int l= relation_forward_index[i];
	//     int h = relation_forward_index[i+1];
	
	//    if(debug) printf("forward testl=%d,h=%d\n",l,h);
	//     for(int j=l; j<h; j++){
	//      if(debug) printf("****states %d, successor %d\n",i, relation_forward[j]);
	//     }
	//   }
	
	//   for(int i=0; i<n_states; i++){
	//     int l= relation_backward_index[i];
	//     int h = relation_backward_index[i+1];
	
	//    if(debug) printf("backward testl=%d,h=%d\n",l,h);
	//     for(int j=l; j<h; j++){
	//      if(debug) printf("****states %d, predecessor %d\n",i, relation_backward[j]);
	//     }
	//   }
	
	delete [] help_index_forward;
	delete [] help_index_backward;
}
	
/**
* Corresponds the SimulationRelation
*/
void
Relation::SimulationRelation(){
	//index of the new relations
	std::vector<int> nr;
	
	int iterations = 1;
	
	while(true){
		//sort the relation. (only for the purpose of accelerations) later ...
		//SortRelation();
		
		//everytime we must calculate it, because the relation changes as the program runs.
		ConstructRelationStructure();
		
		if(debug) {
			printf("\n\nStart of iteration %d, the relation size is %d\n The relation: ", iterations, size_of_relation);
			for(int i=0; i<size_of_relation; i++)  printf("(%d,%d) ", relation[i].x, relation[i].y); 
			printf("\n");
		}
		
		for(int i=0; i<size_of_relation; i++){
			if(debug) printf("[PAIR (%d,%d)]: ", relation[i].x, relation[i].y);
			if(WeakSimulation(relation[i].x, relation[i].y)) //the most important line
				{nr.push_back(i);}
		}
		
		if(debug) printf("End of iteration %d, the new relation size is %d\n", iterations, nr.size());
		
		iterations ++; 
			
		if(nr.size() == size_of_relation) 	break;
		
		size_of_relation = nr.size();
		Pair* new_relation = new Pair[size_of_relation];
		for(int i=0; i<size_of_relation; i++){
			new_relation[i] = relation[nr[i]];
			nr.clear();
		}
	
		delete [] relation;
		relation = new_relation;
	}//end while loop
	
	
	FILE* f = fopen("result.txt", "a");
	fprintf(f, "Simplex re. size: %d\t\t", size_of_relation);

	// divide and conque for lp.	
	//fprintf(f, "LP relation size: %d\t\t", size_of_relation);
	
	fclose(f);
}
	
	
//Parameter pair denotes the index of relation, for which we should check
bool 
Relation::WeakSimulation(int s1, int s2){		
	//a state can weak simulates itself
	if(s1 == s2) {
		if(debug) printf("It can simulates itsself, returns true\n");
			return true;
	}	
	
	//s2 weak simulates s1 if s1 is absorbing. 
	if(row_starts[s1] == row_starts[s1+1]-1) {
		//if the only successor of s1 is itself, s1 is a absorbing state.
		if( s1 == (cols[row_starts[s1]]) ){
			if(debug) printf("%d is absorbing, returns true\n", s1);
			return true;
		}
	}
	
	
	//Line 1, pure stutter step for state s1. We test whether post(s1) \subseteq R^{-1}[s2]
	bool pure_s1 = true;
	for(int i=row_starts[s1]; i<row_starts[s1+1]; i++){
		if( !( in_rel[cols[i]][s2] ) ){
			pure_s1 = false;
			break;
		}
	}
	if(pure_s1){
		if(debug) printf("Pure stutter step for state %d\n", s1);
			return true;
	}//end of the case pure stutter step s1.
	
	// to do, put it later.....
	// Line 2--6. This case consumes much more time then others. Since the reachable states of one states can involve a very large state space. So treate this later if possible ...
	bool pure_s2 = true;
	for(int i=row_starts[s2]; i<row_starts[s2+1]; i++){
		if( !(in_rel[s1][cols[i]]) ){
			pure_s2 = false;
			break;
		}
	}
	if(pure_s2){
		if(debug) printf("pure stutter step of %d, ", s2);
		int size_of_reachable_states = 0;
		int *reachable = ReachableStates(s2, size_of_reachable_states);
			
		//line 3 of WeakSimulation. If we can find a reachable state from s2, which does not simulate s1, the part returns true. 
		for(int i=0; i<size_of_reachable_states; i++){
			if( !in_rel[s1][reachable[i]] ) {
				if(debug) printf("not all reachable states can simulate s1, returns true (line 3)\n");
				return true;
			}
		} //end line 3
		
		//Now line 4 and 5 together. First we assume that line 5 returns true. If exists one state ui, and there does not  exists s \in reach(s2), such that s\in R[u1], it will be set to false.
		bool everything_ok = true;
		for (int ui_index = row_starts[s1]; ui_index < row_starts[s1+1]; ui_index++) {
			int ui = cols[ui_index];
			if( in_rel[ui][s2] ) continue;
			
			//exists s \in reach(s2), such that s\in R[u1]
			bool exists_reachable_simulate_u1 = false;
			for(int re=0; re<size_of_reachable_states; re++){
				if( in_rel[ui][reachable[re]] ){
					exists_reachable_simulate_u1 = true;
					break;
				}
			}
			if(!exists_reachable_simulate_u1){
				everything_ok = false;
				break;
			}//end of the for loop
		}//end of the main for loop
		
		delete [] reachable;
		if(debug) cout<<"pure stutter steps of s2, returns "<<everything_ok<<endl;
		
		return everything_ok;//end of lines 4,5
	}//end line 2--6, end of the pure stutter s2 case.
		
	//if(DIVIDE_AND_CONQUER){ return DivideAndConquer(s1, s2); }
	return DirectSimplex(s1, s2);
}
	
bool
Relation:: DirectSimplex(int s1, int s2){
	//for the headers
	int left = 0, right = 0;
	
	
	FILE *file=fopen("input_simplex.txt","w");  
	string partition_arcs = ConstructTransitionArcs(s1, s2);
	fprintf(file, "%s", partition_arcs.c_str() );
	fprintf(file, "* end of main transitions, i.e., without source and sink.\n");
	
	if(debug) printf(partition_arcs.c_str());
	
	int l_s1 = row_starts[s1],
	h_s1 = row_starts[s1+1];
	int l_s2 = row_starts[s2],
	h_s2 = row_starts[s2+1];     
	left = h_s1 - l_s1;
	right = h_s2 - l_s2;
	
	
	for(int i=l_s1; i<h_s1; i++){
		float cap = non_zeros[i];
		int in_other = 0;
		if( !SimulatedByS2(s2, cols[i]) ) in_other = 1;
		fprintf(file, "s * %d %f %d\n", cols[i], cap, in_other);
	}// end of the for loop
	
	for(int i=l_s2; i<h_s2; i++){
		float cap = non_zeros[i];
		int in_other = 0;
		if( !SimulatesS1(s1, cols[i]) ) in_other = 1;
		fprintf(file, "t %d * %f %d\n", cols[i], cap, in_other);          
	}// end of the for loop


	fprintf(file, "# end of the whole transitions.\n");     
	fprintf(file, "c states s1,s2 (%d, %d), run simplex: number of successors s1: %d; of s2: %d\n", s1, s2, left, right);
	//we need later calculate the number of transitions
	fprintf(file, "p max %d\n", left+right+2);
	//   fprintf(file, "n 1 s\n");
	//   fprintf(file, "n %d t\n", left+right+2);


	fclose(file);


	LpSolver* lp = new LpSolver(debug);      
	bool is_feasible = lp->DirectSimplexFeasible("input_simplex.txt", left, right);
	delete lp;
	
	
	if(is_feasible){   if(debug) printf("returns true(^_*)\n");   return true;  }
	else           {   if(debug) printf("returns false\n");       return false; }
}
	

bool 
Relation::DivideAndConquer(int s1, int s2){
	//construct the sets A1, ... An. where n is the number of partitions.
	int number_partitions = 0;
	set<int> state_s;
	merge(cols+row_starts[s1], cols+row_starts[s1+1], cols+row_starts[s2], cols+row_starts[s2+1], inserter(state_s,state_s.begin()));
	int size_part = state_s.size();
	int* par_array = new int[size_part];
	
	//We define it large enough. The reason that we must initilaize the pointer is that if we do this in the function ConstructStatePartition, we will get segmentation fault later. It seems that you must initialize the pointer, then give it as parameter to a function.
	int* par_array_index = new int[size_part];
	
	//note actually the function ConstructStatePartition always return true, so pls update this line later
	ConstructStatePartition(s1, s2, size_part, number_partitions, par_array, par_array_index);	
	state_s.clear();
	
	// only for debugging  ************************** start
	if(debug) {
		printf("n=%d,", number_partitions);
		for(int i=0; i<number_partitions; i++){
			printf("A%d=(",i);
			for(int j=par_array_index[i]; j<par_array_index[i+1]; j++)
				printf("%d;",par_array[j]);
			printf("),\n");
		}
	}
	// only for debugging  ************************** end	
	
	
	bool return_value;
	//The case n==1, line 17. We have either maxflow algorithm or lp algorithm
	if(number_partitions == 1) {    return_value = ReducedWS(s1, s2);     }
	else{//
		int p1=0, p2=0;
		//in function CheckMaxFlow we either call maxflow or lp
		return_value = CheckMaxFlow(s1, s2, number_partitions, par_array, par_array_index, p1, p2 );
		if(return_value){
			if(FLOW)				return_value = CheckFeasibleFlow(s1, s2, par_array, par_array_index, p1, p2);
			else  				return_value = CheckFeasibleFlowViaLpSolver(s1, s2, par_array, par_array_index, p1, p2);
		}
	}
	
	delete [] par_array;
	delete [] par_array_index;
	
	return return_value;
}
	
	
/** constructs the A1, A2, ... An.
*  

For technical reasons, the states in from par_array[par_array_index[0]] to par_array[par_array_index[1]] are those states, which can reach a state, which can either simulate s1, or be simulated by s2. Note that in the paper this special partition is called An. 

In this function I did not construct the set H as described in the paper. 


the function returns false, if any partition except the first one contains only one element, which indicates that the s2 does not simulate s1. Otherwise, the algorithm returns true. 

this is not the case, for exampe if s1->s and s2->s and s->s, in this case the algorithm should return true. So currently I just comments the corresponding line,

***: later I should update this thing 


(s1,s2) is the states pair
number_partitions: the number of partitions will be saved in this variable and can be used in the callinf function
par_array: The set of states (successor states of s1 and s2)
par_array_index: indicates the start and end of every partition
successor_of_par_array: indicates which successor states it is, 1 indicates of s1, 2 indicates of s2.

for example:

s1-->1,2,3    s2--> 5,6,7,8. Arcs: 1-->6, 1-->5, 2-->7, 3-->7, 3-->8.
number_partitions: 2
par_array: 			1 5 6 2 3 7 8
par_array_index:		0     3      7
*/

//update this function!! 1. the return value.  2. The efficiency. This fuctions consumes 20% of the whole running time. Please check and find a efficient method. Try the sets structures.
void
Relation:: ConstructStatePartition(int s1, int s2,
		int size_part,
		int& number_partitions, int* par_array, 
		int* par_array_index)
{  
	int* explored = new int[n_states]; //indicates whether the state (in state_s) is already explored.
	queue<int> qstates;	// used for the BFS algoriths.
	int processing = 0; // count the number of processing states. Increased after qstates.push().
	for(int i=0; i<n_states; i++){//set the values back to 0.
		//simulated_by_s2[i] = 0;
		//simulating_s1[i] = 0;
		suc_s1[i] = 0;
		suc_s2[i] = 0;
		explored[i] = 0;
	}
	//for (int i = relation_backward_index[s2]; i< relation_backward_index[s2+1];i++) simulated_by_s2[relation_backward[i]] = 1;
	//for (int i = relation_forward_index[s1]; i < relation_forward_index[s1+1]; i++) simulating_s1[relation_forward[i]] = 1;
	for (int i = row_starts[s1]; i < row_starts[s1+1]; i++)  { //insert V1 into the qstates
		suc_s1[cols[i]] = i;
		if(in_rel[cols[i]][s2]){
		//if(simulated_by_s2[cols[i]]){
			qstates.push(cols[i]);
			explored[cols[i]] = 1;
			processing ++;
		}
	}
	for (int i = row_starts[s2]; i < row_starts[s2+1]; i++)  { //insert V2, if not already inserted, into the qstates
		suc_s2[cols[i]] = i;
		//if( simulating_s1[cols[i]] & (!explored[cols[i]]) ){
		if( in_rel[s1][cols[i]] & (!explored[cols[i]]) ){
			qstates.push(cols[i]);
			explored[cols[i]] = 1;
			processing ++;
		}
	}


	par_array_index[0] = 0;
	number_partitions = 0;
	int pi = 0; //is the index for par_array. Increased if a new value is inserted into par_array
	//if qstates is empty, the set An is empty. The special partition is saved at the beginning. In the paper is the last partition.
	if(qstates.empty()) {
		number_partitions ++;
		par_array_index[number_partitions] = 0;
	}

	bool finished = false;
	while( !finished ){
		if(!qstates.empty()) number_partitions++;
		if(debug){
			cout<<"******* before while,\t par_array=(";
			for(int i=0; i<pi; i++){   cout<<";"<<par_array[i];    }
			cout<<"),"<<endl;
		}
		while(!qstates.empty()){
			
			int last = qstates.front();
			qstates.pop();
			par_array[pi] = last;
			pi++;
			if(suc_s1[last]){
				for(int i=relation_forward_index[last]; i<relation_forward_index[last+1]; i++){
					int s = relation_forward[i];
					if( suc_s2[s] && (!explored[s]) ){
						qstates.push(s);
						explored[s] = 1;
						processing ++;
					}
				}
			}
			if(suc_s2[last]){
				for(int i=relation_backward_index[last]; i<relation_backward_index[last+1]; i++){
					int s = relation_backward[i];
					if( suc_s1[s] && (!explored[s]) ){
						qstates.push(s);
						explored[s] = 1;
						processing ++;
					}
				}
			}
		}//end of the first loop
		if(debug){
			cout<<"******* after while,\t par_array=(";
			for(int i=0; i<pi; i++){   cout<<";"<<par_array[i];    }
			cout<<"),"<<endl;
		}

		//now we need still to set the par_array_index
		par_array_index[number_partitions] = pi;

		for (int i = row_starts[s1]; i < row_starts[s1+1]; i++)  { //insert state into the qstates
			if(!(explored[cols[i]])){
				qstates.push(cols[i]);
				explored[cols[i]] = 1;
				processing ++;
				break;
			}
		}
		if(qstates.empty()){
			for (int i = row_starts[s2]; i < row_starts[s2+1]; i++)  { //insert state into the qstates
				if(!(explored[cols[i]])){
					qstates.push(cols[i]);
					explored[cols[i]] = 1;
					processing ++;
					break;
				}
			}
		}
		
		if(qstates.empty()) finished = true;
	}//end of while
	
	delete [] explored;
	
	if(!qstates.empty()){
		cout<<"irgendwas stimmt nicht"<<endl;
		exit(0);
	}
}
		
//return true if state simulates s1, i.e., (s1,state)\in R for the current relation
bool
Relation:: SimulatesS1(int s1, int state){
	//we save successor states of post(s2), which can simulate s1
	int l_frel = relation_forward_index[s1];
	int h_frel = relation_forward_index[s1+1];
	
	for (int i = l_frel; i < h_frel; i++){
		if(state == relation_forward[i]){
			return true;
		}
	}
	
	return false;
}
	
	
//Returns true if s2 can simulates state, i.e., (state, s2)\in R for the current relation
bool
Relation:: SimulatedByS2(int s2, int state){
	
	int l_brel = relation_backward_index[s2];
	int h_brel = relation_backward_index[s2+1];
	
	for (int i = l_brel; i < h_brel; i++){
		if(state == relation_backward[i]){
			return true;
		}
	}
	
	return false;
}
	
	
bool
Relation:: IsSuccessorOf(int suc, int state){
	int l_tra = row_starts[state];
	int h_tra = row_starts[state+1];
	bool is_successor = false;
	for (int i = l_tra; i < h_tra; i++){
		if(cols[i] == suc){
			is_successor = true;
			break;
		}
	}
	return is_successor;
}
	
bool
Relation:: IsSuccessorOf(int suc, int state, int& successor_index){
	int l_tra = row_starts[state];
	int h_tra = row_starts[state+1];
	bool is_successor = false;
	for (int i = l_tra; i < h_tra; i++){
		if(cols[i] == suc){
			is_successor = true;
			successor_index = i;
			break;
		}
	}
	return is_successor;
}
	
	
/*
the parameter index indicates which partition
*/
string
Relation:: ConstructTransitionArcs(int s1, int s2,
				int* par_array, 
				int* par_array_index, 
				int index){
	int l = par_array_index[index];
	int h = par_array_index[index+1];  
	
	string transition_arcs = "";
	
	for(int i=l; i<h; i++){
		if(IsSuccessorOf(par_array[i], s1) ){
			for(int j=l; j<h; j++){
				if(IsSuccessorOf(par_array[j], s2)){
					if( SimulatesS1(par_array[i], par_array[j]) ){                          
					char s[20];
					sprintf(s, "a %d %d\n", par_array[i], par_array[j]);
					string c_plus_s(s);
					(transition_arcs) += c_plus_s;
					}
				}
			}
		}
	
	}//end of for
	
	return transition_arcs;
}
	
/*
this is the transition arcs for all successors of s1 and s2.
*/
string
Relation:: ConstructTransitionArcs(int s1, int s2){
	int l_s1 = row_starts[s1];
	int h_s1 = row_starts[s1+1];
	int l_s2 = row_starts[s2];
	int h_s2 = row_starts[s2+1];     
	
	string transition_arcs = "";
	
	for(int i=l_s1; i<h_s1; i++){
		for(int j=l_s2; j<h_s2; j++){
			if( SimulatesS1(cols[i], cols[j]) ){                          
				char s[20];
				sprintf(s, "a %d %d\n", cols[i], cols[j]);
				string c_plus_s(s);
				(transition_arcs) += c_plus_s;
			}
		}
	}//end of for
	
	return transition_arcs;
}
	
	
	
	
/*zhang  the following code must updated ... */

bool
Relation:: CheckMaxFlow(int s1, int s2, int number_partitions, 
				int* par_array, int* par_array_index, 
				int& p1, int& p2)
{
	//dealing with the partitions except the frist which will be dealt with later. Recall that I put all s1 s2 related states in the first partition. Be careful: How to differenciate 1 (one) and l (letter)
	for(int i=1; i<number_partitions; i++){// number i denotes the i-th partition
		if(debug) printf("checking partition A%d, ", i);
		int l = par_array_index[i];
		int h = par_array_index[i+1];  
		int prob1 = 0, prob2 = 0; // prob1 is the probabilities from s1 to post_i(s1)
		
		/* The states is matched as follows. 
			left denotes the size of post_i(s1), 
			and right denotes post_i(s2). 
			1--> source
			2--> sink
			3 ... 3+left-1 --> left_states[0], ... left_states[left-1]
			3+left, ... 3+left+right --> right_states[0], ... right_states[right-1]
			Example left: 1,2,3, right 1,2,3,4,
			1 (source), 2(sink), 3(1 left), 4(2 left), 5(3 left), 6(1 right), 7(2 right), 8(3 right), 9(4 right).
		*/
		int* left_states = new int[h-l];  //The state maching for post_i(s1).
		int* right_states = new int[h-l]; //The state maching for post_i(s1).
		int* left_cap = new int[h-l];     // save the corresponding capacities to left_states
		int* right_cap = new int[h-l];    // save the corresponding capacities to right_states
		
		int left=0, right=0;
		//by the construction, states in Partition have the property: either P(s1,state)>0 or P(s2,state)>0    
		for(int j=l; j<h; j++){
			int suc_s1_index = suc_s1[par_array[j]]; 
			if(suc_s1_index){     
				int cap = (int) (non_zeros[suc_s1_index] * INPUT_PRECISION);        
				left_states[left] = par_array[j];
				left_cap[left] = cap;
				left++;
				prob1 += cap;
			}      
			int suc_s2_index = suc_s2[par_array[j]];         
			if(suc_s2_index){       
				int cap = (int) (non_zeros[suc_s2_index] * INPUT_PRECISION);
				right_states[right] = par_array[j];
				right_cap[right] = cap;
				right++;
				prob2 += cap;
			}
		}// end of the for loop
		
		if((prob1 == 0) && (prob2 == 0)){
			if(debug) printf("p1=p2=0, error not possbile, please check where are wrong, Relation.cc->CheckMaxFlow()\n\n\n");
			exit(0);
		}else if(prob1 == 0 || prob2==0){
			if(debug) printf("prob1 (%d) or prob2 (%d) equals 0, please check. returns false \n", prob1, prob2);
			return false;
		}
		
		int MAX= 66;
		//(h-l)^2 is the number of maximal transitions from post(s1) to post(s2). (h-l)*2 is the transitions from source or to sink. 5 is the header lines.
		int rows = (h-l)*(h-l) + (h-l)*2 + 5;
		char** inp_input = new char*[rows];
		for(int m=0; m<rows; m++) 
			inp_input[m] = new char[MAX];
		int product = prob1 * prob2;		
		char* arc = new char[MAX];
		int index = 5;
		
		//This time we insert the arcs from the sources and to the sink t in the problem description. Note I use the vector structure here. Therefore I can use the same arc array. The variable arc is pushed in the vector, everytime before it is changed. Otherwise it is probalimatic, since arc is a pointer, and it points everytime to a new thing.
		for(int j=0; j<left; j++){
			sprintf(arc, "a 1 %d %d", j+3, left_cap[j] * prob2 );
			strcpy(inp_input[index], arc);
			index++;
			for(int k=0; k<right; k++){
				//if( SimulatesS1(left_states[j], right_states[k]) ){                          
				if( in_rel[left_states[j]][right_states[k]] ){
					sprintf(arc, "a %d %d %d", j+3, k+3+left, product);
					strcpy(inp_input[index], arc);
					index++;
				}
			}
		}
		
		for(int j=0; j<right; j++){
			sprintf(arc, "a %d 2 %d", j+3+left, right_cap[j] * prob1);
			strcpy(inp_input[index], arc);
			index++;
		}// end of the arcs
		
		bool is_maxflow = false;
		if(FLOW){  
			if(debug) printf("^_*: c states s1,s2 (%d,%d), %d\n", s1 ,s2, i);
			strcpy(inp_input[0], "c set debug to see the states after ^_*:");
			if(debug) printf("^_*: c prob_1=%d , prob_2=%d", prob1, prob2);
			strcpy(inp_input[1], "c set debug to see the probss after ^_*:");
			sprintf(arc, "p max %d %d", left + right + 2, index-5);//5 is the number of head lines.
			strcpy(inp_input[2], arc);	
			strcpy(inp_input[3], "n 1 s");
			strcpy(inp_input[4], "n 2 t");
			
			if(debug){
				for(int i=0; i<index; i++){
					cout<<inp_input[i]<<endl;
				}
			}
			
			//if(product == GetFlow(inp_input, index)) 
				is_maxflow = true;
		}else{/* update later the part for lp solver.
			LpSolver* lp = new LpSolver(debug);
			is_maxflow = lp->CalculateMaxFlow("input_flow.txt", prob1, prob2, left, right);
			delete lp;*/
		}
		
		delete [] arc;		
		for(int m=0; m<rows; m++) { delete [] inp_input[m]; }
		delete [] inp_input;
		delete [] left_states;
		delete [] right_states;
		delete [] left_cap;
		delete [] right_cap;

		if(!is_maxflow){
			if(debug) printf("returns false\n");
			return false;
		}else{
			p1 = prob1;
			p2 = prob2;
		}
	}//end of the special partitions
	
	if(debug) printf("CheckMaxFlow returns true. ");
	return true;
}



bool
Relation:: CheckFeasibleFlow(int s1, int s2,  
				int* par_array, 
				int* par_array_index, 
				int prob1, 
				int prob2){

	//Now with help of the value gamma, we can deal with the first partition.
	int l = par_array_index[0];
	int h = par_array_index[1];  
	
	if(h-l <= 1) {
		if(debug) printf(" A1 contains at most one element. returns true\n");
		return true;
	}  
	
	/* The states is matched as follows. 
		left denotes the size of post_i(s1), 
		and right denotes post_i(s2). 
		1--> SOURCE
		2--> SINK
		3--> source
		4--> sink
		5 ... 5+left-1 --> left_states[0], ... left_states[left-1]
		5+left, ... 5+left+right --> right_states[0], ... right_states[right-1]
		Example left: 1,2,3, right 1,2,3,4,
		1 (SOURCE) 2(SINK) 3(source), 4(sink), 5(1 left), 6(2 left), 7(3 left), 8(1 right), 9(2 right), 10(3 right), 11(4 right).
	*/
	int* left_states = new int[h-l];  //The state maching for post_i(s1).
	int* right_states = new int[h-l]; //The state maching for post_i(s1).
	int* left_cap = new int[h-l];     // save the corresponding capacities to left_states
	int* right_cap = new int[h-l];    // save the corresponding capacities to right_states
	
	int left=0, right=0;
	int prob_other1 = 0, prob_other2 = 0;
	//by the construction, states in Partition have the property: either P(s1,state)>0 or P(s2,state)>0    
	for(int j=l; j<h; j++){
		int suc_s1_index = suc_s1[par_array[j]]; 
		if(suc_s1_index){     
			int cap = (int) (non_zeros[suc_s1_index] * INPUT_PRECISION);        
			left_states[left] = par_array[j];
			left_cap[left] = cap;
			left++;
			if( !(in_rel[par_array[j]][s2]) ) prob_other1 += cap; //if it belongs to MU1
		}      
		int suc_s2_index = suc_s2[par_array[j]];         
		if(suc_s2_index){       
			int cap = (int) (non_zeros[suc_s2_index] * INPUT_PRECISION);
			right_states[right] = par_array[j];
			right_cap[right] = cap;
			right++;
			if( !(in_rel[s1][par_array[j]]) ) prob_other2 += cap;  //if it belongs to MU2
		}
	}// end of the for loop
		
	
	if(prob_other1 == 0 && prob_other2== 0){//which means all states in partition[0] belongs to either PV1 or PV2
		if(debug) printf("p_other1=p_other2=0, returns true.\n");
		return true;
	}

	int MAX= 66, num_header = 5, additional_arcs = 3; //additional arcs for the construction of the feasibility networks.
	//(h-l)^2 is the number of maximal transitions from post(s1) to post(s2). (h-l)*2 is the transitions from source or to sink. 5 is the header lines.
	int rows = (h-l)*(h-l) + (h-l)*2 + num_header + additional_arcs;
	char** inp_input = new char*[rows];
	for(int m=0; m<rows; m++) 
		inp_input[m] = new char[MAX];
	
	//		int product = prob1 * prob2;		
	/* check the reason for this number and write docu!!*/
	//the maximum capacity
	int max_cap = prob_other1*prob2 + prob_other2*prob1;		
	char* arc = new char[MAX];
	int index = num_header;
	
	//This time we insert the arcs from the sources and to the sink t in the problem description. Note I use the vector structure here. Therefore I can use the same arc array. The variable arc is pushed in the vector, everytime before it is changed. Otherwise it is probalimatic, since arc is a pointer, and it points everytime to a new thing.
	for(int j=0; j<left; j++){
		//s is the original source, S is the source of the transforemd graph                
		if( in_rel[left_states[j]][s2] ){ //s
			sprintf(arc, "a 3 %d %d", j+5, left_cap[j] * prob2 );//5 since 4 additional states S,T,s,t and the face that the states in the input inp starts with 0 instead of 0.
			strcpy(inp_input[index], arc);
			index++;
		}else{//S
			sprintf(arc, "a 1 %d %d", j+5, left_cap[j] * prob2 );
			strcpy(inp_input[index], arc);
			index++;
		}
		for(int k=0; k<right; k++){
			//if( SimulatesS1(left_states[j], right_states[k]) ){                          
			if( in_rel[left_states[j]][right_states[k]] ){
				sprintf(arc, "a %d %d %d", j+5, k+5+left, max_cap);
				strcpy(inp_input[index], arc);
				index++;
			}
		}
	}
	
	for(int j=0; j<right; j++){
		if( in_rel[s1][right_states[j]] ){ //t
			sprintf(arc, "a %d 4 %d", j+5+left, right_cap[j] * prob1);
			strcpy(inp_input[index], arc);
			index++;
		}else{//T
			sprintf(arc, "a %d 2 %d", j+5+left, right_cap[j] * prob1);
			strcpy(inp_input[index], arc);
			index++;
		}
	}// end of the arcs

	sprintf(arc, "a 1 4 %d",  prob_other2 * prob1);
	strcpy(inp_input[index], arc);
	index++;
	sprintf(arc, "a 4 3 %d", max_cap);
	strcpy(inp_input[index], arc);
	index++;
	sprintf(arc, "a 3 2 %d", prob_other1 * prob2);
	strcpy(inp_input[index], arc);
	index++;
	

	bool is_maxflow = false;
	if(debug) printf("^_*: c states s1,s2 (%d,%d)\n", s1 ,s2);
	strcpy(inp_input[0], "c set debug to see ^_*:");
	if(debug) printf("^_*: c prob_1=%d , prob_2=%d, prob_other1=%d, prob_other2=%d", prob1, prob2, prob_other1, prob_other2);
	strcpy(inp_input[1], "c set debug to see the probss after ^_*:");
	sprintf(arc, "p max %d %d", left + right + 4, index-num_header);//is 5, the number of head lines.
	strcpy(inp_input[2], arc);	
	strcpy(inp_input[3], "n 1 s");
	strcpy(inp_input[4], "n 2 t");
	
	if(debug){
		for(int i=0; i<index; i++){
			cout<<inp_input[i]<<endl;
		}
	}
	
	//if(max_cap== GetFlow(inp_input, index)) 
		is_maxflow = true;
	
	delete [] arc;		
	for(int m=0; m<rows; m++) { delete [] inp_input[m]; }
	delete [] inp_input;
	delete [] left_states;
	delete [] right_states;
	delete [] left_cap;
	delete [] right_cap;
	
	if(is_maxflow){   if(debug) printf("returns true(^_*)\n");   return true;  }
	else          {   if(debug) printf("returns false\n");       return false; }
}
	

/*bool
Relation:: CheckFeasibleFlow(int s1, int s2,  
				int* par_array, 
				int* par_array_index, 
				int prob1, 
				int prob2){
	//Now with help of the value gamma, we can deal with the first partition.
	int l = par_array_index[0];
	int h = par_array_index[1];  
	
	if(h-l <= 1) {
	if(debug) printf(" A1 contains at most one element. returns true\n");
	return true;
	}  
	
	//for the headers
	int left = 0, right = 0;
	int suc_s1_index = 0, suc_s2_index = 0;         
	int prob_other1 = 0, prob_other2 = 0;
	
	
	FILE *file=fopen("input_flow.txt","w");  
	string partition_arcs = ConstructTransitionArcs(s1, s2, par_array, par_array_index, 0);
	fprintf(file, "%s", partition_arcs.c_str() );
	fprintf(file, "* end of main transitions, i.e., without source and sink.\n");
	
	for(int i=l; i<h; i++){
	if(IsSuccessorOf(par_array[i], s1, suc_s1_index)){  
	int cap = (int) (non_zeros[suc_s1_index] * INPUT_PRECISION);
	
	//s is the original source, S is the source of the transforemd graph                
	if( SimulatedByS2(s2, par_array[i]) ){ // indicates par_array[i] \in PV1
		fprintf(file, "s * %d %d\n", par_array[i], cap);
	}else{ //par_array[i] \in Other1
		fprintf(file, "S * %d %d\n", par_array[i], cap);          
		prob_other1 += cap;
	}
	left++;
	}
	
	if(IsSuccessorOf(par_array[i], s2, suc_s2_index)){ 
	int cap = (int) (non_zeros[suc_s2_index] * INPUT_PRECISION);
	
	if( SimulatesS1(s1, par_array[i]) ){
		fprintf(file, "t %d * %d\n", par_array[i], cap);
	}else{
		fprintf(file, "T %d * %d\n", par_array[i], cap);
		prob_other2 += cap;
	}
	
	right++;
	}
	}// end of the for loop
	
	//we will add three additional transitions in the function CalculateFeasibleFlow (in the transformed graph)
	fprintf(file, "# end of the whole transitions.\n");     
	fprintf(file, "c states s1,s2 (%d, %d), special partition: number of successors s1: %d; of s2: %d\n", s1, s2, left, right);
	fprintf(file, "c prob_1=%d, prob_2=%d, prob_o1=%d, prob_o2=%d\n", prob1,prob2,prob_other1,prob_other2);
	//we need later calculate the number of transitions
	fprintf(file, "p max %d\n", left+right+4);
	fprintf(file, "n %d s\n", left+right+1);
	fprintf(file, "n %d t\n", left+right+4);
	
	fclose(file);
	
	
	if(prob_other1 == 0 && prob_other2== 0){//which means all states in partition[0] belongs to either PV1 or PV2
		if(debug) printf("p_other1=p_other2=0, returns true.\n");
		return true;
	}else{
		MaxFlow* mf = new MaxFlow(debug);      
		bool is_maxflow = mf->CalculateFeasibleFlow("input_flow.txt", prob1, prob2, left, right, 
								prob_other1, prob_other2);
		delete mf;
		
		if(is_maxflow){   if(debug) printf("returns true(^_*)\n");   return true;  }
		else          {   if(debug) printf("returns false\n");       return false; }
		
	}
	
}*/
	

/*bool
Relation:: CheckMaxFlow(int pair, int number_partitions, 
				int* par_array, 
				int* par_array_index, 
				int& p1, 
				int& p2)
{
	int s1 = relation[pair].x;
	int s2 = relation[pair].y;
	
	//dealing with the partitions except the frist which will be dealt with later. Recall that I put all s1 s2 related states in the first partition. 
	//It is not very good to use the variable l. How to differenciate 1 (one) and l (letter): by the colours!!!
	for(int i=1; i<number_partitions; i++){// number i denotes the i-th partition
	if(debug) printf("checking partition A%d, ", i);
	int prob1 = 0, prob2 = 0;
	
	int l = par_array_index[i];
	int h = par_array_index[i+1];  
	
	//the headers
	int left = 0, right = 0;
	int suc_s1_index, suc_s2_index = 0;         
	
	
	FILE *file=fopen("input_flow.txt","w");  
	string partition_arcs = ConstructTransitionArcs(pair, par_array, par_array_index, i);
	
	if(debug) printf("constructed the  main transitions\n");
	fprintf(file, "%s", (partition_arcs).c_str() );
	fprintf(file, "* end of main transitions, i.e., without source and sink.\n");
	
	//by the construction, states in Partition have the property: either P(s1,state)>0 or P(s2,state)>0    
	for(int j=l; j<h; j++){
	if(IsSuccessorOf(par_array[j], s1, suc_s1_index)){      //calculating prob1
		int cap = (int) (non_zeros[suc_s1_index] * INPUT_PRECISION);
		
		prob1 += cap;
		
		fprintf(file, "s * %d %d\n", par_array[j], cap);
		
		left++;
	}
	
	if(IsSuccessorOf(par_array[j], s2, suc_s2_index)){       //calculating prob2
		int cap = (int) (non_zeros[suc_s2_index] * INPUT_PRECISION);
		
		prob2 += cap;
		
		fprintf(file, "t %d * %d\n", par_array[j], cap);
		
		right++;
	}
	}// end of the for loop
	
	
	fprintf(file, "# end of the whole transitions.\n");     
	if(debug) printf("constructed the whole transitions. ");
	
	fprintf(file, "c states s1,s2 (%d, %d),%d-th partition: number of successors s1: %d; of s2: %d\n", s1, s2, i, left, right);
	fprintf(file, "c prob_1=%d, prob_2=%d\n", prob1,prob2);    //we need later calculate the number of transitions
	fprintf(file, "p max %d\n", left+right+2);
	fprintf(file, "n %d s\n", left+right+1);
	fprintf(file, "n %d t\n", left+right+2);
	
	fclose(file);
	if(debug) printf("complete the max flow input file\n");
	
	if((prob1 == 0) && (prob2 == 0)){
	if(debug) printf("p1=p2=0, error not possbile, if this lines appears, please check where are wrong and I just exit here, ups...file Relation.cc, function CheckMaxFlow() \n\n\n");
	exit(0);
	}else if(prob1 == 0){
	if(debug) printf("p1=0, p2=%d, returns false \n", prob2);
	return false;
	}else if(prob2==0){
	if(debug) printf("p1=%d, p2=0, returns false \n", prob1);
	return false;
	}else{
	bool is_maxflow = false;
	if(FLOW){  
		MaxFlow* mf = new MaxFlow(debug);      
		is_maxflow = mf->CalculateMaxFlow("input_flow.txt", prob1, prob2, left, right);
		delete mf;
	}else{
		LpSolver* lp = new LpSolver(debug);
		is_maxflow = lp->CalculateMaxFlow("input_flow.txt", prob1, prob2, left, right);
		delete lp;
	}
	
	if(!is_maxflow){
	if(debug) printf("returns false\n");
		return false;
	}else{
	if(debug) printf("ok for this partition. ");
		p1 = prob1;
		p2 = prob2;
	}
	
	}
	
	}//end of the special partitions
	
	if(debug) printf("CheckMaxFlow returns true. ");
	return true;
}*/

bool
Relation:: CheckFeasibleFlowViaLpSolver(int s1, int s2,  
			int* par_array, 
			int* par_array_index, 
			int prob1, 
			int prob2){
	//Now with help of the value gamma, we can deal with the first partition.
	int l = par_array_index[0];
	int h = par_array_index[1];  
	
	if(h-l <= 1) {
	if(debug) printf(" A1 contains at most one element. returns true\n");
	return true;
	}  
	
	
	
	//for the headers
	int left = 0, right = 0;
	int suc_s1_index = 0, suc_s2_index = 0;         
	int prob_other1 = 0, prob_other2 = 0;
	
	
	FILE *file=fopen("input_lp.txt","w");  
	string partition_arcs = ConstructTransitionArcs(s1, s2, par_array, par_array_index, 0);
	fprintf(file, "%s", partition_arcs.c_str() );
	fprintf(file, "* end of main transitions, i.e., without source and sink.\n");
	
	
	//the purpose of the copy is just to put first all transitions starts with s, then S. This will be needed later, to give an order for the transitions (LpSolver).  
	for(int i=l; i<h; i++){
	if(IsSuccessorOf(par_array[i], s1, suc_s1_index)){  
	int cap = (int) (non_zeros[suc_s1_index] * INPUT_PRECISION);
	
	//s is the original source, S is the source of the transforemd graph                
	if( SimulatedByS2(s2, par_array[i]) ){ // indicates par_array[i] \in PV1
		fprintf(file, "s * %d %d\n", par_array[i], cap);
	}
	left++;
	}
	
	if(IsSuccessorOf(par_array[i], s2, suc_s2_index)){ 
	int cap = (int) (non_zeros[suc_s2_index] * INPUT_PRECISION);
	
	if( SimulatesS1(s1, par_array[i]) ){
		fprintf(file, "t %d * %d\n", par_array[i], cap);
	}
	
	right++;
	}
	}// end of the for loop
	
	for(int i=l; i<h; i++){
		if(IsSuccessorOf(par_array[i], s1, suc_s1_index)){  
			int cap = (int) (non_zeros[suc_s1_index] * INPUT_PRECISION);
			
			//s is the original source, S is the source of the transforemd graph                
			if( !SimulatedByS2(s2, par_array[i]) ){ //par_array[i] \in Other1
				fprintf(file, "S * %d %d\n", par_array[i], cap);          
				prob_other1 += cap;
			}
		}
		
		if(IsSuccessorOf(par_array[i], s2, suc_s2_index)){ 
			int cap = (int) (non_zeros[suc_s2_index] * INPUT_PRECISION);
			
			if( !SimulatesS1(s1, par_array[i]) ){
				fprintf(file, "T %d * %d\n", par_array[i], cap);
				prob_other2 += cap;
			}
		}
	}// end of the for loop
	
	
	//we will add three additional transitions in the function CalculateFeasibleFlow (in the transformed graph)
	fprintf(file, "# end of the whole transitions.\n");     
	fprintf(file, "c states s1,s2 (%d, %d), special partition: number of successors s1: %d; of s2: %d\n", s1, s2, left, right);
	fprintf(file, "c prob_1=%d, prob_2=%d, prob_o1=%d, prob_o2=%d\n", prob1,prob2,prob_other1,prob_other2);
	//we need later calculate the number of transitions
	fprintf(file, "p max %d\n", left+right+4);
	//   fprintf(file, "n %d s\n", left+right+1);
	//   fprintf(file, "n %d t\n", left+right+4);
	
	fclose(file);
	
	
	if(prob_other1 == 0 && prob_other2== 0){//which means all states in partition[0] belongs to either PV1 or PV2
	if(debug) printf("p_other1=p_other2=0, returns true.\n");
	return true;
	}else{
	LpSolver* lp = new LpSolver(debug);      
	bool is_maxflow = lp->CalculateFeasibleFlow("input_lp.txt", prob1, prob2, left, right, 
							prob_other1, prob_other2);
	delete lp;
	
	if(is_maxflow){   if(debug) printf("returns true(^_*)\n");   return true;  }
	else          {   if(debug) printf("returns false\n");       return false; }
	}
	
	}

	
bool 
Relation:: ReducedWS(int s1, int s2){
	if(debug) printf("in function ReducedWS. Now we are dealing the case n==1! I would not suppose to see it!! And if it  does appears, either we deal with it later (the first option), or write it now.\n");
	exit(0);
}
	
bool 
Relation::PureStutterOne(int s1, int s2){
	//Test whether post(s1) \subseteq R^{-1}[s2]
	//l_tra: lower index for transition
	//l_rel: lower index for relation
	int l_tra = row_starts[s1];
	int h_tra = row_starts[s1+1];
	int l_brel= relation_backward_index[s2];
	int h_brel = relation_backward_index[s2+1];
	
	bool flag = false;
	
	for (int i = l_tra; i < h_tra; i++) {
		for(int j=l_brel; j<h_brel; j++){
			if(cols[i] == relation_backward[j]){
				flag = true;
				break;
			}
		}
		
		if (flag) flag = false;
		else return false;
	}
	
	return true;
}
	
bool
Relation:: PureStutterTwo(int s1, int s2){
	//Test whether post(s2) \subseteq R[s1]
	//l_tra: lower index for transition
	//l_rel: lower index for relation
	int l_tra = row_starts[s2];
	int h_tra = row_starts[s2+1];
	int l_frel= relation_forward_index[s1];
	int h_frel = relation_forward_index[s1+1];
	
	bool flag = false;
	
	for (int i = l_tra; i < h_tra; i++) {
		for(int j=l_frel; j<h_frel; j++){
			if(cols[i] == relation_backward[j]){
				flag = true;
				break;
			}
		}
		
		if (flag) flag = false;
		else return false;
	}
	
	return true;
}
	
int*
Relation:: ReachableStates(int state, int& size){
	int* explored = new int[n_states]; //indicates whether the state (in state_s) is already explored.
	for(int i=0; i<n_states; i++) explored[i] = 0;
	explored[state] = 1; //state is reachable from state itself.

	queue<int> qstates;	// used for the BFS algoriths.
	qstates.push(state);
	int processing = 1; // count the number of processed states. Increased after qstates.push().

	//BFS	
	while( !qstates.empty() ){	
		int last = qstates.front();
		qstates.pop();
		
		for (int i = row_starts[last]; i < row_starts[last+1]; i++) {
			if( !explored[cols[i]] ){
				explored[cols[i]] = 1;
				qstates.push(cols[i]);
				processing ++ ;
			}
		}
	}//end while loop
	
	int* r = new int[processing];
	size = processing;
	
	int index = 0;
	for(int i=0; i<n_states; i++){
		if(explored[i]){
			r[index] = i;
			index++;
		}
	}
	
	delete [] explored;
	
	return r;
}
