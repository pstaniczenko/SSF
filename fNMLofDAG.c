// Compile:
// gcc fNMLofDAG.c -o fNMLGet -lgsl -lgslcblas -O3 -Wall  doesn't work, so try:
// gcc fNMLofDAG.c -o fNMLGet -L/usr/local/lib -lgsl -lgslcblas -O3 -Wall -I/usr/local/include

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_short.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_short.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>

// print a matrix of integers
int PrintMatrixShort(FILE * F, gsl_matrix_short * MAT){
  int NR=MAT->size1;
  int NC=MAT->size2;
  int i,j;
  for (i=0;i<NR;i++){
    for (j=0;j<NC;j++){
      fprintf(F, "%d ",gsl_matrix_short_get(MAT,i,j));
    }
    fprintf(F, "\n");
  }
  fprintf(F, "\n");
  return 0;
}

int PrintVector(FILE * F, gsl_vector * VEC){
  int NR=VEC->size;
  int i;
  for (i=0;i<NR;i++){
      fprintf(F, "%1.2f ",gsl_vector_get(VEC,i));
  }
  fprintf(F, "\n");
  return 0;
}

int PrintVectorShort(FILE * F, gsl_vector_short * VEC){
  int NR=VEC->size;
  int i;
  for (i=0;i<NR;i++){
      fprintf(F, "%d ",gsl_vector_short_get(VEC,i));
  }
  fprintf(F, "\n");
  return 0;
}

// return the complexity of a random graph of size X
double complexity(int X){
	if (X == 0) return 0.0;
	// for small X, use the exact value
	double XX = (double) X;
	if (XX < 100.) return log(exp(XX + log(gsl_sf_gamma_inc(XX,XX)) -
						(XX - 1.0) * log(XX)) + 1.0)/log(2);
	// for large X, use the approximation
	return(log(1. +
				sqrt(M_PI * XX/2.) - 1./3. +
				sqrt(2. * M_PI) / (24. * sqrt(XX)) 
				- 4. / (135. * XX) +
				sqrt(2. * M_PI) / (576. * sqrt(XX * XX * XX)) +
				8. / (2835. * XX * XX))/log(2));
}

// NML for a slice of matrix containing Ones 1s and Zeros 0s
double NMLSlice(int Ones, int Zeros){
	double NML = -complexity(Ones + Zeros);
	double p = 0.0;
	// if either Ones or Zeros = 0, the log_likelihood is 0
	if ((Ones > 0) && (Zeros > 0)){
		// compute maximum log likelihood
		p = (float) Ones / (float)(Ones + Zeros);
		NML += ((float) Ones) * log(p)/log(2) + ((float) Zeros) * (log(1. - p)/log(2));
	}
	return NML;
}

int build_random_DAG(int num_nodes, int num_links, gsl_vector_short * myDAG, gsl_rng * r){
	int hier[num_nodes];
	int tot_links = num_nodes * (num_nodes - 1) / 2;
	int lower_tri[tot_links];
	int i;
	// fill the hierarchy
	for (i = 0; i < num_nodes; i++){
		hier[i] = i;
	}
	// shuffle the hierarchy
	gsl_ran_shuffle (r, hier, num_nodes, sizeof (int));
	// fill the links
	for (i = 0; i < tot_links; i++){
		if (i < num_links){
			lower_tri[i] = 1;
		}
		else{
			lower_tri[i] = 0;
		}
	}
	// shuffle the links
	gsl_ran_shuffle (r, lower_tri, tot_links, sizeof (int));
	// fill in the DAG
	for (i = 0; i < (num_nodes + tot_links); i++){
		if (i < num_nodes){
			gsl_vector_short_set(myDAG, i, hier[i]);
		}
		else{
			gsl_vector_short_set(myDAG, i, lower_tri[i - num_nodes]);
		}
	}
	return 0;
}

int mutate_hierarchy(int num_nodes, gsl_vector_short * myDAG, gsl_rng * r){
	int x = gsl_rng_uniform_int(r, num_nodes);
	int y = gsl_rng_uniform_int(r, num_nodes);
	int tmp = gsl_vector_short_get(myDAG, x);
	gsl_vector_short_set(myDAG, x, gsl_vector_short_get(myDAG, y));
	gsl_vector_short_set(myDAG, y, tmp);
	return 0;
}

int mutate_links(int num_nodes, gsl_vector_short * myDAG, gsl_rng * r){
	int tot_links = num_nodes * (num_nodes - 1) / 2;
	int x = gsl_rng_uniform_int(r, tot_links);
	int y = gsl_rng_uniform_int(r, tot_links);
	int tmp = gsl_vector_short_get(myDAG, x + num_nodes);
	gsl_vector_short_set(myDAG, x + num_nodes, gsl_vector_short_get(myDAG, y + num_nodes));
	gsl_vector_short_set(myDAG, y + num_nodes, tmp);
	return 0;
}

// reconstruct DAG matrix from DAG vector
int build_adj(int num_nodes, gsl_vector_short * myDAG, gsl_matrix_short * A){
	int col, row, k;
	k = 0;
	for (col = 0; col < (num_nodes - 1); col++){
		for (row = (col + 1); row < num_nodes; row++){
			gsl_matrix_short_set(A, row, col, gsl_vector_short_get(myDAG, k + num_nodes));
			k++;
		}
	}
	return 0;
}

double getSwitches(int onesA, int zerosA,
				   int onesB, int zerosB,
				   double pA, double pB,
				   int codeA, int codeB,
				   int model){
	if (model < 3){
		// is A a "superset" of B?
		if (codeA > codeB){
			if (pA < pB){
				// we need switches
				return ceil(pB * ((double) onesA + zerosA) - (double) onesA);
			}
		}
		else{
			// B is a "superset" of A
			if (pB < pA){
				// we need switches
				return ceil(pA * ((double) onesB + zerosB) - (double) onesB);
			}
		}
	}
	else{
		// is A a superset of B?
		if ((codeA & codeB) == codeB){
			if (pA < pB){
				return ceil(pB * ((double) onesA + zerosA) - (double) onesA);
			}
		}
		else{
			// B is a superset of A
			if ((codeA & codeB) == codeA){
				if (pB < pA){
					// we need switches
					return ceil(pA * ((double) onesB + zerosB) - (double) onesB);
				}
			}
		}
		
	}
	return 0.0;
}

double fNML(int num_nodes,
			gsl_vector_short * DAG,
			gsl_matrix_short * mDAG,
			int model,
			gsl_matrix_short * Adj,
			gsl_vector * Count,
			gsl_vector * vfNML, 
			gsl_vector * vErrors, 
			gsl_vector * vSwitches, 
			int myvindex){
	// build the matrix corresponding to the DAG
	build_adj(num_nodes, DAG, mDAG);
	// count the number of parents for each node
	int col, row, k, k2, par;
	int col_label;
	int parent_label;
	int current_parent_code;
	int unique_states = 0;
	int num_rows = Adj->size1;

	int ones[num_rows];
	int zeros[num_rows];
	double ps[num_rows];
	
	unsigned int parentcode[num_rows];
	unsigned int parents_state[num_rows];

	double tot_fNML = 0.0;
	double col_fNML = 0.0;
	
	double col_errors = 0.0;
	double tot_errors = 0.0;
	
	double col_switches = 0.0;
	double tot_switches = 0.0;

	// for each column
	for (col = 0; col < num_nodes; col++){
	   	// find the label
		col_label = gsl_vector_short_get(DAG, col);
		// initialize things
		for (k = 0; k < num_rows; k++){
			ones[k] = 0;
			zeros[k] = 0;
			ps[k] = 0.0;
			parentcode[k] = -1;
			if (model == 0){
				parents_state[k] = 1;
			}
			else{
				parents_state[k] = 0;
			}
		}
		// for each possible parent
		for (par = 0; par < num_nodes; par++){
			if (gsl_matrix_short_get(mDAG, par, col) > 0){
				// par is a parent
				parent_label = gsl_vector_short_get(DAG, par);
				// update the state of the parents
				for (k = 0; k < num_rows; k++){
					if (model == 0) {
						parents_state[k] *= gsl_matrix_short_get(Adj, k, parent_label);
					}
					if (model == 1) {
						parents_state[k] += gsl_matrix_short_get(Adj, k, parent_label);
						if (parents_state[k] > 1) parents_state[k] = 1;
					}
					if (model == 2){
						parents_state[k] += gsl_matrix_short_get(Adj, k, parent_label);
					}
					if (model == 3){
						// MAX 32 COLUMNS!!!!!
						parents_state[k] |= ((1 << parent_label) * gsl_matrix_short_get(Adj, k, parent_label)) ;
					}
				}
			}
		}
		// now that each row has a certain parent state, find unique
		// states and set the zeros and ones
		int found;
		unique_states = 0;
		for (row = 0; row < num_rows; row++){
			current_parent_code = parents_state[row];
			// try to find the parent code
			found = 0;
			for (k = 0; k < unique_states; k++){
				if (parentcode[k] == current_parent_code){
					// found the code!
					found++;
					if (gsl_matrix_short_get(Adj, row, col_label) > 0){
						ones[k] += (int) gsl_vector_get(Count, row);
					}
					else{
						zeros[k] += (int) gsl_vector_get(Count, row);
					}
				}
			}
			if (found == 0){
				// never seen the code before
				parentcode[unique_states] = current_parent_code;
				if (gsl_matrix_short_get(Adj, row, col_label) > 0){
					ones[unique_states] += (int) gsl_vector_get(Count, row);
				}
				else{
					zeros[unique_states] += (int) gsl_vector_get(Count, row);
				}
				unique_states++;
			}
		}
		// compute the fNML for the col
		col_fNML = 0.0;
		for (k = 0; k < unique_states; k++){
			ps[k] = ((double) ones[k]) / ((double) (ones[k] + zeros[k]));
			col_fNML += NMLSlice(ones[k], zeros[k]);
		}
		// compute errors and switches
		col_errors = 0.0;
		col_switches = 0.0;
		for (k = 0; k < unique_states; k++){
			if (parentcode[k] == 0){
				// ones[k] are "errors"
				if (unique_states > 1){
					col_errors += ones[k];
				}
			}
			for (k2 = k; k2 < unique_states; k2++){
				if (k != k2){
					// get number of "switches"
					col_switches += getSwitches(ones[k], zeros[k],
												ones[k2], zeros[k2],
												ps[k], ps[k2],
												parentcode[k], parentcode[k2],
												model);
				}
			}
		}
		tot_errors += col_errors;
		tot_switches += col_switches;
		tot_fNML += col_fNML;
	}
	// set the values
	gsl_vector_set(vfNML, myvindex, tot_fNML);
	gsl_vector_set(vErrors, myvindex, tot_errors);
	gsl_vector_set(vSwitches, myvindex, tot_switches);
	return tot_fNML;
}

int main(int argc, char *argv[]){
	// Command arguments
	// ./fNMLGet 128 7 Pew2008-Q5-22429-7-Corrected.txt-z-128-7 myDAG.txt 0 5 10
	// Input parameters
	// 128 number of rows zipped matrix
	// 7 number of cols
	// Pew2008-Q5-22429-7-Corrected.txt-z-128-7 zipped matrix
    // myDAG.txt
	// 0 model: 0 -> And, 1 -> Or, 2 -> Sum, 3 -> Full
    // 123 random seed
	// 5% tolerance for errors
	// 10% tolerance for nonpositive

	int numR = atoi(argv[1]); // Num rows
	int numC = atoi(argv[2]); // Num cols
	char * ZipFileName = argv[3]; // File Name
    char * DAGFileName = argv[4]; // File Name
	int model = atoi(argv[5]); // Model
    int seed = atoi(argv[6]); // Random seed
	double tolerance_errors = atof(argv[7]) / 100.0; // do not penalize if
											 // errors are less than
											 // tolerance % of the ones
	double tolerance_switches = atof(argv[8]) / 100.0;
	int i, j;
	
    int second = 1; // only considering one DAG
    int DAGreducedlength = numC + numC * (numC - 1) / 2; // number of elements in DAG reduced form

	// Set up the random number generator and seed it
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, seed);
	
	// read zipped matrix
	fprintf(stderr, "Reading Data....\n");
	gsl_matrix * Z = gsl_matrix_calloc(numR, numC + 1);
	FILE * F;
	F = fopen(ZipFileName, "rb");
	gsl_matrix_fscanf(F, Z);
	fclose(F);
    
    // read DAG matrix
    fprintf(stderr, "Reading DAG....\n");
	gsl_vector * DAGphillip = gsl_vector_calloc(DAGreducedlength);
	FILE * FF;
	FF = fopen(DAGFileName, "r");
	gsl_vector_fscanf(FF, DAGphillip);
	fclose(FF);
	
	fprintf(stderr, "Build Matrices....\n");
	// build adjacency matrix and count vector
	gsl_matrix_short * Adj = gsl_matrix_short_calloc(numR, numC);
	gsl_vector * Count = gsl_vector_calloc(numR);
	double TotOnes = 0.0;
	double TotZeros = 0.0;
	for (i = 0; i < numR; i++){
		gsl_vector_set(Count, i, gsl_matrix_get(Z, i, numC));
		for (j=0; j < numC; j++){
			gsl_matrix_short_set(Adj, i, j, (int) gsl_matrix_get(Z, i, j));
			if (gsl_matrix_short_get(Adj, i, j) > 0){
				TotOnes += gsl_vector_get(Count, i);
			}
			else{
				TotZeros += gsl_vector_get(Count, i);
			}
		}
	}
    
    // build DAG matrix
    gsl_vector_short * DAGphillipshort = gsl_vector_short_calloc(DAGreducedlength);
    for (i = 0; i < DAGreducedlength; i++){
        gsl_vector_short_set(DAGphillipshort, i, (int) gsl_vector_get(DAGphillip, i));
    }


	gsl_vector * vfNML = gsl_vector_calloc(second);
	gsl_vector * vErrors = gsl_vector_calloc(second);
	gsl_vector * vSwitches = gsl_vector_calloc(second);
	
	fprintf(stderr, "DAG\n");
	PrintVectorShort(stderr, DAGphillipshort);

	gsl_matrix_short * matDAG = gsl_matrix_short_calloc(numC, numC);
	build_adj(numC, DAGphillipshort, matDAG);
	fprintf(stderr, "Corresponding matrix DAG\n");
	PrintMatrixShort(stderr, matDAG);
    
    double MYPENALIZE = -500.0;
    
    double finalfNML = fNML(numC, DAGphillipshort, matDAG, model, Adj, Count, vfNML, vErrors, vSwitches, 0);
    double bestErrors = gsl_vector_get(vErrors, 0);
    double bestSwitches = gsl_vector_get(vSwitches, 0);
    double bestAccept = 1;
    
    if ((bestErrors / TotOnes) > tolerance_errors){
		finalfNML += bestErrors * MYPENALIZE;
		bestAccept = 0;
	}
	// check if switches are tolerable
	if ((bestSwitches / TotZeros) > tolerance_switches){
		finalfNML += bestSwitches * MYPENALIZE;
		bestAccept = 0;
	}
	
	char outfile[1000];
	sprintf(outfile, "%s_knownDAGfNML_%f_Model_%d_E_%d_S_%d_Accept_%1.0f_%1.3f_%1.3f.cout", ZipFileName, finalfNML, model,
			(int) bestErrors, (int) bestSwitches, bestAccept, tolerance_errors, tolerance_switches);
	F = fopen(outfile, "w");
	PrintVectorShort(F, DAGphillipshort);
	fclose(F);
	// free memory
	gsl_vector_short_free(DAGphillipshort);
	gsl_matrix_short_free(Adj);
	gsl_matrix_free(Z);
	gsl_matrix_short_free(matDAG);
	gsl_vector_free(Count);
	gsl_rng_free (r);
	return 0;
}
