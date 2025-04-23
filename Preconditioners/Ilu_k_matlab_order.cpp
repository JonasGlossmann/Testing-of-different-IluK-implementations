//
// Created by jonas on 11.03.2025.
//

#include "Ilu_k_matlab_order.h"

#include "Ilu_k.h"
#include "../misc/TimeControl.h"

DenseVector Ilu_k_matlab_order::apply(DenseVector residual) {

    std::vector<double> result = residual.data;

    //foward substitution
    int num_rows = L_row_ptr.size() -1;
    for (int row = 0; row < num_rows; ++row) {
        int next_row =L_row_ptr[row+1];

        for(int idx = L_row_ptr[row]; idx< next_row; ++idx) {
            result[row] -= L[idx]*result[L_col[idx]];

        }
    }

    //backward substitution
    num_rows = U_row_ptr.size() -2;
    for(int row = num_rows; row >-1 ; --row) {
        int prev_row = U_row_ptr[row]-1;

        for(int idx = U_row_ptr[row+1]-1; idx > prev_row; --idx) {
            int col = U_col[idx];

            if(row != col) {
                result[row] -= U[idx] * result[col];
            } else {
                result[row] = result[row]/U[idx];
            }
        }
    }

    DenseVector result2(result); //todo in later versions we will reduce this implementation to the basics to get better adaptivity
	return result2;
}

bool comp_first_greater(size_t val, const std::pair<size_t, int>& p)
{
	return (val < p.first);
}

Ilu_k_matlab_order::Ilu_k_matlab_order(COOMatrixSorted A_ext, int leve) {

	auto A_tmp = CSRMatrix(A_ext);
    /*
	 *  declare work variables
	 */
	int lfil;
	size_t N, nnzU, nnzL, nnzUold, nnzLold;
	static const int realloc_thresh = 0;
	typedef std::vector<std::pair<size_t,int> >::iterator vitr;
	typedef std::vector<std::pair<size_t,int> >::const_iterator c_vitr;

	// lfil is the level of fill parameter
	// N is the order of the matrix A
	// realloc_thresh (>= 0) is a constant that causes reallocation of
	// the L and U factors if nzmax(.) - nnz(.) > realloc_thresh

	/*
	 *  declare variables for sparse matrix A
	 */
	double *A;
	int *prowA, *pcolindexA;

    /*
	 *  declare variables for output
	 */
    int *prowL, *pcolindexL, *prowU, *pcolindexU;
	double *L, *U;


	A = &A_tmp.val[0];		   // real vector for A
	prowA = &A_tmp.colIndex[0];	   // row indices for elements of A
	pcolindexA = &A_tmp.rowPtr[0]; // index into the columns
	N = A_tmp.cols;           // number of columns of A

	// fill level lfil
	lfil = leve;
	++lfil; // use a shifted level of fill

	// level data structure: 0 <= level[i][j] <= lfil is the level of fill for
	// permitted fill-in entries
	// initially level[i][j] = 1 for all (i,j) such that A(i,j) ~= 0
	std::vector<std::vector<std::pair<size_t,int> > > level(N);
	level[0].reserve(pcolindexA[1]);
	for(size_t k = 0; k != pcolindexA[1]; ++k) {
		level[0].push_back(std::pair<size_t,int>(prowA[k], 1));
	}
	// symbolic phase
	// determine where fill-ins are permitted to occur based on the level of fill
	// count the number of nonzeros in L and U
	nnzU = pcolindexA[1];
	nnzL = 0;
	{
	int nzmax = N;
	// optimize memory usage in the case of no fill-ins, i.e., ILU(0)
	if(lfil == 1) {
		nzmax = pcolindexA[1];
		for(size_t j = 1; j != N; ++j) {
			nzmax = std::max(nzmax, pcolindexA[j + 1] - pcolindexA[j]);
		}
	}
	std::vector<int> work_lev(N, 0); // work vector for level values
	std::vector<size_t> work_idx(nzmax); // index array for work_lev
	for(size_t j = 1; j != N; ++j) {
		size_t nzmax = 0;
		// initialize level[j] work array
		for(size_t k = pcolindexA[j]; k != pcolindexA[j + 1]; ++k) {
			work_lev[prowA[k]] = 1;
			work_idx[nzmax++] = prowA[k];
		}
		size_t ii = 0;
		while(ii < nzmax && work_idx[ii] < j) {
			size_t i = work_idx[ii];
			++nnzL;
			// itrk->first points to the first column index in row i that is > i
			c_vitr itrk = upper_bound(level[i].begin(), level[i].end(), i, comp_first_greater);
			size_t nzmax_old = nzmax;
			while(itrk != level[i].end()) {
				int weight = work_lev[i] + itrk->second;
				// is (i, itrk->first) already a fill-in?
				if(work_lev[itrk->first] != 0) {
					work_lev[itrk->first] = std::min(weight, work_lev[itrk->first]);
				}
				else if(weight <= lfil) {
					work_lev[itrk->first] = weight;
					work_idx[nzmax++] = itrk->first;
				}
				++itrk;
			}
			// merge sorted sublists [0, nzmax_old) and [nzmax_old, nzmax) in work_idx
			// note that elements in positions 0,...,ii are already in their sorted
			// position
			if(nzmax > nzmax_old) {
				inplace_merge(work_idx.begin() + ii + 1, work_idx.begin() + nzmax_old, work_idx.begin() + nzmax);
			}
			++ii;
		}
		nnzU += nzmax;
		// fill level[j] of data structure
		level[j].reserve(nzmax);
		for(size_t k = 0; k != nzmax; ++k) {
			level[j].push_back(std::pair<size_t,int>(work_idx[k], work_lev[work_idx[k]]));
			work_lev[work_idx[k]] = 0;
		}
	}
	} // artificial scope to force destruction of work objects
	nnzU -= nnzL;
	nnzL += N; // include unit diagonal elements


	// create output matrices
	this->L.resize(nnzL,0);
	this->L_row_ptr.resize(N+1,0);
	this->L_col.resize(nnzL,0);
	L = &this->L[0];
	prowL = &this->L_col[0];
	pcolindexL = &this->L_row_ptr[0];

	this->U.resize(nnzU,0);
	this->U_row_ptr.resize(N+1,0);
	this->U_col.resize(nnzU,0);
	U = &this->U[0];
	prowU = &this->U_col[0];
	pcolindexU = &this->U_row_ptr[0];

	// first rows of L and U are known
	pcolindexU[0] = 0;
	pcolindexU[1] = pcolindexA[1];
	for(size_t j = pcolindexA[0]; j != pcolindexA[1]; ++j) {
		prowU[j] = prowA[j];
		U[j] = A[j];
	}


	pcolindexL[0] = 0;
	prowL[0] = 0;

	// numeric phase
	// compute the values of permitted fill-in entries and build L and U
	nnzLold = nnzL; nnzUold = nnzU;
	nnzL = 0; nnzU = pcolindexU[1];
	{
	std::vector<double> work(N, 0.0); // work array for rows of A
	std::vector<char> work_lev(N, 0); // work array for levels
	for(size_t j = 1; j != N; ++j) {
		// unpack jth row of A into the work array
		for(size_t k = pcolindexA[j]; k != pcolindexA[j + 1]; ++k) {
			work[prowA[k]] = A[k];
		}
		// unpack jth row of level into a work array
		for(c_vitr itr = level[j].begin(); itr != level[j].end(); ++itr) {
			work_lev[itr->first] = 1;
		}
		c_vitr itr = level[j].begin();
		while(itr != level[j].end() && itr->first < j) {
			double alpha = work[itr->first] / U[pcolindexU[itr->first]];
			work[itr->first] = alpha;
			for(size_t k = pcolindexU[itr->first] + 1; k != pcolindexU[itr->first + 1]; ++k) {
				// is prowU[k] a permitted fill-in?
				if(work_lev[prowU[k]] == 1) {
					/*todo debug code entfernen
					std::cout << prowU[k] << std::endl;
					if(prowU[k] == 102) {
						std::cout<<std::setprecision (53) <<work[prowU[k]] << "<- a_ij \n   aik-> "<< alpha << "  akj->" << U[k]<< std::endl;
					}*/
					work[prowU[k]] -= alpha * U[k];
				}
				// if memory is an issue you can avoid the work_lev array
				// by using a binary search, however, doing so may
				// significantly increase computing time, especially when
				// the level of fill is large
				// if(b_search(level[j], prowU[k])) {
				//	work[prowU[k]] -= alpha * U[k];
				// }
			}
			++itr;
		}
		// update row j of L and U and zero the work array
		itr = level[j].begin();
		while(itr != level[j].end()) {
			if(work[itr->first] != 0.0) {
				if(itr->first < j) {
					L[nnzL] = work[itr->first];
					prowL[nnzL++] = itr->first;
				}
				else {
					U[nnzU] = work[itr->first];
					prowU[nnzU++] = itr->first;
				}
				work[itr->first] = 0.0;
			}
			work_lev[itr->first] = 0;
			++itr;
		}
		// check that first nonzero in row j of U is a diagonal element
		prowL[nnzL] = j;
		pcolindexL[j + 1] = nnzL;
		pcolindexU[j + 1] = nnzU;
	}
	} // artificial scope to force destruction of work object


	// decrease the maximum number of nonzeros in L and U to nnzL and nnzU if necessary
	if(nnzLold - nnzL > realloc_thresh) {
		this->L.resize(nnzL);
		this->L_col.resize(nnzL);
	}
	if(nnzUold - nnzU > realloc_thresh) {
		this->U.resize(nnzU);
		this->U_col.resize(nnzU);
	}
}
