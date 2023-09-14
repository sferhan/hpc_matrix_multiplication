#include <iostream>

using namespace std;

/*
Applies col major transform.
Takes the row-index @row and column-index @col of a n x n matrix and returns
 the position of this value in a 1-D array assuming col-major scheme.
*/
int col_major_transform(int row, int col, int n) {
    return row + (col * n);
}