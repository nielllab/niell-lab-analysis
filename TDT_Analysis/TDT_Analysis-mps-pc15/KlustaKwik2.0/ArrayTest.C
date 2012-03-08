#include <fstream.h>
#include <iostream.h>
#include "Array.h"


main() {
	Array<int> a(10);
	Array<int> b(10);
	Array<int> c(20);
	int i;
	
	for (i=0; i<10; i++) a[i] = i;
	b=a;
	for (i=0; i<10; i++) cout << a[i] << endl;
	for (i=0; i<10; i++) cout << b[i] << endl;
	
	c=a;
}
	

/*
main() {
	::Array2<int> a2(3,4);
	int i, j;
	
	Array<int> a = a2[1];
	a[1] = 100;
	cout << a[1] << "\n\n\n";
	
	for(i=0; i<3; i++) for(j=0; j<4; j++) a2[i][j] = i + j;
	output();
}

void output() {
	int i,j;
	for(i=0; i<3; i++) {
		for(j=0; j<4 + (i==2); j++) cout << a2[i][j] << " ";
		cout << "\n";
	}
}
*/