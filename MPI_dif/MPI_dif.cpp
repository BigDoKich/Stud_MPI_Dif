#include<iostream>
#include"mpi.h"
#include<fstream>
#include<cmath>
using namespace std;

void Fun(double* x, double* y, int n) {
	for (int i = 0; i < n; i++) {
		y[i] = x[i] * x[i] * x[i] * x[i];
	}
}

double Posled(long long n) {
	double* x = new double[n];
	double* y = new double[n];
	double p1, p2;
	for (int i = 0; i < n; i++) {
		x[i] = i / 10.;
	}
	Fun(x, y, n);

	double t1 = clock();
	double h = fabs(x[1] - x[0]);
	p1 = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);
	p2 = (2 * y[0] - 5 * y[1] + 4 * y[2] - y[3]) / (h * h);
	for (int i = 1; i < n - 1; i++) {
		p1 = (y[i + 1] - y[i - 1]) / (2 * h);
		p2 = (y[i - 1] - 2 * y[i] + y[i + 1]) / (h * h);
	}
	p1 = (y[n - 3] - 4 * y[n - 2] + 3 * y[n - 1]) / (2 * h);
	p2 = (-y[n - 4] + 4 * y[n - 3] - 5 * y[n - 2] + 2 * y[n - 1]) / (h * h);
	double t2 = clock();
	delete[] x;
	delete[] y;
	return (t2 - t1) / CLOCKS_PER_SEC;
}

int main(int argc, char** argv)
{
	const long long n = 10000000;
	double* x = new double[n];
	double* y = new double[n];
	double time_start, time_finish;
	double timeMPI = 0, timePosled = 0;
	ofstream fout{ "result.txt" };
	for (int i = 0; i < n; i++) {
		x[i] = i/10.;
	}
	Fun(x, y, n);
	MPI_Init(&argc, &argv);
	for (int p = 0; p < 100; p++)
	{
		double* p1 = new double[n];
		double* p2 = new double[n];
		double* p1H = new double[n];
		double* p2H = new double[n];
		timeMPI = 0;
		timePosled = 0;
		for (int k = 0; k < 10; k++)
		{
			time_start = MPI_Wtime();
			int rank;
			MPI_Status status;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if (rank == 0) {
				MPI_Recv(&p1[0], n / 2, MPI_DOUBLE, 2, 3, MPI_COMM_WORLD, &status);
				MPI_Recv(&p2[0], n / 2, MPI_DOUBLE, 2, 4, MPI_COMM_WORLD, &status);
				MPI_Recv(&p1H[(n - 1) / 2], n / 2, MPI_DOUBLE, 3, 5, MPI_COMM_WORLD, &status);
				MPI_Recv(&p2H[(n - 1) / 2], n / 2, MPI_DOUBLE, 3, 6, MPI_COMM_WORLD, &status);
				/*for (int i = (n - 1) / 2; i < n; i++) {
					p1[i] = p1H[i];
					p2[i] = p2H[i];
				}*/
				MPI_Recv(&p1[0], 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(&p2[0], 1, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, &status);
				MPI_Recv(&p1H[n - 1], 1, MPI_DOUBLE, 1, 7, MPI_COMM_WORLD, &status);
				MPI_Recv(&p2H[n - 1], 1, MPI_DOUBLE, 1, 8, MPI_COMM_WORLD, &status);
				/*for (int i = 0; i < n / 2; i++) {
					cout << x[i] << "\t" << y[i] << "\t" << p1[i] << "\t" << p2[i] << endl;
				}
				for (int i = n/2; i < n; i++) {
					cout << x[i] << "\t" << y[i] << "\t" << p1H[i] << "\t" << p2H[i] << endl;
				}*/
				//exit(0);
			}
			if (rank == 1) {
				double h = fabs(x[1] - x[0]);
				p1[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);
				p2[0] = (2 * y[0] - 5 * y[1] + 4 * y[2] - y[3]) / (h * h);
				MPI_Send(&p1[0], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
				MPI_Send(&p2[0], 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
				p1H[n - 1] = (y[n - 3] - 4 * y[n - 2] + 3 * y[n - 1]) / (2 * h);
				p2H[n - 1] = (-y[n - 4] + 4 * y[n - 3] - 5 * y[n - 2] + 2 * y[n - 1]) / (h * h);
				MPI_Send(&p1H[n - 1], 1, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);
				MPI_Send(&p2H[n - 1], 1, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
			}
			if (rank == 2) {
				double h = fabs(x[1] - x[0]);
				for (int i = 1; i <= (n - 1) / 2; i++) {
					p1[i] = (y[i + 1] - y[i - 1]) / (2 * h);
					p2[i] = (y[i - 1] - 2 * y[i] + y[i + 1]) / (h * h);
				}
				MPI_Send(&p1[0], n / 2, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
				MPI_Send(&p2[0], n / 2, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
			}
			if (rank == 3) {
				double h = fabs(x[1] - x[0]);
				for (int i = (n - 1) / 2; i < n - 1; i++) {
					p1H[i] = (y[i + 1] - y[i - 1]) / (2 * h);
					p2H[i] = (y[i - 1] - 2 * y[i] + y[i + 1]) / (h * h);
				}
				MPI_Send(&p1H[(n - 1) / 2], n / 2, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
				MPI_Send(&p2H[(n - 1) / 2], n / 2, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
			}
			time_finish = MPI_Wtime();
			if (rank == 0 ) {
				timeMPI += time_finish - time_start;
				timePosled += Posled(n);
			}
			if (rank == 0 && k == 9) {
				cout << "MPI time = " << timeMPI / 10. << endl;
				cout << "Posled time = " << timePosled / 10. << endl;
				cout << "Raznica MPI ot Posled = " << ((timePosled / 10.) / ((timeMPI / 10.) / 100)) - 100 << "%" << endl;
				fout << p + 1 << "\t" << ((timePosled / 10.) / ((timeMPI / 10.) / 100)) - 100 << endl;
			}
		}
		delete[]p1;
		delete[]p2;
		delete[]p1H;
		delete[]p2H;
	}
	MPI_Finalize();
	delete[]x;
	delete[]y;
	return 0;
}