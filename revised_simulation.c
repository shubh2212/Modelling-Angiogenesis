#include <stdio.h>
#include <math.h>
#include <float.h>
#define row 101
#define col 101
#define tsteps 2
#define time_ 1
#define h (1/((double)(row-1)))
#define dt (0.0010)
#define pi 3.14159265359
// For below constants refer Pg. 24 of anderson et al.
#define D  (0.00035)
#define alpha 0.6
#define chi 0.38
#define rho 0.34
#define beta 0.05
#define gamma 0.1
#define neta 0.1
#define epsilon1 0.45
#define epsilon2 0.45
#define epsilon3 0.001
#define k 0.75
#define lambda dt / ((double)(h*h))
#define SCALE 1.0
/*
 * Define any other constants if applicable
 *
*/

/* Global Variables */

double density[row][col][tsteps];	//n in our eqn
double concentration[row][col][tsteps];	//c in our eqn
double f[row][col][tsteps];		//f in our eqn

double c_ij[row][col];
double d_ij[row][col];
double c_ij_time_n[row][col];
double d_ij_time_n[row][col];
double epsilon_ij[row][col];
double epsilon_ij_time_n[row][col];

double nine_pts[row][col];

double Q2[row][col]; //5 points terms
double Q4[row][col]; //5 points terms
double Q5[row][col]; //5 points terms
double Q6[row][col]; //5 points terms
double Q8[row][col]; //5 points terms
double err = 0.5E-4;

void  intialize() {
	int i,j,k1;
	for(i=0;i<row;i++) {
		for(j=0;j<col;j++) {
			c_ij[i][j] = 0.0;
			d_ij[i][j] = 0.0;
			c_ij_time_n[i][j] = 0.0;
			d_ij_time_n[i][j] = 0.0;
			epsilon_ij[i][j] = 0.0;
			epsilon_ij_time_n[i][j] = 0.0;
			nine_pts[i][j]=0.0;
			Q2[i][j] = 0.0;
			Q4[i][j] = 0.0;
			Q5[i][j] = 0.0;
			Q6[i][j] = 0.0;
			Q8[i][j] = 0.0;
			for(k1=0;k1<tsteps;k1++) {
				density[i][j][k1] = 0;
				concentration[i][j][k1] = 0;
				f[i][j][k1] = 0;
			}
		}
	}
	return;	
}
void intialize_coeff() {
			int x,y;
		for(y=0;y<col;y++){
			c_ij[0][y] = ( ( chi/(1+alpha*concentration[0][y][time_]  ) )
						 * ((concentration[1][y][time_] - concentration[0][y][time_])/h) 
						 + rho * ((f[1][y][time_] - f[0][y][time_])/h) 
						 )/D;
			c_ij[row-1][y] = ( ( chi/(1+alpha*concentration[row-1][y][time_]  ) )
						 * ((concentration[row-1][y][time_] - concentration[row-2][y][time_])/h) 
						 + rho * ((f[row-1][y][time_] - f[row-2][y][time_])/h) 
						 )/D;
			c_ij_time_n[0][y] = ( ( chi/(1+alpha*concentration[0][y][time_-1]  ) )
								 * ((concentration[1][y][time_-1] - concentration[0][y][time_-1])/h) 
								 + rho * ((f[1][y][time_-1] - f[0][y][time_-1])/h) 
								 )/D;
			c_ij_time_n[row-1][y] = ( ( chi/(1+alpha*concentration[row-1][y][time_-1]  ) )
								 * ((concentration[row-1][y][time_-1] - concentration[row-2][y][time_-1])/h) 
								 + rho * ((f[row-1][y][time_-1] - f[row-2][y][time_-1])/h) 
								 )/D;	
		}
		for(x=1;x<row-1;x++){
			c_ij[x][0] = ( ( chi/(1+alpha*concentration[x][0][time_]  ) )
						 * ((concentration[x+1][0][time_] - concentration[x-1][0][time_])/(2*h) ) 
						 + rho * ((f[x+1][0][time_] - f[x-1][0][time_])/(2*h)) 
						 )/D;
			c_ij[x][col-1] = ( ( chi/(1+alpha*concentration[x][col-1][time_]  ) )
						 * ((concentration[x+1][col-1][time_] - concentration[x-1][col-1][time_])/(2*h)) 
						 + rho * ((f[x+1][col-1][time_] - f[x-1][col-1][time_])/(2*h)) 
						 )/D;
			c_ij_time_n[x][0] = ( ( chi/(1+alpha*concentration[x][0][time_-1]  ) )
								 * ((concentration[x+1][0][time_-1] - concentration[x-1][0][time_-1])/(2*h)) 
								 + rho * ((f[x+1][1][time_-1] - f[x-1][0][time_-1])/(2*h)) 
								 )/D;
			c_ij_time_n[x][col-1] = ( ( chi/(1+alpha*concentration[x][col-1][time_-1]  ) )
								 * ((concentration[x+1][col-1][time_-1] - concentration[x-1][col-1][time_-1])/(2*h)) 
								 + rho * ((f[x+1][col-1][time_-1] - f[x-1][col-1][time_-1])/(2*h)) 
								 )/D;	
		}

		for(x=0;x<row;x++){
			d_ij[x][0] = ( ( chi/(1+alpha*concentration[x][0][time_]  ) )
						 * ((concentration[x][1][time_] - concentration[x][0][time_])/h) 
						 + rho * ((f[x][1][time_] - f[x][0][time_])/h) 
						 )/D;
			d_ij[x][col-1] = ( ( chi/(1+alpha*concentration[x][col-1][time_]  ) )
						 * ((concentration[x][col-1][time_] - concentration[x][col-2][time_])/h) 
						 + rho * ((f[x][col-1][time_] - f[x][col-2][time_])/h) 
						 )/D;
			d_ij_time_n[x][0] = ( ( chi/(1+alpha*concentration[x][0][time_-1]  ) )
								 * ((concentration[x][1][time_-1] - concentration[x][0][time_-1])/h) 
								 + rho * ((f[x][1][time_-1] - f[x][0][time_-1])/h) 
								 )/D;
			d_ij_time_n[x][col-1] = ( ( chi/(1+alpha*concentration[x][col-1][time_-1]  ) )
								 * ((concentration[x][col-1][time_-1] - concentration[x][col-2][time_-1])/h) 
								 + rho * ((f[x][col-1][time_-1] - f[x][col-2][time_-1])/h) 
								 )/D;	
		}

		for(y=1;y<col-1;y++){
			d_ij[0][y] = ( ( chi/(1+alpha*concentration[0][y][time_]  ) )
						 * ((concentration[0][y+1][time_] - concentration[0][y-1][time_])/(2*h)) 
						 + rho * ((f[0][y+1][time_] - f[0][y-1][time_])/(2*h)) 
						 )/D;
			d_ij[row-1][y] = ( ( chi/(1+alpha*concentration[row-1][y][time_]  ) )
						 * ((concentration[row-1][y+1][time_] - concentration[row-1][y-1][time_])/(2*h)) 
						 + rho * ((f[row-1][y+1][time_] - f[row-1][y-1][time_])/(2*h)) 
						 )/D;
			d_ij_time_n[0][y] = ( ( chi/(1+alpha*concentration[0][y][time_-1]  ) )
								 * ((concentration[0][y+1][time_-1] - concentration[0][y-1][time_-1])/(2*h)) 
								 + rho * ((f[0][y+1][time_-1] - f[0][y-1][time_-1])/(2*h)) 
								 )/D;
			d_ij_time_n[row-1][y] = ( ( chi/(1+alpha*concentration[row-1][y][time_-1]  ) )
								 * ((concentration[row-1][y+1][time_-1] - concentration[row-1][y-1][time_-1])/(2*h)) 
								 + rho * ((f[row-1][y+1][time_-1] - f[row-1][y-1][time_-1])/(2*h)) 
								 )/D;	
		}

		//Finally computing
		for(x=1;x<row-1;x++) {
			for(y=1;y<col-1;y++){
				c_ij[x][y] = ( ( chi/(1+alpha*concentration[x][y][time_]  ) ) 
							* ((concentration[x+1][y][time_] - concentration[x-1][y][time_] )/(2*h)) 
							+ rho*( (f[x+1][y][time_] - f[x-1][y][time_])/(2*h) ) 
							)/D;
				c_ij_time_n[x][y] = ( ( chi/(1+alpha*concentration[x][y][time_-1]  ) ) 
							* ((concentration[x+1][y][time_-1] - concentration[x-1][y][time_-1] )/(2*h)) 
							+ rho*( (f[x+1][y][time_-1] - f[x-1][y][time_-1])/(2*h) ) 
							)/D;

				d_ij[x][y] = ( ( chi/(1+alpha*concentration[x][y][time_]  ) ) 
						  	* ((concentration[x][y+1][time_] - concentration[x][y-1][time_] )/(2*h)) 
						  	+ rho*( (f[x][y+1][time_]- f[x][y-1][time_])/(2*h) ) 
						  	)/D;
				d_ij_time_n[x][y] = ( ( chi/(1+alpha*concentration[x][y][time_-1]  ) ) 
						  	* ((concentration[x][y+1][time_-1] - concentration[x][y-1][time_-1] )/(2*h)) 
						  	+ rho*( (f[x][y+1][time_-1]- f[x][y-1][time_-1])/(2*h) ) 
						  	)/D;
			}
		}


		//Boundary for epsilon ij

		epsilon_ij[0][0] =  ( ( c_ij[1][0] - c_ij[0][0] ) / h) + ( ( d_ij[0][1] - d_ij[0][0] ) / h);
		epsilon_ij[row-1][0] = ( ( c_ij[row-1][0] - c_ij[row-2][0] ) / h) + ( (d_ij[row-1][1] - d_ij[row-1][0] ) / h);
		epsilon_ij[0][col-1] = ( (c_ij[1][col-1] - c_ij[0][col-1] ) / h) + ( (d_ij[0][col-1] - d_ij[0][col-2] ) / h); 
		epsilon_ij[row-1][col-1] =  ( (c_ij[row-1][col-1] - c_ij[row-2][col-1] ) / h) + ( (d_ij[row-1][col-1] - d_ij[row-1][col-2] ) / h);  

		for(y=1;y<col-1;y++){
			double dho_c_x = ( (c_ij[1][y] - c_ij[0][y] ) / h);
			double dho_d_by_y = ( d_ij[0][y+1] - d_ij[0][y-1] ) / (2*h);
			double dho_c_x_time_n = ( (c_ij_time_n[1][y] - c_ij_time_n[0][y] ) / h);
			double dho_d_by_y_time_n = ( d_ij_time_n[0][y+1] - d_ij_time_n[0][y-1] ) / (2*h);
			epsilon_ij[0][y] = dho_c_x + dho_d_by_y;
			epsilon_ij_time_n[0][y] = dho_c_x_time_n + dho_d_by_y_time_n;

			dho_c_x = ( (c_ij[row-1][y] - c_ij[row-2][y] ) / h);
			dho_d_by_y = ( d_ij[row-1][y+1] - d_ij[row-1][y-1] ) / (2*h);
			dho_c_x_time_n = ( (c_ij_time_n[row-1][y] - c_ij_time_n[row-2][y] ) / h);
			dho_d_by_y_time_n = ( d_ij_time_n[row-1][y+1] - d_ij_time_n[row-1][y-1] ) / (2*h);
			epsilon_ij[row-1][y] = dho_c_x + dho_d_by_y;
			epsilon_ij_time_n[row-1][y] = dho_c_x_time_n + dho_d_by_y_time_n;
		}
		for(x=1;x<row-1;x++){
			double dho_c_x = ( c_ij[x+1][0] - c_ij[x-1][0] ) / (2*h);
			double dho_d_by_y = ( (d_ij[x][1] - d_ij[x][0] ) / h);
			double dho_c_x_time_n = ( c_ij_time_n[x+1][0] - c_ij_time_n[x-1][0] ) / (2*h);
			double dho_d_by_y_time_n = ( (d_ij_time_n[x][1] - d_ij_time_n[x][0] ) / h);
			epsilon_ij[x][0] = dho_c_x + dho_d_by_y;
			epsilon_ij_time_n[x][0] = dho_c_x_time_n + dho_d_by_y_time_n;

			dho_c_x = ( c_ij[x+1][col-1] - c_ij[x-1][col-1] ) / (2*h);
			dho_d_by_y = ( d_ij[x][col-1] - d_ij[x][col-2] ) / (2*h);
			dho_c_x_time_n = ( c_ij_time_n[x+1][col-1] - c_ij_time_n[x-1][col-1] ) / (2*h);
			dho_d_by_y_time_n = ( d_ij_time_n[x][col-1] - d_ij_time_n[x][col-2] ) / (2*h);
			epsilon_ij[x][col-1] = dho_c_x + dho_d_by_y;
			epsilon_ij_time_n[x][col-1] = dho_c_x_time_n + dho_d_by_y_time_n;
		}
		for(x=1;x<row-1;x++){
			for(y=1;y<col-1;y++){
				double dho_c_x = ( c_ij[x+1][y] - c_ij[x-1][y] ) / (2*h);
				double dho_d_by_y = ( d_ij[x][y+1] - d_ij[x][y-1] ) / (2*h);
				double dho_c_x_time_n = ( c_ij_time_n[x+1][y] - c_ij_time_n[x-1][y] ) / (2*h);
				double dho_d_by_y_time_n = ( d_ij_time_n[x][y+1] - d_ij_time_n[x][y-1] ) / (2*h);
				epsilon_ij[x][y] = dho_c_x + dho_d_by_y;
				epsilon_ij_time_n[x][y] = dho_c_x_time_n + dho_d_by_y_time_n;
			}
		}
}
void initiate_boundary_conditions() {
	int i,j;
	for(i=0;i<row;i++)	{
		for(j=0;j<col;j++){
			concentration[i][j][0] = exp( -((1.0-(double)i*h) * (1.0-(double)i*h))/epsilon1 );
			f[i][j][0] = k * exp( -(((double)i*h) * ((double)i*h))/epsilon2 );
			density[i][j][0] = sin(6*pi*j*h)*sin(6*pi*j*h) * exp( -(((double)i*h) * ((double)i*h))/epsilon3 );
		}
	}
}

/*void initiate_no_flux_condition() {
	density[0][0][time_] = density[0][0][time_-1];
	density[row-1][0][time_] = density[row-1][0][time_-1];
	int i,j;
	for(j=0;j<col-1;j++){
		double temp_1 = ( ( chi/(1+alpha*concentration[0][j][time_]  ) )
						 * ((concentration[0][j+1][time_] - concentration[0][j][time_])) 
						 + rho * ((f[0][j+1][time_] - f[0][j][time_]))
				  )/D;
		density[0][j+1][time_] = density[0][j][time_] + temp_1*density[0][j][time_];
		double temp_2 = ( ( chi/(1+alpha*concentration[row-1][j][time_]  ) )
						 * ((concentration[row-1][j+1][time_] - concentration[row-1][j][time_])) 
						 + rho * ((f[row-1][j+1][time_] - f[row-1][j][time_]))
				  )/D;
		density[row-1][j+1][time_] = density[row-1][j][time_] + temp_2*density[row-1][j][time_];
	}
	for(i=0;i<row-1;i++){
		double temp_1 = ( ( chi/(1+alpha*concentration[i][0][time_]  ) )
						 * ((concentration[i+1][0][time_] - concentration[i][0][time_])) 
						 + rho * ((f[i+1][0][time_] - f[i][0][time_]))
				  )/D;
		density[i+1][0][time_] = density[i][0][time_] + temp_1*density[i][0][time_];
		double temp_2 = ( ( chi/(1+alpha*concentration[i][col-1][time_]  ) )
						 * ((concentration[i+1][col-1][time_] - concentration[i][col-1][time_])) 
						 + rho * ((f[i+1][col-1][time_] - f[i][col-1][time_]))
				  )/D;
		density[i+1][col-1][time_] = density[i][col-1][time_] + temp_2*density[i][col-1][time_];
	}
}*/

void update_concentration( int i,int j ) {
	//forward difference method to calculate concentration at (n+1)th time level
	 concentration[i][j][time_] = -dt*neta*(density[i][j][time_]*concentration[i][j][time_-1]) + concentration[i][j][time_-1];
}


void update_f( int i, int j){
	//forward difference method to calculate "f" at (n+1)th time level
	f[i][j][time_] = f[i][j][time_-1] + dt*( beta*density[i][j][time_] - gamma*density[i][j][time_]*f[i][j][time_-1] );
}



void update_95_pts(int i,int j ){
	if((j==0||j==col-1) && i!=row-1){
		double temp_1 = ( ( chi/(1+alpha*concentration[i][j][time_]  ) )
						 * ((concentration[i+1][j][time_] - concentration[i][j][time_])) 
						 + rho * ((f[i+1][j][time_] - f[i][j][time_]))
				  )/D;
		nine_pts[i][j] = 0;
		Q2[i][j] = 0.0;
		Q4[i][j] = -(1 + temp_1);//0.0;
		Q5[i][j] = 1.0;
		Q6[i][j] = 0.0;
		Q8[i][j] = 0.0;
		return;
	} 
	else if( (j==0||j==col-1) && i==row-1 ) {
		double temp_1 = ( ( chi/(1+alpha*concentration[i][j][time_]  ) )
						 * ((concentration[i][j][time_] - concentration[i-1][j][time_])) 
						 + rho * ((f[i][j][time_] - f[i-1][j][time_]))
				  )/D;
		nine_pts[i][j] = 0;
		Q2[i][j] = 0.0;
		Q4[i][j] = -(1 + temp_1);
		Q5[i][j] = 1.0;
		Q6[i][j] = 0.0;
		Q8[i][j] = 0.0;
		return;
	}
	else if((i==0||i==row-1) && j!=col-1){
		double temp_1 = ( ( chi/(1+alpha*concentration[i][j][time_]  ) )
						 * ((concentration[i][j+1][time_] - concentration[i][j][time_])) 
						 + rho * ((f[i][j+1][time_] - f[i][j][time_]))
				  )/D;
		nine_pts[i][j] = 0;
		Q2[i][j] = -(1 + temp_1);//0.0;
		Q4[i][j] = 0.0;
		Q5[i][j] = 1.0;
		Q6[i][j] = 0.0;
		Q8[i][j] = 0.0;
		return;
	} 
	else if( (i==0||i==row-1) && j==col-1 ) {
		nine_pts[i][j] = 0;
		double temp_1 = ( ( chi/(1+alpha*concentration[i][j][time_]  ) )
						 * ((concentration[i][j][time_] - concentration[i][j-1][time_])) 
						 + rho * ((f[i][j][time_] - f[i][j-1][time_]))
				  )/D;
		Q2[i][j] = -(1 + temp_1);
		Q4[i][j] = 0.0;
		Q5[i][j] = 1.0;
		Q6[i][j] = 0.0;
		Q8[i][j] = 0.0;
		return;
	}

	double dho_c_ij_by_x = ( c_ij[i+1][j] - c_ij[i-1][j] ) / (2*h);
	double dho_c_ij_by_y = ( c_ij[i][j+1] - c_ij[i][j-1] ) / (2*h);
	double dho2_c_ij_by_x2 = ( c_ij[i+1][j]  -2* c_ij[i][j] + c_ij[i-1][j]) /( h*h );
	double dho2_c_ij_by_y2 = ( c_ij[i][j+1]  -2* c_ij[i][j] + c_ij[i][j-1]) /( h*h );

	double dho_d_ij_by_x = ( d_ij[i+1][j] - d_ij[i-1][j] ) / (2*h);
	double dho_d_ij_by_y = ( d_ij[i][j+1] - d_ij[i][j-1] ) / (2*h);
	double dho2_d_ij_by_x2 = ( d_ij[i+1][j]  -2* d_ij[i][j] + d_ij[i-1][j]) /( h*h );
	double dho2_d_ij_by_y2 = ( d_ij[i][j+1]  -2* d_ij[i][j] + d_ij[i][j-1]) /( h*h );

	double dho_epsilon_ij_by_x = ( epsilon_ij[i+1][j] - epsilon_ij[i-1][j] ) / (2*h);
	double dho_epsilon_ij_by_y = ( epsilon_ij[i][j+1] - epsilon_ij[i][j-1] ) / (2*h);
	double dho2_epsilon_ij_by_x2 = ( epsilon_ij[i+1][j]  -2* epsilon_ij[i][j] + epsilon_ij[i-1][j]) /( h*h );
	double dho2_epsilon_ij_by_y2 = ( epsilon_ij[i][j+1]  -2* epsilon_ij[i][j] + epsilon_ij[i][j-1]) /( h*h );

	double alpha_eqn = 1 + ( h*h*( c_ij[i][j] * c_ij[i][j] - 2 * dho_c_ij_by_x - epsilon_ij[i][j] )/12.0);
	double beta_eqn  = 1 + ( h*h*( d_ij[i][j] * d_ij[i][j] - 2 * dho_d_ij_by_y - epsilon_ij[i][j] )/12.0);
	double gamma_eqn = dho_d_ij_by_x - dho_c_ij_by_x - (c_ij[i][j] * d_ij[i][j]);

	double C_eqn_part1 = (dho2_c_ij_by_x2  - c_ij[i][j] * dho_c_ij_by_x ) * (h*h/12.0);
	double C_eqn_part2 = (dho2_c_ij_by_y2 - d_ij[i][j] * dho_c_ij_by_y ) * (h*h/12.0);
	double C_eqn_part3 = (c_ij[i][j]  - c_ij_time_n[i][j])/2.0; 												//time term
	double C_eqn_part4 =  (2*dho_epsilon_ij_by_x - c_ij[i][j] * epsilon_ij[i][j]) * (h*h/12.0);
	double C_eqn = c_ij[i][j] + C_eqn_part1 + C_eqn_part2 + C_eqn_part3 + C_eqn_part4;

	double D_eqn_part1 = (dho2_d_ij_by_x2  - c_ij[i][j] * dho_d_ij_by_x ) * (h*h/12.0);
	double D_eqn_part2 = (dho2_d_ij_by_y2 - d_ij[i][j] * dho_d_ij_by_y ) * (h*h/12.0);
	double D_eqn_part3 = (d_ij[i][j]  - d_ij_time_n[i][j])/2.0; 												//time term
	double D_eqn_part4 =  (2*dho_epsilon_ij_by_y - d_ij[i][j] * epsilon_ij[i][j]) * (h*h/12.0);
	double D_eqn = d_ij[i][j] + D_eqn_part1 + D_eqn_part2 + D_eqn_part3 + D_eqn_part4;

	double E_eqn_part1 = (dho2_epsilon_ij_by_x2  - c_ij[i][j] * dho_epsilon_ij_by_x ) * (h*h/12.0);

	double E_eqn_part2 = (dho2_epsilon_ij_by_y2  - d_ij[i][j] * dho_epsilon_ij_by_y ) * (h*h/12.0);
	double E_eqn_part3 = (epsilon_ij[i][j] - epsilon_ij_time_n[i][j])/2.0;
	double E_eqn = epsilon_ij[i][j] +E_eqn_part1 + E_eqn_part2 + E_eqn_part3;

	double G_eqn = 0.0;

	double A = 1.0 / (double)D;
	// 9 Pts. of P
	double P1 = lambda * ( -4.0 - 2.0*c_ij[i][j]*h -2.0*d_ij[i][j]*h + gamma_eqn*h*h ) / 24.00;
	double P2 = lambda * ( -12.0*beta_eqn -6.0*D_eqn*h + 4.0 + 2.0*d_ij[i][j]*h)/12.00 -(1/(double)(2.0*h*h)) -(d_ij[i][j]/(double)(4.0*h*h));
	double P3 = lambda * ( -4.0 + 2.0*c_ij[i][j]*h -2.0*d_ij[i][j]*h - gamma_eqn*h*h ) / 24.00;
	double P4 = lambda * ( -12.0*alpha_eqn -6.0*C_eqn*h + 4.0 + 2.0*c_ij[i][j]*h)/12.00 -(1/(double)(2.0*h*h)) -(c_ij[i][j]/(double)(4.0*h*h));
	double P5 = lambda * ( 2.0 * alpha_eqn + 2.0 * beta_eqn - 2.0/(double)3.0) + (2.0/(double)(h*h)) + (epsilon_ij[i][j]/2.0) - E_eqn;
	double P6 = lambda * ( -12.0*alpha_eqn +6.0*C_eqn*h + 4.0 - 2.0*c_ij[i][j]*h)/12.00 -(1/(double)(2.0*h*h)) +(c_ij[i][j]/(double)(4.0*h*h));
	double P7 = lambda * ( -4.0 - 2.0*c_ij[i][j]*h +2.0*d_ij[i][j]*h - gamma_eqn*h*h ) / 24.00;
	double P8 = lambda * ( -12.0*beta_eqn +6.0*D_eqn*h + 4.0 - 2.0*d_ij[i][j]*h)/12.00 -(1/(double)(2.0*h*h)) +(d_ij[i][j]/(double)(4.0*h*h));
	double P9 = lambda * ( -4.0 + 2.0*c_ij[i][j]*h + 2.0*d_ij[i][j]*h + gamma_eqn*h*h ) / 24.00;
	// 5 pts of Q
	Q2[i][j] = A * ( 2.0 + d_ij[i][j]*h )/2.0 -(1/(double)(2.0*h*h)) -(d_ij[i][j]/(double)(4.0*h*h));
	Q4[i][j] = A * ( 2.0 + c_ij[i][j]*h )/2.0 -(1/(double)(2.0*h*h)) -(c_ij[i][j]/(double)(4.0*h*h));
	Q5[i][j] = 8.0 * A + (2.0/(double)(h*h)) + (epsilon_ij[i][j]/2.0);
	Q6[i][j] = A * ( 2.0 - c_ij[i][j]*h )/2.0 -(1/(double)(2.0*h*h)) +(c_ij[i][j]/(double)(4.0*h*h));
	Q8[i][j] = A * ( 2.0 - d_ij[i][j]*h )/2.0 -(1/(double)(2.0*h*h)) +(d_ij[i][j]/(double)(4.0*h*h));


	nine_pts[i][j] = -12*P1*density[i-1][j-1][time_-1] + ( -12* P2 + Q2[i][j] ) * density[i][j-1][time_-1] +
					 -12*P3*density[i+1][j-1][time_-1] + ( -12* P4 + Q4[i][j] ) * density[i-1][j][time_-1] +
					 (-12*P5+Q5[i][j])*density[i][j][time_-1] + 
					 (-12*P6+Q6[i][j])*density[i+1][j][time_-1] + -12* P7 * density[i-1][j+1][time_-1] +
					 (-12*P8+Q8[i][j])*density[i][j+1][time_-1] + -12* P9 * density[i+1][j+1][time_-1];  //Summation of all 9 pts

	return;

}
double norm(double arr[][col]){
	double res = 0.0;
	int i,j;
	for(j=0;j<col;j++){
		for(i=0;i<row;i++){
			res+= fabs(arr[i][j] * arr[i][j]);
		}
	}
	return sqrt(res);
}

double vecmul(double arr1[][col], double arr2[][col]){
	double res = 0.0;
	int i,j;
	for(j=0;j<col;j++){
		for(i=0;i<row;i++){	
			res+= arr1[i][j]*arr2[i][j];
		}
	}
	return res;

}

void matmul2(double df2[][col],double df4[][col], double df5[][col],double df6[][col], double df8[][col], double t[][col], double at[][col]){
	int i,j;

	for(j=0;j<col;j++){
		at[0][j] = t[1][j]-df4[0][j]*t[0][j];
		at[row-1][j] = df4[row-1][j]*t[row-1][j]-t[row-2][j];
	}

	for(i=1;i<row-1;i++){
		at[i][0] = t[i][1]-df2[i][0]*t[i][0];
		at[i][col-1] = df2[i][col-1]*t[i][col-1]-t[i][col-2];
	}
	for(j=1;j<col-1;j++){
		for(i=1;i<row-1;i++){
			at[i][j] = df2[i][j] * t[i][j-1];
			at[i][j] += df4[i][j] * t[i-1][j];
			at[i][j] += df5[i][j] * t[i][j];
			at[i][j] += df6[i][j] * t[i+1][j];
			at[i][j] += df8[i][j] * t[i][j+1];
 		}
	}
}
void bicg(){
	
	double bn[row][col],bm[row][col];
	double at[row][col];
	int i,j;
	for(j=0;j<col;j++){
		for(i=0;i<row;i++) {
			bn[i][j]=density[i][j][time_-1];
			bm[i][j]=0.0;
		}
	}

	matmul2(Q2,Q4,Q5,Q6,Q8,bn,at);
	double res[row][col],res_old[row][col];
	double v1[row][col];double p[row][col];
	double s[row][col];
	for(j=0;j<col;j++){
		for(i=0;i<row;i++){
			res_old[i][j] = res[i][j] = 0.0;
			s[i][j] = 0.0; v1[i][j] =0.0;
			p[i][j]=0.0;
		}
	}
	for(j=0;j<col;j++){
		for(i=0;i<row;i++){
			res_old[i][j] = res[i][j] = ( nine_pts[i][j] - at[i][j] );
			v1[i][j] = 0.0, p[i][j] = 0.0;
		}
	}

	double r0=1,a1=1,w1=1;
	int iter = 0;
	while(norm(res) > err){
		r0 = -w1*r0;
		double rn = vecmul(res,res_old);
		double b1 = ((rn)*(a1))/r0;
		for(j=0;j<col;j++){
			for(i=0;i<row;i++) {
				bm[i][j] = bn[i][j];
				p[i][j] = res[i][j] - b1 * (p[i][j] - w1*v1[i][j]);
			}
		}
		 
		matmul2(Q2,Q4,Q5,Q6,Q8,p,v1);
		a1=rn/vecmul(res_old,v1);
		for(j=0;j<col;j++){
        	for(i=0;i<row;i++){
        		s[i][j] = res[i][j]- a1*v1[i][j];
        	}
        }
 		if(norm(s)<err) {
 			for(j=0;j<col;j++){
        		for(i=0;i<row;i++){
        			bn[i][j] = bn[i][j]+ a1*p[i][j];
        		}
        	}

    break;
 		}
       r0=rn;
       matmul2(Q2,Q4,Q5,Q6,Q8,s,at);
       w1=vecmul(at,s)/vecmul(at,at);

       if(rn==0)
       		break;
       double calc = 0.0;
       for(j=0;j<col;j++){
       		for(i=0;i<row;i++){
          		bn[i][j] += a1 * p[i][j] + w1*s[i][j];
          		calc += (a1 * p[i][j] + w1*s[i][j]) * (a1 * p[i][j] + w1*s[i][j]);
          		res[i][j]=s[i][j] - w1 * at[i][j];
          		bn[i][j]=bm[i][j] + 0.005*(bn[i][j]-bm[i][j]);
        }}
		calc = sqrt(calc);
    	iter++;
    	if(fabs(calc/norm(bn))<=err){
			break;
		}
   if(iter>10000)  break;
	}
	for(j=0;j<col;j++){
		for(i=0;i<row;i++){
			density[i][j][time_] =(bn[i][j]);
		}
	}

}

void output(FILE* f_) {
	int i,j,k1;
	
    	for(j=0;j<col;j++){
		for(i=0;i<row;i++){
			fprintf(f_,"%lf\t%lf\t%lf\t%lf\t%lf\n",i*h,j*h,concentration[i][j][time_-1],
			f[i][j][time_-1],density[i][j][time_-1]);
		}
	}
}

double max(){
	double max_diff = 0.0;
	int i,j;
	for(i=0; i<row; i++){
		for(j=0; j<col; j++){
			if(fabs(density[i][j][time_] - density[i][j][time_-1]) > max_diff)
				max_diff = fabs(density[i][j][time_] - density[i][j][time_-1]);
		}
	}
	return max_diff;
}

int main() {
	int l = dt*(row-1)*(row-1);
	int tol=0.5E-4; 
	int i,j; 	//Loop vars
	int iter_ = 0;
	intialize();
	initiate_boundary_conditions();
	FILE* f_ = fopen("cts_10iter_0.001.dat","w");
	fprintf(f_,"title=\"\"\n");
	fprintf(f_,"variables=x,y,c,f,n\n");
	fprintf(f_,"zone T=\"\",i=%d,j=%d\n",row,col);
	double max_diff = 0.0;
	do {
		for(i=0;i<row;i++){
			for(j=0;j<col;j++){
				concentration[i][j][time_] = concentration[i][j][time_-1];
				f[i][j][time_] = f[i][j][time_-1];
			}
		}

		//ij boundary condns
		for(i=0;i<row;i++) {
			for(j=0;j<col;j++){	
					update_95_pts(i,j);
			}
		}

		bicg();
		for(i=0;i<row;i++){
			for(j=0;j<col;j++){
				update_concentration(i,j);
				update_f(i,j);
			}
		}
		
		
		max_diff = max();
		for(i=0;i<row;i++){
			for(j=0;j<col;j++){
				density[i][j][time_-1] = density[i][j][time_];
				density[i][j][time_] = 0.0;
				concentration[i][j][time_-1] = concentration[i][j][time_];
				concentration[i][j][time_] = 0.0;
				f[i][j][time_-1] = f[i][j][time_];
				f[i][j][time_] = 0.0;
				c_ij[i][j] = 0.0;
				d_ij[i][j] = 0.0;
				epsilon_ij[i][j] = 0.0;
				c_ij_time_n[i][j] = 0.0;
				d_ij_time_n[i][j] = 0.0;
				epsilon_ij_time_n[i][j] = 0.0;
			}
		}
		//printf("At %d iteration Value of density at %d,%d is %lf \n",iter_,row/2-1,col/2-1,density[row/2-1][col/2-1][time_]);
		
	} while(iter_++ < 1000);
	output(f_);
	fclose(f_);
	return 0;
}
