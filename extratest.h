// ==== from CIM paper
inline double getfPaper(double *x, double thesign, PBData &pb){
	double value;
	
	if(thesign > 0){
		value = -pb.epsilonp*(2*x[0] + 12*pow(x[0],2) + 12*pow(x[1],2) - 4*cos(2*x[0] + pow(x[1],2) + pow(x[2],3)) - 4*pow(x[1],2)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3)) - 9*pow(x[2],4)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3)) - 2*sin(2*x[0] + pow(x[1],2) + pow(x[2],3)) - 6*x[2]*sin(2*x[0] + pow(x[1],2) + pow(x[2],3)));
	} else{
		value = -pb.epsilonm*(8*x[0] + 6*x[1] + 12*pow(x[2],2) + 12*cos(3*(pow(x[0],2) + pow(x[1],2))) - 36*pow(x[0],2)*sin(3*(pow(x[0],2) + pow(x[1],2))) - 36*pow(x[1],2)*sin(3*(pow(x[0],2) + pow(x[1],2))));
	}
	return value;
}

inline double getuPaper(double *x, double thesign){
	double value;
	if(thesign > 0){
		value = pow(x[0],4) + x[0]*x[1] + pow(x[1],4) + x[0]*pow(x[2],2) + cos(2*x[0] \
							+ pow(x[1],2) + pow(x[2],3));
	} else {
		value = pow(x[0],3) + x[0]*pow(x[1],2) + pow(x[1],3) + pow(x[2],4) + \
						sin(3*(pow(x[0],2) + pow(x[1],2)));
	}
	return value;
}


inline double getDuPaper(double *x, int s, double thesign){
	double value;
	if(thesign>0){
		if(s==0)
			value = 4*pow(x[0],3) + x[1] + pow(x[2],2) - 2*sin(2*x[0] + pow(x[1],2) + pow(x[2],3));
		if(s==1)
			value = x[0] + 4*pow(x[1],3) - 2*x[1]*sin(2*x[0] + pow(x[1],2) + pow(x[2],3));
		if(s==2)
			value = 2*x[0]*x[2] - 3*pow(x[2],2)*sin(2*x[0] + pow(x[1],2) + pow(x[2],3));
	}else{
		if(s==0)
			value = 3*pow(x[0],2) + pow(x[1],2) + 6*x[0]*cos(3*(pow(x[0],2) + pow(x[1],2)));
		if(s==1)
			value = 2*x[0]*x[1] + 3*pow(x[1],2) + 6*x[1]*cos(3*(pow(x[0],2) + pow(x[1],2)));
		if(s==2)
			value = 4*pow(x[2],3);
	}
	return value;
}

inline double getD2uPaper(double *x, int r, int s, double thesign){
	double value;
	if (thesign>0){
		if(r==0&& s == 0)
			value = 12*pow(x[0],2) - 4*cos(2*x[0] + pow(x[1],2) + pow(x[2],3));
		if(r==0&& s == 1)
			value = 1 - 4*x[1]*cos(2*x[0] + pow(x[1],2) + pow(x[2],3));
		if(r==0&& s == 2)
			value = 2*x[2] - 6*pow(x[2],2)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3));
		if(r==1&& s == 0)
			value = 1 - 4*x[1]*cos(2*x[0] + pow(x[1],2) + pow(x[2],3));
		if(r==1&& s == 1)
			value = 12*pow(x[1],2) - 4*pow(x[1],2)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3)) - 2*sin(2*x[0] + pow(x[1],2) + pow(x[2],3));
		if(r==1&& s == 2)
			value = -6*x[1]*pow(x[2],2)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3));
		if(r==2&& s == 0)
			value = 2*x[2] - 6*pow(x[2],2)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3));
		if(r==2&& s == 1)
			value = -6*x[1]*pow(x[2],2)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3));
		if(r==2&& s == 2)
			value = 2*x[0] - 9*pow(x[2],4)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3)) - 6*x[2]*sin(2*x[0] + pow(x[1],2) + pow(x[2],3));
	}else{
		if(r==0&& s == 0)
			value = 6*x[0] + 6*cos(3*(pow(x[0],2) + pow(x[1],2))) - 36*pow(x[0],2)*sin(3*(pow(x[0],2) + pow(x[1],2)));
		if(r==0&& s == 1)
			value = 2*x[1] - 36*x[0]*x[1]*sin(3*(pow(x[0],2) + pow(x[1],2)));
		if(r==0&& s == 2)
			value = 0;
		if(r==1&& s == 0)
			value = 2*x[1] - 36*x[0]*x[1]*sin(3*(pow(x[0],2) + pow(x[1],2)));
		if(r==1&& s == 1)
			value = 2*x[0] + 6*x[1] + 6*cos(3*(pow(x[0],2) + pow(x[1],2))) - 36*pow(x[1],2)*sin(3*(pow(x[0],2) + pow(x[1],2)));
		if(r==1&& s == 2)
			value = 0;
		if(r==2&& s == 0)
			value = 0;
		if(r==2&& s == 1)
			value = 0;
		if(r==2&& s == 2)
			value = 12*pow(x[2],2);
	}
	return value;
}


inline double gettauPaper(double *x){
	return -pow(x[0],3) + pow(x[0],4) + x[0]*x[1] - x[0]*pow(x[1],2) - \
pow(x[1],3) + pow(x[1],4) + x[0]*pow(x[2],2) - pow(x[2],4) + \
cos(2*x[0] + pow(x[1],2) + pow(x[2],3)) - sin(3*(pow(x[0],2) + \
pow(x[1],2)));
}


inline void getDtauPaper(double *Dtau, double *x){
	Dtau[0] = -3*pow(x[0],2) + 4*pow(x[0],3) + x[1] - pow(x[1],2) + pow(x[2],2) - 6*x[0]*cos(3*(pow(x[0],2) + pow(x[1],2))) - 2*sin(2*x[0] + pow(x[1],2) + pow(x[2],3));
	Dtau[1] = x[0] - 2*x[0]*x[1] - 3*pow(x[1],2) + 4*pow(x[1],3) - 6*x[1]*cos(3*(pow(x[0],2) + pow(x[1],2))) - 2*x[1]*sin(2*x[0] + pow(x[1],2) + pow(x[2],3));
	Dtau[2] = 2*x[0]*x[2] - 4*pow(x[2],3) - 3*pow(x[2],2)*sin(2*x[0] + pow(x[1],2) + pow(x[2],3));
}


inline void getD2tauPaper(double **D2tau, double *x){
	D2tau[0][0] = -6*x[0] + 12*pow(x[0],2) - 6*cos(3*(pow(x[0],2) + pow(x[1],2))) - 4*cos(2*x[0] + pow(x[1],2) + pow(x[2],3)) + 36*pow(x[0],2)*sin(3*(pow(x[0],2) + pow(x[1],2)));
	D2tau[0][1] = 1 - 2*x[1] - 4*x[1]*cos(2*x[0] + pow(x[1],2) + pow(x[2],3)) + 36*x[0]*x[1]*sin(3*(pow(x[0],2) + pow(x[1],2)));
	D2tau[0][2] = 2*x[2] - 6*pow(x[2],2)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3));
	D2tau[1][0] = 1 - 2*x[1] - 4*x[1]*cos(2*x[0] + pow(x[1],2) + pow(x[2],3)) + 36*x[0]*x[1]*sin(3*(pow(x[0],2) + pow(x[1],2)));
	D2tau[1][1] = -2*x[0] - 6*x[1] + 12*pow(x[1],2) - 6*cos(3*(pow(x[0],2) + pow(x[1],2))) - 4*pow(x[1],2)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3)) + 36*pow(x[1],2)*sin(3*(pow(x[0],2) + pow(x[1],2))) - 2*sin(2*x[0] + pow(x[1],2) + pow(x[2],3));
	D2tau[1][2] = -6*x[1]*pow(x[2],2)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3));
	D2tau[2][0] = 2*x[2] - 6*pow(x[2],2)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3));
	D2tau[2][1] = -6*x[1]*pow(x[2],2)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3));
	D2tau[2][2] = 2*x[0] - 12*pow(x[2],2) - 9*pow(x[2],4)*cos(2*x[0] + pow(x[1],2) + pow(x[2],3)) - 6*x[2]*sin(2*x[0] + pow(x[1],2) + pow(x[2],3));
}

//dot( epsp * Dup - epsm * Dum,normal)
inline double getsigmaPaper(double *x, double *normal, PBData &pb){
	double sigma = 0;
	for(int s = 0; s<3; s++){
		sigma += (pb.epsilonp * getDuPaper(x,s,1) - pb.epsilonm * getDuPaper(x,s,-1)) * normal[s]; 
	}
	return sigma;
}

inline void getDsigmaPaper(double *Dsigma, double *x, double *normal, double **Dnormal, PBData &pb){
		double gradp[3], gradm[3], D2p[3][3],D2m[3][3];
		for (int r = 0; r < 3; r++)
		{
			 gradp[r] = getDuPaper(x,r,1);//grad(u+)
			 gradm[r] = getDuPaper(x,r,-1);//grad(u-)
			 for (int s = 0; s < 3; s++)
			 {
					D2p[r][s] = getD2uPaper(x, r, s, 1);
					D2m[r][s] = getD2uPaper(x, r, s, -1);
			 }
		}
		for (int r = 0; r < 3; r++)
		{
			 Dsigma[r] = 0.0;
			 for (int s = 0; s < 3; s++)
					Dsigma[r] += (pb.epsilonp*gradp[s]-pb.epsilonm*gradm[s])*Dnormal[s][r]+
											 (pb.epsilonp*D2p[r][s]-pb.epsilonm*D2m[r][s])*normal[s];
		}

}

// === Torus from Li-tien
inline double getfTorus(double *x, double thesign, PBData &pb){
	double value;
	
	if(thesign > 0){
		value = 3*pb.epsilonp*cos(x[0])*cos(x[1])*cos(x[2]);

	} else{
		value = 3*pb.epsilonm*sin(x[0])*sin(x[1])*sin(x[2]);
	}	
	
	return value;
}

inline double getuTorus(double *x, double thesign){
	double value;
	if(thesign > 0){
		value = cos(x[0]) * cos(x[1]) * cos(x[2]);
	} else{
		value = sin(x[0]) * sin(x[1]) * sin(x[2]);
	}
	return value;
}

inline double getDuTorus(double *x, int s, double thesign){
	double value;
	if(thesign>0){
		if(s==0)
			value = -(cos(x[1])*cos(x[2])*sin(x[0]));
		if(s==1)
			value = -(cos(x[0])*cos(x[2])*sin(x[1]));
		if(s==2)
			value = -(cos(x[0])*cos(x[1])*sin(x[2]));
	}else{
		if(s==0)
			value = cos(x[0])*sin(x[1])*sin(x[2]);
		if(s==1)
			value = cos(x[1])*sin(x[0])*sin(x[2]);
		if(s==2)
			value = cos(x[2])*sin(x[0])*sin(x[1]);
	}
	return value;
}

inline double getD2uTorus(double *x, int r, int s, double thesign){
	double value;
	if (thesign>0){
		if(r==0&& s == 0)
			value = -(cos(x[0])*cos(x[1])*cos(x[2]));
		if(r==0&& s == 1)
			value = cos(x[2])*sin(x[0])*sin(x[1]);
		if(r==0&& s == 2)
			value = cos(x[1])*sin(x[0])*sin(x[2]);
		if(r==1&& s == 0)
			value = cos(x[2])*sin(x[0])*sin(x[1]);
		if(r==1&& s == 1)
			value = -(cos(x[0])*cos(x[1])*cos(x[2]));
		if(r==1&& s == 2)
			value = cos(x[0])*sin(x[1])*sin(x[2]);
		if(r==2&& s == 0)
			value = cos(x[1])*sin(x[0])*sin(x[2]);
		if(r==2&& s == 1)
			value = cos(x[0])*sin(x[1])*sin(x[2]);
		if(r==2&& s == 2)
			value = -(cos(x[0])*cos(x[1])*cos(x[2]));
	}else{
		if(r==0&& s == 0)
			value = -(sin(x[0])*sin(x[1])*sin(x[2]));
		if(r==0&& s == 1)
			value = cos(x[0])*cos(x[1])*sin(x[2]);
		if(r==0&& s == 2)
			value = cos(x[0])*cos(x[2])*sin(x[1]);
		if(r==1&& s == 0)
			value = cos(x[0])*cos(x[1])*sin(x[2]);
		if(r==1&& s == 1)
			value = -(sin(x[0])*sin(x[1])*sin(x[2]));
		if(r==1&& s == 2)
			value = cos(x[1])*cos(x[2])*sin(x[0]);
		if(r==2&& s == 0)
			value = cos(x[0])*cos(x[2])*sin(x[1]);
		if(r==2&& s == 1)
			value = cos(x[1])*cos(x[2])*sin(x[0]);
		if(r==2&& s == 2)
			value = -(sin(x[0])*sin(x[1])*sin(x[2]));
		}
	return value;
}

inline double gettauTorus(double *x){
	return cos(x[0])*cos(x[1])*cos(x[2]) - sin(x[0])*sin(x[1])*sin(x[2]);
}

inline void getDtauTorus(double *Dtau, double *x){
	 Dtau[0] = -(cos(x[1])*cos(x[2])*sin(x[0])) - cos(x[0])*sin(x[1])*sin(x[2]);

	 Dtau[1] = -(cos(x[0])*cos(x[2])*sin(x[1])) - cos(x[1])*sin(x[0])*sin(x[2]);

	 Dtau[2] = -(cos(x[2])*sin(x[0])*sin(x[1])) - cos(x[0])*cos(x[1])*sin(x[2]);
}


inline void getD2tauTorus(double **D2tau, double *x){
	D2tau[0][0] = -(cos(x[0])*cos(x[1])*cos(x[2])) + sin(x[0])*sin(x[1])*sin(x[2]);
	D2tau[0][1] = cos(x[2])*sin(x[0])*sin(x[1]) - cos(x[0])*cos(x[1])*sin(x[2]);
	D2tau[0][2] = -(cos(x[0])*cos(x[2])*sin(x[1])) + cos(x[1])*sin(x[0])*sin(x[2]);
	D2tau[1][0] = cos(x[2])*sin(x[0])*sin(x[1]) - cos(x[0])*cos(x[1])*sin(x[2]);
	D2tau[1][1] = -(cos(x[0])*cos(x[1])*cos(x[2])) + sin(x[0])*sin(x[1])*sin(x[2]);
	D2tau[1][2] = -(cos(x[1])*cos(x[2])*sin(x[0])) + cos(x[0])*sin(x[1])*sin(x[2]);
	D2tau[2][0] = -(cos(x[0])*cos(x[2])*sin(x[1])) + cos(x[1])*sin(x[0])*sin(x[2]);
	D2tau[2][1] = -(cos(x[1])*cos(x[2])*sin(x[0])) + cos(x[0])*sin(x[1])*sin(x[2]);
	D2tau[2][2] = -(cos(x[0])*cos(x[1])*cos(x[2])) + sin(x[0])*sin(x[1])*sin(x[2]);
}

inline double getsigmaTorus(double *x, double *normal, PBData &pb){
	double sigma = 0;
	for(int s = 0; s<3; s++){
		sigma += (pb.epsilonp * getDuTorus(x,s,1) - pb.epsilonm * getDuTorus(x,s,-1)) * normal[s]; 
	}
	return sigma;
}


inline void getDsigmaTorus(double *Dsigma, double *x, double *normal, double **Dnormal, PBData &pb){
		double gradp[3], gradm[3], D2p[3][3],D2m[3][3];
		for (int r = 0; r < 3; r++)
		{
			 gradp[r] = getDuTorus(x,r,1);//grad(u+)
			 gradm[r] = getDuTorus(x,r,-1);//grad(u-)
			 for (int s = 0; s < 3; s++)
			 {
					D2p[r][s] = getD2uTorus(x, r, s, 1);
					D2m[r][s] = getD2uTorus(x, r, s, -1);
			 }
		}
		for (int r = 0; r < 3; r++)
		{
			 Dsigma[r] = 0;
			 for (int s = 0; s < 3; s++)
					Dsigma[r] += (pb.epsilonp*gradp[s]-pb.epsilonm*gradm[s])*Dnormal[s][r]+
											 (pb.epsilonp*D2p[r][s]-pb.epsilonm*D2m[r][s])*normal[s];
      }
}

// === simple function 1
// #define SIMPLE1
#ifdef SIMPLE1
// up(x,y,z) = x^2, um(x,y,z) = -x^2
inline double getuSimple(double *x, double thesign){
	return x[0] * x[0] * ((thesign>0)?1:-1);
}

// lap up = 2, lap(um) = -2;
inline double getfSimple(double *x, double thesign, PBData &pb){
	if (thesign>0){
		return - pb.epsilonp * 2;
	}else{
		return - pb.epsilonm * (-2);
	}

}

// Dup(x,y,z) = [2x,0,0], Dum(x,y,z) = [-2x,0,0]
inline double getDuSimple(double *x, int s, double thesign){
	if (s != 0){
		return 0;
	}else{
		return 2*x[0] * ((thesign>0)?1:-1);	
	}
	
}
// D2up(x,y,z)
inline double getD2uSimple(double *x, int r, int s, double thesign){
	if(r == 0 && s == 0){
		return 2  * ((thesign>0)?1:-1);
	}else{
		return 0;
	}
}

// tau = 2x^2
inline double gettauSimple(double *x){
	return 2 * x[0] * x[0];
}

inline void getDtauSimple(double *Dtau, double *x){
	Dtau[0] = 4*x[0];
	Dtau[1] = 0;
	Dtau[2] = 0;
}


inline void getD2tauSimple(double **D2tau, double *x){
	D2tau[0][0] = 4.0;
	D2tau[0][1] = 0.0;
	D2tau[0][2] = 0.0;
	D2tau[1][0] = 0.0;
	D2tau[1][1] = 0.0;
	D2tau[1][2] = 0.0;
	D2tau[2][0] = 0.0;
	D2tau[2][1] = 0.0;
	D2tau[2][2] = 0.0;	
}

#else // use simple2

// simple function 2 : x^2 + y^2 + z^2 + x y + xz + yz (sign)
inline double getuSimple(double *x, double thesign){
	double sign =  (thesign>0)? 1.0:-1.0;
	return sign*(pow(x[0],2) + x[0]*x[1] + pow(x[1],2) + x[0]*x[2] + x[1]*x[2] + pow(x[2],2)) ;
}

inline double gettauSimple(double *x){
	return 2*pow(x[0],2) + 2*x[0]*x[1] + 2*pow(x[1],2) + 2*x[0]*x[2] + 2*x[1]*x[2] + 2*pow(x[2],2);
}

inline double getfSimple(double *x, double thesign, PBData &pb){
	if (thesign>0){
		return  - 6 * pb.epsilonp;
	}else{
		return  6 * pb.epsilonm;
	}

}

// Dup(x,y,z) = [2 x y^2,2 y x^2,0], Dum(x,y,z) = - Dup
inline double getDuSimple(double *x, int s, double thesign){
	
	if (thesign>0){
		if(s==0)
			return 2*x[0] + x[1] + x[2];
		if(s==1)
			return x[0] + 2*x[1] + x[2];
		if(s==2)
			return x[0] + x[1] + 2*x[2];

	}else{
		if(s==0)
			return -2*x[0] - x[1] - x[2];
		if(s==1)
			return -x[0] - 2*x[1] - x[2];
		if(s==2)
			return -x[0] - x[1] - 2*x[2];
	}
}

// D2up(x,y,z) 
inline double getD2uSimple(double *x, int r, int s, double thesign){
	double d2u[3][3] = {{2,1,1},{1,2,1},{1,1,2}};
	return ((thesign>0)?1.0:-1.0)*d2u[r][s];
}

inline void getDtauSimple(double *Dtau, double *x){
	Dtau[0] = 4*x[0] + 2*x[1] + 2*x[2];
	Dtau[1] = 2*x[0] + 4*x[1] + 2*x[2];
	Dtau[2] = 2*x[0] + 2*x[1] + 4*x[2];
}


inline void getD2tauSimple(double **D2tau, double *x){
	D2tau[0][0] = 4;
	D2tau[0][1] = 2;
	D2tau[0][2] = 2;
	D2tau[1][0] = 2;
	D2tau[1][1] = 4;
	D2tau[1][2] = 2;
	D2tau[2][0] = 2;
	D2tau[2][1] = 2;
	D2tau[2][2] = 4;
}

#endif


// this two is defined with getDuSimple and getD2uSimple, no need to hard code
inline double getsigmaSimple(double *x, double *normal, PBData &pb){
	double sigma = 0;
	for(int s = 0; s<3; s++){
		sigma += (pb.epsilonp * getDuSimple(x,s,1) - pb.epsilonm * getDuSimple(x,s,-1)) * normal[s]; 
	}
	return sigma;
}


inline void getDsigmaSimple(double *Dsigma, double *x, double *normal, double **Dnormal, PBData &pb){
	double gradp[3], gradm[3], D2p[3][3],D2m[3][3];
	for (int r = 0; r < 3; r++)
	{
		 gradp[r] = getDuSimple(x,r,1);//grad(u+)
		 gradm[r] = getDuSimple(x,r,-1);//grad(u-)
		 for (int s = 0; s < 3; s++)
		 {
				D2p[r][s] = getD2uSimple(x, r, s, 1);
				D2m[r][s] = getD2uSimple(x, r, s, -1);
		 }
	}
	for (int r = 0; r < 3; r++)
	{
		 Dsigma[r] = 0;
		 for (int s = 0; s < 3; s++)
				Dsigma[r] += (pb.epsilonp*gradp[s]-pb.epsilonm*gradm[s])*Dnormal[s][r]+
										 (pb.epsilonp*D2p[r][s]-pb.epsilonm*D2m[r][s])*normal[s];
	}
}



// === interface to tryvn.C
inline double getfTest(double *x, double thesign, PBData &pb){
	if (globtestnum == 11){
		return getfPaper(x, thesign, pb);
	}
	if (globtestnum == 10){
		return getfTorus(x, thesign, pb);
	}

	if (globtestnum == 12){
		return getfSimple(x, thesign, pb);
	}
}

inline double getuTest(double *x, double thesign){
	if (globtestnum == 11){
		return getuPaper(x, thesign);
	}
	if (globtestnum == 10){
		return getuTorus(x, thesign);
	}

	if (globtestnum == 12){
		return getuSimple(x, thesign);
	}
}


inline double getDuTest(double *x, int s, double thesign){
	if (globtestnum == 11){
		return getDuPaper(x, s, thesign);
	}
	if (globtestnum == 10){
		return getDuTorus(x, s, thesign);
	}

	if (globtestnum == 12){
		return getDuSimple(x, s, thesign);
	}
}

inline double getD2uTest(double *x, int r, int s, double thesign){
	if (globtestnum == 11){
		return getD2uPaper(x, r, s, thesign);
	}
	if (globtestnum == 10){
		return getD2uTorus(x, r, s, thesign);
	}

	if (globtestnum == 12){
		return getD2uSimple(x, r, s, thesign);
	}
}


inline double gettauTest(double *x){
	if (globtestnum == 11){
		return gettauPaper(x);
	}
	if (globtestnum == 10){
		return gettauTorus(x);
	}
	if (globtestnum == 12){
		return gettauSimple(x);
	}
}


inline void getDtauTest(double *Dtau, double *x){
	if (globtestnum == 11){
		getDtauPaper(Dtau, x);
	}
	if (globtestnum == 10){
		getDtauTorus(Dtau, x);
	}
	if (globtestnum == 12){
		getDtauSimple(Dtau, x);
	}
}


inline void getD2tauTest(double **D2tau, double *x){
	if (globtestnum == 11){
			getD2tauPaper(D2tau, x);
	}
	if (globtestnum == 10){
			getD2tauTorus(D2tau, x);
	}
	if (globtestnum == 12){
			getD2tauSimple(D2tau, x);
	}
}

//dot( epsp * Dup - epsm * Dum,normal)
inline double getsigmaTest(double *x, double *normal, PBData &pb){
	if (globtestnum == 11){
			return getsigmaPaper(x, normal, pb);
	}
	if (globtestnum == 10){
			return getsigmaTorus(x, normal, pb);
	}
	if (globtestnum == 12){
			return getsigmaSimple(x, normal, pb);
	}
}

inline void getDsigmaTest(double *Dsigma, double *x, double *normal, double **Dnormal, PBData &pb){
	if (globtestnum == 11){
		getDsigmaPaper(Dsigma, x, normal, Dnormal, pb);
	}
	if (globtestnum == 10){
		getDsigmaTorus(Dsigma, x, normal, Dnormal, pb);
	}
	if (globtestnum == 12){
		getDsigmaSimple(Dsigma, x, normal, Dnormal, pb);
	}
}

//get exact normal in direction r
double getexactnormalTest(double *x, int r){
	double value;
	if(SURFOPT == 11){//ellipsoid
		if(r==0)
			value = (4*x[0])/sqrt(16*pow(x[0],2) + 36*pow(x[1],2) + 144*pow(x[2],2));
		if(r==1)
			value = (6*x[1])/sqrt(16*pow(x[0],2) + 36*pow(x[1],2) + 144*pow(x[2],2));
		if(r==2)
			value = (12*x[2])/sqrt(16*pow(x[0],2) + 36*pow(x[1],2) + 144*pow(x[2],2));
	}else if(SURFOPT == 10){
			double temp = 0.0;
			for (int m = 0; m < 3; m++)
				temp += x[m]*x[m];
			temp = sqrt(temp);
			value = x[r]/temp;

	}else{
		cerr<<"exact normal not defined"<<endl;
		exit(1);
	}
	return value;
}
