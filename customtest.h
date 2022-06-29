// ==== from CIM paper

double dbc_custom(double *x, int thedim, double thetime){
	double value = 0.0;

	for (int j = 0; j < thedim; j++)
		value += x[j]*x[j];
	value = -1.0/(1.0+value);

	return value;
}

double getf_custom(double *x, double thesign, PBData &pb, GridData &grid){
	double value = 0.0;
	double radius2 = 0.0;
	for (int r = 0; r < grid.dim; r++)
		radius2 += x[r]*x[r];
	if (thesign < 0.0)
	{
		 for (int r = 0; r < grid.dim; r++)
			value -= pb.epsilonm*
					 (8.0*x[r]*x[r]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2))-
					  2.0/((1.0+radius2)*(1.0+radius2)));
	}
	else
	{
	 for (int r = 0; r < grid.dim; r++)
		value -= pb.epsilonp*
				 (-8.0*x[r]*x[r]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2))+
				  2.0/((1.0+radius2)*(1.0+radius2)));
	}
	return value;
}


double gettau_custom(double *x, GridData &grid){
	int r;
	double valuep, valuem, radius2;

	radius2 = 0.0;
	for (int r = 0; r < grid.dim; r++)
		radius2 += x[r]*x[r];

	valuep = -1.0/(1.0+radius2);
	valuem = 1.0/(1.0+radius2);
	return valuep-valuem;
}


void getDtau_custom(double *Dtau, double *x, GridData &grid){
	double valuep, valuem, radius2;

	radius2 = 0.0;
	for (int r = 0; r < grid.dim; r++)
		radius2 += x[r]*x[r];

	for (int s = 0; s < grid.dim; s++)
	{
		valuep = 2.0*x[s]/((1.0+radius2)*(1.0+radius2));
		valuem = -2.0*x[s]/((1.0+radius2)*(1.0+radius2));
		Dtau[s] = valuep-valuem;
	}
}


void getD2tau_custom(double **D2tau, double *x, GridData &grid){
	double valuep, valuem, radius2;

	radius2 = 0.0;
	for (int r = 0; r < grid.dim; r++)
		radius2 += x[r]*x[r];
	for (int r = 0; r < grid.dim; r++)
		for (int s = r; s < grid.dim; s++)
		{
			valuep = -8.0*x[r]*x[s]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2));
			valuem = 8.0*x[r]*x[s]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2));
			if (s == r)
			{
				valuep += 2.0/((1.0+radius2)*(1.0+radius2));
				valuem += -2.0/((1.0+radius2)*(1.0+radius2));
			}
			D2tau[r][s] = valuep-valuem;
			D2tau[s][r] = D2tau[r][s];
		}
}

double getsigma_custom(double *x, double *normal, PBData &pb, GridData &grid){
	int r, s;
	double sigma, gradp[grid.dim], gradm[grid.dim], radius2;

	sigma = 0.0;
	radius2 = 0.0;
	for (r = 0; r < grid.dim; r++)
		radius2 += x[r]*x[r];
	for (r = 0; r < grid.dim; r++)
	{
		gradp[r] = 2.0*x[r]/((1.0+radius2)*(1.0+radius2));
		gradm[r] = -2.0*x[r]/((1.0+radius2)*(1.0+radius2));
		sigma += (pb.epsilonp*gradp[r]-pb.epsilonm*gradm[r])*normal[r];
	}

	return sigma;
}

void getDsigma_custom(double *Dsigma, double *x, double *normal, double **Dnormal, PBData &pb, GridData &grid)
{
	int r, s;
	double gradp[grid.dim], gradm[grid.dim], D2p[grid.dim][grid.dim],
             D2m[grid.dim][grid.dim], radius2;

      radius2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];
      for (r = 0; r < grid.dim; r++)
      {
         gradp[r] = 2.0*x[r]/((1.0+radius2)*(1.0+radius2));//grad(u+)
         gradm[r] = -2.0*x[r]/((1.0+radius2)*(1.0+radius2));//grad(u-)
         for (s = 0; s < grid.dim; s++)
         {
            D2p[r][s] = -8.0*x[r]*x[s]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2));
            D2m[r][s] = 8.0*x[r]*x[s]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2));
            if (s == r)
            {
               D2p[r][s] += 2.0/((1.0+radius2)*(1.0+radius2));
               D2m[r][s] += -2.0/((1.0+radius2)*(1.0+radius2));
            }
         }
      }
      for (r = 0; r < grid.dim; r++)
      {
         Dsigma[r] = 0.0;
         for (s = 0; s < grid.dim; s++)
            Dsigma[r] += (pb.epsilonp*gradp[s]-pb.epsilonm*gradm[s])*Dnormal[s][r]+
                         (pb.epsilonp*D2p[r][s]-pb.epsilonm*D2m[r][s])*normal[s];
      }

}

