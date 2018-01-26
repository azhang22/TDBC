//++++++++++++++++++++++++++filename: air.cpp +++++++++++++++++++++++++++++++++//

//-----------------------air acoustic wave equation--------------------------//
/*******************************************************************************/
/* the program just apply Navier-Stokes equation to calculate acoustic pressure
/* when sound propagation with Eulerian time-domain model
/*******************************************************************************/

//-----------------------scheme from Victor W.Sparrow--------------------------//
/*******************************************************************************/
/* Victor W.Sparrow
/* calculate pressure first, and then calculate velocity
/* apply non-staggered grid(colocated), velocity and pressure are at the same grid point
/* point pieces at y direction(air): 0,1,...IMAX-1;
/* point pieces at z direction(air): 0,1,...JMAX-1;
/********************************************************************************/

//---------------------------scheme from Erik M.Salomons-----------------------//
/*******************************************************************************/
/* Eulerian Time-Domain Model for sound propagation over
/* a finite impedance ground surface.
/* Comparison with frequency-domain models.
/* Author: Erik M.Salomons.
/* ACTA Vol.88(2002) 483-492
/* calculate velocity first, and then calculate pressure
/* apply staggered grid, velocity and pressure are at the different grid point
/* point pieces at y direction(air): 0,1,...IMAX-1;
/* point pieces at z direction(air): 0,1,...JMAX-1;
/* east and west boundary point of v is v[1][j] and v[IMAX-1][j]
/* north and south boundary point of w is w[i][1] and w[i][JMAX-1]
/* the four boundary point v[1][j], v[IMAX-1][j], w[i][1], w[i][JMAX-1]
/* are calculated through equation, not from interpolation. This is different from
/* above colocated scheme
/*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "air.h"
#include "mygeneral.h"

SpeedProfile::SpeedProfile(double velocity_coef,int VM) //double velocity_coef,double sound_speed,int VM
{

	alpha=1.256431;
	b=velocity_coef;
	velocity_method=VM;
	/*
	c=sound_speed;velocity_method=VM;
	switch(velocity_method){
	case 2:
		{
			y0=0;z0=0;L=0.03;
			vorticity=b*(1+2*alpha)*c/L;
			circulation=vorticity*L*L*PI/alpha;
		}
		break;
	case 3:
		{
			y0=0.66;z0=0;L=0.03;
			circulation=b;
		}
		break;
	case 6:
		{
			y0=3.5;z0=4;L=1;
			vorticity=b*(1+2*alpha)*c/L;
			circulation=vorticity*L*L*PI/alpha;
		}
		break;
	default:
		;
	}
	*/
  }

double SpeedProfile::Vav(double x,double y,double z)
{
	switch(velocity_method){
	case 1://"Eulerian time-domain ...",Acta Acustica united with Acustica,Vol.88(2002) 483-492
		{
			//////only for H=4////
			if(z>=5.0){
				return b*5.0;
			}else{
				return b*z;
			}
		}
		break;
	case 2://"study of the sound-vortex interaction ...", Eur.Phys.J.B 32,237-242(2003)
		{
			/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if((r/L)<=0.000001) r=1e-20;
			Vr=circulation/(2*PI*r)*(1-exp(-alpha*r*r/L/L));
			return (z-z0)/r*Vr;
			*/
			double Rr,y_center,z_center;
			y_center=20.0;
			z_center=80.0;
			Rr=sqrt(pow(y-y_center,2)+pow(z-z_center,2));
			if (abs(y-y_center)<=1.0e-6 && abs(z-z_center)<=1.0e-6) return 0.0;
			else if (y-y_center>=0.0) alpha=asin((z-z_center)/Rr);
			else if (y-y_center<=0.0) alpha=-asin((z-z_center)/Rr)+3.1415926;
			return -238.096/Rr*(1-exp(-1.26*pow(Rr,2)))*sin(alpha);
		  }
		break;
	case 3://"transmission of sound through a single vortex", Eur.Phys.J.B 37,229-239(2004)
		{
			/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if(r<=L){//Vr=circulation/(2*PI*L)*r/L;Vav=(z-z0)/r*Vr
				return circulation/(2*PI*L)/L*(z-z0);
			}
			Vr=circulation/(2*PI)/r;
			return (z-z0)/r*Vr;
			*/

			//anticlockwise dirction 
			double Rr,y_center,z_center,vv1,vv2;
			y_center=20.0;
			z_center=81.0;
			Rr=sqrt(pow(y-y_center,2)+pow(z-z_center,2));
			if (abs(y-y_center)<=1.0e-6 && abs(z-z_center)<=1.0e-6) vv1=0.0;
			else{
				if (y-y_center>=0.0) alpha=asin((z-z_center)/Rr);
				else if (y-y_center<=0.0) alpha=-asin((z-z_center)/Rr)+3.1415926;
				vv1=-238.0958/Rr*(1-exp(-1.26*pow(Rr,2)))*sin(alpha);
			}
			
			//clockwise dirction 
			y_center=20.0;
			z_center=79.0;
			Rr=sqrt(pow(y-y_center,2)+pow(z-z_center,2));
			if (abs(y-y_center)<=1.0e-6 && abs(z-z_center)<=1.0e-6) vv2=0.0;
			else{ 
				if (y-y_center>=0.0) alpha=asin((z-z_center)/Rr);
				else if (y-y_center<=0.0) alpha=-asin((z-z_center)/Rr)+3.1415926;
				vv2=238.0958/Rr*(1-exp(-1.26*pow(Rr,2)))*sin(alpha);
			}
			return vv1+vv2;
		 }
		break;
	case 4://reverse coordinate from case 1
		{
			return 0;
		}
		break;
	case 5://book of Solomans: terrain
		{
			return b*log(fabs(z)/0.1+1);
		}
		break;
	case 6:
		{
		/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if((r/L)<=0.00001) return 0;
			Vr=circulation/(2*PI*r)*(1-exp(-alpha*r*r/L/L));
			return (z-z0)/r*Vr;
		*/
		}
		break;
	default:
		return 0;
	}
}

double SpeedProfile::Wav(double x,double y,double z)
{
	switch(velocity_method){
	case 1:
		{
			return 0;
		}
		break;
	case 2:
		{
		/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if((r/L)<=0.000001) r=1e-10;
			Vr=circulation/(2*PI*r)*(1-exp(-alpha*r*r/L/L));
			return -(y-y0)/r*Vr;
		*/
			double Rr,y_center,z_center;
			y_center=20.0;
			z_center=80.0;
			Rr=sqrt(pow(y-y_center,2)+pow(z-z_center,2));
			if (abs(y-y_center)<=1.0e-6 && abs(z-z_center)<=1.0e-6) return 0.0;
			else if (y-y_center>=0.0) alpha=asin((z-z_center)/Rr);
			else if (y-y_center<=0.0) alpha=-asin((z-z_center)/Rr)+3.1415926;
			return 238.096/Rr*(1-exp(-1.26*pow(Rr,2)))*cos(alpha);
		 }
		break;
	case 3:
		{
		/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if(r<=L){//Vr=circulation/(2*PI*L)*r/L;Wav=-(y-y0)/r*Vr;
				return circulation/(2*PI*L)/L*(-(y-y0));
			}
			Vr=circulation/(2*PI)/r;
			return -(y-y0)/r*Vr;
		*/
			// Anticlockwise dirction 
			double Rr,y_center,z_center,ww1,ww2;
			y_center=20.0;
			z_center=81.0;
			Rr=sqrt(pow(y-y_center,2)+pow(z-z_center,2));
			if (abs(y-y_center)<=1.0e-6 && abs(z-z_center)<=1.0e-6) ww1=0.0;
			else{
				if (y-y_center>=0.0) alpha=asin((z-z_center)/Rr);
				else if (y-y_center<=0.0) alpha=-asin((z-z_center)/Rr)+3.1415926;
				ww1=238.0958/Rr*(1-exp(-1.26*pow(Rr,2)))*cos(alpha);
			}

			// Clockwise dirction 
			y_center=20.0;
			z_center=79.0;
			Rr=sqrt(pow(y-y_center,2)+pow(z-z_center,2));
			if (abs(y-y_center)<=1.0e-6 && abs(z-z_center)<=1.0e-6) ww2=0.0;
			else{
				if (y-y_center>=0.0) alpha=asin((z-z_center)/Rr);
				else if (y-y_center<=0.0) alpha=-asin((z-z_center)/Rr)+3.1415926;
				ww2=-238.0958/Rr*(1-exp(-1.26*pow(Rr,2)))*cos(alpha);
			}
			return ww1+ww2;
		}
		break;
	case 4:
		{
			return 0;
		}
		break;
	case 5:
		{
			return b*y;
		}
		break;
	case 6:
		{
		/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if((r/L)<=0.00001) return 0;
			Vr=circulation/(2*PI*r)*(1-exp(-alpha*r*r/L/L));
			return -(y-y0)/r*Vr;
		*/
		}
		break;
	default:
		return 0;//b*log(fabs(27-y)/0.1+1)
	}
}
SpeedProfile::~SpeedProfile(){}
//---------------------------------air member function---------------------------------------//

air::air(scheme_list scheme1,DifferenceStep Step,char *coordi,MovingFrame MF,const double gauss_width,
		 double velocity_coef,int velocity_method1,AirStruct AirPara,double AbsorbZmax,PoreStruct PorePara,PoreStruct PoreParb,
		 hillpore hill1,vortex vortex1,int mpi_rank1,int mpi_size1,int mpi_yarea1,int mpi_zarea1,int mpi_porous1)
   :calculation_2D(scheme1,Step,coordi,MF,gauss_width,AirPara,vortex1,mpi_rank1,mpi_size1,mpi_yarea1,mpi_zarea1,mpi_porous1)
{
	int i,j;
	adiabatic_coef=AirPara.adiabatic_coef;
	Pav=AirPara.Pav; //Aav=AirPara.Aav,sound_speed=AirPara.sound_speed;
	velocity_method=velocity_method1;
	SP=new SpeedProfile(velocity_coef,velocity_method); //velocity_coef,AirPara.sound_speed,velocity_method
	//calculate partial differential of wind speed
	Vav=new double* [mpi_IMAX[mpi_iindex]];Wav=new double* [mpi_IMAX[mpi_iindex]];
	Vav_over_y=new double* [mpi_IMAX[mpi_iindex]];Vav_over_z=new double* [mpi_IMAX[mpi_iindex]];
	Wav_over_y=new double* [mpi_IMAX[mpi_iindex]];Wav_over_z=new double* [mpi_IMAX[mpi_iindex]];
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		Vav[i]=new double [mpi_JMAX[mpi_jindex]];Wav[i]=new double [mpi_JMAX[mpi_jindex]];
		Vav_over_y[i]=new double [mpi_JMAX[mpi_jindex]];Vav_over_z[i]=new double [mpi_JMAX[mpi_jindex]];
		Wav_over_y[i]=new double [mpi_JMAX[mpi_jindex]];Wav_over_z[i]=new double [mpi_JMAX[mpi_jindex]];
	}

	Qv=new double*[mpi_IMAX[mpi_iindex]];Qw=new double *[mpi_IMAX[mpi_iindex]];Qp=new double *[mpi_IMAX[mpi_iindex]];
	for (i=0;i<mpi_IMAX[mpi_iindex];i++){
		Qv[i]=new double[mpi_JMAX[mpi_jindex]];
		Qw[i]=new double[mpi_JMAX[mpi_jindex]];
		Qp[i]=new double[mpi_JMAX[mpi_jindex]];
	}
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			Qv[i][j]=0;Qw[i][j]=0;Qp[i][j]=0;
		}
	}
	whole_Qv=new double*[IMAX];whole_Qw=new double*[IMAX];whole_Qp=new double*[IMAX];
	for (i=0;i<IMAX;i++){
		whole_Qv[i]=new double[mpi_JMAX[mpi_jindex]];
		whole_Qw[i]=new double[mpi_JMAX[mpi_jindex]];
		whole_Qp[i]=new double[mpi_JMAX[mpi_jindex]];
	}
	for(i=0;i<IMAX;i++){
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			whole_Qv[i][j]=0;whole_Qw[i][j]=0;whole_Qp[i][j]=0;
		}
	}

	PML_AbsorbZmax=AbsorbZmax;PML_alpha=2.0;
	//if (mpi_rank == 0) cout << "PML_Height1 = " << PML_Height1 << ", PML_z1 = " << PML_z1 << ", PML_Width1 = " << PML_Width1 << ", PML_y1 = " << PML_y1 << endl;
	if(JMAX1!=0){
		PML_Height1=fabs(leftbound_Z[JMAX1-1]-leftbound_Z[0]);
		PML_z1=leftbound_Z[JMAX1-1];
	}
	if(JMAX2!=0){
		PML_Height2=fabs(leftbound_Z[JMAX-1]-leftbound_Z[JMAX-JMAX2]);
		PML_z2=leftbound_Z[JMAX-JMAX2];
	}
	// PML for y direction
	PML_AbsorbYmax = AbsorbZmax; PML_alphay = PML_alpha;
	/*if (IMAX2 != 0){ ////////////////// PML RH
		PML_y1 = bottombound_Y[IMAX2 - 1];//bottombound_Y[0];
		PML_y2 = bottombound_Y[IMAX - IMAX2];
		PML_Width1 = fabs(bottombound_Y[IMAX2 - 1] - bottombound_Y[0]);// fabs(bottombound_Y[0] - bottombound_Y[0]);
		PML_Width2 = fabs(bottombound_Y[IMAX - 1] - bottombound_Y[IMAX - IMAX2]);
		//PML_AbsorbYmax = 1e3;
		//PML_alphay = PML_alpha;
		if (mpi_rank == 0) cout << "PML_y1 = " << PML_y1 << ", PML_y2 = " << PML_y2 << ", PML_Width1 = " << PML_Width1 << ", PML_Width2 = " << PML_Width2 << endl;
	}*/
	PML_y1 = 0;
	PML_y2 = 0;
	PML_Width1 = 0;
	PML_Width2 = 0;
	if (IMAX01 != 0){ // PML LH
		if (mpi_rank == 0) cout << "IMAX01 = " << IMAX01 << ", bottombound_Y[IMAX01 - 1] = " << bottombound_Y[IMAX01 - 1] << endl;		
		PML_Width1 = fabs(bottombound_Y[IMAX01 - 1] - bottombound_Y[0]);
		PML_y1 = bottombound_Y[IMAX01 - 1];
	}
	if (IMAX02 != 0){ // PML RH
		PML_Width2 = fabs(bottombound_Y[IMAX - 1] - bottombound_Y[IMAX - IMAX02]);
		PML_y2 = bottombound_Y[IMAX - IMAX02];
		if (mpi_rank == 0) cout << "PML_y1 = " << PML_y1 << ", PML_y2 = " << PML_y2 << ", PML_Width1 = " << PML_Width1 << ", PML_Width2 = " << PML_Width2 << endl;
	}
	// if cr_judge is 1, porous media hill exists there
	cr_judge=hill1.curve_judge;
	ystar=hill1.curve_ystar;
	ystop=hill1.curve_ystop;
	zstar=hill1.curve_zstar;
	zstop=hill1.curve_zstop;
	cr_a=hill1.curve_coefa;
	cr_b=hill1.curve_coefb;
	cr_c=hill1.curve_coefc;
	cr_d=hill1.curve_coefd;
	cr_e=hill1.curve_coefe;
	cr_f=hill1.curve_coeff;
	cr_g=hill1.curve_coefg;
	cr_h=hill1.curve_coefh;
	cr_i=hill1.curve_coefi;
	cr_j=hill1.curve_coefj;
	cr_k=hill1.curve_coefk;
	cr_l=hill1.curve_coefl;
	cr_points=hill1.curve_points;
/*
	crpo_eff_density=PorePara.eff_density;
	crpo_porosity=PorePara.porosity;
	crpo_resistivity=PorePara.resistivity;
	crpo_Kp=PorePara.Kp;
// crpo_reistivity should be a variable to be chose
	crpo_Beta = crpo_resistivity*diff_t / crpo_eff_density; // 0 -> ground 3e5 1.5e5
	crpo_Beta1 = 4e8*diff_t / crpo_eff_density; // 1 -> wall 4e8
	crpo_Beta2 = 5e3*diff_t / crpo_eff_density; // 2 -> ground 1.5e5
	crpo_Gama=diff_t/crpo_Kp/crpo_porosity;*/
	
	crpo_eff_density_a = PorePara.eff_density;
	crpo_porosity_a = PorePara.porosity;
	crpo_resistivity_a = PorePara.resistivity;
	crpo_Kp_a = PorePara.Kp;
	crpo_Beta_a = crpo_resistivity_a*diff_t / crpo_eff_density_a;
	crpo_Gama_a = diff_t/crpo_Kp_a/crpo_porosity_a;
	
	crpo_eff_density_b = PoreParb.eff_density;
	crpo_porosity_b = PoreParb.porosity;
	crpo_resistivity_b = PoreParb.resistivity;
	crpo_Kp_b = PoreParb.Kp;
	crpo_Beta_b = crpo_resistivity_b*diff_t / crpo_eff_density_b;
	crpo_Gama_b = diff_t/crpo_Kp_b/crpo_porosity_b;

	cr_y=new double [cr_points]; cr_z=new double [cr_points]; cr_uc=new int *[mpi_IMAX[mpi_iindex]];
	for (i=0;i<mpi_IMAX[mpi_iindex];i++){
		cr_uc[i]=new int [mpi_JMAX[mpi_jindex]];
	}
	// for vortex
	ystar_vor=vortex1.ystar;
	ystop_vor=vortex1.ystop;
	zstar_vor=vortex1.zstar;
	zstop_vor=vortex1.zstop;

	/*double temp_vav, temp_wav;
	ifstream Vav_file("vy_vz.dat",ios::in|ios::binary);
	if(!Vav_file){
		cout<<"cannot open file:" << "vy_vz.dat" <<" for read!!!"<<endl;
		return;
	}
	int in,jn;
	for(i=0;i<whole_IMAX;i++)
	{
		for(j=0;j<JMAX;j++)
		{
			if (i>=500 && i<=1000 &&j>=100 && j<=1100) //adding area
			{
				Vav_file>>temp_vav>>temp_wav;
				if (j>=mpi_j1[mpi_jindex] && j<=mpi_j2[mpi_jindex]&& i>=mpi_i1[mpi_iindex] && i<=mpi_i2[mpi_iindex])
				{
					in=i-mpi_i1[mpi_iindex];
					jn=j-mpi_j1[mpi_jindex];
					Vav[in][jn]=temp_vav;
					Wav[in][jn]=temp_wav;
					//if(mpi_rank==10)cout << "temp_vav = " << temp_vav << endl;
				}
			}
			else
			{
				if (j>=mpi_j1[mpi_jindex] && j<=mpi_j2[mpi_jindex]&& i>=mpi_i1[mpi_iindex] && i<=mpi_i2[mpi_iindex])
				{
					in=i-mpi_i1[mpi_iindex];
					jn=j-mpi_j1[mpi_jindex];
					Vav[in][jn]=0.0;
					Wav[in][jn]=0.0;
				}
			}
		}
	}
	Vav_file.close();*/
}

air::~air()
{
	for(int i=0;i<mpi_IMAX[mpi_iindex];i++){
		delete[] Vav[i];delete[] Wav[i];
		delete[] Vav_over_y[i];delete[] Vav_over_z[i];
		delete[] Wav_over_y[i];delete[] Wav_over_z[i];
		delete[] Qv[i];delete[] Qw[i];delete[] Qp[i];
		delete[] cr_uc[i];
	}
	for (int i=0;i<IMAX;i++){
		delete[] whole_Qv[i];
		delete[] whole_Qw[i];
		delete[] whole_Qp[i];
	}
	delete[] whole_Qv;delete whole_Qw; delete whole_Qp;
	delete[] Vav;delete[] Wav;
	delete[] Vav_over_y;delete[] Vav_over_z;
	delete[] Wav_over_y;delete[] Wav_over_z;
	delete[] Qv; delete[] Qw;delete[] Qp;delete[] cr_uc;
	delete[] cr_y;delete[] cr_z;
	delete SP;
}

void air::SetWindProfile(int N)
{
	int i,j;
	//calculate Z0 on local ground level
	/*
	for(j=0;j<mpi_JMAX[mpi_jindex];j++){
		for(i=0;i<mpi_IMAX[mpi_iindex];i++){
			double Z0;
			Z0=0;
			Vav[i][j]=SP->Vav(0,Y[i][j],Z[i][j]-Z0);
			Wav[i][j]=SP->Wav(0,Y[i][j],Z[i][j]-Z0);
		}
	}
	*/
	for(j=jj1;j<mpi_JMAX[mpi_jindex]-jj2;j++){
		for(i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){

			Vav_over_y[i][j]=(Z_over_z[i][j]*(Vav[i+1][j]-Vav[i-1][j])/diff_y/2-
					Z_over_y[i][j]*(Vav[i][j+1]-Vav[i][j-1])/diff_z/2)/Jacobian[i][j];
			Wav_over_y[i][j]=(Z_over_z[i][j]*(Wav[i+1][j]-Wav[i-1][j])/diff_y/2-
					Z_over_y[i][j]*(Wav[i][j+1]-Wav[i][j-1])/diff_z/2)/Jacobian[i][j];
			Vav_over_z[i][j]=(Y_over_y[i][j]*(Vav[i][j+1]-Vav[i][j-1])/diff_z/2-
					Y_over_z[i][j]*(Vav[i+1][j]-Vav[i-1][j])/diff_y/2)/Jacobian[i][j];
			Wav_over_z[i][j]=(Y_over_y[i][j]*(Wav[i][j+1]-Wav[i][j-1])/diff_z/2-
					Y_over_z[i][j]*(Wav[i+1][j]-Wav[i-1][j])/diff_y/2)/Jacobian[i][j];

		}
	}
}
void air::save_restart_air(char *restartfile)
{
	int i,j;
	ofstream outfile11(restartfile,ios::app|ios::binary);
	outfile11.setf(ios::scientific,ios::floatfield);
	outfile11.setf(ios::showpos);
	outfile11.precision(6);
	outfile11.width(14);
	for (i=0;i<IMAX;i++){
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			outfile11<<whole_Qv[i][j]<<" ";
			outfile11<<whole_Qw[i][j]<<" ";
			outfile11<<whole_Qp[i][j]<<" ";
		}
	}
	outfile11<<endl;
	outfile11.close();
}

void air::input_restart_air(char *restart_infile,int *point_pois)
{
	int i,j,int_pois;
	int_pois=*point_pois;
	ifstream infile(restart_infile,ios::in|ios::binary);
	infile.seekg(int_pois);
	for (i=0;i<IMAX;i++){
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			infile>>whole_Qv[i][j];
			infile>>whole_Qw[i][j];
			infile>>whole_Qp[i][j];
		}
	}
	infile.ignore(100,'\n');
	int_pois=infile.tellg();
	*point_pois=int_pois;
	infile.close();
}
//--------------------------calculate velocity of V and W for new time (nn)----------------//
void air::cal_fvw(double **temp_v,double **temp_w,
			double **temp_p)
{
	int i,in,j,jn;
	double absorb_coefz1, absorb_coefz2, temp_z, absorb_coefz3, absorb_coefz4, absorb_coefy1;
	double vav_a,vav_b,wav_a,wav_b,vav_qva,vav_qvb,vav_qwa,vav_qwb;
	//  PML boundary for y direction
	double absorb_coefy,temp_y;
	for(i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++)
	{
		for(j=jj1;j<mpi_JMAX[mpi_jindex]-jj2;j++)
		{
			in=mpi_i1[mpi_iindex]+i;
			jn=mpi_j1[mpi_jindex]+j;
			if (Vav[i][j]<=0.0)
			{
				vav_a=temp_v[i+1][j]-temp_v[i][j];
				vav_b=temp_v[i][j+1]-temp_v[i][j];
				vav_qva=Qv[i+1][j]-Qv[i][j];
				vav_qvb=Qv[i][j+1]-Qv[i][j];
				vav_qwa=Qw[i+1][j]-Qw[i][j];
				vav_qwb=Qw[i][j+1]-Qw[i][j];
			}
			else
			{
				vav_a=temp_v[i][j]-temp_v[i-1][j];
				vav_b=temp_v[i][j]-temp_v[i][j-1];
				vav_qva=Qv[i][j]-Qv[i-1][j];
				vav_qvb=Qv[i][j]-Qv[i][j-1];
				vav_qwa=Qw[i][j]-Qw[i-1][j];
				vav_qwb=Qw[i][j]-Qw[i][j-1];
			}
			if (Wav[i][j]<=0.0)
			{
				wav_a=temp_v[i+1][j]-temp_v[i][j];
				wav_b=temp_v[i][j+1]-temp_v[i][j];
			}
			else
			{
				wav_a=temp_v[i][j]-temp_v[i-1][j];
				wav_b=temp_v[i][j]-temp_v[i][j-1];
			}

			if (n_order==2 && Y[i][j]>=ystar_vor-1.0e-6 && Y[i][j]<=ystop_vor+1.0e-6 
						&& Z[i][j]>=zstar_vor-1.0e-6 && Z[i][j]<=zstop_vor+1.0e-6)
			{
				if (Vav[i][j]<=0.0)
				{
					vav_a=(-1.0*temp_v[i+2][j]+4.0*temp_v[i+1][j]-3.0*temp_v[i][j])*0.5;
					vav_b=(-1.0*temp_v[i][j+2]+4.0*temp_v[i][j+1]-3.0*temp_v[i][j])*0.5;
					vav_qva=(-1.0*Qv[i+2][j]+4.0*Qv[i+1][j]-3.0*Qv[i][j])*0.5;
					vav_qvb=(-1.0*Qv[i][j+2]+4.0*Qv[i][j+1]-3.0*Qv[i][j])*0.5;
					vav_qwa=(-1.0*Qw[i+2][j]+4.0*Qw[i+1][j]-3.0*Qw[i][j])*0.5;
					vav_qwb=(-1.0*Qw[i][j+2]+4.0*Qw[i][j+1]-3.0*Qw[i][j])*0.5;
				}
				else
				{
					vav_a=(3.0*temp_v[i][j]-4.0*temp_v[i-1][j]+temp_v[i-2][j])*0.5;
					vav_b=(3.0*temp_v[i][j]-4.0*temp_v[i][j-1]+temp_v[i][j-2])*0.5;
					vav_qva=(3.0*Qv[i][j]-4.0*Qv[i-1][j]+Qv[i-2][j])*0.5;
					vav_qvb=(3.0*Qv[i][j]-4.0*Qv[i][j-1]+Qv[i][j-2])*0.5;
					vav_qwa=(3.0*Qw[i][j]-4.0*Qw[i-1][j]+Qw[i-2][j])*0.5;
					vav_qwb=(3.0*Qw[i][j]-4.0*Qw[i][j-1]+Qw[i][j-2])*0.5;
				}
				if (Wav[i][j]<=0.0)
				{
					wav_a=(-1.0*temp_v[i+2][j]+4.0*temp_v[i+1][j]-3.0*temp_v[i][j])*0.5;
					wav_b=(-1.0*temp_v[i][j+2]+4.0*temp_v[i][j+1]-3.0*temp_v[i][j])*0.5;
				}
				else
				{
					wav_a=(3.0*temp_v[i][j]-4.0*temp_v[i-1][j]+temp_v[i-2][j])*0.5;
					wav_b=(3.0*temp_v[i][j]-4.0*temp_v[i][j-1]+temp_v[i][j-2])*0.5;
				}
			}
			else if (n_order==3 && Y[i][j]>=ystar_vor-1.0e-6 && Y[i][j]<=ystop_vor+1.0e-6 
						&& Z[i][j]>=zstar_vor-1.0e-6 && Z[i][j]<=zstop_vor+1.0e-6)
			{
				if (Vav[i][j]<=0.0)
				{
					vav_a=(-1.0*temp_v[i+2][j]+6.0*temp_v[i+1][j]-3.0*temp_v[i][j]-2.0*temp_v[i-1][j])*0.166667;
					vav_b=(-1.0*temp_v[i][j+2]+6.0*temp_v[i][j+1]-3.0*temp_v[i][j]-2.0*temp_v[i][j-1])*0.166667;
					vav_qva=(-1.0*Qv[i+2][j]+6.0*Qv[i+1][j]-3.0*Qv[i][j]-2.0*Qv[i-1][j])*0.166667;
					vav_qvb=(-1.0*Qv[i][j+2]+6.0*Qv[i][j+1]-3.0*Qv[i][j]-2.0*Qv[i][j-1])*0.166667;
					vav_qwa=(-1.0*Qw[i+2][j]+6.0*Qw[i+1][j]-3.0*Qw[i][j]-2.0*Qw[i-1][j])*0.166667;
					vav_qwb=(-1.0*Qw[i][j+2]+6.0*Qw[i][j+1]-3.0*Qw[i][j]-2.0*Qw[i][j-1])*0.166667;
				}
				else
				{
					vav_a=(2.0*temp_v[i+1][j]+3.0*temp_v[i][j]-6.0*temp_v[i-1][j]+temp_v[i-2][j])*0.166667;
					vav_b=(2.0*temp_v[i][j+1]+3.0*temp_v[i][j]-6.0*temp_v[i][j-1]+temp_v[i][j-2])*0.166667;
					vav_qva=(2.0*Qv[i+1][j]+3.0*Qv[i][j]-6.0*Qv[i-1][j]+Qv[i-2][j])*0.166667;
					vav_qvb=(2.0*Qv[i][j+1]+3.0*Qv[i][j]-6.0*Qv[i][j-1]+Qv[i][j-2])*0.166667;
					vav_qwa=(2.0*Qw[i+1][j]+3.0*Qw[i][j]-6.0*Qw[i-1][j]+Qw[i-2][j])*0.166667;
					vav_qwb=(2.0*Qw[i][j+1]+3.0*Qw[i][j]-6.0*Qw[i][j-1]+Qw[i][j-2])*0.166667;
				}
				if (Wav[i][j]<=0.0)
				{
					wav_a=(-1.0*temp_v[i+2][j]+6.0*temp_v[i+1][j]-3.0*temp_v[i][j]-2.0*temp_v[i-1][j])*0.166667;
					wav_b=(-1.0*temp_v[i][j+2]+6.0*temp_v[i][j+1]-3.0*temp_v[i][j]-2.0*temp_v[i][j-1])*0.166667;
				}
				else
				{
					wav_a=(2.0*temp_v[i+1][j]+3.0*temp_v[i][j]-6.0*temp_v[i-1][j]+temp_v[i-2][j])*0.166667;
					wav_b=(2.0*temp_v[i][j+1]+3.0*temp_v[i][j]-6.0*temp_v[i][j-1]+temp_v[i][j-2])*0.166667;
				}
			}
			else if (n_order==5 && Y[i][j]>=ystar_vor-1.0e-6 && Y[i][j]<=ystop_vor+1.0e-6 
						&& Z[i][j]>=zstar_vor-1.0e-6 && Z[i][j]<=zstop_vor+1.0e-6)
			{
				if (Vav[i][j]<=0.0)
				{
					vav_a=-1.0*weno_scheme(temp_v[i+3][j],temp_v[i+2][j],temp_v[i+1][j],
								temp_v[i][j],temp_v[i-1][j],temp_v[i-2][j],diff_y);
					vav_b=-1.0*weno_scheme(temp_v[i][j+3],temp_v[i][j+2],temp_v[i][j+1],
								temp_v[i][j],temp_v[i][j-1],temp_v[i][j-2],diff_z);
					vav_qva=-1.0*weno_scheme(Qv[i+3][j],Qv[i+2][j],Qv[i+1][j],
								Qv[i][j],Qv[i-1][j],Qv[i-2][j],diff_y);
					vav_qvb=-1.0*weno_scheme(Qv[i][j+3],Qv[i][j+2],Qv[i][j+1],
								Qv[i][j],Qv[i][j-1],Qv[i][j-2],diff_z);
					vav_qwa=-1.0*weno_scheme(Qw[i+3][j],Qw[i+2][j],Qw[i+1][j],
								Qw[i][j],Qw[i-1][j],Qw[i-2][j],diff_y);
					vav_qwb=-1.0*weno_scheme(Qw[i][j+3],Qw[i][j+2],Qw[i][j+1],
								Qw[i][j],Qw[i][j-1],Qw[i][j-2],diff_z);
				}
				else
				{

					vav_a=weno_scheme(temp_v[i-3][j],temp_v[i-2][j],temp_v[i-1][j],
								temp_v[i][j],temp_v[i+1][j],temp_v[i+2][j],diff_y);
					vav_b=weno_scheme(temp_v[i][j-3],temp_v[i][j-2],temp_v[i][j-1],
								temp_v[i][j],temp_v[i][j+1],temp_v[i][j+2],diff_z);
					vav_qva=weno_scheme(Qv[i-3][j],Qv[i-2][j],Qv[i-1][j],
								Qv[i][j],Qv[i+1][j],Qv[i+2][j],diff_y);
					vav_qvb=weno_scheme(Qv[i][j-3],Qv[i][j-2],Qv[i][j-1],
								Qv[i][j],Qv[i][j+1],Qv[i][j+2],diff_z);
					vav_qwa=weno_scheme(Qw[i-3][j],Qw[i-2][j],Qw[i-1][j],
								Qw[i][j],Qw[i+1][j],Qw[i+2][j],diff_y);
					vav_qwb=weno_scheme(Qw[i][j-3],Qw[i][j-2],Qw[i][j-1],
								Qw[i][j],Qw[i][j+1],Qw[i][j+2],diff_z);

				}
				if (Wav[i][j]<=0.0)
				{
					wav_a=-1.0*weno_scheme(temp_v[i+3][j],temp_v[i+2][j],temp_v[i+1][j],
								temp_v[i][j],temp_v[i-1][j],temp_v[i-2][j],diff_y);
					wav_b=-1.0*weno_scheme(temp_v[i][j+3],temp_v[i][j+2],temp_v[i][j+1],
								temp_v[i][j],temp_v[i][j-1],temp_v[i][j-2],diff_z);
				}
				else
				{
					wav_a=weno_scheme(temp_v[i-3][j],temp_v[i-2][j],temp_v[i-1][j],
								temp_v[i][j],temp_v[i+1][j],temp_v[i+2][j],diff_y);
					wav_b=weno_scheme(temp_v[i][j-3],temp_v[i][j-2],temp_v[i][j-1],
								temp_v[i][j],temp_v[i][j+1],temp_v[i][j+2],diff_z);
				}
			}			
			v_nn[i][j]=-Vav[i][j]*(Z_over_z[i][j]*diff_t/diff_y*vav_a-Z_over_y[i][j]*diff_t/diff_z*vav_b)/Jacobian[i][j]
						-Wav[i][j]*(Y_over_y[i][j]*diff_t/diff_z*wav_b-Y_over_z[i][j]*diff_t/diff_y*wav_a)/Jacobian[i][j]
				-Aav[j]/Jacobian[i][j]*(Z_over_z[i][j]*diff_t/diff_y*(temp_p[i][j]-temp_p[i-1][j])-
				Z_over_y[i][j]*diff_t/diff_z*(temp_p[i][j]-temp_p[i][j-1]))-
				temp_v[i][j]*diff_t*Vav_over_y[i][j]-temp_w[i][j]*diff_t*Vav_over_z[i][j];

			if (Vav[i][j]<=0.0)
			{
				vav_a=temp_w[i+1][j]-temp_w[i][j];
				vav_b=temp_w[i][j+1]-temp_w[i][j];
			}
			else
			{
				vav_a=temp_w[i][j]-temp_w[i-1][j];
				vav_b=temp_w[i][j]-temp_w[i][j-1];
			}
			if (Wav[i][j]<=0.0)
			{
				wav_a=temp_w[i+1][j]-temp_w[i][j];
				wav_b=temp_w[i][j+1]-temp_w[i][j];
			}
			else
			{
				wav_a=temp_w[i][j]-temp_w[i-1][j];
				wav_b=temp_w[i][j]-temp_w[i][j-1];
			}
			if(n_order==2 && Y[i][j]>=ystar_vor-1.0e-6 && Y[i][j]<=ystop_vor+1.0e-6 
						&& Z[i][j]>=zstar_vor-1.0e-6&& Z[i][j]<=zstop_vor+1.0e-6)
			{
				if (Vav[i][j]<=0.0)
				{
					vav_a=(-1.0*temp_w[i+2][j]+4.0*temp_w[i+1][j]-3.0*temp_w[i][j])*0.5;
					vav_b=(-1.0*temp_w[i][j+2]+4.0*temp_w[i][j+1]-3.0*temp_w[i][j])*0.5;
				}
				else
				{
					vav_a=(3.0*temp_w[i][j]-4.0*temp_w[i-1][j]+temp_w[i-2][j])*0.5;
					vav_b=(3.0*temp_w[i][j]-4.0*temp_w[i][j-1]+temp_w[i][j-2])*0.5;
				}
				if (Wav[i][j]<=0.0)
				{
					wav_a=(-1.0*temp_w[i+2][j]+4.0*temp_w[i+1][j]-3.0*temp_w[i][j])*0.5;
					wav_b=(-1.0*temp_w[i][j+2]+4.0*temp_w[i][j+1]-3.0*temp_w[i][j])*0.5;
				}
				else
				{
					wav_a=(3.0*temp_w[i][j]-4.0*temp_w[i-1][j]+temp_w[i-2][j])*0.5;
					wav_b=(3.0*temp_w[i][j]-4.0*temp_w[i][j-1]+temp_w[i][j-2])*0.5;
				}
			}
			else if(n_order==3 && Y[i][j]>=ystar_vor-1.0e-6 && Y[i][j]<=ystop_vor+1.0e-6 
						&& Z[i][j]>=zstar_vor-1.0e-6&& Z[i][j]<=zstop_vor+1.0e-6)
			{
				if (Vav[i][j]<=0.0)
				{
					vav_a=(-1.0*temp_w[i+2][j]+6.0*temp_w[i+1][j]-3.0*temp_w[i][j]-2.0*temp_w[i-1][j])*0.166667;
					vav_b=(-1.0*temp_w[i][j+2]+6.0*temp_w[i][j+1]-3.0*temp_w[i][j]-2.0*temp_w[i][j-1])*0.166667;
				}
				else
				{
					vav_a=(2.0*temp_w[i+1][j]+3.0*temp_w[i][j]-6.0*temp_w[i-1][j]+temp_w[i-2][j])*0.166667;
					vav_b=(2.0*temp_w[i][j+1]+3.0*temp_w[i][j]-6.0*temp_w[i][j-1]+temp_w[i][j-2])*0.166667;
				}
				if (Wav[i][j]<=0.0)
				{
					wav_a=(-1.0*temp_w[i+2][j]+6.0*temp_w[i+1][j]-3.0*temp_w[i][j]-2.0*temp_w[i-1][j])*0.166667;
					wav_b=(-1.0*temp_w[i][j+2]+6.0*temp_w[i][j+1]-3.0*temp_w[i][j]-2.0*temp_w[i][j-1])*0.166667;
				}
				else
				{
					wav_a=(2.0*temp_w[i+1][j]+3.0*temp_w[i][j]-6.0*temp_w[i-1][j]+temp_w[i-2][j])*0.166667;
					wav_b=(2.0*temp_w[i][j+1]+3.0*temp_w[i][j]-6.0*temp_w[i][j-1]+temp_w[i][j-2])*0.166667;
				}
			}
			else if(n_order==5 && Y[i][j]>=ystar_vor-1.0e-6 && Y[i][j]<=ystop_vor+1.0e-6 
						&& Z[i][j]>=zstar_vor-1.0e-6&& Z[i][j]<=zstop_vor+1.0e-6)
			{
				if (Vav[i][j]<=0.0)
				{
					vav_a=-1.0*weno_scheme(temp_w[i+3][j],temp_w[i+2][j],temp_w[i+1][j],
								temp_w[i][j],temp_w[i-1][j],temp_w[i-2][j],diff_y);
					vav_b=-1.0*weno_scheme(temp_w[i][j+3],temp_w[i][j+2],temp_w[i][j+1],
								temp_w[i][j],temp_w[i][j-1],temp_w[i][j-2],diff_z);
				}
				else
				{
					vav_a=weno_scheme(temp_w[i-3][j],temp_w[i-2][j],temp_w[i-1][j],
								temp_w[i][j],temp_w[i+1][j],temp_w[i+2][j],diff_y);
					vav_b=weno_scheme(temp_w[i][j-3],temp_w[i][j-2],temp_w[i][j-1],
								temp_w[i][j],temp_w[i][j+1],temp_w[i][j+2],diff_z);
				}
				if (Wav[i][j]<=0.0)
				{
					wav_a=-1.0*weno_scheme(temp_w[i+3][j],temp_w[i+2][j],temp_w[i+1][j],
								temp_w[i][j],temp_w[i-1][j],temp_w[i-2][j],diff_y);
					wav_b=-1.0*weno_scheme(temp_w[i][j+3],temp_w[i][j+2],temp_w[i][j+1],
								temp_w[i][j],temp_w[i][j-1],temp_w[i][j-2],diff_z);
				}
				else
				{
					wav_a=weno_scheme(temp_w[i-3][j],temp_w[i-2][j],temp_w[i-1][j],
								temp_w[i][j],temp_w[i+1][j],temp_w[i+2][j],diff_y);
					wav_b=weno_scheme(temp_w[i][j-3],temp_w[i][j-2],temp_w[i][j-1],
								temp_w[i][j],temp_w[i][j+1],temp_w[i][j+2],diff_z);
				}
			}			

			w_nn[i][j]=-Vav[i][j]*(Z_over_z[i][j]*diff_t/diff_y*vav_a-
					Z_over_y[i][j]*diff_t/diff_z*vav_b)/Jacobian[i][j]
					-Wav[i][j]*(Y_over_y[i][j]*diff_t/diff_z*wav_b-
					Y_over_z[i][j]*diff_t/diff_y*wav_a)/Jacobian[i][j]
				-Aav[j]/Jacobian[i][j]*(Y_over_y[i][j]*diff_t/diff_z*(temp_p[i][j]-temp_p[i][j-1])-
				Y_over_z[i][j]*diff_t/diff_y*(temp_p[i][j]-temp_p[i-1][j]))-
				temp_v[i][j]*diff_t*Wav_over_y[i][j]-temp_w[i][j]*diff_t*Wav_over_z[i][j];

			if (jn<JMAX1)
			{
				absorb_coefz1=PML_AbsorbZmax*pow((PML_z1-Z[i][j])/PML_Height1,PML_alpha);
				absorb_coefz3=1.0+25.0*pow((PML_z1-Z[i][j])/PML_Height1,PML_alpha);
				temp_z=Vav[i][j]*(vav_qva/diff_y*Z_over_z[i][j]-
									vav_qvb/diff_z*Z_over_y[i][j])/Jacobian[i][j];
				temp_z += Aav[j]*((Qp[i][j]-Qp[i-1][j])/diff_y*Z_over_z[i][j]-
					(Qp[i][j]-Qp[i][j-1])/diff_z*Z_over_y[i][j])/Jacobian[i][j];
				temp_z +=Qw[i][j]*Vav_over_z[i][j];
				v_nn[i][j] -= absorb_coefz1*(temp_z+temp_v[i][j])*diff_t;

				temp_z=Vav[i][j]*(vav_qwa/diff_y*Z_over_z[i][j]-
									vav_qwb/diff_z*Z_over_y[i][j])/Jacobian[i][j];
				w_nn[i][j] -= absorb_coefz1*(temp_z+temp_w[i][j])*diff_t+(1.0/absorb_coefz3-1.0)*
						(Y_over_y[i][j]*diff_t/diff_z*(temp_p[i][j]-temp_p[i][j-1])-
						Y_over_z[i][j]*diff_t/diff_y*(temp_p[i][j]-temp_p[i-1][j]))*Aav[j]/Jacobian[i][j];
				/*v_nn[i][j] = -1.0 / (1.0 + crpo_Beta*Aav[j])*(diff_t / crpo_eff_density / Jacobian[i][j] * Aav[j]
					* (Z_over_z[i][j] / diff_y*(temp_p[i][j] - temp_p[i - 1][j]) -
					Z_over_y[i][j] / diff_z*(temp_p[i][j] - temp_p[i][j - 1]))) - temp_v[i][j]
					+ temp_v[i][j] / (1.0 + crpo_Beta*Aav[j]);
				w_nn[i][j] = -1.0 / (1.0 + crpo_Beta*Aav[j])*(diff_t / crpo_eff_density / Jacobian[i][j] * Aav[j]
					* (Y_over_y[i][j] / diff_z*(temp_p[i][j] - temp_p[i][j - 1]) -
					Y_over_z[i][j] / diff_y*(temp_p[i][j] - temp_p[i - 1][j]))) - temp_w[i][j]
					+ temp_w[i][j] / (1.0 + crpo_Beta*Aav[j]);*/
			}

			if(jn>=JMAX-JMAX2)
			{
				absorb_coefz2=PML_AbsorbZmax*pow((Z[i][j]-PML_z2)/PML_Height2,PML_alpha);
				absorb_coefz4=1.0+25.0*pow((Z[i][j]-PML_z2)/PML_Height2,PML_alpha);	
				temp_z=Vav[i][j]*(vav_qva/diff_y*Z_over_z[i][j]-
								vav_qvb/diff_z*Z_over_y[i][j])/Jacobian[i][j];
				temp_z += Aav[j]*((Qp[i][j]-Qp[i-1][j])/diff_y*Z_over_z[i][j]-
						(Qp[i][j]-Qp[i][j-1])/diff_z*Z_over_y[i][j])/Jacobian[i][j];
				temp_z +=Qw[i][j]*Vav_over_z[i][j];
				v_nn[i][j] -= absorb_coefz2*(temp_z+temp_v[i][j])*diff_t;

				temp_z=Vav[i][j]*(vav_qwa/diff_y*Z_over_z[i][j]-
									vav_qwb/diff_z*Z_over_y[i][j])/Jacobian[i][j];
				w_nn[i][j] -= absorb_coefz2*(temp_z+temp_w[i][j])*diff_t+(1.0/absorb_coefz4-1.0)*
						(Y_over_y[i][j]*diff_t/diff_z*(temp_p[i][j]-temp_p[i][j-1])-
						Y_over_z[i][j]*diff_t/diff_y*(temp_p[i][j]-temp_p[i-1][j]))*Aav[j]/Jacobian[i][j];
			}
			if (in<IMAX01 && IMAX1 >= 9){

				absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
				absorb_coefy1 = 1.0 + 25.0*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
				absorb_coefz1 = PML_AbsorbZmax*pow((PML_z1 - Z[i][j]) / PML_Height1, PML_alpha);
				absorb_coefz2 = PML_AbsorbZmax*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);

				temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
				v_nn[i][j] -= absorb_coefy * temp_y * diff_t + (1.0 / absorb_coefy1 - 1.0)*	Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] / diff_y * (temp_p[i][j] - temp_p[i - 1][j]) -
					Z_over_y[i][j] / diff_z * (temp_p[i][j] - temp_p[i][j - 1])) * diff_t ;

				temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
					Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) + temp_w[i][j];

				if (jn >= JMAX - JMAX2 && jn<JMAX){
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
					temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
					v_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * absorb_coefz2 * Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] / diff_y * (Qp[i][j] - Qp[i - 1][j]) -
						Z_over_y[i][j] / diff_z * (Qp[i][j] - Qp[i][j - 1])) * diff_t + absorb_coefy * absorb_coefz2 * Qv[i][j] * diff_t;//
					temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
						Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) / absorb_coefz4 + temp_w[i][j] + absorb_coefz2*Qw[i][j];
				}
				else if (jn<JMAX1){
					absorb_coefz4 = 1.0 + 25.0*pow((PML_z1 - Z[i][j]) / PML_Height1, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
					temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
					v_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * absorb_coefz1 * Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] / diff_y * (Qp[i][j] - Qp[i - 1][j]) -
						Z_over_y[i][j] / diff_z * (Qp[i][j] - Qp[i][j - 1])) * diff_t + absorb_coefy * absorb_coefz1 * Qv[i][j] * diff_t;//
					temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
						Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) / absorb_coefz4 + temp_w[i][j] + absorb_coefz1*Qw[i][j];
				}
				w_nn[i][j] -= absorb_coefy*temp_y*diff_t;
			}
			if (in >= IMAX - IMAX02){//

				absorb_coefy = PML_AbsorbYmax*pow((Y[i][j] - PML_y2) / PML_Width2, PML_alphay);
				absorb_coefy1 = 1.0 + 25.0*pow((Y[i][j] - PML_y2) / PML_Width2, PML_alphay);
				absorb_coefz1 = PML_AbsorbZmax*pow((PML_z1 - Z[i][j]) / PML_Height1, PML_alpha);
				absorb_coefz2 = PML_AbsorbZmax*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);

				temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
				v_nn[i][j] -= absorb_coefy*temp_y*diff_t
					+ (1.0 / absorb_coefy1 - 1.0)*
					Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] * diff_t / diff_y*(temp_p[i][j] - temp_p[i - 1][j]) -
					Z_over_y[i][j] * diff_t / diff_z*(temp_p[i][j] - temp_p[i][j - 1]));

				temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
					Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) + temp_w[i][j];

				if (jn >= JMAX - JMAX2 && jn<JMAX){
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((Y[i][j] - PML_y2) / PML_Width2, PML_alphay);
					temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
					v_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * absorb_coefz2 * Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] / diff_y * (Qp[i][j] - Qp[i - 1][j]) -
						Z_over_y[i][j] / diff_z * (Qp[i][j] - Qp[i][j - 1])) * diff_t + absorb_coefy * absorb_coefz2 * Qv[i][j] * diff_t;//
					temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
						Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) / absorb_coefz4 + temp_w[i][j] + absorb_coefz2*Qw[i][j];
				}
				else if (jn<JMAX1){
					absorb_coefz4 = 1.0 + 25.0*pow((PML_z1 - Z[i][j]) / PML_Height1, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y2 - Y[i][j]) / PML_Width2, PML_alphay);
					temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
					v_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * absorb_coefz1 * Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] / diff_y * (Qp[i][j] - Qp[i - 1][j]) -
						Z_over_y[i][j] / diff_z * (Qp[i][j] - Qp[i][j - 1])) * diff_t + absorb_coefy * absorb_coefz1 * Qv[i][j] * diff_t;//
					temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
						Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) / absorb_coefz4 + temp_w[i][j] + absorb_coefz1*Qw[i][j];
				}
				w_nn[i][j] -= absorb_coefy*temp_y*diff_t;
			}
		}
	}
	if (cr_judge!=8 && cr_judge!=0) 
	{
		for(i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){
			for(j=jj1;j<mpi_JMAX[mpi_jindex]-jj2;j++){
				if (cr_uc[i][j]==1){
					v_nn[i][j]=-1.0/(1.0+crpo_Beta_a*Aav[j])*(diff_t/crpo_eff_density_a/Jacobian[i][j]*Aav[j]
						*(Z_over_z[i][j]/diff_y*(temp_p[i][j]-temp_p[i-1][j])-
						Z_over_y[i][j]/diff_z*(temp_p[i][j]-temp_p[i][j-1])))-temp_v[i][j]
						+temp_v[i][j]/(1.0+crpo_Beta_a*Aav[j]);
					w_nn[i][j]=-1.0/(1.0+crpo_Beta_a*Aav[j])*(diff_t/crpo_eff_density_a/Jacobian[i][j]*Aav[j]
						*(Y_over_y[i][j]/diff_z*(temp_p[i][j]-temp_p[i][j-1])-
						Y_over_z[i][j]/diff_y*(temp_p[i][j]-temp_p[i-1][j])))-temp_w[i][j]
						+temp_w[i][j]/(1.0+crpo_Beta_a*Aav[j]);
				}
				if (cr_uc[i][j] == 2){
					v_nn[i][j] = -1.0 / (1.0 + crpo_Beta_b*Aav[j])*(diff_t / crpo_eff_density_b / Jacobian[i][j] * Aav[j]
						* (Z_over_z[i][j] / diff_y*(temp_p[i][j] - temp_p[i - 1][j]) -
						Z_over_y[i][j] / diff_z*(temp_p[i][j] - temp_p[i][j - 1]))) - temp_v[i][j]
						+ temp_v[i][j] / (1.0 + crpo_Beta_b*Aav[j]);
					w_nn[i][j] = -1.0 / (1.0 + crpo_Beta_b*Aav[j])*(diff_t / crpo_eff_density_b / Jacobian[i][j] * Aav[j]
						* (Y_over_y[i][j] / diff_z*(temp_p[i][j] - temp_p[i][j - 1]) -
						Y_over_z[i][j] / diff_y*(temp_p[i][j] - temp_p[i - 1][j]))) - temp_w[i][j]
						+ temp_w[i][j] / (1.0 + crpo_Beta_b*Aav[j]);
				}
			}
		}
	}
	else if (cr_judge==8){
		double Aav1,speed_sound1;
		Aav1=1.0/1.0e10;
		for(i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){
			for(j=jj1;j<mpi_JMAX[mpi_jindex]-jj2;j++){
				if (cr_uc[i][j]==1){
					v_nn[i][j]=-Aav1/Jacobian[i][j]*(Z_over_z[i][j]*diff_t/diff_y*(temp_p[i][j]-temp_p[i-1][j])-
							Z_over_y[i][j]*diff_t/diff_z*(temp_p[i][j]-temp_p[i][j-1]));
					w_nn[i][j]=-Aav1/Jacobian[i][j]*(Y_over_y[i][j]*diff_t/diff_z*(temp_p[i][j]-temp_p[i][j-1])-
							Y_over_z[i][j]*diff_t/diff_y*(temp_p[i][j]-temp_p[i-1][j]));
				}
			}
		}
	}
}
//------------------------------calculate pressure (p) for new time(nn)------------//
void air::cal_fp(double **temp_v,double **temp_w,
			double **temp_p)
{
	int i,in,j,jn,k,ii,jj;
	double absorb_coefz1, absorb_coefz2, temp_z, absorb_coefz3, absorb_coefz4;
	//  PML boundary for y direction
	double absorb_coefy, absorb_coefy1, temp_y;
	double cr_height,cr_length;
	double k_slope1,k_slope2,k_slope3,k_slope4,k_slope5,k_slope6;
// specify the initial condition for circle
	double time_length,k_slope;
	double vav_a,vav_b,wav_a,wav_b,vav_qpa,vav_qpb;

	IMAX1=IMAX1+1;
	/*if (IMAX2>1)
	{
		PML_Width2=fabs(whole_Y[IMAX2-1][0]-whole_Y[0][0]);
		PML_y2=Y_start;
	}*/

	if (cr_judge!=0){
	cr_height=-999.0;
		switch(cr_judge){
		case 1: case 8:
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){
						cr_length=pow((Y[i][j]-cr_c),2)+pow((Z[i][j]-cr_d),2);
						if (sqrt(cr_length)<cr_b)cr_uc[i][j]=1;
					}
				}
			}
			break;
		case 2:
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){
						cr_height=cr_a*pow(Y[i][j],3)+cr_b*pow(Y[i][j],2)+
										cr_c*Y[i][j]+cr_d;
						cr_height=sqrt(cr_height)-4.0;
						if (cr_height>Z[i][j])cr_uc[i][j]=1;
					}
				}
			}
			break;
		case 3:// triangular
			k_slope1=cr_c/(cr_b-cr_a);
			k_slope2=cr_c/(cr_b-cr_d);
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){
						if (Y[i][j]>=cr_a && Y[i][j]<=cr_b){
							cr_height=k_slope1*(Y[i][j]-cr_a);
							if (cr_height>=Z[i][j])cr_uc[i][j]=1;
						}else if(Y[i][j]>=cr_b && Y[i][j]<=cr_d){
							cr_height=k_slope2*(Y[i][j]-cr_d);
							if (cr_height>=Z[i][j])cr_uc[i][j]=1;
						}
					}
				}
			}
			break;
		case 4:// rectangular
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){
							if (Y[i][j]>=cr_a && Y[i][j]<=cr_b &&Z[i][j] >= cr_c &&Z[i][j] <= cr_d) 
								cr_uc[i][j] = 1;
					}
				}
			}
			break;
		case 5:// two wedges
			k_slope1=cr_c/(cr_b-cr_a);
			k_slope2=cr_c/(cr_b-cr_d);
			k_slope3=cr_g/(cr_f-cr_e);
			k_slope4=cr_g/(cr_f-cr_h);
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){

						if (Y[i][j]>=cr_a && Y[i][j]<=cr_b){
							cr_height=k_slope1*(Y[i][j]-cr_a);
							if (cr_height>=Z[i][j])cr_uc[i][j]=1;
						}else if(Y[i][j]>=cr_b && Y[i][j]<=cr_d){
							cr_height=k_slope2*(Y[i][j]-cr_d);
							if (cr_height>=Z[i][j])cr_uc[i][j]=1;
						}else if(Y[i][j]>=cr_e && Y[i][j]<=cr_f){
							cr_height=k_slope3*(Y[i][j]-cr_e);
							if (cr_height>=Z[i][j])cr_uc[i][j]=1;
						}else if(Y[i][j]>=cr_f && Y[i][j]<=cr_h){
							cr_height=k_slope4*(Y[i][j]-cr_h);
							if (cr_height>=Z[i][j])cr_uc[i][j]=1;
						}
					}
				}
			}
			break;
		case 6: // two rectangulars
//	        cr_height=cr_c;
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){
						if (Y[i][j]>=cr_b && Y[i][j]<=cr_d &&
								Z[i][j]<=cr_c)cr_uc[i][j]=1;
						else if (Y[i][j]>=cr_f && Y[i][j]<=cr_h &&
								Z[i][j]<=cr_g)cr_uc[i][j]=1;
					}
				}
			}
			break;
		case 7:
			k_slope1=cr_c/(cr_b-cr_a);
			k_slope2=cr_c/(cr_b-cr_d);
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){
						if (Y[i][j]>=cr_a && Y[i][j]<=cr_b){
							cr_height=k_slope1*(Y[i][j]-cr_a);
							if (cr_height>=Z[i][j])cr_uc[i][j]=1;
						}else if(Y[i][j]>=cr_b && Y[i][j]<=cr_d){
							cr_height=k_slope2*(Y[i][j]-cr_d);
							if (cr_height>=Z[i][j])cr_uc[i][j]=1;
						}else if (Y[i][j]>=cr_f && Y[i][j]<=cr_h &&
							Z[i][j]<=cr_g)cr_uc[i][j]=1;

					}
				}
			}
			break;
		case 9: //water and air
			cr_height=cr_a;
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){
						if(Z[i][j]<=cr_height)cr_uc[i][j]=1;
					}
				}
			}
			break;
		case 10: //Three rectangulars
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){
						if (Y[i][j]>=cr_b && Y[i][j]<=cr_d &&
								Z[i][j]<=cr_c)cr_uc[i][j]=1;
						else if (Y[i][j]>=cr_f && Y[i][j]<=cr_h &&
								Z[i][j]<=cr_g)cr_uc[i][j]=1;
						else if (Y[i][j]>=cr_j && Y[i][j]<=cr_l &&
								Z[i][j]<=cr_k)cr_uc[i][j]=1;
					}
				}
			}
			break;
		case 11: //One rectangular and circular
			cr_height=cr_c;
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){
							if (Y[i][j]>=cr_b && Y[i][j]<=cr_d
								&&Z[i][j]<=cr_height)cr_uc[i][j]=1;
								cr_length=pow((Y[i][j]-cr_g),2)+pow((Z[i][j]-cr_h),2);
								if (sqrt(cr_length)<cr_f)cr_uc[i][j]=1;
					}
				}
			}
			break;
		case 12: //initial condition on the circle
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){
						cr_length=pow((Y[i][j]-cr_c),2)+pow((Z[i][j]-cr_d),2);
						if (sqrt(cr_length<=cr_b))cr_uc[i][j]=1;
					}
				}
			}
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					if (cr_uc[i][j]==1){
						if (cr_uc[i+1][j]==0|| cr_uc[i][j-1]==0||cr_uc[i][j+1]==0){
							cr_uc[i][j]=2;
						}
					}
				}
			}
			break;
		case 13: // two cycles
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					cr_uc[i][j]=0;
					if (Y[i][j]>=ystar-1.0e-6&&Y[i][j]<=ystop+1.0e-6
						&& Z[i][j]>=zstar-1.0e-6&& Z[i][j]<=zstop+1.0e-6){
						cr_length=pow((Y[i][j]-cr_c),2)+pow((Z[i][j]-cr_d),2);
						if (sqrt(cr_length)<cr_b)cr_uc[i][j]=1;
						cr_length=pow((Y[i][j]-cr_g),2)+pow((Z[i][j]-cr_h),2);
						if (sqrt(cr_length)<cr_f)cr_uc[i][j]=1;
					}
				}
			}
			break;
		case 14: // three cylces
			for (i = 0; i<mpi_IMAX[mpi_iindex]; i++){
				for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
					cr_uc[i][j] = 0;
					if (Y[i][j] >= ystar - 1.0e-6&&Y[i][j] <= ystop + 1.0e-6
						&& Z[i][j] >= zstar - 1.0e-6&& Z[i][j] <= zstop + 1.0e-6){
						cr_length = pow((Y[i][j] - cr_c), 2) + pow((Z[i][j] - cr_d), 2);
						if (sqrt(cr_length)<cr_b)cr_uc[i][j] = 1;
						cr_length = pow((Y[i][j] - cr_g), 2) + pow((Z[i][j] - cr_h), 2);
						if (sqrt(cr_length)<cr_f)cr_uc[i][j] = 1;
						cr_length = pow((Y[i][j] - cr_k), 2) + pow((Z[i][j] - cr_l), 2);
						if (sqrt(cr_length)<cr_j)cr_uc[i][j] = 1;
					}
				}
			}
			break;		
		case 15: // cylces array
				for (i = 0; i<mpi_IMAX[mpi_iindex]; i++){
					for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
						cr_uc[i][j] = 0;
						if (Y[i][j] >= ystar - 1.0e-6 &&Y[i][j] <= ystop + 1.0e-6 && Z[i][j] >= zstar - 1.0e-6&& Z[i][j] <= zstop + 1.0e-6){
							for (ii = 0; ii < cr_f; ii++){ // y 
								for (jj = 0; jj < cr_g; jj++){ // z
									cr_length = pow((Y[i][j] - (cr_b + ii*cr_d)), 2) + pow((Z[i][j] - (cr_c + jj*cr_e)), 2);
									if (sqrt(cr_length) < cr_a) cr_uc[i][j] = 1;
								}
							}
						}
					}
				}
				break;
		case 16: // TDBC floor + wall			
				for (i = 0; i<mpi_IMAX[mpi_iindex]; i++){
					for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
						cr_uc[i][j] = 0;
						if (Y[i][j] >= ystar - 1.0e-6 &&Y[i][j] <= ystop + 1.0e-6 && Z[i][j] >= zstar - 1.0e-6&& Z[i][j] <= zstop + 1.0e-6){
							if (Y[i][j]>=cr_a && Y[i][j]<=cr_b &&Z[i][j] >= cr_c &&Z[i][j] <= cr_d) cr_uc[i][j] = 1;
							if (Y[i][j]>=cr_e && Y[i][j]<=cr_f &&Z[i][j] >= cr_g &&Z[i][j] <= cr_h) cr_uc[i][j] = 2;
						}
					}
				}
			break;
		case 17:// roughness with trangulars
			k_slope1 = cr_c / (cr_b - cr_a);
			k_slope2 = cr_c / (cr_b - cr_d);
			for (i = 0; i<mpi_IMAX[mpi_iindex]; i++){
				for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
					cr_uc[i][j] = 0;
					for (k = 0; k<cr_points; k++){
						if (Y[i][j] >= ystar - 1.0e-6 && Y[i][j] <= ystop + 1.0e-6
							&& Z[i][j] >= zstar - 1.0e-6 && Z[i][j] <= zstop + 1.0e-6){
							if (Y[i][j] >= cr_a + k*cr_e && Y[i][j] <= cr_b + k*cr_e){
								cr_height = k_slope1*(Y[i][j] - cr_b - k*cr_e) + 1.1;
								if (cr_height >= Z[i][j]) cr_uc[i][j] = 1;
							}
							else if (Y[i][j] >= cr_b + k*cr_e && Y[i][j] <= cr_d + k*cr_e){
								cr_height = k_slope2*(Y[i][j] - cr_b - k*cr_e) + 1.1;
								if (cr_height >= Z[i][j]) cr_uc[i][j] = 1;
							}
						}
					}
				}
			}
			cr_height = cr_f;
			for (i = 0; i<mpi_IMAX[mpi_iindex]; i++){
				for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
					if (Y[i][j] >= ystar - 1.0e-6 && Y[i][j] <= ystop + 1.0e-6
						&& Z[i][j] >= zstar - 1.0e-6 && Z[i][j] <= zstop + 1.0e-6){
						if (Y[i][j] >= cr_g && Y[i][j] <= cr_h
							&& Z[i][j] <= cr_height) cr_uc[i][j] = 1;
					}
				}
			}
			break;
		default:
			cout<<"it is wrong to entrance cr_judge";
		}
	}

	for (i = ii1; i < mpi_IMAX[mpi_iindex] - ii2; i++){
		for (j = jj1; j < mpi_JMAX[mpi_jindex] - jj2; j++){
			in = mpi_i1[mpi_iindex] + i;
			jn = mpi_j1[mpi_jindex] + j;
			//if(mpi_rank == 0 && i == ii1 && j == jj1) cout << ", jj1 = " << jj1 << ", ii1 = " << ii1;
			if (Vav[i][j] <= 0.0)
			{
				vav_a = temp_p[i + 1][j] - temp_p[i][j];
				vav_b = temp_p[i][j + 1] - temp_p[i][j];
				vav_qpa = Qp[i + 1][j] - Qp[i][j];
				vav_qpb = Qp[i][j + 1] - Qp[i][j];
			}
			else
			{
				vav_a = temp_p[i][j] - temp_p[i - 1][j];
				vav_b = temp_p[i][j] - temp_p[i][j - 1];
				vav_qpa = Qp[i][j] - Qp[i - 1][j];
				vav_qpb = Qp[i][j] - Qp[i][j - 1];
			}
			if (Wav[i][j] <= 0.0)
			{
				wav_a = temp_p[i + 1][j] - temp_p[i][j];
				wav_b = temp_p[i][j + 1] - temp_p[i][j];
			}
			else
			{
				wav_a = temp_p[i][j] - temp_p[i - 1][j];
				wav_b = temp_p[i][j] - temp_p[i][j - 1];
			}
			if (n_order == 2 && Y[i][j] >= ystar_vor - 1.0e-6 && Y[i][j] <= ystop_vor + 1.0e-6
				&& Z[i][j] >= zstar_vor - 1.0e-6 && Z[i][j] <= zstop_vor + 1.0e-6){
				if (Vav[i][j] <= 0.0){
					vav_a = (-1.0*temp_p[i + 2][j] + 4.0*temp_p[i + 1][j] - 3.0*temp_p[i][j])*0.5;
					vav_b = (-1.0*temp_p[i][j + 2] + 4.0*temp_p[i][j + 1] - 3.0*temp_p[i][j])*0.5;
					vav_qpa = (-1.0*Qp[i + 2][j] + 4.0*Qp[i + 1][j] - 3.0*Qp[i][j])*0.5;
					vav_qpb = (-1.0*Qp[i][j + 2] + 4.0*Qp[i][j + 1] - 3.0*Qp[i][j])*0.5;
				}
				else{
					vav_a = (3.0*temp_p[i][j] - 4.0*temp_p[i - 1][j] + temp_p[i - 2][j])*0.5;
					vav_b = (3.0*temp_p[i][j] - 4.0*temp_p[i][j - 1] + temp_p[i][j - 2])*0.5;
					vav_qpa = (3.0*Qp[i][j] - 4.0*Qp[i - 1][j] + Qp[i - 2][j])*0.5;
					vav_qpb = (3.0*Qp[i][j] - 4.0*Qp[i][j - 1] + Qp[i][j - 2])*0.5;
				}
				if (Wav[i][j] <= 0.0){
					wav_a = (-1.0*temp_p[i + 2][j] + 4.0*temp_p[i + 1][j] - 3.0*temp_p[i][j])*0.5;
					wav_b = (-1.0*temp_p[i][j + 2] + 4.0*temp_p[i][j + 1] - 3.0*temp_p[i][j])*0.5;
				}
				else{
					wav_a = (3.0*temp_p[i][j] - 4.0*temp_p[i - 1][j] + temp_p[i - 2][j])*0.5;
					wav_b = (3.0*temp_p[i][j] - 4.0*temp_p[i][j - 1] + temp_p[i][j - 2])*0.5;
				}
			}
			else if (n_order == 3 && Y[i][j] >= ystar_vor - 1.0e-6 && Y[i][j] <= ystop_vor + 1.0e-6
				&& Z[i][j] >= zstar_vor - 1.0e-6&& Z[i][j] <= zstop_vor + 1.0e-6){
				if (Vav[i][j] <= 0.0){
					vav_a = (-1.0*temp_p[i + 2][j] + 6.0*temp_p[i + 1][j] - 3.0*temp_p[i][j] - 2.0*temp_p[i - 1][j])*0.166667;
					vav_b = (-1.0*temp_p[i][j + 2] + 6.0*temp_p[i][j + 1] - 3.0*temp_p[i][j] - 2.0*temp_p[i][j - 1])*0.166667;
					vav_qpa = (-1.0*Qp[i + 2][j] + 6.0*Qp[i + 1][j] - 3.0*Qp[i][j] - 2.0*Qp[i - 1][j])*0.166667;
					vav_qpb = (-1.0*Qp[i][j + 2] + 6.0*Qp[i][j + 1] - 3.0*Qp[i][j] - 2.0*Qp[i][j - 1])*0.166667;
				}
				else{
					vav_a = (2.0*temp_p[i + 1][j] + 3.0*temp_p[i][j] - 6.0*temp_p[i - 1][j] + temp_p[i - 2][j])*0.166667;
					vav_b = (2.0*temp_p[i][j + 1] + 3.0*temp_p[i][j] - 6.0*temp_p[i][j - 1] + temp_p[i][j - 2])*0.166667;
					vav_qpa = (2.0*Qp[i + 1][j] + 3.0*Qp[i][j] - 6.0*Qp[i - 1][j] + Qp[i - 2][j])*0.166667;
					vav_qpb = (2.0*Qp[i][j + 1] + 3.0*Qp[i][j] - 6.0*Qp[i][j - 1] + Qp[i][j - 2])*0.166667;
				}
				if (Wav[i][j] <= 0.0){
					wav_a = (-1.0*temp_p[i + 2][j] + 6.0*temp_p[i + 1][j] - 3.0*temp_p[i][j] - 2.0*temp_p[i - 1][j])*0.166667;
					wav_b = (-1.0*temp_p[i][j + 2] + 6.0*temp_p[i][j + 1] - 3.0*temp_p[i][j] - 2.0*temp_p[i][j - 1])*0.166667;
				}
				else{
					wav_a = (2.0*temp_p[i + 1][j] + 3.0*temp_p[i][j] - 6.0*temp_p[i - 1][j] + temp_p[i - 2][j])*0.166667;
					wav_b = (2.0*temp_p[i][j + 1] + 3.0*temp_p[i][j] - 6.0*temp_p[i][j - 1] + temp_p[i][j - 2])*0.166667;
				}
			}
			else if (n_order == 5 && Y[i][j] >= ystar_vor - 1.0e-6 && Y[i][j] <= ystop_vor + 1.0e-6
				&& Z[i][j] >= zstar_vor - 1.0e-6&& Z[i][j] <= zstop_vor + 1.0e-6)
			{
				if (Vav[i][j] <= 0.0)
				{
					vav_a = -1.0*weno_scheme(temp_p[i + 3][j], temp_p[i + 2][j], temp_p[i + 1][j],
						temp_p[i][j], temp_p[i - 1][j], temp_p[i - 2][j], diff_y);
					vav_b = -1.0*weno_scheme(temp_p[i][j + 3], temp_p[i][j + 2], temp_p[i][j + 1],
						temp_p[i][j], temp_p[i][j - 1], temp_p[i][j - 2], diff_z);
					vav_qpa = -1.0*weno_scheme(Qp[i + 3][j], Qp[i + 2][j], Qp[i + 1][j],
						Qp[i][j], Qp[i - 1][j], Qp[i - 2][j], diff_y);
					vav_qpb = -1.0*weno_scheme(Qp[i][j + 3], Qp[i][j + 2], Qp[i][j + 1],
						Qp[i][j], Qp[i][j - 1], Qp[i][j - 2], diff_z);
				}
				else
				{
					vav_a = weno_scheme(temp_p[i - 3][j], temp_p[i - 2][j], temp_p[i - 1][j],
						temp_p[i][j], temp_p[i + 1][j], temp_p[i + 2][j], diff_y);
					vav_b = weno_scheme(temp_p[i][j - 3], temp_p[i][j - 2], temp_p[i][j - 1],
						temp_p[i][j], temp_p[i][j + 1], temp_p[i][j + 2], diff_z);
					vav_qpa = weno_scheme(Qp[i - 3][j], Qp[i - 2][j], Qp[i - 1][j],
						Qp[i][j], Qp[i + 1][j], Qp[i + 2][j], diff_y);
					vav_qpb = weno_scheme(Qp[i][j - 3], Qp[i][j - 2], Qp[i][j - 1],
						Qp[i][j], Qp[i][j + 1], Qp[i][j + 2], diff_z);
				}
				if (Wav[i][j] <= 0.0)
				{
					wav_a = -1.0*weno_scheme(temp_p[i + 3][j], temp_p[i + 2][j], temp_p[i + 1][j],
						temp_p[i][j], temp_p[i - 1][j], temp_p[i - 2][j], diff_y);
					wav_b = -1.0*weno_scheme(temp_p[i][j + 3], temp_p[i][j + 2], temp_p[i][j + 1],
						temp_p[i][j], temp_p[i][j - 1], temp_p[i][j - 2], diff_z);
				}
				else
				{
					wav_a = weno_scheme(temp_p[i - 3][j], temp_p[i - 2][j], temp_p[i - 1][j],
						temp_p[i][j], temp_p[i + 1][j], temp_p[i + 2][j], diff_y);
					wav_b = weno_scheme(temp_p[i][j - 3], temp_p[i][j - 2], temp_p[i][j - 1],
						temp_p[i][j], temp_p[i][j + 1], temp_p[i][j + 2], diff_z);

				}
			}
			//if (mpi_rank == 0 && i == 876 && j == 1) cout << ", temp_w[i][j] = " << temp_w[i][j] << ", temp_w[i][j+1] = " << temp_w[i][j+1];
			p_nn[i][j] = -Vav[i][j] * (Z_over_z[i][j] * diff_t / diff_y*vav_a -
				Z_over_y[i][j] * diff_t / diff_z*vav_b) / Jacobian[i][j]
				- Wav[i][j] * (Y_over_y[i][j] * diff_t / diff_z*wav_b -
				Y_over_z[i][j] * diff_t / diff_y*wav_a) / Jacobian[i][j]
				- adiabatic_coef*Pav / Jacobian[i][j] *
				(Z_over_z[i][j] * diff_t / diff_y*(temp_v[i + 1][j] - temp_v[i][j]) -
				Z_over_y[i][j] * diff_t / diff_z*(temp_v[i][j + 1] - temp_v[i][j]) +
				Y_over_y[i][j] * diff_t / diff_z*(temp_w[i][j + 1] - temp_w[i][j]) -
				Y_over_z[i][j] * diff_t / diff_y*(temp_w[i + 1][j] - temp_w[i][j])) -
				adiabatic_coef*temp_p[i][j] * diff_t*(Vav_over_y[i][j] + Wav_over_z[i][j]);
			if (jn < JMAX1)
			{
				absorb_coefz1=PML_AbsorbZmax*pow((PML_z1-Z[i][j])/PML_Height1,PML_alpha);
				absorb_coefz3=1.0+25.0*pow((PML_z1-Z[i][j])/PML_Height1,PML_alpha);
				temp_z=adiabatic_coef*Pav*((Qv[i+1][j]-Qv[i][j])/diff_y*Z_over_z[i][j]-
				(Qv[i][j+1]-Qv[i][j])/diff_z*Z_over_y[i][j]);
				temp_z += Vav[i][j]*(vav_qpa/diff_y*Z_over_z[i][j]-
				vav_qpb/diff_z*Z_over_y[i][j]);
				p_nn[i][j] -= absorb_coefz1*(temp_z/Jacobian[i][j]+temp_p[i][j])*diff_t+(1.0/absorb_coefz3-1.0)*(Y_over_y[i][j]*diff_t/diff_z
				*(temp_w[i][j+1]-temp_w[i][j])-Y_over_z[i][j]*diff_t/diff_y*(temp_w[i+1][j]-temp_w[i][j]))
				*adiabatic_coef*Pav/Jacobian[i][j];
				/*p_nn[i][j] = -crpo_Gama / Jacobian[i][j] *
					(Z_over_z[i][j] / diff_y*(temp_v[i + 1][j] - temp_v[i][j]) -
					Z_over_y[i][j] / diff_z*(temp_v[i][j + 1] - temp_v[i][j])
					+ Y_over_y[i][j] / diff_z*(temp_w[i][j + 1] - temp_w[i][j]) -
					Y_over_z[i][j] / diff_y*(temp_w[i + 1][j] - temp_w[i][j]));*/
			}
			if (jn >= JMAX - JMAX2)
			{
				absorb_coefz2 = PML_AbsorbZmax*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);
				absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);
				temp_z = adiabatic_coef*Pav*((Qv[i + 1][j] - Qv[i][j]) / diff_y*Z_over_z[i][j] -
					(Qv[i][j + 1] - Qv[i][j]) / diff_z*Z_over_y[i][j]);
				temp_z += Vav[i][j] * (vav_qpa / diff_y*Z_over_z[i][j] -
					vav_qpb / diff_z*Z_over_y[i][j]);
				p_nn[i][j] -= absorb_coefz2*(temp_z / Jacobian[i][j] + temp_p[i][j])*diff_t + (1.0 / absorb_coefz4 - 1.0)*(Y_over_y[i][j] * diff_t / diff_z
					*(temp_w[i][j + 1] - temp_w[i][j]) - Y_over_z[i][j] * diff_t / diff_y*(temp_w[i + 1][j] - temp_w[i][j]))
					*adiabatic_coef*Pav / Jacobian[i][j];

			}
			if (in < IMAX01 && IMAX1 >= 9){

				absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
				absorb_coefy1 = 1.0 + 25.0*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
				absorb_coefz1 = PML_AbsorbZmax*pow((PML_z1 - Z[i][j]) / PML_Height1, PML_alpha);
				absorb_coefz2 = PML_AbsorbZmax*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);

				temp_y = 1.0 / Jacobian[i][j] * adiabatic_coef*Pav*(Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
					Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) + temp_p[i][j];

				if (jn >= JMAX - JMAX2 && jn < JMAX)
				{
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
					temp_y = 1.0 / Jacobian[i][j] * adiabatic_coef*Pav*(Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
						Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) / absorb_coefz4 + temp_p[i][j];
				}
				else if (jn < JMAX1)
				{
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z1) / PML_Height1, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
					temp_y = 1.0 / Jacobian[i][j] * adiabatic_coef*Pav*(Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
						Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) / absorb_coefz4 + temp_p[i][j];
				}
				p_nn[i][j] -= absorb_coefy*temp_y*diff_t
					+ (1.0 / absorb_coefy1 - 1.0)*(Z_over_z[i][j] * diff_t / diff_y*(temp_v[i + 1][j] - temp_v[i][j]) - Z_over_y[i][j] * diff_t / diff_z*(temp_v[i][j + 1] - temp_v[i][j]))*adiabatic_coef*Pav / Jacobian[i][j];
				if (jn >= JMAX - JMAX2 && jn < JMAX)
				{
					p_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * adiabatic_coef * Pav * absorb_coefz2 / Jacobian[i][j] * (Z_over_z[i][j] / diff_y*(Qv[i + 1][j] - Qv[i][j]) 
						- Z_over_y[i][j] / diff_z*(Qv[i][j + 1] - Qv[i][j])) * diff_t + absorb_coefy * absorb_coefz2 * Qp[i][j] * diff_t;
				}
				else if (jn < JMAX1)
				{
					p_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * adiabatic_coef * Pav * absorb_coefz1 / Jacobian[i][j] * (Z_over_z[i][j] / diff_y*(Qv[i + 1][j] - Qv[i][j])
						- Z_over_y[i][j] / diff_z*(Qv[i][j + 1] - Qv[i][j])) * diff_t + absorb_coefy*absorb_coefz1*Qp[i][j] * diff_t;
				}/**/
			}
			if (in >= IMAX - IMAX02){// && jn>=JMAX1

				absorb_coefy = PML_AbsorbYmax*pow((Y[i][j] - PML_y2) / PML_Width2, PML_alphay);
				absorb_coefy1 = 1.0 + 25.0*pow((Y[i][j] - PML_y2) / PML_Width2, PML_alphay);
				absorb_coefz1 = PML_AbsorbZmax*pow((PML_z1 - Z[i][j]) / PML_Height1, PML_alpha);
				absorb_coefz2 = PML_AbsorbZmax*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);

				temp_y = 1.0 / Jacobian[i][j] * adiabatic_coef*Pav*(Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
					Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) + temp_p[i][j];

				if (jn >= JMAX - JMAX2 && jn<JMAX)
				{
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y2 - Y[i][j]) / PML_Width2, PML_alphay);
					temp_y = 1.0 / Jacobian[i][j] * adiabatic_coef*Pav*(Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
						Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) / absorb_coefz4 + temp_p[i][j];
				}
				else if (jn<JMAX1)
				{
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z1) / PML_Height1, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y2 - Y[i][j]) / PML_Width2, PML_alphay);
					temp_y = 1.0 / Jacobian[i][j] * adiabatic_coef*Pav*(Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
						Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) / absorb_coefz4 + temp_p[i][j];
				}
				p_nn[i][j] -= absorb_coefy*temp_y*diff_t
					+ (1.0 / absorb_coefy1 - 1.0)*(Z_over_z[i][j] * diff_t / diff_y*(temp_v[i + 1][j] - temp_v[i][j]) - Z_over_y[i][j] * diff_t / diff_z*(temp_v[i][j + 1] - temp_v[i][j]))*adiabatic_coef*Pav / Jacobian[i][j];
				if (jn >= JMAX - JMAX2 && jn < JMAX)
				{
					p_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * adiabatic_coef * Pav * absorb_coefz2 / Jacobian[i][j] * (Z_over_z[i][j] / diff_y*(Qv[i + 1][j] - Qv[i][j])
						- Z_over_y[i][j] / diff_z*(Qv[i][j + 1] - Qv[i][j])) * diff_t + absorb_coefy*absorb_coefz2*Qp[i][j] * diff_t;
				}
				else if (jn < JMAX1)
				{
					p_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * adiabatic_coef * Pav * absorb_coefz1 / Jacobian[i][j] * (Z_over_z[i][j] / diff_y*(Qv[i + 1][j] - Qv[i][j])
						- Z_over_y[i][j] / diff_z*(Qv[i][j + 1] - Qv[i][j])) * diff_t + absorb_coefy*absorb_coefz1*Qp[i][j] * diff_t;
				}/**/
			}
		}
	}
	if (cr_judge!=8 && cr_judge!=0) {
		for(i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){
			for (j=jj1;j<mpi_JMAX[mpi_jindex]-jj2;j++){
				if (cr_uc[i][j] == 1 ){
					p_nn[i][j]=-crpo_Gama_a/Jacobian[i][j]*
						(Z_over_z[i][j]/diff_y*(temp_v[i+1][j]-temp_v[i][j])-
						Z_over_y[i][j]/diff_z*(temp_v[i][j+1]-temp_v[i][j])
						+Y_over_y[i][j]/diff_z*(temp_w[i][j+1]-temp_w[i][j])-
						Y_over_z[i][j]/diff_y*(temp_w[i+1][j]-temp_w[i][j]));
				}
				if (cr_uc[i][j] == 2 ){
					p_nn[i][j]=-crpo_Gama_b/Jacobian[i][j]*
						(Z_over_z[i][j]/diff_y*(temp_v[i+1][j]-temp_v[i][j])-
						Z_over_y[i][j]/diff_z*(temp_v[i][j+1]-temp_v[i][j])
						+Y_over_y[i][j]/diff_z*(temp_w[i][j+1]-temp_w[i][j])-
						Y_over_z[i][j]/diff_y*(temp_w[i+1][j]-temp_w[i][j]));
				}
			}
		}
	}
	else if (cr_judge==8){
		double Aav1,speed_sound1;
		Aav1=1.0/1.0e10;
		speed_sound1=340.0;
		for(i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){
			for (j=jj1;j<mpi_JMAX[mpi_jindex]-jj2;j++){
				if (cr_uc[i][j]==1){
					p_nn[i][j]=-speed_sound1*speed_sound1/Aav1/Jacobian[i][j]*
						(Z_over_z[i][j]*diff_t/diff_y*(temp_v[i+1][j]-temp_v[i][j])-
						Z_over_y[i][j]*diff_t/diff_z*(temp_v[i][j+1]-temp_v[i][j])+
						Y_over_y[i][j]*diff_t/diff_z*(temp_w[i][j+1]-temp_w[i][j])-
						Y_over_z[i][j]*diff_t/diff_y*(temp_w[i+1][j]-temp_w[i][j]));
				}
			}
		}
	}
}

double air::weno_scheme(double wen_u1,double wen_u2,double wen_u3,double wen_u4,double wen_u5,double wen_u6,double ddx)
{ 
	double q1,q2,q3,q4,q5;
	double wegh1,wegh2,wegh3,wegh4,wegh5;
	double s1,s2,s3,alpha1,alpha2,alpha3;
	double ux1,ux2,ux3,uxx;
	q1=(wen_u2-wen_u1)/ddx;
	q2=(wen_u3-wen_u2)/ddx;
	q3=(wen_u4-wen_u3)/ddx;
	q4=(wen_u5-wen_u4)/ddx;
	q5=(wen_u6-wen_u5)/ddx;
	s1=1.083333*(q1-2.0*q2+q3)*(q1-2.0*q2+q3)+
		0.25*(q1-4.0*q2+3.0*q3)*(q1-4.0*q2+3.0*q3);
	s2=1.083333*(q2-2.0*q3+q4)*(q2-2.0*q3+q4)+
		0.25*(q2-q4)*(q2-q4);
	s3=1.083333*(q3-2.0*q4+q5)*(q3-2.0*q4+q5)+
		0.25*(3.0*q3-4.0*q4+q5)*(3.0*q3-4.0*q4+q5);
	alpha1=0.1/(s1+1.0e-6); alpha2=0.6/(s2+1.0e-6); alpha3=0.3/(s3+1.0e-6);
	wegh1=alpha1/(alpha1+alpha2+alpha3);
	wegh2=alpha2/(alpha1+alpha2+alpha3);
	wegh3=alpha3/(alpha1+alpha2+alpha3);
	ux1=q1/3.0-7.0/6.0*q2+11.0/6.0*q3;
	ux2=-q2/6.0+5.0/6.0*q3+1.0/3.0*q4;
	ux3=q3/3.0+5.0/6.0*q4-1.0/6.0*q5;
	uxx=wegh1*ux1+wegh2*ux2+wegh3*ux3;
	return uxx*ddx;
}
//-----------------------------update boundary condtions of pressure--------------------//
void air::UpdateBC_pressure(boundary_location BC,int time_judge,int time_current)
{
	int i,j,k,in,jn,rank_send,rank_recv,tag;
	int no,nr1,nr2;
	double time_length,k_slope;
	// specify sender matrix in each block's boundary
	nr1=0;
	nr2=0;
	for(i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){
		if (mpi_jindex!=0){
			for(no=0;no<nno;no++){ 
				pss[nr1]=p_nn[i][no+nno];
				nr1=nr1+1;
			}
		}
		if (mpi_jindex!=mpi_zarea-1){
			for(no=0;no<nno;no++){ 
				jn=mpi_JMAX[mpi_jindex]-nno*2;
				pns[nr2]=p_nn[i][jn+no];
				nr2=nr2+1;
			}
		}
	}
	nr1=0;
	nr2=0;
	for (j=jj1;j<mpi_JMAX[mpi_jindex]-jj2;j++){
		if (mpi_iindex!=0){
			for(no=0;no<nno;no++){ 
				pws[nr1]=p_nn[no+nno][j];
				nr1=nr1+1;
			}
		}
		if (mpi_iindex!=mpi_yarea-1){
			for(no=0;no<nno;no++){ 
				in=mpi_IMAX[mpi_iindex]-nno*2;
				pes[nr2]=p_nn[in+no][j];
				nr2=nr2+1;
			}
		}
	}
	switch(BC){
	case WestBC:
	{
		if (mpi_iindex==0){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++) {
				p_nn[0][j]=p_nn[1][j];////left
			}
			// time series boundary condition
			if (jjmax2!=0 && time_judge==1){
				time_length=time_current*diff_t;
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					jn=mpi_j1[mpi_jindex]+j;
					if (jn>=jjmax1 && jn<jjmax2){
						p_nn[0][j] = sin(2 * 3.1415926 * 500 * time_length); // 2*3.14159*85 = 534.070742
						/*
						for (k=0;k<num_data-1;k++){
							if ((time_step[k]-1.e-6)<=time_length && (time_step[k+1]+1.e-6)>= time_length) {
								k_slope=(init_pressure[k+1]-init_pressure[k])/(time_step[k+1]-time_step[k]);
								p_nn[0][j]=k_slope*(time_length-time_step[k])+init_pressure[k];
								k=num_data;
							}
						}
						*/
					}
				}
			}
		}else{
			tag=1;
			if (mpi_porous!=0)rank_recv=(mpi_iindex-1)*mpi_porous+mpi_jindex;
			else rank_recv=(mpi_iindex-1)*mpi_zarea+mpi_jindex;
			MPI_Send(pws,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
			tag=2;
			rank_send=rank_recv;
			MPI_Recv(pwr,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
		}

	}
	break;
	case EastBC:
	{
		if (mpi_iindex==mpi_yarea-1){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				in=mpi_IMAX[mpi_iindex]-1;
				p_nn[in][j]=p_nn[in-1][j];//right
			}
		}else{
			tag=1;
			if (mpi_porous!=0)rank_send=(mpi_iindex+1)*mpi_porous+mpi_jindex;
			else rank_send=(mpi_iindex+1)*mpi_zarea+mpi_jindex;
			MPI_Recv(per,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
			tag=2;
			rank_recv=rank_send;
			MPI_Send(pes,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
		}
	}
	break;
	case SouthBC:
	{
		if (mpi_jindex==0){
			//cout << "mpi_jindex = " << mpi_jindex << ", mpi_iindex = " << mpi_iindex << endl;
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				p_nn[i][0]=p_nn[i][1]; //Bottom
				/*if (i == 876) 
					cout << ", w_nn[i][2] = " << w_nn[i][2] << ", p_nn[i][2] = " << p_nn[i][2];*/
			}
		}else{
			tag=3;
			if (mpi_porous!=0)rank_recv=mpi_iindex*mpi_porous+mpi_jindex-1;
			else rank_recv=mpi_iindex*mpi_zarea+mpi_jindex-1;
			MPI_Send(pss,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
			tag=4;
			rank_send=rank_recv;
			MPI_Recv(psr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
		}
	}
	break;
	case NorthBC://characteristic boundary condition
	{
		if (mpi_jindex==mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-1;
			for(i=0;i<mpi_IMAX[mpi_iindex];i++) p_nn[i][jn]=p_nn[i][jn-1];//upper
		}else{
			tag=3;
			if (mpi_porous!=0)rank_send=mpi_iindex*mpi_porous+mpi_jindex+1;
			else rank_send=mpi_iindex*mpi_zarea+mpi_jindex+1;
			MPI_Recv(pnr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
			tag=4;
			rank_recv=rank_send;
			MPI_Send(pns,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
		}
	}
	break;
	default:
		cout<<"it is wrong to entrance boundary";
	}

	nr1=0;
	nr2=0;
	for (i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){
		if(mpi_jindex!=0){
			for(no=0;no<nno;no++){
				p_nn[i][no]=psr[nr1];
				nr1=nr1+1;
			}
		}
		if(mpi_jindex!=mpi_zarea-1){
			for(no=0;no<nno;no++){
				jn=mpi_JMAX[mpi_jindex]-nno;
				p_nn[i][jn+no]=pnr[nr2];
				nr2=nr2+1;
			}
		}
	}
	nr1=0;
	nr2=0;
	for (j=jj1;j<mpi_JMAX[mpi_jindex]-jj2;j++){
		if (mpi_iindex!=0){
			for(no=0;no<nno;no++){
				p_nn[no][j]=pwr[nr1];
				nr1=nr1+1;
			}
		}
		if (mpi_iindex!=mpi_yarea-1){
			for(no=0;no<nno;no++){
				in=mpi_IMAX[mpi_iindex]-nno;
				p_nn[in+no][j]=per[nr2];
				nr2=nr2+1;
			}
		}
	}
}

//--------------------------update boundary conditions of velocity(V and W)------------------//
void air::UpdateBC_velocity(boundary_location BC)
{
	//v_nn[1][j] and w_nn[i][1] for Salomons scheme is real boundary value,
	//each time it should be reset.
	// specify sender matrix in each block's boundary
	int i,j,in,jn,rank_send,rank_recv,tag;
	int no,nr1,nr2;
	nr1=0;
	nr2=0;
	for(i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){
		if (mpi_jindex!=0) {
			for(no=0;no<nno;no++){
				vss[nr1]=v_nn[i][no+nno];
				wss[nr1]=w_nn[i][no+nno];
				nr1=nr1+1;
			}
		}
		if (mpi_jindex!=mpi_zarea-1){
			for(no=0;no<nno;no++){
				jn=mpi_JMAX[mpi_jindex]-nno*2;
				vns[nr2]=v_nn[i][jn+no];
				wns[nr2]=w_nn[i][jn+no];
				nr2=nr2+1;
			}
		}
	}
	nr1=0;
	nr2=0;
	for (j=jj1;j<mpi_JMAX[mpi_jindex]-jj2;j++){
		if (mpi_iindex!=0){
			for(no=0;no<nno;no++){
				vws[nr1]=v_nn[no+nno][j];
				wws[nr1]=w_nn[no+nno][j];
				nr1=nr1+1;
			}
		}
		if (mpi_iindex!=mpi_yarea-1){
			for(no=0;no<nno;no++){
				in=mpi_IMAX[mpi_iindex]-nno*2;
				ves[nr2]=v_nn[in+no][j];
				wes[nr2]=w_nn[in+no][j];
				nr2=nr2+1;
			}
		}
	}
	switch(BC){
	case WestBC:
	{
		if (mpi_iindex==0){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++) {
				v_nn[0][j]=0.0;w_nn[0][j]=w_nn[1][j];////left
			}
		}else{
			tag=5;
			if (mpi_porous!=0)rank_recv=(mpi_iindex-1)*mpi_porous+mpi_jindex;
			else rank_recv=(mpi_iindex-1)*mpi_zarea+mpi_jindex;
			MPI_Send(vws,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
			tag=6;
			MPI_Send(wws,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
			tag=7;
			rank_send=rank_recv;
			MPI_Recv(vwr,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
			tag=8;
			MPI_Recv(wwr,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
		}
	}
		break;
	case EastBC:
	{
		if (mpi_iindex==mpi_yarea-1){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				in=mpi_IMAX[mpi_iindex]-1;
				v_nn[in][j]=0.0;//right
				w_nn[in][j]=w_nn[in-1][j];
			}
		}else{
			tag=5;
			if (mpi_porous!=0)rank_send=(mpi_iindex+1)*mpi_porous+mpi_jindex;
			else rank_send=(mpi_iindex+1)*mpi_zarea+mpi_jindex;
			MPI_Recv(ver,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
			tag=6;
			MPI_Recv(wer,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
			tag=7;
			rank_recv=rank_send;
			MPI_Send(ves,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
			tag=8;
			MPI_Send(wes,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
		}
	}
		break;
	case SouthBC:
	{
		if (mpi_jindex==0){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				w_nn[i][0]=0.0;v_nn[i][0]=v_nn[i][1]; //Bottom
			}
		}else{
			tag=9;
			if (mpi_porous!=0)rank_recv=mpi_iindex*mpi_porous+mpi_jindex-1;
			else rank_recv=mpi_iindex*mpi_zarea+mpi_jindex-1;
			MPI_Send(vss,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
			tag=10;
			MPI_Send(wss,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
			tag=11;
			rank_send=rank_recv;
			MPI_Recv(vsr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
			tag=12;
			MPI_Recv(wsr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
		}
	}
		break;
	case NorthBC://characteristic boundary condition
	{
		if (mpi_jindex==mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-1;
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				w_nn[i][jn]=0.0;v_nn[i][jn]=v_nn[i][jn-1];//upper
			}
		}else{
			tag=9;
			if (mpi_porous!=0)rank_send=mpi_iindex*mpi_porous+mpi_jindex+1;
			else rank_send=mpi_iindex*mpi_zarea+mpi_jindex+1;
			MPI_Recv(vnr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
			tag=10;
			MPI_Recv(wnr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&status);
			tag=11;
			rank_recv=rank_send;
			MPI_Send(vns,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
			tag=12;
			MPI_Send(wns,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
		}
	}
	break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
	nr1=0;
	nr2=0;
	for (i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){
		if(mpi_jindex!=0) {
			for(no=0;no<nno;no++){
				v_nn[i][no]=vsr[nr1];
				w_nn[i][no]=wsr[nr1];
				nr1=nr1+1;
			}
		}
		if(mpi_jindex!=mpi_zarea-1){
			for(no=0;no<nno;no++){
				jn=mpi_JMAX[mpi_jindex]-nno;
				v_nn[i][jn+no]=vnr[nr2];
				w_nn[i][jn+no]=wnr[nr2];
				nr2=nr2+1;
			}
		}
	}

	nr1=0;
	nr2=0;
	for (j=jj1;j<mpi_JMAX[mpi_jindex]-jj2;j++){
		if (mpi_iindex!=0){
			for(no=0;no<nno;no++){
				v_nn[no][j]=vwr[nr1];
				w_nn[no][j]=wwr[nr1];
				nr1=nr1+1;
			}
		}
		if (mpi_iindex!=mpi_yarea-1){
			for(no=0;no<nno;no++){
				in=mpi_IMAX[mpi_iindex]-nno;
				v_nn[in+no][j]=ver[nr2];
				w_nn[in+no][j]=wer[nr2];
				nr2=nr2+1;
			}
		}
	}
}
void air::Update_PML_Qvw()
{
	int i, j, in, jn, jn1, i1, i2;
	int mpi_index, rank_num, rank_send, rank_recv, tag;

	for (i = 0; i<mpi_IMAX[mpi_iindex]; i++){
		for (j = 0; j < mpi_JMAX[mpi_jindex]; j++){
			in = mpi_i1[mpi_iindex] + i;
			jn = mpi_j1[mpi_jindex] + j;
			if (jn < JMAX1 || jn >= JMAX - JMAX2){//
				Qv[i][j] += diff_t*v_nn[i][j];
				Qw[i][j] += diff_t*w_nn[i][j];
			}
			if (in < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){//
				Qv[i][j] += diff_t*v_nn[i][j];
				Qw[i][j] += diff_t*w_nn[i][j];	
			}
			if (in >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1){
				Qv[i][j] += diff_t*v_nn[i][j];
				Qw[i][j] += diff_t*w_nn[i][j];
			}
		}
	}
	if (mpi_iindex != 0){
		rank_num = mpi_IMAX[mpi_iindex] * mpi_JMAX[mpi_jindex];
		rank_recv = mpi_jindex;

		if (mpi_iindex != 0)i1 = 1;
		else i1 = 0;
		if (mpi_iindex != mpi_yarea - 1)i2 = 1;
		else i2 = 0;
		mpi_index = 0;
		for (i = i1; i<mpi_IMAX[mpi_iindex] - i2; i++){
			for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
				in = mpi_i1[mpi_iindex] + i;
				jn = mpi_j1[mpi_jindex] + j;
				if (jn<JMAX1 ||jn >= JMAX - JMAX2){// 
					mpi_fs[mpi_index] = Qv[i][j];
					mpi_index = mpi_index + 1;
				}
				if (in < IMAX01 && jn <= JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
					mpi_fs[mpi_index] = Qv[i][j];
					mpi_index = mpi_index + 1;
				}
				if (in >= IMAX - IMAX02 && jn<=JMAX - JMAX2 && jn >= JMAX1){// && IMAX1 >= 0
					mpi_fs[mpi_index] = Qv[i][j];
					mpi_index = mpi_index + 1;
				}
			}
		}
		tag = 24;
		MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);

		mpi_index = 0;
		for (i = i1; i<mpi_IMAX[mpi_iindex] - i2; i++){
			for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
				in = mpi_i1[mpi_iindex] + i;
				jn = mpi_j1[mpi_jindex] + j;
				if (jn<JMAX1 ||jn >= JMAX - JMAX2){// 
					mpi_fs[mpi_index] = Qw[i][j];
					mpi_index = mpi_index + 1;
				}
				if (in < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
					mpi_fs[mpi_index] = Qw[i][j];
					mpi_index = mpi_index + 1;
				}
				if (in >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1){// && IMAX1 >= 0
					mpi_fs[mpi_index] = Qw[i][j];
					mpi_index = mpi_index + 1;
				}
			}
		}
		tag = 25;
		MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);

		rank_send = rank_recv;
		rank_num = IMAX*mpi_JMAX[mpi_jindex];
		tag = 26;
		MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		mpi_index = 0;
		for (i = 0; i<IMAX; i++){
			for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
				jn = mpi_j1[mpi_jindex] + j;
				if (jn<JMAX1 ||jn >= JMAX - JMAX2){// 
					whole_Qv[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
				if (i < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
					whole_Qv[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
				if (i >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1){// && IMAX1 >= 0
					whole_Qv[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
			}
		}

		tag = 27;
		MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		mpi_index = 0;
		for (i = 0; i<IMAX; i++){
			for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
				jn = mpi_j1[mpi_jindex] + j;
				if (jn<JMAX1 ||jn >= JMAX - JMAX2){// 
					whole_Qw[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
				if (i < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
					whole_Qw[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
				if (i >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1){// && IMAX1 >= 0
					whole_Qw[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
			}
		}
	}
	else{
		for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
			for (i = 0; i<mpi_IMAX[mpi_iindex]; i++){
				in = mpi_i1[mpi_iindex] + i;
				jn = mpi_j1[mpi_jindex] + j;
				if (jn<JMAX1 ||jn >= JMAX - JMAX2){// 
					whole_Qv[in][j] = Qv[i][j];
					whole_Qw[in][j] = Qw[i][j];
				}
				if (in < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
					whole_Qv[in][j] = Qv[i][j];
					whole_Qw[in][j] = Qw[i][j];
				}
				if (in >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1){// && IMAX1 >= 0
					whole_Qv[in][j] = Qv[i][j];
					whole_Qw[in][j] = Qw[i][j];
				}
			}
		}
		for (in = 1; in<mpi_yarea; in++){
			if (mpi_porous != 0)rank_send = (mpi_iindex + in)*mpi_porous + mpi_jindex;
			else rank_send = (mpi_iindex + in)*mpi_zarea + mpi_jindex;
			jn = mpi_jindex;
			rank_num = mpi_IMAX[in] * mpi_JMAX[jn];
			if (in != 0)i1 = 1;
			else i1 = 0;
			if (in != mpi_yarea - 1)i2 = 1;
			else i2 = 0;
			mpi_index = 0;
			tag = 24;
			MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			for (i = mpi_i1[in] + i1; i<mpi_i2[in] + 1 - i2; i++){
				for (j = 0; j<mpi_JMAX[jn]; j++){
					jn1 = mpi_j1[mpi_jindex] + j;
					if (jn1<JMAX1 || jn1 >= JMAX - JMAX2){//
						whole_Qv[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
					if (i < IMAX01 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
						whole_Qv[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX02 && jn1<JMAX - JMAX2 && jn1 >= JMAX1){// && IMAX1 >= 0
						whole_Qv[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
				}
			}
			mpi_index = 0;
			tag = 25;
			MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			for (i = mpi_i1[in] + i1; i<mpi_i2[in] + 1 - i2; i++){
				for (j = 0; j<mpi_JMAX[jn]; j++){
					jn1 = mpi_j1[mpi_jindex] + j;
					if (jn1<JMAX1 || jn1 >= JMAX - JMAX2){//
						whole_Qw[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
					if (i < IMAX01 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
						whole_Qw[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX02 && jn1<JMAX - JMAX2 && jn1 >= JMAX1){// && IMAX1 >= 0
						whole_Qw[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
				}
			}
		}
		for (in = 1; in<mpi_yarea; in++){
			jn = mpi_jindex;
			rank_num = IMAX*mpi_JMAX[jn];
			if (mpi_porous != 0) rank_recv = (mpi_iindex + in)*mpi_porous + mpi_jindex;
			else rank_recv = (mpi_iindex + in)*mpi_zarea + mpi_jindex;

			mpi_index = 0;
			for (i = 0; i<IMAX; i++){
				for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
					jn1 = mpi_j1[jn] + j;
					if (jn1<JMAX1 || jn1 >= JMAX - JMAX2){//
						mpi_fs[mpi_index] = whole_Qv[i][j];
						mpi_index = mpi_index + 1;
					}
					if (i < IMAX01 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
						mpi_fs[mpi_index] = whole_Qv[i][j];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX02 && jn1<JMAX - JMAX2 && jn1 >= JMAX1){// && IMAX1 >= 0
						mpi_fs[mpi_index] = whole_Qv[i][j];
						mpi_index = mpi_index + 1;
					}
				}
			}
			tag = 26;
			MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);

			mpi_index = 0;
			for (i = 0; i<IMAX; i++){
				for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
					jn1 = mpi_j1[mpi_jindex] + j;
					if (jn1<JMAX1 || jn1 >= JMAX - JMAX2){//
						mpi_fs[mpi_index] = whole_Qw[i][j];
						mpi_index = mpi_index + 1;
					}
					if (i < IMAX01 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
						mpi_fs[mpi_index] = whole_Qw[i][j];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX02 && jn1<JMAX - JMAX2 && jn1 >= JMAX1){// && IMAX1 >= 0
						mpi_fs[mpi_index] = whole_Qw[i][j];
						mpi_index = mpi_index + 1;
					}
				}
			}
			tag = 27;
			MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
	}
}

void air::Update_PML_Qp()
{
	int i, j, in, jn, jn1, i1, i2;
	int mpi_index, rank_num, rank_send, rank_recv, tag;
	for (i = 0; i < mpi_IMAX[mpi_iindex]; i++){
		for (j = 0; j < mpi_JMAX[mpi_jindex]; j++){
			in = mpi_i1[mpi_iindex] + i;
			jn = mpi_j1[mpi_jindex] + j;
			if (jn < JMAX1)Qp[i][j] += diff_t*p_nn[i][j];//
			if (jn >= JMAX - JMAX2) Qp[i][j] += diff_t*p_nn[i][j];
			if (in < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){//jn<JMAX - JMAX2 && jn >= JMAX1 &&
				Qp[i][j] += diff_t*p_nn[i][j];
			}
			if (in >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1){ // && jn < JMAX - JMAX2 && jn >= JMAX1
				Qp[i][j] += diff_t*p_nn[i][j];
			}
		}
	}
	if (mpi_iindex != 0){
		rank_num = mpi_IMAX[mpi_iindex] * mpi_JMAX[mpi_jindex];
		rank_recv = mpi_jindex;
		if (mpi_iindex != 0)i1 = 1;
		else i1 = 0;
		if (mpi_iindex != mpi_yarea - 1)i2 = 1;
		else i2 = 0;
		mpi_index = 0;
		for (i = i1; i<mpi_IMAX[mpi_iindex] - i2; i++){
			for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
				in = mpi_i1[mpi_iindex] + i;
				jn = mpi_j1[mpi_jindex] + j;
				if (jn<JMAX1 || jn >= JMAX - JMAX2){//
					mpi_fs[mpi_index] = Qp[i][j];
					mpi_index = mpi_index + 1;
				}
				if (in < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
					mpi_fs[mpi_index] = Qp[i][j];
					mpi_index = mpi_index + 1;
				}
				if (in >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1){// && IMAX1 >= 0
					mpi_fs[mpi_index] = Qp[i][j];
					mpi_index = mpi_index + 1;
				}
			}
		}
		tag = 28;
		MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);

		rank_send = rank_recv;
		rank_num = IMAX*mpi_JMAX[mpi_jindex];
		tag = 29;
		MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		mpi_index = 0;
		for (i = 0; i<IMAX; i++){
			for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
				jn = mpi_j1[mpi_jindex] + j;
				if (jn<JMAX1 || jn >= JMAX - JMAX2){//
					whole_Qp[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
				if (i < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
					whole_Qp[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
				if (i >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1){// && IMAX1 >= 0
					whole_Qp[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
			}
		}
	}
	else{
		for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
			for (i = 0; i<mpi_IMAX[mpi_iindex]; i++){
				in = mpi_i1[mpi_iindex] + i;
				jn = mpi_j1[mpi_jindex] + j;
				if (jn<JMAX1 || jn >= JMAX - JMAX2){//
					whole_Qp[in][j] = Qp[i][j];
				}
				if (in < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
					whole_Qp[in][j] = Qp[i][j];
				}
				if (in >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1){// && IMAX1 >= 0
					whole_Qp[in][j] = Qp[i][j];
				}
			}
		}
		for (in = 1; in<mpi_yarea; in++){
			if (mpi_porous != 0) rank_send = (mpi_iindex + in)*mpi_porous + mpi_jindex;
			else rank_send = (mpi_iindex + in)*mpi_zarea + mpi_jindex;

			jn = mpi_jindex;
			rank_num = mpi_IMAX[in] * mpi_JMAX[jn];
			if (in != 0)i1 = 1;
			else i1 = 0;
			if (in != mpi_yarea - 1)i2 = 1;
			else i2 = 0;
			mpi_index = 0;
			tag = 28;
			MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			for (i = mpi_i1[in] + i1; i<mpi_i2[in] + 1 - i2; i++){
				for (j = 0; j<mpi_JMAX[jn]; j++){
					jn1 = mpi_j1[mpi_jindex] + j;
					if (jn1<JMAX1 || jn1 >= JMAX - JMAX2){//
						whole_Qp[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
					if (i < IMAX01 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
						whole_Qp[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX02 && jn1<JMAX - JMAX2 && jn1 >= JMAX1){// && IMAX1 >= 0
						whole_Qp[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
				}
			}
		}
		for (in = 1; in<mpi_yarea; in++){
			jn = mpi_jindex;
			rank_num = IMAX*mpi_JMAX[jn];
			if (mpi_porous != 0) rank_recv = (mpi_iindex + in)*mpi_porous + mpi_jindex;
			else rank_recv = (mpi_iindex + in)*mpi_zarea + mpi_jindex;
			mpi_index = 0;
			for (i = 0; i<IMAX; i++){
				for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
					jn1 = mpi_j1[mpi_jindex] + j;
					if (jn1<JMAX1 || jn1 >= JMAX - JMAX2){//
						mpi_fs[mpi_index] = whole_Qp[i][j];
						mpi_index = mpi_index + 1;
					}
					if (i < IMAX01 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
						mpi_fs[mpi_index] = whole_Qp[i][j];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX02 && jn1<JMAX - JMAX2 && jn1 >= JMAX1){// && IMAX1 >= 0
						mpi_fs[mpi_index] = whole_Qp[i][j];
						mpi_index = mpi_index + 1;
					}
				}
			}
			tag = 29;
			MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
	}
}



void air::SetQMove(int N)
{
	int i, j, in, jn;
	int temp_I;
	if (N == 0){
		for (i = 0; i<mpi_IMAX[mpi_iindex]; i++){
			for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
				in = mpi_i1[mpi_iindex] + i;
				jn = mpi_j1[mpi_jindex] + j;
				if (jn<JMAX1 || jn >= JMAX - JMAX2){//
					Qv[i][j] = 0.0;
					Qw[i][j] = 0.0;
					Qp[i][j] = 0.0;
				}
				if (in < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 9
					Qv[i][j] = 0.0;
					Qw[i][j] = 0.0;
					Qp[i][j] = 0.0;
				}
				if (in >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1 ){//&& IMAX1 >= 0
					Qv[i][j] = 0.0;
					Qw[i][j] = 0.0;
					Qp[i][j] = 0.0;
				}
				temp_I = 0;
			}
		}
	}
	else{
		temp_I = move_frame.IMAX - move_frame.lead_DI;
	}

	if (temp_I != 0){
		for (i = 0; i<mpi_IMAX[mpi_iindex]; i++){
			for (j = 0; j<mpi_JMAX[mpi_jindex]; j++){
				in = mpi_i1[mpi_iindex] + i;
				jn = mpi_j1[mpi_jindex] + j;
				if (in<temp_I){
					if (jn<JMAX1 || jn >= JMAX - JMAX2){//
						Qv[i][j] = whole_Qv[in + move_frame.trail_DI - 1][j];
						Qw[i][j] = whole_Qw[in + move_frame.trail_DI - 1][j];
						Qp[i][j] = whole_Qp[in + move_frame.trail_DI - 1][j];
					}
					if (in >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 0
						Qv[i][j] = 0.0;
						Qw[i][j] = 0.0;
						Qp[i][j] = 0.0;
					}
					if (in < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1){//&& IMAX1 >= 9
						Qv[i][j] = 0.0;
						Qw[i][j] = 0.0;
						Qp[i][j] = 0.0;
					}
				}
				else{
					if (jn<JMAX1 || jn >= JMAX - JMAX2){//
						Qv[i][j] = 0.0;
						Qw[i][j] = 0.0;
						Qp[i][j] = 0.0;
					}
					if (in >= IMAX - IMAX02 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 9){// && IMAX1 >= 0
						Qv[i][j] = 0.0;
						Qw[i][j] = 0.0;
						Qp[i][j] = 0.0;
					}
					if (in < IMAX01 && jn<JMAX - JMAX2 && jn >= JMAX1){// && IMAX1 >= 9
						Qv[i][j] = 0.0;
						Qw[i][j] = 0.0;
						Qp[i][j] = 0.0;
					}
				}
			}
		}
	}
}