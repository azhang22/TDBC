//+++++++++++++++++++++++filename: calculation_2D.cpp +++++++++++++++++++++++++++++//

//-----------------------fathe class for air and porous class------------------//
/*******************************************************************************/
/* air and porous class will inherit from the class calculation_2D.
/*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "calculation_2D.h"

calculation_2D::calculation_2D(scheme_list scheme1,DifferenceStep Step,char *coordi,
					MovingFrame MF,const double gauss_width,AirStruct AirPara,vortex vortex1,
					int mpi_rank1,int mpi_size1,int mpi_yarea1,int mpi_zarea1,int mpi_porous1)
{
	int mpi_ny,mpi_nz;
	mpi_rank=mpi_rank1;
	mpi_size=mpi_size1;
	mpi_yarea=mpi_yarea1;
	mpi_zarea=mpi_zarea1;
	mpi_porous=mpi_porous1;
	scheme=scheme1;
	move_frame=MF;
	diff_t=Step.diff_t;
	diff_y=Step.diff_y;
	diff_z=Step.diff_z;
	gaussian_coef=4.0*log(2.0)/pow(gauss_width,2);
	//gaussian_coef=4.0*log(2.0)/pow(7.5e-3,2);
	//gaussian_coef=4.0*log(2.0)/pow(0.6,2);
	//gaussian_coef=4.0*log(2.0)/pow(5.0e-3,2);
	//gaussian_coef=40*7*7/pow(550*diff_y,2);
	//gaussian_coef=40*7*7/pow(400*diff_y,2);

	int i,j,k;
	//read data from files
	ifstream myfile(coordi,ios::in|ios::binary);
	if(!myfile){
		cout<<"cannot open file:"<<coordi<<" for read!!!"<<endl;
		return;
	}
	
	myfile >> whole_IMAX >> whole_JMAX >> whole_JMAX1 >> whole_JMAX2 >> whole_KMAX1 >> whole_KMAX2 >> whole_IJKMAX;
	//the following setting is for one moving frame
	int N_limit;
	if(MF.Judge==1){
		if ((whole_IJKMAX-move_frame.IMAX)%(move_frame.lead_DI-1)!=0){
			N_limit=int((whole_IJKMAX-move_frame.IMAX)/(move_frame.lead_DI-1)+1.0e-10);
			whole_IJKMAX=move_frame.IMAX+(N_limit+1)*(move_frame.lead_DI-1);
		}
		IMAX=move_frame.IMAX;
		JMAX=whole_JMAX;
		JMAX1=whole_JMAX1;
		JMAX2=whole_JMAX2;
		IMAX2=whole_KMAX2;
		IMAX1 = 0;
		IMAX02 = whole_KMAX2;
		IMAX01 = whole_KMAX1;

	}else{
		IMAX=whole_IMAX;JMAX=whole_JMAX;
		JMAX1=whole_JMAX1;JMAX2=whole_JMAX2;
		IMAX2=whole_KMAX2;
		IMAX1=0;
		IMAX02 = whole_KMAX2;
		IMAX01 = whole_KMAX1;
		whole_IJKMAX=IMAX;
	}
	// mpi:specify the domain decomposition
	mpi_ny=int(IMAX/mpi_yarea);
	mpi_nz=int(JMAX/mpi_zarea);
	mpi_iindex=mpi_rank/mpi_zarea;
	mpi_jindex=mpi_rank%mpi_zarea;
	if (mpi_porous!=0){
		mpi_iindex=mpi_rank/mpi_porous;
		mpi_jindex=mpi_rank%mpi_porous;
	}
	mpi_i1=new int[mpi_yarea];
	mpi_i2=new int[mpi_yarea];
	mpi_j1=new int[mpi_zarea];
	mpi_j2=new int[mpi_zarea];
	mpi_IMAX=new int[mpi_yarea];
	mpi_JMAX=new int[mpi_zarea];
	//overlap information
	n_order=vortex1.n_order;
	switch(n_order){
	case 1:
		nno=1;
		break;
	case 2:
		nno=2;
		break;
	case 3:
		nno=2;
		break;
	case 5:
		nno=3;
		break;
	default:
		nno=1;
	}

	for (i=0;i<mpi_yarea;i++){
		if (i!=0)mpi_i1[i]=i*mpi_ny+1-nno;
		else mpi_i1[i]=0;
		if (i!=mpi_yarea-1) mpi_i2[i]=(i+1)*mpi_ny+nno;
		else mpi_i2[i]=IMAX-1;
		mpi_IMAX[i]=mpi_i2[i]-mpi_i1[i]+1;
	}
	for (j=0;j<mpi_zarea;j++){
		if (j!=0)mpi_j1[j]=j*mpi_nz+1-nno;
		else mpi_j1[j]=0; 
		if (j!=mpi_zarea-1) mpi_j2[j]=(j+1)*mpi_nz+nno;
		else mpi_j2[j]=JMAX-1;
		mpi_JMAX[j]=mpi_j2[j]-mpi_j1[j]+1;
	}
	whole_Y=new double* [IMAX];whole_Z=new double* [IMAX];
	for(i=0;i<IMAX;i++){
		whole_Y[i]=new double [mpi_JMAX[mpi_jindex]];
		whole_Z[i]=new double [mpi_JMAX[mpi_jindex]];
	}
	leftbound_Y=new double[JMAX];leftbound_Z=new double[JMAX];
	bottombound_Y=new double[whole_IMAX];
	int jn;
	double temp_Y,temp_Z;
	for(j=0;j<JMAX;j++){
		for(i=0;i<whole_IMAX;i++){
			myfile>>temp_Y;
			if (j>=mpi_j1[mpi_jindex]&&j<=mpi_j2[mpi_jindex]){
				jn=j-mpi_j1[mpi_jindex];
				whole_Y[i][jn]=temp_Y;
			}
			if (i==0) leftbound_Y[j]=temp_Y;
			if (j==0) bottombound_Y[i]=temp_Y;
		}
	}
	for(j=0;j<JMAX;j++){
		for(i=0;i<whole_IMAX;i++){
			myfile>>temp_Z;
			if (j>=mpi_j1[mpi_jindex]&&j<=mpi_j2[mpi_jindex]){
				jn=j-mpi_j1[mpi_jindex];
				whole_Z[i][jn]=temp_Z;
			}
			if (i==0)leftbound_Z[j]=temp_Z;
		}
	}
	myfile.close();
	whole_v=new double* [IMAX];
	whole_w=new double* [IMAX];
	whole_p=new double* [IMAX];
	for(i=0;i<IMAX;i++){
		whole_v[i]=new double [mpi_JMAX[mpi_jindex]];
		whole_w[i]=new double [mpi_JMAX[mpi_jindex]];
		whole_p[i]=new double [mpi_JMAX[mpi_jindex]];
	}
	// establish the start of computational domain
	if (mpi_iindex!=0) ii1=nno;
	else ii1=1;
	if (mpi_iindex!=mpi_yarea-1) ii2=nno;
	else ii2=1;

	if (mpi_jindex!=0) jj1=nno;
	else jj1=1;
	if (mpi_jindex!=mpi_zarea-1)jj2=nno;
	else jj2=1;
	
	// mpi: The specify sender and receiver matrix for pressure;
	pss=new double[mpi_IMAX[mpi_iindex]*nno];psr=new double[mpi_IMAX[mpi_iindex]*nno];
	pns=new double[mpi_IMAX[mpi_iindex]*nno];pnr=new double[mpi_IMAX[mpi_iindex]*nno];
	pes=new double[mpi_JMAX[mpi_jindex]*nno];per=new double[mpi_JMAX[mpi_jindex]*nno];
	pws=new double[mpi_JMAX[mpi_jindex]*nno];pwr=new double[mpi_JMAX[mpi_jindex]*nno];
	// v,w velocity;
	vss=new double[mpi_IMAX[mpi_iindex]*nno];vsr=new double[mpi_IMAX[mpi_iindex]*nno];
	vns=new double[mpi_IMAX[mpi_iindex]*nno];vnr=new double[mpi_IMAX[mpi_iindex]*nno];
	ves=new double[mpi_JMAX[mpi_jindex]*nno];ver=new double[mpi_JMAX[mpi_jindex]*nno];
	vws=new double[mpi_JMAX[mpi_jindex]*nno];vwr=new double[mpi_JMAX[mpi_jindex]*nno];
	wss=new double[mpi_IMAX[mpi_iindex]*nno];wsr=new double[mpi_IMAX[mpi_iindex]*nno];
	wns=new double[mpi_IMAX[mpi_iindex]*nno];wnr=new double[mpi_IMAX[mpi_iindex]*nno];
	wes=new double[mpi_JMAX[mpi_jindex]*nno];wer=new double[mpi_JMAX[mpi_jindex]*nno];
	wws=new double[mpi_JMAX[mpi_jindex]*nno];wwr=new double[mpi_JMAX[mpi_jindex]*nno];
	//mpi_fps;mpi_fpr;mpi_fvs;mpi_fvr;mpi_fws;mpi_fwr
	mpi_fs=new double[IMAX*(mpi_JMAX[mpi_zarea-1]+nno)];
	mpi_fr=new double[IMAX*(mpi_JMAX[mpi_zarea-1]+nno)];
	mpi_fr1=new double[IMAX*(mpi_JMAX[mpi_zarea-1]+nno)];
	mpi_fr2=new double[IMAX*(mpi_JMAX[mpi_zarea-1]+nno)];
	mpi_fr3 = new double[IMAX*(mpi_JMAX[mpi_zarea - 1] + nno)];
	mpi_fr4 = new double[IMAX*(mpi_JMAX[mpi_zarea - 1] + nno)];

	//add the temperature distribution in the air
	double speed_air,temperature_air;
	Aav=new double[mpi_JMAX[mpi_jindex]];
//	double slope_temperature;
	double slope_sound,speed_air1,speed_air2;
	if (leftbound_Z[JMAX-1]<1.e-6){  // temperature at the ground
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			temperature_air=AirPara.temperature1;
			speed_air=331.3*sqrt(1+temperature_air/273.15);
			Aav[j]=speed_air*speed_air/AirPara.adiabatic_coef/AirPara.Pav;
		}
	}else{
		/*
		slope_temperature=(AirPara.temperature2-AirPara.temperature1)/fabs(leftbound_Z[JMAX-1]-leftbound_Z[0]);
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			jn=mpi_j1[mpi_jindex]+j;
			temperature_air=slope_temperature*(leftbound_Z[jn]-leftbound_Z[0])+AirPara.temperature1;
			speed_air=331.3*sqrt(1+temperature_air/273.15);
			Aav[j]=speed_air*speed_air/AirPara.adiabatic_coef/AirPara.Pav;
		}
		*/
		speed_air1=331.3*sqrt(1+AirPara.temperature1/273.15);
		speed_air2=331.3*sqrt(1+AirPara.temperature2/273.15);
		slope_sound=(speed_air2-speed_air1)/fabs(leftbound_Z[JMAX-1]-leftbound_Z[0]);
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			jn=mpi_j1[mpi_jindex]+j;
			speed_air=slope_sound*(leftbound_Z[jn]-leftbound_Z[0])+speed_air1;
			Aav[j]=speed_air*speed_air/AirPara.adiabatic_coef/AirPara.Pav;
		}
	}
	Y=new double* [mpi_IMAX[mpi_iindex]];Z=new double* [mpi_IMAX[mpi_iindex]];
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		Y[i]=new double [mpi_JMAX[mpi_jindex]];Z[i]=new double [mpi_JMAX[mpi_jindex]];
	}
	Y_over_y=new double* [mpi_IMAX[mpi_iindex]];Y_over_z=new double* [mpi_IMAX[mpi_iindex]];
	Z_over_y=new double* [mpi_IMAX[mpi_iindex]];Z_over_z=new double* [mpi_IMAX[mpi_iindex]];
	Jacobian=new double* [mpi_IMAX[mpi_iindex]];
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		Y_over_y[i]=new double [mpi_JMAX[mpi_jindex]];Y_over_z[i]=new double [mpi_JMAX[mpi_jindex]];
		Z_over_y[i]=new double [mpi_JMAX[mpi_jindex]];Z_over_z[i]=new double [mpi_JMAX[mpi_jindex]];
		Jacobian[i]=new double [mpi_JMAX[mpi_jindex]];
	}
	//allocate memory for velocity and pressure
	v_n=new double* [mpi_IMAX[mpi_iindex]];w_n=new double* [mpi_IMAX[mpi_iindex]];
	p_n=new double* [mpi_IMAX[mpi_iindex]];
	v_nn=new double* [mpi_IMAX[mpi_iindex]];w_nn=new double* [mpi_IMAX[mpi_iindex]];
	p_nn=new double* [mpi_IMAX[mpi_iindex]];
	psi_1 = new double* [mpi_IMAX[mpi_iindex]];psi_2 = new double* [mpi_IMAX[mpi_iindex]];
	psi_3 = new double* [mpi_IMAX[mpi_iindex]];psi_4 = new double* [mpi_IMAX[mpi_iindex]];
	psi_5 = new double* [mpi_IMAX[mpi_iindex]];psi_6 = new double* [mpi_IMAX[mpi_iindex]];
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		v_n[i]=new double [mpi_JMAX[mpi_jindex]];w_n[i]=new double [mpi_JMAX[mpi_jindex]];
		p_n[i]=new double [mpi_JMAX[mpi_jindex]];
		v_nn[i]=new double [mpi_JMAX[mpi_jindex]];w_nn[i]=new double [mpi_JMAX[mpi_jindex]];
		p_nn[i]=new double [mpi_JMAX[mpi_jindex]];
		psi_1[i]=new double [mpi_JMAX[mpi_jindex]];psi_2[i]=new double [mpi_JMAX[mpi_jindex]];
		psi_3[i]=new double [mpi_JMAX[mpi_jindex]];psi_4[i]=new double [mpi_JMAX[mpi_jindex]];
		psi_5[i]=new double [mpi_JMAX[mpi_jindex]];psi_6[i]=new double [mpi_JMAX[mpi_jindex]];
	}

	// add the timeseries data
	jjmax1=0;
	jjmax2=0;
	judge_left=1;
	// time_step=new double[100];
	// init_pressure=new double[100];
	if (judge_left==1){
		/*
		char file_interp[100]=" ";
		strcpy(file_interp,"intial_pressure_position.dat");
		ifstream init_file(file_interp,ios::in|ios::binary);
		init_file>>num_data;
		init_file.ignore(100,'\n');
		for (k=0;k<num_data;k++){
			init_file>>time_step[k]>>init_pressure[k];
			init_file.ignore(100,'\n');
		}
		*/
		if(whole_JMAX==3201){
			jjmax1=201;
			jjmax2=3001;
		}
		if(whole_JMAX==1501){
			jjmax1=1;
			jjmax2=1251;
		}
//		init_file.close();
	}
	//a1 = -0.55590;	a2 = 0.13751; a3 = 0.48279; a4 = 0.14898; a5 = 0.28696; a6 = -0.00037;
	//gamma1 = 0.57559; gamma2 = 0.01157; gamma3 = 0.63507; gamma4 = 0.11357; gamma5 = 0.41958; gamma6 = 2.63762;
	a1 = 0.27941;	a2 = -0.07950; a3 = 0.16538; a4 = 0.19467; a5 = 0.17342; a6 = 0.53160;
	gamma1 = 0.82597; gamma2 = 0.04010; gamma3 = 0.33455; gamma4 = 0.01627; gamma5 = 0.10557; gamma6 = 2.31909;
}

calculation_2D::~calculation_2D()
{
	int i;
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		delete[] p_n[i];delete[] w_n[i];delete[] v_n[i];
		delete[] p_nn[i];delete[] w_nn[i];delete[] v_nn[i];
		delete[] Y_over_y[i];delete[] Y_over_z[i];
		delete[] Z_over_y[i];delete[] Z_over_z[i];
		delete[] Y[i];delete[] Z[i];delete[] Jacobian[i];
	}
	delete[] v_n;delete[] w_n;delete[] p_n;
	delete[] v_nn;delete[] w_nn;delete[] p_nn;
	delete[] Y_over_y;delete[] Y_over_z;
	delete[] Z_over_y;delete[] Z_over_z;
	delete[] Y;delete[] Z;delete[] Jacobian;
	for(i=0;i<IMAX;i++){
		delete[] whole_v[i];delete[] whole_w[i];delete[] whole_p[i];
		delete[] whole_Y[i];delete[] whole_Z[i];
	}
	delete[] whole_v;delete[] whole_w;delete[] whole_p;
	delete[] whole_Y;delete[] whole_Z;
	delete[] Aav;
	// delete[] time_step; delete[] init_pressure;
	delete[] pss; delete[] psr; delete[] pns; delete[] pnr;
	delete[] pes; delete[] per; delete[] pws; delete[] pwr;
	delete[] vss; delete[] vsr; delete[] vns; delete[] vnr;
	delete[] ves; delete[] ver; delete[] vws; delete[] vwr;
	delete[] wss; delete[] wsr; delete[] wns; delete[] wnr;
	delete[] wes; delete[] wer; delete[] wws; delete[] wwr;
	delete[] mpi_fs;delete[] mpi_fr;delete[] mpi_fr1;delete[] mpi_fr2;delete[] mpi_fr3;delete[] mpi_fr4;
	delete[] mpi_i1;delete[] mpi_i2;delete[] mpi_j1;
	delete[] mpi_j2;delete[] mpi_IMAX;delete[] mpi_JMAX;
	delete[] leftbound_Y;delete[] leftbound_Z;
	delete[] bottombound_Y;
}
double calculation_2D::cal_InitialPressure(double y,double z)
{
	double r;
	r=sqrt(pow(y,2)+pow(z,2));
	return (1.0*exp(-gaussian_coef*r*r));//this is from Salomons paper(2002) p=exp(-40*r*r)
}

void calculation_2D::Set_InitialCond(Position source)
{

	int i,j;
	for(j=0;j<mpi_JMAX[mpi_jindex];j++){
		for(i=0;i<whole_IMAX;i++){
			whole_v[i][j]=0.0;
			whole_w[i][j]=0.0;
			if (judge_left==1) whole_p[i][j]=0.0;
			else {
				whole_p[i][j]=cal_InitialPressure(whole_Y[i][j]-source.y,
				whole_Z[i][j]-source.z);
			}			
		}
	}
	if (mpi_rank == 0) cout << "Source initiated " << endl;
}

void calculation_2D::UpdateSource(double f0, int N, Position source)
{
	int i, j, k;
	double t0, t, omega0;
	t0 = 1 / (f0 - 1e-20);
	t = N*diff_t;
	omega0 = f0 * 2 * 3.1415926;
	//if (mpi_rank == 0) cout << "Start to set source." << endl;
	for (int j = 0; j < mpi_JMAX[mpi_jindex]; j++){
		for (int i = 0; i < (whole_IMAX); i++){
			if (whole_Y[i][j] == source.y && whole_Z[i][j] == source.z)
				p_nn[i][j] = -1 * (t - t0)*omega0*omega0*exp(-1 * (t - t0)*(t - t0)*omega0*omega0 / 2);
		}
	}
	//if (mpi_rank == 0) cout << "Source updated." << endl;
}

void calculation_2D::get_position(int &ii,int &jj,Position receiver)
{
	int i,j,jn;
/*
	for (i=0;i<whole_IJKMAX;i++){
		for (j=0;j<JMAX;j++){
			if (((receiver.y-(leftbound_Y[j]+i*diff_y))<(diff_y/8.0))&&
				((receiver.z-leftbound_Z[j])<(diff_z/8.0))){
				ii=i;jj=j;
				return;
			}
		}
	}
*/
// for fixed domain we can use the following method
	ii=0;
	jj=0;
	for (j=0;j<JMAX;j++){
		if((receiver.z-leftbound_Z[j])<(diff_z/8.0)) {
			jj=j;
			break;
		}
	}
	for (i=0;i<IMAX;i++){
		if ((receiver.y-bottombound_Y[i])<(diff_y/8.0)){
			ii=i;
			return;
		}
	}
}
void calculation_2D::get_circle(double y_center, double z_center, double radius_cir,
	int num_rec, int *temp_i1, int *temp_j1, double *areay_a, double *areaz_b)
{
	// for fixed domain
	int i, j;
	double angle, yy, zz;
	for (i = 0; i<num_rec; i++){
		angle = i * 2 * 3.1415926 / double(num_rec - 1);
		yy = cos(angle)*radius_cir + y_center;
		zz = sin(angle)*radius_cir + z_center;
		temp_i1[i] = int(yy / diff_y);
		temp_j1[i] = int(zz / diff_z);
		areay_a[i] = yy / diff_y - int(yy / diff_y);
		areaz_b[i] = zz / diff_z - int(zz / diff_z);
		if (temp_i1[i] == 0) temp_i1[i] = 1;
	}
}

void calculation_2D::get_line(double y_start, double z_start, double y_end, double z_end,
	int num_rec, int *temp_i1, int *temp_j1, double *areay_a, double *areaz_b)
{
	// for fixed domain
	int i, j;
	double space, yy, zz, slope, intercept, co;
	for (i = 0; i < num_rec; i++){
		space = sqrt((y_start - y_end)*(y_start - y_end) + (z_start - z_end) *(z_start - z_end)) / double(num_rec - 1);
		if (y_start - y_end != 0 && z_start - z_end != 0){
			slope = (z_start - z_end) / (y_start - y_end);
			intercept = z_start - slope*y_start;
			co = (y_start - y_end) / sqrt((y_start - y_end)*(y_start - y_end) + (z_start - z_end) *(z_start - z_end));
			yy = y_start + i*space*co;
			zz = slope*yy + intercept;
		}
		if (y_start == y_end){
			yy = y_start;
			zz = z_start + i * space;
			if (mpi_rank == 0) cout << "zz = " << zz << endl;
		}
		if (z_start == z_end){
			yy = y_start + i * space;
			zz = z_start;
			if (mpi_rank == 0) cout << "yy = " << yy << endl;
		}

		temp_i1[i] = int(yy / diff_y);
		temp_j1[i] = int(zz / diff_z);
		if (mpi_rank == 0) cout << "temp_j1[i]" << temp_j1[i] << endl;
		areay_a[i] = yy / diff_y - int(yy / diff_y);
		areaz_b[i] = zz / diff_z - int(zz / diff_z);
		//if (temp_i1[i] == 0) temp_i1[i] = 1;
	}
}

void calculation_2D::UpdateInitialCond(int N)
{
	//transfer data under time n+1 to data under time n
	int i,j,in,jn,tag,rank_send,rank_recv,rank_num;
	int i1,i2,j1,j2,mpi_index;
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			v_n[i][j]=v_nn[i][j];
			w_n[i][j]=w_nn[i][j];
			p_n[i][j]=p_nn[i][j];
		}
		//if (mpi_rank == 0 && i == 876) cout << "After each step, w_n[i][0] = " << w_n[i][0]<< endl;
	}
	//update v,w,p for whole domain
	if (mpi_iindex!=0){
		rank_num=mpi_IMAX[mpi_iindex]*mpi_JMAX[mpi_jindex];
		rank_recv=mpi_jindex;
		if (mpi_iindex!=0)i1=nno;
		else i1=0;
		if (mpi_iindex!=mpi_yarea-1)i2=nno;
		else i2=0;
		if (N!=0){
			mpi_index=0;
			for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
				for(j=0;j<mpi_JMAX[mpi_jindex];j++){
					mpi_fs[mpi_index]=Y[i][j];
					mpi_index=mpi_index+1;
				}
			}
			tag=16;
			MPI_Send(mpi_fs,rank_num,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);

			mpi_index=0;
			for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
				for(j=0;j<mpi_JMAX[mpi_jindex];j++){
					mpi_fs[mpi_index]=Z[i][j];
					mpi_index=mpi_index+1;
				}
			}
			tag=17;
			MPI_Send(mpi_fs,rank_num,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
		}
		mpi_index=0;
		for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				mpi_fs[mpi_index]=v_nn[i][j];
				mpi_index=mpi_index+1;
			}
		}
		tag=18;
		MPI_Send(mpi_fs,rank_num,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);

		mpi_index=0;
		for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				mpi_fs[mpi_index]=w_nn[i][j];
				mpi_index=mpi_index+1;
			}
		}
		tag=19;
		MPI_Send(mpi_fs,rank_num,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);

		mpi_index=0;
		for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				mpi_fs[mpi_index]=p_nn[i][j];
				mpi_index=mpi_index+1;
			}
		}
		tag=20;
		MPI_Send(mpi_fs,rank_num,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);

		rank_send=rank_recv;
		rank_num=IMAX*mpi_JMAX[mpi_jindex];
		tag=21;
		MPI_Recv(mpi_fr,rank_num,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
		mpi_index=0;
		for (i=0;i<IMAX;i++){
			for (j=0;j<mpi_JMAX[mpi_jindex];j++){
				whole_v[i][j]=mpi_fr[mpi_index];
				mpi_index=mpi_index+1;
			}
		}

		tag=22;
		MPI_Recv(mpi_fr,rank_num,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
		mpi_index=0;
		for (i=0;i<IMAX;i++){
			for (j=0;j<mpi_JMAX[mpi_jindex];j++){
				whole_w[i][j]=mpi_fr[mpi_index];
				mpi_index=mpi_index+1;
			}
		}

		tag=23;
		MPI_Recv(mpi_fr,rank_num,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
		mpi_index=0;
		for (i=0;i<IMAX;i++){
			for (j=0;j<mpi_JMAX[mpi_jindex];j++){
				whole_p[i][j]=mpi_fr[mpi_index];
				mpi_index=mpi_index+1;
			}
		}
	}else{
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				whole_v[in][j]=v_nn[i][j];
				whole_w[in][j]=w_nn[i][j];
				whole_p[in][j]=p_nn[i][j];
				if (N!=0){
					whole_Y[in][j]=Y[i][j];
					whole_Z[in][j]=Z[i][j];
				}
			}
		}
		for (in=1;in<mpi_yarea;in++){
			if (mpi_porous!=0)rank_send=(mpi_iindex+in)*mpi_porous+mpi_jindex;
			else rank_send=(mpi_iindex+in)*mpi_zarea+mpi_jindex;
			jn=mpi_jindex;
			rank_num=mpi_IMAX[in]*mpi_JMAX[jn];
			if (in!=0)i1=nno;
			else i1=0;
			if (in!=mpi_yarea-1)i2=nno;
			else i2=0;
			if (N!=0){
				mpi_index=0;
				tag=16;
				MPI_Recv(mpi_fr,rank_num,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
				for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
					for(j=0;j<mpi_JMAX[jn];j++){
						whole_Y[i][j]=mpi_fr[mpi_index];
						mpi_index=mpi_index+1;
					}
				}
				mpi_index=0;
				tag=17;
				MPI_Recv(mpi_fr,rank_num,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
				for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
					for(j=0;j<mpi_JMAX[jn];j++){
						whole_Z[i][j]=mpi_fr[mpi_index];
						mpi_index=mpi_index+1;
					}
				}
			}
			mpi_index=0;
			tag=18;
			MPI_Recv(mpi_fr,rank_num,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
			for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
				for(j=0;j<mpi_JMAX[jn];j++){
					whole_v[i][j]=mpi_fr[mpi_index];
					mpi_index=mpi_index+1;
				}
			}
			mpi_index=0;
			tag=19;
			MPI_Recv(mpi_fr,rank_num,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
			for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
				for(j=0;j<mpi_JMAX[jn];j++){
					whole_w[i][j]=mpi_fr[mpi_index];
					mpi_index=mpi_index+1;
				}
			}
			mpi_index=0;
			tag=20;
			MPI_Recv(mpi_fr,rank_num,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
			for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
				for(j=0;j<mpi_JMAX[jn];j++){
					whole_p[i][j]=mpi_fr[mpi_index];
					mpi_index=mpi_index+1;
				}
			}
		}
		for (in=1;in<mpi_yarea;in++){
			jn=mpi_jindex;
			rank_num=IMAX*mpi_JMAX[jn];
			if (mpi_porous!=0)rank_recv=(mpi_iindex+in)*mpi_porous+mpi_jindex;
			else rank_recv=(mpi_iindex+in)*mpi_zarea+mpi_jindex;
			mpi_index=0;
			for (i=0;i<IMAX;i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					mpi_fs[mpi_index]=whole_v[i][j];
					mpi_index=mpi_index+1;
				}
			}
			tag=21;
			MPI_Send(mpi_fs,rank_num,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);

			mpi_index=0;
			for (i=0;i<IMAX;i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					mpi_fs[mpi_index]=whole_w[i][j];
					mpi_index=mpi_index+1;
				}
			}
			tag=22;
			MPI_Send(mpi_fs,rank_num,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);

			mpi_index=0;
			for (i=0;i<IMAX;i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					mpi_fs[mpi_index]=whole_p[i][j];
					mpi_index=mpi_index+1;
				}
			}
			tag=23;
			MPI_Send(mpi_fs,rank_num,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
		}
	}
}


double calculation_2D::get_pressure(int ii,int jj,int N)
{
	int i,il,in,jn,rank_send,rank_recv,tag;
	double temp_pvalue,temp_pvalue1;
	double p_rec;
	if ((ii>N*(move_frame.trail_DI-1))&&(ii<IMAX-1+N*(move_frame.trail_DI-1))){
		il=ii-N*(move_frame.trail_DI-1);
		if(il>=mpi_i1[mpi_iindex]+ii1 && il<=mpi_i2[mpi_iindex]-ii2&&
			jj>=mpi_j1[mpi_jindex]+jj1 && jj<=mpi_j2[mpi_jindex]-jj2){
			in=il-mpi_i1[mpi_iindex];
			jn=jj-mpi_j1[mpi_jindex];
			temp_pvalue = p_n[in][jn];
		}else temp_pvalue = 0.0;
	}else return 0.0;
	MPI_Allreduce(&temp_pvalue,&temp_pvalue1,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
	return temp_pvalue1;
}

void calculation_2D::SetMovingFrame(int N)
//N is total frame number starting from 0
{
	int mpi_ny,i;
	if(N==0){
		IMAX=whole_IMAX;
	}else{
		IMAX=move_frame.IMAX+N*(move_frame.lead_DI-move_frame.trail_DI);
	}
// Redefine the number of grids in y drection
	mpi_ny=int(IMAX/mpi_yarea);
	for (i=0;i<mpi_yarea;i++){
		mpi_i1[i]=i*mpi_ny;
		if (i!=mpi_yarea-1)mpi_i2[i]=(i+1)*mpi_ny+nno;
		else mpi_i2[i]=IMAX-1;
		mpi_IMAX[i]=mpi_i2[i]-mpi_i1[i]+1;
	};
	SetLeftRightMoving(N);
}

void calculation_2D::SetLeftRightMoving(int N)
{
	int i,j,in,jn;
	//transfer coordinate
	if (N==0){
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				Y[i][j]=whole_Y[in][j];
				Z[i][j]=whole_Z[in][j];
			}
		}
		Y_start=bottombound_Y[0];
	}else{
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				jn=mpi_j1[mpi_jindex]+j;
				Y[i][j]=bottombound_Y[in]+N*(move_frame.trail_DI-1)*diff_y;
				Z[i][j]=leftbound_Z[jn];
			}
		}
		Y_start=bottombound_Y[0]+N*(move_frame.trail_DI-1)*diff_y;
	}

	//calculate partial differential of coordinate transformation using central difference
		for(j=jj1;j<mpi_JMAX[mpi_jindex]-jj2;j++){
			for(i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){
				Y_over_y[i][j]=(Y[i+1][j]-Y[i-1][j])/diff_y/2;
				Z_over_y[i][j]=(Z[i+1][j]-Z[i-1][j])/diff_y/2;
				Y_over_z[i][j]=(Y[i][j+1]-Y[i][j-1])/diff_z/2;
				Z_over_z[i][j]=(Z[i][j+1]-Z[i][j-1])/diff_z/2;
				Jacobian[i][j]=Y_over_y[i][j]*Z_over_z[i][j]-Y_over_z[i][j]*Z_over_y[i][j];
		}
	}
	//update v,w,p for moving frame
	int temp_I;
	if (N==0){
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				v_n[i][j]=whole_v[in][j];
				w_n[i][j]=whole_w[in][j];
				p_n[i][j]=whole_p[in][j];
			}
		}
		temp_I=0;
   	 }else temp_I=move_frame.IMAX-move_frame.lead_DI;
	// temp_I is the overlapping part after frame is moved in y direction
	if (temp_I!=0){
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				if (in<temp_I){
					v_n[i][j]=whole_v[in+move_frame.trail_DI-1][j];
					w_n[i][j]=whole_w[in+move_frame.trail_DI-1][j];
					p_n[i][j]=whole_p[in+move_frame.trail_DI-1][j];
				}else{
					v_n[i][j]=0.0;
					w_n[i][j]=0.0;
					p_n[i][j]=0.0;
				}
			}
		}
		// specfy the left boundary condition which is used to avoid numerical error in the case wind shear.
		if (mpi_iindex==0){
			for (j=0;j<mpi_JMAX[mpi_jindex];j++)p_n[0][j]=p_n[1][j];
		}
	}
}
void calculation_2D::save_restart_cal(char *restartfile)
{
	int i,j;
	ofstream outfile11(restartfile,ios::app|ios::binary);
	outfile11.setf(ios::scientific,ios::floatfield);
	outfile11.setf(ios::showpos);
	outfile11.precision(6);
	outfile11.width(14);
	for(j=0;j<mpi_JMAX[mpi_jindex];j++){
		for(i=0;i<IMAX;i++){
			outfile11<<whole_Y[i][j]<<" ";
			outfile11<<whole_Z[i][j]<<" ";
			outfile11<<whole_v[i][j]<<" ";
			outfile11<<whole_w[i][j]<<" ";
			outfile11<<whole_p[i][j]<<" ";
		}
	}
	outfile11<<endl;
	outfile11.close();
}

void calculation_2D::input_restart_cal(char *restart_infile,int *point_pois,int N)
{
	int i,j,int_pois;
	N=N-1;
	if(N==(int)((whole_IJKMAX-move_frame.IMAX)/(move_frame.lead_DI-1)+1.0e-6)){
		IMAX=whole_IJKMAX-N*(move_frame.trail_DI-1);
	}else{
		IMAX=move_frame.IMAX+N*(move_frame.lead_DI-move_frame.trail_DI);
	}
	int_pois=*point_pois;
	ifstream infile(restart_infile,ios::in|ios::binary);
	infile.seekg(int_pois);
	for (j=0;j<mpi_JMAX[mpi_jindex];j++){
		for (i=0;i<IMAX;i++){
			infile>>whole_Y[i][j];
			infile>>whole_Z[i][j];
			infile>>whole_v[i][j];
			infile>>whole_w[i][j];
			infile>>whole_p[i][j];
		}
	}
	infile.ignore(100,'\n');
	int_pois=infile.tellg();
	*point_pois=int_pois;
	infile.close();
}

//-----------------------------update boundary condtions of pressure--------------------//
void calculation_2D::UpdateBC_pressure(boundary_location BC,int time_judge,int time_current)
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
						p_nn[0][j] = sin(2 * 3.1415926 * 500 * time_length);// 2*3.1415926*85 = 534.070742
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
			MPI_Recv(pwr,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
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
			MPI_Recv(per,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
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
			}
		}else{
			tag=3;
			if (mpi_porous!=0)rank_recv=mpi_iindex*mpi_porous+mpi_jindex-1;
			else rank_recv=mpi_iindex*mpi_zarea+mpi_jindex-1;
			MPI_Send(pss,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
			tag=4;
			rank_send=rank_recv;
			MPI_Recv(psr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
		}
	}
	break;
	case NorthBC: //characteristic boundary condition
	{
		if (mpi_jindex==mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-1;
			for(i=0;i<mpi_IMAX[mpi_iindex];i++) p_nn[i][jn]=p_nn[i][jn-1];//upper
		}else{
			tag=3;
			if (mpi_porous!=0)rank_send=mpi_iindex*mpi_porous+mpi_jindex+1;
			else rank_send=mpi_iindex*mpi_zarea+mpi_jindex+1;
			MPI_Recv(pnr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
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
/*
void calculation_2D::Cal_TDBC(double Z_inf, double tau)
{
	int i;
	double dt = diff_t / tau;
	double psi;
	if (mpi_jindex==0){
		for(i=0;i<mpi_IMAX[mpi_iindex];i++)
		{
			psi_1[i][3] = exp(-1 * gamma1 * dt) * psi_1[i][3] + a1 * w_nn[i][3] * dt;
			psi_2[i][3] = exp(-1 * gamma2 * dt) * psi_2[i][3] + a2 * w_nn[i][3] * dt;
			psi_3[i][3] = exp(-1 * gamma3 * dt) * psi_3[i][3] + a3 * w_nn[i][3] * dt;
			psi_4[i][3] = exp(-1 * gamma4 * dt) * psi_4[i][3] + a4 * w_nn[i][3] * dt;
			psi_5[i][3] = exp(-1 * gamma5 * dt) * psi_5[i][3] + a5 * w_nn[i][3] * dt;
			psi_6[i][3] = exp(-1 * gamma6 * dt) * psi_6[i][3] + a6 * w_nn[i][3] * dt;
			psi = psi_1[i][3] + psi_2[i][3] + psi_3[i][3] + psi_4[i][3] + psi_5[i][3] + psi_6[i][3];
			p_nn[i][3] = Z_inf * (w_nn[i][3] + psi);
			if (i == 876) cout << ", psi = " << psi << ", w_nn[i][3] = " << w_nn[i][3] << ", p_nn[i][3] = " << p_nn[i][3]; 
		}
	}	
}*/

void calculation_2D::Cal_TDBC(double Z_inf, double tau)
{
	int i;
	double dt = diff_t / tau;
	double sum;
	if (mpi_jindex==0){
		//cout << ", Z_inf = " << Z_inf << ", tau = " << tau ;
		for(i=0;i<mpi_IMAX[mpi_iindex];i++)
		{			
			sum = exp(-1 * gamma1 * dt) * psi_1[i][0];
			sum += exp(-1 * gamma2 * dt) * psi_2[i][0];
			sum += exp(-1 * gamma3 * dt) * psi_3[i][0];
			sum += exp(-1 * gamma4 * dt) * psi_4[i][0];
			sum += exp(-1 * gamma5 * dt) * psi_5[i][0];
			sum += exp(-1 * gamma6 * dt) * psi_6[i][0];
			
			w_nn[i][0] = p_nn[i][0] / Z_inf - sum;
			w_nn[i][0] = w_nn[i][0] / (1.0 + (a1 + a2 + a3 + a4 + a5 + a6) * dt);
			//if (i == 876) cout << ", sum = " << sum << ", w_nn[i][0] = " << w_nn[i][0] << ", p_nn[i][0] = " << p_nn[i][0]; 
			
			psi_1[i][0] = exp(-1 * gamma1 * dt) * psi_1[i][0] + a1 * w_nn[i][0] * dt;
			psi_2[i][0] = exp(-1 * gamma2 * dt) * psi_2[i][0] + a2 * w_nn[i][0] * dt;
			psi_3[i][0] = exp(-1 * gamma3 * dt) * psi_3[i][0] + a3 * w_nn[i][0] * dt;
			psi_4[i][0] = exp(-1 * gamma4 * dt) * psi_4[i][0] + a4 * w_nn[i][0] * dt;
			psi_5[i][0] = exp(-1 * gamma5 * dt) * psi_5[i][0] + a5 * w_nn[i][0] * dt;
			psi_6[i][0] = exp(-1 * gamma6 * dt) * psi_6[i][0] + a6 * w_nn[i][0] * dt;
		}
	}	
}

void calculation_2D::UpdateBC_pressure(boundary_location BC, double Z_inf, double tau)
{
	int i,j,k,in,jn,rank_send,rank_recv,tag;
	int no,nr1,nr2;
	double time_length,k_slope;
	// specify sender matrix in each block's boundary
	nr1=0;;
	for(i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){
		if (mpi_jindex!=0){
			for(no=0;no<nno;no++){ 
				pss[nr1]=p_nn[i][no+nno];
				nr1=nr1+1;
			}
		}
	}
	switch(BC){
	case SouthBC:
	{
		if (mpi_jindex==0){
			double dt = diff_t / tau;
			double psi, w;
			for(i=0;i<mpi_IMAX[mpi_iindex];i++)
			{
				w = sqrt(w_nn[i][0] * w_nn[i][0] + v_nn[i][0] * v_nn[i][0]);
				/*if (i == 876) cout << "before, psi_1 = " << psi_1[i][0] << ", psi_2 = " << psi_2[i][0] << ", psi_3 = " << psi_3[i][0] << ", psi_4 = " << psi_4[i][0] << ", psi_5 = " << psi_5[i][0] << ", psi_6 = " <<psi_6[i][0] << endl;*/
				psi_1[i][0] = exp(-1 * gamma1 * dt) * psi_1[i][0] + a1 * w * dt;
				psi_2[i][0] = exp(-1 * gamma2 * dt) * psi_2[i][0] + a2 * w * dt;
				psi_3[i][0] = exp(-1 * gamma3 * dt) * psi_3[i][0] + a3 * w * dt;
				psi_4[i][0] = exp(-1 * gamma4 * dt) * psi_4[i][0] + a4 * w * dt;
				psi_5[i][0] = exp(-1 * gamma5 * dt) * psi_5[i][0] + a5 * w * dt;
				psi_6[i][0] = exp(-1 * gamma6 * dt) * psi_6[i][0] + a6 * w * dt;
				psi = psi_1[i][0] + psi_2[i][0] + psi_3[i][0] + psi_4[i][0] + psi_5[i][0] + psi_6[i][0];
				/*if (i == 876) cout << "after, psi_1 = " << psi_1[i][0] << ", psi_2 = " << psi_2[i][0] << ", psi_3 = " << psi_3[i][0] << ", psi_4 = " << psi_4[i][0] << ", psi_5 = " << psi_5[i][0] << ", psi_6 = " <<psi_6[i][0] << endl;*/
				p_nn[i][0] = Z_inf * (w + psi);
				if (i == 876) cout << ", psi = " << psi <<", w = " << w << ", p_nn[i][0] = " << p_nn[i][0];
			}
		}else{
			tag=3;
			if (mpi_porous!=0)rank_recv=mpi_iindex*mpi_porous+mpi_jindex-1;
			else rank_recv=mpi_iindex*mpi_zarea+mpi_jindex-1;
			MPI_Send(pss,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);
			tag=4;
			rank_send=rank_recv;
			MPI_Recv(psr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
		}
	}
	break;
	default:
		cout<<"it is wrong to entrance boundary";
	}

	nr1=0;
	for (i=ii1;i<mpi_IMAX[mpi_iindex]-ii2;i++){
		if(mpi_jindex!=0){
			for(no=0;no<nno;no++){
				p_nn[i][no]=psr[nr1];
				nr1=nr1+1;
			}
		}
	}
}


//--------------------------update boundary conditions of velocity(V and W)------------------//
void calculation_2D::UpdateBC_velocity(boundary_location BC)
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
			MPI_Recv(vwr,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
			tag=8;
			MPI_Recv(wwr,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
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
			MPI_Recv(ver,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
			tag=6;
			MPI_Recv(wer,(mpi_JMAX[mpi_jindex]-jj1-jj2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
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
			MPI_Recv(vsr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
			tag=12;
			MPI_Recv(wsr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
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
			MPI_Recv(vnr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
			tag=10;
			MPI_Recv(wnr,(mpi_IMAX[mpi_iindex]-ii1-ii2)*nno,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD, &status);
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

void calculation_2D::UpdateBC_velocity(boundary_location BC,calculation_2D& porous)
{
	//update absorbing boundary
	int AirDy_CouplingDy,AirDz_CouplingDz;
	AirDy_CouplingDy=(int)(diff_y/porous.diff_y);
	AirDz_CouplingDz=(int)(diff_z/porous.diff_z);
	int i,j;
	switch(BC){
	case SouthBC:
		if(AirDz_CouplingDz<1.5){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				v_nn[i][0]=porous.v_nn[i][porous.mpi_JMAX[mpi_jindex]-2];
				porous.v_nn[i][porous.mpi_JMAX[mpi_jindex]-1]=v_nn[i][1];
				w_nn[i][0]=porous.w_nn[i][porous.mpi_JMAX[mpi_jindex]-2];
				porous.w_nn[i][porous.mpi_JMAX[mpi_jindex]-1]=w_nn[i][1];
			}
		}else{//transient interface
			for(i=0;i<(mpi_IMAX[mpi_iindex]-1);i++){
				for(int k=0;k<(AirDy_CouplingDy);k++){
					double temp_v;
					temp_v=(v_nn[i+1][1]-v_nn[i][1])/(AirDy_CouplingDy)*k+v_nn[i][1];
					porous.v_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-1]=
						(temp_v-porous.v_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2])/
						(AirDz_CouplingDz+1)+porous.v_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2];
					double temp_w;
					temp_w=(w_nn[i+1][1]-w_nn[i][1])/(AirDy_CouplingDy)*k+w_nn[i][1];
					porous.w_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-1]=
						(temp_w-porous.w_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2])/
						(AirDz_CouplingDz+1)+porous.w_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2];
				}
			}
			for(i=ii1;i<(mpi_IMAX[mpi_iindex]-ii2);i++){
				v_nn[i][0]=porous.v_nn[i*AirDy_CouplingDy][porous.mpi_JMAX[mpi_jindex]-1];
				w_nn[i][0]=porous.w_nn[i*AirDy_CouplingDy][porous.mpi_JMAX[mpi_jindex]-1];
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void calculation_2D::UpdateBC_pressure(boundary_location BC,calculation_2D& porous)
{
	//update absorbing boundary
	int AirDy_CouplingDy,AirDz_CouplingDz;
	AirDy_CouplingDy=(int)(diff_y/porous.diff_y);
	AirDz_CouplingDz=(int)(diff_z/porous.diff_z);
	int i,j;
	switch(BC){
	case SouthBC:
	{
		if(AirDz_CouplingDz<1.5){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				p_nn[i][0]=porous.p_nn[i][porous.mpi_JMAX[mpi_jindex]-2];
				porous.p_nn[i][porous.mpi_JMAX[mpi_jindex]-1]=p_nn[i][1];
			}
		}else{//transient interface
			for(i=0;i<mpi_IMAX[mpi_iindex]-1;i++){
				for(int k=0;k<AirDy_CouplingDy;k++){
					double temp_p;
					temp_p=(p_nn[i+1][1]-p_nn[i][1])/(AirDy_CouplingDy)*k+p_nn[i][1];
					porous.p_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-1]=
						(temp_p-porous.p_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2])/
						(AirDz_CouplingDz+1)+porous.p_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2];
				}
			}
			for(i=ii1;i<(mpi_IMAX[mpi_iindex]-ii2);i++){
				p_nn[i][0]=porous.p_nn[i*AirDy_CouplingDy][porous.mpi_JMAX[mpi_jindex]-1];
			}
		}
	}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}


//cal_fvw() and cal_fp() to get v_nn=fvw*diff_t;p_nn=fp*diff_t//
void calculation_2D::Cal_velocity(void)
{
	int i, j;
	switch (scheme){
	case FB_v_p://forward - backward scheme (1st)
	{
		cal_fvw(v_n, w_n, p_n);
		for (j = jj1; j < mpi_JMAX[mpi_jindex] - jj2; j++){
			for (i = ii1; i < mpi_IMAX[mpi_iindex] - ii2; i++){
				v_nn[i][j] = v_n[i][j] + v_nn[i][j];
				w_nn[i][j] = w_n[i][j] + w_nn[i][j];
				if (fabs(v_nn[i][j]) < 1.0e-200) v_nn[i][j] = 0.0;
				if (fabs(w_nn[i][j]) < 1.0e-200) w_nn[i][j] = 0.0; 
			}
		}
	}
	break;
	case FB_p_v://forward - backward scheme (1st)

		cal_fvw(v_n, w_n, p_nn);
		for (j = jj1; j < mpi_JMAX[mpi_jindex] - jj2; j++){
			for (i = ii1; i < mpi_IMAX[mpi_iindex] - ii2; i++){
				v_nn[i][j] = v_n[i][j] + v_nn[i][j];
				w_nn[i][j] = w_n[i][j] + w_nn[i][j];
				if (fabs(v_nn[i][j]) < 1e-200) v_nn[i][j] = 0.0;
				if (fabs(w_nn[i][j]) < 1e-200) w_nn[i][j] = 0.0;
			}
		}
		break;
	default:
		printf("wrong with scheme name");
	}//end switch
}

void calculation_2D::Cal_pressure(void)
{
	int i, j;
	//double temp_p;
	switch (scheme){
	case FB_v_p://forward - backward scheme (1st)
	{
		cal_fp(v_nn, w_nn, p_n);
		for (j = jj1; j < mpi_JMAX[mpi_jindex] - jj2; j++){
			for (i = ii1; i < mpi_IMAX[mpi_iindex] - ii2; i++){
				p_nn[i][j] = p_n[i][j] + p_nn[i][j];
				if (fabs(p_nn[i][j]) < 1e-200) p_nn[i][j] = 0.0;
			}
		}
	}
	break;
	case FB_p_v://forward - backward scheme (1st)
	{
		cal_fp(v_n, w_n, p_n);
		for (j = jj1; j < mpi_JMAX[mpi_jindex] - jj2; j++){
			for (i = ii1; i < mpi_IMAX[mpi_iindex] - ii2; i++){
				p_nn[i][j] = p_n[i][j] + p_nn[i][j];
				if (fabs(p_nn[i][j]) < 1e-200) p_nn[i][j] = 0.0;
			}
		}
	}
	break;
	default:
		printf("wrong with scheme name");
	}//end switch
}
void calculation_2D::mpi_send_data(void)
{
	int i,j,i1,i2,j1,j2,in,jn;
	int rank_num,mpi_index,tag,rank_send,rank_recv;
	int save_ii;
	save_ii=7;
	mpi_index=0;
	rank_recv=0;
	rank_num=IMAX*mpi_JMAX[mpi_jindex];
	if (mpi_jindex!=0)j1=nno;
	else j1=0;
	if (mpi_jindex!=mpi_zarea-1)j2=nno;
	else j2=0;
	for (j = j1; j < mpi_JMAX[mpi_jindex] - j2; j++){
		for (i = 0; i < IMAX; i++){
			if (i%save_ii == 0 || i == IMAX - 1){
				mpi_fs[mpi_index] = whole_Y[i][j];
				mpi_index = mpi_index + 1;
			}
		}
	}
	tag=13;
	MPI_Send(mpi_fs,rank_num,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);

	mpi_index=0;
	for (j = j1; j < mpi_JMAX[mpi_jindex] - j2; j++){
		for (i = 0; i < IMAX; i++){
			if (i%save_ii == 0 || i == IMAX - 1){
				mpi_fs[mpi_index] = whole_Z[i][j];
				mpi_index = mpi_index + 1;
			}
		}
	}
	tag=14;
	MPI_Send(mpi_fs,rank_num,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);

	mpi_index=0;
	for (j = j1; j < mpi_JMAX[mpi_jindex] - j2; j++){
		for (i = 0; i < IMAX; i++){
			if (i%save_ii == 0 || i == IMAX - 1){
				mpi_fs[mpi_index] = whole_p[i][j];
				mpi_index = mpi_index + 1;
			}
		}
	}
	tag=15;
	MPI_Send(mpi_fs,rank_num,MPI_DOUBLE,rank_recv,tag, MPI_COMM_WORLD);

	mpi_index = 0;
	for (j = j1; j < mpi_JMAX[mpi_jindex] - j2; j++){
		for (i = 0; i < IMAX; i++){
			if (i%save_ii == 0 || i == IMAX - 1){
				mpi_fs[mpi_index] = whole_v[i][j];
				mpi_index = mpi_index + 1;
			}
		}
	}
	tag = 151;
	MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);

	mpi_index = 0;
	for (j = j1; j < mpi_JMAX[mpi_jindex] - j2; j++){
		for (i = 0; i < IMAX; i++){
			if (i%save_ii == 0 || i == IMAX - 1){
				mpi_fs[mpi_index] = whole_w[i][j];
				mpi_index = mpi_index + 1;
			}
		}
	}
	tag = 152;
	MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);/**/
}
void getOstream(ostream& output_stream, calculation_2D& my_object)
{
	int IMAX1,IMAX,JMAX,i,j,jn,j1,j2;
	int mpi_index,rank_num,rank_send,tag;
	double save_idd;int save_ii;
	save_idd=7.0;
	save_ii=7;
	IMAX=my_object.IMAX;JMAX=my_object.JMAX;
	if ((IMAX-1)%save_ii==0)IMAX1=int((IMAX-1+1.0e-6)/save_idd)+1;
	else IMAX1=int((IMAX-1+1.0e-6)/save_idd)+2;
	output_stream<<"I="<<IMAX1<<",J="<<JMAX<<",F=POINT"<<endl;
	if (my_object.mpi_zarea-1!=0)j2=my_object.nno;//nno=j2=2
	else j2=0;
	for(j=0;j<my_object.mpi_JMAX[0]-j2;j++){
		for(i=0;i<IMAX;i++){
			if (i%save_ii==0||i==IMAX-1){
				output_stream << (my_object.whole_Y[i][j]) << " " << (my_object.whole_Z[i][j]) << " "
					<< my_object.whole_p[i][j] << " " << my_object.whole_v[i][j]  //endl;
					<< " " << my_object.whole_w[i][j] << endl;
			}
		}
	}
	for (rank_send = 1; rank_send<my_object.mpi_zarea; rank_send++){//my_object.mpi_zarea=10
		//cout << "rank_send = " << rank_send << ", jn = " << my_object.mpi_JMAX[rank_send] << endl;
		jn=my_object.mpi_JMAX[rank_send];
		if (rank_send != 0)j1 = my_object.nno;//j1 = 2
		else j1=0;
		if (rank_send != my_object.mpi_zarea - 1)j2 = my_object.nno;
		else j2=0;
		mpi_index=0;
		rank_num=IMAX*jn;
		tag = 13;
		MPI_Recv(my_object.mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &my_object.status);
		tag = 14;
		MPI_Recv(my_object.mpi_fr1, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &my_object.status);
		tag = 15;
		MPI_Recv(my_object.mpi_fr2, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &my_object.status);
		tag = 151;
		MPI_Recv(my_object.mpi_fr3, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &my_object.status);
		tag = 152;
		MPI_Recv(my_object.mpi_fr4, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &my_object.status);
		for(j=j1;j<jn-j2;j++){
			for(i=0;i<IMAX;i++){
				if (i%save_ii==0||i==IMAX-1){
					output_stream << (my_object.mpi_fr[mpi_index]) << " " << (my_object.mpi_fr1[mpi_index]) << " "
						<< my_object.mpi_fr2[mpi_index] << " " << my_object.mpi_fr3[mpi_index]//<< endl;
						 << " " << my_object.mpi_fr4[mpi_index] << endl;
					mpi_index=mpi_index+1;
				}
			}
		}
		//cout << ""
	}
}

ostream& operator<<(ostream& output_stream, calculation_2D& my_object)
{
	getOstream(output_stream,my_object);
	return output_stream;
}