//++++++++++++++++++++++++++++filename: airpore.cpp ++++++++++++++++++++++++++//

//-------------------------------airpore class---------------------------//
/*******************************************************************************/
/* the class construct computational environment to calculate
/* the coupling of air and porous
/* Note: AirStep.diff_y and AirStep.diff_z which is read through input.txt is the space step
/* of the domain above transient domain. Space step of whole air domain is calculated again
/* inside the fuction void airpore::get_coordi().
/* move_frame.IMAX and move_frame.diff_I is the grid points of moving frame and moving
/* grid points each time
/*******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "airpore.h"
#include "air.h"
#include "porous.h"
#include "time.h"
#include <direct.h>

//****************************airpore code*******************************//
airpore::airpore(char *input_file)
{
	ifstream infile(input_file);
	infile>>gauss_width;infile.ignore(100,'\n');
	infile>>move_frame.IMAX>>move_frame.lead_DI>>move_frame.trail_DI>>move_frame.Judge;infile.ignore(100,'\n');
	int scheme_type;
	infile>>scheme_type;infile.ignore(100,'\n');
	switch(scheme_type)
	{
	case 1:
		scheme=FB_v_p;break;
	case 2:
		scheme=FB_p_v;break;
	case 3:
		scheme=FB_vp;break;
	case 4:
		scheme=LeapTrap;break;
	case 5:
		scheme=RK4;break;
	case 6:
		scheme=RK2;break;
	case 7:
		scheme=ABM;break;
	default:
		scheme=FB_p_v;
	}

	int boundary_choice[4];
	int i;
	for(i=0;i<4;i++) infile>>boundary_choice[i];
	infile.ignore(100,'\n');
	infile.ignore(100,'\n');
	for(i=0;i<4;i++){
		switch(boundary_choice[i]){
		case 1://rigid
			{
				switch(i){
		case 0:
			air_boundary.air_north=rigid;break;
		case 1:
			air_boundary.air_south=rigid;break;
		case 2:
			air_boundary.air_west=rigid;break;
		case 3:
			air_boundary.air_east=rigid;break;
		default:
			cout<<"please entrance your choice again!";
				}
			}
			break;
		case 2://absorbing layer
			{
				switch(i){
		case 0:
			air_boundary.air_north=absorbing;break;
		case 1:
			air_boundary.air_south=absorbing;break;
		case 2:
			air_boundary.air_west=absorbing;break;
		case 3:
			air_boundary.air_east=absorbing;break;
		default:
			cout<<"please entrance your choice again!";
				}
			}
			break;
		case 3://porous media
			{
				switch(i){
		case 0:
			air_boundary.air_north=porous_media;break;
		case 1:
			air_boundary.air_south=porous_media;break;
		case 2:
			air_boundary.air_west=porous_media;break;
		case 3:
			air_boundary.air_east=porous_media;break;
		default:
			cout<<"please entrance your choice again!";
				}
			}
			break;
		case 4://porous media
			{
				switch(i){
		case 0:
			air_boundary.air_north=TDBC;break;
		case 1:
			air_boundary.air_south=TDBC;break;
		case 2:
			air_boundary.air_west=TDBC;break;
		case 3:
			air_boundary.air_east=TDBC;break;
		default:
			cout<<"please entrance your choice again!";
				}
			}
			break;
		default:
			cout<<"please entrance your choice again!";
		}
	}
	if(air_boundary.air_west==rigid)
	{
		if(air_boundary.air_east==rigid)
		{
			if(air_boundary.air_north==rigid)
			{
				CaseNo=1;
			}
			else
			{//absorbing
				CaseNo=2;
			}
		}
		else
		{
			if(air_boundary.air_south==absorbing)
			{
				CaseNo=6;
			}
			else
			{//absorbing
				CaseNo=5;
			}
		}
	}
	else
	{//absorbing
		if(air_boundary.air_south==absorbing)
		{
			CaseNo=4;
		}
		else
		{
			CaseNo=3;
		}
	}
	infile>>time_domain>>out_difft>>FFT_m>>frequency;infile.ignore(100,'\n');
	FFT_N=1;
	for (i=0;i<FFT_m;i++) FFT_N *= 2;
	//air properties
	infile>>velocity_coef;infile>>velocity_method;infile.ignore(100,'\n');
	//infile>>AirPara.adiabatic_coef>>AirPara.Pav>>AirPara.sound_speed;infile.ignore(100,'\n');
	infile>>AirPara.adiabatic_coef>>AirPara.Pav>>AirPara.temperature1>>AirPara.temperature2;
	infile.ignore(100,'\n');infile.ignore(100,'\n');

	speed_moveframe=331.3*sqrt(1+AirPara.temperature1/273.15);
	//AirPara.Aav=pow(AirPara.sound_speed,2)/AirPara.adiabatic_coef/AirPara.Pav;
	//absorbing layer structure
	infile>>AbsorbPara.Cs>>AbsorbPara.porosity>>AbsorbPara.resistivity;infile.ignore(100,'\n');
//	AbsorbPara.eff_density=AbsorbPara.Cs/AbsorbPara.porosity;  ///AirPara.Aav;
	AbsorbPara.Kp=1.0/AirPara.adiabatic_coef/AirPara.Pav;
	SideAbsorbPara=AbsorbPara;SideAbsorbPara.resistivity=AbsorbPara.resistivity/4.0;
	//porous media structure
	infile>>PorePara.q>>PorePara.porosity>>PorePara.resistivity;infile.ignore(100,'\n');
	//PorePara.eff_density=PorePara.Cs/PorePara.porosity/AirPara.Aav;
	PorePara.Cs = PorePara.q*PorePara.q;
	PorePara.eff_density = PorePara.Cs/PorePara.porosity;
	PorePara.Kp = 1.0/AirPara.adiabatic_coef/AirPara.Pav;
	double rho = AirPara.adiabatic_coef * AirPara.Pav / 340.29 / 340.29;
	PorePara.Z_inf = rho * 340.29 * PorePara.q / PorePara.porosity;
	PorePara.tau = AirPara.adiabatic_coef * rho * PorePara.q * PorePara.q / PorePara.resistivity / PorePara.porosity;
	
	infile>>PoreParb.q>>PoreParb.porosity>>PoreParb.resistivity;infile.ignore(100,'\n');
	//PorePara.eff_density=PorePara.Cs/PorePara.porosity/AirPara.Aav;
	PoreParb.Cs = PoreParb.q*PoreParb.q;
	PoreParb.eff_density = PoreParb.Cs/PoreParb.porosity;
	PoreParb.Kp = 1.0/AirPara.adiabatic_coef/AirPara.Pav;
	PoreParb.Z_inf = rho * 340.29 * PoreParb.q / PoreParb.porosity;
	PoreParb.tau = AirPara.adiabatic_coef * rho * PoreParb.q * PoreParb.q / PoreParb.resistivity / PoreParb.porosity;

	//source and receiver
	infile>>source.y>>source.z>>receiver.y>>receiver.z;infile.ignore(100,'\n');

	//air space dimension
	infile>>AirStep.diff_t;infile.ignore(100,'\n');
	infile>>AirStep.diff_y>>AirStep.diff_z;infile.ignore(100,'\n');
	infile>>AirFrame.left>>AirFrame.right>>AirFrame.upper>>AirFrame.lower;infile.ignore(100,'\n');

	//absorbing boundary setting
	infile>>NorthStep.diff_y>>NorthStep.diff_z;infile.ignore(100,'\n');
	infile>>NorthFrame.left>>NorthFrame.right>>NorthFrame.upper>>NorthFrame.lower;
	infile.ignore(100,'\n');
	infile>>SouthFrame.left>>SouthFrame.right>>SouthFrame.upper>>SouthFrame.lower;
	infile.ignore(100,'\n');

	//west frame absorbing boundary
	infile>>WestStep.diff_y>>WestStep.diff_z;infile.ignore(100,'\n');
	infile>>WestFrame.left>>WestFrame.right>>WestFrame.upper>>WestFrame.lower;
	infile.ignore(100,'\n');
	infile>>WestSouthFrame.left>>WestSouthFrame.right>>WestSouthFrame.upper>>WestSouthFrame.lower;
	infile.ignore(100,'\n');
	infile>>WestNorthFrame.left>>WestNorthFrame.right>>WestNorthFrame.upper>>WestNorthFrame.lower;
	infile.ignore(100,'\n');

	//east frame absorbing boundary
	infile>>EastStep.diff_y>>EastStep.diff_z;infile.ignore(100,'\n');
	infile>>EastFrame.left>>EastFrame.right>>EastFrame.upper>>EastFrame.lower;
	infile.ignore(100,'\n');
	infile>>EastSouthFrame.left>>EastSouthFrame.right>>EastSouthFrame.upper>>EastSouthFrame.lower;
	infile.ignore(100,'\n');
	infile>>EastNorthFrame.left>>EastNorthFrame.right>>EastNorthFrame.upper>>EastNorthFrame.lower;
	infile.ignore(100,'\n');

	SouthStep.diff_t=WestStep.diff_t=NorthStep.diff_t=
	WestSouthStep.diff_t=WestNorthStep.diff_t=EastStep.diff_t=
	EastSouthStep.diff_t=EastNorthStep.diff_t=AirStep.diff_t;
	SouthStep=NorthStep;
	WestNorthStep.diff_y=WestStep.diff_y;WestNorthStep.diff_z=NorthStep.diff_y;
	EastNorthStep.diff_y=EastStep.diff_y;EastNorthStep.diff_z=NorthStep.diff_y;
	EastSouthStep.diff_y=EastStep.diff_y;EastSouthStep.diff_z=SouthStep.diff_z;

	// MPI: Domain blocks yarea*zarea
	infile>>mpi_yarea>>mpi_zarea;infile.ignore(100,'\n');

	// output.data
    infile.ignore(100,'\n');infile.ignore(100,'\n');
    infile>>restart>>restart_out>>restart_infile;

	// curve setting,adds the porous mediea hill,
	/* curve_judge: judge the curve or not; curve_coefa,curve_coefb, curve_coefc
	and curve_coefd: The coefficient of equations:y=a*x^3+b*x^2+c*x+d;
	*/
	infile.ignore(100,'\n');infile.ignore(100,'\n');
	infile>>hill1.curve_judge>>hill1.curve_coefa>>hill1.curve_coefb>>
			  hill1.curve_coefc>>hill1.curve_coefd>>hill1.curve_coefe>>
			  hill1.curve_coeff>>hill1.curve_coefg>>hill1.curve_coefh>>
			  hill1.curve_coefi>>hill1.curve_coefj>>hill1.curve_coefk>>
			  hill1.curve_coefl>>hill1.curve_points;
	infile.ignore(100,'\n');
	// The domain we need to calculate for immersed boundary
	infile>>hill1.curve_ystar>>hill1.curve_ystop>>hill1.curve_zstar>>hill1.curve_zstop;
	infile.ignore(100,'\n');
	infile>>vortex1.n_order>>vortex1.ystar>>vortex1.ystop>>vortex1.zstar>>vortex1.zstop;
	infile.ignore(100,'\n');
	infile>>rec_circ1.judge>>rec_circ1.num>>rec_circ1.radius>>
												rec_circ1.center_y>>rec_circ1.center_z>>rec_circ1.lineend_y>>rec_circ1.lineend_z;
	infile.ignore(100,'\n');
	infile>> time_source >> f_0;
	
	//porous media setting
	infile.ignore(100,'\n');infile.ignore(100,'\n');
	if(air_boundary.air_south==porous_media)
	{
		infile>>SouthFrame.left>>SouthFrame.right>>SouthFrame.upper>>SouthFrame.lower;
		infile.ignore(100,'\n');
		infile>>SouthStep.diff_y>>SouthStep.diff_z;
		infile.ignore(100,'\n');
		EastSouthStep.diff_z=SouthStep.diff_z;
	}
	else
	{
		infile.ignore(100,'\n');infile.ignore(100,'\n');
		infile.ignore(100,'\n');infile.ignore(100,'\n');
	}
	infile.close();
	//test time step for moving frame according to cp*dx/dt=m
	//AirStep.diff_t=dx*m in input.txt
	/*
	if((AirStep.diff_t/AirStep.diff_z)>0.05){
		AirStep.diff_t=AirStep.diff_t/AirPara.sound_speed;
		SouthStep.diff_t=WestStep.diff_t=NorthStep.diff_t=
			WestSouthStep.diff_t=WestNorthStep.diff_t=EastStep.diff_t=AirStep.diff_t;
	}
	*/
}


airpore::~airpore()
{
}


void airpore::Cal_velocity()
{
	AirMedia->Cal_velocity();
	switch(CaseNo){
	case 1:
		{
			AirMedia->UpdateBC_velocity(NorthBC);
			AirMedia->UpdateBC_velocity(EastBC);
			AirMedia->UpdateBC_velocity(WestBC);
			if(air_boundary.air_south==rigid || air_boundary.air_south == TDBC){
				AirMedia->UpdateBC_velocity(SouthBC);
				if(air_boundary.air_south == TDBC){				
					double Z_inf = PorePara.Z_inf;
					double tau = PorePara.tau;
					//double tau = 1.1e-3;
					AirMedia->Cal_TDBC(Z_inf,tau);
				}
			}
			else{
				if (mpi_jindex!=0)AirMedia->UpdateBC_velocity(SouthBC);
				else{
					south_pore->Cal_velocity();
					south_pore->UpdateBC_velocity(SouthBC);
					south_pore->UpdateBC_velocity(EastBC);
					south_pore->UpdateBC_velocity(WestBC);
					AirMedia->UpdateBC_velocity(SouthBC,*south_pore);
				}
			}
			AirMedia->Update_PML_Qvw();
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void airpore::Cal_pressure(int time_judge,int temp_time)
{
	AirMedia->Cal_pressure();
	switch(CaseNo){
	case 1:
		{
			AirMedia->UpdateBC_pressure(NorthBC,time_judge,temp_time);
			AirMedia->UpdateBC_pressure(EastBC,time_judge,temp_time);
			AirMedia->UpdateBC_pressure(WestBC,time_judge,temp_time);
			if(air_boundary.air_south==rigid || air_boundary.air_south == TDBC){
				AirMedia->UpdateBC_pressure(SouthBC,time_judge,temp_time);
			}else{
				if (mpi_jindex!=0)AirMedia->UpdateBC_pressure(SouthBC,time_judge,temp_time);
				else{
					south_pore->Cal_pressure();
					south_pore->UpdateBC_pressure(SouthBC,time_judge,temp_time);
					south_pore->UpdateBC_pressure(EastBC,time_judge,temp_time);
					south_pore->UpdateBC_pressure(WestBC,time_judge,temp_time);
					AirMedia->UpdateBC_pressure(SouthBC,*south_pore);
				}
			}
			AirMedia -> Update_PML_Qp();

		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void airpore::get_output(int mpi_rank1,int mpi_size1)
{
	//initialize all computational domains
	mpi_rank=mpi_rank1;
	mpi_size=mpi_size1;
	mpi_iindex=mpi_rank/mpi_zarea;
	mpi_jindex=mpi_rank%mpi_zarea;
	SetInitialCond();
	//output coordinate (combine all domains into one coordinate system)
	//output_coordi();//for plot3d in Tecplot

	//main calculation
	int *temp_i,*temp_j,IMAX,JMAX,NO_Rec_y;
	int i,j,mm,index_pressure;
	double **pressure_temp;
	clock_t begin_calp, begin_calv, end_calv, end_temp;
	Rec_count=10;
	pressure_temp=new double *[Rec_count];
	for(i=0;i<Rec_count;i++) pressure_temp[i]=new double[time_domain];
	for (i=0;i<Rec_count;i++){
		for (j=0;j<time_domain;j++){
			pressure_temp[i][j]=0.0;
		}
	}
	IMAX=AirMedia->get_whole_IMAX();
	JMAX=AirMedia->get_whole_JMAX();
	NO_Rec_y=Rec_count;
	temp_i=new int [NO_Rec_y];temp_j=new int [NO_Rec_y];
	
	for (i = 0; i < NO_Rec_y; i++){
		int ii, jj;
		struct Position rec;
		if (i <= 10) {
			rec.y = 35;
			rec.z = i * 0.1;
		}
		AirMedia->get_position(ii, jj, rec);
		if (ii == 0 && jj == 0){ ii = 1; jj = 1; }
		temp_i[i] = ii; temp_j[i] = jj;
		if (ii > IMAX || jj > JMAX ){ ii = 1; jj = 1; cout << "receiver # " << i << " out of domain! rest!" << endl; }
		temp_i[i]=ii;temp_j[i]=jj;
	}
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank == 0) cout << "receiver set" << endl;
		_mkdir("solution");
	//output of p_t(y) -- z is fixed z=temp_j
	//if pty results files exists, remove them
	if (mpi_rank == 0){
		for (i = 0; i < NO_Rec_y; i++){
			char tempnameid[10] = "";
			char ptyfilename[100] = ("solution\\pty");

			if (i < 10 && i >= 0) strcat_s(ptyfilename, "00");
			else if (i >= 10 && i < 100) strcat_s(ptyfilename, "0");
			sprintf_s(tempnameid, "%d", i);
			strcat_s(ptyfilename, tempnameid);// pty001
			strcat_s(ptyfilename, ".dat");	// pty001.dat

			ofstream ofout(ptyfilename, ios::out | ios::binary);
			if (ofout.is_open()) remove(ptyfilename);
			ofout.close();
		}
		
		//_mkdir("/panfs/pfs.local/work/zczheng/j391z772/Results/2D/20171020_2D_TDBC_study");
	}
	if (restart==0){
		index_pressure=0;
	}

	int MF_count=0,temp_time=0;//MF_count to count how many frames when moving
	int time_series1,time_series2,time_series3;
	int time_series4,time_series5,time_series6;
	time_series1=3100000;time_series2=3100000;time_series3=3100000;
	time_series4=3100000;time_series5=3100000;time_series6=3100000;

	int int_pois,*point_pois;
	if(restart==1){
		char restart_temp1[20];
		sprintf_s(restart_temp1,"%d",mpi_rank);
		ifstream infile;
                strcat_s(restart_infile,restart_temp1);

		strcat_s(restart_infile,".dat");
		infile.open(restart_infile,ios::in|ios::binary);
		infile>>MF_count>>temp_time>>index_pressure;

		infile.ignore(100,'\n');
		int_pois=infile.tellg();
		point_pois=&int_pois;
		infile.close();

		ifstream infile1("pty1.dat",ios::in|ios::binary);
		ifstream infile2("pty2.dat",ios::in|ios::binary);
		ifstream infile3("pty3.dat",ios::in|ios::binary);
		ifstream infile4("pty4.dat",ios::in|ios::binary);
		ifstream infile5("pty5.dat",ios::in|ios::binary);
		ifstream infile6("pty6.dat",ios::in|ios::binary);
		ifstream infile7("pty7.dat",ios::in|ios::binary);

		if (mpi_rank==0){
			for(j=0;j<index_pressure;j++){
				if(j<time_series1){
					infile1>>pressure_temp[0][j];
					infile1.ignore(100,'\n');
				}
				if (j<time_series2){
					infile2>>pressure_temp[1][j];
					infile2.ignore(100,'\n');
				}
				if (j<time_series3){
					infile3>>pressure_temp[2][j];
					infile3.ignore(100,'\n');
				}

				if(j<time_series4){
					infile4>>pressure_temp[3][j];
					infile4.ignore(100,'\n');
				}

				if (j<time_series5){
	           			infile5>>pressure_temp[4][j];
			        	infile5.ignore(100,'\n');
            			}
				if (j<time_series6){
					infile6>>pressure_temp[5][j];
					infile6.ignore(100,'\n');
				}
	        	infile7>>pressure_temp[6][j];
		    	infile7.ignore(100,'\n');
			}
		}
		infile1.close();infile2.close();infile3.close();
                infile4.close();infile5.close();infile6.close();
		infile7.close();
	}

	ofstream outfile[10];
	if (mpi_rank == 0){
		for (i = 0; i < NO_Rec_y; i++){
			char tempnameid1[10] = "";
			char pty_infile[100] = ("solution\\pty");

			if (i < 10 && i >= 0) strcat_s(pty_infile, "00");
			else if (i >= 10 && i < 100) strcat_s(pty_infile, "0");
			sprintf_s(tempnameid1, "%d", i);
			strcat_s(pty_infile, tempnameid1);	// pty001
			strcat_s(pty_infile, ".dat");			// pty001.dat
			if (mpi_rank != 0)sprintf_s(pty_infile, "%d", mpi_rank);

			outfile[i].open(pty_infile, ios::app);		
			outfile[i].setf(ios::scientific, ios::floatfield);
			outfile[i].precision(6);
			if (mpi_rank != 0){
				outfile[i].close();
				remove(pty_infile);
			}
		}
		cout << "Receiver files all set!" << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	double temp,temp1,temp2,temp3,temp4;
	int n;
	if (temp_time>=time_domain){  //judge to use the time boundary
		time_judge=1;
	}else{
		time_judge=0;
	}
	do {
		//when count=0, computation will cover whole frame; next, computation will cover half of frame
		int MF_limit;
		double MF_lefttime1;
		double MF_lefttime2;
		if(MF_count==(int)((IMAX-move_frame.IMAX)/(move_frame.lead_DI-1)+1.e-6)||(IMAX==move_frame.IMAX)){
			MF_limit=time_domain-temp_time;
		}else if(MF_count==0){
			MF_limit=(int)((move_frame.IMAX-1)*AirStep.diff_y/speed_moveframe/AirStep.diff_t*3.5/4.0);
		}else{
			MF_limit=(int)((move_frame.trail_DI-1)*AirStep.diff_y/speed_moveframe/AirStep.diff_t);
			MF_lefttime1=0.0;
			MF_lefttime2=0.0;
			if (MF_count%10==0){
				MF_lefttime1=(move_frame.trail_DI-1)*AirStep.diff_y/speed_moveframe/AirStep.diff_t-MF_limit;
				MF_lefttime1=MF_lefttime1*10.0;
			}
			if(MF_count%100==0){
				MF_lefttime2=(move_frame.trail_DI-1)*AirStep.diff_y/speed_moveframe/AirStep.diff_t-MF_limit;
				MF_lefttime2=MF_lefttime2*100-int(MF_lefttime2*10)*10;
			}
			MF_limit=MF_limit+int(MF_lefttime1)+int(MF_lefttime2);
			if ((temp_time+MF_limit)>time_domain)MF_limit=time_domain-temp_time;
		}

		if (restart==1){
			//input_restartfile(point_pois,MF_count);
			restart=0;
		}
		//update Y,Z and v,w,p for moving frame(MF)
		SetMovingFrame(MF_count);

		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank == 0) cout << "Moing frame set" << ", MF_limit = " << MF_limit << endl;

		for(n=0;n<MF_limit;n++){
			//print out data
			index_pressure=index_pressure+1;
			
			for (i = 0; i < NO_Rec_y; i++){
				temp = AirMedia->get_pressure(temp_i[i], temp_j[i], MF_count);
				if (mpi_rank == 0){
					if (fabs(temp)<1.0e-100) temp = 0.0;
					outfile[i] << temp << endl;
					pressure_temp[i][index_pressure - 1] = temp;
				}
			}
			
			if (mpi_rank == 0) cout << "Time step = " << index_pressure << ", receiver #2 reading is " << pressure_temp[1][index_pressure - 1];
			
			MPI_Barrier(MPI_COMM_WORLD);
			begin_calv = clock();
			//calculate pressure and velocity to advance
			if((scheme==FB_p_v)||(scheme==LeapTrap)){
				
				Cal_pressure(time_judge,temp_time+n);
				Cal_velocity();

			}else if(scheme==FB_v_p){
				Cal_velocity();
				Cal_pressure(time_judge,temp_time+n);
			}else{
				Cal_pressure(time_judge,temp_time+n);
				Cal_velocity();
				UpdateInitialCond(MF_count);
				Cal_pressure(time_judge,temp_time+n);
				Cal_velocity();
			}
			//update source 
				end_calv = clock();
				if (mpi_rank == 0) 	cout << ", one time step using time: " << ((double(end_calv - begin_calv)) / CLOCKS_PER_SEC) << "s." << endl;
			UpdateSource(f_0, index_pressure, source);


			UpdateInitialCond(MF_count);
			if ((n + temp_time + 1) % out_difft == 0){//output of contour
				if (mpi_iindex == 0 && mpi_jindex != 0)mpisend_data_contour();
				if (mpi_rank == 0) get_data_contour(n + temp_time);
			}
			if (temp_time+n+1>=time_domain && time_judge==1) time_judge=0;
		}
		MF_count=MF_count+1;
		temp_time=MF_limit+temp_time;
	}while((MF_count<=(int)((IMAX-move_frame.IMAX)/(move_frame.lead_DI-1)+1.e-6))&&
		(temp_time<time_domain));
	if (mpi_rank==0){
		for (i = 0; i < NO_Rec_y; i++) outfile[i].close();
	}
	//delete new variabe to vacate memory
	delete AirMedia;
	switch(CaseNo){
	case 1:
		if (mpi_jindex==0){
			if(air_boundary.air_south==porous_media) delete south_pore;
		}
		break;
	default:
		cout<<"it is wrong to delete object";
	}
	delete[] temp_i;delete[] temp_j;
	for (i = 0; i < Rec_count;i++)delete[] pressure_temp[i];
	delete[] pressure_temp;
}


void airpore::UpdateInitialCond(int MF_count)
{
	AirMedia->UpdateInitialCond(MF_count);
	switch(CaseNo){
			case 1:
			if (mpi_jindex==0){
				if(air_boundary.air_south==porous_media) south_pore->UpdateInitialCond(MF_count);
			}
				break;
			default:
				cout<<"it is wrong to entrance boundary";
	}
}
void airpore::save_restartfile(char *restartfile)
{
	AirMedia->save_restart_cal(restartfile);
	AirMedia->save_restart_air(restartfile);
	switch(CaseNo){
	case 1:
		if (mpi_jindex==0){
			if (air_boundary.air_south==porous_media) south_pore->save_restart_cal(restartfile);
		}
		break;
	case 2:
		{
			north_pore->save_restart_cal(restartfile);
			if(AbsorbPara.Cs<0){
				north_pore->save_restart_air(restartfile);
				north_pore->save_restart_pml(restartfile);
			}
			if(air_boundary.air_south!=rigid) south_pore->save_restart_cal(restartfile);
			if((AbsorbPara.Cs<0)&&(air_boundary.air_south==absorbing)){
				south_pore->save_restart_air(restartfile);
				south_pore->save_restart_pml(restartfile);
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}
void airpore::input_restartfile(int *point_pois,int MF_count)
{
	AirMedia->input_restart_cal(restart_infile,point_pois,MF_count);
	AirMedia->input_restart_air(restart_infile,point_pois);
	switch(CaseNo){
	case 1:
		if (mpi_jindex==0){
			if (air_boundary.air_south==porous_media) south_pore->input_restart_cal(restart_infile,point_pois,MF_count);
		}
		break;
	case 2:
		{
			north_pore->input_restart_cal(restart_infile,point_pois,MF_count);
			if(AbsorbPara.Cs<0){
				north_pore->input_restart_air(restart_infile,point_pois);
				north_pore->input_restart_pml(restart_infile,point_pois);
			}
			if(air_boundary.air_south!=rigid) south_pore->input_restart_cal(restart_infile,point_pois,MF_count);
			if((AbsorbPara.Cs<0)&&(air_boundary.air_south==absorbing)){
				south_pore->input_restart_air(restart_infile,point_pois);
				south_pore->input_restart_pml(restart_infile,point_pois);
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}

}


void airpore::SetMovingFrame(int MF_count)
{
	AirMedia->SetMovingFrame(MF_count);
	AirMedia->SetWindProfile(MF_count);
	AirMedia->SetQMove(MF_count);
	switch(CaseNo){
	case 1:
		if (mpi_jindex==0){
			if(air_boundary.air_south==porous_media) south_pore->SetMovingFrame(MF_count);
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void airpore::UpdateSource(double f0, int N, Position)
{
	AirMedia->UpdateSource(f0, N, source);
}

void airpore::SetInitialCond()
{
	char filename[100];
	strcpy_s(filename,"coordi_air.dat");
	AirMedia= new air(scheme,AirStep,filename,move_frame,gauss_width,
		velocity_coef,velocity_method,AirPara,AbsorbPara.resistivity,PorePara,PoreParb,hill1,
		vortex1,mpi_rank,mpi_size,mpi_yarea,mpi_zarea,0);
	AirMedia->Set_InitialCond(source);
	switch(CaseNo){
	case 1:
		if (mpi_jindex==0){
			if(air_boundary.air_south==porous_media){
				strcpy_s(filename,"coordi_pore.dat");
				south_pore=new porous(scheme,SouthStep,filename,move_frame,gauss_width,PorePara,AirPara,
														vortex1,mpi_rank,mpi_size,mpi_yarea,1,mpi_zarea);
				south_pore->Set_InitialCond(source);
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void airpore::get_FFT_y1(int out_type)
{
	int i, j, num;
	double dt;
	num = 51;
	dt = AirStep.diff_t;

	for (j = 1; j <= num; j++){
		char pty_file[100] = "solution\\pty";
		char p_f_file[100] = "solution\\p_f2_";
		char temp[100] = "";
		sprintf_s(temp, "%d", j);
		strcat_s(temp, ".dat");
		strcat_s(pty_file, temp);
		cout << "start to read " << pty_file << endl;
		ifstream infile(pty_file, ios::in | ios::binary);
		if (!infile) { cout << "Cannot open file " << pty_file << " for read!!!" << endl; return; }
		double *ppr, *pr, *pi;
		ppr = new double[time_domain];
		pi = new double[FFT_N];
		pr = new double[FFT_N];
		for (i = 0; i < time_domain; i++)ppr[i] = 0.0;
		for (i = 0; i < FFT_N; i++){
			pr[i] = 0.0;
			pi[i] = 0.0;
		}
		for (i = 0; i < time_domain; i++){
			infile >> ppr[i];
			infile.ignore(100, '\n');
		}
		infile.close();
		strcat_s(p_f_file, temp);
		cout << "start to write " << p_f_file << endl;
		if (receiver.y != 0){
			for (i = 0; i < time_domain; i++) { pr[i] = ppr[i]; pi[i] = 0; }
			for (i = time_domain; i < FFT_N; i++) { pr[i] = 0; pi[i] = 0; }
			FFT(1, FFT_m, pr, pi);
			FFT_output(FFT_m, pr, pi, p_f_file, dt);
		}
		delete[] pr; delete[] pi; delete[] ppr;
	}
}

void airpore::mpisend_data_contour(){
	AirMedia->mpi_send_data();
}

//the following is the output format for the command "load data files" in Tecplot
void airpore::get_data_contour(int n)
{
	char p_contour[200]="pt",temp[20];
	int nn;
	nn=(n+1)/out_difft;
	if(nn<10) strcat_s(p_contour,"00");
	if((nn>=10)&&(nn<100)) strcat_s(p_contour,"0");
	int ii;
	ii=sprintf_s(temp,"%d",nn);
	strcat_s(p_contour,temp);strcat_s(p_contour,".dat");

	//output
	ofstream fp_p_t(p_contour,ios::out | ios::trunc | ios::binary);
	fp_p_t.setf(ios::scientific,ios::floatfield);
	fp_p_t.precision(6);

	switch(CaseNo){
	case 1:
		{
			if(air_boundary.air_south==rigid || air_boundary.air_south==TDBC){
				fp_p_t << "VARIABLES = \"Y\", \"Z\", \"P\", \"V\", \"W\"" << endl;
				fp_p_t << "ZONE T = \"Air Zone\",";
				fp_p_t<<*AirMedia;
			}else{

				fp_p_t<<"ZONE T = \"Air Zone\",";
				fp_p_t<<*AirMedia;
				fp_p_t<<"ZONE T = \"south pore Zone\",";
				fp_p_t<<*south_pore;
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
	fp_p_t.close();
}
