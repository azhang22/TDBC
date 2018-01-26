//++++++++++++++++++++++++++++++filename: air.h +++++++++++++++++++++++++++++++++++++++//

//----------------------------description of variable----------------------------------//
/***************************************************************************************/
/* Vav,Wav is average speed. Vav_over_y means patial (Vav)/patial (y)
/* Pav is average atmospheric pressure
/* Aav is a reversal of average density
/* adiabatic_coef is adiabatic coefficient
/* Profile_Vav,Profile_Wav point to the average wind speed profile
/**************************************************************************************/
#ifndef _AIR_H
#define _AIR_H
#include "mygeneral.h"
#include "calculation_2D.h"

//**************************wind speed profile*****************************//
//z is with respect to local ground level
class SpeedProfile{
private:
	int velocity_method;
	double y0,z0,L,alpha;//m
	double b,c,circulation,vorticity;//b is Mach number for vortex back flow
public:
	SpeedProfile(){};
	SpeedProfile(double velocity_coef,int VM);
	~SpeedProfile();
	double Vav(double x,double y,double z);
	double Wav(double x,double y,double z);
};

class air:public calculation_2D
{
	protected:
		double **Vav_over_y,**Vav_over_z,**Wav_over_y,**Wav_over_z,**Vav,**Wav;
		double Pav,adiabatic_coef;//sound_speed, Aav;
		int velocity_method;
		SpeedProfile *SP;

		double **Qv,**Qw,**Qp;
		double PML_AbsorbZmax,PML_alpha;
		double PML_Height1,PML_z1,PML_Height2,PML_z2;

		// PML boundary prameters for the y direction
		double PML_AbsorbYmax,PML_alphay;
		double PML_Width1,PML_y1,PML_Width2,PML_y2;

		double **whole_Qv,**whole_Qw,**whole_Qp;

		// this is added for porous media hill
		int cr_judge,cr_points;
		double ystar,ystop,zstar,zstop;
		double cr_a,cr_b,cr_c,cr_d,cr_e,cr_f,cr_g,cr_h,cr_i,cr_j,cr_k,cr_l;

		// for the hill
		double *cr_y,*cr_z;
		int **cr_uc;

		double crpo_eff_density,crpo_porosity,crpo_resistivity;
		double crpo_eff_density_a,crpo_porosity_a,crpo_resistivity_a;
		double crpo_eff_density_b,crpo_porosity_b,crpo_resistivity_b;
		double crpo_Kp,crpo_Beta,crpo_Beta1,crpo_Beta2,crpo_Gama;
		double crpo_Kp_a,crpo_Beta_a,crpo_Gama_a;
		double crpo_Kp_b,crpo_Beta_b,crpo_Gama_b;
		
		//vortex
		double ystar_vor,ystop_vor,zstar_vor,zstop_vor;	
	public:
		air(){};
		air(scheme_list,DifferenceStep,char *coordi,MovingFrame MF,const double,
			double velocity_coef,int velocity_method1,AirStruct AirPara,double AbsorbZmax,
			PoreStruct PorePara,PoreStruct PoreParb,hillpore hill1,vortex vortex1,int mpi_rank1,int mpi_size1,int mpi_yarea1,
			int mpi_zarea1,int mpi_porous1);
		~air();
		virtual void cal_fvw(double **temp_v,double **temp_w,
			double **temp_p);
		virtual void cal_fp(double **temp_v,double **temp_w,
			double **temp_p);
		virtual void save_restart_air(char *restartfile);
		virtual void input_restart_air(char *restartfile,int *point_pois);
		virtual void SetWindProfile(int N);
		virtual void UpdateBC_pressure(boundary_location BC,int time_judge,int time_current);
		virtual void UpdateBC_velocity(boundary_location BC);
		//virtual void Cal_VWav(int index_pressure);
		//virtual void cal_vwav(int T);

		virtual void Update_PML_Qvw();
		virtual void Update_PML_Qp();
		virtual void SetQMove(int N);
		double weno_scheme(double wen_u1,double wen_u2,double wen_u3,
								double wen_u4,double wen_u5,double wen_u6,double ddx);
};

#endif
