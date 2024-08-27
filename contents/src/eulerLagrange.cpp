// Configured==true
#include "eulerLagrange.h"

const cpp_dec_float_25 epsilon_0{ 1. / (4.*M_PI) }; // Geometrised-Gaussian units


Precalculated::Precalculated(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda):
	sigma_p2{ sigma_p*sigma_p },
	sigma_p3{ sigma_p*sigma_p2 },
	sigma_p4{ sigma_p2*sigma_p2 },
	delta_p2{ delta_p*delta_p },
	delta_p3{ delta_p*delta_p2 },
	delta_p4{ delta_p2*delta_p2 },
	chi_p2{ chi_p*chi_p },
	chi_p3{ chi_p*chi_p2 },
	chi_p4{ chi_p2*chi_p2 },

	tdot{ (p->t - p->t_prev) / dlambda },
	rdot{ (p->r - p->r_prev) / dlambda },
	phidot{ (p->phi - p->phi_prev) / dlambda },
	thetadot{ (p->theta - p->theta_prev) / dlambda },

	tsdot{ (ps->t - ps->t_prev) / dlambda },
	rsdot{ (ps->r - ps->r_prev) / dlambda },
	phisdot{ (ps->phi - ps->phi_prev) / dlambda },
	thetasdot{ (ps->theta - ps->theta_prev) / dlambda },
	tsddot{
		 (ps->t - 2.*ps->t_prev + ps->t_prevprev) / (dlambda*dlambda)},
	rsddot{
		(ps->r - 2.*ps->r_prev + ps->r_prevprev) / (dlambda*dlambda)},
	phisddot{
		(ps->phi - 2.*ps->phi_prev + ps->phi_prevprev) / (dlambda*dlambda)},
	thetasddot{
		(ps->theta - 2.*ps->theta_prev + ps->theta_prevprev) / (dlambda*dlambda)},

	nx{p->r*sin(p->theta)*cos(p->phi) - ps->r*sin(ps->theta)*cos(ps->phi)},
	ny{p->r*sin(p->theta)*sin(p->phi) - ps->r*sin(ps->theta)*sin(ps->phi)},
	nz{p->r*cos(p->theta) - ps->r*cos(ps->theta)},
	nmag{sqrt(nx*nx + ny*ny + nz*nz)},
	nmag2{nmag*nmag},
	nmag5{nmag*nmag2*nmag2},
	nmag6{nmag5*nmag},
	nmag7{nmag6*nmag},
	nmag8{nmag6*nmag2},
	nmag12{nmag6*nmag6},
	nmag13{nmag6*nmag7},
	nmag14{nmag7*nmag7},

	ral_sqr{BH->a2 + BH->l2 + p->r2},
	sintheta{sin(p->theta)},
	sintheta2{sintheta*sintheta},
	sintheta3{sintheta2*sintheta},
	sintheta4{sintheta2*sintheta2},
	sintheta5{sintheta2*sintheta3},
	sintheta6{sintheta3*sintheta3},
		sintheta7{sintheta3*sintheta4},

	costheta{cos(p->theta)},
	sin2theta{sin(2.*p->theta)},
	sinphi{sin(p->phi)},
	cosphi{cos(p->phi)},

	sinthetas{sin(ps->theta)},
	costhetas{cos(ps->theta)},
	sinphis{sin(ps->phi)},
	cosphis{cos(ps->phi)},
	sinthetasdot{sin(thetasdot)},
	costhetasdot{cos(thetasdot)},
	sinphisdot{sin(phisdot)},
	cosphisdot{cos(phisdot)},
	sinthetas2{sinthetas*sinthetas},
	costhetas2{costhetas*costhetas},
	sinphis2{sinphis*sinphis},
	cosphis2{cosphis*cosphis},
	rsdot2{rsdot*rsdot},
	phisdot2{phisdot*phisdot}
{
	/*
	Initalise the struct.
	*/
	return;
}
cpp_dec_float_25 calc_t_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const Precalculated* precalculated,
	const cpp_dec_float_25 dlambda)
{
	/*
	Calculate the next value of t.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_25 sigma_p: Precalculated sigma value.
		- const cpp_dec_float_25 delta_p: Precalculated delta value.
		- const cpp_dec_float_25 chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_25 dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_25 t_new: New coordinate.
	*/

	cpp_dec_float_25 t_new{};
	cpp_dec_float_25 temp[40];
	temp[0] = ps->epsilon*p->r;
	temp[1] = -BH->a+precalculated->phisdot*(BH->a2+p->r2);
	temp[2] = p->mass*precalculated->nmag12*precalculated->sigma_p2;
	temp[3] = -precalculated->nmag6+ps->sigma6;
	temp[4] = precalculated->ny*p->r*precalculated->cosphi-precalculated->nx*p->r*precalculated->sinphi;
	temp[5] = (precalculated->rdot*(precalculated->nx*precalculated->cosphi+precalculated->ny*precalculated->sinphi))/precalculated->nmag+(precalculated->phidot*temp[4])/precalculated->nmag+(precalculated->rsdot*(-(precalculated->nx*precalculated->cosphis)-precalculated->ny*precalculated->sinphis))/precalculated->nmag+(precalculated->phisdot*(-(precalculated->ny*ps->r*precalculated->cosphis)+precalculated->nx*ps->r*precalculated->sinphis))/precalculated->nmag;
	temp[6] = precalculated->rdot*ps->sigma6;
	temp[7] = -2.*BH->mass+2.*p->r;
	temp[8] = temp[6]*temp[3];
	temp[9] = p->mass*precalculated->nmag12*delta_p*sigma_p;
	temp[10] = ps->epsilon*ps->sigma6;
	temp[11] = temp[10]*temp[3];
	temp[12] = p->mass*precalculated->nmag13*sigma_p;
	temp[13] = -8.*temp[0]*(BH->a2+p->r2)*temp[1]*temp[6]*temp[3];
	temp[14] = precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot;
	temp[15] = p->mass*precalculated->nmag12*sigma_p;
	temp[16] = 0.+precalculated->phidot*(BH->a*(BH->a2+p->r2)+precalculated->phisdot*(-BH->a4-2.*BH->a2*p->r2-p->r4))+BH->a*temp[1]*precalculated->tdot;
	temp[17] = precalculated->nmag14*delta_p;
	temp[18] = temp[1]*ps->sigma6;
	temp[19] = 48.*ps->epsilon*(BH->a2+p->r2)*temp[18]*temp[3]*temp[5];
	temp[20] = BH->charge*p->r2;
	temp[21] = p->mass*precalculated->sigma_p2;
	temp[22] = -precalculated->phidot*(BH->a2+p->r2)+BH->a*precalculated->tdot;
	temp[23] = (8.*temp[0]*temp[8]*(BH->a*temp[1]+delta_p))/temp[2];
	temp[24] = (2.*temp[7]*precalculated->rdot*precalculated->tdot)/sigma_p;
	temp[25] = (4.*temp[11]*(BH->a*precalculated->phisddot*(BH->a2-p->r2)*delta_p-temp[7]*precalculated->rdot*(BH->a*temp[1]-2.*delta_p)))/temp[9];
	temp[26] = (2.*BH->a2)/sigma_p-(2.*delta_p)/sigma_p;
	temp[27] = (4.*p->r*(BH->a2-p->r2)*precalculated->rdot*temp[14])/precalculated->sigma_p2;
	temp[28] = 48.*temp[10]*(-0.5*precalculated->nmag6+ps->sigma6)*(temp[16]*delta_p+precalculated->tdot*precalculated->delta_p2-precalculated->rdot*precalculated->rsdot*precalculated->sigma_p2)*temp[4];
	temp[29] = (4.*BH->a*p->r*precalculated->rdot*temp[22])/precalculated->sigma_p2;
	temp[30] = temp[23]-(p->charge*BH->charge*precalculated->rdot)/(p->mass*sigma_p);
	temp[31] = temp[25]-(24.*temp[10]*(BH->a*temp[1]-delta_p)*temp[5])/(p->mass*precalculated->nmag7*sigma_p);
	temp[32] = temp[27]-(4.*ps->epsilon*precalculated->phisddot*pow(BH->a2-p->r2,2)*ps->sigma6*temp[3])/temp[15];
	temp[33] = 2.*p->charge*temp[20]*precalculated->rdot;
	temp[34] = temp[30]-temp[24];
	temp[35] = temp[7]*temp[8];
	temp[36] = 48.*temp[11]*(BH->a*temp[1]+delta_p)*temp[5];
	temp[37] = -2.*BH->a*(BH->a2+p->r2)*(temp[13]/temp[2]-temp[32]-temp[28]/(p->mass*temp[17]*sigma_p)-(24.*ps->epsilon*(BH->a2+p->r2)*temp[18]*temp[5])/(p->mass*precalculated->nmag7*sigma_p)-temp[19]/temp[12]);
	temp[38] = temp[29]-(4.*p->r*precalculated->rdot*precalculated->tdot*delta_p)/precalculated->sigma_p2;
	temp[39] = temp[33]/temp[21]-temp[38]+temp[34]+(4.*ps->epsilon*temp[35]*(BH->a*temp[1]+delta_p))/temp[9]-temp[31]+temp[36]/temp[12];


	t_new = 2. * p->t - p->t_prev + dlambda*dlambda*(-((temp[37]/sigma_p-(2.*pow(BH->a2+p->r2,2)*temp[39])/sigma_p)/((4.*BH->a2*pow(BH->a2+p->r2,2))/precalculated->sigma_p2-(2.*pow(BH->a2+p->r2,2)*temp[26])/sigma_p)));
	return t_new;
}

cpp_dec_float_25 calc_r_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const Precalculated* precalculated,
	const cpp_dec_float_25 dlambda)
{
	/*
	Calculate the next value of r.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_25 sigma_p: Precalculated sigma value.
		- const cpp_dec_float_25 delta_p: Precalculated delta value.
		- const cpp_dec_float_25 chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_25 dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_25 r_new: New coordinate.
	*/

	cpp_dec_float_25 r_new{};
	cpp_dec_float_25 temp[66];
	temp[0] = precalculated->rdot*precalculated->rdot;
	temp[1] = -2.*BH->mass+2.*p->r;
	temp[2] = ps->epsilon*ps->sigma6;
	temp[3] = precalculated->nmag12*delta_p;
	temp[4] = p->r*precalculated->rdot;
	temp[5] = temp[1]*precalculated->rdot;
	temp[6] = temp[4]*precalculated->rsdot;
	temp[7] = -precalculated->nmag6+ps->sigma6;
	temp[8] = BH->a-precalculated->phisdot*(BH->a2+p->r2);
	temp[9] = precalculated->rdot*precalculated->rsdot;
	temp[10] = ps->sigma6*temp[7];
	temp[11] = -BH->a+precalculated->phisdot*(BH->a2+p->r2);
	temp[12] = (precalculated->phidot*(BH->a2+p->r2)*temp[11]+BH->a*temp[8]*precalculated->tdot)*delta_p-precalculated->tdot*precalculated->delta_p2+temp[9]*precalculated->sigma_p2;
	temp[13] = precalculated->ny*p->r*precalculated->cosphi-precalculated->nx*p->r*precalculated->sinphi;
	temp[14] = -(precalculated->ny*ps->r*precalculated->cosphis)+precalculated->nx*ps->r*precalculated->sinphis;
	temp[15] = precalculated->nx*precalculated->cosphi+precalculated->ny*precalculated->sinphi;
	temp[16] = 2.*p->r*temp[0];
	temp[17] = p->charge*BH->charge;
	temp[18] = p->mass*precalculated->sigma_p2;
	temp[19] = (2.*p->r*pow(precalculated->tdot,2)*delta_p)/precalculated->sigma_p2;
	temp[20] = (temp[1]*temp[0]*sigma_p)/precalculated->delta_p2;
	temp[21] = -BH->a3-BH->a*p->r2+precalculated->phisdot*(BH->a4+2.*BH->a2*p->r2+p->r4);
	temp[22] = 0.+2.*precalculated->rdot*precalculated->rsdot;
	temp[23] = -8.*temp[6]*delta_p+(4.*temp[5]*precalculated->rsdot-4.*precalculated->rsddot*delta_p)*sigma_p;
	temp[24] = -4.*temp[5]*precalculated->rsdot+4.*precalculated->rsddot*delta_p;
	temp[25] = (8.*ps->epsilon*p->r*temp[10]*temp[12])/(p->mass*temp[3]*precalculated->sigma_p2);
	temp[26] = (24.*temp[2]*temp[12]*temp[15])/(p->mass*precalculated->nmag8*delta_p*sigma_p);
	temp[27] = precalculated->rsdot*ps->sigma12;
	temp[28] = (48.*precalculated->phidot*temp[13])/precalculated->nmag;
	temp[29] = 48.*precalculated->phisdot*temp[14];
	temp[30] = precalculated->rsdot*ps->sigma6;
	temp[31] = (24.*precalculated->phidot*temp[13])/precalculated->nmag;
	temp[32] = 24.*precalculated->phisdot*temp[14];
	temp[33] = (2.*temp[17]*p->r2*precalculated->tdot)/temp[18];
	temp[34] = temp[19]-(p->charge*BH->charge*precalculated->tdot)/(p->mass*sigma_p);
	temp[35] = temp[1]*(precalculated->phidot*temp[21]+precalculated->tdot*(BH->a*(BH->a+precalculated->phisdot*(-BH->a2-p->r2))-2.*delta_p))+2.*p->r*temp[22]*sigma_p;
	temp[36] = p->mass*precalculated->nmag12*precalculated->delta_p2;
	temp[37] = p->mass*precalculated->nmag12*precalculated->delta_p2*sigma_p;
	temp[38] = precalculated->nmag14*delta_p;
	temp[39] = -48.*precalculated->rdot*temp[15];
	temp[40] = (48.*precalculated->rsdot*(+(precalculated->nx*precalculated->cosphis)+precalculated->ny*precalculated->sinphis))/precalculated->nmag;
	temp[41] = temp[31]+(24.*precalculated->rsdot*(-(precalculated->nx*precalculated->cosphis)-precalculated->ny*precalculated->sinphis))/precalculated->nmag;
	temp[42] = temp[33]-(2.*p->r*pow(precalculated->phidot*(BH->a2-p->r2)+BH->a*precalculated->tdot,2))/precalculated->sigma_p2;
	temp[43] = (ps->epsilon*ps->sigma6*temp[23])/(p->mass*precalculated->nmag6*precalculated->delta_p2);
	temp[44] = (4.*ps->epsilon*temp[1]*temp[10]*temp[12])/temp[37];
	temp[45] = p->mass*temp[38]*sigma_p;
	temp[46] = temp[39]/precalculated->nmag-temp[28]-temp[40]-temp[29]/precalculated->nmag;
	temp[47] = (24.*precalculated->rdot*temp[15])/precalculated->nmag+temp[41]+temp[32]/precalculated->nmag;
	temp[48] = temp[42]+temp[34];
	temp[49] = temp[20]+(4.*temp[2]*temp[7]*temp[35])/(p->mass*temp[3]*sigma_p);
	temp[50] = temp[25]+temp[44];
	temp[51] = temp[2]*temp[7];
	temp[52] = (ps->epsilon*temp[27]*sigma_p*temp[46])/(p->mass*precalculated->nmag13*delta_p);
	temp[53] = temp[48]-(temp[1]*pow(precalculated->tdot,2))/sigma_p;
	temp[54] = temp[6]*delta_p;
	temp[55] = temp[50]+temp[26];
	temp[56] = ps->epsilon*temp[30]*sigma_p*temp[47];
	temp[57] = temp[53]+temp[49];
	temp[58] = 8.*temp[54]+temp[24]*sigma_p;
	temp[59] = temp[51]*temp[12];
	temp[60] = p->mass*precalculated->nmag7*delta_p;
	temp[61] = temp[57]-temp[43];
	temp[62] = temp[55]+(48.*temp[59]*temp[15])/temp[45];
	temp[63] = temp[61]-(ps->epsilon*ps->sigma12*temp[58])/temp[36];
	temp[64] = temp[63]-temp[62];
	temp[65] = temp[64]-temp[52];


	r_new = 2. * p->r - p->r_prev + dlambda*dlambda*((-0.5*delta_p*(temp[16]/delta_p-temp[65]+temp[56]/temp[60]))/sigma_p);
	return r_new;
}

cpp_dec_float_25 calc_phi_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const Precalculated* precalculated,
	const cpp_dec_float_25 dlambda)
{
	/*
	Calculate the next value of phi.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_25 sigma_p: Precalculated sigma value.
		- const cpp_dec_float_25 delta_p: Precalculated delta value.
		- const cpp_dec_float_25 chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_25 dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_25 phi_new: New coordinate.
	*/

	cpp_dec_float_25 phi_new{};
	cpp_dec_float_25 temp[491];
	temp[0] = BH->a5*precalculated->nmag13;
	temp[1] = precalculated->rdot*delta_p;
	temp[2] = precalculated->nmag13*p->charge;
	temp[3] = BH->a*temp[2];
	temp[4] = 0.08333333333333333*BH->a6;
	temp[5] = p->r*precalculated->rdot;
	temp[6] = p->mass*precalculated->nmag13;
	temp[7] = temp[6]*precalculated->phidot;
	temp[8] = p->r7*precalculated->rdot;
	temp[9] = BH->a6*ps->epsilon;
	temp[10] = precalculated->phisdot*temp[5];
	temp[11] = 0.5*BH->a4*ps->epsilon;
	temp[12] = precalculated->phisdot*p->r3;
	temp[13] = 0.5*BH->a2*ps->epsilon;
	temp[14] = precalculated->phisdot*p->r5;
	temp[15] = ps->epsilon*precalculated->nmag7;
	temp[16] = temp[9]*precalculated->nmag;
	temp[17] = BH->a4*ps->epsilon;
	temp[18] = temp[12]*precalculated->rdot;
	temp[19] = ps->epsilon*precalculated->nmag;
	temp[20] = ps->sigma12*precalculated->delta_p2;
	temp[21] = temp[0]*p->charge;
	temp[22] = delta_p*sigma_p;
	temp[23] = BH->a3*temp[2];
	temp[24] = precalculated->rdot*temp[22];
	temp[25] = 0.16666666666666666*BH->a5;
	temp[26] = precalculated->rdot*ps->sigma6;
	temp[27] = 0.3333333333333333*BH->a3;
	temp[28] = p->r2*temp[26];
	temp[29] = 0.3333333333333333*BH->a3;
	temp[30] = 0.16666666666666666*BH->a;
	temp[31] = p->r4*temp[26];
	temp[32] = 0.16666666666666666*BH->a;
	temp[33] = delta_p*sigma_p;
	temp[34] = BH->mass*precalculated->nmag;
	temp[35] = ps->sigma12*temp[33];
	temp[36] = p->r2*precalculated->rdot;
	temp[37] = p->r3*precalculated->rdot;
	temp[38] = p->r4*precalculated->rdot;
	temp[39] = p->r5*precalculated->rdot;
	temp[40] = p->mass*BH->mass;
	temp[41] = precalculated->tdot*temp[22];
	temp[42] = precalculated->tdot*temp[33];
	temp[43] = temp[40]*precalculated->nmag13;
	temp[44] = temp[15]*precalculated->phisddot;
	temp[45] = p->r2*ps->sigma6;
	temp[46] = p->r4*ps->sigma6;
	temp[47] = p->r6*ps->sigma6;
	temp[48] = ps->sigma12*precalculated->delta_p2;
	temp[49] = p->r2*temp[48];
	temp[50] = BH->a2*temp[19];
	temp[51] = p->r4*temp[48];
	temp[52] = p->r6*temp[48];
	temp[53] = 0.5*temp[9]*precalculated->nmag5;
	temp[54] = precalculated->delta_p2*sigma_p;
	temp[55] = 1.5*temp[17]*precalculated->nmag5;
	temp[56] = precalculated->nx*precalculated->cosphi+precalculated->ny*precalculated->sinphi;
	temp[57] = temp[26]*temp[54];
	temp[58] = precalculated->rdot*temp[48];
	temp[59] = temp[17]*precalculated->phisdot;
	temp[60] = BH->a2*ps->epsilon;
	temp[61] = temp[48]*sigma_p;
	temp[62] = precalculated->phisdot*p->r6;
	temp[63] = 0.5*BH->a7*ps->epsilon;
	temp[64] = precalculated->phidot*ps->sigma6;
	temp[65] = p->r*precalculated->cosphi;
	temp[66] = ps->epsilon*precalculated->nmag5;
	temp[67] = delta_p*sigma_p;
	temp[68] = precalculated->nx*p->r;
	temp[69] = temp[66]*precalculated->phidot;
	temp[70] = precalculated->nmag5*precalculated->phidot;
	temp[71] = precalculated->ny*temp[65]-temp[68]*precalculated->sinphi;
	temp[72] = ps->sigma12*temp[67];
	temp[73] = ps->epsilon*precalculated->phidot;
	temp[74] = p->r2*temp[72];
	temp[75] = precalculated->phidot*precalculated->phisdot;
	temp[76] = p->r4*temp[72];
	temp[77] = ps->sigma6*precalculated->tdot;
	temp[78] = precalculated->nmag5*precalculated->phisdot;
	temp[79] = temp[67]*temp[71];
	temp[80] = temp[45]*precalculated->tdot;
	temp[81] = temp[46]*precalculated->tdot;
	temp[82] = ps->sigma12*precalculated->tdot;
	temp[83] = ps->epsilon*precalculated->phisdot;
	temp[84] = p->r2*temp[82];
	temp[85] = temp[83]*temp[84];
	temp[86] = temp[83]*p->r4;
	temp[87] = 0.5*BH->a5*temp[66];
	temp[88] = temp[69]*temp[45];
	temp[89] = temp[69]*temp[46];
	temp[90] = temp[49]*sigma_p;
	temp[91] = temp[51]*sigma_p;
	temp[92] = temp[77]*temp[54];
	temp[93] = temp[66]*temp[80];
	temp[94] = temp[66]*precalculated->phisdot;
	temp[95] = temp[82]*temp[54];
	temp[96] = temp[13]*precalculated->nmag5;
	temp[97] = precalculated->delta_p3*sigma_p;
	temp[98] = temp[11]*precalculated->nmag5;
	temp[99] = precalculated->sigma_p3*temp[71];
	temp[100] = precalculated->sigma_p3*temp[71];
	temp[101] = precalculated->rsdot*ps->sigma12;
	temp[102] = temp[36]*temp[101];
	temp[103] = temp[60]*precalculated->nmag5;
	temp[104] = ps->sigma6*delta_p;
	temp[105] = temp[36]*precalculated->rsdot;
	temp[106] = delta_p*precalculated->sigma_p3;
	temp[107] = temp[53]*precalculated->phisdot;
	temp[108] = temp[55]*precalculated->phisdot;
	temp[109] = -(precalculated->nx*precalculated->cosphis)-precalculated->ny*precalculated->sinphis;
	temp[110] = precalculated->rsdot*temp[61];
	temp[111] = temp[60]*precalculated->phisdot;
	temp[112] = ps->sigma6*temp[54];
	temp[113] = -(precalculated->ny*ps->r*precalculated->cosphis)+precalculated->nx*ps->r*precalculated->sinphis;
	temp[114] = 0.041666666666666664*temp[21];
	temp[115] = 0.08333333333333333*temp[23];
	temp[116] = 0.041666666666666664*temp[3];
	temp[117] = temp[4]*temp[7];
	temp[118] = 0.25*BH->a4*temp[7];
	temp[119] = 0.25*BH->a2*temp[7];
	temp[120] = 0.08333333333333333*temp[7];
	temp[121] = 0.16666666666666666*temp[9];
	temp[122] = temp[10]*ps->sigma6;
	temp[123] = temp[11]*precalculated->nmag7;
	temp[124] = temp[13]*precalculated->nmag7;
	temp[125] = 0.16666666666666666*temp[15];
	temp[126] = temp[8]*ps->sigma6;
	temp[127] = 0.16666666666666666*temp[16];
	temp[128] = temp[17]*precalculated->nmag;
	temp[129] = temp[50]*temp[14];
	temp[130] = 0.16666666666666666*temp[19];
	temp[131] = temp[8]*ps->sigma12;
	temp[132] = 0.020833333333333332*temp[21];
	temp[133] = temp[23]*BH->charge;
	temp[134] = 0.020833333333333332*temp[3];
	temp[135] = temp[25]*ps->epsilon;
	temp[136] = temp[26]*temp[33];
	temp[137] = BH->a5*temp[15];
	temp[138] = ps->sigma6*temp[22];
	temp[139] = BH->mass*precalculated->nmag7;
	temp[140] = temp[26]*temp[22];
	temp[141] = temp[139]*temp[31];
	temp[142] = temp[15]*p->r5;
	temp[143] = 0.16666666666666666*BH->a5;
	temp[144] = precalculated->rdot*ps->sigma12;
	temp[145] = temp[5]*temp[35];
	temp[146] = ps->epsilon*temp[34];
	temp[147] = temp[27]*temp[19];
	temp[148] = temp[32]*temp[146];
	temp[149] = temp[30]*temp[19];
	temp[150] = 0.08333333333333333*BH->a5;
	temp[151] = 0.08333333333333333*BH->a5;
	temp[152] = 0.16666666666666666*BH->a3;
	temp[153] = 0.16666666666666666*BH->a3;
	temp[154] = 0.08333333333333333*BH->a;
	temp[155] = 0.08333333333333333*BH->a;
	temp[156] = temp[4]*temp[44];
	temp[157] = precalculated->delta_p2*sigma_p;
	temp[158] = 0.25*BH->a2*temp[44];
	temp[159] = 0.08333333333333333*temp[44];
	temp[160] = 0.08333333333333333*temp[16];
	temp[161] = temp[48]*sigma_p;
	temp[162] = precalculated->phisddot*temp[49];
	temp[163] = temp[51]*sigma_p;
	temp[164] = temp[19]*precalculated->phisddot;
	temp[165] = temp[57]*temp[56];
	temp[166] = temp[28]*temp[54];
	temp[167] = temp[60]*temp[78];
	temp[168] = 0.5*temp[66]*temp[62];
	temp[169] = (BH->a6*temp[83]*temp[58]*sigma_p*temp[56])/precalculated->nmag;
	temp[170] = (3.*temp[111]*temp[38]*temp[61]*temp[56])/precalculated->nmag;
	temp[171] = temp[63]*precalculated->nmag5;
	temp[172] = 0.5*BH->a8*temp[69];
	temp[173] = ps->sigma6*temp[79];
	temp[174] = temp[67]*temp[71];
	temp[175] = precalculated->phisdot*temp[45];
	temp[176] = temp[55]*temp[75];
	temp[177] = temp[13]*temp[70];
	temp[178] = temp[47]*temp[174];
	temp[179] = (BH->a8*temp[73]*precalculated->phisdot*temp[72]*temp[71])/precalculated->nmag;
	temp[180] = temp[75]*temp[74];
	temp[181] = (BH->a3*temp[73]*temp[76]*temp[71])/precalculated->nmag;
	temp[182] = (BH->a2*temp[73]*temp[62]*temp[72]*temp[71])/precalculated->nmag;
	temp[183] = temp[98]*temp[80];
	temp[184] = temp[94]*temp[80];
	temp[185] = BH->a3*temp[94];
	temp[186] = (BH->a6*ps->epsilon*temp[82]*temp[67]*temp[71])/precalculated->nmag;
	temp[187] = (BH->a4*ps->epsilon*temp[84]*temp[67]*temp[71])/precalculated->nmag;
	temp[188] = BH->a3*temp[86]*temp[82]*temp[67]*temp[71];
	temp[189] = temp[54]*temp[71];
	temp[190] = 0.5*BH->a*temp[89];
	temp[191] = (BH->a5*temp[73]*temp[61]*temp[71])/precalculated->nmag;
	temp[192] = (BH->a*temp[73]*temp[91]*temp[71])/precalculated->nmag;
	temp[193] = temp[92]*temp[71];
	temp[194] = BH->a3*temp[184];
	temp[195] = 0.5*BH->a*temp[94];
	temp[196] = (2.*temp[17]*temp[95]*temp[71])/precalculated->nmag;
	temp[197] = (2.*temp[60]*temp[84]*temp[54]*temp[71])/precalculated->nmag;
	temp[198] = (BH->a*temp[86]*temp[95]*temp[71])/precalculated->nmag;
	temp[199] = temp[97]*temp[71];
	temp[200] = temp[82]*temp[97];
	temp[201] = (ps->epsilon*temp[84]*temp[97]*temp[71])/precalculated->nmag;
	temp[202] = temp[96]*temp[105];
	temp[203] = (BH->a4*ps->epsilon*precalculated->rdot*temp[101]*precalculated->sigma_p3*temp[71])/precalculated->nmag;
	temp[204] = precalculated->rdot*precalculated->rsdot;
	temp[205] = ps->epsilon*precalculated->rdot;
	temp[206] = (ps->epsilon*temp[102]*temp[106]*temp[71])/precalculated->nmag;
	temp[207] = p->r2*precalculated->rsdot;
	temp[208] = 1.5*temp[167];
	temp[209] = temp[112]*temp[109];
	temp[210] = (BH->a6*temp[83]*temp[110]*temp[109])/precalculated->nmag;
	temp[211] = (3.*temp[111]*p->r4*temp[110]*temp[109])/precalculated->nmag;
	temp[212] = temp[53]*precalculated->phisdot2;
	temp[213] = precalculated->phisdot2*temp[45];
	temp[214] = precalculated->phisdot2*temp[46];
	temp[215] = precalculated->phisdot2*temp[47];
	temp[216] = ps->epsilon*precalculated->phisdot2;
	temp[217] = (3.*temp[17]*precalculated->phisdot2*temp[90]*temp[113])/precalculated->nmag;
	temp[218] = ps->epsilon*precalculated->phisdot2*temp[52]*sigma_p*temp[113];
	temp[219] = BH->a4+2.*BH->a2*p->r2+p->r4;
	temp[220] = temp[114]*BH->charge;
	temp[221] = temp[115]*BH->charge;
	temp[222] = temp[116]*BH->charge;
	temp[223] = temp[117]*temp[5];
	temp[224] = temp[118]*temp[37];
	temp[225] = temp[119]*temp[39];
	temp[226] = temp[120]*temp[8];
	temp[227] = temp[121]*precalculated->nmag7;
	temp[228] = temp[123]*temp[18];
	temp[229] = temp[124]*temp[14];
	temp[230] = temp[125]*precalculated->phisdot;
	temp[231] = temp[127]*temp[10];
	temp[232] = temp[128]*temp[18];
	temp[233] = temp[129]*precalculated->rdot;
	temp[234] = temp[131]*precalculated->delta_p2;
	temp[235] = 0.041666666666666664*temp[133];
	temp[236] = BH->charge*p->r4;
	temp[237] = temp[135]*temp[139];
	temp[238] = temp[137]*temp[5];
	temp[239] = temp[27]*ps->epsilon;
	temp[240] = temp[29]*temp[15];
	temp[241] = temp[30]*ps->epsilon;
	temp[242] = temp[32]*temp[142];
	temp[243] = temp[146]*temp[144];
	temp[244] = temp[19]*temp[145];
	temp[245] = temp[36]*ps->sigma12;
	temp[246] = temp[37]*temp[35];
	temp[247] = temp[38]*ps->sigma12;
	temp[248] = temp[39]*temp[35];
	temp[249] = temp[43]*precalculated->rdot;
	temp[250] = temp[5]*temp[42];
	temp[251] = temp[43]*temp[36];
	temp[252] = temp[6]*temp[37];
	temp[253] = temp[154]*temp[43];
	temp[254] = temp[155]*temp[6];
	temp[255] = temp[156]*ps->sigma6;
	temp[256] = BH->a4*temp[44];
	temp[257] = temp[158]*temp[46];
	temp[258] = temp[47]*precalculated->delta_p2;
	temp[259] = 0.25*temp[128];
	temp[260] = 0.25*temp[50];
	temp[261] = 0.08333333333333333*temp[164];
	temp[262] = temp[108]*temp[166];
	temp[263] = temp[31]*temp[54];
	temp[264] = temp[57]*temp[56];
	temp[265] = (3.*temp[59]*temp[36]*temp[61]*temp[56])/precalculated->nmag;
	temp[266] = temp[171]*temp[64];
	temp[267] = precalculated->phisdot*temp[173];
	temp[268] = 1.5*temp[9]*temp[70];
	temp[269] = 0.5*BH->a3*temp[89];
	temp[270] = temp[46]*temp[174];
	temp[271] = (BH->a7*temp[73]*temp[72]*temp[71])/precalculated->nmag;
	temp[272] = (3.*temp[9]*temp[180]*temp[71])/precalculated->nmag;
	temp[273] = temp[182]-temp[53]*temp[77]*temp[79];
	temp[274] = temp[183]*temp[79];
	temp[275] = 0.5*temp[185];
	temp[276] = temp[186]-(BH->a7*temp[83]*temp[82]*temp[67]*temp[71])/precalculated->nmag;
	temp[277] = temp[188]/precalculated->nmag;
	temp[278] = BH->a3*temp[88];
	temp[279] = temp[190]*temp[54];
	temp[280] = (2.*BH->a3*temp[73]*temp[90]*temp[71])/precalculated->nmag;
	temp[281] = precalculated->phisdot*temp[193];
	temp[282] = temp[194]*temp[189];
	temp[283] = temp[196]-(BH->a5*temp[83]*temp[95]*temp[71])/precalculated->nmag;
	temp[284] = temp[198]-temp[96]*temp[77]*temp[97]*temp[71];
	temp[285] = temp[201]-temp[98]*temp[204]*ps->sigma6*temp[99];
	temp[286] = (BH->a2*ps->epsilon*temp[102]*precalculated->sigma_p3*temp[71])/precalculated->nmag;
	temp[287] = (BH->a2*temp[205]*temp[101]*temp[106]*temp[71])/precalculated->nmag;
	temp[288] = temp[207]*temp[209];
	temp[289] = precalculated->rsdot*temp[209];
	temp[290] = temp[210]+(3.*temp[59]*p->r2*temp[110]*temp[109])/precalculated->nmag;
	temp[291] = temp[212]*temp[112];
	temp[292] = temp[213]*temp[54];
	temp[293] = temp[103]*temp[214];
	temp[294] = 0.5*temp[66]*temp[215];
	temp[295] = (BH->a6*temp[216]*temp[61]*temp[113])/precalculated->nmag;
	temp[296] = p->mass*precalculated->nmag13*(BH->a2+p->r2)*temp[219]*precalculated->delta_p2*sigma_p;
	temp[297] = temp[220]*p->r2;
	temp[298] = temp[221]*p->r4;
	temp[299] = temp[222]*p->r6;
	temp[300] = temp[223]*precalculated->delta_p2;
	temp[301] = temp[225]*precalculated->delta_p2;
	temp[302] = temp[227]*temp[122];
	temp[303] = temp[228]*ps->sigma6;
	temp[304] = temp[229]*temp[26];
	temp[305] = temp[230]*temp[126];
	temp[306] = temp[231]*temp[20];
	temp[307] = 0.5*temp[233];
	temp[308] = temp[130]*precalculated->phisdot;
	temp[309] = temp[235]*p->r2;
	temp[310] = temp[134]*temp[236];
	temp[311] = 0.16666666666666666*temp[238];
	temp[312] = temp[28]*temp[33];
	temp[313] = p->r3*temp[140];
	temp[314] = temp[141]*temp[33];
	temp[315] = temp[143]*temp[243];
	temp[316] = temp[29]*temp[146];
	temp[317] = temp[147]*temp[246];
	temp[318] = temp[149]*temp[248];
	temp[319] = temp[151]*temp[6];
	temp[320] = temp[152]*temp[251];
	temp[321] = temp[252]*temp[42];
	temp[322] = temp[254]*temp[39];
	temp[323] = 0.25*temp[256];
	temp[324] = temp[257]*temp[157];
	temp[325] = temp[160]*precalculated->phisddot;
	temp[326] = temp[260]*precalculated->phisddot;
	temp[327] = temp[107]*temp[165];
	temp[328] = temp[208]*temp[263];
	temp[329] = temp[169]+temp[265];
	temp[330] = temp[62]*temp[58];
	temp[331] = temp[266]*temp[174];
	temp[332] = BH->a5*temp[88];
	temp[333] = temp[268]*temp[175];
	temp[334] = temp[176]*temp[270];
	temp[335] = temp[271]-temp[179];
	temp[336] = temp[272]-temp[181];
	temp[337] = temp[75]*temp[76];
	temp[338] = temp[273]+temp[63]*temp[78]*temp[77]*temp[174];
	temp[339] = temp[81]*temp[79];
	temp[340] = temp[187]-(2.*BH->a5*temp[85]*temp[67]*temp[71])/precalculated->nmag;
	temp[341] = temp[191]+temp[280];
	temp[342] = temp[66]*temp[92];
	temp[343] = BH->a2*temp[93];
	temp[344] = temp[282]+temp[195]*temp[81]*temp[189];
	temp[345] = temp[85]*temp[54];
	temp[346] = temp[284]-0.5*temp[93]*temp[199];
	temp[347] = temp[285]-temp[202]*ps->sigma6*temp[100];
	temp[348] = temp[204]*temp[104];
	temp[349] = temp[66]*temp[105];
	temp[350] = temp[287]+temp[206];
	temp[351] = temp[108]*temp[288];
	temp[352] = temp[168]*precalculated->rsdot;
	temp[353] = temp[211]+(ps->epsilon*temp[62]*temp[110]*temp[109])/precalculated->nmag;
	temp[354] = temp[294]*temp[54];
	temp[355] = temp[217]+(3.*temp[60]*precalculated->phisdot2*temp[91]*temp[113])/precalculated->nmag;
	temp[356] = temp[297]*temp[1];
	temp[357] = temp[299]*temp[1];
	temp[358] = temp[224]*precalculated->delta_p2;
	temp[359] = temp[302]*precalculated->delta_p2;
	temp[360] = temp[304]*precalculated->delta_p2;
	temp[361] = temp[306]+0.5*temp[232]*temp[20];
	temp[362] = temp[132]*BH->charge;
	temp[363] = temp[310]*temp[24];
	temp[364] = temp[311]*temp[138];
	temp[365] = temp[240]*temp[313];
	temp[366] = temp[242]*temp[136];
	temp[367] = temp[25]*temp[244];
	temp[368] = temp[317]-temp[148]*temp[247]*temp[22];
	temp[369] = temp[319]*temp[250];
	temp[370] = temp[153]*temp[321];
	temp[371] = temp[322]*temp[42];
	temp[372] = temp[323]*temp[45];
	temp[373] = temp[159]*temp[258];
	temp[374] = temp[259]*temp[162];
	temp[375] = temp[261]*temp[52];
	temp[376] = temp[262]*temp[56];
	temp[377] = temp[168]*temp[264];
	temp[378] = (ps->epsilon*temp[330]*sigma_p*temp[56])/precalculated->nmag;
	temp[379] = temp[269]*temp[174];
	temp[380] = precalculated->phisdot*temp[178];
	temp[381] = temp[336]+(3.*temp[17]*temp[337]*temp[71])/precalculated->nmag;
	temp[382] = temp[276]+temp[340];
	temp[383] = temp[64]*temp[189];
	temp[384] = temp[279]*temp[71];
	temp[385] = BH->a4*temp[342];
	temp[386] = temp[87]*temp[281];
	temp[387] = temp[344]+temp[283];
	temp[388] = BH->a3*temp[345];
	temp[389] = temp[346]+(BH->a2*ps->epsilon*temp[200]*temp[71])/precalculated->nmag;
	temp[390] = 0.5*temp[349];
	temp[391] = temp[350]+temp[107]*temp[289];
	temp[392] = temp[352]*temp[112];
	temp[393] = temp[353]-temp[291]*temp[113];
	temp[394] = 1.5*temp[293];
	temp[395] = temp[354]*temp[113];
	temp[396] = temp[356]+temp[298]*temp[1];
	temp[397] = temp[358]+temp[301];
	temp[398] = temp[359]+temp[303]*precalculated->delta_p2;
	temp[399] = temp[361]+temp[307]*temp[20];
	temp[400] = temp[309]*temp[24];
	temp[401] = temp[364]-temp[239]*temp[139]*temp[312];
	temp[402] = temp[315]*temp[22];
	temp[403] = temp[245]*temp[22];
	temp[404] = temp[150]*temp[249];
	temp[405] = temp[320]*temp[41];
	temp[406] = temp[38]*temp[41];
	temp[407] = temp[255]*temp[157];
	temp[408] = temp[324]+temp[373]*sigma_p;
	temp[409] = temp[374]*sigma_p;
	temp[410] = temp[375]*sigma_p;
	temp[411] = temp[376]+temp[328]*temp[56];
	temp[412] = temp[378]+temp[331];
	temp[413] = temp[332]*temp[174];
	temp[414] = temp[379]-temp[334];
	temp[415] = temp[335]+(2.*BH->a5*temp[73]*temp[74]*temp[71])/precalculated->nmag;
	temp[416] = temp[382]-temp[277];
	temp[417] = temp[278]*temp[189];
	temp[418] = temp[192]+temp[385]*temp[71];
	temp[419] = temp[387]+temp[197];
	temp[420] = temp[389]+temp[347];
	temp[421] = 0.5*temp[103];
	temp[422] = temp[390]*temp[104];
	temp[423] = temp[351]+temp[208]*p->r4*temp[289];
	temp[424] = temp[393]-temp[55]*temp[292]*temp[113];
	temp[425] = temp[295]+temp[355];
	temp[426] = temp[396]+temp[357];
	temp[427] = temp[226]*precalculated->delta_p2;
	temp[428] = temp[399]+temp[308]*temp[234];
	temp[429] = temp[363]+temp[237]*temp[136];
	temp[430] = temp[366]+temp[402];
	temp[431] = temp[368]+temp[318];
	temp[432] = temp[369]-temp[405];
	temp[433] = temp[371]-temp[407];
	temp[434] = temp[408]-temp[325]*temp[161];
	temp[435] = temp[410]+temp[327];
	temp[436] = temp[329]+temp[170];
	temp[437] = temp[413]-temp[333]*temp[79];
	temp[438] = temp[415]-temp[381];
	temp[439] = BH->a5*temp[184];
	temp[440] = temp[275]*temp[339];
	temp[441] = temp[417]+temp[384];
	temp[442] = temp[386]-temp[343]*temp[189];
	temp[443] = temp[420]+temp[203];
	temp[444] = temp[348]*temp[100];
	temp[445] = temp[391]+temp[423];
	temp[446] = temp[290]+temp[424];
	temp[447] = temp[395]-temp[425];
	temp[448] = temp[426]+temp[300];
	temp[449] = temp[398]+temp[360];
	temp[450] = temp[428]-temp[362]*temp[24];
	temp[451] = temp[365]-temp[241]*temp[314];
	temp[452] = temp[431]-temp[404]*temp[41];
	temp[453] = temp[433]-temp[372]*temp[157];
	temp[454] = temp[435]+temp[411];
	temp[455] = temp[412]-temp[172]*temp[267];
	temp[456] = temp[438]-temp[338];
	temp[457] = temp[440]+temp[416];
	temp[458] = temp[441]-temp[341];
	temp[459] = temp[419]-(2.*temp[388]*temp[71])/precalculated->nmag;
	temp[460] = temp[445]+temp[392]*temp[109];
	temp[461] = temp[448]+temp[397];
	temp[462] = temp[305]*precalculated->delta_p2;
	temp[463] = temp[401]+temp[451];
	temp[464] = temp[316]*temp[403];
	temp[465] = temp[370]-temp[253]*temp[406];
	temp[466] = temp[326]*temp[163];
	temp[467] = temp[436]+temp[455];
	temp[468] = temp[177]*temp[380];
	temp[469] = temp[439]*temp[79];
	temp[470] = temp[458]-temp[418];
	temp[471] = temp[443]+temp[286];
	temp[472] = temp[422]*temp[99];
	temp[473] = temp[394]*temp[54];
	temp[474] = temp[461]+temp[427];
	temp[475] = temp[450]-temp[400];
	temp[476] = temp[430]-temp[367];
	temp[477] = temp[432]+temp[465];
	temp[478] = temp[409]+temp[466];
	temp[479] = temp[467]+temp[437];
	temp[480] = temp[456]+temp[274];
	temp[481] = temp[87]*temp[383];
	temp[482] = temp[459]-temp[471];
	temp[483] = temp[472]-temp[460];
	temp[484] = temp[474]-temp[449];
	temp[485] = temp[429]-temp[463];
	temp[486] = temp[452]+temp[477];
	temp[487] = temp[478]+temp[454];
	temp[488] = temp[414]-temp[468];
	temp[489] = temp[457]-temp[481];
	temp[490] = temp[482]-temp[421]*temp[444];


	phi_new = 2. * p->phi - p->phi_prev + dlambda*dlambda*((24.*(0.+temp[484]-temp[462]+temp[475]-temp[485]+temp[476]+temp[464]-temp[486]-temp[453]+temp[434]-temp[487]-temp[377]+temp[479]+temp[488]-temp[480]+temp[469]+temp[489]-temp[470]-temp[442]-temp[490]+temp[483]+temp[446]-temp[473]*temp[113]-temp[447]+temp[218]/precalculated->nmag))/temp[296]);
	return phi_new;
}

cpp_dec_float_25 calc_theta_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const Precalculated* precalculated,
	const cpp_dec_float_25 dlambda)
{
	/*
	Calculate the next value of theta.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_25 sigma_p: Precalculated sigma value.
		- const cpp_dec_float_25 delta_p: Precalculated delta value.
		- const cpp_dec_float_25 chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_25 dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_25 theta_new: New coordinate.
	*/

	cpp_dec_float_25 theta_new{};
	cpp_dec_float_25 temp[23];
	temp[0] = p->r*precalculated->rdot;
	temp[1] = precalculated->thetasdot*sigma_p;
	temp[2] = precalculated->nx*precalculated->cosphi+precalculated->ny*precalculated->sinphi;
	temp[3] = precalculated->nx*p->r;
	temp[4] = -(precalculated->nx*precalculated->cosphis)-precalculated->ny*precalculated->sinphis;
	temp[5] = precalculated->nx*ps->r;
	temp[6] = precalculated->ny*p->r*precalculated->cosphi-temp[3]*precalculated->sinphi;
	temp[7] = (ps->epsilon*ps->sigma6*(-8.*temp[0]*precalculated->thetasdot-4.*precalculated->thetasddot*sigma_p))/(p->mass*precalculated->nmag6);
	temp[8] = ps->sigma12*temp[1];
	temp[9] = (48.*precalculated->phidot*temp[6])/precalculated->nmag;
	temp[10] = 48.*precalculated->phisdot*(-(precalculated->ny*ps->r*precalculated->cosphis)+temp[5]*precalculated->sinphis);
	temp[11] = (24.*precalculated->phidot*temp[6])/precalculated->nmag;
	temp[12] = 24.*precalculated->phisdot*(-(precalculated->ny*ps->r*precalculated->cosphis)+temp[5]*precalculated->sinphis);
	temp[13] = temp[7]+(ps->epsilon*ps->sigma12*(8.*temp[0]*precalculated->thetasdot+4.*precalculated->thetasddot*sigma_p))/(p->mass*precalculated->nmag12);
	temp[14] = p->mass*precalculated->nmag13;
	temp[15] = (24.*precalculated->rdot*temp[2])/precalculated->nmag+temp[11]+(24.*precalculated->rsdot*temp[4])/precalculated->nmag+temp[12]/precalculated->nmag;
	temp[16] = -48.*precalculated->rdot*temp[2];
	temp[17] = (48.*precalculated->rsdot*temp[4])/precalculated->nmag;
	temp[18] = ps->sigma6*temp[1];
	temp[19] = ps->epsilon*temp[8]*(temp[16]/precalculated->nmag-temp[9]-temp[17]-temp[10]/precalculated->nmag);
	temp[20] = temp[13]+temp[19]/temp[14];
	temp[21] = p->mass*precalculated->nmag7;
	temp[22] = ps->epsilon*temp[18]*temp[15];


	theta_new = 2. * p->theta - p->theta_prev + dlambda*dlambda*((-0.5*(0.+temp[20]+temp[22]/temp[21]))/sigma_p);
	return theta_new;
}

std::tuple<cpp_dec_float_25, cpp_dec_float_25,
	cpp_dec_float_25, cpp_dec_float_25> eulerMoveMathematica(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda)
{
	/*
	Calculate the next value of theta.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_25 sigma_p: Precalculated sigma value.
		- const cpp_dec_float_25 delta_p: Precalculated delta value.
		- const cpp_dec_float_25 chi_p: Precalculated chi value.
		- const cpp_dec_float_25 dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_25 t_new: New t coordinate.
		- cpp_dec_float_25 r_new: New r coordinate.
		- cpp_dec_float_25 phi_new: New phi coordinate.
		- cpp_dec_float_25 theta_new: New theta coordinate.
	*/

	cpp_dec_float_25 t_new{};
	cpp_dec_float_25 r_new{};
	cpp_dec_float_25 phi_new{};
	cpp_dec_float_25 theta_new{};

	Precalculated precalculated{BH, p, ps, sigma_p, delta_p, chi_p, dlambda};

	t_new = calc_t_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda);

	r_new = calc_r_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda);

	phi_new = calc_phi_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda);

	theta_new = calc_theta_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda);

	return {t_new, r_new, phi_new, theta_new};
}
