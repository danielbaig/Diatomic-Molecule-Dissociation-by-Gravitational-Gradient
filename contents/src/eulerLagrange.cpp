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
	cpp_dec_float_25 temp[47];
	temp[0] = -2.*precalculated->chi_p2*delta_p;
	temp[1] = (0.+2.*p->r*precalculated->rdot);
	temp[2] = -2.*BH->mass+2.*p->r;
	temp[3] = -precalculated->nmag6+ps->sigma6;
	temp[4] = temp[2]*precalculated->rdot;
	temp[5] = ps->epsilon*temp[1];
	temp[6] = delta_p-BH->a2*precalculated->sintheta2;
	temp[7] = ps->epsilon*ps->sigma6;
	temp[8] = precalculated->nz*precalculated->costheta+precalculated->nx*precalculated->cosphi*precalculated->sintheta+precalculated->ny*precalculated->sinphi*precalculated->sintheta;
	temp[9] = (precalculated->thetasdot*(-(precalculated->nx*ps->r*precalculated->cosphis*precalculated->costhetas)-precalculated->ny*ps->r*precalculated->costhetas*precalculated->sinphis+precalculated->nz*ps->r*precalculated->sinthetas))/precalculated->nmag;
	temp[10] = -(precalculated->ny*ps->r*precalculated->cosphis*precalculated->sinthetas)+precalculated->nx*ps->r*precalculated->sinphis*precalculated->sinthetas;
	temp[11] = (precalculated->phidot*(precalculated->ny*p->r*precalculated->cosphi*precalculated->sintheta-precalculated->nx*p->r*precalculated->sinphi*precalculated->sintheta))/precalculated->nmag;
	temp[12] = 2.*BH->a*(BH->a2+p->r2)*precalculated->sintheta2;
	temp[13] = precalculated->tdot-precalculated->phidot*chi_p;
	temp[14] = temp[5]*ps->sigma6;
	temp[15] = 2.*BH->a*(BH->a2+p->r2);
	temp[16] = p->mass*precalculated->nmag12*precalculated->sigma_p2;
	temp[17] = (2.*chi_p*delta_p-temp[15]*precalculated->sintheta2);
	temp[18] = temp[11]+temp[9];
	temp[19] = precalculated->nx*precalculated->cosphis;
	temp[20] = precalculated->phisdot*temp[10];
	temp[21] = (precalculated->rdot*temp[8])/precalculated->nmag+temp[18]+(precalculated->rsdot*(-(precalculated->nz*precalculated->costhetas)-temp[19]*precalculated->sinthetas-precalculated->ny*precalculated->sinphis*precalculated->sinthetas))/precalculated->nmag+temp[20]/precalculated->nmag;
	temp[22] = (temp[0]/sigma_p+(2.*pow(BH->a2+p->r2,2)*precalculated->sintheta2)/sigma_p);
	temp[23] = (0.+temp[2]*precalculated->rdot);
	temp[24] = (2.*temp[4]*(-precalculated->tdot+precalculated->phidot*chi_p))/sigma_p;
	temp[25] = (4.*temp[14]*temp[3]*temp[6])/temp[16];
	temp[26] = p->mass*precalculated->nmag7*sigma_p;
	temp[27] = p->mass*precalculated->nmag13*sigma_p;
	temp[28] = (2.*temp[1]*chi_p*(-precalculated->tdot+precalculated->phidot*chi_p)*delta_p)/precalculated->sigma_p2;
	temp[29] = (2.*(BH->a2+p->r2)*temp[1]*(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot)*precalculated->sintheta2)/precalculated->sigma_p2;
	temp[30] = (precalculated->ny*p->r*precalculated->cosphi*precalculated->sintheta-precalculated->nx*p->r*precalculated->sinphi*precalculated->sintheta);
	temp[31] = p->mass*precalculated->nmag14*sigma_p;
	temp[32] = 48.*temp[7]*temp[3]*temp[17]*temp[21];
	temp[33] = ((-2.*delta_p)/sigma_p+(2.*BH->a2*precalculated->sintheta2)/sigma_p);
	temp[34] = (2.*temp[1]*temp[13]*delta_p)/precalculated->sigma_p2;
	temp[35] = temp[24]-(2.*BH->a*temp[1]*(-precalculated->phidot*(BH->a2+p->r2)+BH->a*precalculated->tdot)*precalculated->sintheta2)/precalculated->sigma_p2;
	temp[36] = temp[28]+(2.*temp[4]*chi_p*temp[13])/sigma_p;
	temp[37] = p->mass*precalculated->nmag12*sigma_p;
	temp[38] = temp[7]*(-0.5*precalculated->nmag6+ps->sigma6);
	temp[39] = BH->a*precalculated->phidot;
	temp[40] = -0.5*delta_p+0.5*BH->a2*precalculated->sintheta2;
	temp[41] = 2.*chi_p*delta_p;
	temp[42] = temp[22]*(0.+temp[34]-(4.*ps->epsilon*temp[23]*ps->sigma6*temp[3])/temp[37]+temp[35]+temp[25]+(24.*temp[7]*temp[6]*temp[21])/temp[26]+(48.*temp[7]*temp[3]*temp[6]*temp[21])/temp[27]);
	temp[43] = (96.*temp[38]*temp[30]*(precalculated->phidot*chi_p*delta_p+temp[39]*(-BH->a2-p->r2)*precalculated->sintheta2+precalculated->tdot*temp[40]))/temp[31];
	temp[44] = temp[41]/sigma_p-temp[12]/sigma_p;
	temp[45] = 0.+2.*temp[4]*chi_p;
	temp[46] = (4.*temp[14]*temp[3]*(2.*chi_p*delta_p-temp[15]*precalculated->sintheta2))/temp[16];


	t_new = 2. * p->t - p->t_prev + dlambda*dlambda*(-((-temp[42]+temp[44]*(0.+temp[36]+(4.*temp[7]*temp[3]*temp[45])/temp[37]-temp[29]-temp[46]+temp[43]-(24.*temp[7]*temp[17]*temp[21])/temp[26]-temp[32]/temp[27]))/(pow(temp[41]/sigma_p-temp[12]/sigma_p,2)-temp[33]*(temp[0]/sigma_p+(2.*pow(BH->a2+p->r2,2)*precalculated->sintheta2)/sigma_p))));
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
	cpp_dec_float_25 temp[21];
	temp[0] = precalculated->rdot*precalculated->rdot;
	temp[1] = -2.*BH->mass+2.*p->r;
	temp[2] = ps->sigma6*(-precalculated->nmag6+ps->sigma6);
	temp[3] = 2.*precalculated->phidot*chi_p*delta_p-2.*BH->a*precalculated->phidot*(BH->a2+p->r2)*precalculated->sintheta2+precalculated->tdot*(-delta_p+BH->a2*precalculated->sintheta2);
	temp[4] = precalculated->sintheta+precalculated->ny*precalculated->sinphi;
	temp[5] = (precalculated->nz*precalculated->costheta+precalculated->nx*precalculated->cosphi*temp[4]*precalculated->sintheta);
	temp[6] = -2.*p->r*temp[0];
	temp[7] = 0.+2.*p->r*precalculated->rdot;
	temp[8] = (temp[1]*pow(precalculated->tdot-precalculated->phidot*chi_p,2))/sigma_p;
	temp[9] = p->mass*precalculated->nmag12*sigma_p;
	temp[10] = (2.*p->r*pow(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot,2)*precalculated->sintheta2)/precalculated->sigma_p2;
	temp[11] = p->mass*precalculated->nmag8*sigma_p;
	temp[12] = p->mass*precalculated->nmag14*sigma_p;
	temp[13] = (2.*precalculated->rdot*temp[7])/delta_p;
	temp[14] = temp[8]-(8.*ps->epsilon*temp[1]*temp[2]*(-0.5*precalculated->tdot+precalculated->phidot*chi_p))/temp[9];
	temp[15] = p->mass*precalculated->nmag12*precalculated->sigma_p2;
	temp[16] = 48.*ps->epsilon*temp[2]*temp[5]*temp[3];
	temp[17] = temp[13]-(2.*p->r*pow(precalculated->tdot-precalculated->phidot*chi_p,2)*delta_p)/precalculated->sigma_p2;
	temp[18] = (24.*ps->epsilon*ps->sigma6*temp[5]*temp[3])/temp[11];
	temp[19] = temp[17]+temp[14];
	temp[20] = temp[10]+(8.*ps->epsilon*p->r*temp[2]*temp[3])/temp[15];


	r_new = 2. * p->r - p->r_prev + dlambda*dlambda*((-0.5*delta_p*(temp[6]/delta_p+temp[19]-(temp[1]*temp[0]*sigma_p)/precalculated->delta_p2+temp[20]+temp[18]+temp[16]/temp[12]))/sigma_p);
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
	cpp_dec_float_25 temp[700];
	temp[0] = 2.*ps->epsilon;
	temp[1] = p->r*precalculated->rdot;
	temp[2] = 2.*ps->epsilon;
	temp[3] = temp[1]*ps->sigma12;
	temp[4] = BH->mass*precalculated->nmag7;
	temp[5] = precalculated->chi_p2*precalculated->delta_p2;
	temp[6] = ps->sigma6*temp[5];
	temp[7] = precalculated->nmag*precalculated->rdot;
	temp[8] = precalculated->nmag*temp[3];
	temp[9] = precalculated->phidot*temp[1];
	temp[10] = precalculated->sintheta2-2.*BH->a2;
	temp[11] = precalculated->phidot*p->r3;
	temp[12] = precalculated->delta_p2*precalculated->sintheta2;
	temp[13] = p->r5*precalculated->rdot;
	temp[14] = ps->epsilon*precalculated->nmag7;
	temp[15] = temp[14]*p->r3;
	temp[16] = ps->sigma6*chi_p;
	temp[17] = 4.*BH->a3*ps->epsilon;
	temp[18] = BH->a*ps->epsilon;
	temp[19] = p->r3*precalculated->rdot;
	temp[20] = BH->a3*p->mass;
	temp[21] = p->mass*precalculated->nmag13;
	temp[22] = precalculated->sintheta2+2.*BH->a2;
	temp[23] = precalculated->sintheta2-BH->a2*temp[21];
	temp[24] = BH->a3*p->mass;
	temp[25] = precalculated->nmag13*precalculated->rdot;
	temp[26] = BH->a3*temp[21];
	temp[27] = chi_p*delta_p;
	temp[28] = p->mass*BH->mass;
	temp[29] = p->r2*precalculated->rdot;
	temp[30] = temp[27]*sigma_p;
	temp[31] = temp[21]*temp[19];
	temp[32] = temp[30]*precalculated->sintheta2;
	temp[33] = precalculated->chi_p2*delta_p;
	temp[34] = BH->a3*temp[21];
	temp[35] = sigma_p*precalculated->sintheta2;
	temp[36] = precalculated->phidot*temp[29];
	temp[37] = temp[21]*temp[11];
	temp[38] = sigma_p*precalculated->sintheta2;
	temp[39] = BH->a2*temp[14];
	temp[40] = temp[33]*temp[35];
	temp[41] = temp[7]*ps->sigma12;
	temp[42] = ps->epsilon*temp[8];
	temp[43] = temp[25]*precalculated->tdot;
	temp[44] = temp[33]*temp[38];
	temp[45] = precalculated->chi_p3*delta_p;
	temp[46] = temp[45]*temp[35];
	temp[47] = delta_p*precalculated->sintheta4;
	temp[48] = temp[20]*precalculated->nmag13;
	temp[49] = temp[21]*precalculated->phidot;
	temp[50] = 2.*BH->a6*temp[14];
	temp[51] = 4.*BH->a4*temp[15];
	temp[52] = ps->sigma6*temp[47];
	temp[53] = delta_p*precalculated->sintheta4;
	temp[54] = 4.*BH->a4*ps->epsilon;
	temp[55] = temp[19]*ps->sigma12;
	temp[56] = precalculated->nmag*temp[13];
	temp[57] = temp[21]*temp[9];
	temp[58] = BH->a4*temp[37];
	temp[59] = temp[27]*precalculated->sintheta4;
	temp[60] = BH->a5*temp[14];
	temp[61] = BH->a3*temp[15];
	temp[62] = BH->a5*temp[42];
	temp[63] = 4.*BH->a3*ps->epsilon;
	temp[64] = temp[27]*precalculated->sintheta4;
	temp[65] = temp[11]*precalculated->rdot;
	temp[66] = BH->a6*ps->epsilon;
	temp[67] = ps->sigma6*sigma_p;
	temp[68] = temp[50]*temp[1];
	temp[69] = temp[4]*temp[29];
	temp[70] = temp[51]*precalculated->rdot;
	temp[71] = ps->epsilon*temp[4];
	temp[72] = temp[67]*precalculated->sintheta4;
	temp[73] = 2.*BH->a6*temp[42];
	temp[74] = BH->a4*ps->epsilon;
	temp[75] = temp[29]*ps->sigma12;
	temp[76] = temp[54]*precalculated->nmag;
	temp[77] = BH->a2*ps->epsilon;
	temp[78] = precalculated->rdot*ps->sigma12;
	temp[79] = 2.*temp[77]*temp[56];
	temp[80] = temp[28]*temp[43];
	temp[81] = temp[21]*temp[1];
	temp[82] = sigma_p*precalculated->sintheta4;
	temp[83] = temp[31]*precalculated->tdot;
	temp[84] = p->r4*precalculated->rdot;
	temp[85] = precalculated->tdot*sigma_p;
	temp[86] = precalculated->phidot*precalculated->rdot;
	temp[87] = chi_p*temp[82];
	temp[88] = chi_p*sigma_p;
	temp[89] = BH->a2*temp[28];
	temp[90] = temp[49]*temp[13];
	temp[91] = precalculated->rdot*temp[16];
	temp[92] = sigma_p*precalculated->sintheta4;
	temp[93] = precalculated->nmag*temp[55];
	temp[94] = temp[81]*precalculated->tdot;
	temp[95] = precalculated->nmag13*temp[29];
	temp[96] = temp[28]*precalculated->nmag13;
	temp[97] = precalculated->chi_p2*temp[92];
	temp[98] = precalculated->rdot*ps->sigma6;
	temp[99] = 6.*ps->epsilon;
	temp[100] = precalculated->chi_p2*precalculated->delta_p3;
	temp[101] = precalculated->nx*precalculated->cosphi;
	temp[102] = ps->epsilon*temp[78];
	temp[103] = precalculated->nz*precalculated->costheta+temp[101]*precalculated->sintheta+precalculated->ny*precalculated->sinphi*precalculated->sintheta;
	temp[104] = chi_p*precalculated->delta_p2;
	temp[105] = temp[104]*sigma_p;
	temp[106] = 6.*temp[77]*precalculated->nmag5;
	temp[107] = sigma_p*precalculated->sintheta2;
	temp[108] = 6.*temp[66]*precalculated->nmag5;
	temp[109] = precalculated->sintheta4*temp[103];
	temp[110] = ps->sigma6*delta_p;
	temp[111] = 6.*temp[77]*precalculated->nmag5;
	temp[112] = precalculated->sintheta4*temp[103];
	temp[113] = sigma_p*precalculated->sintheta4;
	temp[114] = delta_p*temp[113];
	temp[115] = temp[77]*p->r4;
	temp[116] = 12.*BH->a5*ps->epsilon;
	temp[117] = delta_p*sigma_p;
	temp[118] = BH->a3*ps->epsilon;
	temp[119] = temp[29]*temp[16];
	temp[120] = temp[30]*precalculated->sintheta4;
	temp[121] = 6.*BH->a8*ps->epsilon;
	temp[122] = temp[67]*precalculated->sintheta6;
	temp[123] = temp[122]*temp[103];
	temp[124] = sigma_p*precalculated->sintheta6;
	temp[125] = temp[74]*p->r4;
	temp[126] = temp[99]*precalculated->nmag5;
	temp[127] = (precalculated->ny*p->r*precalculated->cosphi*precalculated->sintheta-precalculated->nx*p->r*precalculated->sinphi*precalculated->sintheta);
	temp[128] = precalculated->sintheta-precalculated->nx*p->r;
	temp[129] = ps->epsilon*precalculated->nmag5;
	temp[130] = temp[100]*sigma_p;
	temp[131] = temp[128]*precalculated->sinphi;
	temp[132] = ps->sigma12*temp[130];
	temp[133] = p->r*precalculated->cosphi;
	temp[134] = temp[118]*precalculated->nmag5;
	temp[135] = temp[18]*precalculated->nmag5;
	temp[136] = temp[107]*(precalculated->ny*temp[133]*temp[131]*precalculated->sintheta);
	temp[137] = precalculated->ny*temp[133]*temp[131]*precalculated->sintheta;
	temp[138] = precalculated->phidot*p->r2;
	temp[139] = temp[107]*temp[127];
	temp[140] = ps->sigma12*temp[105];
	temp[141] = temp[140]*precalculated->sintheta2;
	temp[142] = ps->sigma6*precalculated->tdot;
	temp[143] = temp[105]*precalculated->sintheta2;
	temp[144] = ps->sigma12*temp[5];
	temp[145] = temp[66]*precalculated->nmag5;
	temp[146] = temp[74]*precalculated->nmag5;
	temp[147] = p->r4*temp[110];
	temp[148] = precalculated->phidot*ps->sigma12;
	temp[149] = temp[138]*ps->sigma12;
	temp[150] = temp[77]*precalculated->phidot;
	temp[151] = temp[116]*precalculated->nmag5;
	temp[152] = p->r2*temp[142];
	temp[153] = BH->a5*ps->epsilon;
	temp[154] = temp[118]*p->r2;
	temp[155] = precalculated->tdot*temp[114];
	temp[156] = temp[129]*precalculated->phidot;
	temp[157] = temp[16]*temp[114];
	temp[158] = ps->sigma12*precalculated->tdot;
	temp[159] = precalculated->phidot*p->r4;
	temp[160] = BH->a8*ps->epsilon;
	temp[161] = ps->sigma12*temp[124];
	temp[162] = ps->sigma6*temp[85];
	temp[163] = temp[162]*precalculated->sintheta6;
	temp[164] = p->r2*ps->sigma12;
	temp[165] = ps->sigma6*precalculated->thetasdot;
	temp[166] = ps->r*precalculated->cosphis;
	temp[167] = precalculated->ny*ps->r;
	temp[168] = precalculated->sinphis+precalculated->nz*ps->r;
	temp[169] = ps->sigma12*precalculated->thetasdot;
	temp[170] = temp[167]*precalculated->costhetas;
	temp[171] = 12.*temp[134];
	temp[172] = (-(precalculated->nx*temp[166]*precalculated->costhetas)-temp[170]*temp[168]*precalculated->sinthetas);
	temp[173] = temp[170]*temp[168];
	temp[174] = temp[164]*precalculated->thetasdot;
	temp[175] = temp[5]*temp[107];
	temp[176] = -(precalculated->nx*temp[166]*precalculated->costhetas)-temp[173]*precalculated->sinthetas;
	temp[177] = temp[169]*temp[120];
	temp[178] = temp[121]*precalculated->nmag5;
	temp[179] = p->r2*temp[165];
	temp[180] = p->r4*temp[165];
	temp[181] = temp[169]*temp[124];
	temp[182] = ps->sigma6*temp[130];
	temp[183] = precalculated->nx*precalculated->cosphis;
	temp[184] = -(precalculated->nz*precalculated->costhetas)-temp[183]*precalculated->sinthetas-precalculated->ny*precalculated->sinphis*precalculated->sinthetas;
	temp[185] = temp[18]*p->r2;
	temp[186] = temp[6]*temp[107];
	temp[187] = temp[144]*temp[107];
	temp[188] = precalculated->rsdot*temp[110];
	temp[189] = temp[188]*temp[113];
	temp[190] = precalculated->rsdot*ps->sigma12;
	temp[191] = temp[190]*temp[114];
	temp[192] = precalculated->rsdot*temp[157];
	temp[193] = temp[190]*temp[120];
	temp[194] = precalculated->rsdot*temp[122];
	temp[195] = temp[194]*temp[184];
	temp[196] = precalculated->rsdot*temp[161];
	temp[197] = (-(precalculated->ny*temp[166]*precalculated->sinthetas)+precalculated->nx*ps->r*precalculated->sinphis*precalculated->sinthetas);
	temp[198] = ps->r*precalculated->sinphis;
	temp[199] = temp[16]*precalculated->delta_p2;
	temp[200] = precalculated->nx*temp[198];
	temp[201] = -(precalculated->ny*temp[166]*precalculated->sinthetas)+temp[200]*precalculated->sinthetas;
	temp[202] = temp[110]*temp[113];
	temp[203] = temp[202]*temp[201];
	temp[204] = temp[147]*temp[113];
	temp[205] = precalculated->phisdot*ps->sigma12;
	temp[206] = temp[74]*precalculated->phisdot;
	temp[207] = temp[77]*precalculated->phisdot;
	temp[208] = 12.*temp[134];
	temp[209] = temp[118]*precalculated->phisdot;
	temp[210] = temp[122]*temp[197];
	temp[211] = precalculated->phisdot*p->r2;
	temp[212] = temp[0]*precalculated->nmag7;
	temp[213] = precalculated->delta_p3-temp[2]*temp[8];
	temp[214] = temp[0]*temp[4];
	temp[215] = temp[6]*sigma_p;
	temp[216] = precalculated->nmag7*temp[1];
	temp[217] = temp[41]*temp[5];
	temp[218] = temp[8]*temp[5];
	temp[219] = BH->a4*temp[57];
	temp[220] = temp[10]*temp[37];
	temp[221] = chi_p*temp[12];
	temp[222] = precalculated->nmag13*precalculated->phidot;
	temp[223] = temp[1]*temp[16];
	temp[224] = 4.*BH->a*temp[15];
	temp[225] = precalculated->sintheta2+temp[17]*temp[8];
	temp[226] = 4.*temp[18]*temp[93];
	temp[227] = 2.*temp[48]*temp[9];
	temp[228] = 2.*BH->a*temp[37];
	temp[229] = temp[5]*temp[10];
	temp[230] = temp[1]*temp[6];
	temp[231] = temp[42]*temp[5];
	temp[232] = temp[9]*precalculated->chi_p3;
	temp[233] = temp[43]*temp[32];
	temp[234] = precalculated->tdot*temp[30];
	temp[235] = BH->a*temp[28];
	temp[236] = precalculated->tdot*temp[32];
	temp[237] = BH->a3*temp[96];
	temp[238] = temp[44]-temp[34]*temp[9];
	temp[239] = temp[36]*temp[44];
	temp[240] = precalculated->rdot*temp[44];
	temp[241] = temp[40]+4.*temp[39];
	temp[242] = temp[40]+4.*temp[77];
	temp[243] = temp[44]-4.*BH->a2;
	temp[244] = BH->a2*temp[80];
	temp[245] = sigma_p*temp[23];
	temp[246] = temp[44]-temp[89]*precalculated->nmag13;
	temp[247] = temp[46]+BH->a7*temp[57];
	temp[248] = temp[37]*precalculated->rdot;
	temp[249] = precalculated->phidot*temp[13];
	temp[250] = p->r7*precalculated->rdot;
	temp[251] = temp[52]+temp[70]*temp[52];
	temp[252] = ps->sigma6*temp[53];
	temp[253] = temp[76]*temp[55];
	temp[254] = ps->sigma12*temp[53];
	temp[255] = temp[59]-4.*temp[58];
	temp[256] = temp[59]-2.*BH->a2;
	temp[257] = 4.*temp[60]*temp[223];
	temp[258] = temp[91]*temp[53];
	temp[259] = temp[59]-temp[63]*temp[93];
	temp[260] = temp[33]*precalculated->sintheta4;
	temp[261] = precalculated->sintheta4-2.*temp[66];
	temp[262] = temp[68]*temp[67];
	temp[263] = temp[54]*temp[69];
	temp[264] = temp[67]*precalculated->sintheta4;
	temp[265] = temp[72]+2.*temp[39];
	temp[266] = 2.*temp[66]*BH->mass;
	temp[267] = temp[73]*temp[82];
	temp[268] = BH->mass*precalculated->nmag;
	temp[269] = temp[82]+2.*temp[77];
	temp[270] = temp[78]*temp[92];
	temp[271] = temp[82]+BH->a6*temp[80];
	temp[272] = temp[82]+2.*BH->a4;
	temp[273] = temp[85]*precalculated->sintheta4;
	temp[274] = BH->a2*temp[96];
	temp[275] = temp[273]-BH->a2*temp[21];
	temp[276] = BH->a6*temp[96];
	temp[277] = temp[87]+BH->a6*temp[57];
	temp[278] = 2.*BH->a4*temp[96];
	temp[279] = 2.*temp[58]*precalculated->rdot;
	temp[280] = temp[89]*temp[222];
	temp[281] = BH->a2*temp[90];
	temp[282] = 4.*BH->a5*temp[71];
	temp[283] = 4.*temp[60]*temp[223];
	temp[284] = temp[69]*temp[16];
	temp[285] = temp[61]*temp[91];
	temp[286] = temp[153]*BH->mass;
	temp[287] = 4.*temp[62]*temp[88];
	temp[288] = temp[63]*temp[268];
	temp[289] = temp[17]*temp[93];
	temp[290] = BH->a5*temp[80];
	temp[291] = BH->a5*temp[94];
	temp[292] = precalculated->sintheta4-temp[24]*BH->mass;
	temp[293] = temp[26]*temp[19];
	temp[294] = temp[87]+BH->a5*temp[96];
	temp[295] = BH->a5*temp[57];
	temp[296] = temp[82]+temp[237]*temp[36];
	temp[297] = temp[97]-2.*BH->a8;
	temp[298] = ps->sigma6*precalculated->sintheta6;
	temp[299] = precalculated->sintheta6-2.*BH->a4;
	temp[300] = precalculated->sintheta6+2.*BH->a8;
	temp[301] = 4.*temp[66]*temp[93];
	temp[302] = 2.*temp[74]*temp[56];
	temp[303] = temp[126]*temp[98];
	temp[304] = (12.*temp[102]*temp[130]*temp[103])/precalculated->nmag;
	temp[305] = temp[103]-12.*temp[135];
	temp[306] = temp[107]*temp[103];
	temp[307] = temp[102]*temp[143];
	temp[308] = (24.*temp[18]*temp[75]*temp[143]*temp[103])/precalculated->nmag;
	temp[309] = temp[78]*temp[175];
	temp[310] = temp[108]*temp[98];
	temp[311] = 12.*temp[146];
	temp[312] = temp[110]*sigma_p;
	temp[313] = temp[111]*temp[84];
	temp[314] = (12.*temp[66]*temp[78]*temp[114]*temp[103])/precalculated->nmag;
	temp[315] = (12.*temp[115]*temp[78]*temp[114]*temp[103])/precalculated->nmag;
	temp[316] = temp[117]*temp[112];
	temp[317] = temp[102]*temp[120];
	temp[318] = (24.*temp[118]*temp[75]*temp[120]*temp[103])/precalculated->nmag;
	temp[319] = temp[123]-6.*temp[146];
	temp[320] = temp[103]+(12.*BH->a8*temp[102]*temp[124]*temp[103])/precalculated->nmag;
	temp[321] = 12.*temp[125]*temp[78]*temp[124]*temp[103];
	temp[322] = chi_p*precalculated->delta_p3;
	temp[323] = ps->epsilon*temp[158];
	temp[324] = 6.*temp[156]*temp[182];
	temp[325] = ps->epsilon*precalculated->phidot;
	temp[326] = 6.*temp[134]*temp[142];
	temp[327] = temp[139]-6.*temp[135];
	temp[328] = temp[136]+(12.*temp[118]*temp[158]*precalculated->delta_p2*temp[107]*temp[137])/precalculated->nmag;
	temp[329] = precalculated->phidot*temp[199];
	temp[330] = temp[138]*temp[199];
	temp[331] = temp[118]*precalculated->phidot;
	temp[332] = (24.*temp[18]*temp[138]*temp[141]*temp[137])/precalculated->nmag;
	temp[333] = (24.*temp[77]*temp[158]*temp[143]*temp[137])/precalculated->nmag;
	temp[334] = 6.*temp[145]*precalculated->phidot;
	temp[335] = temp[138]*temp[202];
	temp[336] = precalculated->phidot*temp[204];
	temp[337] = temp[148]*temp[114];
	temp[338] = (24.*temp[74]*temp[149]*temp[114]*temp[137])/precalculated->nmag;
	temp[339] = temp[151]*temp[142];
	temp[340] = temp[208]*temp[152];
	temp[341] = (24.*temp[153]*ps->sigma12*temp[155]*temp[137])/precalculated->nmag;
	temp[342] = 12.*BH->a5*temp[156];
	temp[343] = temp[171]*temp[138];
	temp[344] = (24.*temp[153]*temp[148]*temp[120]*temp[137])/precalculated->nmag;
	temp[345] = 6.*temp[146]*temp[142];
	temp[346] = (12.*temp[74]*temp[158]*temp[120]*temp[137])/precalculated->nmag;
	temp[347] = temp[122]*temp[137];
	temp[348] = temp[122]*temp[127];
	temp[349] = temp[148]*temp[124];
	temp[350] = (24.*temp[66]*temp[149]*temp[124]*temp[137])/precalculated->nmag;
	temp[351] = 6.*BH->a7*temp[129];
	temp[352] = 6.*BH->a5*temp[129];
	temp[353] = temp[137]+(12.*BH->a7*ps->epsilon*ps->sigma12*temp[85]*precalculated->sintheta6*temp[137])/precalculated->nmag;
	temp[354] = temp[130]*temp[172];
	temp[355] = temp[169]*temp[130];
	temp[356] = temp[171]*temp[165];
	temp[357] = 12.*temp[135];
	temp[358] = temp[143]*temp[176];
	temp[359] = temp[169]*temp[143];
	temp[360] = (24.*temp[18]*temp[174]*temp[143]*temp[176])/precalculated->nmag;
	temp[361] = temp[108]*temp[165];
	temp[362] = temp[311]*temp[179];
	temp[363] = temp[111]*temp[180];
	temp[364] = (12.*temp[66]*temp[169]*temp[114]*temp[176])/precalculated->nmag;
	temp[365] = (12.*temp[115]*temp[169]*temp[114]*temp[176])/precalculated->nmag;
	temp[366] = temp[120]*temp[172];
	temp[367] = (24.*temp[154]*temp[177]*temp[176])/precalculated->nmag;
	temp[368] = temp[179]*temp[124];
	temp[369] = temp[146]*temp[180];
	temp[370] = (12.*temp[160]*temp[181]*temp[176])/precalculated->nmag;
	temp[371] = (12.*temp[125]*temp[181]*temp[176])/precalculated->nmag;
	temp[372] = (12.*ps->epsilon*precalculated->rsdot*temp[132]*temp[184])/precalculated->nmag;
	temp[373] = p->r2*precalculated->rsdot;
	temp[374] = temp[184]+(24.*temp[118]*precalculated->rsdot*temp[141]*temp[184])/precalculated->nmag;
	temp[375] = precalculated->rsdot*temp[186];
	temp[376] = precalculated->rsdot*temp[187];
	temp[377] = temp[108]*temp[189];
	temp[378] = p->r2*temp[189];
	temp[379] = temp[111]*p->r4;
	temp[380] = temp[184]-(12.*temp[66]*temp[191]*temp[184])/precalculated->nmag;
	temp[381] = (12.*temp[115]*temp[191]*temp[184])/precalculated->nmag;
	temp[382] = temp[184]-(24.*temp[153]*temp[193]*temp[184])/precalculated->nmag;
	temp[383] = temp[195]-12.*temp[145];
	temp[384] = 6.*temp[146]*p->r4;
	temp[385] = (12.*temp[160]*temp[196]*temp[184])/precalculated->nmag;
	temp[386] = (12.*temp[125]*temp[196]*temp[184])/precalculated->nmag;
	temp[387] = (12.*ps->epsilon*precalculated->phisdot*temp[132]*temp[201])/precalculated->nmag;
	temp[388] = temp[211]*temp[199];
	temp[389] = (24.*temp[209]*temp[141]*temp[201])/precalculated->nmag;
	temp[390] = temp[106]*precalculated->phisdot;
	temp[391] = temp[207]*temp[187];
	temp[392] = temp[108]*precalculated->phisdot;
	temp[393] = temp[203]+temp[111]*precalculated->phisdot;
	temp[394] = (24.*temp[206]*temp[164]*temp[114]*temp[201])/precalculated->nmag;
	temp[395] = temp[151]*precalculated->phisdot;
	temp[396] = temp[211]*temp[157];
	temp[397] = temp[153]*temp[205];
	temp[398] = (24.*temp[209]*temp[164]*temp[120]*temp[201])/precalculated->nmag;
	temp[399] = temp[210]-6.*temp[146];
	temp[400] = p->r4*temp[122];
	temp[401] = (12.*temp[160]*precalculated->phisdot*temp[161]*temp[201])/precalculated->nmag;
	temp[402] = p->r4*temp[161];
	temp[403] = precalculated->nmag13*sigma_p;
	temp[404] = precalculated->sintheta2-BH->a*p->r2;
	temp[405] = delta_p*temp[22];
	temp[406] = precalculated->sintheta2+p->r4*delta_p;
	temp[407] = temp[27]*precalculated->sintheta2;
	temp[408] = BH->a2*temp[33];
	temp[409] = temp[212]*temp[1];
	temp[410] = precalculated->chi_p2*temp[213];
	temp[411] = temp[214]*precalculated->rdot;
	temp[412] = temp[215]-temp[2]*BH->mass;
	temp[413] = temp[218]*sigma_p;
	temp[414] = temp[104]*temp[220];
	temp[415] = temp[221]-p->mass*temp[222];
	temp[416] = BH->a3*temp[14];
	temp[417] = temp[12]-temp[224]*temp[91];
	temp[418] = temp[225]*temp[104];
	temp[419] = temp[226]*temp[104];
	temp[420] = temp[227]*temp[5];
	temp[421] = temp[228]*precalculated->rdot;
	temp[422] = temp[22]*temp[231];
	temp[423] = temp[12]-temp[24]*BH->mass;
	temp[424] = temp[1]*temp[234];
	temp[425] = temp[235]*temp[95];
	temp[426] = temp[83]*temp[32];
	temp[427] = temp[86]*temp[238];
	temp[428] = temp[96]*temp[239];
	temp[429] = temp[240]-4.*BH->a2;
	temp[430] = temp[241]*temp[1];
	temp[431] = temp[242]*BH->mass;
	temp[432] = temp[42]*temp[40];
	temp[433] = temp[33]*temp[245];
	temp[434] = temp[246]*temp[86];
	temp[435] = temp[57]*temp[247];
	temp[436] = BH->a5*temp[248];
	temp[437] = 3.*temp[48]*temp[249];
	temp[438] = temp[49]*temp[250];
	temp[439] = temp[251]+2.*temp[39];
	temp[440] = temp[73]*temp[53];
	temp[441] = temp[53]-temp[79]*temp[254];
	temp[442] = precalculated->rdot*temp[256];
	temp[443] = temp[47]+4.*temp[61];
	temp[444] = temp[62]*temp[259];
	temp[445] = temp[57]*temp[260];
	temp[446] = temp[33]*temp[261];
	temp[447] = temp[72]+temp[262]*precalculated->sintheta4;
	temp[448] = 2.*BH->a2*temp[71];
	temp[449] = temp[13]*temp[72];
	temp[450] = temp[41]*temp[92];
	temp[451] = 4.*temp[74]*temp[268];
	temp[452] = temp[253]*temp[269];
	temp[453] = temp[270]-temp[79]*ps->sigma12;
	temp[454] = temp[94]*temp[272];
	temp[455] = temp[273]-2.*BH->a4;
	temp[456] = temp[274]*temp[84];
	temp[457] = temp[273]-temp[276]*temp[86];
	temp[458] = temp[278]*temp[36];
	temp[459] = temp[88]*precalculated->sintheta4;
	temp[460] = temp[91]*temp[92];
	temp[461] = temp[82]+temp[17]*temp[284];
	temp[462] = temp[92]-4.*temp[286];
	temp[463] = temp[287]*precalculated->sintheta4;
	temp[464] = temp[289]*temp[459];
	temp[465] = temp[291]*temp[88];
	temp[466] = precalculated->tdot*temp[87];
	temp[467] = temp[294]*temp[86];
	temp[468] = precalculated->chi_p2*temp[296];
	temp[469] = temp[297]*temp[14];
	temp[470] = 4.*BH->a6*temp[15];
	temp[471] = temp[14]*temp[13];
	temp[472] = temp[300]*temp[42];
	temp[473] = temp[301]*precalculated->sintheta6;
	temp[474] = temp[303]*temp[130];
	temp[475] = temp[171]*temp[91];
	temp[476] = temp[107]*temp[305];
	temp[477] = temp[306]+(24.*BH->a3*temp[307]*temp[103])/precalculated->nmag;
	temp[478] = (12.*temp[77]*temp[309]*temp[103])/precalculated->nmag;
	temp[479] = temp[312]*temp[109];
	temp[480] = temp[112]-temp[314]-(24.*temp[74]*temp[75]*temp[114]*temp[103])/precalculated->nmag;
	temp[481] = (24.*BH->a5*temp[317]*temp[103])/precalculated->nmag;
	temp[482] = 12.*temp[145];
	temp[483] = temp[319]*temp[84];
	temp[484] = (24.*temp[66]*temp[75]*temp[124]*temp[103])/precalculated->nmag;
	temp[485] = sigma_p*temp[127];
	temp[486] = temp[322]*sigma_p;
	temp[487] = temp[324]*temp[137];
	temp[488] = temp[326]*precalculated->delta_p2;
	temp[489] = temp[328]+(12.*temp[185]*temp[158]*precalculated->delta_p2*temp[107]*temp[137])/precalculated->nmag;
	temp[490] = temp[332]-12.*temp[77];
	temp[491] = temp[142]*temp[143];
	temp[492] = temp[111]*precalculated->phidot;
	temp[493] = temp[334]*temp[202];
	temp[494] = temp[146]*temp[335];
	temp[495] = temp[336]*temp[137];
	temp[496] = temp[338]+(12.*temp[150]*p->r4*ps->sigma12*temp[114]*temp[137])/precalculated->nmag;
	temp[497] = ps->sigma12*temp[155];
	temp[498] = temp[342]*temp[157];
	temp[499] = temp[157]*temp[137];
	temp[500] = temp[118]*temp[149];
	temp[501] = temp[345]*temp[120];
	temp[502] = 6.*BH->a8*temp[156];
	temp[503] = temp[145]*temp[138];
	temp[504] = temp[146]*temp[159];
	temp[505] = temp[160]*temp[349];
	temp[506] = temp[350]-(12.*temp[74]*temp[159]*temp[161]*temp[137])/precalculated->nmag;
	temp[507] = 12.*temp[153]*temp[164]*temp[85]*precalculated->sintheta6*temp[137];
	temp[508] = temp[356]*temp[143];
	temp[509] = temp[179]*temp[358];
	temp[510] = temp[360]-temp[106]*temp[165];
	temp[511] = temp[77]*temp[169];
	temp[512] = temp[361]*temp[114];
	temp[513] = temp[114]*temp[176];
	temp[514] = temp[172]-temp[364]-(24.*temp[74]*temp[174]*temp[114]*temp[176])/precalculated->nmag;
	temp[515] = (24.*temp[153]*temp[177]*temp[176])/precalculated->nmag;
	temp[516] = temp[482]*temp[368];
	temp[517] = temp[369]*temp[124];
	temp[518] = (24.*temp[66]*temp[174]*temp[124]*temp[176])/precalculated->nmag;
	temp[519] = temp[372]-temp[171]*precalculated->rsdot;
	temp[520] = temp[373]*temp[199];
	temp[521] = (24.*temp[185]*precalculated->rsdot*temp[141]*temp[184])/precalculated->nmag;
	temp[522] = temp[377]*temp[184];
	temp[523] = temp[184]+temp[379]*temp[189];
	temp[524] = p->r2*temp[191];
	temp[525] = temp[381]+temp[151]*temp[192];
	temp[526] = temp[192]*temp[382];
	temp[527] = temp[178]*temp[383];
	temp[528] = temp[384]*temp[194];
	temp[529] = (24.*temp[66]*p->r2*temp[196]*temp[184])/precalculated->nmag;
	temp[530] = temp[197]-temp[387]-temp[171]*precalculated->phisdot;
	temp[531] = temp[107]*temp[201];
	temp[532] = temp[18]*temp[211];
	temp[533] = temp[390]*temp[186];
	temp[534] = temp[392]*temp[203];
	temp[535] = temp[393]*temp[204];
	temp[536] = temp[66]*temp[205];
	temp[537] = temp[394]-(12.*temp[207]*p->r4*ps->sigma12*temp[114]*temp[201])/precalculated->nmag;
	temp[538] = temp[398]-temp[178]*precalculated->phisdot;
	temp[539] = precalculated->phisdot*temp[400];
	temp[540] = temp[66]*temp[211];
	temp[541] = 12.*temp[206]*temp[402]*temp[201];
	temp[542] = temp[403]*(chi_p*delta_p-BH->a3*temp[404]*precalculated->sintheta2);
	temp[543] = 2.*BH->a3*temp[407];
	temp[544] = temp[27]*precalculated->sintheta2;
	temp[545] = temp[409]*ps->sigma6;
	temp[546] = precalculated->delta_p3+temp[411]*temp[215];
	temp[547] = sigma_p+temp[0]*temp[413];
	temp[548] = precalculated->rdot*temp[415];
	temp[549] = temp[416]*temp[223];
	temp[550] = temp[418]*precalculated->sintheta2;
	temp[551] = temp[420]*precalculated->sintheta2;
	temp[552] = temp[230]*temp[422];
	temp[553] = temp[423]*temp[233];
	temp[554] = precalculated->sintheta2-temp[425]*temp[236];
	temp[555] = BH->a*temp[428];
	temp[556] = temp[429]*temp[71];
	temp[557] = ps->sigma6*temp[431];
	temp[558] = temp[432]+temp[244]*temp[433];
	temp[559] = temp[434]*temp[46];
	temp[560] = temp[47]+3.*temp[436];
	temp[561] = temp[47]+BH->a*temp[438];
	temp[562] = temp[439]*temp[13];
	temp[563] = temp[253]*temp[441];
	temp[564] = temp[255]*temp[442];
	temp[565] = temp[257]*temp[443];
	temp[566] = temp[444]*temp[64];
	temp[567] = temp[26]*temp[65];
	temp[568] = temp[4]*precalculated->rdot;
	temp[569] = temp[72]+temp[70]*temp[264];
	temp[570] = temp[449]+temp[266]*temp[450];
	temp[571] = temp[92]-temp[452]*temp[268];
	temp[572] = temp[92]-BH->a6*temp[454];
	temp[573] = temp[455]*temp[83];
	temp[574] = temp[275]*temp[13];
	temp[575] = temp[459]-temp[458]*temp[87];
	temp[576] = temp[84]*temp[87];
	temp[577] = temp[87]+temp[282]*temp[460];
	temp[578] = 4.*temp[285]*temp[462];
	temp[579] = temp[463]-temp[288]*temp[75];
	temp[580] = temp[87]+temp[465]*temp[292];
	temp[581] = precalculated->tdot*temp[467];
	temp[582] = temp[97]-temp[34]*temp[65];
	temp[583] = temp[470]*temp[98];
	temp[584] = ps->sigma6*temp[472];
	temp[585] = temp[473]+temp[302]*ps->sigma12;
	temp[586] = temp[474]*temp[103];
	temp[587] = precalculated->delta_p2*temp[476];
	temp[588] = temp[477]+temp[308]-temp[106]*precalculated->rdot;
	temp[589] = temp[109]+temp[311]*temp[29];
	temp[590] = temp[480]-temp[315]+temp[151]*temp[91];
	temp[591] = temp[316]-temp[481]-temp[318]-temp[178]*precalculated->rdot;
	temp[592] = temp[484]+temp[321]/precalculated->nmag;
	temp[593] = (12.*temp[323]*temp[486]*temp[137])/precalculated->nmag;
	temp[594] = temp[488]*temp[327];
	temp[595] = temp[489]+temp[208]*temp[329];
	temp[596] = temp[330]*temp[139];
	temp[597] = temp[490]*precalculated->nmag5;
	temp[598] = temp[492]*temp[6];
	temp[599] = (12.*temp[150]*temp[187]*temp[137])/precalculated->nmag;
	temp[600] = temp[495]+(12.*temp[66]*temp[337]*temp[137])/precalculated->nmag;
	temp[601] = temp[127]-temp[341]-(24.*temp[154]*temp[497]*temp[137])/precalculated->nmag;
	temp[602] = temp[501]*temp[127];
	temp[603] = temp[347]+12.*temp[503];
	temp[604] = temp[348]-(12.*temp[505]*temp[137])/precalculated->nmag;
	temp[605] = p->r2*temp[163];
	temp[606] = temp[507]/precalculated->nmag;
	temp[607] = (12.*ps->epsilon*temp[355]*temp[176])/precalculated->nmag;
	temp[608] = temp[510]*temp[175];
	temp[609] = temp[511]*temp[175];
	temp[610] = temp[512]*temp[176];
	temp[611] = temp[363]*temp[114];
	temp[612] = temp[151]*temp[165];
	temp[613] = temp[208]*temp[179];
	temp[614] = temp[367]-temp[178]*temp[165];
	temp[615] = temp[172]-6.*temp[517];
	temp[616] = temp[518]+temp[371];
	temp[617] = temp[182]*temp[184];
	temp[618] = temp[107]*temp[184];
	temp[619] = temp[107]*temp[374];
	temp[620] = temp[375]*temp[184];
	temp[621] = temp[522]+temp[311]*temp[378];
	temp[622] = temp[74]*temp[524];
	temp[623] = temp[525]*temp[184];
	temp[624] = temp[526]-(24.*temp[154]*temp[193]*temp[184])/precalculated->nmag;
	temp[625] = temp[386]+temp[126]*precalculated->phisdot;
	temp[626] = temp[197]-temp[357]*temp[388];
	temp[627] = temp[532]*temp[141];
	temp[628] = temp[533]*temp[201];
	temp[629] = temp[534]+temp[311]*temp[211];
	temp[630] = temp[536]*temp[114];
	temp[631] = temp[537]+temp[395]*temp[157];
	temp[632] = temp[197]-(24.*temp[397]*temp[120]*temp[201])/precalculated->nmag;
	temp[633] = temp[201]+temp[401]+(24.*temp[540]*temp[161]*temp[201])/precalculated->nmag;
	temp[634] = p->r2*delta_p;
	temp[635] = precalculated->sintheta2-temp[543]-2.*BH->a;
	temp[636] = temp[545]*temp[410];
	temp[637] = temp[546]-temp[2]*temp[216];
	temp[638] = temp[547]-temp[219]*temp[414];
	temp[639] = 4.*temp[549]*temp[417];
	temp[640] = temp[550]+temp[419]*precalculated->sintheta2;
	temp[641] = temp[552]*temp[23];
	temp[642] = temp[26]*temp[424];
	temp[643] = temp[426]+temp[237]*temp[427];
	temp[644] = temp[37]*temp[556];
	temp[645] = temp[557]*temp[41];
	temp[646] = temp[1]*precalculated->tdot;
	temp[647] = temp[560]*temp[47];
	temp[648] = temp[47]+temp[68]*temp[562];
	temp[649] = temp[563]-2.*BH->a6;
	temp[650] = temp[90]*temp[64];
	temp[651] = temp[258]-4.*temp[566];
	temp[652] = temp[567]*temp[446];
	temp[653] = temp[263]*temp[569];
	temp[654] = temp[265]*temp[570];
	temp[655] = temp[75]*temp[571];
	temp[656] = temp[271]*temp[572];
	temp[657] = temp[573]*temp[82];
	temp[658] = temp[457]*temp[277];
	temp[659] = temp[459]-temp[280]*temp[576];
	temp[660] = temp[461]*temp[92];
	temp[661] = temp[87]+temp[579]*temp[87];
	temp[662] = temp[580]*temp[95];
	temp[663] = temp[581]*temp[97];
	temp[664] = temp[582]*temp[469];
	temp[665] = temp[583]*temp[299];
	temp[666] = precalculated->sintheta6+temp[585]*precalculated->sintheta6;
	temp[667] = temp[588]*temp[186];
	temp[668] = temp[310]*temp[117];
	temp[669] = temp[313]*temp[312];
	temp[670] = temp[109]+temp[208]*temp[119];
	temp[671] = temp[29]*temp[483];
	temp[672] = temp[592]+temp[126]*temp[142];
	temp[673] = temp[487]+(12.*temp[325]*temp[132]*temp[137])/precalculated->nmag;
	temp[674] = 12.*temp[135];
	temp[675] = (24.*temp[331]*temp[141]*temp[137])/precalculated->nmag;
	temp[676] = temp[139]-temp[599]-temp[493]*temp[127];
	temp[677] = temp[600]+temp[496]+temp[339]*temp[114];
	temp[678] = temp[498]*temp[127];
	temp[679] = temp[344]+(24.*temp[500]*temp[120]*temp[137])/precalculated->nmag;
	temp[680] = temp[604]-temp[506]-temp[351]*temp[163];
	temp[681] = temp[606]+temp[126]*temp[165];
	temp[682] = temp[172]-temp[357]*temp[509];
	temp[683] = temp[608]*temp[176];
	temp[684] = temp[610]+temp[362]*temp[513];
	temp[685] = temp[612]*temp[120];
	temp[686] = temp[366]-temp[515]-temp[614]*temp[124];
	temp[687] = temp[370]+temp[616];
	temp[688] = temp[617]-temp[519]*temp[199];
	temp[689] = temp[619]+temp[521]-temp[106]*temp[620];
	temp[690] = temp[621]*temp[523];
	temp[691] = temp[623]+temp[208]*p->r2;
	temp[692] = p->r2*temp[195];
	temp[693] = temp[184]+temp[385]+temp[529]+temp[625]*temp[182];
	temp[694] = temp[389]+(24.*temp[627]*temp[201])/precalculated->nmag;
	temp[695] = temp[629]*temp[535];
	temp[696] = temp[631]*temp[201];
	temp[697] = temp[632]-temp[538]*temp[210];
	temp[698] = temp[539]*temp[633];
	temp[699] = p->mass*temp[542]*(BH->a4*temp[405]*temp[634]*temp[406]*temp[635]*p->r2*temp[544]+temp[408]*precalculated->sintheta2);


	phi_new = 2. * p->phi - p->phi_prev + dlambda*dlambda*((-2.*(0.+temp[636]*precalculated->chi_p2*temp[637]*temp[412]*temp[217]*temp[638]*temp[548]*temp[13]*temp[221]-temp[639]*precalculated->delta_p2*temp[640]+temp[551]+temp[421]*temp[229]*temp[14]*temp[641]*temp[232]*temp[553]+temp[642]*temp[554]+BH->a*temp[643]*temp[40]+temp[555]-BH->a*temp[644]*temp[98]*temp[430]*temp[645]*temp[243]*temp[558]*temp[646]*temp[559]+BH->a2*temp[435]*temp[647]+temp[437]*temp[561]*temp[648]*temp[252]-temp[440]-temp[649]*temp[57]*temp[564]*temp[650]+temp[565]*temp[651]+BH->a5*temp[445]+temp[652]*temp[568]*temp[447]-temp[653]-temp[448]*temp[84]*temp[654]-temp[267]+temp[451]*temp[655]*p->r4*temp[453]*temp[656]*temp[28]*temp[95]*temp[657]+temp[456]*temp[574]*temp[658]*temp[575]+temp[279]*temp[659]+temp[281]*temp[577]-temp[283]*temp[660]-temp[578]*temp[41]*temp[661]+temp[464]-temp[290]*temp[662]*temp[466]+temp[293]*temp[663]-temp[295]*temp[468]*temp[664]*temp[1]*temp[298]-temp[665]*temp[471]*temp[584]*temp[666]+temp[586]-temp[304]-temp[475]*temp[587]*temp[119]*precalculated->delta_p2*temp[667]*temp[103]+temp[478]+temp[668]*temp[589]*temp[479]+temp[669]*temp[590]*temp[117]*temp[670]*temp[591]*temp[123]-temp[482]*temp[671]*temp[122]*temp[320]+temp[672]*temp[322]*temp[485]-temp[593]-temp[673]-temp[594]*temp[152]*precalculated->delta_p2*temp[595]*temp[136]+temp[674]*temp[596]-temp[675]-temp[597]*temp[491]*temp[137]+temp[333]+temp[598]*temp[676]-12.*temp[494]*temp[127]-temp[106]*temp[677]*temp[137]+temp[340]*temp[114]*temp[601]-temp[678]-temp[343]*temp[499]+temp[679]+temp[602]-temp[346]+temp[502]*temp[603]*temp[347]+6.*temp[504]*temp[680]*temp[127]-temp[352]*temp[605]*temp[353]+temp[681]*temp[354]-temp[607]-temp[508]*temp[682]+(24.*temp[118]*temp[359]*temp[176])/precalculated->nmag+temp[683]+(12.*temp[609]*temp[176])/precalculated->nmag+temp[684]+temp[611]*temp[514]-temp[365]+temp[685]*temp[176]+temp[613]*temp[686]*temp[172]-temp[516]*temp[615]*temp[176]+temp[687]+temp[126]*precalculated->rsdot*temp[688]*temp[618]-temp[357]*temp[520]*temp[689]+(12.*temp[77]*temp[376]*temp[184])/precalculated->nmag+temp[690]*temp[380]-(24.*temp[622]*temp[184])/precalculated->nmag-temp[691]*temp[624]-temp[527]*temp[692]-temp[528]*temp[693]*temp[530]*temp[199]*temp[107]*temp[626]*temp[531]+temp[694]-temp[628]+(12.*temp[391]*temp[201])/precalculated->nmag+temp[695]*temp[197]-(12.*temp[630]*temp[201])/precalculated->nmag-temp[696]+temp[208]*temp[396]*temp[697]-temp[482]*temp[211]*temp[399]*temp[698]+temp[541]/precalculated->nmag))/temp[699]);
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
	cpp_dec_float_25 temp[24];
	temp[0] = (BH->l+BH->a*precalculated->costheta);
	temp[1] = ps->epsilon*ps->sigma6;
	temp[2] = 2.*precalculated->phidot*chi_p*delta_p-2.*BH->a*precalculated->phidot*(BH->a2+p->r2)*precalculated->sintheta2+precalculated->tdot*(-delta_p+BH->a2*precalculated->sintheta2);
	temp[3] = precalculated->sinphi-precalculated->nz*p->r;
	temp[4] = temp[1]*(-precalculated->nmag6+ps->sigma6);
	temp[5] = p->r*precalculated->costheta;
	temp[6] = 2.*BH->l*precalculated->sintheta+BH->a*precalculated->sin2theta;
	temp[7] = (2.*(precalculated->rdot*precalculated->rdot)*temp[0]*precalculated->sintheta)/delta_p;
	temp[8] = (2.*pow(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot,2)*temp[0]*precalculated->sintheta3)/precalculated->sigma_p2;
	temp[9] = p->r*precalculated->cosphi;
	temp[10] = precalculated->ny*temp[5];
	temp[11] = p->mass*precalculated->nmag8*sigma_p;
	temp[12] = temp[9]*precalculated->costheta;
	temp[13] = p->mass*precalculated->nmag14*sigma_p;
	temp[14] = (2.*precalculated->phidot*(precalculated->tdot-precalculated->phidot*chi_p)*delta_p*temp[6])/sigma_p;
	temp[15] = (BH->a2+p->r2)*precalculated->sin2theta;
	temp[16] = temp[7]+(2.*pow(precalculated->tdot-precalculated->phidot*chi_p,2)*temp[0]*delta_p*precalculated->sintheta)/precalculated->sigma_p2;
	temp[17] = (precalculated->nx*temp[12]+temp[10]*temp[3]*precalculated->sintheta);
	temp[18] = 4.*temp[4]*(BH->a2*precalculated->tdot*precalculated->sin2theta-2.*precalculated->phidot*(BH->a*temp[15]-delta_p*temp[6]));
	temp[19] = temp[16]-temp[8];
	temp[20] = temp[0]*precalculated->sintheta;
	temp[21] = (24.*temp[1]*temp[17]*temp[2])/temp[11];
	temp[22] = (pow(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot,2)*precalculated->sin2theta)/sigma_p;
	temp[23] = temp[19]-(8.*temp[4]*temp[20]*temp[2])/(p->mass*precalculated->nmag12*precalculated->sigma_p2);


	theta_new = 2. * p->theta - p->theta_prev + dlambda*dlambda*((-0.5*(0.+temp[23]+temp[21]+(48.*temp[4]*temp[17]*temp[2])/temp[13]-temp[22]-temp[14]-temp[18]/(p->mass*precalculated->nmag12*sigma_p)))/sigma_p);
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
