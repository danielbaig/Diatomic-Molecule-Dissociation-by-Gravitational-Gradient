// Configured==true
#include "eulerLagrange.h"

const cpp_dec_float_n epsilon_0{ 1. / (4.*M_PI) }; // Geometrised-Gaussian units


Precalculated::Precalculated(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const cpp_dec_float_n dlambda):
	sigma_p2{ sigma_p*sigma_p },
	sigma_p3{ sigma_p*sigma_p2 },
	sigma_p4{ sigma_p2*sigma_p2 },
	delta_p2{ delta_p*delta_p },
	delta_p3{ delta_p*delta_p2 },
	delta_p4{ delta_p2*delta_p2 },
	chi_p2{ chi_p*chi_p },
	chi_p3{ chi_p*chi_p2 },
	chi_p4{ chi_p2*chi_p2 },

	tdot{ (p->t_prev - p->t_prevprev) / (dlambda) },
	rdot{ (p->r_prev - p->r_prevprev) / (dlambda) },
	phidot{ (p->phi_prev - p->phi_prevprev) / (dlambda) },
	thetadot{ (p->theta_prev - p->theta_prevprev) / (dlambda) },

	tsdot{ (ps->t_prev - ps->t_prevprev) / (dlambda) },
	rsdot{ (ps->r_prev - ps->r_prevprev) / (dlambda) },
	phisdot{ (ps->phi_prev - ps->phi_prevprev) / (dlambda) },
	thetasdot{ (ps->theta_prev - ps->theta_prevprev) / (dlambda) },
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


Precalculated::~Precalculated()
{
	/*
	Delete the struct.
	*/
delete[] temp0;
temp0 = nullptr;
delete[] temp1;
temp1 = nullptr;
delete[] temp2;
temp2 = nullptr;
}


cpp_dec_float_n calc_t_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda, cpp_dec_float_n* temp)
{
	/*
	Calculate the next value of t.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_n sigma_p: Precalculated sigma value.
		- const cpp_dec_float_n delta_p: Precalculated delta value.
		- const cpp_dec_float_n chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_n dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_n t_new: New coordinate.
	*/

	cpp_dec_float_n t_new{};
	temp[0] = (2.*BH->a2*(BH->a2+p->r2))/sigma_p+(2.*p->r2*(BH->a2+p->r2))/sigma_p-(2.*precalculated->chi_p2*delta_p)/sigma_p; // (...)
	temp[1] = -1.*precalculated->phidot*(BH->a2+p->r2)+BH->a*precalculated->tdot; // (...)
	temp[2] = (8.*ps->epsilon*p->r*precalculated->rdot*ps->sigma6*(-1.*precalculated->nmag6+ps->sigma6)*(-1.*BH->a2+delta_p))/(p->mass*precalculated->nmag12*precalculated->sigma_p2); // +...-
	temp[3] = precalculated->rdot*ps->sigma6; // *...*
	temp[4] = (2.*(-2.*BH->mass+2.*p->r)*precalculated->rdot*(-1.*precalculated->tdot+precalculated->phidot*chi_p))/sigma_p; // +...+
	temp[5] = (precalculated->phidot*(precalculated->ny*p->r*precalculated->cosphi-precalculated->nx*p->r*precalculated->sinphi))/precalculated->nmag; // +...+
	temp[6] = precalculated->phisdot*(-(precalculated->ny*ps->r*precalculated->cosphis)+precalculated->nx*ps->r*precalculated->sinphis); // (...)
	temp[7] = -1.*BH->a2+delta_p; // (...)
	temp[8] = temp[5]+(precalculated->rsdot*(-(precalculated->nx*precalculated->cosphis)-precalculated->ny*precalculated->sinphis))/precalculated->nmag; // +...+
	temp[9] = -(4.*p->charge*BH->charge*p->r2*precalculated->rdot*chi_p)/(p->mass*precalculated->sigma_p2); // -...+
	temp[10] = (8.*ps->epsilon*p->r*temp[3]*(-1.*precalculated->nmag6+ps->sigma6)*(BH->a3+BH->a*p->r2-chi_p*delta_p))/(p->mass*precalculated->nmag12*precalculated->sigma_p2); // +...+
	temp[11] = (2.*(-2.*BH->mass+2.*p->r)*precalculated->rdot*chi_p*(precalculated->tdot-precalculated->phidot*chi_p))/sigma_p; // +...-
	temp[12] = ps->epsilon*ps->sigma6; // *...*
	temp[13] = BH->a3+BH->a*p->r2-chi_p*delta_p; // (...)
	temp[14] = -precalculated->nx*p->r; // -...*
	temp[15] = 24.*temp[12]*temp[13]*((precalculated->rdot*(precalculated->nx*precalculated->cosphi+precalculated->ny*precalculated->sinphi))/precalculated->nmag+temp[8]+temp[6]/precalculated->nmag); // (...)
	temp[16] = -2.*BH->a*(BH->a2+p->r2); // (...)
	temp[17] = (2.*BH->a2)/sigma_p-(2.*delta_p)/sigma_p; // (...)
	temp[18] = 4.*p->charge*BH->charge*p->r2*precalculated->rdot; // (...)
	temp[19] = p->r*precalculated->rdot; // *...*
	temp[20] = (4.*temp[19]*(precalculated->tdot-precalculated->phidot*chi_p)*delta_p)/precalculated->sigma_p2; // +...+
	temp[21] = -(4.*ps->epsilon*(-2.*BH->mass+2.*p->r)*temp[3]*(-1.*precalculated->nmag6+ps->sigma6))/(p->mass*precalculated->nmag12*sigma_p); // -...+
	temp[22] = -1.*precalculated->nmag6+ps->sigma6; // (...)
	temp[23] = precalculated->nx*precalculated->cosphi+precalculated->ny*precalculated->sinphi; // (...)
	temp[24] = -4.*p->r*(BH->a2+p->r2)*precalculated->rdot*(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot); // (...)
	temp[25] = temp[10]+(4.*precalculated->phidot*p->r*(BH->a2+p->r2)*precalculated->rdot)/sigma_p; // +...-
	temp[26] = temp[11]-(8.*ps->epsilon*temp[3]*temp[22]*(BH->a*p->r-0.5*(-2.*BH->mass+2.*p->r)*chi_p))/(p->mass*precalculated->nmag12*sigma_p); // +...-
	temp[27] = p->r*precalculated->cosphi; // *...+
	temp[28] = temp[15]/(p->mass*precalculated->nmag7*sigma_p); // +...+
	temp[29] = 2.*chi_p*delta_p; // (...)
	temp[30] = -(4.*BH->a*temp[19]*temp[1])/precalculated->sigma_p2; // -...+
	temp[31] = -(4.*BH->a*precalculated->phidot*p->r*precalculated->rdot)/sigma_p; // -...+
	temp[32] = 48.*temp[12]*temp[22]*temp[7]*((precalculated->rdot*temp[23])/precalculated->nmag+temp[8]+temp[6]/precalculated->nmag); // (...)
	temp[33] = (4.*temp[19]*chi_p*(-1.*precalculated->tdot+precalculated->phidot*chi_p)*delta_p)/precalculated->sigma_p2; // +...+
	temp[34] = 48.*temp[12]*(-0.5*precalculated->nmag6+ps->sigma6)*(precalculated->tdot*temp[7]+precalculated->phidot*temp[13])*(precalculated->ny*temp[27]+temp[14]*precalculated->sinphi); // (...)
	temp[35] = temp[30]+temp[20]; // +...+
	temp[36] = p->charge*BH->charge; // *...*
	temp[37] = temp[31]+temp[21]; // +...+
	temp[38] = temp[12]*temp[7]; // *...*
	temp[39] = temp[16]/sigma_p+temp[29]/sigma_p; // (...)
	temp[40] = temp[9]+temp[33]; // +...+
	temp[41] = 4.*temp[19]*temp[1]; // (...)
	temp[42] = temp[26]-temp[34]/(p->mass*precalculated->nmag14*sigma_p); // +...+
	temp[43] = temp[35]+temp[2]; // +...-
	temp[44] = temp[37]+temp[4]; // +...+
	temp[45] = (precalculated->rdot*temp[23])/precalculated->nmag+temp[8]+temp[6]/precalculated->nmag; // (...)
	temp[46] = temp[40]+temp[25]; // +...-
	temp[47] = (2.*temp[36]*precalculated->rdot*chi_p)/(p->mass*sigma_p); // +...+
	temp[48] = temp[43]-(2.*temp[36]*precalculated->rdot)/(p->mass*sigma_p); // +...+
	temp[49] = temp[24]/precalculated->sigma_p2+temp[46]-temp[41]/sigma_p+temp[47]+temp[42]+temp[28]+(48.*temp[12]*temp[22]*temp[13]*temp[45])/(p->mass*precalculated->nmag13*sigma_p); // (...)
	temp[50] = temp[48]+temp[44]; // +...+
	temp[51] = temp[50]+(24.*temp[38]*temp[45])/(p->mass*precalculated->nmag7*sigma_p); // +...+


	t_new = 2. * p->t - p->t_prev + dlambda*dlambda*(-((-1.*temp[0]*(temp[18]/(p->mass*precalculated->sigma_p2)+temp[51]+temp[32]/(p->mass*precalculated->nmag13*sigma_p))+temp[39]*temp[49])/(pow(temp[16]/sigma_p+temp[29]/sigma_p,2)-temp[17]*temp[0])));
	return t_new;
}

cpp_dec_float_n calc_r_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda, cpp_dec_float_n* temp)
{
	/*
	Calculate the next value of r.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_n sigma_p: Precalculated sigma value.
		- const cpp_dec_float_n delta_p: Precalculated delta value.
		- const cpp_dec_float_n chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_n dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_n r_new: New coordinate.
	*/

	cpp_dec_float_n r_new{};
	temp[0] = 2.*p->r*(precalculated->rdot*precalculated->rdot); // (...)
	temp[1] = (4.*p->charge*BH->charge*p->r2*(-1.*precalculated->tdot+precalculated->phidot*chi_p))/(p->mass*precalculated->sigma_p2); // +...-
	temp[2] = -1.*precalculated->nmag6+ps->sigma6; // (...)
	temp[3] = BH->a3+BH->a*p->r2-chi_p*delta_p; // (...)
	temp[4] = -(4.*precalculated->phidot*p->r*(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot))/sigma_p; // -...+
	temp[5] = -1.*precalculated->tdot+precalculated->phidot*chi_p; // (...)
	temp[6] = 2.*BH->a*precalculated->phidot*p->r+(-2.*BH->mass+2.*p->r)*(precalculated->tdot-precalculated->phidot*chi_p); // (...)
	temp[7] = 24.*ps->epsilon*ps->sigma6*(precalculated->tdot*(-1.*BH->a2+delta_p)+precalculated->phidot*temp[3])*(precalculated->nx*precalculated->cosphi+precalculated->ny*precalculated->sinphi); // (...)
	temp[8] = (2.*p->r*pow(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot,2))/precalculated->sigma_p2; // +...+
	temp[9] = -(8.*ps->epsilon*p->r*ps->sigma6*temp[2]*(precalculated->tdot*(-1.*BH->a2+delta_p)+precalculated->phidot*temp[3]))/(p->mass*precalculated->nmag12*precalculated->sigma_p2); // -...+
	temp[10] = ps->sigma6*temp[2]; // *...*
	temp[11] = -((-2.*BH->mass+2.*p->r)*(precalculated->rdot*precalculated->rdot)*sigma_p)/precalculated->delta_p2; // -...-
	temp[12] = -1.*BH->a2+delta_p; // (...)
	temp[13] = precalculated->nx*precalculated->cosphi+precalculated->ny*precalculated->sinphi; // (...)
	temp[14] = temp[8]+temp[1]; // +...-
	temp[15] = temp[9]+temp[4]; // +...+
	temp[16] = -(2.*p->charge*BH->charge*temp[5])/(p->mass*sigma_p); // -...+
	temp[17] = temp[11]-temp[7]/(p->mass*precalculated->nmag8*sigma_p); // +...-
	temp[18] = temp[14]-(2.*p->r*pow(precalculated->tdot-precalculated->phidot*chi_p,2)*delta_p)/precalculated->sigma_p2; // +...+
	temp[19] = (4.*ps->epsilon*temp[10]*temp[6])/(p->mass*precalculated->nmag12*sigma_p); // +...+
	temp[20] = temp[18]+temp[15]; // +...+
	temp[21] = temp[16]+temp[19]; // +...+
	temp[22] = ps->epsilon*temp[10]; // *...*
	temp[23] = temp[20]+((-2.*BH->mass+2.*p->r)*pow(precalculated->tdot-precalculated->phidot*chi_p,2))/sigma_p; // +...+
	temp[24] = temp[23]+temp[21]; // +...+
	temp[25] = precalculated->tdot*temp[12]+precalculated->phidot*temp[3]; // (...)
	temp[26] = temp[24]+temp[17]; // +...-
	temp[27] = 48.*temp[22]*temp[25]*temp[13]; // (...)


	r_new = 2. * p->r - p->r_prev + dlambda*dlambda*((-0.5*delta_p*(temp[0]/delta_p+temp[26]-temp[27]/(p->mass*precalculated->nmag14*sigma_p)))/sigma_p);
	return r_new;
}

cpp_dec_float_n calc_phi_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda, cpp_dec_float_n* temp)
{
	/*
	Calculate the next value of phi.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_n sigma_p: Precalculated sigma value.
		- const cpp_dec_float_n delta_p: Precalculated delta value.
		- const cpp_dec_float_n chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_n dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_n phi_new: New coordinate.
	*/

	cpp_dec_float_n phi_new{};
	temp[0] = -1.*BH->a6*precalculated->nmag13; // -...*
	temp[1] = -2.*BH->a4*precalculated->nmag13; // -...*
	temp[2] = -BH->a2*precalculated->nmag13; // -...*
	temp[3] = BH->a5*precalculated->nmag13; // +...*
	temp[4] = precalculated->rdot*chi_p; // *...+
	temp[5] = p->charge*BH->charge; // *...*
	temp[6] = precalculated->phidot*p->r; // *...*
	temp[7] = -3.*BH->a5*p->mass; // -...*
	temp[8] = precalculated->phidot*p->r3; // *...*
	temp[9] = -3.*BH->a3*p->mass; // -...*
	temp[10] = precalculated->phidot*p->r5; // *...*
	temp[11] = -BH->a*p->mass; // -...*
	temp[12] = precalculated->phidot*p->r7; // *...*
	temp[13] = 2.*BH->a6*p->mass; // +...*
	temp[14] = temp[6]*precalculated->rdot; // *...*
	temp[15] = temp[5]*p->r2; // *...*
	temp[16] = chi_p*delta_p; // *...+
	temp[17] = precalculated->nmag13*temp[8]; // *...*
	temp[18] = temp[5]*p->r4; // *...*
	temp[19] = 2.*BH->a2*p->mass; // +...*
	temp[20] = temp[10]*precalculated->rdot; // *...*
	temp[21] = temp[14]*precalculated->chi_p2; // *...*
	temp[22] = precalculated->rdot*precalculated->chi_p2; // *...*
	temp[23] = BH->a4*p->mass; // +...*
	temp[24] = temp[14]*chi_p; // *...*
	temp[25] = temp[19]*temp[17]; // +...*
	temp[26] = chi_p*precalculated->delta_p2; // *...+
	temp[27] = temp[20]*chi_p; // *...*
	temp[28] = -2.*BH->a3*p->mass; // -...*
	temp[29] = temp[21]*precalculated->delta_p2; // *...-
	temp[30] = temp[22]*precalculated->delta_p2; // *...+
	temp[31] = temp[14]*precalculated->chi_p3; // *...*
	temp[32] = 0.5*BH->a6*precalculated->nmag13; // +...*
	temp[33] = -BH->a7*p->mass; // -...*
	temp[34] = temp[14]*sigma_p; // *...+
	temp[35] = temp[15]*precalculated->rdot; // *...*
	temp[36] = 0.5*BH->a2*precalculated->nmag13; // +...*
	temp[37] = -BH->a3*p->mass; // -...*
	temp[38] = temp[20]*sigma_p; // *...+
	temp[39] = ps->epsilon*BH->mass; // *...*
	temp[40] = precalculated->rdot*ps->sigma6; // *...*
	temp[41] = precalculated->nmag7*p->r; // *...*
	temp[42] = BH->a4*temp[39]; // *...*
	temp[43] = p->r2*temp[40]; // *...*
	temp[44] = -6.*BH->a4*ps->epsilon; // -...*
	temp[45] = p->r3*temp[40]; // *...*
	temp[46] = 2.*BH->a2*temp[39]; // +...*
	temp[47] = p->r4*temp[40]; // *...*
	temp[48] = -2.*BH->a2*ps->epsilon; // -...*
	temp[49] = p->r5*temp[40]; // *...*
	temp[50] = -2.*BH->a6*temp[39]; // -...*
	temp[51] = precalculated->rdot*ps->sigma12; // *...*
	temp[52] = precalculated->nmag*p->r; // *...*
	temp[53] = temp[42]*precalculated->nmag; // *...*
	temp[54] = 6.*BH->a4*ps->epsilon; // +...*
	temp[55] = p->r3*temp[51]; // *...*
	temp[56] = -2.*BH->a2*temp[39]; // -...*
	temp[57] = p->r4*temp[51]; // *...*
	temp[58] = 2.*BH->a2*ps->epsilon; // +...*
	temp[59] = p->r5*temp[51]; // *...*
	temp[60] = -BH->a6*p->mass; // -...*
	temp[61] = precalculated->nmag13*precalculated->rdot; // *...*
	temp[62] = precalculated->nmag13*p->r; // *...*
	temp[63] = precalculated->tdot*sigma_p; // *...-
	temp[64] = precalculated->nmag13*p->r2; // *...*
	temp[65] = 3.*BH->a4*p->mass; // +...*
	temp[66] = p->r3*precalculated->rdot; // *...*
	temp[67] = BH->mass*precalculated->nmag13; // *...*
	temp[68] = precalculated->tdot*sigma_p; // *...+
	temp[69] = BH->a6*p->mass; // +...*
	temp[70] = precalculated->phidot*precalculated->rdot; // *...*
	temp[71] = temp[5]*precalculated->rdot; // *...*
	temp[72] = temp[24]*sigma_p; // *...+
	temp[73] = p->mass*temp[67]; // *...*
	temp[74] = p->r2*precalculated->rdot; // *...*
	temp[75] = precalculated->nmag13*temp[35]; // *...*
	temp[76] = BH->a2*temp[73]; // +...*
	temp[77] = p->r4*precalculated->rdot; // *...*
	temp[78] = temp[27]*sigma_p; // *...-
	temp[79] = temp[39]*precalculated->nmag7; // *...*
	temp[80] = temp[41]*temp[40]; // *...*
	temp[81] = temp[79]*temp[43]; // *...*
	temp[82] = 2.*BH->a3*ps->epsilon; // +...*
	temp[83] = temp[45]*chi_p; // *...*
	temp[84] = 2.*BH->a5*temp[39]; // +...*
	temp[85] = temp[51]*chi_p; // *...*
	temp[86] = -2.*BH->a5*ps->epsilon; // -...*
	temp[87] = 2.*BH->a3*temp[39]; // +...*
	temp[88] = p->r2*temp[85]; // *...*
	temp[89] = -2.*BH->a3*ps->epsilon; // -...*
	temp[90] = temp[55]*chi_p; // *...*
	temp[91] = BH->a5*p->mass; // +...*
	temp[92] = temp[61]*precalculated->tdot; // *...*
	temp[93] = temp[62]*precalculated->rdot; // *...*
	temp[94] = BH->mass*temp[64]; // *...*
	temp[95] = precalculated->tdot*chi_p; // *...*
	temp[96] = temp[66]*temp[95]; // *...*
	temp[97] = temp[73]*precalculated->phidot; // *...*
	temp[98] = temp[21]*sigma_p; // *...+
	temp[99] = precalculated->phidot*p->r2; // *...*
	temp[100] = p->mass*temp[17]; // *...*
	temp[101] = temp[14]*delta_p; // *...*
	temp[102] = BH->a3*temp[100]; // *...*
	temp[103] = delta_p*sigma_p; // *...+
	temp[104] = p->mass*precalculated->nmag13; // *...*
	temp[105] = ps->epsilon*temp[52]; // *...*
	temp[106] = precalculated->nmag*temp[55]; // *...*
	temp[107] = precalculated->tdot*delta_p; // *...*
	temp[108] = temp[66]*temp[107]; // *...*
	temp[109] = BH->a3*precalculated->nmag13; // *...*
	temp[110] = -0.5*BH->a*temp[75]; // -...*
	temp[111] = -2.*BH->a3*temp[79]; // -...*
	temp[112] = 4.*BH->a3*ps->epsilon; // +...*
	temp[113] = delta_p*sigma_p; // *...-
	temp[114] = 2.*BH->a*ps->epsilon; // +...*
	temp[115] = temp[83]*temp[103]; // *...+
	temp[116] = temp[85]*temp[113]; // *...-
	temp[117] = temp[85]*temp[103]; // *...+
	temp[118] = precalculated->nmag*temp[88]; // *...*
	temp[119] = temp[90]*temp[103]; // *...+
	temp[120] = BH->mass*temp[92]; // *...*
	temp[121] = temp[93]*temp[95]; // *...*
	temp[122] = temp[94]*precalculated->rdot; // *...*
	temp[123] = temp[37]*temp[67]; // +...*
	temp[124] = temp[22]*temp[103]; // *...+
	temp[125] = BH->a3*temp[104]; // +...*
	temp[126] = temp[11]*temp[67]; // +...*
	temp[127] = BH->a*temp[100]; // +...*
	temp[128] = temp[46]*precalculated->nmag7; // +...*
	temp[129] = temp[80]*precalculated->chi_p2; // *...*
	temp[130] = temp[51]*precalculated->chi_p2; // *...*
	temp[131] = temp[130]*temp[113]; // *...-
	temp[132] = temp[120]*precalculated->chi_p2; // *...*
	temp[133] = temp[93]*precalculated->tdot; // *...*
	temp[134] = temp[76]*temp[70]; // +...*
	temp[135] = -BH->a2*temp[104]; // -...*
	temp[136] = temp[104]*temp[24]; // *...*
	temp[137] = -2.*temp[100]*precalculated->rdot; // -...*
	temp[138] = -2.*BH->a*ps->epsilon; // -...*
	temp[139] = 2.*BH->a*temp[105]; // +...*
	temp[140] = BH->a*p->mass; // +...*
	temp[141] = precalculated->delta_p2*sigma_p; // *...+
	temp[142] = ps->sigma6*sigma_p; // *...*
	temp[143] = p->r*precalculated->cosphi; // *...-
	temp[144] = -12.*BH->a6*ps->epsilon; // -...*
	temp[145] = temp[99]*temp[142]; // *...*
	temp[146] = -precalculated->nx*p->r; // -...*
	temp[147] = precalculated->nmag5*precalculated->phidot; // *...*
	temp[148] = (12.*BH->a8*ps->epsilon*precalculated->phidot*ps->sigma12*sigma_p*(precalculated->ny*temp[143]+temp[146]*precalculated->sinphi))/precalculated->nmag; // +...+
	temp[149] = BH->a4*ps->epsilon; // *...*
	temp[150] = p->r4*ps->sigma12; // *...*
	temp[151] = 6.*BH->a7*ps->epsilon; // +...*
	temp[152] = ps->sigma6*precalculated->tdot; // *...*
	temp[153] = 6.*BH->a5*ps->epsilon; // +...*
	temp[154] = p->r2*temp[152]; // *...*
	temp[155] = precalculated->ny*temp[143]+temp[146]*precalculated->sinphi; // (...)
	temp[156] = temp[147]*p->r4; // *...*
	temp[157] = sigma_p*temp[155]; // *...-
	temp[158] = precalculated->phidot*ps->sigma12; // *...*
	temp[159] = -(24.*temp[149]*temp[99]*ps->sigma12*delta_p*sigma_p*temp[155])/precalculated->nmag; // -...-
	temp[160] = -12.*BH->a5*ps->epsilon; // -...*
	temp[161] = ps->sigma6*temp[107]; // *...*
	temp[162] = BH->a3*ps->epsilon; // *...*
	temp[163] = p->r2*temp[161]; // *...*
	temp[164] = (24.*BH->a5*ps->epsilon*ps->sigma12*temp[107]*sigma_p*temp[155])/precalculated->nmag; // +...+
	temp[165] = 12.*BH->a5*ps->epsilon; // +...*
	temp[166] = chi_p*delta_p; // *...*
	temp[167] = 12.*temp[162]; // +...*
	temp[168] = temp[99]*ps->sigma6; // *...*
	temp[169] = BH->a5*ps->epsilon; // *...*
	temp[170] = -(24.*temp[162]*temp[99]*ps->sigma12*temp[166]*sigma_p*temp[155])/precalculated->nmag; // -...+
	temp[171] = 12.*temp[149]*ps->sigma12*temp[95]*delta_p*sigma_p*temp[155]; // (...)
	temp[172] = temp[152]*precalculated->delta_p2; // *...*
	temp[173] = ps->epsilon*precalculated->nmag5; // *...*
	temp[174] = -(12.*temp[162]*ps->sigma12*precalculated->tdot*precalculated->delta_p2*sigma_p*temp[155])/precalculated->nmag; // -...-
	temp[175] = -12.*temp[162]; // -...*
	temp[176] = ps->sigma6*chi_p; // *...*
	temp[177] = -12.*BH->a*temp[173]; // -...*
	temp[178] = precalculated->delta_p2*sigma_p; // *...*
	temp[179] = (24.*BH->a*ps->epsilon*temp[99]*ps->sigma12*chi_p*temp[178]*temp[155])/precalculated->nmag; // +...+
	temp[180] = ps->epsilon*ps->sigma12; // *...*
	temp[181] = -6.*BH->a2*ps->epsilon; // -...*
	temp[182] = precalculated->chi_p2*temp[178]; // *...*
	temp[183] = temp[158]*temp[182]; // *...*
	temp[184] = -6.*temp[173]*ps->sigma6; // -...*
	temp[185] = sigma_p*temp[155]; // *...+
	temp[186] = 6.*ps->epsilon; // +...*
	temp[187] = ps->sigma6*precalculated->chi_p2; // *...*
	temp[188] = 12.*ps->epsilon*temp[158]*precalculated->chi_p2*precalculated->delta_p3*sigma_p*temp[155]; // (...)
	temp[189] = BH->a4*delta_p+2.*BH->a2*p->r2*delta_p+p->r4*delta_p-2.*BH->a3*chi_p*delta_p-2.*BH->a*p->r2*temp[16]+BH->a2*precalculated->chi_p2*delta_p; // (...)
	temp[190] = temp[0]*temp[15]; // +...*
	temp[191] = temp[1]*temp[18]; // +...*
	temp[192] = temp[2]*temp[5]; // +...*
	temp[193] = temp[3]*temp[15]; // +...*
	temp[194] = precalculated->nmag13*temp[18]; // *...*
	temp[195] = temp[33]*precalculated->nmag13; // +...*
	temp[196] = temp[17]*precalculated->rdot; // *...*
	temp[197] = temp[20]*delta_p; // *...+
	temp[198] = precalculated->nmag13*temp[12]; // *...*
	temp[199] = temp[13]*precalculated->nmag13; // +...*
	temp[200] = temp[75]*temp[16]; // *...+
	temp[201] = temp[100]*precalculated->rdot; // *...*
	temp[202] = precalculated->rdot*temp[16]; // *...+
	temp[203] = temp[27]*delta_p; // *...-
	temp[204] = temp[21]*delta_p; // *...+
	temp[205] = temp[22]*delta_p; // *...+
	temp[206] = temp[17]*temp[205]; // *...+
	temp[207] = temp[24]*precalculated->delta_p2; // *...+
	temp[208] = p->mass*precalculated->nmag13; // +...*
	temp[209] = temp[28]*precalculated->nmag13; // +...*
	temp[210] = BH->a2*temp[104]; // +...*
	temp[211] = temp[32]*temp[71]; // +...*
	temp[212] = BH->a4*temp[75]; // +...*
	temp[213] = -2.*BH->a5*temp[201]; // -...*
	temp[214] = temp[18]*precalculated->rdot; // *...*
	temp[215] = 2.*BH->a6*temp[79]; // +...*
	temp[216] = -4.*BH->a6*ps->epsilon; // -...*
	temp[217] = 4.*temp[42]*precalculated->nmag7; // +...*
	temp[218] = precalculated->nmag7*temp[45]; // *...*
	temp[219] = temp[48]*precalculated->nmag7; // +...*
	temp[220] = precalculated->nmag*temp[51]; // *...*
	temp[221] = temp[51]*sigma_p; // *...-
	temp[222] = p->r2*temp[51]; // *...*
	temp[223] = temp[54]*temp[106]; // +...*
	temp[224] = precalculated->nmag*temp[57]; // *...*
	temp[225] = temp[59]*sigma_p; // *...+
	temp[226] = BH->mass*temp[61]; // *...*
	temp[227] = temp[93]*temp[63]; // *...-
	temp[228] = p->mass*temp[122]; // *...*
	temp[229] = precalculated->nmag13*temp[66]; // *...*
	temp[230] = temp[77]*temp[68]; // *...+
	temp[231] = p->r5*precalculated->rdot; // *...*
	temp[232] = temp[67]*temp[70]; // *...*
	temp[233] = -0.5*BH->a5*precalculated->nmag13; // -...*
	temp[234] = precalculated->nmag13*temp[72]; // *...+
	temp[235] = chi_p*sigma_p; // *...-
	temp[236] = temp[75]*temp[235]; // *...-
	temp[237] = chi_p*sigma_p; // *...+
	temp[238] = temp[77]*temp[237]; // *...+
	temp[239] = -2.*BH->a5*temp[79]; // -...*
	temp[240] = 2.*temp[169]*temp[80]; // +...*
	temp[241] = temp[81]*temp[237]; // *...+
	temp[242] = temp[83]*sigma_p; // *...+
	temp[243] = precalculated->nmag*temp[85]; // *...*
	temp[244] = temp[85]*sigma_p; // *...+
	temp[245] = temp[118]*sigma_p; // *...+
	temp[246] = temp[90]*sigma_p; // *...+
	temp[247] = temp[120]*temp[235]; // *...-
	temp[248] = temp[121]*sigma_p; // *...+
	temp[249] = temp[95]*sigma_p; // *...+
	temp[250] = precalculated->nmag13*temp[96]; // *...*
	temp[251] = temp[22]*sigma_p; // *...+
	temp[252] = precalculated->nmag13*temp[98]; // *...+
	temp[253] = BH->a3*temp[100]; // +...*
	temp[254] = 2.*BH->a5*temp[104]; // +...*
	temp[255] = 4.*temp[102]*precalculated->rdot; // +...*
	temp[256] = temp[20]*temp[103]; // *...+
	temp[257] = temp[80]*temp[103]; // *...+
	temp[258] = -2.*BH->a4*temp[105]; // -...*
	temp[259] = temp[48]*temp[106]; // +...*
	temp[260] = p->mass*temp[93]; // *...*
	temp[261] = temp[135]*temp[108]; // +...*
	temp[262] = temp[109]*temp[71]; // *...*
	temp[263] = temp[110]*temp[166]; // +...*
	temp[264] = temp[40]*chi_p; // *...*
	temp[265] = temp[112]*temp[80]; // +...*
	temp[266] = -2.*BH->a*temp[81]; // -...*
	temp[267] = temp[114]*precalculated->nmag7; // +...*
	temp[268] = -4.*BH->a3*temp[105]; // -...*
	temp[269] = BH->a*temp[39]; // *...*
	temp[270] = temp[138]*precalculated->nmag; // +...*
	temp[271] = temp[120]*chi_p; // *...*
	temp[272] = temp[28]*temp[121]; // +...*
	temp[273] = temp[122]*temp[95]; // *...*
	temp[274] = temp[250]*temp[103]; // *...+
	temp[275] = temp[36]*temp[5]; // +...*
	temp[276] = temp[125]*temp[21]; // +...*
	temp[277] = temp[99]*temp[124]; // *...+
	temp[278] = temp[128]*temp[40]; // +...*
	temp[279] = temp[48]*temp[129]; // +...*
	temp[280] = precalculated->nmag*temp[130]; // *...*
	temp[281] = -BH->a2*p->mass; // -...*
	temp[282] = BH->a2*p->mass; // +...*
	temp[283] = precalculated->chi_p2*temp[103]; // *...+
	temp[284] = temp[31]*temp[113]; // *...-
	temp[285] = temp[137]*chi_p; // +...*
	temp[286] = temp[138]*temp[80]; // +...*
	temp[287] = temp[139]*temp[85]; // +...*
	temp[288] = temp[121]*temp[141]; // *...+
	temp[289] = temp[21]*precalculated->delta_p2; // *...*
	temp[290] = temp[147]*temp[142]; // *...*
	temp[291] = precalculated->nmag5*temp[145]; // *...*
	temp[292] = temp[142]*temp[155]; // *...+
	temp[293] = BH->a6*ps->epsilon; // *...*
	temp[294] = (12.*temp[149]*precalculated->phidot*temp[150]*sigma_p*temp[155])/precalculated->nmag; // +...+
	temp[295] = temp[154]*temp[157]; // *...-
	temp[296] = temp[180]*precalculated->tdot; // *...*
	temp[297] = -(12.*temp[169]*p->r2*ps->sigma12*precalculated->tdot*sigma_p*temp[155])/precalculated->nmag; // -...+
	temp[298] = temp[168]*delta_p; // *...*
	temp[299] = 6.*BH->a2*ps->epsilon; // +...*
	temp[300] = delta_p*temp[157]; // *...-
	temp[301] = temp[158]*delta_p; // *...*
	temp[302] = temp[159]-(12.*BH->a2*ps->epsilon*precalculated->phidot*temp[150]*delta_p*sigma_p*temp[155])/precalculated->nmag; // +...+
	temp[303] = (24.*temp[162]*p->r2*ps->sigma12*temp[107]*sigma_p*temp[155])/precalculated->nmag; // +...+
	temp[304] = temp[168]*temp[166]; // *...*
	temp[305] = temp[169]*temp[158]; // *...*
	temp[306] = temp[170]+temp[44]*precalculated->nmag5*ps->sigma6*temp[95]*delta_p*temp[185]; // +...+
	temp[307] = 6.*BH->a*temp[173]; // +...*
	temp[308] = temp[174]-(12.*BH->a*ps->epsilon*p->r2*ps->sigma12*precalculated->tdot*temp[178]*temp[155])/precalculated->nmag; // +...+
	temp[309] = (24.*temp[162]*temp[158]*chi_p*temp[178]*temp[155])/precalculated->nmag; // +...+
	temp[310] = -(24.*BH->a2*temp[180]*temp[95]*temp[178]*temp[155])/precalculated->nmag; // -...+
	temp[311] = temp[184]*temp[95]; // +...*
	temp[312] = (12.*temp[180]*temp[95]*precalculated->delta_p3*sigma_p*temp[155])/precalculated->nmag; // +...+
	temp[313] = -1.*BH->a3-BH->a*p->r2+chi_p*delta_p; // (...)
	temp[314] = temp[190]*precalculated->rdot; // +...+
	temp[315] = temp[192]*p->r6; // +...*
	temp[316] = temp[193]*temp[4]; // +...+
	temp[317] = precalculated->rdot*chi_p; // *...+
	temp[318] = temp[7]*temp[196]; // +...*
	temp[319] = precalculated->nmag13*temp[197]; // *...+
	temp[320] = temp[199]*temp[14]; // +...*
	temp[321] = 4.*BH->a4*temp[201]; // +...*
	temp[322] = temp[194]*temp[202]; // *...+
	temp[323] = -BH->a5*temp[104]; // -...*
	temp[324] = temp[2]*temp[15]; // +...*
	temp[325] = temp[37]*temp[206]; // +...+
	temp[326] = temp[25]*precalculated->rdot; // +...*
	temp[327] = temp[209]*temp[29]; // +...-
	temp[328] = temp[210]*temp[31]; // +...*
	temp[329] = temp[211]*sigma_p; // +...+
	temp[330] = temp[212]*sigma_p; // +...+
	temp[331] = temp[36]*temp[214]; // +...*
	temp[332] = precalculated->nmag13*temp[38]; // *...+
	temp[333] = temp[216]*temp[80]; // +...*
	temp[334] = temp[43]*sigma_p; // *...+
	temp[335] = temp[218]*sigma_p; // *...+
	temp[336] = temp[47]*sigma_p; // *...+
	temp[337] = temp[49]*sigma_p; // *...+
	temp[338] = temp[220]*sigma_p; // *...+
	temp[339] = temp[105]*temp[221]; // *...-
	temp[340] = temp[222]*sigma_p; // *...+
	temp[341] = temp[56]*temp[224]; // +...*
	temp[342] = precalculated->nmag*temp[225]; // *...+
	temp[343] = temp[13]*temp[227]; // +...-
	temp[344] = temp[65]*temp[229]; // +...*
	temp[345] = temp[73]*temp[230]; // *...+
	temp[346] = temp[69]*temp[232]; // +...*
	temp[347] = temp[71]*temp[237]; // *...+
	temp[348] = 2.*BH->a4*temp[97]; // +...*
	temp[349] = -0.5*BH->a3*temp[236]; // -...-
	temp[350] = temp[76]*precalculated->phidot; // +...*
	temp[351] = temp[239]*temp[40]; // +...*
	temp[352] = -2.*BH->a3*temp[241]; // -...+
	temp[353] = temp[84]*temp[243]; // +...*
	temp[354] = temp[52]*temp[244]; // *...+
	temp[355] = temp[89]*precalculated->nmag; // +...*
	temp[356] = -BH->a5*p->mass; // -...*
	temp[357] = BH->a3*temp[228]; // +...*
	temp[358] = temp[37]*temp[250]; // +...*
	temp[359] = temp[97]*temp[251]; // *...+
	temp[360] = temp[123]*temp[99]; // +...*
	temp[361] = temp[254]*temp[101]; // +...*
	temp[362] = 2.*BH->a*temp[104]; // +...*
	temp[363] = temp[149]*temp[257]; // *...+
	temp[364] = temp[258]*temp[51]; // +...*
	temp[365] = -BH->a4*temp[260]; // -...*
	temp[366] = temp[261]*sigma_p; // +...-
	temp[367] = temp[166]*sigma_p; // *...+
	temp[368] = temp[111]*temp[264]; // +...*
	temp[369] = chi_p*temp[113]; // *...+
	temp[370] = chi_p*temp[103]; // *...+
	temp[371] = temp[87]*precalculated->nmag; // +...*
	temp[372] = 2.*temp[269]*temp[118]; // +...*
	temp[373] = BH->a3*p->mass; // +...*
	temp[374] = temp[272]*temp[103]; // +...+
	temp[375] = temp[11]*temp[274]; // +...+
	temp[376] = temp[275]*temp[124]; // +...+
	temp[377] = temp[126]*temp[277]; // +...+
	temp[378] = temp[278]*temp[283]; // +...+
	temp[379] = temp[56]*temp[280]; // +...*
	temp[380] = temp[52]*temp[131]; // *...+
	temp[381] = temp[282]*temp[133]; // +...*
	temp[382] = precalculated->chi_p3*temp[113]; // *...+
	temp[383] = temp[136]*temp[141]; // *...+
	temp[384] = temp[286]*chi_p; // +...*
	temp[385] = temp[287]*temp[141]; // +...+
	temp[386] = BH->a*temp[104]; // +...*
	temp[387] = -6.*BH->a8*ps->epsilon; // -...*
	temp[388] = temp[144]*temp[291]; // +...*
	temp[389] = temp[156]*temp[292]; // *...+
	temp[390] = temp[293]*temp[99]; // *...*
	temp[391] = temp[294]+temp[151]*precalculated->nmag5*temp[152]*temp[185]; // +...+
	temp[392] = temp[296]*sigma_p; // *...*
	temp[393] = temp[297]+6.*temp[293]*temp[147]*ps->sigma6*delta_p*temp[185]; // +...+
	temp[394] = temp[299]*temp[156]; // +...*
	temp[395] = -(12.*temp[293]*temp[301]*sigma_p*temp[155])/precalculated->nmag; // -...+
	temp[396] = temp[175]*precalculated->nmag5; // +...*
	temp[397] = temp[303]+temp[165]*temp[147]*ps->sigma6*temp[166]*temp[185]; // +...+
	temp[398] = temp[305]*temp[166]; // *...*
	temp[399] = temp[306]+temp[171]/precalculated->nmag; // +...+
	temp[400] = temp[172]*temp[185]; // *...+
	temp[401] = precalculated->delta_p2*temp[157]; // *...+
	temp[402] = temp[176]*temp[401]; // *...+
	temp[403] = chi_p*temp[178]; // *...*
	temp[404] = temp[309]+temp[179]; // +...+
	temp[405] = ps->sigma6*temp[95]; // *...*
	temp[406] = temp[181]*temp[147]; // +...*
	temp[407] = (12.*BH->a2*ps->epsilon*temp[183]*temp[155])/precalculated->nmag; // +...+
	temp[408] = temp[147]*temp[187]; // *...*
	temp[409] = precalculated->nmag13*temp[313]; // *...*
	temp[410] = temp[314]+temp[191]*precalculated->rdot; // +...+
	temp[411] = BH->a3*temp[194]; // +...*
	temp[412] = temp[195]*temp[14]; // +...*
	temp[413] = temp[9]*temp[319]; // +...+
	temp[414] = precalculated->rdot*delta_p; // *...+
	temp[415] = temp[321]*temp[16]; // +...+
	temp[416] = temp[19]*precalculated->nmag13; // +...*
	temp[417] = temp[324]*temp[205]; // +...+
	temp[418] = precalculated->nmag13*temp[207]; // *...+
	temp[419] = temp[27]*precalculated->delta_p2; // *...+
	temp[420] = temp[328]*precalculated->delta_p2; // +...+
	temp[421] = temp[330]+temp[213]*sigma_p; // +...+
	temp[422] = temp[215]*temp[40]; // +...*
	temp[423] = temp[217]*temp[334]; // +...+
	temp[424] = temp[128]*temp[336]; // +...+
	temp[425] = temp[50]*temp[338]; // +...+
	temp[426] = -4.*temp[53]*temp[340]; // -...+
	temp[427] = temp[341]*sigma_p; // +...+
	temp[428] = temp[60]*temp[226]; // +...*
	temp[429] = -2.*BH->a4*temp[228]; // -...*
	temp[430] = -BH->a2*temp[345]; // -...+
	temp[431] = temp[231]*temp[68]; // *...+
	temp[432] = temp[233]*temp[347]; // +...+
	temp[433] = temp[348]*temp[74]; // +...*
	temp[434] = -2.*BH->a4*temp[201]; // -...*
	temp[435] = temp[135]*temp[78]; // +...+
	temp[436] = temp[240]*temp[235]; // +...+
	temp[437] = precalculated->nmag7*temp[242]; // *...+
	temp[438] = temp[87]*temp[245]; // +...+
	temp[439] = temp[91]*temp[247]; // +...+
	temp[440] = temp[357]*temp[249]; // +...+
	temp[441] = -BH->a5*temp[359]; // -...+
	temp[442] = temp[360]*temp[251]; // +...+
	temp[443] = temp[361]*sigma_p; // +...+
	temp[444] = temp[362]*temp[256]; // +...+
	temp[445] = temp[58]*temp[218]; // +...*
	temp[446] = temp[259]*temp[113]; // +...+
	temp[447] = temp[366]-0.5*temp[262]*temp[367]; // +...+
	temp[448] = temp[265]*temp[369]; // +...+
	temp[449] = temp[267]*temp[115]; // +...+
	temp[450] = temp[268]*temp[117]; // +...+
	temp[451] = temp[270]*temp[119]; // +...+
	temp[452] = temp[374]+temp[140]*temp[273]*temp[103]; // +...+
	temp[453] = temp[376]+temp[276]*temp[103]; // +...+
	temp[454] = temp[378]+temp[279]*temp[103]; // +...+
	temp[455] = temp[281]*temp[132]; // +...*
	temp[456] = temp[134]*temp[382]; // +...+
	temp[457] = -2.*BH->a2*temp[383]; // -...+
	temp[458] = temp[384]*temp[141]; // +...+
	temp[459] = temp[386]*temp[289]; // +...*
	temp[460] = temp[290]*temp[155]; // *...+
	temp[461] = temp[44]*temp[389]; // +...+
	temp[462] = 24.*temp[390]*ps->sigma12*sigma_p*temp[155]; // (...)
	temp[463] = precalculated->nmag5*temp[295]; // *...-
	temp[464] = temp[393]+12.*temp[149]*precalculated->nmag5*temp[298]*temp[185]; // +...+
	temp[465] = precalculated->nmag5*temp[161]; // *...*
	temp[466] = temp[164]+temp[397]; // +...+
	temp[467] = temp[304]*temp[157]; // *...-
	temp[468] = temp[399]+6.*temp[162]*precalculated->nmag5*temp[400]; // +...+
	temp[469] = temp[147]*temp[402]; // *...+
	temp[470] = temp[403]*temp[155]; // *...+
	temp[471] = BH->a2*temp[173]; // *...*
	temp[472] = temp[310]+temp[406]*ps->sigma6*temp[182]*temp[155]; // +...+
	temp[473] = temp[312]+temp[186]*temp[408]*precalculated->delta_p3*temp[157]; // +...-
	temp[474] = temp[410]+temp[315]*precalculated->rdot; // +...+
	temp[475] = temp[318]*delta_p; // +...+
	temp[476] = temp[11]*temp[198]; // +...*
	temp[477] = BH->a3*temp[200]; // +...+
	temp[478] = BH->a*temp[322]; // +...+
	temp[479] = temp[323]*temp[204]; // +...+
	temp[480] = temp[23]*temp[418]; // +...+
	temp[481] = temp[208]*temp[419]; // +...+
	temp[482] = BH->a*temp[100]; // *...*
	temp[483] = temp[420]+temp[329]; // +...+
	temp[484] = temp[421]+temp[331]*sigma_p; // +...+
	temp[485] = temp[333]*sigma_p; // +...+
	temp[486] = temp[44]*temp[335]; // +...+
	temp[487] = temp[425]+4.*BH->a6*temp[339]; // +...+
	temp[488] = temp[427]+temp[58]*temp[342]; // +...+
	temp[489] = temp[429]*temp[68]; // +...+
	temp[490] = temp[430]+temp[210]*temp[431]; // +...+
	temp[491] = temp[60]*temp[234]; // +...+
	temp[492] = temp[349]+temp[434]*temp[237]; // +...+
	temp[493] = temp[351]*temp[237]; // +...+
	temp[494] = temp[82]*temp[437]; // +...+
	temp[495] = temp[86]*temp[354]; // +...+
	temp[496] = temp[439]+temp[356]*temp[248]; // +...+
	temp[497] = temp[441]+temp[91]*temp[252]; // +...+
	temp[498] = temp[443]+temp[255]*temp[103]; // +...+
	temp[499] = temp[445]*temp[113]; // +...+
	temp[500] = temp[446]+temp[365]*temp[107]*sigma_p; // +...+
	temp[501] = temp[448]+temp[266]*temp[370]; // +...+
	temp[502] = temp[450]+temp[372]*temp[113]; // +...+
	temp[503] = temp[452]+temp[375]; // +...+
	temp[504] = temp[453]+temp[377]; // +...+
	temp[505] = temp[454]+temp[379]*temp[103]; // +...+
	temp[506] = temp[381]*temp[283]; // +...+
	temp[507] = temp[457]+temp[285]*temp[141]; // +...+
	temp[508] = temp[459]*sigma_p; // +...+
	temp[509] = temp[388]*temp[155]; // +...+
	temp[510] = temp[462]/precalculated->nmag; // +...+
	temp[511] = -(12.*BH->a7*temp[392]*temp[155])/precalculated->nmag; // -...+
	temp[512] = temp[302]+temp[160]*temp[465]*temp[157]; // +...+
	temp[513] = temp[167]*precalculated->nmag5; // +...*
	temp[514] = temp[468]+temp[307]*temp[154]*temp[401]; // +...+
	temp[515] = temp[168]*temp[470]; // *...+
	temp[516] = temp[471]*temp[405]; // *...*
	temp[517] = temp[407]+temp[311]*precalculated->delta_p3*temp[185]; // +...+
	temp[518] = temp[409]*temp[189]; // *...*
	temp[519] = temp[474]+temp[316]; // +...+
	temp[520] = temp[412]*delta_p; // +...+
	temp[521] = temp[413]+temp[476]*temp[414]; // +...+
	temp[522] = temp[415]+temp[478]; // +...+
	temp[523] = temp[479]+temp[417]; // +...+
	temp[524] = temp[326]*temp[26]; // +...+
	temp[525] = -2.*temp[482]*temp[30]; // -...+
	temp[526] = temp[484]+temp[37]*temp[332]; // +...+
	temp[527] = temp[423]+temp[486]; // +...+
	temp[528] = temp[487]+temp[426]; // +...+
	temp[529] = temp[488]+temp[428]*temp[68]; // +...+
	temp[530] = temp[490]+temp[346]*temp[237]; // +...+
	temp[531] = temp[492]+temp[350]*temp[238]; // +...+
	temp[532] = temp[352]+temp[494]; // +...+
	temp[533] = temp[495]+temp[438]; // +...+
	temp[534] = temp[496]+temp[440]; // +...+
	temp[535] = temp[497]+temp[442]; // +...+
	temp[536] = temp[498]+temp[444]; // +...+
	temp[537] = temp[499]+temp[364]*temp[103]; // +...+
	temp[538] = temp[368]*temp[103]; // +...+
	temp[539] = temp[371]*temp[116]; // +...+
	temp[540] = temp[373]*temp[271]; // +...*
	temp[541] = temp[123]*precalculated->phidot; // +...*
	temp[542] = temp[505]+temp[58]*temp[380]; // +...+
	temp[543] = temp[456]+temp[135]*temp[284]; // +...+
	temp[544] = temp[140]*temp[288]; // +...+
	temp[545] = temp[509]+temp[461]; // +...+
	temp[546] = temp[391]+temp[153]*temp[463]; // +...+
	temp[547] = ps->sigma6*temp[300]; // *...+
	temp[548] = temp[396]*temp[163]; // +...*
	temp[549] = temp[513]*temp[467]; // +...-
	temp[550] = temp[514]+temp[308]; // +...+
	temp[551] = temp[177]*temp[515]; // +...+
	temp[552] = temp[516]*temp[401]; // *...+
	temp[553] = temp[519]+temp[411]*temp[317]; // +...+
	temp[554] = temp[320]*temp[16]; // +...+
	temp[555] = temp[416]*temp[203]; // +...+
	temp[556] = temp[480]+temp[524]; // +...+
	temp[557] = temp[525]+temp[483]; // +...+
	temp[558] = temp[526]+temp[422]*sigma_p; // +...+
	temp[559] = temp[424]+temp[219]*temp[337]; // +...+
	temp[560] = temp[529]+temp[343]; // +...+
	temp[561] = temp[530]+temp[432]; // +...+
	temp[562] = temp[531]+temp[435]; // +...+
	temp[563] = temp[532]+temp[353]*sigma_p; // +...+
	temp[564] = temp[534]+temp[358]*sigma_p; // +...+
	temp[565] = temp[536]+2.*temp[363]; // +...+
	temp[566] = temp[447]+temp[263]*sigma_p; // +...+
	temp[567] = temp[449]+temp[539]; // +...+
	temp[568] = temp[540]*temp[103]; // +...+
	temp[569] = temp[504]+temp[127]*temp[124]; // +...+
	temp[570] = temp[506]+temp[543]; // +...+
	temp[571] = temp[385]+temp[544]; // +...+
	temp[572] = temp[545]+temp[148]; // +...+
	temp[573] = temp[511]+temp[464]; // +...+
	temp[574] = temp[395]+temp[512]; // +...+
	temp[575] = temp[466]+temp[549]; // +...-
	temp[576] = temp[550]+temp[175]*temp[469]; // +...+
	temp[577] = temp[472]+temp[517]; // +...+
	temp[578] = temp[553]+temp[520]; // +...+
	temp[579] = temp[554]+temp[477]; // +...+
	temp[580] = temp[523]+temp[325]; // +...+
	temp[581] = temp[327]+temp[557]; // +...+
	temp[582] = temp[558]+temp[485]; // +...+
	temp[583] = temp[528]+temp[223]*sigma_p; // +...+
	temp[584] = temp[344]*temp[63]; // +...+
	temp[585] = temp[433]*temp[235]; // +...+
	temp[586] = temp[436]+temp[563]; // +...+
	temp[587] = temp[564]+temp[535]; // +...+
	temp[588] = temp[565]+temp[537]; // +...+
	temp[589] = temp[538]+temp[501]; // +...+
	temp[590] = temp[451]+temp[568]; // +...+
	temp[591] = temp[569]+temp[542]; // +...+
	temp[592] = temp[570]+temp[507]; // +...+
	temp[593] = temp[508]+temp[387]*temp[460]; // +...+
	temp[594] = temp[573]+temp[394]*temp[547]; // +...+
	temp[595] = temp[575]-(24.*temp[398]*sigma_p*temp[155])/precalculated->nmag; // +...+
	temp[596] = temp[578]+temp[475]; // +...+
	temp[597] = temp[522]+temp[555]; // +...+
	temp[598] = temp[481]+temp[581]; // +...+
	temp[599] = temp[582]+temp[527]; // +...+
	temp[600] = temp[560]+temp[489]; // +...+
	temp[601] = temp[491]+temp[585]; // +...+
	temp[602] = temp[586]+temp[533]; // +...+
	temp[603] = temp[587]+temp[253]*temp[251]; // +...+
	temp[604] = temp[589]+temp[567]; // +...+
	temp[605] = temp[503]+temp[541]*temp[124]; // +...+
	temp[606] = temp[592]+temp[458]; // +...+
	temp[607] = temp[572]+temp[510]; // +...+
	temp[608] = temp[574]+temp[548]*temp[185]; // +...+
	temp[609] = temp[404]+12.*temp[552]; // +...+
	temp[610] = temp[596]+temp[521]; // +...+
	temp[611] = temp[580]+temp[556]; // +...+
	temp[612] = temp[599]+temp[559]; // +...+
	temp[613] = temp[584]+temp[561]; // +...+
	temp[614] = temp[493]+temp[602]; // +...+
	temp[615] = temp[603]+temp[588]; // +...+
	temp[616] = temp[604]+temp[502]; // +...+
	temp[617] = temp[591]+temp[455]*temp[103]; // +...+
	temp[618] = temp[607]+temp[546]; // +...+
	temp[619] = temp[595]+temp[576]; // +...+
	temp[620] = temp[577]+temp[473]; // +...-
	temp[621] = temp[610]+temp[579]; // +...+
	temp[622] = temp[598]+temp[195]*temp[34]; // +...+
	temp[623] = temp[613]+temp[601]; // +...+
	temp[624] = temp[355]*temp[246]; // +...+
	temp[625] = temp[566]+temp[616]; // +...+
	temp[626] = temp[617]+temp[606]; // +...+
	temp[627] = temp[618]+temp[594]; // +...+
	temp[628] = temp[551]+temp[609]; // +...+
	temp[629] = temp[621]+temp[597]; // +...+
	temp[630] = temp[612]+temp[583]; // +...+
	temp[631] = temp[562]+temp[614]; // +...+
	temp[632] = temp[500]+temp[625]; // +...+
	temp[633] = temp[626]+temp[571]; // +...+
	temp[634] = temp[608]+temp[619]; // +...+
	temp[635] = temp[629]+temp[611]; // +...+
	temp[636] = temp[600]+temp[623]; // +...+
	temp[637] = temp[615]+temp[632]; // +...+
	temp[638] = temp[633]+temp[593]; // +...+
	temp[639] = temp[628]+temp[620]; // +...-
	temp[640] = temp[635]+temp[622]; // +...+
	temp[641] = temp[631]+temp[624]; // +...+
	temp[642] = temp[605]+temp[638]; // +...+
	temp[643] = temp[640]+temp[630]; // +...+
	temp[644] = temp[637]+temp[590]; // +...+
	temp[645] = temp[634]+temp[639]; // +...-
	temp[646] = temp[643]+temp[636]; // +...+
	temp[647] = temp[642]+temp[627]; // +...+
	temp[648] = temp[646]+temp[641]; // +...+
	temp[649] = temp[648]+temp[644]; // +...+
	temp[650] = temp[649]+temp[647]; // +...+
	temp[651] = temp[650]+temp[645]; // +...-


	phi_new = 2. * p->phi - p->phi_prev + dlambda*dlambda*((2.*(+temp[651]-temp[188]/precalculated->nmag))/(p->mass*temp[518]*sigma_p));
	return phi_new;
}

cpp_dec_float_n calc_theta_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda, cpp_dec_float_n* temp)
{
	/*
	Calculate the next value of theta.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_n sigma_p: Precalculated sigma value.
		- const cpp_dec_float_n delta_p: Precalculated delta value.
		- const cpp_dec_float_n chi_p: Precalculated chi value.
		- const Precalculated* precalculated: Instance of the precalculated variables struct.
		- const cpp_dec_float_n dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_n theta_new: New coordinate.
	*/

	cpp_dec_float_n theta_new{};


	theta_new = 2. * p->theta - p->theta_prev + dlambda*dlambda*(0.);
	return theta_new;
}

std::tuple<cpp_dec_float_n, cpp_dec_float_n,
	cpp_dec_float_n, cpp_dec_float_n> eulerMoveMathematica(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const cpp_dec_float_n dlambda)
{
	/*
	Calculate the next value of theta.

	Inputs:
		- const BlackHole* BH: Instance of the black hole.
		- Particle* p: Instance of the particle for which the next coord is being calculated.
		- Particle* ps: Instance of the source particle.
		- const cpp_dec_float_n sigma_p: Precalculated sigma value.
		- const cpp_dec_float_n delta_p: Precalculated delta value.
		- const cpp_dec_float_n chi_p: Precalculated chi value.
		- const cpp_dec_float_n dlambda: Lambda step.

	Outputs:
		- cpp_dec_float_n t_new: New t coordinate.
		- cpp_dec_float_n r_new: New r coordinate.
		- cpp_dec_float_n phi_new: New phi coordinate.
		- cpp_dec_float_n theta_new: New theta coordinate.
	*/

	cpp_dec_float_n t_new{};
	cpp_dec_float_n r_new{};
	cpp_dec_float_n phi_new{};
	cpp_dec_float_n theta_new{};

	Precalculated precalculated{BH, p, ps, sigma_p, delta_p, chi_p, dlambda};

	t_new = calc_t_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda, precalculated.temp0);

	r_new = calc_r_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda, precalculated.temp1);

	phi_new = calc_phi_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda, precalculated.temp2);

	theta_new = calc_theta_next(BH,
		p, ps,
		sigma_p, delta_p,
		chi_p, &precalculated, dlambda);

	return {t_new, r_new, phi_new, theta_new};
}
