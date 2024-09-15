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
	temp[0] = 2.*BH->a2*delta_p; // (...)
	temp[1] = p->charge*BH->charge; // *...*
	temp[2] = -(4.*BH->a*p->r*precalculated->rdot*(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot))/precalculated->sigma_p2; // -...+
	temp[3] = (4.*p->r*precalculated->rdot*(BH->a*precalculated->phidot-precalculated->tdot)*delta_p)/precalculated->sigma_p2; // +...-
	temp[4] = (4.*ps->epsilon*(-2.*BH->mass+2.*p->r)*precalculated->rdot*ps->sigma6*(-1.*precalculated->nmag6+ps->sigma6))/(p->mass*precalculated->nmag12*sigma_p); // +...+
	temp[5] = (precalculated->rdot*(precalculated->nx*precalculated->cosphi+precalculated->ny*precalculated->sinphi))/precalculated->nmag+(precalculated->phidot*(precalculated->ny*p->r*precalculated->cosphi-precalculated->nx*p->r*precalculated->sinphi))/precalculated->nmag+(precalculated->rsdot*(-(precalculated->nx*precalculated->cosphis)-precalculated->ny*precalculated->sinphis))/precalculated->nmag+(precalculated->phisdot*(-(precalculated->ny*ps->r*precalculated->cosphis)+precalculated->nx*ps->r*precalculated->sinphis))/precalculated->nmag; // (...)
	temp[6] = BH->a2+p->r2-delta_p; // (...)
	temp[7] = (4.*BH->a*p->r*precalculated->rdot*(-1.*BH->a*precalculated->phidot+precalculated->tdot)*delta_p)/precalculated->sigma_p2; // +...-
	temp[8] = precalculated->rdot*ps->sigma6; // *...*
	temp[9] = (2.*BH->a*(-2.*BH->mass+2.*p->r)*precalculated->rdot*(BH->a*precalculated->phidot-precalculated->tdot))/sigma_p; // +...-
	temp[10] = temp[1]*precalculated->rdot; // *...*
	temp[11] = ps->epsilon*ps->sigma6; // *...*
	temp[12] = precalculated->phidot*temp[6]; // *...+
	temp[13] = precalculated->ny*p->r*precalculated->cosphi-precalculated->nx*p->r*precalculated->sinphi; // (...)
	temp[14] = 48.*BH->a*temp[11]*(-1.*precalculated->nmag6+ps->sigma6)*temp[6]*temp[5]; // (...)
	temp[15] = 2.*BH->a*delta_p; // (...)
	temp[16] = 4.*temp[1]*p->r2*precalculated->rdot; // (...)
	temp[17] = temp[2]+(8.*ps->epsilon*p->r*temp[8]*(-1.*precalculated->nmag6+ps->sigma6)*(BH->a2-delta_p))/(p->mass*precalculated->nmag12*precalculated->sigma_p2); // +...+
	temp[18] = -1.*BH->a*precalculated->phidot+precalculated->tdot; // (...)
	temp[19] = 48.*temp[11]*(-1.*precalculated->nmag6+ps->sigma6)*(BH->a2-delta_p)*temp[5]; // (...)
	temp[20] = 4.*p->r*(BH->a2+p->r2)*precalculated->rdot*(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot); // (...)
	temp[21] = -(8.*BH->a*ps->epsilon*p->r*temp[8]*(-1.*precalculated->nmag6+ps->sigma6)*temp[6])/(p->mass*precalculated->nmag12*precalculated->sigma_p2); // -...+
	temp[22] = -1.*precalculated->nmag6+ps->sigma6; // (...)
	temp[23] = temp[9]-(4.*p->r*precalculated->rdot*(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot))/sigma_p; // +...+
	temp[24] = -1.*BH->a2+delta_p; // (...)
	temp[25] = -(24.*BH->a*temp[11]*temp[6]*temp[5])/(p->mass*precalculated->nmag7*sigma_p); // -...-
	temp[26] = (-2.*BH->a2)/sigma_p+(2.*delta_p)/sigma_p; // (...)
	temp[27] = temp[17]+temp[3]; // +...-
	temp[28] = (4.*BH->a*precalculated->phidot*p->r*precalculated->rdot)/sigma_p; // +...+
	temp[29] = (24.*temp[11]*(BH->a2-delta_p)*temp[5])/(p->mass*precalculated->nmag7*sigma_p); // +...+
	temp[30] = -(4.*temp[1]*p->r2*precalculated->rdot*chi_p)/(p->mass*precalculated->sigma_p2); // -...+
	temp[31] = (8.*BH->a*ps->epsilon*(p->r-0.5*(-2.*BH->mass+2.*p->r))*temp[8]*temp[22])/(p->mass*precalculated->nmag12*sigma_p); // +...+
	temp[32] = 2.*BH->a*(BH->a2+p->r2); // (...)
	temp[33] = temp[27]-(2.*temp[1]*precalculated->rdot)/(p->mass*sigma_p); // +...+
	temp[34] = temp[32]/sigma_p-temp[15]/sigma_p; // (...)
	temp[35] = temp[30]+temp[21]; // +...+
	temp[36] = precalculated->phidot*p->r; // *...*
	temp[37] = temp[31]+temp[23]; // +...+
	temp[38] = (48.*temp[11]*(-0.5*precalculated->nmag6+ps->sigma6)*(BH->a*temp[12]+precalculated->tdot*temp[24])*temp[13])/(p->mass*precalculated->nmag14*sigma_p); // +...+
	temp[39] = temp[16]/(p->mass*precalculated->sigma_p2)+temp[33]+temp[28]+temp[4]+(2.*(-2.*BH->mass+2.*p->r)*precalculated->rdot*temp[18])/sigma_p+temp[29]+temp[19]/(p->mass*precalculated->nmag13*sigma_p); // (...)
	temp[40] = -1.*((-2.*pow(BH->a2+p->r2,2))/sigma_p+temp[0]/sigma_p)*temp[39]; // -...+
	temp[41] = -(4.*temp[36]*(BH->a2+p->r2)*precalculated->rdot)/sigma_p; // -...+
	temp[42] = temp[35]+temp[7]; // +...+
	temp[43] = (2.*temp[10]*chi_p)/(p->mass*sigma_p); // +...+
	temp[44] = temp[42]+temp[41]; // +...+
	temp[45] = temp[38]+temp[25]; // +...-
	temp[46] = temp[20]/precalculated->sigma_p2+temp[44]+temp[37]+temp[43]+temp[45]-temp[14]/(p->mass*precalculated->nmag13*sigma_p); // (...)
	temp[47] = +temp[40]+temp[34]*temp[46]; // (...)


	t_new = 2. * p->t - p->t_prev + dlambda*dlambda*(-(temp[47]/(pow(temp[32]/sigma_p-temp[15]/sigma_p,2)-temp[26]*((-2.*pow(BH->a2+p->r2,2))/sigma_p+temp[0]/sigma_p))));
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
	temp[0] = -2.*p->r*(precalculated->rdot*precalculated->rdot); // (...)
	temp[1] = (4.*p->charge*BH->charge*p->r2*(-1.*precalculated->tdot+precalculated->phidot*chi_p))/(p->mass*precalculated->sigma_p2); // +...+
	temp[2] = p->r*ps->sigma6; // *...*
	temp[3] = BH->a*precalculated->phidot*(BH->a2+p->r2-delta_p)+precalculated->tdot*(-1.*BH->a2+delta_p); // (...)
	temp[4] = BH->a*precalculated->phidot; // *...+
	temp[5] = precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot; // (...)
	temp[6] = 2.*BH->a*precalculated->phidot*p->r+(-2.*BH->mass+2.*p->r)*(-1.*temp[4]+precalculated->tdot); // (...)
	temp[7] = ((-2.*BH->mass+2.*p->r)*(precalculated->rdot*precalculated->rdot)*sigma_p)/precalculated->delta_p2; // +...+
	temp[8] = 48.*ps->epsilon*ps->sigma6*(-1.*precalculated->nmag6+ps->sigma6)*temp[3]*(precalculated->nx*precalculated->cosphi+precalculated->ny*precalculated->sinphi); // (...)
	temp[9] = -(2.*p->r*pow(precalculated->phidot*(BH->a2+p->r2)-BH->a*precalculated->tdot,2))/precalculated->sigma_p2; // -...+
	temp[10] = -1.*precalculated->nmag6+ps->sigma6; // (...)
	temp[11] = -((-2.*BH->mass+2.*p->r)*pow(-1.*temp[4]+precalculated->tdot,2))/sigma_p; // -...+
	temp[12] = -(2.*p->charge*BH->charge*(-1.*precalculated->tdot+precalculated->phidot*chi_p))/(p->mass*sigma_p); // -...+
	temp[13] = temp[9]+temp[1]; // +...+
	temp[14] = (8.*ps->epsilon*temp[2]*temp[10]*temp[3])/(p->mass*precalculated->nmag12*precalculated->sigma_p2); // +...+
	temp[15] = temp[12]+temp[7]; // +...+
	temp[16] = ps->sigma6*temp[3]; // *...*
	temp[17] = temp[13]+(2.*p->r*pow(-1.*temp[4]+precalculated->tdot,2)*delta_p)/precalculated->sigma_p2; // +...+
	temp[18] = temp[15]+(24.*ps->epsilon*temp[16]*(precalculated->nx*precalculated->cosphi+precalculated->ny*precalculated->sinphi))/(p->mass*precalculated->nmag8*sigma_p); // +...+
	temp[19] = temp[17]+temp[14]; // +...+
	temp[20] = precalculated->phidot*p->r; // *...*
	temp[21] = ps->epsilon*ps->sigma6; // *...*
	temp[22] = temp[19]+temp[11]; // +...+
	temp[23] = -(4.*temp[21]*temp[10]*temp[6])/(p->mass*precalculated->nmag12*sigma_p); // -...+
	temp[24] = temp[22]+(4.*temp[20]*temp[5])/sigma_p; // +...+
	temp[25] = temp[24]+temp[23]; // +...+
	temp[26] = temp[25]+temp[18]; // +...+


	r_new = 2. * p->r - p->r_prev + dlambda*dlambda*((0.5*delta_p*(temp[0]/delta_p+temp[26]+temp[8]/(p->mass*precalculated->nmag14*sigma_p)))/sigma_p);
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
	temp[0] = BH->a5*precalculated->nmag13; // *...*
	temp[1] = 0.16666666666666666*BH->a3; // +...*
	temp[2] = p->charge*BH->charge; // *...*
	temp[3] = BH->a*precalculated->nmag13; // *...*
	temp[4] = -0.08333333333333333*BH->a4; // -...*
	temp[5] = temp[2]*p->r2; // *...*
	temp[6] = -0.08333333333333333*BH->a2; // -...*
	temp[7] = temp[2]*p->r4; // *...*
	temp[8] = -0.16666666666666666*BH->a3; // -...*
	temp[9] = temp[5]*precalculated->rdot; // *...*
	temp[10] = precalculated->rdot*delta_p; // *...+
	temp[11] = precalculated->phidot*p->r5; // *...*
	temp[12] = -0.08333333333333333*p->mass; // -...*
	temp[13] = precalculated->phidot*p->r7; // *...*
	temp[14] = BH->a2*precalculated->nmag13; // *...*
	temp[15] = 0.08333333333333333*precalculated->nmag13; // +...*
	temp[16] = 0.08333333333333333*temp[3]; // +...*
	temp[17] = p->mass*precalculated->nmag13; // *...*
	temp[18] = -0.08333333333333333*precalculated->nmag13; // -...*
	temp[19] = -0.041666666666666664*temp[0]; // -...*
	temp[20] = -0.08333333333333333*BH->a6; // -...*
	temp[21] = p->r*precalculated->rdot; // *...*
	temp[22] = BH->a4*p->mass; // *...*
	temp[23] = precalculated->nmag13*precalculated->phidot; // *...*
	temp[24] = BH->a3*precalculated->nmag13; // *...*
	temp[25] = temp[17]*precalculated->phidot; // *...*
	temp[26] = BH->a2*p->mass; // *...*
	temp[27] = temp[23]*p->r4; // *...*
	temp[28] = -0.041666666666666664*temp[3]; // -...*
	temp[29] = -0.16666666666666666*BH->a2; // -...*
	temp[30] = -0.16666666666666666*BH->a5; // -...*
	temp[31] = temp[21]*ps->sigma6; // *...*
	temp[32] = ps->epsilon*BH->mass; // *...*
	temp[33] = p->r2*precalculated->rdot; // *...*
	temp[34] = BH->a3*ps->epsilon; // *...*
	temp[35] = p->r3*precalculated->rdot; // *...*
	temp[36] = BH->a*temp[32]; // *...*
	temp[37] = p->r4*precalculated->rdot; // *...*
	temp[38] = -0.16666666666666666*BH->a; // -...*
	temp[39] = p->r5*precalculated->rdot; // *...*
	temp[40] = 0.16666666666666666*BH->a5; // +...*
	temp[41] = temp[21]*ps->sigma12; // *...*
	temp[42] = temp[32]*precalculated->nmag; // *...*
	temp[43] = 0.3333333333333333*temp[34]; // +...*
	temp[44] = ps->sigma12*sigma_p; // *...-
	temp[45] = temp[36]*precalculated->nmag; // *...*
	temp[46] = 0.16666666666666666*BH->a; // +...*
	temp[47] = temp[39]*ps->sigma12; // *...*
	temp[48] = BH->a5*temp[17]; // *...*
	temp[49] = precalculated->tdot*sigma_p; // *...-
	temp[50] = p->mass*BH->mass; // *...*
	temp[51] = temp[33]*precalculated->tdot; // *...*
	temp[52] = temp[35]*temp[49]; // *...-
	temp[53] = BH->a*temp[50]; // *...*
	temp[54] = temp[37]*precalculated->tdot; // *...*
	temp[55] = BH->a*temp[17]; // *...*
	temp[56] = precalculated->tdot*sigma_p; // *...+
	temp[57] = BH->a4*precalculated->nmag13; // *...*
	temp[58] = chi_p*sigma_p; // *...+
	temp[59] = temp[14]*temp[9]; // *...*
	temp[60] = temp[24]*temp[2]; // *...*
	temp[61] = delta_p*sigma_p; // *...+
	temp[62] = BH->a4*temp[25]; // *...*
	temp[63] = temp[6]*temp[50]; // +...*
	temp[64] = temp[16]*temp[9]; // +...*
	temp[65] = BH->a2*temp[25]; // *...*
	temp[66] = 0.16666666666666666*temp[17]; // +...*
	temp[67] = temp[43]*precalculated->nmag7; // +...*
	temp[68] = -0.16666666666666666*temp[36]; // -...*
	temp[69] = ps->sigma6*temp[61]; // *...+
	temp[70] = BH->a*ps->epsilon; // *...*
	temp[71] = temp[35]*ps->sigma6; // *...*
	temp[72] = temp[34]*precalculated->nmag; // *...*
	temp[73] = temp[45]*temp[33]; // *...*
	temp[74] = delta_p*sigma_p; // *...-
	temp[75] = temp[70]*precalculated->nmag; // *...*
	temp[76] = temp[8]*temp[17]; // +...*
	temp[77] = 0.08333333333333333*temp[53]; // +...*
	temp[78] = temp[51]*temp[61]; // *...+
	temp[79] = temp[35]*precalculated->tdot; // *...*
	temp[80] = temp[2]*precalculated->rdot; // *...*
	temp[81] = chi_p*temp[61]; // *...+
	temp[82] = temp[80]*precalculated->delta_p2; // *...*
	temp[83] = precalculated->delta_p2*sigma_p; // *...-
	temp[84] = temp[35]*precalculated->delta_p2; // *...*
	temp[85] = precalculated->nmag7*temp[31]; // *...*
	temp[86] = precalculated->nmag*temp[41]; // *...*
	temp[87] = temp[55]*temp[21]; // *...*
	temp[88] = precalculated->delta_p2*sigma_p; // *...+
	temp[89] = BH->a7*ps->epsilon; // *...*
	temp[90] = precalculated->phidot*ps->sigma6; // *...*
	temp[91] = p->r*precalculated->cosphi; // *...-
	temp[92] = -BH->a5*ps->epsilon; // -...*
	temp[93] = precalculated->phidot*p->r2; // *...*
	temp[94] = -precalculated->nx*p->r; // -...*
	temp[95] = precalculated->nmag5*precalculated->phidot; // *...*
	temp[96] = (BH->a7*ps->epsilon*precalculated->phidot*ps->sigma12*sigma_p*(precalculated->ny*temp[91]+temp[94]*precalculated->sinphi))/precalculated->nmag; // +...+
	temp[97] = ps->epsilon*precalculated->phidot; // *...*
	temp[98] = precalculated->ny*temp[91]+temp[94]*precalculated->sinphi; // (...)
	temp[99] = -(BH->a4*ps->epsilon*p->r2*ps->sigma12*precalculated->tdot*sigma_p*temp[98])/precalculated->nmag; // -...+
	temp[100] = precalculated->nmag5*temp[93]; // *...*
	temp[101] = temp[95]*p->r4; // *...*
	temp[102] = sigma_p*temp[98]; // *...-
	temp[103] = delta_p*sigma_p; // *...*
	temp[104] = -(4.*temp[34]*temp[93]*ps->sigma12*temp[103]*temp[98])/precalculated->nmag; // -...-
	temp[105] = -1.5*BH->a4*ps->epsilon; // -...*
	temp[106] = ps->sigma6*precalculated->tdot; // *...*
	temp[107] = temp[103]*temp[98]; // *...+
	temp[108] = precalculated->tdot*temp[103]; // *...*
	temp[109] = (2.*BH->a2*ps->epsilon*p->r2*ps->sigma12*temp[108]*temp[98])/precalculated->nmag; // +...-
	temp[110] = -BH->a*ps->epsilon; // -...*
	temp[111] = ps->sigma6*precalculated->delta_p2; // *...*
	temp[112] = precalculated->phidot*ps->sigma12; // *...*
	temp[113] = (2.*temp[70]*temp[93]*ps->sigma12*precalculated->delta_p2*sigma_p*temp[98])/precalculated->nmag; // +...+
	temp[114] = 0.5*ps->epsilon; // +...*
	temp[115] = p->r2*temp[106]; // *...*
	temp[116] = -(3.*BH->a2*ps->epsilon*ps->sigma12*precalculated->tdot*precalculated->delta_p2*sigma_p*temp[98])/precalculated->nmag; // -...-
	temp[117] = 0.5*temp[70]*precalculated->nmag5; // +...*
	temp[118] = -(BH->a*temp[97]*ps->sigma12*precalculated->delta_p3*sigma_p*temp[98])/precalculated->nmag; // -...-
	temp[119] = ps->epsilon*ps->sigma12*precalculated->tdot*precalculated->delta_p3*sigma_p*temp[98]; // (...)
	temp[120] = temp[0]*temp[5]; // *...*
	temp[121] = temp[1]*precalculated->nmag13; // +...*
	temp[122] = temp[16]*temp[2]; // +...*
	temp[123] = temp[4]*precalculated->nmag13; // +...*
	temp[124] = precalculated->nmag13*temp[7]; // *...*
	temp[125] = temp[8]*precalculated->nmag13; // +...*
	temp[126] = temp[3]*temp[7]; // *...*
	temp[127] = temp[6]*temp[17]; // +...*
	temp[128] = precalculated->rdot*delta_p; // *...+
	temp[129] = temp[13]*temp[10]; // *...+
	temp[130] = temp[59]*chi_p; // *...*
	temp[131] = temp[15]*temp[7]; // +...*
	temp[132] = chi_p*delta_p; // *...+
	temp[133] = 0.08333333333333333*temp[17]; // +...*
	temp[134] = temp[18]*temp[9]; // +...*
	temp[135] = temp[19]*temp[80]; // +...*
	temp[136] = temp[25]*temp[21]; // *...*
	temp[137] = temp[22]*BH->mass; // *...*
	temp[138] = -0.08333333333333333*temp[24]; // -...*
	temp[139] = temp[62]*temp[35]; // *...*
	temp[140] = temp[26]*BH->mass; // *...*
	temp[141] = temp[28]*temp[7]; // +...*
	temp[142] = temp[29]*temp[17]; // +...*
	temp[143] = temp[30]*ps->epsilon; // +...*
	temp[144] = temp[1]*temp[32]; // +...*
	temp[145] = temp[33]*ps->sigma6; // *...*
	temp[146] = temp[34]*precalculated->nmag7; // *...*
	temp[147] = temp[36]*precalculated->nmag7; // *...*
	temp[148] = temp[38]*ps->epsilon; // +...*
	temp[149] = temp[39]*ps->sigma6; // *...*
	temp[150] = ps->epsilon*temp[86]; // *...*
	temp[151] = temp[42]*temp[33]; // *...*
	temp[152] = temp[43]*precalculated->nmag; // +...*
	temp[153] = temp[45]*temp[37]; // *...*
	temp[154] = temp[46]*ps->epsilon; // +...*
	temp[155] = temp[47]*sigma_p; // *...+
	temp[156] = temp[48]*temp[21]; // *...*
	temp[157] = BH->a3*temp[50]; // *...*
	temp[158] = temp[51]*sigma_p; // *...+
	temp[159] = -0.08333333333333333*temp[53]; // -...*
	temp[160] = temp[54]*sigma_p; // *...+
	temp[161] = temp[55]*temp[39]; // *...*
	temp[162] = temp[57]*temp[80]; // *...*
	temp[163] = temp[59]*temp[58]; // *...+
	temp[164] = temp[60]*precalculated->rdot; // *...*
	temp[165] = temp[62]*temp[21]; // *...*
	temp[166] = temp[23]*temp[33]; // *...*
	temp[167] = 0.4166666666666667*temp[65]; // +...*
	temp[168] = temp[11]*precalculated->rdot; // *...*
	temp[169] = temp[33]*temp[69]; // *...+
	temp[170] = temp[70]*precalculated->nmag7; // *...*
	temp[171] = temp[72]*temp[41]; // *...*
	temp[172] = temp[73]*ps->sigma12; // *...*
	temp[173] = temp[75]*temp[35]; // *...*
	temp[174] = temp[76]*temp[21]; // +...*
	temp[175] = temp[77]*precalculated->nmag13; // +...*
	temp[176] = temp[79]*temp[61]; // *...+
	temp[177] = temp[80]*chi_p; // *...*
	temp[178] = -0.041666666666666664*precalculated->nmag13; // -...*
	temp[179] = temp[6]*temp[136]; // +...*
	temp[180] = temp[25]*temp[84]; // *...*
	temp[181] = temp[85]*temp[88]; // *...+
	temp[182] = 0.08333333333333333*temp[87]; // +...*
	temp[183] = 0.041666666666666664*precalculated->nmag13; // +...*
	temp[184] = temp[90]*sigma_p; // *...*
	temp[185] = temp[92]*temp[100]; // +...*
	temp[186] = -0.5*temp[34]*temp[101]; // -...*
	temp[187] = temp[96]+(2.*BH->a5*ps->epsilon*temp[93]*ps->sigma12*sigma_p*temp[98])/precalculated->nmag; // +...+
	temp[188] = ps->epsilon*precalculated->nmag5; // *...*
	temp[189] = 0.5*BH->a4*temp[188]; // +...*
	temp[190] = -(BH->a6*ps->epsilon*ps->sigma12*precalculated->tdot*sigma_p*temp[98])/precalculated->nmag; // -...+
	temp[191] = 0.5*temp[70]*temp[101]; // +...*
	temp[192] = -(3.*BH->a5*temp[97]*ps->sigma12*temp[103]*temp[98])/precalculated->nmag; // -...+
	temp[193] = temp[105]*precalculated->nmag5; // +...*
	temp[194] = -BH->a2*temp[188]; // -...*
	temp[195] = (3.*BH->a4*ps->epsilon*ps->sigma12*temp[108]*temp[98])/precalculated->nmag; // +...+
	temp[196] = temp[110]*temp[100]; // +...*
	temp[197] = (3.*temp[34]*temp[112]*precalculated->delta_p2*sigma_p*temp[98])/precalculated->nmag; // +...+
	temp[198] = temp[114]*precalculated->nmag5; // +...*
	temp[199] = temp[116]-(ps->epsilon*p->r2*ps->sigma12*precalculated->tdot*precalculated->delta_p2*sigma_p*temp[98])/precalculated->nmag; // +...+
	temp[200] = precalculated->nmag13*p->r4; // *...*
	temp[201] = temp[120]*precalculated->rdot; // *...+
	temp[202] = temp[122]*p->r6; // +...*
	temp[203] = temp[123]*temp[9]; // +...*
	temp[204] = temp[124]*precalculated->rdot; // *...*
	temp[205] = -0.16666666666666666*temp[126]; // -...*
	temp[206] = temp[12]*precalculated->nmag13; // +...*
	temp[207] = temp[130]*delta_p; // *...+
	temp[208] = precalculated->rdot*temp[132]; // *...+
	temp[209] = temp[133]*temp[168]; // +...*
	temp[210] = temp[134]*chi_p; // +...*
	temp[211] = temp[135]*sigma_p; // +...+
	temp[212] = 0.08333333333333333*temp[137]; // +...*
	temp[213] = temp[9]*sigma_p; // *...-
	temp[214] = 0.08333333333333333*temp[140]; // +...*
	temp[215] = temp[141]*precalculated->rdot; // +...*
	temp[216] = temp[143]*temp[85]; // +...*
	temp[217] = precalculated->nmag7*temp[145]; // *...*
	temp[218] = temp[146]*temp[71]; // *...*
	temp[219] = temp[147]*temp[37]; // *...*
	temp[220] = temp[148]*precalculated->nmag7; // +...*
	temp[221] = temp[150]*sigma_p; // *...+
	temp[222] = ps->sigma12*sigma_p; // *...+
	temp[223] = -0.16666666666666666*temp[153]; // -...*
	temp[224] = 0.08333333333333333*temp[156]; // +...*
	temp[225] = temp[157]*precalculated->nmag13; // *...*
	temp[226] = temp[159]*precalculated->nmag13; // +...*
	temp[227] = temp[161]*temp[56]; // *...+
	temp[228] = temp[162]*temp[58]; // *...+
	temp[229] = 0.08333333333333333*temp[164]; // +...*
	temp[230] = temp[165]*temp[61]; // *...+
	temp[231] = temp[64]*temp[61]; // +...+
	temp[232] = temp[35]*temp[61]; // *...+
	temp[233] = temp[67]*temp[31]; // +...*
	temp[234] = temp[68]*precalculated->nmag7; // +...*
	temp[235] = temp[170]*temp[71]; // *...*
	temp[236] = temp[171]*temp[61]; // *...+
	temp[237] = temp[172]*temp[74]; // *...-
	temp[238] = temp[173]*ps->sigma12; // *...*
	temp[239] = precalculated->tdot*temp[61]; // *...+
	temp[240] = temp[17]*temp[176]; // *...+
	temp[241] = temp[177]*temp[74]; // *...+
	temp[242] = temp[28]*temp[82]; // +...*
	temp[243] = -0.16666666666666666*temp[180]; // -...*
	temp[244] = temp[46]*temp[150]; // +...*
	temp[245] = precalculated->tdot*temp[88]; // *...+
	temp[246] = -0.5*temp[89]*precalculated->nmag5; // -...*
	temp[247] = ps->sigma6*temp[102]; // *...+
	temp[248] = sigma_p*temp[98]; // *...+
	temp[249] = (BH->a3*temp[97]*p->r4*ps->sigma12*sigma_p*temp[98])/precalculated->nmag; // +...+
	temp[250] = temp[190]+temp[99]; // +...+
	temp[251] = temp[188]*temp[90]; // *...*
	temp[252] = temp[34]*temp[100]; // *...*
	temp[253] = temp[191]*ps->sigma6; // +...*
	temp[254] = temp[192]+temp[104]; // +...-
	temp[255] = p->r4*ps->sigma12; // *...*
	temp[256] = temp[193]*temp[106]; // +...*
	temp[257] = temp[194]*temp[115]; // +...*
	temp[258] = temp[109]-1.5*temp[34]*precalculated->nmag5*temp[90]*precalculated->delta_p2*temp[102]; // +...+
	temp[259] = precalculated->delta_p2*temp[248]; // *...+
	temp[260] = temp[199]+temp[117]*temp[90]*precalculated->delta_p3*temp[102]; // +...+
	temp[261] = BH->a2+p->r2-delta_p; // (...)
	temp[262] = temp[121]*temp[7]; // +...*
	temp[263] = temp[202]*precalculated->rdot; // +...+
	temp[264] = temp[204]*chi_p; // *...+
	temp[265] = temp[9]*delta_p; // *...+
	temp[266] = temp[127]*temp[11]; // +...*
	temp[267] = 0.16666666666666666*temp[207]; // +...+
	temp[268] = temp[209]*precalculated->delta_p2; // +...+
	temp[269] = temp[211]+temp[20]*temp[136]*sigma_p; // +...+
	temp[270] = -0.25*temp[139]; // -...*
	temp[271] = temp[214]*temp[27]; // +...*
	temp[272] = temp[215]*sigma_p; // +...+
	temp[273] = temp[168]*sigma_p; // *...+
	temp[274] = temp[144]*temp[217]; // +...*
	temp[275] = temp[218]*sigma_p; // *...+
	temp[276] = temp[219]*ps->sigma6; // *...*
	temp[277] = temp[149]*sigma_p; // *...+
	temp[278] = temp[8]*temp[151]; // +...*
	temp[279] = temp[152]*temp[35]; // +...*
	temp[280] = temp[154]*precalculated->nmag; // +...*
	temp[281] = -0.08333333333333333*temp[225]; // -...*
	temp[282] = temp[226]*temp[160]; // +...+
	temp[283] = 0.041666666666666664*temp[228]; // +...+
	temp[284] = temp[229]*temp[61]; // +...+
	temp[285] = temp[63]*temp[166]; // +...*
	temp[286] = temp[167]*temp[232]; // +...+
	temp[287] = temp[233]*delta_p; // +...*
	temp[288] = 0.3333333333333333*temp[235]; // +...*
	temp[289] = 0.16666666666666666*temp[237]; // +...-
	temp[290] = temp[174]*temp[239]; // +...+
	temp[291] = temp[38]*temp[240]; // +...+
	temp[292] = temp[178]*temp[9]; // +...*
	temp[293] = temp[179]*temp[83]; // +...+
	temp[294] = temp[148]*temp[181]; // +...+
	temp[295] = temp[182]*temp[245]; // +...+
	temp[296] = temp[246]*temp[184]; // +...*
	temp[297] = temp[186]*ps->sigma6; // +...*
	temp[298] = temp[249]+0.5*BH->a6*temp[188]*temp[106]*temp[248]; // +...+
	temp[299] = BH->a5*temp[251]; // *...*
	temp[300] = 2.*temp[252]*ps->sigma6; // +...*
	temp[301] = delta_p*temp[102]; // *...+
	temp[302] = -(BH->a*temp[97]*temp[255]*temp[103]*temp[98])/precalculated->nmag; // -...+
	temp[303] = temp[196]*temp[111]; // +...*
	temp[304] = temp[113]+1.5*BH->a2*temp[188]*temp[106]*temp[259]; // +...+
	temp[305] = temp[260]+temp[118]; // +...-
	temp[306] = temp[106]*precalculated->delta_p3; // *...*
	temp[307] = temp[200]*temp[261]; // *...*
	temp[308] = temp[262]*precalculated->rdot; // +...+
	temp[309] = temp[6]*temp[264]; // +...+
	temp[310] = temp[205]*temp[10]; // +...+
	temp[311] = temp[206]*temp[129]; // +...+
	temp[312] = temp[64]*precalculated->delta_p2; // +...+
	temp[313] = temp[269]+temp[212]*temp[166]*sigma_p; // +...+
	temp[314] = temp[271]*precalculated->rdot; // +...*
	temp[315] = temp[216]*sigma_p; // +...+
	temp[316] = -0.3333333333333333*temp[275]; // -...+
	temp[317] = temp[220]*temp[277]; // +...+
	temp[318] = temp[278]*temp[222]; // +...+
	temp[319] = temp[223]*temp[222]; // +...+
	temp[320] = temp[224]*temp[49]; // +...+
	temp[321] = temp[1]*temp[17]; // +...*
	temp[322] = temp[282]+0.08333333333333333*temp[227]; // +...+
	temp[323] = temp[284]+0.16666666666666666*temp[230]; // +...+
	temp[324] = temp[66]*temp[168]; // +...*
	temp[325] = temp[234]*temp[169]; // +...+
	temp[326] = -0.3333333333333333*temp[236]; // -...+
	temp[327] = temp[238]*temp[61]; // *...+
	temp[328] = temp[291]+temp[6]*precalculated->nmag13*temp[241]; // +...+
	temp[329] = temp[243]*sigma_p; // +...+
	temp[330] = temp[244]*temp[88]; // +...+
	temp[331] = temp[177]*temp[83]; // *...+
	temp[332] = temp[185]*temp[247]; // +...+
	temp[333] = temp[187]+temp[298]; // +...+
	temp[334] = temp[250]+1.5*temp[299]*temp[107]; // +...+
	temp[335] = temp[254]+temp[302]; // +...+
	temp[336] = temp[257]*temp[107]; // +...+
	temp[337] = temp[303]*temp[248]; // +...+
	temp[338] = temp[198]*temp[115]; // +...*
	temp[339] = temp[305]-0.5*temp[188]*temp[306]*temp[248]; // +...+
	temp[340] = temp[307]*delta_p; // *...*
	temp[341] = temp[308]+temp[263]; // +...+
	temp[342] = temp[309]+temp[125]*temp[265]; // +...+
	temp[343] = temp[311]+temp[267]; // +...+
	temp[344] = temp[312]+temp[268]; // +...+
	temp[345] = temp[313]+temp[138]*temp[213]; // +...+
	temp[346] = temp[272]+temp[142]*temp[273]; // +...+
	temp[347] = temp[316]+0.16666666666666666*temp[276]*sigma_p; // +...+
	temp[348] = temp[319]+temp[280]*temp[155]; // +...+
	temp[349] = temp[321]*temp[52]; // +...+
	temp[350] = 0.041666666666666664*temp[163]; // +...+
	temp[351] = temp[231]+temp[286]; // +...+
	temp[352] = temp[287]*sigma_p; // +...+
	temp[353] = temp[288]*temp[74]; // +...+
	temp[354] = -0.3333333333333333*temp[327]; // -...+
	temp[355] = temp[328]+temp[292]*temp[81]; // +...+
	temp[356] = temp[329]+temp[294]; // +...+
	temp[357] = temp[183]*temp[331]; // +...+
	temp[358] = temp[332]+temp[297]*temp[248]; // +...+
	temp[359] = temp[334]+temp[300]*temp[107]; // +...+
	temp[360] = temp[256]*temp[301]; // +...+
	temp[361] = temp[258]+temp[337]; // +...+
	temp[362] = temp[338]*precalculated->delta_p2; // +...*
	temp[363] = temp[341]+temp[203]*chi_p; // +...+
	temp[364] = temp[266]*temp[128]; // +...+
	temp[365] = temp[344]+temp[210]*precalculated->delta_p2; // +...+
	temp[366] = temp[346]+temp[315]; // +...+
	temp[367] = temp[347]+temp[317]; // +...+
	temp[368] = temp[318]+temp[279]*temp[44]; // +...+
	temp[369] = temp[349]+temp[322]; // +...+
	temp[370] = temp[323]+temp[285]*temp[61]; // +...+
	temp[371] = temp[352]+temp[325]; // +...+
	temp[372] = temp[289]+temp[354]; // +...+
	temp[373] = temp[355]+temp[242]*sigma_p; // +...+
	temp[374] = temp[330]+temp[295]; // +...+
	temp[375] = temp[358]+temp[333]; // +...+
	temp[376] = temp[359]+temp[253]*temp[301]; // +...+
	temp[377] = temp[195]+temp[361]; // +...+
	temp[378] = temp[362]*temp[102]; // +...+
	temp[379] = temp[363]+temp[342]; // +...+
	temp[380] = temp[343]+temp[131]*temp[208]; // +...+
	temp[381] = temp[314]*sigma_p; // +...+
	temp[382] = temp[274]*sigma_p; // +...+
	temp[383] = temp[40]*temp[221]; // +...+
	temp[384] = temp[320]+temp[281]*temp[158]; // +...+
	temp[385] = temp[370]+temp[351]; // +...+
	temp[386] = temp[371]+temp[353]; // +...+
	temp[387] = temp[290]+temp[175]*temp[78]; // +...+
	temp[388] = temp[374]+temp[357]; // +...+
	temp[389] = temp[375]+temp[189]*temp[115]*temp[102]; // +...+
	temp[390] = temp[377]+temp[197]; // +...+
	temp[391] = temp[379]+temp[310]; // +...+
	temp[392] = temp[365]+temp[345]; // +...+
	temp[393] = temp[381]+temp[366]; // +...+
	temp[394] = temp[383]+temp[368]; // +...+
	temp[395] = temp[369]+temp[283]; // +...+
	temp[396] = temp[324]*temp[61]; // +...+
	temp[397] = temp[372]+temp[387]; // +...+
	temp[398] = temp[356]+temp[388]; // +...+
	temp[399] = temp[389]+temp[376]; // +...+
	temp[400] = temp[336]+temp[390]; // +...+
	temp[401] = temp[391]+temp[364]; // +...+
	temp[402] = temp[270]*sigma_p; // +...+
	temp[403] = temp[382]+temp[367]; // +...+
	temp[404] = temp[384]+temp[395]; // +...+
	temp[405] = temp[396]+temp[386]; // +...+
	temp[406] = temp[373]+temp[293]; // +...+
	temp[407] = temp[399]+temp[335]; // +...+
	temp[408] = temp[304]+temp[378]; // +...+
	temp[409] = temp[401]+temp[380]; // +...+
	temp[410] = temp[393]+temp[403]; // +...+
	temp[411] = temp[404]+temp[350]; // +...+
	temp[412] = temp[326]+temp[397]; // +...+
	temp[413] = temp[296]*temp[98]; // +...+
	temp[414] = temp[400]+temp[408]; // +...+
	temp[415] = temp[409]+temp[392]; // +...+
	temp[416] = temp[394]+temp[348]; // +...+
	temp[417] = temp[405]+temp[412]; // +...+
	temp[418] = temp[413]+temp[407]; // +...+
	temp[419] = temp[415]+temp[402]; // +...+
	temp[420] = temp[411]+temp[385]; // +...+
	temp[421] = temp[398]+temp[418]; // +...+
	temp[422] = temp[419]+temp[410]; // +...+
	temp[423] = temp[417]+temp[406]; // +...+
	temp[424] = temp[414]+temp[339]; // +...+
	temp[425] = temp[422]+temp[416]; // +...+
	temp[426] = temp[421]+temp[360]; // +...+
	temp[427] = temp[425]+temp[420]; // +...+
	temp[428] = temp[427]+temp[423]; // +...+
	temp[429] = temp[428]+temp[426]; // +...+
	temp[430] = temp[429]+temp[424]; // +...+


	phi_new = 2. * p->phi - p->phi_prev + dlambda*dlambda*((-24.*(0.08333333333333333*temp[201]+temp[430]+temp[119]/precalculated->nmag))/(p->mass*temp[340]*sigma_p));
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
